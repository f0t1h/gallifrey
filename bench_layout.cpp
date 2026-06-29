#include <CLI/CLI.hpp>
#include <inttypes.h>
#include <stdint.h>

#include <algorithm>
#include <chrono>
#include <cstdio>
#include <cstring>
#include <numeric>
#include <random>
#include <string>
#include <vector>

#include "asm.h"
#include "bench_layout.h"
#include "ext/lcptools/lcptools_ho.h"
#include "lcmer.h"
#include "syncmer.h"

namespace {

// --- Minimal FASTA reader (same shape as dump_cores.cpp) ---
struct FastaReader {
    FILE *fp = nullptr;
    std::string name, seq;
    explicit FastaReader(const char *path) { fp = fopen(path, "r"); }
    ~FastaReader() {
        if (fp) fclose(fp);
    }
    bool next() {
        if (!fp) return false;
        name.clear();
        seq.clear();
        int c;
        while ((c = fgetc(fp)) != EOF && c != '>') {
        }
        if (c == EOF) return false;
        while ((c = fgetc(fp)) != EOF && c != '\n' && c != '\r')
            name.push_back((char)c);
        while ((c = fgetc(fp)) != EOF) {
            if (c == '>') {
                ungetc(c, fp);
                break;
            }
            if (c != '\n' && c != '\r' && c != ' ') seq.push_back((char)c);
        }
        return !seq.empty();
    }
};

inline void ascii_to_bam_nibbles(const std::string &seq,
                                 std::vector<uint8_t> &out) {
    size_t len = seq.size();
    out.assign((len + 1) / 2, 0);
    for (size_t i = 0; i < len; ++i) {
        uint8_t nib = seq_nt16_table[(uint8_t)seq[i]];
        if (i & 1)
            out[i / 2] |= nib;
        else
            out[i / 2] = nib << 4;
    }
}

inline uint64_t code2(char c) {
    switch (c) {
        case 'c':
        case 'C':
            return 1;
        case 'g':
        case 'G':
            return 2;
        case 'u':
        case 'U':
        case 't':
        case 'T':
            return 3;
        default:
            return 0;  // A and anything else
    }
}

// 112 B inline union = 14 words = 448 bases at 2 bits/base.
constexpr size_t SSO_WORDS = 14;
constexpr uint32_t SSO_BASES = SSO_WORDS * 32;

// Layout A: today's design — a pointer to a separately heap-allocated buffer.
struct SubA {
    uint32_t len;
    uint8_t skip;
    uint64_t origin;
    uint64_t *data;
};

// Layout B: SSO-128 — data inline when it fits, else a heap pointer.
struct SubB {
    uint32_t len;
    uint8_t skip;
    bool inl;
    uint64_t origin;
    union {
        uint64_t *heap;
        uint64_t words[SSO_WORDS];
    } u;
    const uint64_t *words() const { return inl ? u.words : u.heap; }
};
static_assert(sizeof(SubB) == 128, "SubB must be one cache-line pair (128 B)");

// Layout C: arena — struct holds an 8-byte offset into one contiguous blob.
struct SubC {
    uint32_t len;
    uint8_t skip;
    uint64_t origin;
    uint64_t offset;
};

// Layout D: arena with a 4-byte (word) offset, packed to 16 B so the metadata
// vector is denser. skip normalized to 0. offset addresses uint64_t words, so
// uint32_t spans 32 GB of packed sequence — fine per arena chunk.
struct SubD {
    uint64_t origin;
    uint32_t len;
    uint32_t offset;
};
static_assert(sizeof(SubD) == 16, "SubD must be 16 B");

// The representative hot inner loop: read every base's 2-bit code (the access
// pattern of convert_to_ascii / reverse_complement / syncmer scanning) and mix
// into an accumulator so the compiler can't elide the loads.
inline uint64_t read_all_bases(const uint64_t *w, uint32_t len, uint8_t skip,
                               uint64_t acc) {
    const size_t base_bit = (size_t)skip * 2;
    for (uint32_t i = 0; i < len; ++i) {
        const size_t bit = base_bit + (size_t)i * 2;
        acc = acc * 1099511628211ull + ((w[bit >> 6] >> (bit & 63)) & 3);
    }
    return acc;
}

template <class Get>
double time_sweep(size_t n, const std::vector<uint32_t> &order, int reps,
                  Get get, uint64_t &sink) {
    double best = 1e300;
    for (int r = 0; r < reps; ++r) {
        auto t0 = std::chrono::steady_clock::now();
        uint64_t acc = 1469598103934665603ull;
        for (size_t k = 0; k < n; ++k) {
            const uint32_t idx = order[k];
            const uint64_t *w;
            uint32_t len;
            uint8_t skip;
            get(idx, w, len, skip);
            acc = read_all_bases(w, len, skip, acc);
        }
        auto t1 = std::chrono::steady_clock::now();
        sink += acc;
        best = std::min(
            best, std::chrono::duration<double>(t1 - t0).count());
    }
    return best;
}

}  // namespace

static int run_bench_layout(const std::string &input, int syncmer_k,
                            int submer_s, int submer_t, int lcp_level,
                            uint64_t seed, uint64_t max_reads, int reps) {
    syncmer_params spara;
    spara.k = syncmer_k;
    spara.mincode.s = submer_s;
    spara.mincode.t = submer_t;
    using SyncmerType = syncmer::GenericSyncmerRule<uint32_t, syncmer::MurMur>;
    SyncmerType syncmer_functor(spara, (unsigned)seed);

    FastaReader reader(input.c_str());
    if (!reader.fp) {
        fprintf(stderr, "Cannot open input: %s\n", input.c_str());
        return 1;
    }

    // --- Extraction: pack every core subseq into one master blob (offset+len).
    std::vector<uint64_t> master;     // contiguous packed words (becomes arena)
    std::vector<SubC> meta;           // len + offset into master
    std::vector<uint8_t> nibbles;
    lps str = {0, 0, 0, nullptr};
    charsi chs;
    pos_scratch pos_buf;
    uint64_t reads = 0;
    unsigned __int128 total_bases = 0;

    while (reader.next() && reads < max_reads) {
        const std::string &s = reader.seq;
        const size_t len = s.size();
        ascii_to_bam_nibbles(s, nibbles);
        lcp_read lr;
        parse_cores_from_bam_record(nibbles.data(), len, lcp_level, &str,
                                    syncmer_functor, chs, lr, pos_buf);
        uint32_t prev_end = 0;
        for (size_t i = 0; i < lr.cores.size(); ++i) {
            const uint32_t cur_end = lr.pos[i].second;
            const uint32_t sub_len = cur_end - prev_end;
            const uint32_t from = prev_end;
            prev_end = cur_end;
            if (sub_len == 0) continue;
            const uint32_t words = sub_len / 32 + 1;
            const uint64_t off = master.size();
            master.resize(off + words, 0);
            for (uint32_t b = 0; b < sub_len; ++b) {
                const uint64_t code = code2(s[from + b]);
                const size_t bit = (size_t)b * 2;
                master[off + (bit >> 6)] |= code << (bit & 63);
            }
            meta.push_back(SubC{sub_len, 0, /*origin*/ reads, off});
            total_bases += sub_len;
        }
        ++reads;
    }
    free_lps(&str);

    const size_t n = meta.size();
    if (n == 0) {
        fprintf(stderr, "No subseqs extracted.\n");
        return 1;
    }

    // --- Build Layout A (scattered heap), with churn so the allocator can't
    // hand back one contiguous run (approximates the real, long-running,
    // interleaved allocation pattern).
    std::mt19937_64 rng(seed ? seed : 0xC0FFEE);
    std::vector<SubA> A(n);
    std::vector<uint8_t *> junk;
    junk.reserve(n / 2 + 1);
    for (size_t i = 0; i < n; ++i) {
        const uint32_t words = meta[i].len / 32 + 1;
        A[i].len = meta[i].len;
        A[i].skip = 0;
        A[i].origin = meta[i].origin;
        A[i].data = new uint64_t[words];
        std::memcpy(A[i].data, master.data() + meta[i].offset,
                    words * sizeof(uint64_t));
        // churn: allocate a small block; keep ~half (the rest become holes)
        uint8_t *j = new uint8_t[16 + (rng() & 0xFF)];
        if (rng() & 1) {
            junk.push_back(j);
        } else {
            delete[] j;
        }
    }

    // --- Build Layout B (SSO-128).
    std::vector<SubB> B(n);
    for (size_t i = 0; i < n; ++i) {
        const uint32_t words = meta[i].len / 32 + 1;
        B[i].len = meta[i].len;
        B[i].skip = 0;
        B[i].origin = meta[i].origin;
        if (meta[i].len <= SSO_BASES) {
            B[i].inl = true;
            std::memset(B[i].u.words, 0, sizeof(B[i].u.words));
            std::memcpy(B[i].u.words, master.data() + meta[i].offset,
                        words * sizeof(uint64_t));
        } else {
            B[i].inl = false;
            B[i].u.heap = new uint64_t[words];
            std::memcpy(B[i].u.heap, master.data() + meta[i].offset,
                        words * sizeof(uint64_t));
        }
    }

    // --- Layout C (arena): master is already the contiguous blob.
    std::vector<SubC> &C = meta;
    const uint64_t *arena = master.data();

    // --- Layout D (arena, 4-byte offset): same blob, denser metadata vector.
    std::vector<SubD> D(n);
    for (size_t i = 0; i < n; ++i) {
        D[i].origin = meta[i].origin;
        D[i].len = meta[i].len;
        D[i].offset = (uint32_t)meta[i].offset;
    }

    // --- Orders.
    std::vector<uint32_t> seq_order(n);
    std::iota(seq_order.begin(), seq_order.end(), 0u);
    std::vector<uint32_t> shuf_order = seq_order;
    std::shuffle(shuf_order.begin(), shuf_order.end(), rng);

    uint64_t sink = 0;
    const double GB = (double)(uint64_t)total_bases;  // bases (for ns/base)

    auto getA = [&](uint32_t i, const uint64_t *&w, uint32_t &l, uint8_t &s) {
        w = A[i].data;
        l = A[i].len;
        s = A[i].skip;
    };
    auto getB = [&](uint32_t i, const uint64_t *&w, uint32_t &l, uint8_t &s) {
        w = B[i].words();
        l = B[i].len;
        s = B[i].skip;
    };
    auto getC = [&](uint32_t i, const uint64_t *&w, uint32_t &l, uint8_t &s) {
        w = arena + C[i].offset;
        l = C[i].len;
        s = C[i].skip;
    };
    auto getD = [&](uint32_t i, const uint64_t *&w, uint32_t &l, uint8_t &s) {
        w = arena + D[i].offset;
        l = D[i].len;
        s = 0;
    };

    struct Row {
        const char *name;
        double seqt, shuft;
    };
    Row rows[4];
    rows[0] = {"A scattered-heap (24B + ptr)",
               time_sweep(n, seq_order, reps, getA, sink),
               time_sweep(n, shuf_order, reps, getA, sink)};
    rows[1] = {"B SSO-128 (inline<=448b)",
               time_sweep(n, seq_order, reps, getB, sink),
               time_sweep(n, shuf_order, reps, getB, sink)};
    rows[2] = {"C arena (8B offset, 24B)",
               time_sweep(n, seq_order, reps, getC, sink),
               time_sweep(n, shuf_order, reps, getC, sink)};
    rows[3] = {"D arena (4B offset, 16B)",
               time_sweep(n, seq_order, reps, getD, sink),
               time_sweep(n, shuf_order, reps, getD, sink)};

    // --- Report.
    fprintf(stderr,
            "\n=== layout access bench (k=%d s=%d t=%d l=%d) ===\n", syncmer_k,
            submer_s, submer_t, lcp_level);
    fprintf(stderr, "reads=%" PRIu64 "  subseqs=%zu  total_bases=%.0f\n", reads,
            n, GB);
    fprintf(stderr, "struct sizes: A=%zu B=%zu C=%zu bytes\n", sizeof(SubA),
            sizeof(SubB), sizeof(SubC));
    fprintf(stderr,
            "vector footprint: A=%.1f MB  B=%.1f MB  C=%.1f MB(+arena %.1f "
            "MB)\n",
            sizeof(SubA) * n / 1e6, sizeof(SubB) * n / 1e6, sizeof(SubC) * n / 1e6,
            master.size() * 8 / 1e6);
    fprintf(stderr, "%-32s  %12s  %12s\n", "layout", "seq ns/base",
            "shuffled ns/base");
    for (auto &r : rows) {
        fprintf(stderr, "%-32s  %12.3f  %12.3f\n", r.name, r.seqt / GB * 1e9,
                r.shuft / GB * 1e9);
    }
    fprintf(stderr, "(checksum %llu)\n", (unsigned long long)sink);

    // cleanup
    for (auto &a : A) delete[] a.data;
    for (auto &b : B)
        if (!b.inl) delete[] b.u.heap;
    for (auto *j : junk) delete[] j;
    return 0;
}

BenchLayoutProgram::BenchLayoutProgram(CLI::App &parent) {
    auto &app = *parent.add_subcommand(
        "bench-layout", "Benchmark CrumbledSubSeq64 storage layouts");
    app.add_option("input", input, "Input FASTA file")->required();
    app.add_option("--syncmer-k", syncmer_k, "K-mer length")
        ->capture_default_str();
    app.add_option("--submer-s", submer_s, "Submer size")->capture_default_str();
    app.add_option("--submer-t", submer_t, "Submer pos")->capture_default_str();
    app.add_option("-l,--level", lcp_level, "LCP level")->capture_default_str();
    app.add_option("--seed", seed, "Random seed")->capture_default_str();
    app.add_option("--max-reads", max_reads, "Cap reads processed")
        ->capture_default_str();
    app.add_option("--reps", reps, "Timing repetitions (min taken)")
        ->capture_default_str();
    app.callback([this]() {
        return run_bench_layout(input, syncmer_k, submer_s, submer_t, lcp_level,
                                seed, max_reads, reps);
    });
}
