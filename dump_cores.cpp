#include <CLI/CLI.hpp>
#include <inttypes.h>
#include <stdint.h>

#include <algorithm>
#include <cstdio>
#include <map>
#include <string>
#include <vector>

#include "asm.h"
#include "dump_cores.h"
#include "ext/lcptools/lcptools_ho.h"
#include "lcmer.h"
#include "syncmer.h"

struct FastaReader {
    FILE *fp = nullptr;
    std::string name;
    std::string seq;

    explicit FastaReader(const char *path) { fp = fopen(path, "r"); }
    ~FastaReader() {
        if (fp) fclose(fp);
    }
    bool next() {
        if (!fp) return false;
        name.clear();
        seq.clear();
        int c;
        while ((c = fgetc(fp)) != EOF) {
            if (c == '>') break;
        }
        if (c == EOF) return false;
        while ((c = fgetc(fp)) != EOF && c != '\n' && c != '\r') {
            name.push_back(static_cast<char>(c));
        }
        while ((c = fgetc(fp)) != EOF) {
            if (c == '>') {
                ungetc(c, fp);
                break;
            }
            if (c != '\n' && c != '\r' && c != ' ') {
                seq.push_back(static_cast<char>(c));
            }
        }
        return !seq.empty();
    }
};

static void ascii_to_bam_nibbles(const std::string &seq,
                                 std::vector<uint8_t> &out) {
    size_t len = seq.size();
    out.resize((len + 1) / 2);
    for (size_t i = 0; i < len; ++i) {
        uint8_t nibble = seq_nt16_table[static_cast<uint8_t>(seq[i])];
        if (i & 1) {
            out[i / 2] |= nibble;
        } else {
            out[i / 2] = nibble << 4;
        }
    }
}

static int run_dump_cores(const std::string &input_path,
                          const std::string &output_path, int syncmer_k,
                          int submer_s, int submer_t, int lcp_level,
                          uint64_t seed) {
    syncmer_params spara;
    spara.k = syncmer_k;
    spara.mincode.s = submer_s;
    spara.mincode.t = submer_t;

    using SYNCMER_HASHER = syncmer::MurMur;
    using SyncmerType = syncmer::GenericSyncmerRule<uint32_t, SYNCMER_HASHER>;
    SyncmerType syncmer_functor(spara, static_cast<unsigned int>(seed));

    FILE *out = fopen(output_path.c_str(), "w");
    if (!out) {
        fprintf(stderr, "Cannot open output: %s\n", output_path.c_str());
        return 1;
    }
    fprintf(out,
            "read_name\tread_len\tcore_idx\tcore_label\tstart\tend\tspan\t"
            "subseq_len\trle_subseq_len\n");

    // subseq_len = end[i] - end[i-1] (with end[-1] = 0) is the length of the
    // CrumbledSubSeq64 actually stored per core in the assembler
    // (build_lr_cores_with_seqs extracts the non-overlapping portion
    // [prev_end, cur_end)).  That is the quantity that drives the SSO inline
    // capacity decision, not the core span (end - start).
    std::map<uint64_t, uint64_t> hist;      // subseq_len (original bases)
    std::map<uint64_t, uint64_t> hist_rle;  // subseq_len if RLE-compressed
    uint64_t total_subseqs = 0;

    FastaReader reader(input_path.c_str());
    if (!reader.fp) {
        fprintf(stderr, "Cannot open input: %s\n", input_path.c_str());
        fclose(out);
        return 1;
    }

    std::vector<uint8_t> nibbles;
    lps str = {0, 0, 0, nullptr};
    charsi chs;
    pos_scratch pos_buf;
    int read_count = 0;

    while (reader.next()) {
        size_t len = reader.seq.size();
        ascii_to_bam_nibbles(reader.seq, nibbles);

        lcp_read lr;
        parse_cores_from_bam_record(nibbles.data(), len, lcp_level, &str,
                                    syncmer_functor, chs, lr, pos_buf);

        const std::string &s = reader.seq;
        uint32_t prev_end = 0;
        for (size_t i = 0; i < lr.cores.size(); ++i) {
            const uint32_t cur_end = lr.pos[i].second;
            const uint32_t subseq_len = cur_end - prev_end;
            // RLE-compressed length of the same span [prev_end, cur_end):
            // count global homopolymer-collapsed bases (compare across the
            // span boundary, matching the whole-read RLE the pipeline uses).
            uint32_t rle_len = 0;
            for (uint32_t j = prev_end; j < cur_end && j < s.size(); ++j) {
                if (j == 0 || s[j] != s[j - 1]) ++rle_len;
            }
            prev_end = cur_end;
            fprintf(out, "%s\t%zu\t%zu\t%" PRIu32 "\t%u\t%u\t%u\t%u\t%u\n",
                    reader.name.c_str(), len, i, lr.cores[i], lr.pos[i].first,
                    cur_end, cur_end - lr.pos[i].first, subseq_len, rle_len);
            ++hist[subseq_len];
            ++hist_rle[rle_len];
            ++total_subseqs;
        }
        ++read_count;
    }

    free_lps(&str);
    fclose(out);

    // Distribution summary of the stored CrumbledSubSeq64 lengths. Reported
    // both as-is (original bases, what the pipeline stores today) and as it
    // would be if the stored sequence were RLE-compressed.
    auto report = [&](const char *label, const std::map<uint64_t, uint64_t> &h) {
        if (total_subseqs == 0) return;
        std::vector<std::pair<uint64_t, uint64_t>> sorted(h.begin(), h.end());
        std::sort(sorted.begin(), sorted.end());

        auto percentile = [&](double p) -> uint64_t {
            const uint64_t target =
                static_cast<uint64_t>(p * (total_subseqs - 1));
            uint64_t cum = 0;
            for (const auto &[len, cnt] : sorted) {
                cum += cnt;
                if (cum > target) return len;
            }
            return sorted.back().first;
        };
        auto frac_le = [&](uint64_t bound) -> double {
            uint64_t cum = 0;
            for (const auto &[len, cnt] : sorted) {
                if (len <= bound) cum += cnt;
            }
            return 100.0 * static_cast<double>(cum) /
                   static_cast<double>(total_subseqs);
        };

        unsigned __int128 sum = 0;
        for (const auto &[len, cnt] : sorted) sum += (unsigned __int128)len * cnt;
        const double mean =
            static_cast<double>((double)(uint64_t)sum) / (double)total_subseqs;

        fprintf(stderr, "\n=== %s (k=%d s=%d t=%d l=%d) ===\n", label, syncmer_k,
                submer_s, submer_t, lcp_level);
        fprintf(stderr, "count=%" PRIu64 "  mean=%.2f  min=%" PRIu64
                        "  max=%" PRIu64 "\n",
                total_subseqs, mean, sorted.front().first, sorted.back().first);
        fprintf(stderr,
                "p50=%" PRIu64 "  p90=%" PRIu64 "  p95=%" PRIu64
                "  p99=%" PRIu64 "\n",
                percentile(0.50), percentile(0.90), percentile(0.95),
                percentile(0.99));
        fprintf(stderr,
                "<=16:%.2f%%  <=32:%.2f%%  <=48:%.2f%%  <=64:%.2f%%  "
                "<=96:%.2f%%\n",
                frac_le(16), frac_le(32), frac_le(48), frac_le(64),
                frac_le(96));
    };
    report("CrumbledSubSeq64 length (original bases)", hist);
    report("CrumbledSubSeq64 length (if RLE-compressed)", hist_rle);

    // SSO economics. Net memory delta vs the current 24 B struct + per-object
    // heap allocation, accounting for BOTH the struct growth (paid on every
    // instance) and the heap allocation saved on inline-eligible instances.
    //   inline_buf_bytes: bytes available for inline 2-bit data
    //   struct_growth   : new_struct_size - 24 (paid on every instance)
    //   inline_bases    = inline_buf_bytes * 4  (8 B -> 32 bases)
    //   heap_now(len)   = storage_size(len)*8 [+ malloc overhead]
    //   per instance: net = struct_growth - (len<=inline_bases ? heap_now : 0)
    // Negative net = memory saving.
    auto sso_report = [&](const char *tag, uint64_t inline_buf_bytes,
                          int64_t struct_growth) {
        if (total_subseqs == 0) return;
        const uint64_t inline_bases = inline_buf_bytes * 4;
        unsigned __int128 saved_payload = 0, saved_real = 0;
        uint64_t inline_count = 0;
        for (const auto &[len, cnt] : hist) {
            if (len > inline_bases) continue;
            const uint64_t payload = (len / 32 + 1) * 8;
            const uint64_t real = std::max<uint64_t>(payload, 16) + 16;
            inline_count += cnt;
            saved_payload += (unsigned __int128)payload * cnt;
            saved_real += (unsigned __int128)real * cnt;
        }
        const double frac = 100.0 * (double)inline_count / (double)total_subseqs;
        // net = struct_growth*N - saved_heap, averaged per instance
        const double net_payload =
            (double)struct_growth -
            (double)(uint64_t)saved_payload / (double)total_subseqs;
        const double net_real =
            (double)struct_growth -
            (double)(uint64_t)saved_real / (double)total_subseqs;
        fprintf(stderr, "\n--- %s (k=%d s=%d t=%d l=%d) ---\n", tag, syncmer_k,
                submer_s, submer_t, lcp_level);
        fprintf(stderr,
                "inline-eligible (len<=%" PRIu64 "): %" PRIu64 " / %" PRIu64
                " = %.3f%%   (struct growth %+" PRId64 " B/inst)\n",
                inline_bases, inline_count, total_subseqs, frac, struct_growth);
        fprintf(stderr,
                "avg NET delta/instance: %+.3f B (payload) | %+.3f B (with "
                "malloc overhead)   [negative = saving]\n",
                net_payload, net_real);
    };
    sso_report("SSO 8B-base, struct stays 24B (inline 32 bases)", 8, 0);
    sso_report("SSO 64B struct, 48B inline union (inline 192 bases)", 48, 40);
    sso_report("SSO 128B struct, 112B inline union (inline 448 bases)", 112,
               104);

    fprintf(stderr, "Processed %d reads, output written to %s\n", read_count,
            output_path.c_str());
    return 0;
}

DumpCoresProgram::DumpCoresProgram(CLI::App &parent) {
    auto &app =
        *parent.add_subcommand("dump-cores", "Read FASTA and dump LCP cores");
    app.add_option("input", input, "Input FASTA file")->required();
    app.add_option("-o,--output", output, "Output TSV file")->required();
    app.add_option("--syncmer-k", syncmer_k, "K-mer length")
        ->capture_default_str();
    app.add_option("--submer-s", submer_s,
                   "Size of submers for mincode syncmers")
        ->capture_default_str();
    app.add_option("--submer-t", submer_t,
                   "Position of submers for mincode syncmers")
        ->capture_default_str();
    app.add_option("-l,--level", lcp_level, "LCP level to compute")
        ->capture_default_str();
    app.add_option("--seed", seed, "Random seed")->capture_default_str();
    app.callback([this]() {
        return run_dump_cores(input, output, syncmer_k, submer_s, submer_t,
                              lcp_level, seed);
    });
}
