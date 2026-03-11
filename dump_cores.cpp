// Standalone tool: read FASTA, compute LCP cores, dump to TSV.
// Build:  (from build dir)
//   cmake --build . --target dump_cores
// Usage:
//   ./dump_cores <input.fa> <output.tsv> [--syncmer-k 11] [--submer-s 4]
//        [--submer-t 5] [--lcp-level 4] [--seed 0]

#define LCPTOOLS_IMPL

#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>

#include "asm.h"
#include "ext/lcptools/lcptools_ho.h"
#include "htslib/sam.h"  // for seq_nt16_table, bam_seqi
#include "lcmer.h"
#include "syncmer.h"

// Simple FASTA reader — yields one record at a time.
struct FastaReader {
    FILE *fp = nullptr;
    std::string name;
    std::string seq;

    FastaReader(const char *path) { fp = fopen(path, "r"); }
    ~FastaReader() {
        if (fp) fclose(fp);
    }
    // Returns false when EOF.
    bool next() {
        if (!fp) return false;
        name.clear();
        seq.clear();
        int c;
        // Skip to next '>'
        while ((c = fgetc(fp)) != EOF) {
            if (c == '>') break;
        }
        if (c == EOF) return false;
        // Read name line
        while ((c = fgetc(fp)) != EOF && c != '\n' && c != '\r') {
            name.push_back(c);
        }
        // Read sequence lines until next '>' or EOF
        while ((c = fgetc(fp)) != EOF) {
            if (c == '>') {
                ungetc(c, fp);
                break;
            }
            if (c != '\n' && c != '\r' && c != ' ') {
                seq.push_back(c);
            }
        }
        return !seq.empty();
    }
};

// Convert ASCII DNA to BAM nibble-packed format.
static void ascii_to_bam_nibbles(const std::string &seq,
                                  std::vector<uint8_t> &out) {
    size_t len = seq.size();
    out.resize((len + 1) / 2);
    for (size_t i = 0; i < len; ++i) {
        uint8_t nibble = seq_nt16_table[(uint8_t)seq[i]];
        if (i & 1) {
            out[i / 2] |= nibble;
        } else {
            out[i / 2] = nibble << 4;
        }
    }
}

int main(int argc, char **argv) {
    // Defaults matching quickrunlasmquast.sh
    int syncmer_k = 11;
    int submer_s = 4;
    int submer_t = 5;
    int lcp_level = 4;
    uint64_t seed = 0;

    if (argc < 3) {
        fprintf(stderr,
                "Usage: %s <input.fa> <output.tsv> [--syncmer-k K] "
                "[--submer-s S] [--submer-t T] [--lcp-level L] [--seed N]\n",
                argv[0]);
        return 1;
    }
    const char *input_path = argv[1];
    const char *output_path = argv[2];

    for (int i = 3; i < argc; i += 2) {
        if (i + 1 >= argc) break;
        if (!strcmp(argv[i], "--syncmer-k"))
            syncmer_k = atoi(argv[i + 1]);
        else if (!strcmp(argv[i], "--submer-s"))
            submer_s = atoi(argv[i + 1]);
        else if (!strcmp(argv[i], "--submer-t"))
            submer_t = atoi(argv[i + 1]);
        else if (!strcmp(argv[i], "--lcp-level"))
            lcp_level = atoi(argv[i + 1]);
        else if (!strcmp(argv[i], "--seed"))
            seed = atol(argv[i + 1]);
    }

    // Build syncmer functor
    syncmer_params spara;
    spara.k = syncmer_k;
    spara.mincode.s = submer_s;
    spara.mincode.t = submer_t;

    using SYNCMER_HASHER = syncmer::MurMur;
    using SyncmerType = syncmer::GenericSyncmerRule<uint32_t, SYNCMER_HASHER>;
    SyncmerType syncmer_functor(spara, static_cast<unsigned int>(seed));

    FILE *out = fopen(output_path, "w");
    if (!out) {
        fprintf(stderr, "Cannot open output: %s\n", output_path);
        return 1;
    }
    fprintf(out, "read_name\tread_len\tcore_idx\tcore_label\tstart\tend\tspan\n");

    FastaReader reader(input_path);
    if (!reader.fp) {
        fprintf(stderr, "Cannot open input: %s\n", input_path);
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

        for (size_t i = 0; i < lr.cores.size(); ++i) {
            fprintf(out, "%s\t%zu\t%zu\t%llu\t%u\t%u\t%u\n", reader.name.c_str(),
                    len, i, (unsigned long long)lr.cores[i], lr.pos[i].first,
                    lr.pos[i].second, lr.pos[i].second - lr.pos[i].first);
        }
        read_count++;
    }

    free_lps(&str);
    fclose(out);
    fprintf(stderr, "Processed %d reads, output written to %s\n", read_count,
            output_path);
    return 0;
}
