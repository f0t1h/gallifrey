#include <CLI/CLI.hpp>
#include <inttypes.h>
#include <stdint.h>

#include <cstdio>
#include <string>
#include <vector>

#include "asm.h"
#include "dump_cores.h"
#include "ext/lcptools/lcptools_ho.h"
#include "htslib/sam.h"
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
    fprintf(out, "read_name\tread_len\tcore_idx\tcore_label\tstart\tend\tspan\n");

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

        for (size_t i = 0; i < lr.cores.size(); ++i) {
            fprintf(out, "%s\t%zu\t%zu\t%" PRIu32 "\t%u\t%u\t%u\n",
                    reader.name.c_str(), len, i, lr.cores[i], lr.pos[i].first,
                    lr.pos[i].second, lr.pos[i].second - lr.pos[i].first);
        }
        ++read_count;
    }

    free_lps(&str);
    fclose(out);
    fprintf(stderr, "Processed %d reads, output written to %s\n", read_count,
            output_path.c_str());
    return 0;
}

int dump_cores_program(CLI::App &parent) {
    struct parameters {
        std::string input;
        std::string output;
        int syncmer_k;
        int submer_s;
        int submer_t;
        int lcp_level;
        uint64_t seed;
    };
    auto *app_ =
        parent.add_subcommand("dump-cores", "Read FASTA and dump LCP cores");
    auto &app = *app_;
    parameters *p = new parameters{"", "", 11, 4, 5, 4, 0};
    app.add_option("input", p->input, "Input FASTA file")->required();
    app.add_option("-o,--output", p->output, "Output TSV file")->required();
    app.add_option("--syncmer-k", p->syncmer_k, "K-mer length")
        ->capture_default_str();
    app.add_option("--submer-s", p->submer_s,
                   "Size of submers for mincode syncmers")
        ->capture_default_str();
    app.add_option("--submer-t", p->submer_t,
                   "Position of submers for mincode syncmers")
        ->capture_default_str();
    app.add_option("-l,--level", p->lcp_level, "LCP level to compute")
        ->capture_default_str();
    app.add_option("--seed", p->seed, "Random seed")->capture_default_str();
    app.callback([=]() {
        int ret = run_dump_cores(p->input, p->output, p->syncmer_k, p->submer_s,
                                 p->submer_t, p->lcp_level, p->seed);
        delete p;
        return ret;
    });
    return 0;
}
