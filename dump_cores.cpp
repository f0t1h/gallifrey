#include "dump_cores.h"

#include <inttypes.h>
#include <stdint.h>

#include <CLI/CLI.hpp>
#include <cstdio>
#include <string>
#include <vector>

#include "asm.h"
#include "ext/lcptools/lcptools_ho.h"
#include "fastx.h"
#include "lcmer.h"
#include "syncmer.h"

static input_format dump_input_format(const std::string &path) {
    const input_format f = detect_input_format(path);
    if (f == input_format::fastx_plain || f == input_format::fastx_compressed) {
        return f;
    }
    return input_format::invalid;
}

static int run_dump_cores(const std::string &input_path,
                          const std::string &output_path, int syncmer_k,
                          int submer_s, int submer_t, int lcp_level,
                          uint64_t seed) {
    syncmer_params spara;
    spara.k = syncmer_k;
    spara.smer.s = submer_s;
    spara.smer.t = submer_t;

    using SYNCMER_HASHER = syncmer::MurMur;
    using SyncmerType = syncmer::GenericSubmerRule<uint32_t, SYNCMER_HASHER>;
    SyncmerType syncmer_functor(spara, static_cast<unsigned int>(seed));

    FILE *out = fopen(output_path.c_str(), "w");
    if (!out) {
        fprintf(stderr, "Cannot open output: %s\n", output_path.c_str());
        return 1;
    }
    fprintf(out, "read_name\tread_len\tcore_idx\tcore_label\tstart\tend\tspan\n");

    const input_format fmt = dump_input_format(input_path);
    if (fmt == input_format::invalid) {
        fprintf(stderr, "Not a FASTA/FASTQ input: %s\n", input_path.c_str());
        fclose(out);
        return 1;
    }
    fastx_reader reader;
    if (!reader.open(input_path, fmt == input_format::fastx_compressed, 1)) {
        fprintf(stderr, "Cannot open input: %s\n", input_path.c_str());
        fclose(out);
        return 1;
    }

    seq_record rec;
    lps str = {0, 0, 0, nullptr};
    charsi chs;
    pos_scratch pos_buf;
    int read_count = 0;

    while (reader.next(rec)) {
        lcp_read lr;
        parse_cores_from_nibbles(rec.nibbles, rec.len, lcp_level, &str,
                                 syncmer_functor, chs, lr, pos_buf);

        for (size_t i = 0; i < lr.cores.size(); ++i) {
            fprintf(out, "%s\t%u\t%zu\t%" PRIu32 "\t%u\t%u\t%u\n",
                    reader.name().c_str(), rec.len, i, lr.cores[i],
                    lr.pos[i].first, lr.pos[i].second,
                    lr.pos[i].second - lr.pos[i].first);
        }
        ++read_count;
    }
    free(rec.nibbles);
    if (reader.failed()) {
        fprintf(stderr, "Malformed FASTA/FASTQ input: %s\n",
                input_path.c_str());
        free_lps(&str);
        fclose(out);
        return 1;
    }

    free_lps(&str);
    fclose(out);
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
                   "s-mer length s for syncmer selection")
        ->capture_default_str();
    app.add_option("--submer-t", submer_t,
                   "Offset t of the smallest s-mer (open syncmer)")
        ->capture_default_str();
    app.add_option("-l,--level", lcp_level, "LCP level to compute")
        ->capture_default_str();
    app.add_option("--seed", seed, "Random seed")->capture_default_str();
    app.callback([this]() {
        // CLI11 callbacks return void; surface a nonzero rc as a
        // CLI::RuntimeError so CLI11_PARSE prints it and exits with that code.
        if (int rc = run_dump_cores(input, output, syncmer_k, submer_s,
                                    submer_t, lcp_level, seed);
            rc != 0) {
            throw CLI::RuntimeError("dump-cores failed", rc);
        }
    });
}
