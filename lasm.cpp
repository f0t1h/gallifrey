#include <CLI/CLI.hpp>
#include <map>
#include <string>
// #include <vector>

#include "CLI/CLI.hpp"

#define LCPTOOLS_IMPL
#include "asm.h"
#include "bamdoctor.h"
#include "dump_cores.h"
#include "ext/lcptools/lcptools_ho.h"
using std::string;

// Each subcommand is a command object that owns its parsed parameters and
// registers itself with the parent CLI::App on construction. The objects
// live for all of main (enclosing CLI11_PARSE), so bound options and the
// `this`-capturing callbacks stay valid. Copy/move is deleted.

class AsmProgram {
    asm_config cfg;
    std::string file_reading_mode{"sequential"};

public:
    explicit AsmProgram(CLI::App &parent) {
        auto &app = *parent.add_subcommand("hifi", "run hifi assembly");
        app.add_option("input", cfg.inputs,
                       "Input read files: BAM, CRAM, SAM, or FASTA/FASTQ "
                       "(plain or gzipped). Format is detected from content; "
                       "formats may be mixed.")
            ->required();
        app.add_option("-o,--output", cfg.output, "Output gfa file")
            ->required();
        app.add_option("-l,--level", cfg.lcp_level, "LCP level to compute")
            ->required()
            ->check(CLI::Range(1, 64));
        app.add_option("--dbg-lcmer-l", cfg.dbg_lcmer_l,
                       "DBG LCMER window length")
            ->check(CLI::Range(1, 64))
            ->capture_default_str();
        app.add_option("--olc-lcmer-l", cfg.olc_lcmer_l,
                       "OLC (overlap) LCMER window length. Overlap detection "
                       "queries the DBG (built at --dbg-lcmer-l), so this must "
                       "equal --dbg-lcmer-l or overlaps are lost")
            ->check(CLI::Range(1, 64))
            ->capture_default_str();
        app.add_option("-t,--threads", cfg.threads.work,
                       "Number of worker threads")
            ->check(CLI::Range(1, MAXIMUM_WORK_THREADS))
            ->capture_default_str();
        app.add_option("-@,--io-threads", cfg.threads.io,
                       "Number of IO threads (Used by htslib for decompressing)")
            ->check(CLI::Range(0, 64))
            ->capture_default_str();
        app.add_option("--random-seed", cfg.seed, "Random seed")
            ->capture_default_str();
        app.add_option("--file-reading-mode", file_reading_mode,
                       "Input reader mode: parallel or sequential")
            ->check(CLI::IsMember({"parallel", "sequential"}))
            ->capture_default_str();

        app.add_flag("--rle", cfg.rle, "Enable Run Length Encoding")
            ->capture_default_str();

        app.add_flag("--print-raw-dbg", cfg.print_precleanup_graph,
                     "Print DBG before cleanup (with .raw.gfa suffix)")
            ->capture_default_str();
        app.add_flag("--print-bi-dbg", cfg.print_prelayout_graph,
                     "Print DBG before layout (with .bi.gfa suffix)")
            ->capture_default_str();
        app.add_flag("--silent", cfg.silent, "Do not print progress to stderr")
            ->capture_default_str();
        app.add_flag("--drop-crumbs,!--keep-crumbs", cfg.drop_crumbs,
                     "Drop crumb nodes (fewer than olc-lcmer-l+1 cores, so "
                     "invisible to overlap detection) that have a single "
                     "predecessor and continue into a longer node, instead of "
                     "emitting them as C segments. Halves the contig count at "
                     "a small genome-fraction cost.")
            ->capture_default_str();
        app.add_option("--spaced-seed-mask", cfg.syncmer.spaced_seed_mask,
                       "Mask for spaced seeding. 0->remove the position from "
                       "selection entirely (effective k = k - #zeros; syncmer "
                       "s-windows span only kept positions), 1->keep, "
                       "2->keep but mask left bit, 3->keep but mask right "
                       "bit. With --rle, positions index RLE runs.")
            ->capture_default_str();
        app.add_option("--syncmer-k", cfg.syncmer.k,
                       "K-mer length (table is 4^k bits; k=15 ⇒ 128 MiB)")
            ->check(CLI::Range(3u, 32u))
            ->capture_default_str();
        app.add_option("--density", cfg.syncmer.density,
                       "Mincode-submer selection (Edgar 2021): keep k-mers whose "
                       "hashed code < H*density (compression c = 1/density). "
                       "'density' is the resulting metric, not a syncmer.")
            ->check(CLI::Range(0.0, 1.0));
        app.add_option("--submer-s", cfg.syncmer.smer.s,
                       "s-mer length s for syncmer selection. Without --submer-t "
                       "this is a closed syncmer (smallest s-mer at the first or "
                       "last position).")
            ->check(CLI::Range(1ul, 31ul))
            ->excludes("--density");
        app.add_option("--submer-t", cfg.syncmer.smer.t,
                       "Offset t of the smallest s-mer: selects an open syncmer "
                       "(smallest s-mer at position t).")
            ->check(CLI::Range(0ul, 31ul))
            ->excludes("--density")
            ->needs("--submer-s");
        app.add_flag("--canonical", cfg.syncmer.canonical,
                     "Strand-canonical selection: apply the submer rule to "
                     "min(kmer, reverse-complement) so a k-mer and its reverse "
                     "complement are always co-selected.")
            ->capture_default_str();
        auto bfg = app.add_option_group("bloomfilter",
                                        "Parameters to setup bloom filter");
        bfg->add_option("--projected-element-count", cfg.bloom.projected_count,
                        "Projected element count (human-scale ~3e9 is fine; the "
                        "filter prints its allocated size on startup)")
            ->check(CLI::Range((int64_t)1, (int64_t)10000000000))
            ->capture_default_str();
        bfg->add_option("--false-positive-probability",
                        cfg.bloom.false_positive_rate,
                        "False positive probability")
            ->capture_default_str();
        auto gcg = app.add_option_group("graphcleaning",
                                        "Parameters to setup graph cleaning");
        gcg->add_option("--iterations", cfg.graph_cleaning.iterations,
                        "Number of iterations of graph cleaning")
            ->check(CLI::Range(0, 100))
            ->capture_default_str();
        gcg->add_option("--max-tip-length", cfg.graph_cleaning.max_tip_length,
                        "Maximum tip length to cut in terms of cores")
            ->check(CLI::Range(0, 1000000))
            ->capture_default_str();
        gcg->add_option(
               "--min-bubble-depth-difference-ratio",
               cfg.graph_cleaning.min_bubble_depth_difference_percentage,
               "Pops the shallower bubble arm when min_depth(shallow) / "
               "min_depth(deep) < arg; arg >= 1 pops regardless of depth.")
            ->capture_default_str();
        gcg->add_option(
               "--min-bubble-length-difference",
               cfg.graph_cleaning.min_bubble_length_difference,
               "Maximum length ratio between bubble arms to pop: pops when "
               "longer_arm <= shorter_arm * arg; more divergent arms are kept.")
            ->capture_default_str();
        gcg->add_option(
               "--max-bubble-arm-length",
               cfg.graph_cleaning.max_bubble_arm_length,
               "Maximum arm length (in cores) of a bubble arm to pop; "
               "longer arms are treated as real copy divergence and kept.")
            ->check(CLI::Range(0, 1000000))
            ->capture_default_str();
        app.callback([this]() {
            // Validate the syncmer configuration (selection mode + mask) up
            // front; a CLI error makes CLI11 print the message and exit nonzero.
            if (auto v = cfg.syncmer.validate(); !v) {
                throw CLI::ValidationError(v.error());
            }
            cfg.parallel_file_reading = (file_reading_mode == "parallel");
            // CLI11 callbacks return void; surface a nonzero rc as a
            // CLI::RuntimeError so the failure reaches the process exit status.
            if (int rc = run_assembly(cfg); rc != 0) {
                throw CLI::RuntimeError("assembly failed", rc);
            }
        });
    }
    AsmProgram(const AsmProgram &) = delete;
    AsmProgram &operator=(const AsmProgram &) = delete;
};

class SamReadTesting {
    string input;
    int lcp_level = 5;
    int threads = 8;
    SamReadTestingMode mode = SamReadTestingMode::BASIC;
    string job = "read_count";
    double density = 0.01;

public:
    explicit SamReadTesting(CLI::App &parent) {
        auto &app = *parent.add_subcommand("bamdoctor", "read bam file fast");
        app.add_option("input", input, "Input bam file")->required();
        app.add_option("-l,--level", lcp_level, "LCP level to compute")
            ->required();
        app.add_option("-t,--threads", threads, "Number of threads")
            ->check(CLI::Range(1, MAXIMUM_WORK_THREADS))
            ->capture_default_str();
        app.add_option("--density", density, "Density of selected k-mers")
            ->check(CLI::Range(0.0, 1.0))
            ->required();

        std::map<std::string, SamReadTestingMode> SamReadTestingMode_map = {
            {"basic", SamReadTestingMode::BASIC},
            {"pool", SamReadTestingMode::POOL},
        };
        std::map<std::string, std::string> job_map = {
            {"read_count", "read_count"},
            {"core_count", "core_count"},
        };
        app.add_option("--job", job, "Job to run")
            ->required()
            ->transform(CLI::CheckedTransformer(job_map, CLI::ignore_case));
        app.add_option("--mode", mode, "Mode to run")
            ->required()
            ->transform(CLI::CheckedTransformer(SamReadTestingMode_map,
                                                CLI::ignore_case));
        app.callback([this]() {
            // CLI11 callbacks return void; re-raise a nonzero rc as a CLI error
            // so it reaches the process exit status.
            if (int rc = run_sam_read_testing(input, threads, lcp_level, mode,
                                              job, density);
                rc != 0) {
                throw CLI::RuntimeError("bamdoctor failed", rc);
            }
        });
    }
    SamReadTesting(const SamReadTesting &) = delete;
    SamReadTesting &operator=(const SamReadTesting &) = delete;
};

int main(int argc, char **argv) {
    CLI::App app("lasm");
    app.require_subcommand(1);
    // app.set_help_all_flag("--help-all", "Expand all help");

    AsmProgram hifi{app};
    SamReadTesting bamdoctor{app};
    DumpCoresProgram dumpcores{app};

    CLI11_PARSE(app, argc, argv);
    return 0;
}
