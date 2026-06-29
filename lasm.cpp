#include <CLI/CLI.hpp>
#include <map>
#include <string>
// #include <vector>

#include "CLI/CLI.hpp"

#define LCPTOOLS_IMPL
#include "asm.h"
#include "bamdoctor.h"
#include "bench_layout.h"
#include "dump_cores.h"
#include "ext/lcptools/lcptools_ho.h"
using std::string;

// Each subcommand is a command object that owns its parsed parameters and
// registers itself with the parent CLI::App on construction. The objects live
// for the duration of main (which encloses CLI11_PARSE), so the options bound
// to their members and the callbacks capturing `this` stay valid through
// parsing. Copy/move is deleted so a registered `this` can never dangle.

class AsmProgram {
    asm_config cfg;
    std::string file_reading_mode{"sequential"};

public:
    explicit AsmProgram(CLI::App &parent) {
        auto &app = *parent.add_subcommand("hifi", "run hifi assembly");
        app.add_option("input", cfg.inputs, "Input bam files")->required();
        app.add_option("-o,--output", cfg.output, "Output gfa file")
            ->required();
        app.add_option("-l,--level", cfg.lcp_level, "LCP level to compute")
            ->required();
        app.add_option("--dbg-lcmer-l", cfg.dbg_lcmer_l,
                       "DBG LCMER window length")
            ->capture_default_str();
        app.add_option("-t,--threads", cfg.threads.work,
                       "Number of worker threads")
            ->capture_default_str();
        app.add_option("-@,--io-threads", cfg.threads.io,
                       "Number of IO threads (Used by htslib for decompressing)")
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
        app.add_option("--spaced-seed-mask", cfg.syncmer.spaced_seed_mask,
                       "Mask for spaced seeding. (0->skip, 1->keep, 2->mask "
                       "left bit, 3->mask right bit.")
            ->capture_default_str();
        app.add_option("--syncmer-k", cfg.syncmer.k, "K-mer length")
            ->capture_default_str();
        app.add_option("--density", cfg.syncmer.density,
                       "Density of selected k-mers")
            ->check(CLI::Range(0.0, 1.0));
        app.add_option("--submer-s", cfg.syncmer.mincode.s,
                       "Size of submers for mincode syncmers")
            ->excludes("--density");
        app.add_option("--submer-t", cfg.syncmer.mincode.t,
                       "Position of submers for mincode syncmers (Closed "
                       "syncmer if omitted)")
            ->excludes("--density")
            ->needs("--submer-s");
        auto bfg = app.add_option_group("bloomfilter",
                                        "Parameters to setup bloom filter");
        bfg->add_option("--projected-element-count", cfg.bloom.projected_count,
                        "Projected element count")
            ->capture_default_str();
        bfg->add_option("--false-positive-probability",
                        cfg.bloom.false_positive_rate,
                        "False positive probability")
            ->capture_default_str();
        auto gcg = app.add_option_group("graphcleaning",
                                        "Parameters to setup graph cleaning");
        gcg->add_option("--iterations", cfg.graph_cleaning.iterations,
                        "Number of iterations of graph cleaning")
            ->capture_default_str();
        gcg->add_option("--max-tip-length", cfg.graph_cleaning.max_tip_length,
                        "Maximum tip length to cut in terms of cores")
            ->capture_default_str();
        gcg->add_option(
               "--min-bubble-depth-difference-ratio",
               cfg.graph_cleaning.min_bubble_depth_difference_percentage,
               "Minimum bubble length to pop in terms of cores")
            ->capture_default_str();
        gcg->add_option(
               "--min-bubble-length-difference",
               cfg.graph_cleaning.min_bubble_length_difference,
               "Minimum difference between the two ends of the bubbles to "
               "pop. if depth(P1) * arg > depth(P2) then pop();")
            ->capture_default_str();
        gcg->add_flag(
               "--open-simple-loops", cfg.graph_cleaning.enable_loop_opening,
               "Enable simple loop opening during graph cleaning. Warning: may "
               "cause misassembly if loop copy counts are mispredicted.")
            ->capture_default_str();
        app.callback([this]() {
            // Validate the syncmer configuration (selection mode + mask) up
            // front. Throwing a CLI error makes CLI11 print the message and
            // exit with a non-zero status.
            if (auto v = cfg.syncmer.validate(); !v) {
                throw CLI::ValidationError(v.error());
            }
            cfg.parallel_file_reading = (file_reading_mode == "parallel");
            // Propagate assembly failures to the process exit status: CLI11
            // callbacks return void, so surface a nonzero result as a
            // CLI::RuntimeError, which CLI11_PARSE prints and exits with.
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
            ->capture_default_str();
        app.add_option("--density", density, "Density of selected k-mers")
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
            return run_sam_read_testing(input, threads, lcp_level, mode, job,
                                        density);
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
    BenchLayoutProgram benchlayout{app};

    CLI11_PARSE(app, argc, argv);
    return 0;
}
