#include <CLI/CLI.hpp>
#include <string>
#include <vector>

#include "CLI/CLI.hpp"

#define LCPTOOLS_IMPL
#include "asm.h"
#include "bamdoctor.h"
#include "ext/lcptools/lcptools_ho.h"
#ifdef BUILD_HISTOGRAM_UTILS
#include "stats.h"
#include "sthist.h"
#endif

using std::string;

// auto lam = []() -> int {
//     CPUType type = nibble2base_resolve();
//     std::cerr << "CPU Type: " << (int)type << std::endl;
//     return 0;
// }();
struct bloom_parameters;
int asm_program(CLI::App &parent) {
    auto *app_ = parent.add_subcommand("hifi", "run hifi assembly");
    auto &app = *app_;
    asm_config *p = new asm_config;
    auto *inputs = new std::vector<std::string>;
    auto *file_reading_mode = new std::string{"sequential"};
    app.add_option("input", *inputs, "Input bam files")->required();
    app.add_option("-o,--output", p->output, "Output gfa file")->required();
    app.add_option("-l,--level", p->lcp_level, "LCP level to compute")
        ->required();
    app.add_option("--dbg-lcmer-l", p->dbg_lcmer_l, "DBG LCMER window length")
        ->capture_default_str();
    app.add_option("-t,--threads", p->threads.work, "Number of worker threads")
        ->capture_default_str();
    app.add_option("-@,--io-threads", p->threads.io,
                   "Number of IO threads (Used by htslib for decompressing)")
        ->capture_default_str();
    app.add_option("--random-seed", p->seed, "Random seed")
        ->capture_default_str();
    app.add_option("--file-reading-mode", *file_reading_mode,
                   "Input reader mode: parallel or sequential")
        ->check(CLI::IsMember({"parallel", "sequential"}))
        ->capture_default_str();

    app.add_flag("--rle", p->rle, "Enable Run Length Encoding")
        ->capture_default_str();

    app.add_flag("--print-raw-dbg", p->print_precleanup_graph,
                 "Print DBG before cleanup (with .raw suffix)")
        ->capture_default_str();
    app.add_flag("--silent", p->silent, "Do not print progress to stderr")
        ->capture_default_str();
    app.add_option("--spaced-seed-mask", p->syncmer.spaced_seed_mask,
                   "Mask for spaced seeding. (0->skip, 1->keep, 2->mask left "
                   "bit, 3->mask right bit.")
        ->capture_default_str();
    app.add_option("--syncmer-k", p->syncmer.k, "K-mer length")
        ->capture_default_str();
    app.add_option("--density", p->syncmer.density,
                   "Density of selected k-mers")
        ->expected(0.0F, 1.0F);
    app.add_option("--submer-s", p->syncmer.mincode.s,
                   "Size of submers for mincode syncmers")
        ->excludes("--density");
    app.add_option("--submer-t", p->syncmer.mincode.t,
                   "Position of submers for mincode syncmers (Closed syncmer "
                   "if omitted)")
        ->excludes("--density")
        ->needs("--submer-s");
    auto bfg =
        app.add_option_group("bloomfilter", "Parameters to setup bloom filter");
    bfg->add_option("--projected-element-count", p->bloom.projected_count,
                    "Projected element count")
        ->capture_default_str();
    bfg->add_option("--false-positive-probability",
                    p->bloom.false_positive_rate, "False positive probability")
        ->capture_default_str();
    auto gcg = app.add_option_group("graphcleaning",
                                    "Parameters to setup graph cleaning");
    gcg->add_option("--iterations", p->graph_cleaning.iterations,
                    "Number of iterations of graph cleaning")
        ->capture_default_str();
    gcg->add_option("--max-tip-length", p->graph_cleaning.max_tip_length,
                    "Maximum tip length to cut in terms of cores")
        ->capture_default_str();
    gcg->add_option("--min-bubble-depth-difference-ratio",
                    p->graph_cleaning.min_bubble_depth_difference_percentage,
                    "Minimum bubble length to pop in terms of cores")
        ->capture_default_str();
    gcg->add_option("--min-bubble-length-difference",
                    p->graph_cleaning.min_bubble_length_difference,
                    "Minimum difference between the two ends of the bubbles to "
                    "pop. if depth(P1) * arg > depth(P2) then pop();")
        ->capture_default_str();
    gcg->add_flag(
           "--open-simple-loops", p->graph_cleaning.enable_loop_opening,
           "Enable simple loop opening during graph cleaning. Warning: may "
           "cause misassembly if loop copy counts are mispredicted.")
        ->capture_default_str();
    app.callback([=]() {
        if (inputs->empty()) {
            delete file_reading_mode;
            delete inputs;
            delete p;
            return -1;
        }
        p->parallel_file_reading = (*file_reading_mode == "parallel");
        p->inputs = *inputs;
        p->input = inputs->front();
        int ret = run_assembly(*p);
        delete file_reading_mode;
        delete inputs;
        delete p;
        return ret;
    });
    return 0;
}
int map_prog(CLI::App &parent) {
    auto *app_ = parent.add_subcommand("map", "run hifi assembly");
    auto &app = *app_;
    asm_config *p = new asm_config;
    auto *inputs = new std::vector<std::string>;
    auto *file_reading_mode = new std::string{"sequential"};
    app.add_option("input", *inputs, "Input bam files")->required();
    app.add_option("-o,--output", p->output,
                   "Second input file (query; kept for current run_mapping semantics)")
        ->required();
    app.add_option("-l,--level", p->lcp_level, "LCP level to compute")
        ->required();
    app.add_option("--dbg-lcmer-l", p->dbg_lcmer_l, "DBG LCMER window length")
        ->capture_default_str();
    app.add_option("-t,--threads", p->threads.work, "Number of worker threads")
        ->capture_default_str();
    app.add_option("-@,--io-threads", p->threads.io, "Number of IO threads")
        ->capture_default_str();
    app.add_option("--random-seed", p->seed, "Random seed")
        ->capture_default_str();
    app.add_option("--file-reading-mode", *file_reading_mode,
                   "Input reader mode: parallel or sequential")
        ->check(CLI::IsMember({"parallel", "sequential"}))
        ->capture_default_str();
    app.add_flag("--rle", p->rle, "Enable Run Length Encoding")
        ->capture_default_str();
    app.add_option("--spaced-seed-mask", p->syncmer.spaced_seed_mask,
                   "Mask for spaced seeding. (0->skip, 1->keep, 2->mask left "
                   "bit, 3->mask right bit.")
        ->capture_default_str();
    app.add_option("--syncmer-k", p->syncmer.k, "K-mer length")
        ->capture_default_str();
    app.add_option("--density", p->syncmer.density,
                   "Density of selected k-mers")
        ->expected(0.0F, 1.0F);
    app.add_option("--submer-s", p->syncmer.mincode.s,
                   "Size of submers for mincode syncmers")
        ->excludes("--density");
    app.add_option("--submer-t", p->syncmer.mincode.t,
                   "Position of submers for mincode syncmers (Closed syncmer "
                   "if omitted)")
        ->excludes("--density")
        ->needs("--submer-s");
    app.callback([=]() {
        if (inputs->empty()) {
            delete file_reading_mode;
            delete inputs;
            delete p;
            return -1;
        }
        p->parallel_file_reading = (*file_reading_mode == "parallel");
        p->inputs = *inputs;
        p->input = inputs->front();
        int ret = run_mapping(*p);
        delete file_reading_mode;
        delete inputs;
        delete p;
        return ret;
    });
    return 0;
}
#ifdef BUILD_HISTOGRAM_UTILS
void stats_program(CLI::App &parent) {
    auto *app_ = parent.add_subcommand("bam-stats", "Print stats on input bam");
    auto &app = *app_;
    struct parameters {
        string input_file;
        int update_frequency;
        int lcmer_level;
        double density;
    };
    parameters *params = new parameters{"", 100, 5, 0.01};
    app.add_option("input", params->input_file, "Input bam file")->required();
    app.add_option("--update-frequency", params->update_frequency,
                   "Update frequency")
        ->capture_default_str();
    app.add_option("-l,--level", params->lcmer_level,
                   "Depth of lcmers to be computed")
        ->capture_default_str();
    app.add_option("--density", params->density, "Density of selected k-mers")
        ->required();

    app.callback([=]() {
        int ret = run_stats(params->input_file, params->update_frequency,
                            params->lcmer_level, params->density);
        delete params;
        return ret;
    });
}

void stdin_hist_prog(CLI::App &parent) {
    auto *app_ = parent.add_subcommand("mephisto",
                                       "Testing HTSLib Multi File movements");
    auto &app = *app_;
    struct parameters {
        int dims;
        int freq;
    };
    // Allocate parameters as heap and set default values
    parameters *params = new parameters{1, 10};
    app.add_option("--dims", params->dims, "screen dimensions")
        ->capture_default_str();
    app.add_option("--freq", params->freq, "Frequency of histogram updates")
        ->capture_default_str();

    app.callback([=]() {
        int ret = run_hist_on_stdin(params->dims, params->freq);
        delete params;
        return ret;
    });
}
#endif

int sam_read_testing(CLI::App &parent) {
    struct parameters {
        string input;
        int lcp_level;
        int threads;
        SamReadTestingMode mode;
        string job;
        double density;
    };
    auto *app_ = parent.add_subcommand("bamdoctor", "read bam file fast");
    auto &app = *app_;
    parameters *p =
        new parameters{"", 5, 8, SamReadTestingMode::BASIC, "read_count", 0.01};
    app.add_option("input", p->input, "Input bam file")->required();
    app.add_option("-l,--level", p->lcp_level, "LCP level to compute")
        ->required();
    app.add_option("-t,--threads", p->threads, "Number of threads")
        ->capture_default_str();
    app.add_option("--density", p->density, "Density of selected k-mers")
        ->required();

    std::map<std::string, SamReadTestingMode> SamReadTestingMode_map = {
        {"basic", SamReadTestingMode::BASIC},
        {"pool", SamReadTestingMode::POOL},
    };

    std::map<std::string, std::string> job_map = {
        {"read_count", "read_count"},
        {"core_count", "core_count"},
    };
    app.add_option("--job", p->job, "Job to run")
        ->required()
        ->transform(CLI::CheckedTransformer(job_map, CLI::ignore_case));
    app.add_option("--mode", p->mode, "Mode to run")
        ->required()
        ->transform(
            CLI::CheckedTransformer(SamReadTestingMode_map, CLI::ignore_case));
    app.callback([=]() {
        int ret = run_sam_read_testing(p->input, p->threads, p->lcp_level,
                                       p->mode, p->job, p->density);
        delete p;
        return ret;
    });
    return 0;
}

#define DOCTEST_CONFIG_IMPLEMENT
#include "doctest.h"

int tests(int argc, char **argv) {
    doctest::Context context;
    context.setOption("abort-after", 5);
    context.setOption("order-by", "name");
    context.applyCommandLine(argc, argv);

    // overrides
    context.setOption("no-breaks", true);

    int res = context.run();

    if (context.shouldExit()) return res;
    return res;
}

int main(int argc, char **argv) {
    if (argc > 1 && string{argv[1]} == "tests") {
        return tests(argc - 1, argv + 1);
    }
    CLI::App app("lasm");
    app.require_subcommand(1);
    // app.set_help_all_flag("--help-all", "Expand all help");

    asm_program(app);
    map_prog(app);
#ifdef BUILD_HISTOGRAM_UTILS
    stats_program(app);
    stdin_hist_prog(app);
#endif
    sam_read_testing(app);
    CLI11_PARSE(app, argc, argv);
    return 0;
}
