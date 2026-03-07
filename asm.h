#pragma once

#ifndef _ASM_H_
#define _ASM_H_
#include <inttypes.h>
#include <stdint.h>

#include <string>
#include <vector>
struct syncmer_params {
    static constexpr unsigned long invalid_param_value = 127127127;
    enum class typeenum { Density, MinCodeOpen, MinCodeClosed, Invalid };
    constexpr typeenum type() const {
        if (invalid_density != invalid_param_value) {
            return typeenum::Density;
        }
        if (mincode.s != invalid_param_value) {
            if (mincode.t != invalid_param_value) {
                return typeenum::MinCodeOpen;
            } else {
                return typeenum::MinCodeClosed;
            }
        }
        return typeenum::Invalid;
    }
    unsigned k = 8;
    union {
        double density;
        unsigned long invalid_density = invalid_param_value;
    };
    struct {
        unsigned long s = invalid_param_value;
        unsigned long t = invalid_param_value;
    } mincode;
    std::string spaced_seed_mask = "11111111111111111111111111111111";
    uint64_t compute_spaced_mask() const {
        if (spaced_seed_mask.size() > 32) {
            fprintf(stderr, "Error: mask string is too long %zu > 32\n",
                    spaced_seed_mask.size());
            exit(EXIT_FAILURE);
        }
        uint64_t acc = -1;
        for (char c : spaced_seed_mask) {
            acc <<= 2;
            switch (c) {
                case '0':
                    break;
                case '1':
                    acc |= 3;
                    break;
                case '2':
                    acc |= 2;
                    break;
                case '3':
                    acc |= 1;
                    break;
                default:
                    fprintf(stderr,
                            "Error: invalid character '%c' in mask string\n",
                            c);
                    exit(EXIT_FAILURE);
            }
        }
        if ((int64_t)acc != -1) {
            fprintf(stderr,
                    "Using spaced seeds with mask %s built into 0x%" PRIX64
                    "!\n",
                    spaced_seed_mask.c_str(), acc);
        }
        return acc;
    }
    void print_cfg(FILE *f) const {
        fprintf(f, "spaced-seed-mask: \"%s\"\n", spaced_seed_mask.c_str());
        fprintf(f, "k: %d\n", k);
        if (invalid_density != invalid_param_value) {
            fprintf(f, "kmer-selection: \"density\"\n");
            fprintf(f, "Density: %f\n", density);
        } else {
            if (mincode.t != invalid_param_value) {
                fprintf(f, "kmer-selection: \"open-syncmers\"\n");
                fprintf(f, "syncmer-s:%lu\nsyncmer-t:%lu\n", mincode.s, mincode.t);
            } else {
                fprintf(f, "kmer-selection: \"closed-syncmers\"\n");
                fprintf(f, "syncmer-s: %lu\n", mincode.s);
            }
        }
    }
};
struct asm_config {
    std::vector<std::string> inputs;
    std::string input;
    std::string output;
    int lcp_level = 3;
    int dbg_lcmer_l = 8;
    int olc_lcmer_l = 8;
    struct {
        int io = 4;
        int work = 8;
    } threads;
    size_t seed = 1453;
    bool rle = true;
    bool parallel_file_reading = false;
    syncmer_params syncmer;
    bool print_precleanup_graph = false;
    bool silent = false;
    struct {
        uint64_t projected_count = 100'000'000;
        double false_positive_rate = 0.001;
    } bloom;
    struct {
        int iterations = 2;
        int max_tip_length = 8 * 2;
        int min_bubble_length_difference = 2;
        double min_bubble_depth_difference_percentage = 2;
        bool enable_loop_opening = false;
    } graph_cleaning;

    void print_cfg(FILE *f) const {
        if (!inputs.empty()) {
            fprintf(f, "Inputs (%zu):\n", inputs.size());
            for (const auto &in : inputs) {
                fprintf(f, "\t%s\n", in.c_str());
            }
        } else {
            fprintf(f, "Input: %s\n", input.c_str());
        }
        fprintf(f, "Output: %s\n", output.c_str());
        fprintf(f, "LCP Level: %d\n", lcp_level);
        fprintf(f, "DBG LCMER Length: %d\n", dbg_lcmer_l);
        fprintf(f, "OLC LCMER Length: %d\n", olc_lcmer_l);
        fprintf(f, "Threads:\n\tIO: %d\n\tWork: %d\n", threads.io, threads.work);
        fprintf(f, "Seed: %zu\n", seed);
        fprintf(f, "RLE: %s\n", rle ? "true" : "false");
        fprintf(f, "File Reading Mode: %s\n",
                parallel_file_reading ? "parallel" : "sequential");
        syncmer.print_cfg(f);
        fprintf(f,
                "Bloom Filter Parameters:\n\tProjected Count: %" PRIu64
                "\n\tFalse Positive "
                "Rate: %.3f\n",
                bloom.projected_count, bloom.false_positive_rate);
        fprintf(f,
                "Graph Cleaning Parameters: \n\tIterations: %d\n\tMax Tip "
                "Length: %d\n\tMin Bubble Length Difference: %d\n\tMin Bubble "
                "Depth Difference Rate: %.2f\n\tOpen Simple Loops: %s\n",
                graph_cleaning.iterations, graph_cleaning.max_tip_length,
                graph_cleaning.min_bubble_length_difference,
                graph_cleaning.min_bubble_depth_difference_percentage,
                graph_cleaning.enable_loop_opening ? "true" : "false");
    }
};


int run_assembly(const asm_config &param);
int run_mapping(const asm_config &param);
#endif
