#pragma once

#ifndef _ASM_H_
#define _ASM_H_

#include <inttypes.h>
#include <stdint.h>

#include <cassert>
#include <expected>
#include <optional>
#include <string>
#include <vector>
struct syncmer_params {
    enum class typeenum {
        Density,
        MinCodeOpen,
        MinCodeClosed,
        LcpSubsample,
        Invalid
    };
    constexpr typeenum type() const {
        // LCP-subsample is an alternative to syncmer selection: compute regular
        // LCP cores to a low level, subsample them by density, then deepen.
        // `density` is shared: with --lcp-low-level it is the LCP-core subsample
        // rate; otherwise it selects density syncmers (checked lcp-low first).
        if (lcp_low_level.has_value()) {
            return typeenum::LcpSubsample;
        }
        if (density.has_value()) {
            return typeenum::Density;
        }
        if (mincode.s.has_value()) {
            if (mincode.t.has_value()) {
                return typeenum::MinCodeOpen;
            } else {
                return typeenum::MinCodeClosed;
            }
        }
        return typeenum::Invalid;
    }
    unsigned k = 8;
    std::optional<double> density;
    struct {
        std::optional<unsigned long> s;
        std::optional<unsigned long> t;
    } mincode;
    // LCP-subsample mode parameter. When set, `density` (above) is reused as the
    // fraction of low-level LCP cores kept, rather than the syncmer density.
    std::optional<int> lcp_low_level;  // LCP level before subsampling
    std::string spaced_seed_mask = "11111111111111111111111111111111";

    // Validating mask builder: returns the packed mask or an error message.
    // This is the single source of truth for mask parsing; the total
    // compute_spaced_mask() below delegates here.
    std::expected<uint64_t, std::string> compute_spaced_mask_checked() const {
        if (spaced_seed_mask.size() > 32) {
            return std::unexpected(
                "mask string is too long " +
                std::to_string(spaced_seed_mask.size()) + " > 32");
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
                    return std::unexpected(
                        std::string("invalid character '") + c +
                        "' in mask string");
            }
        }
        return acc;
    }

    // Total accessor for the hot construction paths (member-init lists). The
    // config is validated once up front via validate(); reaching the error
    // branch here means that contract was broken.
    uint64_t compute_spaced_mask() const {
        auto r = compute_spaced_mask_checked();
        if (!r) {
            assert(false && "mask must be validated before use");
            __builtin_unreachable();
        }
        if ((int64_t)*r != -1) {
            fprintf(stderr,
                    "Using spaced seeds with mask %s built into 0x%" PRIX64
                    "!\n",
                    spaced_seed_mask.c_str(), *r);
        }
        return *r;
    }

    // Validate the whole syncmer configuration once, before any syncmer
    // functor is constructed. Lets callers surface a clean error instead of
    // exit()-ing from deep inside a constructor.
    std::expected<void, std::string> validate() const {
        if (type() == typeenum::Invalid) {
            return std::unexpected(std::string(
                "Provide a core selection: --density, --submer-s "
                "(optionally with --submer-t), or --lcp-low-level (with "
                "--density as the subsample rate)"));
        }
        if (type() == typeenum::LcpSubsample) {
            if (!density.has_value()) {
                return std::unexpected(std::string(
                    "--lcp-low-level requires --density (used as the LCP-core "
                    "subsample rate)"));
            }
            const double d = *density;
            if (d <= 0.0 || d > 1.0) {
                return std::unexpected(
                    std::string("--density must be in (0, 1]"));
            }
            if (*lcp_low_level < 1) {
                return std::unexpected(
                    std::string("--lcp-low-level must be >= 1"));
            }
            return {};  // LCP mode does not use the spaced-seed mask
        }
        if (auto r = compute_spaced_mask_checked(); !r) {
            return std::unexpected(r.error());
        }
        return {};
    }
    void print_cfg(FILE *f) const {
        fprintf(f, "spaced-seed-mask: \"%s\"\n", spaced_seed_mask.c_str());
        fprintf(f, "k: %d\n", k);
        if (lcp_low_level.has_value()) {
            fprintf(f, "core-selection: \"lcp-subsample\"\n");
            fprintf(f, "lcp-low-level: %d\n", *lcp_low_level);
            fprintf(f, "subsample-density: %f\n", density.value_or(0.0));
        } else if (density.has_value()) {
            fprintf(f, "kmer-selection: \"density\"\n");
            fprintf(f, "Density: %f\n", *density);
        } else if (mincode.s.has_value()) {
            if (mincode.t.has_value()) {
                fprintf(f, "kmer-selection: \"open-syncmers\"\n");
                fprintf(f, "syncmer-s:%lu\nsyncmer-t:%lu\n", *mincode.s,
                        *mincode.t);
            } else {
                fprintf(f, "kmer-selection: \"closed-syncmers\"\n");
                fprintf(f, "syncmer-s: %lu\n", *mincode.s);
            }
        } else {
            fprintf(f, "kmer-selection: \"none\"\n");
        }
    }
};
struct asm_config {
    std::vector<std::string> inputs;
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
    bool parallel_file_reading = true;
    syncmer_params syncmer;
    bool print_precleanup_graph = false;
    bool print_prelayout_graph = false;
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
            fprintf(f, "Inputs (0):\n");
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
#endif
