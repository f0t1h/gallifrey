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
// Submer selection parameters (terminology per Edgar 2021, "Syncmers are more
// sensitive than minimizers…"). A "submer" is a k-mer selected by a submer
// rule. The rules here are: a mincode submer (density-based: select if the
// hashed code < H/c), and syncmers (select by the position of the smallest
// s-mer µ(κ,s) — closed = first or last s-mer, open = the s-mer at offset t).
struct syncmer_params {
    enum class typeenum { MincodeSubmer, OpenSyncmer, ClosedSyncmer, Invalid };
    constexpr typeenum type() const {
        if (density.has_value()) {
            return typeenum::MincodeSubmer;
        }
        if (smer.s.has_value()) {
            if (smer.t.has_value()) {
                return typeenum::OpenSyncmer;
            } else {
                return typeenum::ClosedSyncmer;
            }
        }
        return typeenum::Invalid;
    }
    unsigned k = 8;
    // Mincode-submer parameter: target density d (compression factor c = 1/d).
    std::optional<double> density;
    // Syncmer parameters: s-mer length s and (open syncmer) offset t.
    struct {
        std::optional<unsigned long> s;
        std::optional<unsigned long> t;
    } smer;
    // Apply the rule to the strand-canonical k-mer (min of the k-mer and its
    // reverse complement) so a k-mer and its reverse complement are always
    // co-selected. Opt-in; default off preserves the directional behavior.
    bool canonical = false;
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

    // Effective k after spaced-seed compaction: each '0' mask char removes
    // its position, so selection runs on the dense (k - #zeros)-mer. Only
    // positions within the k-mer window count (chars beyond k are shifted
    // above the window and never apply). '2'/'3' keep their position.
    std::expected<unsigned, std::string> effective_k_checked() const {
        auto r = compute_spaced_mask_checked();
        if (!r) {
            return std::unexpected(r.error());
        }
        unsigned kept = 0;
        for (unsigned i = 0; i < k; ++i) {
            kept += ((*r >> (2 * i)) & 3) != 0;
        }
        return kept;
    }

    // Validate the whole syncmer configuration once, before any syncmer
    // functor is constructed. Lets callers surface a clean error instead of
    // exit()-ing from deep inside a constructor.
    std::expected<void, std::string> validate() const {
        if (type() == typeenum::Invalid) {
            return std::unexpected(
                std::string("Provide a k-mer selection: either --density or "
                            "--submer-s (optionally with --submer-t)"));
        }
        auto ek = effective_k_checked();
        if (!ek) {
            return std::unexpected(ek.error());
        }
        const unsigned k_eff = *ek;
        if (k_eff == 0) {
            return std::unexpected(
                "--spaced-seed-mask removes every k-mer position (effective "
                "k = 0)");
        }
        // Cross-field numeric checks: reject configurations that select no
        // k-mers or index an s-window out of the (compacted) k-mer before any
        // functor is constructed (a per-option CLI::Range can't express these
        // relations). s/t bounds are against the EFFECTIVE k: masked positions
        // are non-existent in the selection computation.
        if (density.has_value() && *density == 0.0) {
            return std::unexpected(std::string(
                "--density must be > 0 (0 selects no k-mers -> empty "
                "assembly)"));
        }
        if (smer.s.has_value()) {
            if (*smer.s >= k_eff) {
                return std::unexpected(
                    "--submer-s (" + std::to_string(*smer.s) +
                    ") must be < effective k (" + std::to_string(k_eff) +
                    " = --syncmer-k " + std::to_string(k) + " minus " +
                    std::to_string(k - k_eff) + " masked positions)");
            }
            if (smer.t.has_value() && *smer.t >= k_eff - *smer.s + 1) {
                return std::unexpected(
                    "--submer-t (" + std::to_string(*smer.t) +
                    ") must be < effective-k - s + 1 (" +
                    std::to_string(k_eff - *smer.s + 1) + ")");
            }
        }
        return {};
    }
    void print_cfg(FILE *f) const {
        fprintf(f, "spaced-seed-mask: \"%s\"\n", spaced_seed_mask.c_str());
        fprintf(f, "k: %d\n", k);
        if (auto ek = effective_k_checked(); ek && *ek != k) {
            fprintf(f, "effective-k (after spaced-seed compaction): %u\n",
                    *ek);
        }
        if (density.has_value()) {
            fprintf(f, "submer-rule: \"mincode-submer\"\n");
            fprintf(f, "density: %f\n", *density);
        } else if (smer.s.has_value()) {
            if (smer.t.has_value()) {
                fprintf(f, "submer-rule: \"open-syncmer\"\n");
                fprintf(f, "syncmer-s: %lu\nsyncmer-t: %lu\n", *smer.s,
                        *smer.t);
            } else {
                fprintf(f, "submer-rule: \"closed-syncmer\"\n");
                fprintf(f, "syncmer-s: %lu\n", *smer.s);
            }
        } else {
            fprintf(f, "submer-rule: \"none\"\n");
        }
        fprintf(f, "canonical: %s\n", canonical ? "true" : "false");
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
    // When non-empty, DBG construction classifies each read's lcmer insertions
    // (novel/coverage/branch-open/branch-join) and writes a per-read TSV here.
    // Opt-in analysis only; empty leaves the construction path byte-identical.
    std::string dump_insert_types_path;

    // --thread-bridge: read-threading bridge pass. During construction it retains
    // (streaming) the walks of reads that cross a repeat, capping redundancy per
    // bridge via a counting bloom; after cleaning it votes read-supported unitig
    // adjacencies across tangles and confirms unambiguous, well-supported bridges.
    // Opt-in; off leaves the pipeline byte-identical.
    struct {
        bool enable = false;
        int cap = 8;          // max retained reads per bridge key (counting bloom)
        int min_support = 5;  // reads needed to confirm a bridge adjacency
        int min_anchor = 10;  // min unitig length (cores) to be a bridge anchor
        // When non-empty, also materialize a scaffold FASTA at this path: the
        // bridge edges are walked (orientation-aware) to concatenate the joined
        // supernode sequences (N-gap between); unbridged supernodes are emitted
        // as-is. Opt-in; the GFA output is unaffected.
        std::string fasta_path;
    } thread_bridge;
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
    struct {
        // M-contig bridging through C nodes (see
        // docs/contiguity-and-spaced-seed-plan.md, Task 1). 0 = off,
        // 1 = thread (topology: unique through-path only),
        // 2 = reads (topology + spanning-read evidence).
        int mode = 0;
        // Minimum EC (supporting DBG edges) on each bridge superedge.
        size_t min_edge_support = 1;
        // Max total superedges on the C node: exact 2 in mode 1; looser cap
        // in mode 2 where read evidence arbitrates the tangle.
        int max_c_degree = 2;
        int max_c_degree_reads = 6;
        // Reads that must span both junction anchors (mode 2).
        int min_read_support = 2;
        // Chaining iterations (M-C-M-C-... merges become possible after a
        // merge rewires the supernode graph).
        int max_passes = 4;
    } bridge;

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
