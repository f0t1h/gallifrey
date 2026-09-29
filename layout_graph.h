#ifndef _LAYOUT_GRAPH_H
#define _LAYOUT_GRAPH_H

#include <sys/stat.h>

#include <algorithm>
#include <atomic>
#include <chrono>
#include <cinttypes>
#include <cstdio>
#include <cstring>
#include <gtl/phmap_fwd_decl.hpp>
#include <gtl/vector.hpp>
#include <map>
#include <optional>
#include <queue>
// #include <stack>
#include <string>
#include <thread>
#include <tuple>
#include <utility>

#include "align.h"
#include "graph.h"

namespace layout_graph_algo {

using layout_result_t =
    gtl::vector<gtl::flat_hash_map<u64, std::tuple<i64, i64, bool>>>;

template <class RecigMapT>
inline layout_result_t connected_comp_plus_layout(const RecigMapT &recig) {
    layout_result_t layouts;
    std::queue<std::pair<u64, typename RecigMapT::mapped_type>> q;
    gtl::flat_hash_set<u64> visited;

    // Precompute node → length. A node's length is the query.len of the first
    // overlap keyed by that node (when > 0); otherwise the target_span.len of
    // the first overlap (iteration order) that references it as a target.
    gtl::flat_hash_map<u64, i64> first_qry_len;  // key node → its first query.len
    gtl::flat_hash_map<u64, i64> first_tgt_len;  // target node → first tgt len>0
    for (const auto &kv : recig) {
        const auto &rm = kv.second;
        first_qry_len.try_emplace(kv.first, static_cast<i64>(rm.query.len));
        if (rm.target_span.len > 0)
            first_tgt_len.try_emplace(rm.target,
                                      static_cast<i64>(rm.target_span.len));
    }

    auto infer_len_for_node = [&](u64 id) -> i64 {
        auto qit = first_qry_len.find(id);
        if (qit != first_qry_len.end() && qit->second > 0) return qit->second;
        auto tit = first_tgt_len.find(id);
        if (tit != first_tgt_len.end()) return tit->second;
        // Invariant I1: every participating node carries a core-span length
        // (CoreSpan.len ≥ 1). Warn if violated; treat as 0.
        std::fprintf(stderr,
                     "[layout] WARNING: node %" PRIu64 " carries no core length "
                     "(invariant I1 violated); treating as 0\n",
                     static_cast<uint64_t>(id));
        return 0;
    };

    auto to_end_inclusive = [](i64 start, i64 len) -> i64 {
        return (len <= 0) ? start : start + len - 1;
    };

    for (const auto &nit : recig) {
        if (visited.contains(nit.first)) continue;
        layouts.emplace_back();
        const i64 root_start = 0;
        const i64 root_len = infer_len_for_node(nit.first);
        const i64 root_end = to_end_inclusive(root_start, root_len);
        layouts.back().emplace(nit.first,
                               std::make_tuple(root_start, root_end, true));
        visited.emplace(nit.first);

        auto start_range = recig.equal_range(nit.first);
        for (auto it = start_range.first; it != start_range.second; ++it) {
            q.push(*it);
        }

        while (!q.empty()) {
            std::pair<u64, typename RecigMapT::mapped_type> node = q.front();
            q.pop();

            u64 src = node.first;
            u64 tgt = node.second.target;
            if (visited.contains(tgt)) continue;
            visited.emplace(tgt);

            auto [spos, epos, s_orient] = layouts.back().at(src);
            (void)epos;
            bool target_orient = !s_orient;

            long long qs =
                s_orient ? (node.second.query.len - node.second.query.to)
                         : node.second.query.from;
            long long ts =
                target_orient
                    ? (node.second.target_span.len - node.second.target_span.to)
                    : node.second.target_span.from;

            i64 delta = static_cast<i64>(qs - ts);
            i64 target_start = spos + delta;
            const i64 target_len = infer_len_for_node(tgt);
            const i64 target_end = to_end_inclusive(target_start, target_len);

            layouts.back().emplace(
                tgt, std::make_tuple(target_start, target_end, target_orient));

            auto range = recig.equal_range(tgt);
            for (auto it = range.first; it != range.second; ++it) {
                if (!visited.contains(it->second.target)) q.push(*it);
            }
        }
    }
    return layouts;
}

}  // namespace layout_graph_algo

namespace layout_graph_emit {

// ── Tunable constants ───────────────────────────────────────────────
// Reserve hint: expected seed-chain-DAG nodes (cores) per contig.
inline constexpr size_t SCD_CORES_PER_CONTIG_HINT = 64;
// Seed-chain shift gate (build_scd_graph). A chained seed's layout shift is
// where the layout put the contig's core minus where it put the matched node.
// Consecutive seeds of one true overlap keep their shift within
// SCD_SEED_SHIFT_TOL cores of each other; a run of fewer than
// SCD_MIN_SEED_RUN seeds set off by a larger jump is a core-label match
// between repeat copies and is not merged.
inline constexpr i64 SCD_SEED_SHIFT_TOL = 64;
inline constexpr size_t SCD_MIN_SEED_RUN = 2;
// Default max length (bp) of a post-layout C-node tip to prune.
inline constexpr size_t DEFAULT_POST_LAYOUT_MAX_TIP_LEN = 50000;
// Layout-boundary fuzz: a node counts as touching an end if it is within
// max(lcmer_l * BOUNDARY_FUZZ_LCMER_MULT, span / BOUNDARY_FUZZ_SPAN_DIV).
inline constexpr i64 BOUNDARY_FUZZ_LCMER_MULT = 2;
inline constexpr i64 BOUNDARY_FUZZ_SPAN_DIV = 20;

// Reverse-complement a sequence-fragment vector: reverse fragment order and
// reverse-complement each fragment (so the whole sequence is reverse-complemented).
inline gtl::vector<SequenceContainerType> reverse_complement_seqs(
    const gtl::vector<SequenceContainerType> &seqs) {
    gtl::vector<SequenceContainerType> rev;
    rev.reserve(seqs.size());
    for (auto it = seqs.rbegin(); it != seqs.rend(); ++it)
        rev.emplace_back(it->reverse_complement());
    return rev;
}

// Which end of a layout an inter-supernode edge attaches to.
inline constexpr uint8_t SIDE_START = 1;
inline constexpr uint8_t SIDE_END = 2;
inline const char *side_label(uint8_t side) {
    if (side == SIDE_START) return "start";
    if (side == SIDE_END) return "end";
    return "na";
}
inline char ori_char(bool ori) { return ori ? '+' : '-'; }

struct ultra_layout_node_info {
    u64 id;
    i64 from;
    i64 to;
    bool ori;
};

struct ultra_layout_job {
    u64 layout_idx;
    gtl::vector<ultra_layout_node_info> lnodes;
};

struct lane_node {
    u64 id = 0;
    i64 from = 0;
    i64 to = 0;
    bool is_forward_lane = true;
    const gtl::vector<SequenceContainerType> *raw_seqs = nullptr;
};

struct scd_node {
    uint32_t core_label;
    size_t contig_idx;
    size_t core_idx;
    uint16_t weight;
};

struct comp_path {
    i64 layout_from;
    i64 layout_to;
    int64_t total_weight;
    gtl::vector<uint32_t> path;
};

using scd_adj_t = gtl::vector<std::pair<uint32_t, uint16_t>>;

// 64-bit fingerprint of the homopolymer-compressed sequence of each core's
// span lr.pos[k] = [from, to) on the concatenation of `seqs`: one base per
// homopolymer and the number of homopolymers, not their lengths. Cores are
// built on homopolymer-compressed sequence, so copies of one core whose
// homopolymers differ in length (a common HiFi error) agree, while cores with
// equal labels but different compressed sequence (a 32-bit label collision,
// or repeat copies whose labels agree on the boundary children only) do not.
//
// One pass over the bases keeps the polynomial hash of all homopolymers so
// far; each core's first and last base record it, and the span's hash is the
// difference of the two prefixes. Memory is O(cores), not O(bases).
inline gtl::vector<u64> core_hpc_fingerprints(
    const gtl::vector<SequenceContainerType> &seqs, const lcp_read &lr) {
    constexpr u64 mul = 0x9E3779B97F4A7C15ULL;
    const size_t M = lr.pos.size();
    // (base, 2k) for core k's first base, (base, 2k + 1) for its last.
    gtl::vector<std::pair<uint32_t, uint32_t>> marks;
    marks.reserve(2 * M);
    for (uint32_t k = 0; k < M; ++k) {
        const auto [from, to] = lr.pos[k];
        if (from >= to) continue;
        marks.emplace_back(from, 2 * k);
        marks.emplace_back(to - 1, 2 * k + 1);
    }
    std::sort(marks.begin(), marks.end());

    // At a core's first base: hash of the runs before it and that run's
    // index. At its last base: hash of the runs through it and the run count.
    gtl::vector<u64> h_first(M, 0), h_last(M, 0);
    gtl::vector<uint32_t> r_first(M, 0), r_last(M, 0);
    u64 before = 0, through = 0;  // hash of runs [0, run - 1) / [0, run)
    uint32_t run = 0;             // runs started so far
    int prev = -1;
    uint32_t pos = 0;
    size_t mi = 0;
    for (const auto &f : seqs) {
        for (size_t b = 0; b < f.len && mi < marks.size(); ++b, ++pos) {
            const size_t bit = (static_cast<size_t>(f.skip_first_n) + b) * 2;
            const int c =
                static_cast<int>((f.data[bit >> 6] >> (bit & 63)) & 3);
            if (c != prev) {
                prev = c;
                before = through;
                through = through * mul + static_cast<u64>(c + 1);
                ++run;
            }
            for (; mi < marks.size() && marks[mi].first == pos; ++mi) {
                const uint32_t k = marks[mi].second >> 1;
                if (marks[mi].second & 1) {
                    h_last[k] = through;
                    r_last[k] = run;
                } else {
                    h_first[k] = before;
                    r_first[k] = run - 1;
                }
            }
        }
        if (mi == marks.size()) break;
    }

    gtl::vector<u64> fps(M, 0);
    for (size_t k = 0; k < M; ++k) {
        if (lr.pos[k].first >= lr.pos[k].second) continue;
        const uint32_t runs = r_last[k] - r_first[k];
        u64 p = 1;  // mul^runs
        for (u64 base = mul, e = runs; e; e >>= 1, base *= base)
            if (e & 1) p *= base;
        u64 h = h_last[k] - h_first[k] * p;
        h = (h ^ runs) * mul;
        fps[k] = h ^ (h >> 32);
    }
    return fps;
}

// Label-matched seeds seen by build_scd_graph, and those rejected because the
// two cores' homopolymer-compressed fingerprints differ.
struct scd_seed_stats {
    u64 seeds = 0;
    u64 rejected = 0;
};

// ── Phase 1: Build seed-chain DAG from contig cores ─────────────────
// With `core_fps` (core_hpc_fingerprints per contig, indexed like `nodes`),
// a label match becomes a seed only if both cores' fingerprints agree.
inline scd_seed_stats build_scd_graph(
    const gtl::vector<size_t> &contig_order,
    const gtl::flat_hash_map<size_t, lcp_read> &contig_lr,
    const gtl::vector<lane_node> &nodes,
    gtl::vector<scd_node> &scd,
    gtl::vector<scd_adj_t> &scd_fwd,
    gtl::flat_hash_map<uint32_t, gtl::vector<uint32_t>> &label_nodes,
    gtl::vector<gtl::vector<std::pair<size_t, size_t>>> &node_contigs,
    const gtl::vector<gtl::vector<u64>> *core_fps = nullptr)
{
    scd_seed_stats stats;
    scd.reserve(contig_order.size() * SCD_CORES_PER_CONTIG_HINT);
    scd_fwd.reserve(contig_order.size() * SCD_CORES_PER_CONTIG_HINT);
    node_contigs.reserve(contig_order.size() * SCD_CORES_PER_CONTIG_HINT);

    auto add_edge = [&](int32_t from, uint32_t to) {
        if (from < 0) return;
        for (auto &[t, ew] : scd_fwd[static_cast<size_t>(from)]) {
            if (t == to) { ++ew; return; }
        }
        scd_fwd[static_cast<size_t>(from)].push_back({to, 1});
    };

    for (size_t ci : contig_order) {
        const auto &lr = contig_lr.find(ci)->second;
        const size_t M = lr.cores.size();
        if (M == 0) continue;
        const size_t gn = scd.size();

        if (gn == 0) {
            int32_t last_node = -1;
            for (size_t k = 0; k < M; ++k) {
                uint32_t cur = static_cast<uint32_t>(scd.size());
                scd.push_back({lr.cores[k], ci, k, 1});
                scd_fwd.emplace_back();
                node_contigs.push_back({{ci, k}});
                label_nodes[lr.cores[k]].push_back(cur);
                add_edge(last_node, cur);
                last_node = static_cast<int32_t>(cur);
            }
            continue;
        }

        // Compute topological ranks (component-ordered)
        gtl::vector<uint32_t> topo_rank(gn, 0);
        {
            gtl::vector<uint32_t> in_deg(gn, 0);
            for (size_t u = 0; u < gn; ++u)
                for (auto &[v, w] : scd_fwd[u])
                    if (v < gn) ++in_deg[v];

            gtl::vector<uint32_t> uf(gn);
            for (uint32_t i = 0; i < gn; ++i) uf[i] = i;
            auto uf_find = [&](uint32_t x) -> uint32_t {
                while (uf[x] != x) { uf[x] = uf[uf[x]]; x = uf[x]; }
                return x;
            };
            for (size_t u = 0; u < gn; ++u)
                for (auto &[v, w] : scd_fwd[u])
                    if (v < gn) uf[uf_find(u)] = uf_find(v);

            gtl::flat_hash_map<uint32_t, gtl::vector<uint32_t>> comp_srcs;
            for (uint32_t i = 0; i < gn; ++i)
                if (in_deg[i] == 0) comp_srcs[uf_find(i)].push_back(i);

            gtl::vector<uint32_t> roots;
            roots.reserve(comp_srcs.size());
            for (auto &[r, _] : comp_srcs) roots.push_back(r);
            std::sort(roots.begin(), roots.end(), [&](uint32_t a, uint32_t b) {
                i64 ma = INT64_MAX, mb = INT64_MAX;
                for (uint32_t s : comp_srcs[a])
                    ma = std::min(ma, nodes[scd[s].contig_idx].from);
                for (uint32_t s : comp_srcs[b])
                    mb = std::min(mb, nodes[scd[s].contig_idx].from);
                return ma < mb;
            });

            uint32_t rank = 0;
            for (uint32_t root : roots) {
                std::queue<uint32_t> q;
                for (uint32_t s : comp_srcs[root]) q.push(s);
                while (!q.empty()) {
                    uint32_t u = q.front(); q.pop();
                    topo_rank[u] = rank++;
                    for (auto &[v, w] : scd_fwd[u])
                        if (v < gn && --in_deg[v] == 0) q.push(v);
                }
            }
        }

        // Collect seeds
        struct scd_seed { uint32_t rank; size_t k; uint32_t v; };
        gtl::vector<scd_seed> seeds;
        const u64 *ci_fps = core_fps ? (*core_fps)[ci].data() : nullptr;
        for (size_t k = 0; k < M; ++k) {
            auto it = label_nodes.find(lr.cores[k]);
            if (it == label_nodes.end()) continue;
            for (uint32_t v : it->second) {
                if (v >= gn) continue;
                ++stats.seeds;
                if (ci_fps) {
                    const auto &sv = scd[v];
                    if (ci_fps[k] != (*core_fps)[sv.contig_idx][sv.core_idx]) {
                        ++stats.rejected;
                        continue;
                    }
                }
                seeds.push_back({topo_rank[v], k, v});
            }
        }

        // Sort by rank asc, k desc for LIS
        std::sort(seeds.begin(), seeds.end(),
                  [](const scd_seed &a, const scd_seed &b) {
                      if (a.rank != b.rank) return a.rank < b.rank;
                      return a.k > b.k;
                  });

        // LIS on k with traceback
        const size_t S = seeds.size();
        gtl::vector<size_t> tails;
        gtl::vector<size_t> tail_idx;
        gtl::vector<int32_t> lis_prev(S, -1);
        for (size_t i = 0; i < S; ++i) {
            size_t k = seeds[i].k;
            auto pos = std::lower_bound(tails.begin(), tails.end(), k);
            size_t lp = static_cast<size_t>(pos - tails.begin());
            lis_prev[i] =
                (lp > 0) ? static_cast<int32_t>(tail_idx[lp - 1]) : -1;
            if (pos == tails.end()) {
                tails.push_back(k);
                tail_idx.push_back(i);
            } else {
                *pos = k;
                tail_idx[lp] = i;
            }
        }

        // Traceback
        gtl::vector<int32_t> query_to_node(M, -1);
        if (!tails.empty()) {
            for (int32_t idx = static_cast<int32_t>(tail_idx.back());
                 idx >= 0; idx = lis_prev[idx])
                query_to_node[seeds[idx].k] =
                    static_cast<int32_t>(seeds[idx].v);
        }

        // Unchain isolated seeds. Level-5 core labels hash only a core's
        // boundary children, so diverged repeat copies share labels; one such
        // seed fuses two layout positions into one DAG component, and the
        // heaviest path then drops whichever arm is lighter (Mb-scale losses,
        // sometimes behind a false join). Seeds of a true overlap come in runs
        // of near-constant layout shift; split the chain where the shift jumps
        // by more than SCD_SEED_SHIFT_TOL and drop runs shorter than
        // SCD_MIN_SEED_RUN.
        {
            gtl::vector<size_t> run;
            i64 prev_shift = 0;
            auto close_run = [&]() {
                if (run.size() < SCD_MIN_SEED_RUN)
                    for (size_t k : run) query_to_node[k] = -1;
                run.clear();
            };
            for (size_t k = 0; k < M; ++k) {
                if (query_to_node[k] < 0) continue;
                const auto &sv = scd[static_cast<size_t>(query_to_node[k])];
                const i64 shift = nodes[ci].from + static_cast<i64>(k) -
                                  nodes[sv.contig_idx].from -
                                  static_cast<i64>(sv.core_idx);
                if (!run.empty() &&
                    std::abs(shift - prev_shift) > SCD_SEED_SHIFT_TOL)
                    close_run();
                run.push_back(k);
                prev_shift = shift;
            }
            close_run();
        }

        // Merge matched, create new for unmatched
        int32_t last_node = -1;
        for (size_t k = 0; k < M; ++k) {
            uint32_t cur_node;
            if (query_to_node[k] >= 0) {
                cur_node = static_cast<uint32_t>(query_to_node[k]);
                scd[cur_node].weight++;
                node_contigs[cur_node].push_back({ci, k});
            } else {
                cur_node = static_cast<uint32_t>(scd.size());
                scd.push_back({lr.cores[k], ci, k, 1});
                scd_fwd.emplace_back();
                node_contigs.push_back({{ci, k}});
                label_nodes[lr.cores[k]].push_back(cur_node);
            }
            add_edge(last_node, cur_node);
            last_node = static_cast<int32_t>(cur_node);
        }
    }
    return stats;
}

// Compiled only under -DLASM_LAYOUT_RESOLVE_OVERLAP. Drops any SCD component
// overlapping a heavier one by >50% of the shorter interval, keeping the
// heavier. Off by default; disconnected components are instead emitted as
// separate contigs (see collect_path_sequences).
inline void resolve_overlapping_components(gtl::vector<comp_path> &comp_paths) {
    if (comp_paths.size() < 2) return;

    // Greedy weighted-interval selection: consider components heaviest-first,
    // accept one unless it duplicates (>50% of the shorter interval overlaps)
    // a component already accepted.
    gtl::vector<uint32_t> order(comp_paths.size());
    for (uint32_t i = 0; i < order.size(); ++i) order[i] = i;
    std::sort(order.begin(), order.end(), [&](uint32_t a, uint32_t b) {
        if (comp_paths[a].total_weight != comp_paths[b].total_weight)
            return comp_paths[a].total_weight > comp_paths[b].total_weight;
        return a < b;  // stable tie-break
    });

    gtl::vector<uint32_t> kept_idx;
    for (uint32_t idx : order) {
        const auto &cand = comp_paths[idx];
        const i64 cand_len = std::max<i64>(1, cand.layout_to - cand.layout_from);
        bool dup = false;
        for (uint32_t ki : kept_idx) {
            const auto &k = comp_paths[ki];
            const i64 lo = std::max(cand.layout_from, k.layout_from);
            const i64 hi = std::min(cand.layout_to, k.layout_to);
            const i64 ov = hi > lo ? hi - lo : 0;
            const i64 k_len = std::max<i64>(1, k.layout_to - k.layout_from);
            if (ov * 2 > std::min(cand_len, k_len)) { dup = true; break; }
        }
        if (!dup) kept_idx.push_back(idx);
    }

    if (kept_idx.size() == comp_paths.size()) return;  // nothing dropped

    // Re-emit the survivors in layout order (moving to avoid copying paths).
    std::sort(kept_idx.begin(), kept_idx.end(), [&](uint32_t a, uint32_t b) {
        return comp_paths[a].layout_from < comp_paths[b].layout_from;
    });
    gtl::vector<comp_path> kept;
    kept.reserve(kept_idx.size());
    for (uint32_t idx : kept_idx) kept.push_back(std::move(comp_paths[idx]));
    comp_paths = std::move(kept);
}

// ── Phase 2: Find heaviest path per connected component ─────────────
// Single longest-weight DP with an implicit virtual source (wired to every
// in-degree-0 node) and virtual sink (wired from every out-degree-0 node):
// dp[] realizes the source side, best_end the sink side. One source→sink
// traceback per component.
inline gtl::vector<comp_path> find_heaviest_paths(
    const gtl::vector<scd_node> &scd,
    const gtl::vector<scd_adj_t> &scd_fwd,
    const gtl::vector<lane_node> &nodes,
    u64 layout_idx)
{
    const size_t n = scd.size();

    gtl::vector<uint32_t> topo_order;
    topo_order.reserve(n);
    {
        gtl::vector<uint32_t> in_deg(n, 0);
        for (size_t u = 0; u < n; ++u)
            for (auto &[v, w] : scd_fwd[u]) ++in_deg[v];
        std::queue<uint32_t> q;
        for (uint32_t i = 0; i < n; ++i)
            if (in_deg[i] == 0) q.push(i);
        while (!q.empty()) {
            uint32_t u = q.front(); q.pop();
            topo_order.push_back(u);
            for (auto &[v, w] : scd_fwd[u])
                if (--in_deg[v] == 0) q.push(v);
        }
        if (topo_order.size() < n) {
            std::fprintf(stderr,
                         "[layout-dag] layout=%" PRIu64
                         " WARNING: %zu/%zu nodes stuck in cycles\n",
                         static_cast<uint64_t>(layout_idx),
                         n - topo_order.size(), n);
            gtl::flat_hash_set<uint32_t> visited(topo_order.begin(),
                                                 topo_order.end());
            for (uint32_t i = 0; i < n; ++i)
                if (!visited.contains(i)) topo_order.push_back(i);
        }
    }

    gtl::vector<scd_adj_t> scd_rev(n);
    for (size_t u = 0; u < n; ++u)
        for (auto &[v, w] : scd_fwd[u])
            scd_rev[v].push_back({static_cast<uint32_t>(u), w});

    gtl::vector<int64_t> dp(n, 0);
    gtl::vector<int32_t> parent(n, -1);
    for (uint32_t v : topo_order) {
        dp[v] = scd[v].weight;
        parent[v] = -1;
        for (auto &[u, w] : scd_rev[v]) {
            int64_t cand = dp[u] + scd[v].weight;
            if (cand > dp[v]) {
                dp[v] = cand;
                parent[v] = static_cast<int32_t>(u);
            }
        }
    }

    gtl::vector<int32_t> comp_id(n, -1);
    int32_t num_comps = 0;
    for (size_t i = 0; i < n; ++i) {
        if (comp_id[i] >= 0) continue;
        std::queue<uint32_t> bfs;
        bfs.push(static_cast<uint32_t>(i));
        comp_id[i] = num_comps;
        while (!bfs.empty()) {
            uint32_t u = bfs.front(); bfs.pop();
            for (auto &[v, w] : scd_fwd[u])
                if (comp_id[v] < 0) { comp_id[v] = num_comps; bfs.push(v); }
            for (auto &[u2, w] : scd_rev[u])
                if (comp_id[u2] < 0) { comp_id[u2] = num_comps; bfs.push(u2); }
        }
        ++num_comps;
    }

    gtl::vector<comp_path> comp_paths;
    for (int32_t c = 0; c < num_comps; ++c) {
        int32_t best_end = -1;
        int64_t best_dp = -1;
        for (size_t i = 0; i < n; ++i) {
            if (comp_id[i] != c) continue;
            if (scd_fwd[i].empty() && dp[i] > best_dp) {
                best_dp = dp[i];
                best_end = static_cast<int32_t>(i);
            }
        }
        if (best_end < 0) {
            for (size_t i = 0; i < n; ++i) {
                if (comp_id[i] != c) continue;
                if (dp[i] > best_dp) {
                    best_dp = dp[i];
                    best_end = static_cast<int32_t>(i);
                }
            }
        }
        if (best_end < 0) continue;

        // Traceback the heaviest path ending at best_end.
        gtl::vector<uint32_t> path;
        for (int32_t cur = best_end; cur >= 0; cur = parent[cur])
            path.push_back(static_cast<uint32_t>(cur));
        std::reverse(path.begin(), path.end());

        // Invariant I3: consecutive path nodes share an SCD edge (parent[]
        // points only at a real predecessor). Verified defensively.
        for (size_t pi = 1; pi < path.size(); ++pi) {
            bool linked = false;
            for (auto &[v, w] : scd_fwd[path[pi - 1]])
                if (v == path[pi]) { linked = true; break; }
            if (!linked) [[unlikely]]
                std::fprintf(stderr,
                    "[layout-dag] layout=%" PRIu64 " WARNING: emitted path not "
                    "edge-connected at index %zu (%s:%d)\n",
                    static_cast<uint64_t>(layout_idx), pi, __FILE__, __LINE__);
        }

        i64 min_pos = INT64_MAX, max_pos = INT64_MIN;
        for (size_t i = 0; i < n; ++i) {
            if (comp_id[i] != c) continue;
            min_pos = std::min(min_pos, nodes[scd[i].contig_idx].from);
            max_pos = std::max(max_pos, nodes[scd[i].contig_idx].to);
        }
        comp_paths.push_back({min_pos, max_pos, best_dp, std::move(path)});
    }

    std::sort(comp_paths.begin(), comp_paths.end(),
              [](const auto &a, const auto &b) {
                  return a.layout_from < b.layout_from;
              });
#ifdef LASM_LAYOUT_RESOLVE_OVERLAP
    // Off by default; see the note on resolve_overlapping_components.
    resolve_overlapping_components(comp_paths);
#endif

    return comp_paths;
}

// ── Phase 3: Collect output sequences from heaviest paths ───────────
// Returns one sequence group per SCD component; groups are never concatenated.
// The caller emits each as its own contig.
inline gtl::vector<gtl::vector<SequenceContainerType>> collect_path_sequences(
    const gtl::vector<comp_path> &paths,
    const gtl::vector<gtl::vector<std::pair<size_t, size_t>>> &node_contigs,
    gtl::flat_hash_map<size_t, gtl::vector<SequenceContainerType>> &contig_core_seqs,
    size_t n_nodes)
{
    gtl::vector<gtl::vector<SequenceContainerType>> comp_seqs;
    comp_seqs.reserve(paths.size());

    for (size_t ci = 0; ci < paths.size(); ++ci) {
        gtl::vector<SequenceContainerType> merged_seqs;
        merged_seqs.reserve(n_nodes);

        const auto &path = paths[ci].path;
        size_t emit_contig = SIZE_MAX;

        for (size_t pi = 0; pi < path.size(); ++pi) {
            uint32_t nid = path[pi];
            const auto &nc = node_contigs[nid];

            size_t use_contig = SIZE_MAX;
            size_t use_core = SIZE_MAX;
            if (emit_contig != SIZE_MAX) {
                for (auto &[cid, kid] : nc) {
                    if (cid == emit_contig) {
                        use_contig = cid;
                        use_core = kid;
                        break;
                    }
                }
            }

            if (use_contig == SIZE_MAX) {
                // Prefer a contig shared with the next node so consecutive
                // path nodes stay on the same contig; otherwise pick the
                // longest available core sequence.
                auto pick_contig = [&]() -> std::pair<size_t, size_t> {
                    if (pi + 1 < path.size()) {
                        const auto &nc_next = node_contigs[path[pi + 1]];
                        for (auto &[cid, kid] : nc) {
                            for (auto &[cid2, kid2] : nc_next) {
                                if (cid == cid2) return {cid, kid};
                            }
                        }
                    }
                    size_t best_cid = SIZE_MAX, best_kid = SIZE_MAX,
                           best_len = 0;
                    for (auto &[cid, kid] : nc) {
                        auto it = contig_core_seqs.find(cid);
                        size_t len =
                            (it != contig_core_seqs.end() &&
                             kid < it->second.size())
                                ? it->second[kid].len
                                : 0;
                        if (best_cid == SIZE_MAX || len > best_len) {
                            best_cid = cid;
                            best_kid = kid;
                            best_len = len;
                        }
                    }
                    return {best_cid, best_kid};
                };
                const auto [picked_contig, picked_core] = pick_contig();
                use_contig = picked_contig;
                use_core = picked_core;
            }

            if (use_contig == SIZE_MAX) continue;
            emit_contig = use_contig;

            auto cs_it = contig_core_seqs.find(use_contig);
            if (cs_it != contig_core_seqs.end() &&
                use_core < cs_it->second.size()) {
                auto &core_seq = cs_it->second[use_core];
                if (core_seq.len > 0)
                    merged_seqs.push_back(std::move(core_seq));
            }
        }
        if (!merged_seqs.empty()) comp_seqs.push_back(std::move(merged_seqs));
    }
    return comp_seqs;
}

// ── Main entry: build layout sequence from positioned nodes ─────────
template <class AssemblerT>
// Returns one sequence group per SCD component (see collect_path_sequences):
// disconnected components are never concatenated. Adds this layout's SCD seed
// counts to `seed_stats`.
inline gtl::vector<gtl::vector<SequenceContainerType>>
build_layout_sequence_from_nodes(
    [[maybe_unused]] AssemblerT &as, dbg_index &dbg,
    const gtl::vector<ultra_layout_node_info> &lnodes, u64 layout_idx,
    scd_seed_stats &seed_stats) {
    if (lnodes.empty()) return {};

    // ── Build lane_nodes ────────────────────────────────────────────
    gtl::vector<lane_node> nodes;
    nodes.reserve(lnodes.size());
    for (const auto &n : lnodes) {
        auto it = dbg.dbg.find(n.id);
        if (it == dbg.dbg.end()) continue;
        if (it->second.seq.empty()) continue;
        nodes.push_back({n.id, n.from, n.to, n.ori, &it->second.seq});
    }
    if (nodes.empty()) return {};

    // Sort by layout position
    gtl::vector<size_t> sorted_idx(nodes.size());
    for (size_t i = 0; i < nodes.size(); ++i) sorted_idx[i] = i;
    std::sort(sorted_idx.begin(), sorted_idx.end(), [&](size_t ia, size_t ib) {
        const auto &a = nodes[ia];
        const auto &b = nodes[ib];
        if (a.from != b.from) return a.from < b.from;
        if (a.to != b.to) return a.to < b.to;
        return a.id < b.id;
    });

    // ── Fast path: single contig or 1fw+1rv → return the larger ──
    {
        size_t n_fw = 0, n_rv = 0, fi = SIZE_MAX, ri = SIZE_MAX;
        for (size_t i = 0; i < nodes.size(); ++i) {
            if (nodes[i].is_forward_lane) { ++n_fw; fi = i; }
            else                          { ++n_rv; ri = i; }
        }
        if (n_fw <= 1 && n_rv <= 1 && (n_fw + n_rv) <= 2) {
            size_t fw_len = 0, rv_len = 0;
            if (fi != SIZE_MAX)
                for (const auto &s : *nodes[fi].raw_seqs) fw_len += s.len;
            if (ri != SIZE_MAX)
                for (const auto &s : *nodes[ri].raw_seqs) rv_len += s.len;

            size_t pick = (fw_len >= rv_len && fi != SIZE_MAX) ? fi : ri;
            if (pick == SIZE_MAX) return {};

            gtl::vector<gtl::vector<SequenceContainerType>> out;
            out.reserve(1);
            if (!nodes[pick].is_forward_lane) {
                out.push_back(reverse_complement_seqs(*nodes[pick].raw_seqs));
                return out;
            }
            const auto &raw = *nodes[pick].raw_seqs;
            gtl::vector<SequenceContainerType> result;
            result.reserve(1);
            result.emplace_back(raw.begin(), raw.end());
            out.push_back(std::move(result));
            return out;
        }
    }

    // ── Build per-contig core sequences ─────────────────────────────
    gtl::flat_hash_map<size_t, gtl::vector<SequenceContainerType>> contig_core_seqs;
    gtl::flat_hash_map<size_t, lcp_read> contig_lr;
    gtl::vector<size_t> contig_order;
    contig_order.reserve(sorted_idx.size());
    core_buf lps_buf;
    // --strict-layout-labels: per-contig homopolymer-compressed core
    // fingerprints, indexed like `nodes`.
    gtl::vector<gtl::vector<u64>> core_fps;
    if (as.strict_layout_labels) core_fps.resize(nodes.size());

    for (size_t ci : sorted_idx) {
        const auto &raw = *nodes[ci].raw_seqs;
        gtl::vector<SequenceContainerType> rev_seqs;
        const gtl::vector<SequenceContainerType> *seqs_ptr = &raw;
        if (!nodes[ci].is_forward_lane) {
            rev_seqs = reverse_complement_seqs(raw);
            seqs_ptr = &rev_seqs;
        }

        auto [lr, cseqs] = AssemblerT::build_lr_cores_with_seqs(
            *seqs_ptr, as.syncmer_functor, as.lcp_level, lps_buf);

        if (lr.cores.empty()) {
            const auto &frags = *seqs_ptr;
            size_t total_len = 0;
            for (const auto &f : frags) total_len += f.len;
            if (total_len == 0) continue;
            SequenceContainerType synth(frags.begin(), frags.end());
            gtl::vector<SequenceContainerType> synth_vec;
            synth_vec.push_back(std::move(synth));
            lcp_read synth_lr;
            synth_lr.cores.push_back(0);
            synth_lr.pos.emplace_back(0u, static_cast<uint32_t>(total_len));
            contig_lr[ci] = std::move(synth_lr);
            if (as.strict_layout_labels)
                core_fps[ci] = core_hpc_fingerprints(frags, contig_lr[ci]);
            contig_core_seqs[ci] = std::move(synth_vec);
            contig_order.push_back(ci);
            continue;
        }
        if (as.strict_layout_labels)
            core_fps[ci] = core_hpc_fingerprints(*seqs_ptr, lr);
        contig_lr[ci] = std::move(lr);
        contig_core_seqs[ci] = std::move(cseqs);
        contig_order.push_back(ci);
    }

    if (contig_order.empty()) {
        std::fprintf(stderr, "[layout-dag] layout=%" PRIu64
                     " SKIPPED: no contigs with cores\n",
                     static_cast<uint64_t>(layout_idx));
        return {};
    }

    // ── Phase 1: Build seed-chain DAG ───────────────────────────────
    gtl::vector<scd_node> scd;
    gtl::vector<scd_adj_t> scd_fwd;
    gtl::flat_hash_map<uint32_t, gtl::vector<uint32_t>> label_nodes;
    gtl::vector<gtl::vector<std::pair<size_t, size_t>>> node_contigs;
    const scd_seed_stats stats = build_scd_graph(
        contig_order, contig_lr, nodes, scd, scd_fwd, label_nodes, node_contigs,
        as.strict_layout_labels ? &core_fps : nullptr);
    seed_stats.seeds += stats.seeds;
    seed_stats.rejected += stats.rejected;

    if (scd.empty()) {
        std::fprintf(stderr, "[layout-dag] layout=%" PRIu64
                     " SKIPPED: SCD graph empty\n",
                     static_cast<uint64_t>(layout_idx));
        return {};
    }

    // ── Phase 2: Find heaviest paths ────────────────────────────────
    auto comp_paths = find_heaviest_paths(scd, scd_fwd, nodes, layout_idx);

    // ── Phase 3: Collect output sequences ───────────────────────────
    auto result = collect_path_sequences(
        comp_paths, node_contigs, contig_core_seqs, scd.size());


    return result;
}
// One SCD component's sequences …
using layout_comp_seqs_t = gtl::vector<SequenceContainerType>;
// … and all components of one layout (emitted as separate contigs).
using layout_seq_t = gtl::vector<layout_comp_seqs_t>;

template <class AssemblerT>
inline gtl::vector<layout_seq_t> build_layout_sequences_parallel(
    AssemblerT &as, dbg_index &dbg, const gtl::vector<ultra_layout_job> &jobs) {
    gtl::vector<layout_seq_t> out_seqs(jobs.size());
    if (jobs.empty()) return out_seqs;

    using job_q_t = moodycamel::BlockingConcurrentQueue<size_t>;
    job_q_t job_queue;
    std::atomic<bool> done{false};

    int nthreads = as.threads.work;
    if (nthreads < 1) nthreads = 1;
    nthreads = std::min<int>(nthreads, jobs.size());

    gtl::vector<std::thread> workers;
    workers.reserve(nthreads);
    gtl::vector<scd_seed_stats> worker_stats(static_cast<size_t>(nthreads));
    for (int tid = 0; tid < nthreads; ++tid) {
        workers.emplace_back([&, tid]() {
            scd_seed_stats &stats = worker_stats[static_cast<size_t>(tid)];
            size_t job_idx = 0;
            while (true) {
                if (!job_queue.wait_dequeue_timed(job_idx,
                                                  WAIT_QUEUE_TIMEOUT)) {
                    if (done.load()) break;
                    continue;
                }
                out_seqs[job_idx] = build_layout_sequence_from_nodes(
                    as, dbg, jobs[job_idx].lnodes, jobs[job_idx].layout_idx,
                    stats);
            }
            while (job_queue.try_dequeue(job_idx)) {
                out_seqs[job_idx] = build_layout_sequence_from_nodes(
                    as, dbg, jobs[job_idx].lnodes, jobs[job_idx].layout_idx,
                    stats);
            }
        });
    }

    for (size_t i = 0; i < jobs.size(); ++i) {
        while (!job_queue.try_enqueue(i)) {
            std::this_thread::sleep_for(std::chrono::nanoseconds(10));
        }
    }
    done.store(true);
    for (auto &worker : workers) {
        worker.join();
    }
    if (as.strict_layout_labels) {
        scd_seed_stats total;
        for (const auto &s : worker_stats) {
            total.seeds += s.seeds;
            total.rejected += s.rejected;
        }
        std::fprintf(stderr,
                     "[layout-dag] strict-layout-labels: rejected %" PRIu64
                     " of %" PRIu64 " label-matched SCD seeds\n",
                     total.rejected, total.seeds);
    }

    return out_seqs;
}

// ── Post-layout graph types ─────────────────────────────────────────
struct supernode_out {
    std::string name;
    gtl::vector<SequenceContainerType> owned_seqs;  // layout seqs (owned)
    const gtl::vector<SequenceContainerType> *ref_seqs =
        nullptr;  // unitig seqs (borrowed)
    bool is_layout;
    u64 aux_id;

    const gtl::vector<SequenceContainerType> &seqs() const {
        return ref_seqs ? *ref_seqs : owned_seqs;
    }
};

struct edge_agg {
    u64 src_super = 0;
    u64 dst_super = 0;
    char src_ori = '+';
    char dst_ori = '+';
    uint8_t src_side = 0;
    uint8_t dst_side = 0;
    u64 example_from = 0;
    u64 example_to = 0;
    size_t count = 0;
};

using edge_key = std::tuple<u64, u64, char, char, uint8_t, uint8_t>;
using supernode_vec = gtl::vector<supernode_out>;
using edge_map = std::map<edge_key, edge_agg>;

// ── Post-layout graph cleaning ──────────────────────────────────────
struct PostLayoutCleaning {
    struct adj_entry {
        u64 neighbor;  // supernode id
        edge_key key;  // key into merged_edges
    };

    // Per-supernode adjacency (directed: fwd = outgoing, rev = incoming)
    struct adj_info {
        gtl::vector<adj_entry> fwd;  // outgoing edges
        gtl::vector<adj_entry> rev;  // incoming edges
    };

    using adj_map = gtl::flat_hash_map<u64, adj_info>;

    static void build_adjacency(const edge_map &merged_edges, adj_map &adj) {
        adj.clear();
        for (const auto &[k, e] : merged_edges) {
            adj[e.src_super].fwd.push_back({e.dst_super, k});
            adj[e.dst_super].rev.push_back({e.src_super, k});
        }
    }

    // Drop the reciprocal adjacency entry for a removed edge from a neighbor's
    // fwd/rev list (matched by unique edge key).
    static void drop_adj_edge(gtl::vector<adj_entry> &lst, const edge_key &k) {
        lst.erase(std::remove_if(lst.begin(), lst.end(),
                                 [&](const adj_entry &a) { return a.key == k; }),
                  lst.end());
    }

    // Remove short C-node tips: C-nodes with no successors (or no
    // predecessors) and sequence length below max_tip_len. Each pass decides
    // removals from a frozen snapshot of `adj`, then applies them while
    // patching neighbor adjacency incrementally.
    static size_t cut_tips(supernode_vec &supernodes, edge_map &merged_edges,
                           adj_map &adj, size_t max_tip_len) {
        size_t removed = 0;
        bool changed = true;
        gtl::vector<size_t> to_remove;
        while (changed) {
            changed = false;
            to_remove.clear();

            // Decide removals from the current (frozen) adjacency snapshot.
            for (size_t i = 0; i < supernodes.size(); ++i) {
                auto &sn = supernodes[i];
                if (sn.name.empty()) continue;  // already removed
                if (sn.is_layout) continue;     // never remove M-nodes

                auto adj_it = adj.find(i);
                size_t in_deg =
                    (adj_it != adj.end()) ? adj_it->second.rev.size() : 0;
                size_t out_deg =
                    (adj_it != adj.end()) ? adj_it->second.fwd.size() : 0;

                // Tip: dead-end on at least one side
                if (in_deg != 0 && out_deg != 0) continue;
                if (in_deg == 0 && out_deg == 0) continue;  // isolated — keep

                // Check sequence length
                size_t seq_len = 0;
                for (const auto &s : sn.seqs()) seq_len += s.len;
                if (seq_len >= max_tip_len) continue;

                to_remove.push_back(i);
            }

            // Apply removals: erase edges and patch each neighbor's reciprocal
            // adjacency entry so the next pass sees updated degrees.
            for (size_t i : to_remove) {
                auto adj_it = adj.find(i);
                if (adj_it != adj.end()) {
                    for (const auto &e : adj_it->second.fwd) {  // i → neighbor
                        merged_edges.erase(e.key);
                        if (e.neighbor == i) continue;  // self-loop
                        auto nit = adj.find(e.neighbor);
                        if (nit != adj.end())
                            drop_adj_edge(nit->second.rev, e.key);
                    }
                    for (const auto &e : adj_it->second.rev) {  // neighbor → i
                        merged_edges.erase(e.key);
                        if (e.neighbor == i) continue;  // self-loop
                        auto nit = adj.find(e.neighbor);
                        if (nit != adj.end())
                            drop_adj_edge(nit->second.fwd, e.key);
                    }
                }
                supernodes[i].name.clear();
                adj.erase(i);
                changed = true;
                ++removed;
            }
        }
        return removed;
    }

    // Remove edges pointing to/from removed supernodes.
    static void fix_edges(const supernode_vec &supernodes,
                          edge_map &merged_edges) {
        std::erase_if(merged_edges, [&](const auto &kv) {
            const auto &e = kv.second;
            return e.src_super >= supernodes.size() ||
                   e.dst_super >= supernodes.size() ||
                   supernodes[e.src_super].name.empty() ||
                   supernodes[e.dst_super].name.empty();
        });
    }

    static void run(supernode_vec &supernodes, edge_map &merged_edges,
                    size_t max_tip_len = DEFAULT_POST_LAYOUT_MAX_TIP_LEN) {
        adj_map adj;
        build_adjacency(merged_edges, adj);

        size_t tips = cut_tips(supernodes, merged_edges, adj, max_tip_len);
        fix_edges(supernodes, merged_edges);

        if (tips > 0)
            std::fprintf(stderr, "[post-layout-clean] Removed %zu tips\n",
                         tips);
    }
};

// ── Orientation resolution ──────────────────────────────────────────
// Each layout assigns every DBG node an orientation (+/-). Normally a node
// belongs to exactly one layout, so nothing needs resolving. `*out_conflicts`
// counts any node claimed by two layouts with opposite orientation (a
// violation of the single-owner invariant).
enum class orientation_policy {
    first_wins,     // first layout (in `layouts` order) to claim the node wins
    majority_vote,  // orientation backed by the most layouts wins (tie → '+')
};

inline gtl::flat_hash_map<u64, bool> resolve_layout_orientations(
    const layout_graph_algo::layout_result_t &layouts,
    orientation_policy policy = orientation_policy::first_wins,
    size_t *out_conflicts = nullptr) {
    gtl::flat_hash_map<u64, std::pair<uint32_t, uint32_t>> votes;  // {plus,minus}
    gtl::flat_hash_map<u64, bool> first_ori;                       // first-seen
    for (const auto &layout : layouts) {
        for (const auto &[id, t] : layout) {
            const bool ori = std::get<2>(t);
            auto &v = votes[id];
            if (v.first == 0 && v.second == 0) first_ori.emplace(id, ori);
            if (ori) ++v.first; else ++v.second;
        }
    }

    size_t conflicts = 0;
    gtl::flat_hash_map<u64, bool> result;
    result.reserve(votes.size());
    for (const auto &[id, v] : votes) {
        if (v.first > 0 && v.second > 0) ++conflicts;
        bool ori;
        switch (policy) {
            case orientation_policy::majority_vote:
                ori = (v.first >= v.second);  // tie → '+'
                break;
            case orientation_policy::first_wins:
            default:
                ori = first_ori[id];
                break;
        }
        result[id] = ori;
    }
    if (out_conflicts) *out_conflicts = conflicts;
    return result;
}

// ── Phase: build one layout job per component + per-node boundary masks ──
// A node's boundary mask records whether it sits at the start and/or end of its
// layout (within a fuzz window), which later gates inter-supernode edges.
template <class AssemblerT>
inline void build_layout_jobs(
    AssemblerT &as, const layout_graph_algo::layout_result_t &layouts,
    gtl::vector<ultra_layout_job> &layout_jobs,
    gtl::flat_hash_set<u64> &nodes_in_layouts,
    gtl::flat_hash_map<u64, uint8_t> &node_boundary_mask) {
    u64 layout_idx = 0;
    for (const auto &current_layout : layouts) {
        gtl::vector<ultra_layout_node_info> lnodes;
        lnodes.reserve(current_layout.size());
        i64 min_start = INT64_MAX;
        i64 max_end = INT64_MIN;
        for (const auto &[ctg_id, layout] : current_layout) {
            auto [from, to, ori] = layout;
            lnodes.push_back({ctg_id, from, to, ori});
            min_start = std::min(min_start, from);
            max_end = std::max(max_end, to);
        }
        if (lnodes.empty()) {
            ++layout_idx;
            continue;
        }
        std::sort(lnodes.begin(), lnodes.end(),
                  [](const ultra_layout_node_info &a,
                     const ultra_layout_node_info &b) {
                      if (a.from != b.from) return a.from < b.from;
                      return a.id < b.id;
                  });

        const i64 span = std::max<i64>(1, max_end - min_start);
        const i64 fuzzy = std::max<i64>(as.lcmer_l * BOUNDARY_FUZZ_LCMER_MULT,
                                        span / BOUNDARY_FUZZ_SPAN_DIV);
        for (const auto &n : lnodes) {
            uint8_t mask = 0;
            if (n.from <= min_start + fuzzy) mask |= SIDE_START;
            if (n.to >= max_end - fuzzy) mask |= SIDE_END;
            if (mask == 0) {
                i64 d_start = n.from - min_start;
                if (d_start < 0) d_start = -d_start;
                i64 d_end = max_end - n.to;
                if (d_end < 0) d_end = -d_end;
                mask = (d_start <= d_end) ? SIDE_START : SIDE_END;
            }
            node_boundary_mask[n.id] |= mask;
            nodes_in_layouts.insert(n.id);
        }

        layout_jobs.push_back({layout_idx, std::move(lnodes)});
        ++layout_idx;
    }
}

// ── Phase: find C-node components fully contained inside layouts ─────────
// Populates `candidate_set` with the non-layout ("C") nodes, then returns the
// subset whose connected component (over C-nodes) never touches a layout
// boundary or an external non-layout node.
inline gtl::flat_hash_set<u64> find_contained_cnodes(
    dbg_index &dbg, const gtl::flat_hash_set<u64> &nodes_in_layouts,
    const gtl::flat_hash_map<u64, uint8_t> &node_boundary_mask,
    const gtl::flat_hash_map<u64, bool> &node_layout_ori,
    gtl::flat_hash_set<u64> &candidate_set) {
    // Successor map from prev edges.
    gtl::flat_hash_map<u64, gtl::vector<u64>> successors;
    for (const auto &[to_id, to_node] : dbg.dbg) {
        if (to_node.redirected()) continue;
        for (u64 from_id : to_node.prev) successors[from_id].push_back(to_id);
    }

    for (u64 node_id : dbg.active_set) {
        if (nodes_in_layouts.contains(node_id)) continue;
        auto it = dbg.dbg.find(node_id);
        if (it == dbg.dbg.end()) continue;
        if (it->second.redirected() || it->second.seq.empty()) continue;
        candidate_set.insert(node_id);
    }

    gtl::flat_hash_set<u64> cc_visited;
    gtl::flat_hash_set<u64> contained_nodes;
    size_t standalone_components = 0;

    for (u64 seed : candidate_set) {
        if (cc_visited.contains(seed)) continue;

        gtl::vector<u64> component;
        std::queue<u64> bfs;
        bfs.push(seed);
        cc_visited.insert(seed);
        bool has_boundary_neighbor = false;
        bool has_external_edge = false;

        while (!bfs.empty()) {
            u64 cur = bfs.front();
            bfs.pop();
            component.push_back(cur);
            auto db_it = dbg.dbg.find(cur);

            for (u64 pred_id : db_it->second.prev) {
                if (candidate_set.contains(pred_id)) {
                    if (!cc_visited.contains(pred_id)) {
                        cc_visited.insert(pred_id);
                        bfs.push(pred_id);
                    }
                } else {
                    has_external_edge = true;
                    if (!nodes_in_layouts.contains(pred_id)) {
                        has_boundary_neighbor = true;
                    } else {
                        auto bm_it = node_boundary_mask.find(pred_id);
                        if (bm_it != node_boundary_mask.end()) {
                            bool ori = node_layout_ori.at(pred_id);
                            uint8_t needed = ori ? SIDE_END : SIDE_START;
                            if (bm_it->second & needed)
                                has_boundary_neighbor = true;
                        }
                    }
                }
            }
            for (u64 succ_id : successors[cur]) {
                if (candidate_set.contains(succ_id)) {
                    if (!cc_visited.contains(succ_id)) {
                        cc_visited.insert(succ_id);
                        bfs.push(succ_id);
                    }
                } else {
                    has_external_edge = true;
                    if (!nodes_in_layouts.contains(succ_id)) {
                        has_boundary_neighbor = true;
                    } else {
                        auto bm_it = node_boundary_mask.find(succ_id);
                        if (bm_it != node_boundary_mask.end()) {
                            bool ori = node_layout_ori.at(succ_id);
                            uint8_t needed = ori ? SIDE_START : SIDE_END;
                            if (bm_it->second & needed)
                                has_boundary_neighbor = true;
                        }
                    }
                }
            }
        }

        // Contained: has external adjacency and none reaches a layout
        // boundary (hangs off layout interiors). A component with no external
        // adjacency is standalone and must be kept.
        if (!has_boundary_neighbor && has_external_edge) {
            for (u64 nid : component) contained_nodes.insert(nid);
        } else if (!has_external_edge && !component.empty()) {
            ++standalone_components;
        }
    }

    size_t contained_bp = 0;
    for (u64 node_id : contained_nodes) {
        auto it = dbg.dbg.find(node_id);
        for (const auto &s : it->second.seq) contained_bp += s.len;
    }
    if (!contained_nodes.empty() || standalone_components > 0)
        std::fprintf(stderr,
                     "[layout-emit] Skipped %zu contained C-nodes (%zu bp); "
                     "kept %zu standalone components\n",
                     contained_nodes.size(), contained_bp,
                     standalone_components);
    return contained_nodes;
}

// ── Phase: aggregate DBG edges into inter-supernode edges ───────────────
// Parallel over the DBG; boundary/orientation gate which edges survive. Edges
// between nodes of the same supernode are dropped (except self-loops).
inline edge_map aggregate_superedges(
    dbg_index &dbg, const gtl::flat_hash_map<u64, u64> &node_to_super,
    const gtl::flat_hash_set<u64> &nodes_in_layouts,
    const gtl::flat_hash_map<u64, bool> &node_layout_ori,
    const gtl::flat_hash_map<u64, uint8_t> &node_boundary_mask, int nthreads) {
    if (nthreads < 1) nthreads = 1;
    gtl::vector<std::map<edge_key, edge_agg>> per_thread(nthreads);

    for_each_map(
        dbg.dbg,
        [&](int tid, const auto &entry) {
            const auto &[to_id, to_node] = entry;
            if (to_node.redirected()) return;
            auto to_super_it = node_to_super.find(to_id);
            if (to_super_it == node_to_super.end()) return;

            for (u64 from_id : to_node.prev) {
                auto from_super_it = node_to_super.find(from_id);
                if (from_super_it == node_to_super.end()) continue;
                // Skip intra-supernode edges, but keep self-loops
                if (from_super_it->second == to_super_it->second &&
                    from_id != to_id)
                    continue;

                bool accept = true;
                uint8_t src_side = 0, dst_side = 0;
                char src_ori = '+', dst_ori = '+';

                if (nodes_in_layouts.contains(from_id)) {
                    bool ori = node_layout_ori.at(from_id);
                    src_ori = ori_char(ori);
                    src_side = ori ? SIDE_END : SIDE_START;
                    if ((node_boundary_mask.at(from_id) & src_side) == 0)
                        accept = false;
                }
                if (nodes_in_layouts.contains(to_id)) {
                    bool ori = node_layout_ori.at(to_id);
                    dst_ori = ori_char(ori);
                    dst_side = ori ? SIDE_START : SIDE_END;
                    if ((node_boundary_mask.at(to_id) & dst_side) == 0)
                        accept = false;
                }
                if (!accept) continue;

                edge_key k{from_super_it->second, to_super_it->second,
                           src_ori,             dst_ori,
                           src_side,            dst_side};
                auto &agg = per_thread[tid][k];
                if (agg.count == 0) {
                    agg.src_super = from_super_it->second;
                    agg.dst_super = to_super_it->second;
                    agg.src_ori = src_ori;
                    agg.dst_ori = dst_ori;
                    agg.src_side = src_side;
                    agg.dst_side = dst_side;
                    agg.example_from = from_id;
                    agg.example_to = to_id;
                }
                agg.count++;
            }
        },
        nthreads);

    edge_map merged_edges;
    for (auto &local : per_thread)
        for (auto &[k, agg] : local) {
            auto &dst = merged_edges[k];
            if (dst.count == 0)
                dst = agg;
            else
                dst.count += agg.count;
        }
    return merged_edges;
}

// ── Phase: format supernode S-lines (parallel) and write the GFA ────────
inline void write_gfa(const supernode_vec &supernodes,
                      const edge_map &merged_edges, int nthreads, FILE *out) {
    // Each worker converts a supernode's sequences to ASCII and builds the
    // complete `S` line (prefix + sequence + tags) off the serial write path.
    const size_t n_super = supernodes.size();
    gtl::vector<std::string> super_line_strs(n_super);
    {
        int nthreads_print = nthreads < 1 ? 1 : nthreads;
        nthreads_print = std::min<int>(nthreads_print, n_super);

        std::atomic<size_t> next_idx{0};
        gtl::vector<std::thread> workers;
        workers.reserve(nthreads_print);
        for (int t = 0; t < nthreads_print; ++t) {
            workers.emplace_back([&]() {
                charsi local_chs;
                while (true) {
                    size_t i = next_idx.fetch_add(1);
                    if (i >= n_super) break;
                    const auto &sn = supernodes[i];
                    if (sn.name.empty()) continue;  // removed
                    const auto &sn_seqs = sn.seqs();
                    std::string &dst = super_line_strs[i];
                    size_t total_len = 0;
                    for (const auto &s : sn_seqs) total_len += s.len;

                    dst.reserve(total_len + sn.name.size() + 48);
                    dst.append("S\t");
                    dst.append(sn.name);
                    dst.push_back('\t');

                    const size_t seq_start = dst.size();
                    for (const auto &s : sn_seqs) {
                        // No gap fabrication: skip zero-length entries rather
                        // than emitting a run of 'N'.
                        if (s.len == 0) continue;
                        local_chs.require(s.len + 32);
                        s.convert_to_ascii((char *)local_chs.data);
                        dst.append((const char *)local_chs.data, s.len);
                    }
                    if (dst.size() == seq_start) dst.push_back('*');  // no seq

                    dst.append("\tTP:Z:");
                    dst.append(sn.is_layout ? "layout" : "unitig");
                    dst.push_back('\t');
                    dst.append(sn.is_layout ? "LY" : "UI");
                    dst.append(":i:");
                    dst.append(std::to_string(static_cast<uint64_t>(sn.aux_id)));
                    dst.push_back('\n');
                }
            });
        }
        for (auto &w : workers) w.join();
    }

    // Accumulate output into a 4 MiB buffer, flushing in large fwrites.
    // Oversized lines (megabase S-lines) bypass the buffer with a direct write.
    constexpr size_t FLUSH_AT = size_t{1} << 22;  // 4 MiB
    std::string wbuf;
    wbuf.reserve(FLUSH_AT + 4096);
    auto flush = [&]() {
        if (!wbuf.empty()) {
            std::fwrite(wbuf.data(), 1, wbuf.size(), out);
            wbuf.clear();
        }
    };
    auto emit = [&](const char *data, size_t len) {
        if (len >= FLUSH_AT) {
            flush();
            std::fwrite(data, 1, len, out);
            return;
        }
        if (wbuf.size() + len > FLUSH_AT) flush();
        wbuf.append(data, len);
    };

    for (size_t i = 0; i < n_super; ++i) {
        if (super_line_strs[i].empty()) continue;  // removed
        emit(super_line_strs[i].data(), super_line_strs[i].size());
    }

    char lbuf[1024];
    for (const auto &[_, edge] : merged_edges) {
        const auto &src_name = supernodes[edge.src_super].name;
        const auto &dst_name = supernodes[edge.dst_super].name;
        int n = std::snprintf(
            lbuf, sizeof(lbuf),
            "L\t%s\t%c\t%s\t%c\t0M\tSO:Z:%c\tTO:Z:%c\tSS:Z:%s\tTS:Z:%s\t"
            "EC:i:%zu\tFN:i:%" PRIu64 "\tTN:i:%" PRIu64 "\n",
            src_name.c_str(), edge.src_ori, dst_name.c_str(), edge.dst_ori,
            edge.src_ori, edge.dst_ori, side_label(edge.src_side),
            side_label(edge.dst_side), edge.count,
            static_cast<uint64_t>(edge.example_from),
            static_cast<uint64_t>(edge.example_to));
        if (n > 0)
            emit(lbuf,
                 std::min<size_t>(static_cast<size_t>(n), sizeof(lbuf) - 1));
    }
    flush();
}

// ── End joins: join overlapping dead segment ends ───────────────────────
// When the + and − strand unitigs of one locus stop at different holes, the
// two layouts end in segments whose free ends overlap as reverse complements:
// the last bases of segment A, read toward its free end, are the reverse
// complement of the last bases of segment B, read toward B's free end. Such a
// pair is one sequence: A, then B read backwards.
//
// run() goes through five steps:
//  1. free_ends: every free end (no link on that side), with its tail (the
//     last kEndWindow bases, read toward the end) and the tail's cores one
//     LCP level below the assembly level (denser cores, so overlaps of a few
//     hundred bases hold several).
//  2. find_overlaps: end pairs whose tails share at least two core pairs on
//     one diagonal, confirmed by extending one shared core pair with the
//     ksw2 aligner (align.h) to both free ends (within kBand) at
//     >= kMinIdentity.
//  3. pick_joins: overlaps whose ends overlap nothing else, whose end
//     unitigs are below the locus depth limit (one locus, not a collapsed
//     repeat) and whose tails are not tandem-periodic (collapsed satellite
//     arrays can read below the depth limit).
//  4. merge_chains: every chain of joined segments becomes one segment.
//  5. rewire_links: links of merged segments move to the merged segment.
struct EndJoin {
    static constexpr size_t kEndWindow = 10000;
    // Positional tolerance in bases: the extension band, the diagonal bin
    // that shared core pairs vote in, and how far short of a free end an
    // extension may stop.
    static constexpr int kBand = 64;
    static constexpr double kMinIdentity = 0.99;
    static constexpr size_t kMaxKeyOcc = 4;
    // Separates substantial within-window recurrence (tandem arrays) from
    // occasional duplicate pairs, which stay well below it.
    static constexpr double kMaxRepeatedPairFraction = 0.30;

    // A free segment end. `tail` is the segment's last bases read toward the
    // end, `tail_rc` its reverse complement (read from the end inward);
    // `tail_cores` and `tail_rc_cores` are their cores.
    struct end_t {
        u64 sid;
        bool at_end;  // the segment's end side (else its start side)
        gtl::vector<SequenceContainerType> tail, tail_rc;  // one window each
        lcp_read tail_cores, tail_rc_cores;
    };

    // Ends a and b overlap when tail(a) ends with the start of tail_rc(b).
    // The join point is where their alignment ends; each end loses the
    // bases beyond it at its free end (trim_a, trim_b).
    struct overlap_t {
        size_t a, b;
        size_t trim_a, trim_b;
    };

    // A segment side: sid * 2 + at_end.
    static u64 side(u64 sid, bool at_end) { return sid * 2 + at_end; }

    // The join at one segment side: the partner side and the bases cut off
    // this side.
    struct join_t {
        u64 partner;
        size_t trim;
    };
    using join_map = gtl::flat_hash_map<u64, join_t>;  // side -> join

    // Where a merged segment went: the merged segment's id, and whether the
    // segment is read reversed in it.
    struct moved_t {
        u64 sid;
        bool reversed;
    };
    using moved_map = gtl::flat_hash_map<u64, moved_t>;  // old sid -> moved

    // Counts for the log line.
    struct stats_t {
        size_t ends = 0, candidates = 0, aligned = 0, deep = 0, tandem = 0,
               joined = 0;
    };

    static u64 pair_key(const lcp_read &lr, size_t c) {
        return (static_cast<u64>(lr.cores[c - 1]) << 32) | lr.cores[c];
    }

    static size_t segment_length(const supernode_out &sn) {
        size_t len = 0;
        for (const auto &s : sn.seqs()) len += s.len;
        return len;
    }

    // Fraction of consecutive core-pair occurrences whose key occurs again
    // in the same window (including the first occurrence).
    static double repeated_pair_fraction(const lcp_read &lr) {
        if (lr.cores.size() < 2) return 0.0;
        gtl::flat_hash_map<u64, size_t> counts;
        for (size_t c = 1; c < lr.cores.size(); ++c) ++counts[pair_key(lr, c)];
        size_t repeated = 0;
        for (const auto &[_, n] : counts)
            if (n > 1) repeated += n;
        return static_cast<double>(repeated) / (lr.cores.size() - 1);
    }

    static bool is_tandem(const lcp_read &window) {
        return repeated_pair_fraction(window) >= kMaxRepeatedPairFraction;
    }

    // Per-thread buffers for confirm, reused across candidates: the unpacked
    // bases right and left of the seed in tail(a) and tail_rc(b), and the
    // extensions' CIGARs.
    struct confirm_scratch {
        aligner dp;
        gtl::vector<uint8_t> ra, rb, la, lb;
        cigar right, left;
    };

    // Step 1: the free ends with their tails and the tails' cores.
    template <class AssemblerT>
    static gtl::vector<end_t> free_ends(const AssemblerT &as,
                                        const supernode_vec &supernodes,
                                        const edge_map &edges, int threads) {
        gtl::flat_hash_set<u64> linked;  // linked segment sides
        for (const auto &[_, e] : edges) {
            linked.insert(side(e.src_super, e.src_ori == '+'));
            linked.insert(side(e.dst_super, e.dst_ori != '+'));
        }
        gtl::vector<end_t> ends;
        for (size_t sid = 0; sid < supernodes.size(); ++sid) {
            const auto &sn = supernodes[sid];
            if (sn.name.empty()) continue;
            const size_t len = segment_length(sn);
            const size_t w = std::min(len, kEndWindow);
            for (const bool at_end : {false, true}) {
                if (linked.contains(side(sid, at_end))) continue;
                // The window as stored in the segment; at the start side it
                // is read away from the end, so the tail is its reverse
                // complement.
                end_t e{sid, at_end, {}, {}, {}, {}};
                if (at_end) {
                    SequenceContainerType fwd(sn.seqs(), len - w, len);
                    e.tail_rc.push_back(fwd.reverse_complement());
                    e.tail.push_back(std::move(fwd));
                } else {
                    SequenceContainerType fwd(sn.seqs(), 0, w);
                    e.tail.push_back(fwd.reverse_complement());
                    e.tail_rc.push_back(std::move(fwd));
                }
                ends.push_back(std::move(e));
            }
        }
        const int level = std::max(1, as.lcp_level - 1);
        parallel_for(ends.size(), threads, [&](size_t i) {
            core_buf str;
            ends[i].tail_cores = AssemblerT::build_lr_cores(
                ends[i].tail, as.syncmer_functor, level, str);
            ends[i].tail_rc_cores = AssemblerT::build_lr_cores(
                ends[i].tail_rc, as.syncmer_functor, level, str);
        });
        return ends;
    }

    // A shared core pair: its end in tail(a) and in tail_rc(b).
    struct seed_t {
        size_t a;
        uint32_t pos_a, pos_b;
    };

    // For every end a > b whose tail shares core pairs with tail_rc(b) (keys
    // occurring at most kMaxKeyOcc times among all tails), a seed on the
    // most-voted diagonal (the position in tail(a) where tail_rc(b) starts,
    // voted in kBand bins), kept when >= 2 core pairs vote for it.
    using tail_index =  // core-pair key -> (end, pair end in its tail)
        gtl::flat_hash_map<u64, gtl::vector<std::pair<uint32_t, uint32_t>>>;
    static gtl::vector<seed_t> candidates_for(size_t b,
                                              const gtl::vector<end_t> &ends,
                                              const tail_index &index) {
        struct bin_t {
            uint32_t votes = 0;
            uint32_t pos_a = 0, pos_b = 0;  // the last pair in the bin
        };
        const lcp_read &query = ends[b].tail_rc_cores;
        gtl::flat_hash_map<std::pair<uint32_t, long>, bin_t> bins;
        for (size_t c = 1; c < query.cores.size(); ++c) {
            auto it = index.find(pair_key(query, c));
            if (it == index.end() || it->second.size() > kMaxKeyOcc) continue;
            const uint32_t pos_b = query.pos[c].second;
            for (const auto &[a, pos_a] : it->second) {
                if (a <= b || ends[a].sid == ends[b].sid) continue;
                const long diagonal = (long)pos_a - (long)pos_b;
                if (diagonal < 0) continue;
                bin_t &bin = bins[{a, diagonal / kBand}];
                ++bin.votes;
                bin.pos_a = pos_a;
                bin.pos_b = pos_b;
            }
        }
        gtl::flat_hash_map<uint32_t, bin_t> best;  // a -> its best bin
        for (const auto &[a_bin, bin] : bins)
            if (bin.votes > best[a_bin.first].votes) best[a_bin.first] = bin;
        gtl::vector<seed_t> out;
        for (const auto &[a, bin] : best)
            if (bin.votes >= 2) out.push_back({a, bin.pos_a, bin.pos_b});
        return out;
    }

    // The overlap of ends a and b through `seed`, if extending it reaches
    // the end of tail(a) and the start of tail_rc(b) (both free ends, within
    // kBand) at >= kMinIdentity.
    static std::optional<overlap_t> confirm(const gtl::vector<end_t> &ends,
                                            size_t b, const seed_t &seed) {
        thread_local confirm_scratch sc;
        const SequenceContainerType &ta = ends[seed.a].tail[0];
        const SequenceContainerType &tb = ends[b].tail_rc[0];
        const size_t pa = seed.pos_a, pb = seed.pos_b;
        // Right of the seed, toward the end of tail(a).
        sc.ra.resize(ta.len - pa);
        ta.unpack_bytes(pa, ta.len, sc.ra.data());
        const size_t rb_to = std::min(tb.len, pb + sc.ra.size() + kBand);
        sc.rb.resize(rb_to - pb);
        tb.unpack_bytes(pb, rb_to, sc.rb.data());
        sc.right.clear();
        int right_b = 0, right_a = 0;
        sc.dp.extend(sc.rb.data(), int(sc.rb.size()), sc.ra.data(),
                     int(sc.ra.size()), kBand, sc.right, right_b, right_a);
        if (pa + right_a + kBand < ta.len) return std::nullopt;
        // Left of the seed on the reversed prefixes, toward the start of
        // tail_rc(b).
        size_t la_from = 0;
        if (pa > pb + kBand) la_from = pa - pb - kBand;
        sc.la.resize(pa - la_from);
        ta.unpack_bytes_reversed(la_from, pa, sc.la.data());
        sc.lb.resize(pb);
        tb.unpack_bytes_reversed(0, pb, sc.lb.data());
        sc.left.clear();
        int left_b = 0, left_a = 0;
        sc.dp.extend(sc.lb.data(), int(sc.lb.size()), sc.la.data(),
                     int(sc.la.size()), kBand, sc.left, left_b, left_a);
        if (size_t(left_b + kBand) < pb) return std::nullopt;
        const align_stats r =
            sc.dp.score_cigar(sc.right, sc.rb.data(), sc.ra.data());
        const align_stats l =
            sc.dp.score_cigar(sc.left, sc.lb.data(), sc.la.data());
        const uint32_t matches = r.matches + l.matches;
        if (matches < kMinIdentity * (matches + r.nm() + l.nm()))
            return std::nullopt;
        // The join point is where the right extension ends: tail(a) keeps its
        // bases up to it, b drops its bases before it.
        return overlap_t{seed.a, b, ta.len - (pa + right_a), pb + right_b};
    }

    // Step 2: every overlap, each end pair checked once (b < a; the other
    // direction is its reverse complement).
    static gtl::vector<overlap_t> find_overlaps(const gtl::vector<end_t> &ends,
                                                int threads, stats_t &stats) {
        tail_index index;
        for (size_t a = 0; a < ends.size(); ++a) {
            const lcp_read &tail = ends[a].tail_cores;
            for (size_t c = 1; c < tail.cores.size(); ++c)
                index[pair_key(tail, c)].emplace_back((uint32_t)a,
                                                      tail.pos[c].second);
        }
        gtl::vector<gtl::vector<overlap_t>> found(ends.size());  // by b
        std::atomic<size_t> candidates{0};
        parallel_for(ends.size(), threads, [&](size_t b) {
            for (const seed_t &seed : candidates_for(b, ends, index)) {
                ++candidates;
                if (auto overlap = confirm(ends, b, seed))
                    found[b].push_back(*overlap);
            }
        });
        gtl::vector<overlap_t> overlaps;
        for (const auto &per_end : found)
            overlaps.insert(overlaps.end(), per_end.begin(), per_end.end());
        stats.candidates = candidates;
        stats.aligned = overlaps.size();
        return overlaps;
    }

    // Step 3: the overlaps kept as joins, recorded at both of their sides.
    template <class EndDepth>
    static join_map pick_joins(const gtl::vector<end_t> &ends,
                               const gtl::vector<overlap_t> &overlaps,
                               double locus_limit, const EndDepth &end_depth,
                               stats_t &stats) {
        gtl::vector<size_t> n_overlaps(ends.size(), 0);
        for (const overlap_t &o : overlaps) {
            ++n_overlaps[o.a];
            ++n_overlaps[o.b];
        }
        auto deep = [&](size_t e) {
            return end_depth(ends[e].sid, ends[e].at_end) >= locus_limit;
        };
        auto tandem = [&](size_t e) { return is_tandem(ends[e].tail_cores); };
        auto side_of = [&](size_t e) {
            return side(ends[e].sid, ends[e].at_end);
        };
        join_map joins;
        for (const overlap_t &o : overlaps) {
            if (n_overlaps[o.a] != 1 || n_overlaps[o.b] != 1) continue;
            if (deep(o.a) || deep(o.b)) {
                ++stats.deep;
                continue;
            }
            if (tandem(o.a) || tandem(o.b)) {
                ++stats.tandem;
                continue;
            }
            joins[side_of(o.a)] = {side_of(o.b), o.trim_a};
            joins[side_of(o.b)] = {side_of(o.a), o.trim_b};
        }
        return joins;
    }

    // A chain member: the segment, whether it is read reversed, and the
    // bases cut off its first and last side in that direction.
    struct member_t {
        u64 sid;
        bool reversed;
        size_t trim_first, trim_last;
    };

    // The chain from segment `first`, whose free side is the chain's outer
    // end: read forward when that is its start, reversed when its end.
    static gtl::vector<member_t> walk_chain(const join_map &joins, u64 first,
                                            bool reversed) {
        gtl::vector<member_t> members;
        size_t trim_first = 0;
        for (u64 sid = first;;) {
            const auto out = joins.find(side(sid, !reversed));  // side we leave
            if (out == joins.end()) {
                members.push_back({sid, reversed, trim_first, 0});
                return members;
            }
            members.push_back({sid, reversed, trim_first, out->second.trim});
            const u64 next = out->second.partner;
            trim_first = joins.at(next).trim;
            sid = next / 2;
            reversed = next % 2;  // entered at its end: read reversed
        }
    }

    // One segment spelling the trimmed members in order, named J<new_sid>.
    // Its layout tags come from its first layout member (else its first
    // member).
    static supernode_out merged_segment(const supernode_vec &supernodes,
                                        const gtl::vector<member_t> &members,
                                        size_t new_sid) {
        const supernode_out &tag_src = [&]() -> const supernode_out & {
            for (const member_t &m : members)
                if (supernodes[m.sid].is_layout) return supernodes[m.sid];
            return supernodes[members[0].sid];
        }();
        supernode_out merged;
        merged.name = "J" + std::to_string(new_sid);
        merged.is_layout = tag_src.is_layout;
        merged.aux_id = tag_src.aux_id;
        for (const member_t &m : members) {
            const auto &src = supernodes[m.sid];
            const size_t len = segment_length(src);
            if (m.trim_first + m.trim_last >= len) continue;
            if (m.reversed)  // the kept range, cut first, then reversed
                merged.owned_seqs.push_back(
                    SequenceContainerType(src.seqs(), m.trim_last,
                                          len - m.trim_first)
                        .reverse_complement());
            else
                merged.owned_seqs.emplace_back(src.seqs(), m.trim_first,
                                               len - m.trim_last);
        }
        return merged;
    }

    // Step 4: every chain of joined segments becomes one new segment and its
    // members are emptied. A chain starts at a segment joined on one side
    // only; a cycle has no such segment and is left as it is. Returns where
    // each member went.
    static moved_map merge_chains(supernode_vec &supernodes,
                                  const join_map &joins, stats_t &stats) {
        moved_map moved;
        const size_t n_orig = supernodes.size();
        for (u64 sid = 0; sid < n_orig; ++sid) {
            if (supernodes[sid].name.empty()) continue;  // merged already
            const bool start_joined = joins.contains(side(sid, false));
            const bool end_joined = joins.contains(side(sid, true));
            if (start_joined == end_joined) continue;  // no chain end
            const auto members = walk_chain(joins, sid, start_joined);
            const size_t new_sid = supernodes.size();
            supernodes.push_back(merged_segment(supernodes, members, new_sid));
            for (const member_t &m : members) {
                moved[m.sid] = {new_sid, m.reversed};
                supernodes[m.sid].name.clear();
                gtl::vector<SequenceContainerType>().swap(
                    supernodes[m.sid].owned_seqs);
            }
            stats.joined += members.size() - 1;
        }
        return moved;
    }

    // Step 5: links of merged segments move to the merged segment, with
    // orientation and side flipped for reversed members; links that become
    // identical are combined.
    static void rewire_links(const moved_map &moved, edge_map &edges) {
        auto move_end = [&](u64 &super, char &ori, uint8_t &side) {
            auto m = moved.find(super);
            if (m == moved.end()) return;
            super = m->second.sid;
            if (!m->second.reversed) return;
            if (ori == '+')
                ori = '-';
            else
                ori = '+';
            if (side == SIDE_START)
                side = SIDE_END;
            else if (side == SIDE_END)
                side = SIDE_START;
        };
        edge_map rewired;
        for (const auto &[_, link] : edges) {
            edge_agg e = link;
            move_end(e.src_super, e.src_ori, e.src_side);
            move_end(e.dst_super, e.dst_ori, e.dst_side);
            const edge_key key{e.src_super, e.dst_super, e.src_ori,
                               e.dst_ori,   e.src_side,  e.dst_side};
            if (auto [it, inserted] = rewired.try_emplace(key, e); !inserted)
                it->second.count += e.count;
        }
        edges.swap(rewired);
    }

    // end_depth(sid, at_end): depth of the unitig at that segment end.
    template <class AssemblerT, class EndDepth>
    static size_t run(const AssemblerT &as, supernode_vec &supernodes,
                      edge_map &edges, double locus_limit,
                      const EndDepth &end_depth) {
        if (locus_limit <= 0.0) {
            std::fprintf(stderr,
                         "[end-join] skipped: no locus depth limit (no depth "
                         "model)\n");
            return 0;
        }
        const int threads = std::max(1, as.threads.work);
        stats_t stats;
        const gtl::vector<end_t> ends =
            free_ends(as, supernodes, edges, threads);
        stats.ends = ends.size();
        const auto overlaps = find_overlaps(ends, threads, stats);
        const join_map joins =
            pick_joins(ends, overlaps, locus_limit, end_depth, stats);
        const moved_map moved = merge_chains(supernodes, joins, stats);
        if (!moved.empty()) rewire_links(moved, edges);
        std::fprintf(stderr,
                     "[end-join] free ends %zu, candidates %zu, aligned %zu, "
                     "rejected deep %zu, rejected tandem %zu, joins %zu\n",
                     stats.ends, stats.candidates, stats.aligned, stats.deep,
                     stats.tandem, stats.joined);
        return stats.joined;
    }
};

// Depth of the unitig at a segment end: a C segment's own unitig, or the
// layout member reaching furthest towards that end (0 when unknown). A layout
// whose cores split into several components is emitted as several segments
// (M<idx>, M<idx>_1, ...) but its members are not recorded per component, so
// only the first segment of a layout gets a depth; the others report
// infinity and their ends are not joined.
class segment_end_depth {
   public:
    segment_end_depth(const supernode_vec &supernodes, const dbg_index &dbg,
                      const gtl::vector<ultra_layout_job> &layout_jobs)
        : supernodes_(supernodes), dbg_(dbg) {
        for (const auto &job : layout_jobs) job_of_[job.layout_idx] = &job;
        for (u64 sid = 0; sid < supernodes.size(); ++sid)
            if (supernodes[sid].is_layout)
                first_segment_of_.try_emplace(supernodes[sid].aux_id, sid);
    }

    double operator()(u64 sid, bool at_end) const {
        const auto &sn = supernodes_[sid];
        if (!sn.is_layout) return node_depth(sn.aux_id);
        if (first_segment_of_.at(sn.aux_id) != sid)
            return std::numeric_limits<double>::infinity();
        auto j = job_of_.find(sn.aux_id);
        if (j == job_of_.end() || j->second->lnodes.empty()) return 0.0;
        const auto &ln = j->second->lnodes;
        if (at_end) {
            const auto furthest = std::max_element(
                ln.begin(), ln.end(),
                [](const auto &x, const auto &y) { return x.to < y.to; });
            return node_depth(furthest->id);
        }
        const auto furthest = std::min_element(
            ln.begin(), ln.end(),
            [](const auto &x, const auto &y) { return x.from < y.from; });
        return node_depth(furthest->id);
    }

   private:
    double node_depth(u64 id) const {
        auto it = dbg_.dbg.find(id);
        if (it == dbg_.dbg.end() || it->second.redirected() ||
            it->second.seq.empty())
            return 0.0;
        return it->second.avg_core_count();
    }

    const supernode_vec &supernodes_;
    const dbg_index &dbg_;
    gtl::flat_hash_map<u64, const ultra_layout_job *> job_of_;
    gtl::flat_hash_map<u64, u64> first_segment_of_;  // layout -> segment
};

// Sequences of the output segments (supernodes still named after post-layout
// cleaning and end joins), each concatenated into one owned container: the
// segments a --multipass pass seeds into the next pass.
inline gtl::vector<SequenceContainerType> segment_sequences(
    const supernode_vec &supernodes, int nthreads) {
    gtl::vector<const supernode_out *> kept;
    for (const auto &sn : supernodes) {
        if (!sn.name.empty()) kept.push_back(&sn);
    }
    gtl::vector<SequenceContainerType> out(kept.size());
    std::atomic<size_t> next{0};
    gtl::vector<std::thread> workers;
    for (int t = 0; t < std::max(1, nthreads); ++t) {
        workers.emplace_back([&]() {
            for (size_t i = next++; i < kept.size(); i = next++) {
                const auto &seqs = kept[i]->seqs();
                out[i] = SequenceContainerType(seqs.begin(), seqs.end());
            }
        });
    }
    for (auto &w : workers) w.join();
    return out;
}

// Layout and emission of the cleaned DBG: fold, layout supernodes, edge
// aggregation, post-layout cleaning and end joins, then
// `sink(supernodes, edges)` (`sink_name` labels its timing line).
template <class AssemblerT, class Sink>
inline void build_merged_graph(AssemblerT &as, dbg_index &dbg,
                               const char *sink_name, Sink &&sink) {
    INIT_TIMERS(LayoutEmit);
    auto [recig, layouts, dropped_crumbs] = as.fold_dbg_graph(dbg);
    TIME_CHECKPOINT(LayoutEmit, stderr, "fold_dbg_graph (%zu layouts)",
                    layouts.size());

    gtl::vector<supernode_out> supernodes;
    gtl::flat_hash_map<u64, u64> node_to_super;
    gtl::flat_hash_set<u64> nodes_in_layouts;
    gtl::flat_hash_map<u64, uint8_t> node_boundary_mask;
    node_to_super.reserve(dbg.active_set.size() + 8);
    nodes_in_layouts.reserve(dbg.active_set.size() + 8);

    // Resolve per-node orientation up front (see resolve_layout_orientations).
    // A nonzero conflict count means the single-owner invariant was violated.
    size_t ori_conflicts = 0;
    auto node_layout_ori = resolve_layout_orientations(
        layouts, orientation_policy::first_wins, &ori_conflicts);
    if (ori_conflicts) [[unlikely]]
        std::fprintf(stderr,
                     "[layout-emit] WARNING: %zu node(s) claimed by multiple "
                     "layouts with conflicting orientation (policy=first_wins)\n",
                     ori_conflicts);

    // ── Phase 1: one job per layout component + per-node boundary masks ──
    gtl::vector<ultra_layout_job> layout_jobs;
    layout_jobs.reserve(layouts.size());
    build_layout_jobs(as, layouts, layout_jobs, nodes_in_layouts,
                      node_boundary_mask);

    // ── Phase 2: build layout ("M") supernode sequences in parallel ─────
    auto layout_seqs = build_layout_sequences_parallel(as, dbg, layout_jobs);
    TIME_CHECKPOINT(LayoutEmit, stderr, "Build layout sequences");
    for (size_t i = 0; i < layout_jobs.size(); ++i) {
        const auto &job = layout_jobs[i];
        auto &comps = layout_seqs[i];
        if (comps.empty()) continue;
        // One contig per SCD component. The first keeps the plain "M<idx>"
        // name; any further component gets an "_<c>" suffix for a unique GFA
        // segment name.
        const u64 primary_sid = supernodes.size();
        for (size_t c = 0; c < comps.size(); ++c) {
            std::string name = "M" + std::to_string(job.layout_idx);
            if (c) name += "_" + std::to_string(c);
            supernodes.push_back({std::move(name), std::move(comps[c]), nullptr,
                                  true, job.layout_idx});
        }
        // Edges are aggregated per node; map every node of the job to the
        // primary component. Secondary components carry no claimed adjacency.
        for (const auto &n : job.lnodes) {
            node_to_super[n.id] = primary_sid;
        }
    }

    // ── Phase 3: drop C-nodes already contained in layouts; emit the rest
    //             as unitig ("C") supernodes ─────────────────────────────
    gtl::flat_hash_set<u64> candidate_set;
    auto contained_nodes =
        find_contained_cnodes(dbg, nodes_in_layouts, node_boundary_mask,
                              node_layout_ori, candidate_set);
    for (u64 node_id : candidate_set) {
        if (contained_nodes.contains(node_id)) continue;
        if (dropped_crumbs.contains(node_id)) continue;
        auto it = dbg.dbg.find(node_id);
        if (it == dbg.dbg.end()) continue;
        const u64 sid = supernodes.size();
        supernodes.push_back(
            {"C" + std::to_string(node_id), {}, &it->second.seq, false, node_id});
        node_to_super[node_id] = sid;
    }
    TIME_CHECKPOINT(LayoutEmit, stderr, "Supernode creation");

    // ── Phase 4: aggregate inter-supernode edges ────────────────────────
    edge_map merged_edges =
        aggregate_superedges(dbg, node_to_super, nodes_in_layouts,
                             node_layout_ori, node_boundary_mask,
                             as.threads.work);
    TIME_CHECKPOINT(LayoutEmit, stderr, "Edge aggregation");

    // ── Phase 5: post-layout graph cleaning ─────────────────────────────
    PostLayoutCleaning::run(supernodes, merged_edges);
    TIME_CHECKPOINT(LayoutEmit, stderr, "Post-layout cleaning");

    // ── Phase 5b: join overlapping dead ends of one locus ───────────────
    if (as.end_joins) {
        EndJoin::run(as, supernodes, merged_edges, dbg.locus_limit,
                     segment_end_depth(supernodes, dbg, layout_jobs));
        TIME_CHECKPOINT(LayoutEmit, stderr, "End joins");
    }

    // ── Phase 6: hand the graph to the sink (GFA writing) ───────────────
    std::forward<Sink>(sink)(std::as_const(supernodes),
                             std::as_const(merged_edges));
    TIME_CHECKPOINT(LayoutEmit, stderr, "%s", sink_name);
}

template <class AssemblerT>
inline void emit_ultra_merged_graph(AssemblerT &as, dbg_index &dbg, FILE *out) {
    build_merged_graph(
        as, dbg, "GFA writing",
        [&](const supernode_vec &supernodes, const edge_map &merged_edges) {
            write_gfa(supernodes, merged_edges, as.threads.work, out);
        });
}

template <class AssemblerT>
inline void emit_ultra_merged_graph(AssemblerT &as, dbg_index &dbg) {
    emit_ultra_merged_graph(as, dbg, stdout);
}

}  // namespace layout_graph_emit

#endif  // _LAYOUT_GRAPH_H
