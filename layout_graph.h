#ifndef _LAYOUT_GRAPH_H
#define _LAYOUT_GRAPH_H

#include <sys/stat.h>

#include <algorithm>
#include <atomic>
#include <chrono>
#include <cinttypes>
#include <cstdio>
#include <cstring>
#include <expected>
#include <gtl/phmap_fwd_decl.hpp>
#include <gtl/vector.hpp>
#include <map>
#include <memory>
#include <optional>
#include <queue>
// #include <stack>
#include <string>
#include <thread>
#include <tuple>
#include <utility>

#include "graph.h"

namespace layout_graph_algo {

using layout_result_t =
    gtl::vector<gtl::flat_hash_map<u64, std::tuple<i64, i64, bool>>>;

template <class RecigMapT>
inline layout_result_t connected_comp_plus_layout(const RecigMapT &recig) {
    layout_result_t layouts;
    std::queue<std::pair<u64, typename RecigMapT::mapped_type>> q;
    gtl::flat_hash_set<u64> visited;

    // Precompute node → length once (was an O(N) scan of `recig` per node).
    // Mirrors the original precedence: a node's length is the `query.len` of the
    // first overlap keyed by that node (when > 0); otherwise the `target_span.len`
    // of the first overlap (in iteration order) that references it as a target.
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
        // (CoreSpan.len ≥ 1). fold_dbg_graph only admits nodes with
        // cores.size() ≥ olc_lcmer_l+1, so this branch is unreachable on real
        // input (C11-verified: zero hits on chr22). Warn rather than silently
        // emitting a zero-length node if the invariant is ever violated.
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

// ── Tunable constants (were scattered magic numbers) ────────────────
// Reserve hint: expected seed-chain-DAG nodes (cores) per contig.
inline constexpr size_t SCD_CORES_PER_CONTIG_HINT = 64;
// Default max length (bp) of a post-layout C-node tip to prune.
inline constexpr size_t DEFAULT_POST_LAYOUT_MAX_TIP_LEN = 50000;
// Layout-boundary fuzz: a node counts as touching an end if it is within
// max(olc_lcmer_l * BOUNDARY_FUZZ_LCMER_MULT, span / BOUNDARY_FUZZ_SPAN_DIV).
inline constexpr i64 BOUNDARY_FUZZ_LCMER_MULT = 2;
inline constexpr i64 BOUNDARY_FUZZ_SPAN_DIV = 20;
// (Retired: N_GAP_LEN. Disconnected SCD components used to be concatenated with
// a fabricated 64-'N' gap — root cause R7, docs/layout-graph-improvement-plan.md
// — asserting an adjacency the data does not support. They are now emitted as
// separate contigs; a real join needs P3.3's oriented, validated DBG-edge path.)

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

// ── Phase 1: Build seed-chain DAG from contig cores ─────────────────
inline void build_scd_graph(
    const gtl::vector<size_t> &contig_order,
    const gtl::flat_hash_map<size_t, lcp_read> &contig_lr,
    const gtl::vector<lane_node> &nodes,
    gtl::vector<scd_node> &scd,
    gtl::vector<scd_adj_t> &scd_fwd,
    gtl::flat_hash_map<uint32_t, gtl::vector<uint32_t>> &label_nodes,
    gtl::vector<gtl::vector<std::pair<size_t, size_t>>> &node_contigs)
{
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
        for (size_t k = 0; k < M; ++k) {
            auto it = label_nodes.find(lr.cores[k]);
            if (it == label_nodes.end()) continue;
            for (uint32_t v : it->second) {
                if (v >= gn) continue;
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
}

// DISABLED EXPERIMENT — not compiled into the default build. Enable with
// -DLASM_LAYOUT_RESOLVE_OVERLAP to trial it on a full genome.
//
// Drops any SCD component overlapping a heavier one by >50% of the shorter
// interval, keeping the heavier. It does clear the M79 duplication, but it is
// NOT trustworthy as a general rule:
//   * the 50% threshold is a magic constant with no principled basis;
//   * it acts on layout_from/layout_to (matching-phase coordinates), which are
//     a proxy for duplicated sequence, not proof — on other data it could drop
//     genuinely distinct sequence;
//   * it was derived from, and "validated" on, a single occurrence (layout 79
//     was the only multi-component layout in the whole l5 chr22 run).
// The underlying defect is R7 (docs/layout-graph-improvement-plan.md:37): the
// layout fabricated adjacency between disconnected components. That fabrication
// is removed at the emission site instead; see collect_path_sequences.
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
// traceback per component; no post-hoc span extension.
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

        // Traceback the heaviest path ending at best_end. Because dp[] is the
        // longest-weight path from a virtual source wired to every in-degree-0
        // node, and best_end is the heaviest out-degree-0 node (a virtual sink
        // wired from every real sink), this single traceback IS the heaviest
        // virtual-source→virtual-sink path for the component. No post-hoc span
        // extension exists (the old forward/backward heuristics — C13's bug
        // class — are gone); the path is edge-connected by construction.
        gtl::vector<uint32_t> path;
        for (int32_t cur = best_end; cur >= 0; cur = parent[cur])
            path.push_back(static_cast<uint32_t>(cur));
        std::reverse(path.begin(), path.end());

        // Invariant I3: consecutive path nodes share an SCD edge. parent[] only
        // ever points at a real predecessor, so this holds by construction —
        // verify defensively since NDEBUG makes assert() a no-op on release.
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

    // Per-layout component diagnostic (R3/R4 root-cause work): a layout whose
    // nodes fail to chain into ONE seed-chain path is the condition that used to
    // trigger the fabricated N-gap join. Overlapping spans additionally mean the
    // same region was reconstructed twice.
    if (comp_paths.size() > 1 && getenv("LASM_DEBUG_COMPSPANS")) {
        bool any_overlap = false;
        for (size_t i = 0; i + 1 < comp_paths.size() && !any_overlap; ++i)
            for (size_t j = i + 1; j < comp_paths.size(); ++j)
                if (std::min(comp_paths[i].layout_to, comp_paths[j].layout_to) >
                    std::max(comp_paths[i].layout_from,
                             comp_paths[j].layout_from)) {
                    any_overlap = true;
                    break;
                }
        for (size_t i = 0; i < comp_paths.size(); ++i)
            std::fprintf(stderr,
                         "[layout-dag] layout=%" PRIu64 " ncomp=%zu comp=%zu "
                         "from=%" PRId64 " to=%" PRId64 " weight=%" PRId64
                         " overlapping=%s\n",
                         static_cast<uint64_t>(layout_idx), comp_paths.size(), i,
                         static_cast<int64_t>(comp_paths[i].layout_from),
                         static_cast<int64_t>(comp_paths[i].layout_to),
                         static_cast<int64_t>(comp_paths[i].total_weight),
                         any_overlap ? "yes" : "no");
    }
    return comp_paths;
}

// ── Phase 3: Collect output sequences from heaviest paths ───────────
// Returns ONE sequence group per SCD component. Components are weakly
// disconnected in the seed-chain DAG, so there is no evidence for their
// adjacency, order, or separation: they must not be concatenated. (Previously
// they were joined with a fabricated 64-N gap — root cause R7,
// docs/layout-graph-improvement-plan.md:37 — which asserted an adjacency the
// data does not support and produced the M79 tandem duplication.) The caller
// emits each group as its own contig; a real join needs the oriented, validated
// DBG-edge approach of P3.3, not a synthetic gap.
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
                // path nodes stay on the same contig (real read adjacency);
                // otherwise pick the longest available core sequence — the most
                // complete representative — instead of an arbitrary first contig.
                // (No fallback logging: it also fires legitimately at every path
                // start, so it would be noise rather than signal.)
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
// disconnected components are never concatenated.
inline gtl::vector<gtl::vector<SequenceContainerType>>
build_layout_sequence_from_nodes(
    [[maybe_unused]] AssemblerT &as, dbg_index &dbg,
    const gtl::vector<ultra_layout_node_info> &lnodes,
    [[maybe_unused]] int hash_combine_len, u64 layout_idx = 0) {
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
    lps lps_buf = {0, 0, 0, 0};

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
            contig_core_seqs[ci] = std::move(synth_vec);
            contig_order.push_back(ci);
            continue;
        }
        contig_lr[ci] = std::move(lr);
        contig_core_seqs[ci] = std::move(cseqs);
        contig_order.push_back(ci);
    }

    if (contig_order.empty()) {
        std::fprintf(stderr, "[layout-dag] layout=%" PRIu64
                     " SKIPPED: no contigs with cores\n",
                     static_cast<uint64_t>(layout_idx));
        free_lps(&lps_buf);
        return {};
    }

    // ── Phase 1: Build seed-chain DAG ───────────────────────────────
    gtl::vector<scd_node> scd;
    gtl::vector<scd_adj_t> scd_fwd;
    gtl::flat_hash_map<uint32_t, gtl::vector<uint32_t>> label_nodes;
    gtl::vector<gtl::vector<std::pair<size_t, size_t>>> node_contigs;
    build_scd_graph(contig_order, contig_lr, nodes,
                    scd, scd_fwd, label_nodes, node_contigs);

    if (scd.empty()) {
        std::fprintf(stderr, "[layout-dag] layout=%" PRIu64
                     " SKIPPED: SCD graph empty\n",
                     static_cast<uint64_t>(layout_idx));
        free_lps(&lps_buf);
        return {};
    }

    // ── Phase 2: Find heaviest paths ────────────────────────────────
    auto comp_paths = find_heaviest_paths(scd, scd_fwd, nodes, layout_idx);

    // ── Phase 3: Collect output sequences ───────────────────────────
    auto result = collect_path_sequences(
        comp_paths, node_contigs, contig_core_seqs, scd.size());

    free_lps(&lps_buf);
    return result;
}

// One SCD component's sequences …
using layout_comp_seqs_t = gtl::vector<SequenceContainerType>;
// … and all components of one layout (emitted as separate contigs).
using layout_seq_t = gtl::vector<layout_comp_seqs_t>;

template <class AssemblerT>
inline gtl::vector<layout_seq_t> build_layout_sequences_parallel(
    AssemblerT &as, dbg_index &dbg, const gtl::vector<ultra_layout_job> &jobs,
    int hash_combine_len) {
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
    for (int tid = 0; tid < nthreads; ++tid) {
        workers.emplace_back([&]() {
            size_t job_idx = 0;
            while (true) {
                if (!job_queue.wait_dequeue_timed(job_idx,
                                                  WAIT_QUEUE_TIMEOUT)) {
                    if (done.load()) break;
                    continue;
                }
                out_seqs[job_idx] = build_layout_sequence_from_nodes(
                    as, dbg, jobs[job_idx].lnodes, hash_combine_len,
                    jobs[job_idx].layout_idx);
            }
            while (job_queue.try_dequeue(job_idx)) {
                out_seqs[job_idx] = build_layout_sequence_from_nodes(
                    as, dbg, jobs[job_idx].lnodes, hash_combine_len,
                    jobs[job_idx].layout_idx);
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

enum class bridge_realization : uint8_t {
    gap_safe,
    overlap_unresolved,
    gap_uncertain,
};

inline const char *bridge_realization_label(bridge_realization status) {
    switch (status) {
        case bridge_realization::gap_safe: return "gap_safe";
        case bridge_realization::overlap_unresolved:
            return "overlap_unresolved";
        case bridge_realization::gap_uncertain: return "gap_uncertain";
    }
    return "gap_uncertain";
}

struct signed_gap_summary {
    uint32_t count = 0;
    int64_t median = 0;
    int64_t p10 = 0;
    int64_t p90 = 0;
    int64_t mad = 0;
};

// One read's signed anchor-run separation. Keeping qname attached until final
// superedge projection prevents one read from being counted again when several
// confirmed anchor pairs collapse onto the same superedge.
struct bridge_read_gap_observation {
    u64 qname = 0;  // qname hash (read identity)
    int64_t gap = 0;
};

// Sort `values` by `less`, then collapse each run of `same`-equal records to the
// run's lower empirical median element (index floor((n-1)/2)). Shared by the
// distinct-vote dedup and the per-qname read-gap collapse; both put `.gap` last
// in `less`, so the retained element carries each group's lower median gap.
template <class T, class Less, class Same>
inline gtl::vector<T> collapse_runs_lower_median(gtl::vector<T> values,
                                                 Less less, Same same) {
    std::sort(values.begin(), values.end(), less);
    gtl::vector<T> out;
    out.reserve(values.size());
    for (size_t first = 0; first < values.size();) {
        size_t last = first + 1;
        while (last < values.size() && same(values[last], values[first])) ++last;
        out.push_back(values[first + (last - first - 1) / 2]);
        first = last;
    }
    return out;
}

// Collapse all observations for each qname to its deterministic lower empirical
// median. The returned gaps therefore have exactly one entry per distinct read.
inline gtl::vector<int64_t> collapse_bridge_read_gaps_by_qname(
    gtl::vector<bridge_read_gap_observation> observations) {
    const auto collapsed = collapse_runs_lower_median(
        std::move(observations),
        [](const auto &a, const auto &b) {
            return std::tie(a.qname, a.gap) < std::tie(b.qname, b.gap);
        },
        [](const auto &a, const auto &b) { return a.qname == b.qname; });
    gtl::vector<int64_t> gaps;
    gaps.reserve(collapsed.size());
    for (const auto &o : collapsed) gaps.push_back(o.gap);
    return gaps;
}

// Sorts observations in place. Quantiles use the deterministic lower empirical
// index floor((n - 1) * p), and MAD uses the same median convention.
inline signed_gap_summary summarize_signed_gaps(gtl::vector<int64_t> &gaps) {
    signed_gap_summary result;
    if (gaps.empty()) return result;
    std::sort(gaps.begin(), gaps.end());
    const size_t n = gaps.size();
    auto quantile = [&](size_t numerator) {
        return gaps[((n - 1) * numerator) / 100];
    };
    result.count = static_cast<uint32_t>(n);
    result.median = quantile(50);
    result.p10 = quantile(10);
    result.p90 = quantile(90);
    gtl::vector<int64_t> deviations;
    deviations.reserve(n);
    for (int64_t gap : gaps)
        deviations.push_back(gap >= result.median ? gap - result.median
                                                  : result.median - gap);
    std::sort(deviations.begin(), deviations.end());
    result.mad = deviations[(n - 1) / 2];
    return result;
}

inline bridge_realization classify_bridge_realization(
    const signed_gap_summary &gap) {
    if (gap.median <= 0) return bridge_realization::overlap_unresolved;
    // MAD can remain tiny when a near-zero mode hides a substantial alternate
    // gap population. Keep such conflicting observations graph-only.
    if (gap.median <= 500 && gap.p90 - gap.median > 1000)
        return bridge_realization::gap_uncertain;
    if (gap.count >= 5 &&
        (gap.mad <= 500 || gap.p90 - gap.p10 <= 1000))
        return bridge_realization::gap_safe;
    return bridge_realization::gap_uncertain;
}

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
    // --thread-bridge: set on read-threaded scaffold edges injected after
    // layout cleaning. Tagged distinctly in the GFA so they are separable from
    // layout junction edges. tb_support = confirming read count; tb_gap =
    // signed anchor separation in read bases.
    bool thread_bridge = false;
    uint32_t tb_support = 0;
    int64_t tb_gap = 0;
    signed_gap_summary tb_gap_summary;
    bridge_realization tb_realization = bridge_realization::gap_uncertain;
    bool tb_c_orientation_unresolved = false;
};

// One cleaned-unitig representative traversed inside a bridge's collapsed tangle,
// with the raw DBG nodes at which the representative supporting read entered and
// exited it. Tier-2 gap materialization reconstructs the gap sequence by walking
// these reps through the layout's overlap-aware assembler; entry/exit fix the
// sub-span and orientation the read actually traversed.
struct tangle_segment {
    u64 rep = 0;
    u64 entry_raw = 0;
    u64 exit_raw = 0;
};

// A confirmed read-threaded bridge: DBG unitig u abuts v (directed u->v), backed
// by `support` reads, separated by an estimated `gap` bases of collapsed tangle.
// Produced by run_thread_bridge; consumed at emit to add a scaffold edge between
// the supernodes carrying u and v.
struct bridge_edge {
    u64 u = 0;
    u64 v = 0;
    uint32_t support = 0;
    int64_t gap = 0;
    signed_gap_summary gap_summary;
    bridge_realization realization = bridge_realization::gap_uncertain;
    gtl::vector<bridge_read_gap_observation> read_gap_observations;
    // Representative tangle path (empty unless Tier-2 capture found one). The
    // representative walk is the supporting read whose tangle length is closest
    // to `gap`; ties broken deterministically for -t1 reproducibility.
    gtl::vector<tangle_segment> gap_path;

    bridge_edge() = default;
    bridge_edge(u64 from, u64 to, uint32_t read_support, int64_t signed_gap = 0)
        : u{from}, v{to}, support{read_support}, gap{signed_gap} {}
    bridge_edge(u64 from, u64 to, uint32_t read_support, int64_t signed_gap,
                signed_gap_summary summary, bridge_realization status,
                gtl::vector<bridge_read_gap_observation> observations)
        : u{from},
          v{to},
          support{read_support},
          gap{signed_gap},
          gap_summary{summary},
          realization{status},
          read_gap_observations{std::move(observations)} {}
};

inline bool thread_bridge_is_fasta_safe(const edge_agg &edge) {
    return edge.thread_bridge &&
           edge.tb_realization == bridge_realization::gap_safe &&
           edge.tb_gap > 0 && !edge.tb_c_orientation_unresolved;
}

using edge_key = std::tuple<u64, u64, char, char, uint8_t, uint8_t>;
using supernode_vec = gtl::vector<supernode_out>;
using edge_map = std::map<edge_key, edge_agg>;

// ── Post-layout graph cleaning ──────────────────────────────────────
// Mirrors GraphCleaning in graph.h: separate struct with individual
// cleaning functions operating on the supernode/superedge graph.
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

    // Remove short C-node tips: C-nodes with no successors (or no predecessors)
    // and sequence length below max_tip_len.
    //
    // Each pass decides all removals from a *consistent snapshot* of `adj`
    // (matching the previous full-rebuild-between-passes behavior), then applies
    // them while patching neighbor adjacency incrementally — avoiding an
    // O(edges) rebuild per pass (measured serial hotspot).
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
// Each layout assigns every DBG node it contains an orientation (+/-). Under the
// current pipeline a node belongs to exactly one layout (connected_comp_plus_
// layout uses a single shared `visited` set), so there is normally nothing to
// resolve. A conflict — the same node claimed by two layouts with opposite
// orientation — can only arise if that single-owner invariant is relaxed (e.g.
// future repeat handling that shares a node across layouts). This function makes
// the resolution explicit and deterministic instead of the old arbitrary
// last-write-wins, and `*out_conflicts` surfaces any violation of the invariant.
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
        const i64 fuzzy = std::max<i64>(as.olc_lcmer_l * BOUNDARY_FUZZ_LCMER_MULT,
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
// boundary or an external non-layout node — i.e. their sequence is already
// represented inside a layout, so they can be dropped.
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

        // A component is contained only if it HAS external adjacency and none
        // of it reaches a layout boundary — i.e. it hangs off layout
        // interiors, so its sequence is already represented inside a layout.
        // A component with NO external adjacency is standalone (isolated
        // node, or a closed circular component whose neighbors are all
        // internal): nothing else represents its sequence, so it must be
        // kept. Previously the vacuous "no boundary neighbor" case dropped
        // such components (and e.g. plasmids) from the assembly.
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

// ── Phase: inject read-threaded bridge edges (--thread-bridge) ──────────
// After layout cleaning, add a scaffold superedge for each confirmed bridge
// (u->v) whose anchors resolve to distinct, live supernodes on real boundary
// sides. Reuses the boundary/orientation gate of aggregate_superedges. A strict
// 1:1 filter (each known physical (supernode, side) endpoint used by exactly
// one bridge) keeps M-only additions simple-path only. Standalone C endpoints
// have no resolved side, remain diagnostic GFA links, and are never used for
// FASTA scaffolding. No sequence work is performed.
struct bridge_injection_result {
    size_t added = 0;
    size_t invalid_endpoint = 0;
    size_t intra_supernode = 0;
    size_t aggregated_records = 0;
    size_t existing_edge = 0;
    size_t projection_conflict = 0;
    size_t endpoint_reuse = 0;
};

inline bridge_injection_result inject_thread_bridge_edges(
    const gtl::vector<bridge_edge> &bridges, const supernode_vec &supernodes,
    const gtl::flat_hash_map<u64, u64> &node_to_super,
    const gtl::flat_hash_set<u64> &nodes_in_layouts,
    const gtl::flat_hash_map<u64, bool> &node_layout_ori,
    const gtl::flat_hash_map<u64, uint8_t> &node_boundary_mask,
    edge_map &merged_edges) {
    bridge_injection_result result;
    if (bridges.empty()) return result;

    // Resolve an anchor to (supernode, orientation char, side). is_src selects
    // the exit side (u) vs the entry side (v), matching aggregate_superedges.
    struct resolved {
        bool ok = false;
        bool orientation_known = false;
        u64 super = 0;
        char ori = '+';
        uint8_t side = 0;
    };
    auto resolve = [&](u64 id, bool is_src) -> resolved {
        auto it = node_to_super.find(id);
        if (it == node_to_super.end()) return {};
        const u64 super = it->second;
        if (super >= supernodes.size() || supernodes[super].name.empty())
            return {};  // supernode removed by post-layout cleaning
        char ori = '+';
        uint8_t side = 0;
        if (nodes_in_layouts.contains(id)) {
            auto oit = node_layout_ori.find(id);
            if (oit == node_layout_ori.end()) return {};
            const bool o = oit->second;
            ori = ori_char(o);
            side = is_src ? (o ? SIDE_END : SIDE_START)
                          : (o ? SIDE_START : SIDE_END);
            auto mit = node_boundary_mask.find(id);
            if (mit == node_boundary_mask.end() || (mit->second & side) == 0)
                return {};  // anchor is not on a free contig boundary
        }
        return {true, nodes_in_layouts.contains(id), super, ori, side};
    };

    using physical_endpoint = std::pair<u64, uint8_t>;
    using physical_pair =
        std::tuple<u64, uint8_t, u64, uint8_t>;
    using super_pair = std::pair<u64, u64>;
    struct projection {
        physical_pair endpoints;
        std::pair<std::optional<char>, std::optional<char>> orientations;
        u64 first_anchor = 0;
        u64 second_anchor = 0;
        bool first_endpoint_known = false;
        bool second_endpoint_known = false;
        bool c_orientation_unresolved = false;
    };
    auto flip = [](char ori) { return ori == '+' ? '-' : '+'; };
    auto project = [&](const resolved &src, const resolved &dst, u64 u,
                       u64 v) {
        const physical_endpoint src_ep{src.super, src.side};
        const physical_endpoint dst_ep{dst.super, dst.side};
        const std::optional<char> src_ori =
            src.orientation_known ? std::optional<char>{src.ori} : std::nullopt;
        const std::optional<char> dst_ori =
            dst.orientation_known ? std::optional<char>{dst.ori} : std::nullopt;
        if (src_ep <= dst_ep)
            return projection{{src.super, src.side, dst.super, dst.side},
                              {src_ori, dst_ori}, u, v,
                              src.orientation_known, dst.orientation_known,
                              !src.orientation_known || !dst.orientation_known};
        // Reversing a GFA link flips both segment orientations.
        return projection{{dst.super, dst.side, src.super, src.side},
                           {dst_ori ? std::optional<char>{flip(*dst_ori)}
                                    : std::nullopt,
                            src_ori ? std::optional<char>{flip(*src_ori)}
                                    : std::nullopt},
                           v, u,
                           dst.orientation_known, src.orientation_known,
                           !src.orientation_known || !dst.orientation_known};
    };
    auto projections_agree = [](const projection &a, const projection &b) {
        if (a.endpoints != b.endpoints) return false;
        const auto orientation_agrees = [](const std::optional<char> &x,
                                           const std::optional<char> &y) {
            return !x || !y || *x == *y;
        };
        return orientation_agrees(a.orientations.first,
                                  b.orientations.first) &&
               orientation_agrees(a.orientations.second,
                                  b.orientations.second);
    };
    auto merge_projection_orientation = [](projection &dst,
                                           const projection &src) {
        if (!dst.orientations.first)
            dst.orientations.first = src.orientations.first;
        if (!dst.orientations.second)
            dst.orientations.second = src.orientations.second;
        dst.first_endpoint_known |= src.first_endpoint_known;
        dst.second_endpoint_known |= src.second_endpoint_known;
        dst.c_orientation_unresolved |= src.c_orientation_unresolved;
    };

    struct candidate_group {
        projection projected;
        bool initialized = false;
        bool conflict = false;
        gtl::vector<bridge_read_gap_observation> read_gaps;
    };
    std::map<super_pair, candidate_group> groups;
    size_t bridge_record_index = 0;
    for (const auto &b : bridges) {
        if (b.u == b.v) continue;
        resolved su = resolve(b.u, /*is_src=*/true);
        resolved sv = resolve(b.v, /*is_src=*/false);
        if (!su.ok || !sv.ok) {
            ++result.invalid_endpoint;
            continue;
        }
        if (su.super == sv.super) {
            ++result.intra_supernode;
            continue;
        }
        const super_pair pair{std::min(su.super, sv.super),
                              std::max(su.super, sv.super)};
        const projection p = project(su, sv, b.u, b.v);
        auto &group = groups[pair];
        if (!group.initialized) {
            group.projected = p;
            group.initialized = true;
        } else if (!projections_agree(group.projected, p)) {
            group.conflict = true;
        } else {
            ++result.aggregated_records;
            merge_projection_orientation(group.projected, p);
        }
        if (!b.read_gap_observations.empty()) {
            group.read_gaps.insert(group.read_gaps.end(),
                                   b.read_gap_observations.begin(),
                                   b.read_gap_observations.end());
        } else {
            // Compatibility for API callers that only provide aggregate
            // support: represent each count as a record-local synthetic read via
            // a distinct qname hash so the collapse keeps them separate.
            // Production bridges always take the qname-preserving path above.
            for (uint32_t i = 0; i < b.support; ++i)
                group.read_gaps.push_back(
                    {(static_cast<u64>(bridge_record_index) << 32) | i, b.gap});
        }
        ++bridge_record_index;
    }

    struct candidate {
        edge_key key;
        edge_agg agg;
        physical_endpoint first_endpoint;
        physical_endpoint second_endpoint;
        bool first_endpoint_known = false;
        bool second_endpoint_known = false;
    };
    gtl::vector<candidate> candidates;
    std::map<physical_endpoint, size_t> endpoint_use;
    for (auto &[_, group] : groups) {
        if (group.conflict) {
            ++result.projection_conflict;
            continue;
        }
        auto [ss, sside, ds, dside] = group.projected.endpoints;
        bool equivalent_existing = false;
        bool conflicting_existing = false;
        for (const auto &[existing_key, existing] : merged_edges) {
            (void)existing_key;
            const projection ep = project(
                {true, existing.src_side != 0, existing.src_super,
                 existing.src_ori, existing.src_side},
                {true, existing.dst_side != 0, existing.dst_super,
                 existing.dst_ori, existing.dst_side},
                existing.example_from, existing.example_to);
            if (ep.endpoints != group.projected.endpoints) continue;
            if (projections_agree(ep, group.projected))
                equivalent_existing = true;
            else
                conflicting_existing = true;
        }
        if (equivalent_existing) {
            ++result.existing_edge;
            continue;
        }
        if (conflicting_existing) {
            ++result.projection_conflict;
            continue;
        }

        auto distinct_gaps =
            collapse_bridge_read_gaps_by_qname(std::move(group.read_gaps));
        const signed_gap_summary gap = summarize_signed_gaps(distinct_gaps);
        const uint32_t support = gap.count;
        const char sori = group.projected.orientations.first.value_or('+');
        const char dori = group.projected.orientations.second.value_or('+');
        edge_agg agg;
        agg.src_super = ss;
        agg.dst_super = ds;
        agg.src_ori = sori;
        agg.dst_ori = dori;
        agg.src_side = sside;
        agg.dst_side = dside;
        agg.example_from = group.projected.first_anchor;
        agg.example_to = group.projected.second_anchor;
        agg.count = support;
        agg.thread_bridge = true;
        agg.tb_support = support;
        agg.tb_gap = gap.median;
        agg.tb_gap_summary = gap;
        agg.tb_realization = classify_bridge_realization(gap);
        agg.tb_c_orientation_unresolved =
            group.projected.c_orientation_unresolved;
        edge_key key{ss, ds, sori, dori, sside, dside};
        physical_endpoint first{ss, sside}, second{ds, dside};
        if (group.projected.first_endpoint_known) ++endpoint_use[first];
        if (group.projected.second_endpoint_known) ++endpoint_use[second];
        candidates.push_back({key, agg, first, second,
                              group.projected.first_endpoint_known,
                              group.projected.second_endpoint_known});
    }

    for (const auto &c : candidates) {
        // 1:1 filter: each orientation-resolved endpoint is used by exactly
        // one candidate bridge. C/unknown side 0 is not a physical endpoint.
        if ((c.first_endpoint_known && endpoint_use[c.first_endpoint] != 1) ||
            (c.second_endpoint_known && endpoint_use[c.second_endpoint] != 1)) {
            ++result.endpoint_reuse;
            continue;
        }
        if (merged_edges.count(c.key)) {
            ++result.existing_edge;
            continue;
        }
        merged_edges[c.key] = c.agg;
        ++result.added;
    }
    std::fprintf(stderr,
                 "thread-bridge inject: added=%zu aggregated=%zu "
                 "invalid=%zu intra=%zu existing=%zu conflict=%zu reuse=%zu\n",
                 result.added, result.aggregated_records,
                 result.invalid_endpoint, result.intra_supernode,
                 result.existing_edge, result.projection_conflict,
                 result.endpoint_reuse);
    return result;
}

// ── Phase: format supernode S-lines (parallel) and write the GFA ────────
inline void write_gfa(const supernode_vec &supernodes,
                      const edge_map &merged_edges, int nthreads, FILE *out) {
    // Each worker converts a supernode's sequences to ASCII and builds the
    // complete `S` line (prefix + sequence + tags). Doing the whole line here
    // keeps the megabase `%s` formatting off the serial write path below.
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
                        // No gap fabrication: disconnected components are
                        // emitted as separate contigs (R7), so a zero-length
                        // entry contributes nothing rather than a run of 'N'.
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

    // Accumulate output into a 4 MiB buffer and flush in large fwrites to
    // minimize write() syscalls. Oversized lines (megabase S-lines) bypass the
    // buffer with a direct write.
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
    char tbuf[160];
    for (const auto &[_, edge] : merged_edges) {
        const auto &src_name = supernodes[edge.src_super].name;
        const auto &dst_name = supernodes[edge.dst_super].name;
        // Read-threaded scaffold edges carry TP:Z:bridge + TB:i:<support> so a
        // consumer can separate them from layout junction edges.
        tbuf[0] = '\0';
        if (edge.thread_bridge)
            std::snprintf(tbuf, sizeof(tbuf),
                          "\tTP:Z:bridge\tTB:i:%u\tGL:i:%" PRId64
                          "\tTR:Z:%s%s",
                          edge.tb_support, edge.tb_gap,
                          bridge_realization_label(edge.tb_realization),
                          edge.tb_c_orientation_unresolved
                              ? "\tCO:Z:unresolved"
                              : "");
        int n = std::snprintf(
            lbuf, sizeof(lbuf),
            "L\t%s\t%c\t%s\t%c\t0M\tSO:Z:%c\tTO:Z:%c\tSS:Z:%s\tTS:Z:%s\t"
            "EC:i:%zu\tFN:i:%" PRIu64 "\tTN:i:%" PRIu64 "%s\n",
            src_name.c_str(), edge.src_ori, dst_name.c_str(), edge.dst_ori,
            edge.src_ori, edge.dst_ori, side_label(edge.src_side),
            side_label(edge.dst_side), edge.count,
            static_cast<uint64_t>(edge.example_from),
            static_cast<uint64_t>(edge.example_to), tbuf);
        if (n > 0)
            emit(lbuf,
                 std::min<size_t>(static_cast<size_t>(n), sizeof(lbuf) - 1));
    }
    flush();
}

// ── Phase: materialize a scaffold FASTA from the bridge edges (opt-in) ───
// Walks the thread-bridge subgraph (only TP:Z:bridge superedges) through simple
// paths, concatenating the joined supernode sequences orientation-aware with an
// N-gap between; every unbridged live supernode is emitted as its own record, so
// the FASTA is a complete assembly (usable directly by QUAST). Branchy bridge
// components (a node with >2 bridge ends) are NOT scaffolded — their members are
// emitted standalone. Opt-in: does not touch the GFA path.
inline std::expected<void, std::string> write_thread_bridge_fasta(
    const supernode_vec &supernodes, const edge_map &merged_edges,
    const std::string &path) {
    using file_handle = std::unique_ptr<FILE, decltype(&std::fclose)>;
    file_handle fa{std::fopen(path.c_str(), "w"), &std::fclose};
    if (!fa) {
        return std::unexpected(
            "Cannot open --thread-bridge-fasta output: " + path);
    }
    const std::string paths_path = path + ".paths.tsv";
    file_handle paths{std::fopen(paths_path.c_str(), "w"), &std::fclose};
    if (!paths) {
        return std::unexpected(
            "Cannot open thread-bridge path sidecar: " + paths_path);
    }
    // Both required outputs are open before either receives content.
    std::fprintf(paths.get(),
                 "scaffold\torder\tsegment\torientation\tgap_before\t"
                 "support\tfrom_anchor\tto_anchor\tsegment_start\t"
                 "segment_end\n");
    auto flip = [](char o) { return o == '+' ? '-' : '+'; };

    // Supernode -> forward ASCII sequence (built lazily, reused across scaffolds).
    charsi chs;
    auto ascii_of = [&](u64 s) {
        std::string out;
        const auto &seqs = supernodes[s].seqs();
        size_t tot = 0;
        for (const auto &q : seqs) tot += q.len;
        out.reserve(tot);
        for (const auto &q : seqs) {
            if (q.len == 0) continue;
            chs.require(q.len + 32);
            q.convert_to_ascii((char *)chs.data);
            out.append((const char *)chs.data, q.len);
        }
        return out;
    };
    auto revcomp = [](const std::string &s) {
        std::string r(s.size(), 'N');
        for (size_t i = 0; i < s.size(); ++i) {
            switch (s[s.size() - 1 - i]) {
                case 'A': r[i] = 'T'; break;
                case 'C': r[i] = 'G'; break;
                case 'G': r[i] = 'C'; break;
                case 'T': r[i] = 'A'; break;
                default:  r[i] = 'N'; break;
            }
        }
        return r;
    };

    // Undirected bridge adjacency. GL remains an anchor-run separation rather
    // than an exact supernode-boundary distance. It is realized as an N-gap only
    // for edges that passed thread_bridge_is_fasta_safe's conservative gates.
    struct binc {
        u64 other;
        char xo;
        char oo;
        int64_t gap;
        uint32_t support;
        u64 from_anchor;
        u64 to_anchor;
    };
    gtl::flat_hash_map<u64, gtl::vector<binc>> badj;
    for (const auto &[k, e] : merged_edges) {
        (void)k;
        if (!thread_bridge_is_fasta_safe(e)) continue;
        badj[e.src_super].push_back(
            {e.dst_super, e.src_ori, e.dst_ori, e.tb_gap, e.tb_support,
             e.example_from, e.example_to});
        badj[e.dst_super].push_back(
            {e.src_super, flip(e.dst_ori), flip(e.src_ori), e.tb_gap,
             e.tb_support, e.example_to, e.example_from});
    }

    gtl::flat_hash_set<u64> done;  // supernodes already emitted
    size_t n_scaf = 0;

    // Emit scaffolds: start from each unused node that is a path endpoint
    // (exactly 1 bridge neighbor) in a strictly-simple component.
    auto component_is_simple = [&](u64 seed, gtl::vector<u64> &members) {
        // BFS the component; reject if any node has >2 bridge neighbors.
        gtl::flat_hash_set<u64> seen{seed};
        gtl::vector<u64> stack{seed};
        int endpoints = 0;
        bool simple = true;
        while (!stack.empty()) {
            u64 x = stack.back();
            stack.pop_back();
            members.push_back(x);
            const size_t deg = badj[x].size();
            if (deg > 2) simple = false;
            if (deg == 1) ++endpoints;
            for (const auto &b : badj[x])
                if (seen.insert(b.other).second) stack.push_back(b.other);
        }
        return simple && endpoints == 2;  // a path has exactly two ends
    };

    for (const auto &[start, incs] : badj) {
        if (done.contains(start) || incs.size() != 1) continue;
        gtl::vector<u64> members;
        if (!component_is_simple(start, members)) continue;
        // Walk the path. The start strand is NOT free: it must be the
        // orientation the start's (single) edge requires (xo), so the junction
        // faces the appended neighbor. Seeding '+' unconditionally produces an
        // anti-colinear order when the anchor sits on the start node's 5' side.
        std::string scaffold;
        u64 cur = start;
        char strand = badj[start].front().xo;
        int64_t pending_gap = 0;  // gap bases of the edge that led to `cur`
        uint32_t pending_support = 0;
        u64 pending_from_anchor = 0;
        u64 pending_to_anchor = 0;
        size_t path_order = 0;
        gtl::flat_hash_set<u64> visited;
        while (true) {
            visited.insert(cur);
            done.insert(cur);
            std::string s = ascii_of(cur);
            if (strand == '-') s = revcomp(s);
            if (!scaffold.empty())
                scaffold.append(static_cast<size_t>(pending_gap), 'N');
            const size_t segment_start = scaffold.size();
            scaffold += s;
            std::fprintf(paths.get(),
                         "scaffold_%zu\t%zu\t%s\t%c\t%" PRId64
                         "\t%u\t%" PRIu64 "\t%" PRIu64 "\t%zu\t%zu\n",
                         n_scaf, path_order++, supernodes[cur].name.c_str(),
                         strand, pending_gap, pending_support,
                         pending_from_anchor, pending_to_anchor, segment_start,
                         scaffold.size());
            // next unvisited bridge neighbor
            const binc *nxt = nullptr;
            for (const auto &b : badj[cur])
                if (!visited.contains(b.other)) { nxt = &b; break; }
            if (!nxt) break;
            char next_strand = (strand == nxt->xo) ? nxt->oo : flip(nxt->oo);
            pending_gap = nxt->gap;
            pending_support = nxt->support;
            pending_from_anchor = nxt->from_anchor;
            pending_to_anchor = nxt->to_anchor;
            cur = nxt->other;
            strand = next_strand;
        }
        std::fprintf(fa.get(), ">scaffold_%zu segments=%zu len=%zu\n%s\n",
                     n_scaf++, members.size(), scaffold.size(),
                     scaffold.c_str());
    }

    // Every remaining live supernode (unbridged, or a member of a rejected
    // branchy component) is emitted standalone so the FASTA is complete.
    for (u64 s = 0; s < supernodes.size(); ++s) {
        if (supernodes[s].name.empty() || done.contains(s)) continue;
        std::string seq = ascii_of(s);
        if (seq.empty()) continue;
        std::fprintf(fa.get(), ">%s len=%zu\n%s\n", supernodes[s].name.c_str(),
                     seq.size(), seq.c_str());
    }
    std::fprintf(stderr, "thread-bridge: scaffold FASTA (%zu scaffolds) -> %s\n",
                  n_scaf, path.c_str());
    return {};
}

template <class AssemblerT>
inline std::expected<void, std::string> emit_ultra_merged_graph(
    AssemblerT &as, dbg_index &dbg, FILE *out) {
    INIT_TIMERS(LayoutEmit);
    auto [recig, layouts] = as.fold_dbg_graph(dbg);
    TIME_CHECKPOINT(LayoutEmit, stderr, "fold_dbg_graph (%zu layouts)",
                    layouts.size());

    gtl::vector<supernode_out> supernodes;
    gtl::flat_hash_map<u64, u64> node_to_super;
    gtl::flat_hash_set<u64> nodes_in_layouts;
    gtl::flat_hash_map<u64, uint8_t> node_boundary_mask;
    node_to_super.reserve(dbg.active_set.size() + 8);
    nodes_in_layouts.reserve(dbg.active_set.size() + 8);

    // Resolve per-node orientation up front (deterministic; see
    // resolve_layout_orientations). Normally each node is in exactly one layout,
    // so there are no conflicts; a nonzero count means the single-owner invariant
    // was violated upstream.
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
    const int hash_combine_len = as.olc_lcmer_l;
    auto layout_seqs =
        build_layout_sequences_parallel(as, dbg, layout_jobs, hash_combine_len);
    TIME_CHECKPOINT(LayoutEmit, stderr, "Build layout sequences");
    for (size_t i = 0; i < layout_jobs.size(); ++i) {
        const auto &job = layout_jobs[i];
        auto &comps = layout_seqs[i];
        if (comps.empty()) continue;
        // One contig per SCD component: components are disconnected in the
        // seed-chain DAG, so no adjacency between them may be asserted (R7).
        // The first component keeps the plain "M<idx>" name; any further
        // component gets an "_<c>" suffix to keep GFA segment names unique.
        const u64 primary_sid = supernodes.size();
        for (size_t c = 0; c < comps.size(); ++c) {
            std::string name = "M" + std::to_string(job.layout_idx);
            if (c) name += "_" + std::to_string(c);
            supernodes.push_back({std::move(name), std::move(comps[c]), nullptr,
                                  true, job.layout_idx});
        }
        // Edges are aggregated per node; map every node of the job to the
        // primary component, exactly as when the job produced a single contig.
        // Secondary components carry no claimed adjacency.
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

    // ── Phase 5.5: inject read-threaded bridge scaffold edges ───────────
    // Added after cleaning so cleaning cannot cut them; endpoints are re-checked
    // against the cleaned supernode set inside the injector.
    if (!as.confirmed_bridges.empty()) {
        const bridge_injection_result injected = inject_thread_bridge_edges(
            as.confirmed_bridges, supernodes, node_to_super, nodes_in_layouts,
            node_layout_ori, node_boundary_mask, merged_edges);
        TIME_CHECKPOINT(LayoutEmit, stderr,
                        "Thread-bridge edges injected (%zu of %zu confirmed)",
                        injected.added, as.confirmed_bridges.size());
    }

    // ── Phase 6: format supernode S-lines (parallel) + write GFA ────────
    write_gfa(supernodes, merged_edges, as.threads.work, out);
    TIME_CHECKPOINT(LayoutEmit, stderr, "GFA writing");

    // ── Phase 7: opt-in scaffold FASTA materialized from bridge edges ───
    if (!as.thread_bridge.fasta_path.empty()) {
        if (auto result = write_thread_bridge_fasta(
                supernodes, merged_edges, as.thread_bridge.fasta_path);
            !result)
            return result;
        TIME_CHECKPOINT(LayoutEmit, stderr, "Scaffold FASTA");
    }
    return {};
}

template <class AssemblerT>
inline std::expected<void, std::string> emit_ultra_merged_graph(
    AssemblerT &as, dbg_index &dbg) {
    return emit_ultra_merged_graph(as, dbg, stdout);
}

}  // namespace layout_graph_emit

#endif  // _LAYOUT_GRAPH_H
