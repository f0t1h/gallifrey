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
#include <queue>
// #include <stack>
#include <string>
#include <thread>
#include <tuple>

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
    [[maybe_unused]] int hash_combine_len, u64 layout_idx = 0,
    const gtl::flat_hash_set<u64> *adj_crumbs = nullptr) {
    if (lnodes.empty()) return {};

    // ── Build lane_nodes ────────────────────────────────────────────
    gtl::vector<lane_node> nodes;
    gtl::vector<lane_node> crumbs;  // adj-joined crumbs (emitted via abutment)
    nodes.reserve(lnodes.size());
    for (const auto &n : lnodes) {
        auto it = dbg.dbg.find(n.id);
        if (it == dbg.dbg.end()) continue;
        if (it->second.seq.empty()) continue;
        if (adj_crumbs && adj_crumbs->contains(n.id)) {
            crumbs.push_back({n.id, n.from, n.to, n.ori, &it->second.seq});
            continue;
        }
        nodes.push_back({n.id, n.from, n.to, n.ori, &it->second.seq});
    }
    gtl::vector<size_t> sorted_idx(nodes.size());
    for (size_t i = 0; i < sorted_idx.size(); ++i) sorted_idx[i] = i;
    std::sort(sorted_idx.begin(), sorted_idx.end(),
              [&](size_t a, size_t b) {
                  return nodes[a].from < nodes[b].from;
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


    if (getenv("LASM_ADJ_JOIN") && result.size() > 1 && !crumbs.empty()) {
        // Data-supported abutment: SCD components separated only by
        // adj-joined crumbs are positionally adjacent; the crumb edges are
        // read evidence for the join (unlike the fabricated N-gap of R7).
        std::sort(crumbs.begin(), crumbs.end(),
                  [](const lane_node &a, const lane_node &b) {
                      return a.from < b.from;
                  });
        for (size_t i = 0; i + 1 < comp_paths.size();) {
            const i64 gap_lo = comp_paths[i].layout_to;
            const i64 gap_hi = comp_paths[i + 1].layout_from;
            if (gap_lo >= gap_hi) {
                ++i;
                continue;
            }
            gtl::vector<const lane_node *> gap_crumbs;
            for (const auto &c : crumbs) {
                if (c.from >= gap_lo && c.to <= gap_hi) {
                    gap_crumbs.push_back(&c);
                }
            }
            if (gap_crumbs.empty()) {
                ++i;
                continue;
            }
            for (const lane_node *c : gap_crumbs) {
                if (c->is_forward_lane) {
                    for (const auto &s : *c->raw_seqs) {
                        result[i].push_back(
                            SequenceContainerType::copy_of(s));
                    }
                } else {
                    gtl::vector<SequenceContainerType> rc =
                        reverse_complement_seqs(*c->raw_seqs);
                    for (auto &s : rc) result[i].push_back(std::move(s));
                }
            }
            for (auto &s : result[i + 1]) {
                result[i].push_back(std::move(s));
            }
            comp_paths[i].layout_to = comp_paths[i + 1].layout_to;
            result.erase(result.begin() + i + 1);
            comp_paths.erase(comp_paths.begin() + i + 1);
        }
    }

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
    int hash_combine_len, const gtl::flat_hash_set<u64> *adj_crumbs = nullptr) {
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
                    jobs[job_idx].layout_idx, adj_crumbs);
            }
            while (job_queue.try_dequeue(job_idx)) {
                out_seqs[job_idx] = build_layout_sequence_from_nodes(
                    as, dbg, jobs[job_idx].lnodes, hash_combine_len,
                    jobs[job_idx].layout_idx, adj_crumbs);
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
        } else {
            wbuf.append(data, len);
        }
    };
    auto emit_final = [&]() {
        flush();
    };

    for (size_t i = 0; i < n_super; ++i) {
        if (super_line_strs[i].empty()) continue;  // removed
        emit(super_line_strs[i].data(), super_line_strs[i].size());
    }
    emit_final();


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

template <class AssemblerT>
inline void emit_ultra_merged_graph(AssemblerT &as, dbg_index &dbg, FILE *out) {
    INIT_TIMERS(LayoutEmit);
    auto [recig, layouts] = as.fold_dbg_graph(dbg);
    TIME_CHECKPOINT(LayoutEmit, stderr, "fold_dbg_graph (%zu layouts)",
                    layouts.size());

    // Adj-joined crumbs: reciprocal-map entries whose query span is below
    // the lcmer threshold (fold_dbg_graph's LASM_ADJ_JOIN block). Emitted by
    // abutment between SCD components, not as their own contigs.
    gtl::flat_hash_set<u64> adj_crumbs;
    if (getenv("LASM_ADJ_JOIN")) {
        for (const auto &[k, ov] : recig) {
            if (ov.query.len < static_cast<u64>(as.olc_lcmer_l + 1)) {
                adj_crumbs.insert(k);
            }
        }
    }
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
        build_layout_sequences_parallel(as, dbg, layout_jobs, hash_combine_len,
                                        &adj_crumbs);
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

    // ── Phase 6: format supernode S-lines (parallel) + write GFA ────────
    write_gfa(supernodes, merged_edges, as.threads.work, out);
    TIME_CHECKPOINT(LayoutEmit, stderr, "GFA writing");
}

template <class AssemblerT>
inline void emit_ultra_merged_graph(AssemblerT &as, dbg_index &dbg) {
    emit_ultra_merged_graph(as, dbg, stdout);
}

}  // namespace layout_graph_emit

#endif  // _LAYOUT_GRAPH_H
