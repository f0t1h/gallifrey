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
    // Mirrors the original precedence: a node's length is the `qry.size` of the
    // first entry keyed by that node (when > 0); otherwise the `tgt.size` of the
    // first entry (in iteration order) that references it as a target.
    gtl::flat_hash_map<u64, i64> first_qry_len;  // key node → its first qry.size
    gtl::flat_hash_map<u64, i64> first_tgt_len;  // target node → first tgt.size>0
    for (const auto &kv : recig) {
        const auto &rm = kv.second;
        first_qry_len.try_emplace(kv.first, static_cast<i64>(rm.qry.size));
        if (rm.tgt.size > 0)
            first_tgt_len.try_emplace(rm.node, static_cast<i64>(rm.tgt.size));
    }

    auto infer_len_for_node = [&](u64 id) -> i64 {
        auto qit = first_qry_len.find(id);
        if (qit != first_qry_len.end() && qit->second > 0) return qit->second;
        auto tit = first_tgt_len.find(id);
        if (tit != first_tgt_len.end()) return tit->second;
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
            u64 tgt = node.second.node;
            if (visited.contains(tgt)) continue;
            visited.emplace(tgt);

            auto [spos, epos, s_orient] = layouts.back().at(src);
            (void)epos;
            bool target_orient = !s_orient;

            long long qs = s_orient
                               ? (node.second.qry.size - node.second.qry.to)
                               : node.second.qry.from;
            long long ts = target_orient
                               ? (node.second.tgt.size - node.second.tgt.to)
                               : node.second.tgt.from;

            i64 delta = static_cast<i64>(qs - ts);
            i64 target_start = spos + delta;
            const i64 target_len = infer_len_for_node(tgt);
            const i64 target_end = to_end_inclusive(target_start, target_len);

            layouts.back().emplace(
                tgt, std::make_tuple(target_start, target_end, target_orient));

            auto range = recig.equal_range(tgt);
            for (auto it = range.first; it != range.second; ++it) {
                if (!visited.contains(it->second.node)) q.push(*it);
            }
        }
    }
    return layouts;
}

}  // namespace layout_graph_algo

namespace layout_graph_emit {
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
    scd.reserve(contig_order.size() * 64);
    scd_fwd.reserve(contig_order.size() * 64);
    node_contigs.reserve(contig_order.size() * 64);

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

// ── Phase 2: Find heaviest path per connected component ─────────────
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

        gtl::vector<uint32_t> path;
        for (int32_t cur = best_end; cur >= 0; cur = parent[cur])
            path.push_back(static_cast<uint32_t>(cur));
        std::reverse(path.begin(), path.end());

        // Extend at both ends
        {
            gtl::flat_hash_set<uint32_t> on_path(path.begin(), path.end());
            i64 path_max_to = INT64_MIN;
            for (uint32_t nid : path)
                path_max_to = std::max(path_max_to,
                                       nodes[scd[nid].contig_idx].to);
            i64 path_min_from = INT64_MAX;
            for (uint32_t nid : path)
                path_min_from = std::min(path_min_from,
                                         nodes[scd[nid].contig_idx].from);

            // Forward extension
            int32_t ext_sink = -1;
            i64 ext_to = path_max_to;
            for (size_t i = 0; i < n; ++i) {
                if (comp_id[i] != c) continue;
                if (on_path.contains(static_cast<uint32_t>(i))) continue;
                i64 to = nodes[scd[i].contig_idx].to;
                if (to > ext_to && scd_fwd[i].empty()) {
                    ext_to = to;
                    ext_sink = static_cast<int32_t>(i);
                }
            }
            if (ext_sink >= 0) {
                gtl::vector<uint32_t> ext_path;
                for (int32_t cur = ext_sink; cur >= 0; cur = parent[cur])
                    ext_path.push_back(static_cast<uint32_t>(cur));
                std::reverse(ext_path.begin(), ext_path.end());
                size_t branch_pos = 0;
                bool found_overlap = false;
                for (size_t i = 0; i < ext_path.size(); ++i) {
                    if (on_path.contains(ext_path[i])) {
                        branch_pos = i;
                        found_overlap = true;
                    }
                }
                if (!found_overlap) [[unlikely]] {
                    std::fprintf(stderr,
                        "[layout-dag] WARNING: forward extension has no "
                        "overlap with main path (%s:%d)\n",
                        __FILE__, __LINE__);
                } else {
                    for (size_t i = branch_pos + 1; i < ext_path.size(); ++i)
                        path.push_back(ext_path[i]);
                }
            }

            // Backward extension
            int32_t ext_source = -1;
            i64 ext_from = path_min_from;
            for (size_t i = 0; i < n; ++i) {
                if (comp_id[i] != c) continue;
                if (on_path.contains(static_cast<uint32_t>(i))) continue;
                i64 from = nodes[scd[i].contig_idx].from;
                if (from < ext_from && scd_rev[i].empty()) {
                    ext_from = from;
                    ext_source = static_cast<int32_t>(i);
                }
            }
            if (ext_source >= 0) {
                // Greedily walk max-dp successors from the earlier source until
                // we rejoin the main path. Only prepend if we actually reconnect
                // (mirrors the forward extension's overlap guard) — otherwise the
                // prefix would splice a disconnected chain onto the layout. The
                // `walked` set bounds the walk so a cyclic SCD can't loop forever.
                gtl::vector<uint32_t> prefix;
                gtl::flat_hash_set<uint32_t> walked;
                uint32_t cur = static_cast<uint32_t>(ext_source);
                bool reconnected = true;
                while (!on_path.contains(cur)) {
                    if (!walked.insert(cur).second) {  // revisited → cycle
                        reconnected = false;
                        break;
                    }
                    prefix.push_back(cur);
                    uint32_t best_succ = UINT32_MAX;
                    int64_t best_dp_val = -1;
                    for (auto &[v, w] : scd_fwd[cur]) {
                        if (dp[v] > best_dp_val) {
                            best_dp_val = dp[v];
                            best_succ = v;
                        }
                    }
                    if (best_succ == UINT32_MAX) {  // dead-end before path
                        reconnected = false;
                        break;
                    }
                    cur = best_succ;
                }
                if (reconnected && !prefix.empty()) {
                    prefix.insert(prefix.end(), path.begin(), path.end());
                    path = std::move(prefix);
                }
            }
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
    return comp_paths;
}

// ── Phase 3: Collect output sequences from heaviest paths ───────────
inline gtl::vector<SequenceContainerType> collect_path_sequences(
    const gtl::vector<comp_path> &paths,
    const gtl::vector<gtl::vector<std::pair<size_t, size_t>>> &node_contigs,
    gtl::flat_hash_map<size_t, gtl::vector<SequenceContainerType>> &contig_core_seqs,
    size_t n_nodes)
{
    gtl::vector<SequenceContainerType> merged_seqs;
    merged_seqs.reserve(n_nodes);

    for (size_t ci = 0; ci < paths.size(); ++ci) {
        if (ci > 0 && !merged_seqs.empty())
            merged_seqs.emplace_back(0, 0);

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
                // path nodes stay on the same contig; otherwise fall back to
                // this node's first contig.
                auto pick_contig = [&]() -> std::pair<size_t, size_t> {
                    if (pi + 1 < path.size()) {
                        const auto &nc_next = node_contigs[path[pi + 1]];
                        for (auto &[cid, kid] : nc) {
                            for (auto &[cid2, kid2] : nc_next) {
                                if (cid == cid2) return {cid, kid};
                            }
                        }
                    }
                    if (!nc.empty()) return {nc[0].first, nc[0].second};
                    return {SIZE_MAX, SIZE_MAX};
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
    }
    return merged_seqs;
}

// ── Main entry: build layout sequence from positioned nodes ─────────
template <class AssemblerT>
inline gtl::vector<SequenceContainerType> build_layout_sequence_from_nodes(
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

            if (!nodes[pick].is_forward_lane) {
                const auto &raw = *nodes[pick].raw_seqs;
                gtl::vector<SequenceContainerType> rev;
                rev.reserve(raw.size());
                for (auto rit = raw.rbegin(); rit != raw.rend(); ++rit)
                    rev.emplace_back(rit->reverse_complement());
                return rev;
            }
            const auto &raw = *nodes[pick].raw_seqs;
            gtl::vector<SequenceContainerType> result;
            result.reserve(1);
            result.emplace_back(raw.begin(), raw.end());
            return result;
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
            rev_seqs.reserve(raw.size());
            for (auto rit = raw.rbegin(); rit != raw.rend(); ++rit)
                rev_seqs.emplace_back(rit->reverse_complement());
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

using layout_seq_t = gtl::vector<SequenceContainerType>;

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
                    size_t max_tip_len = 50000) {
        adj_map adj;
        build_adjacency(merged_edges, adj);

        size_t tips = cut_tips(supernodes, merged_edges, adj, max_tip_len);
        fix_edges(supernodes, merged_edges);

        if (tips > 0)
            std::fprintf(stderr, "[post-layout-clean] Removed %zu tips\n",
                         tips);
    }
};

template <class AssemblerT>
inline void emit_ultra_merged_graph(AssemblerT &as, dbg_index &dbg, FILE *out) {
    INIT_TIMERS(LayoutEmit);
    auto [recig, layouts] = as.fold_dbg_graph(dbg);
    TIME_CHECKPOINT(LayoutEmit, stderr, "fold_dbg_graph (%zu layouts)",
                    layouts.size());

    static constexpr uint8_t SIDE_START = 1;
    static constexpr uint8_t SIDE_END = 2;
    auto side_label = [&](uint8_t side) -> const char * {
        if (side == SIDE_START) return "start";
        if (side == SIDE_END) return "end";
        return "na";
    };
    auto ori_char = [](bool ori) -> char { return ori ? '+' : '-'; };

    gtl::vector<supernode_out> supernodes;
    gtl::flat_hash_map<u64, u64> node_to_super;
    gtl::flat_hash_set<u64> nodes_in_layouts;
    gtl::flat_hash_map<u64, bool> node_layout_ori;
    gtl::flat_hash_map<u64, uint8_t> node_boundary_mask;
    node_to_super.reserve(dbg.active_set.size() + 8);
    nodes_in_layouts.reserve(dbg.active_set.size() + 8);

    const int hash_combine_len = as.olc_lcmer_l;
    gtl::vector<ultra_layout_job> layout_jobs;
    layout_jobs.reserve(layouts.size());

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
        const i64 fuzzy = std::max<i64>(as.olc_lcmer_l * 2, span / 20);
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
            auto [ori_it, inserted] = node_layout_ori.try_emplace(n.id, n.ori);
            if (!inserted && ori_it->second != n.ori) [[unlikely]] {
                std::fprintf(
                    stderr,
                    "[layout-emit] WARNING: node %" PRIu64
                    " orientation conflict: was %c, now %c (overwriting)\n",
                    n.id, ori_it->second ? '+' : '-', n.ori ? '+' : '-');
                ori_it->second = n.ori;
            }
            nodes_in_layouts.insert(n.id);
        }

        layout_jobs.push_back({layout_idx, std::move(lnodes)});
        ++layout_idx;
    }

    auto layout_seqs =
        build_layout_sequences_parallel(as, dbg, layout_jobs, hash_combine_len);
    TIME_CHECKPOINT(LayoutEmit, stderr, "Build layout sequences");
    for (size_t i = 0; i < layout_jobs.size(); ++i) {
        const auto &job = layout_jobs[i];
        const u64 sid = supernodes.size();
        supernodes.push_back({"M" + std::to_string(job.layout_idx),
                              std::move(layout_seqs[i]), nullptr, true,
                              job.layout_idx});
        for (const auto &n : job.lnodes) {
            node_to_super[n.id] = sid;
        }
    }

    // ── Skip C-nodes whose DBG neighbors are all inside layouts ────────
    // Build successor map from prev edges
    gtl::flat_hash_map<u64, gtl::vector<u64>> successors;
    for (const auto &[to_id, to_node] : dbg.dbg) {
        if (to_node.redirected()) continue;
        for (u64 from_id : to_node.prev) {
            successors[from_id].push_back(to_id);
        }
    }

    // Remove C-node connected components whose external edges all go to
    // layout-interior nodes. Such components are tips/paths hanging off
    // interior layout unitigs — their sequence is already in the layout.
    gtl::flat_hash_set<u64> candidate_set;
    for (u64 node_id : dbg.active_set) {
        if (nodes_in_layouts.contains(node_id)) continue;
        auto it = dbg.dbg.find(node_id);
        if (it == dbg.dbg.end()) continue;
        if (it->second.redirected() || it->second.seq.empty()) continue;
        candidate_set.insert(node_id);
    }

    // Find connected components among C-nodes
    gtl::flat_hash_set<u64> cc_visited;
    gtl::flat_hash_set<u64> contained_nodes;

    for (u64 seed : candidate_set) {
        if (cc_visited.contains(seed)) continue;

        // BFS within C-nodes
        gtl::vector<u64> component;
        std::queue<u64> bfs;
        bfs.push(seed);
        cc_visited.insert(seed);

        bool has_boundary_neighbor = false;

        while (!bfs.empty()) {
            u64 cur = bfs.front();
            bfs.pop();
            component.push_back(cur);
            auto db_it = dbg.dbg.find(cur);

            // Check all neighbors
            for (u64 pred_id : db_it->second.prev) {
                if (candidate_set.contains(pred_id)) {
                    if (!cc_visited.contains(pred_id)) {
                        cc_visited.insert(pred_id);
                        bfs.push(pred_id);
                    }
                } else if (!nodes_in_layouts.contains(pred_id)) {
                    // Neighbor is outside layouts and not a C-node candidate
                    has_boundary_neighbor = true;
                }
                // If in layout: check if it's at a layout boundary
                if (nodes_in_layouts.contains(pred_id)) {
                    auto bm_it = node_boundary_mask.find(pred_id);
                    if (bm_it != node_boundary_mask.end()) {
                        bool ori = node_layout_ori.at(pred_id);
                        uint8_t needed = ori ? SIDE_END : SIDE_START;
                        if (bm_it->second & needed)
                            has_boundary_neighbor = true;
                    }
                }
            }
            for (u64 succ_id : successors[cur]) {
                if (candidate_set.contains(succ_id)) {
                    if (!cc_visited.contains(succ_id)) {
                        cc_visited.insert(succ_id);
                        bfs.push(succ_id);
                    }
                } else if (!nodes_in_layouts.contains(succ_id)) {
                    has_boundary_neighbor = true;
                }
                if (nodes_in_layouts.contains(succ_id)) {
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

        if (!has_boundary_neighbor && !component.empty()) {
            for (u64 nid : component) contained_nodes.insert(nid);
        }
    }

    size_t contained_bp = 0;
    for (u64 node_id : contained_nodes) {
        auto it = dbg.dbg.find(node_id);
        for (const auto &s : it->second.seq) contained_bp += s.len;
    }
    if (!contained_nodes.empty()) {
        std::fprintf(stderr,
                     "[layout-emit] Skipped %zu contained C-nodes (%zu bp)\n",
                     contained_nodes.size(), contained_bp);
    }

    // Create supernodes for remaining C-nodes
    for (u64 node_id : candidate_set) {
        if (contained_nodes.contains(node_id)) continue;
        auto it = dbg.dbg.find(node_id);
        if (it == dbg.dbg.end()) continue;

        const u64 sid = supernodes.size();
        supernodes.push_back({"C" + std::to_string(node_id),
                              {},
                              &it->second.seq,
                              false,
                              node_id});
        node_to_super[node_id] = sid;
    }

    TIME_CHECKPOINT(LayoutEmit, stderr, "Supernode creation");
    edge_map merged_edges;
    {
        int nthreads_edge = as.threads.work;
        if (nthreads_edge < 1) nthreads_edge = 1;
        gtl::vector<std::map<edge_key, edge_agg>> per_thread(nthreads_edge);

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

                    edge_key k{from_super_it->second,
                               to_super_it->second,
                               src_ori,
                               dst_ori,
                               src_side,
                               dst_side};
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
            nthreads_edge);

        // Merge per-thread maps
        for (auto &local : per_thread) {
            for (auto &[k, agg] : local) {
                auto &dst = merged_edges[k];
                if (dst.count == 0) {
                    dst = agg;
                } else {
                    dst.count += agg.count;
                }
            }
        }
    }

    TIME_CHECKPOINT(LayoutEmit, stderr, "Edge aggregation");

    // ── Post-layout graph cleaning ──────────────────────────────────
    PostLayoutCleaning::run(supernodes, merged_edges);
    TIME_CHECKPOINT(LayoutEmit, stderr, "Post-layout cleaning");

    // ── Parallel S-line formatting ───────────────────────────────────
    // Each worker converts a supernode's sequences to ASCII and builds the
    // complete `S` line (prefix + sequence + tags). Doing the whole line here
    // keeps the megabase `%s` formatting off the serial write path below.
    constexpr size_t N_GAP_LEN = 64;
    const size_t n_super = supernodes.size();
    gtl::vector<std::string> super_line_strs(n_super);
    {
        int nthreads_print = as.threads.work;
        if (nthreads_print < 1) nthreads_print = 1;
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
                    for (const auto &s : sn_seqs)
                        total_len += (s.len == 0) ? N_GAP_LEN : s.len;

                    dst.reserve(total_len + sn.name.size() + 48);
                    dst.append("S\t");
                    dst.append(sn.name);
                    dst.push_back('\t');

                    const size_t seq_start = dst.size();
                    for (const auto &s : sn_seqs) {
                        if (s.len == 0) {
                            dst.append(N_GAP_LEN, 'N');
                            continue;
                        }
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
    TIME_CHECKPOINT(LayoutEmit, stderr, "Format S-lines");

    // ── Buffered write ───────────────────────────────────────────────
    // Accumulate output into a 4 MiB buffer and flush in large fwrites to
    // minimize write() syscalls (the serial write was I/O-bound at ~60% util).
    // Oversized lines (megabase S-lines) bypass the buffer with a direct write.
    {
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
            // Clamp to what actually fit (snprintf returns the untruncated
            // length); names are small so truncation should never happen.
            if (n > 0)
                emit(lbuf, std::min<size_t>(static_cast<size_t>(n),
                                            sizeof(lbuf) - 1));
        }
        flush();
    }
    TIME_CHECKPOINT(LayoutEmit, stderr, "GFA writing");
}

template <class AssemblerT>
inline void emit_ultra_merged_graph(AssemblerT &as, dbg_index &dbg) {
    emit_ultra_merged_graph(as, dbg, stdout);
}

}  // namespace layout_graph_emit

#endif  // _LAYOUT_GRAPH_H
