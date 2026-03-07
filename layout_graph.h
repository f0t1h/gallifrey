#ifndef _LAYOUT_GRAPH_H
#define _LAYOUT_GRAPH_H

#include <algorithm>
#include <atomic>
#include <chrono>
#include <cstdio>
#include <cstring>
#include <sys/stat.h>
#include <functional>
#include <gtl/phmap_fwd_decl.hpp>
#include <gtl/vector.hpp>
#include <Eigen/Dense>
#include <cmath>
#include <map>
#include <optional>
#include <queue>
#include <stack>
#include <stdexcept>
#include <string>
#include <thread>
#include <tuple>

#include "graph.h"

using layout_small_graph = gtl::flat_hash_map<uint64_t, dbg_index::bidbg_val>;

// Break cycles by duplicating the cycling sub-path.
// Each cycle found via DFS is unrolled once: nodes are copied with fresh IDs
// and the back-edge is redirected to the copy, turning the cycle into a linear
// tail.
inline void unroll_cycles(layout_small_graph &lsg) {
    u64 next_id = 0;
    for (const auto &[id, _] : lsg) next_id = std::max(next_id, id);
    ++next_id;

    for (;;) {
        // Build forward adjacency from prev lists.
        gtl::flat_hash_map<u64, gtl::vector<u64>> fwd;
        for (const auto &[id, val] : lsg)
            for (u64 p : val.prev)
                if (lsg.contains(p)) fwd[p].push_back(id);

        // DFS to find one cycle.
        gtl::flat_hash_set<u64> visited, on_stack;
        gtl::flat_hash_map<u64, u64> parent;
        gtl::vector<u64> cycle;
        bool found = false;

        std::function<void(u64)> dfs = [&](u64 u) {
            if (found) return;
            visited.insert(u);
            on_stack.insert(u);
            for (u64 v : fwd[u]) {
                if (found) return;
                if (on_stack.contains(v)) {
                    // Back-edge u→v. Extract cycle: v → ... → u
                    cycle.push_back(v);
                    for (u64 cur = u; cur != v; cur = parent[cur])
                        cycle.push_back(cur);
                    std::reverse(cycle.begin(), cycle.end());
                    found = true;
                    return;
                }
                if (!visited.contains(v)) {
                    parent[v] = u;
                    dfs(v);
                }
            }
            on_stack.erase(u);
        };

        for (const auto &[id, _] : lsg) {
            if (!visited.contains(id)) dfs(id);
            if (found) break;
        }

        if (!found) break;  // DAG — done.

        fprintf(stderr, "Unrolling cycle of %zu nodes (graph size: %zu)\n",
                cycle.size(), lsg.size());

        // Duplicate cycle nodes with fresh IDs.
        // First pass: collect data from originals (before any inserts
        // invalidate refs).
        struct dup_info {
            u64 orig_id, dup_id;
            unsigned core_count;
            gtl::vector<SequenceContainerType> seq;
        };
        gtl::vector<dup_info> dups;
        dups.reserve(cycle.size());
        gtl::flat_hash_map<u64, u64> orig_to_dup;
        for (u64 id : cycle) {
            u64 dup_id = next_id++;
            orig_to_dup[id] = dup_id;
            auto &orig = lsg.find(id)->second;
            gtl::vector<SequenceContainerType> seq_copy;
            seq_copy.reserve(orig.seq.size());
            for (const auto &s : orig.seq) {
                size_t n_words =
                    SequenceContainerType::storage_size(s.len + s.skip_first_n);
                uint64_t *buf = new uint64_t[n_words];
                std::memcpy(buf, s.data, n_words * sizeof(uint64_t));
                seq_copy.emplace_back(s.len, s.skip_first_n, buf,
                                      s.origin_node);
            }
            dups.push_back(
                {id, dup_id, static_cast<unsigned>(orig.counts.core_count),
                 std::move(seq_copy)});
        }
        // Second pass: insert duplicates into the map.
        for (auto &d : dups) {
            auto [it, _ins] = lsg.emplace(d.dup_id, d.core_count);
            it->second.seq = std::move(d.seq);
        }

        // Wire the duplicated chain.
        // dup[0].prev = {cycle.back()}  — the back-edge now leads into the copy
        lsg.find(orig_to_dup[cycle[0]])->second.prev = {cycle.back()};
        for (size_t i = 1; i < cycle.size(); ++i)
            lsg.find(orig_to_dup[cycle[i]])->second.prev = {
                orig_to_dup[cycle[i - 1]]};

        // Remove the back-edge from the original: drop cycle.back() from
        // cycle[0].prev
        auto &head_prev = lsg.find(cycle[0])->second.prev;
        head_prev.erase(
            std::remove(head_prev.begin(), head_prev.end(), cycle.back()),
            head_prev.end());
    }
}

// Find the longest path in the DAG by node count.
// Returns pairs of (node_id, sequences) along the path (source → sink order).
inline gtl::vector<std::pair<u64, gtl::vector<SequenceContainerType>>>
longest_path(layout_small_graph &lsg) {
    if (lsg.empty()) return {};

    // unroll_cycles(lsg);

    // Build forward adjacency list and in-degree map.
    gtl::flat_hash_map<u64, gtl::vector<u64>> fwd;
    gtl::flat_hash_map<u64, size_t> in_degree;

    for (const auto &[id, val] : lsg) {
        in_degree[id];
        fwd[id];
    }

    for (const auto &[id, val] : lsg) {
        for (u64 p : val.prev) {
            if (lsg.contains(p)) {
                fwd[p].push_back(id);
                in_degree[id]++;
            }
        }
    }
    // Kahn's BFS with DP for longest distance.
    gtl::flat_hash_map<u64, size_t> dist;
    gtl::flat_hash_map<u64, u64> pred;

    std::queue<u64> q;
    size_t processed = 0;
    for (const auto &[id, deg] : in_degree) {
        dist[id] = 1;
        if (deg == 0) {
            q.push(id);
        }
    }

    u64 best_node = lsg.begin()->first;
    size_t best_dist = 1;

    while (!q.empty()) {
        u64 u = q.front();
        q.pop();
        ++processed;

        if (dist[u] > best_dist) {
            best_dist = dist[u];
            best_node = u;
        }

        for (u64 v : fwd[u]) {
            if (dist[u] + 1 > dist[v]) {
                dist[v] = dist[u] + 1;
                pred[v] = u;
            }
            if (--in_degree[v] == 0) {
                q.push(v);
            }
        }
    }

    if (processed < lsg.size()) {
        fprintf(stderr,
                "Warning: layout graph contains cycles (%zu of %zu nodes "
                "unreachable)\n",
                lsg.size() - processed, lsg.size());
    }

    // Reconstruct path and collect node_id + sequences (source → sink).
    gtl::vector<u64> path;
    for (u64 cur = best_node;;) {
        path.push_back(cur);
        auto it = pred.find(cur);
        if (it == pred.end()) break;
        cur = it->second;
    }
    std::reverse(path.begin(), path.end());

    gtl::vector<std::pair<u64, gtl::vector<SequenceContainerType>>> result;
    result.reserve(path.size());
    for (u64 id : path) {
        auto it = lsg.find(id);
        result.emplace_back(id, std::move(it->second.seq));
    }
    return result;
}

struct layout_boundary {
    gtl::flat_hash_set<u64>
        incoming_fwd;  // predecessors of forward-oriented layout contigs
    gtl::flat_hash_set<u64>
        incoming_rev;  // predecessors of reverse-oriented layout contigs
    gtl::flat_hash_set<u64>
        outgoing_fwd;  // successors of forward-oriented layout contigs
    gtl::flat_hash_set<u64>
        outgoing_rev;  // successors of reverse-oriented layout contigs
    gtl::vector<u64> contained;  // non-layout nodes fully between layout nodes

    // Convenience: all incoming / all outgoing regardless of orientation.
    gtl::flat_hash_set<u64> all_incoming() const {
        gtl::flat_hash_set<u64> s(incoming_fwd);
        s.insert(incoming_rev.begin(), incoming_rev.end());
        return s;
    }
    gtl::flat_hash_set<u64> all_outgoing() const {
        gtl::flat_hash_set<u64> s(outgoing_fwd);
        s.insert(outgoing_rev.begin(), outgoing_rev.end());
        return s;
    }
};

// Given the original graph, active node set, and a map of layout node IDs to
// their orientation (true = forward, false = reverse), find non-layout nodes
// that border or are contained within the layout.
// Incoming = non-layout nodes in a layout node's prev list (topological
// predecessors). Outgoing = non-layout nodes whose prev list contains a layout
// node (topological successors). Each set is split by the orientation of the
// connected layout contig.
inline layout_boundary find_layout_boundary(
    const dbg_index::dbg_map_t &dbg, const auto &active_nodes,
    const gtl::flat_hash_map<u64, bool> &layout_ori) {
    layout_boundary result;

    // Build a flat set for quick membership checks.
    gtl::flat_hash_set<u64> layout_set;
    layout_set.reserve(layout_ori.size());
    for (const auto &[id, _] : layout_ori) {
        layout_set.insert(id);
    }

    // 1. Incoming: scan layout nodes' prev — active predecessors not in layout.
    //    Split by the orientation of the layout contig they connect to.
    for (const auto &[lid, ori] : layout_ori) {
        auto it = dbg.find(lid);
        if (it == dbg.end()) continue;
        for (u64 p : it->second.prev) {
            if (!layout_set.contains(p) && active_nodes.contains(p)) {
                if (ori) {
                    result.incoming_fwd.insert(p);
                } else {
                    result.incoming_rev.insert(p);
                }
            }
        }
    }

    // 2. Outgoing: active non-layout nodes whose prev contains a layout node.
    //    Split by the orientation of the layout contig they point back to.
    for (u64 id : active_nodes) {
        if (layout_set.contains(id)) continue;
        auto it = dbg.find(id);
        if (it == dbg.end()) continue;
        for (u64 p : it->second.prev) {
            auto oit = layout_ori.find(p);
            if (oit != layout_ori.end()) {
                if (oit->second) {
                    result.outgoing_fwd.insert(id);
                } else {
                    result.outgoing_rev.insert(id);
                }
            }
        }
    }

    // 3. Contained: nodes that appear in both incoming and outgoing (any
    // orientation).
    auto all_in = result.all_incoming();
    auto all_out = result.all_outgoing();
    for (u64 id : all_in) {
        if (all_out.contains(id)) {
            result.contained.push_back(id);
        }
    }

    return result;
}

namespace layout_graph_algo {
using layout_result_t =
    gtl::vector<gtl::flat_hash_map<u64, std::tuple<int, int, bool>>>;

template <class RecigMapT>
inline layout_result_t connected_comp_plus_layout(const RecigMapT &recig) {
    layout_result_t layouts;
    std::queue<std::pair<u64, typename RecigMapT::mapped_type>> q;
    gtl::flat_hash_set<u64> visited;

    auto infer_len_for_node = [&](u64 id) -> int {
        auto out_range = recig.equal_range(id);
        if (out_range.first != out_range.second) {
            const auto &rm = out_range.first->second;
            if (rm.qry.size > 0) return static_cast<int>(rm.qry.size);
        }
        for (const auto &kv : recig) {
            if (kv.second.node != id) continue;
            const auto &rm = kv.second;
            if (rm.tgt.size > 0) return static_cast<int>(rm.tgt.size);
        }
        return 0;
    };

    auto to_end_inclusive = [&](int start, int len) -> int {
        if (len <= 0) return start;
        if (start > std::numeric_limits<int>::max() - (len - 1)) {
            return std::numeric_limits<int>::max();
        }
        if (start < std::numeric_limits<int>::min() + (len - 1)) {
            return std::numeric_limits<int>::min();
        }
        return start + (len - 1);
    };

    for (const auto &nit : recig) {
        if (visited.contains(nit.first)) continue;
        layouts.emplace_back();
        const int root_start = 0;
        const int root_len = infer_len_for_node(nit.first);
        const int root_end = to_end_inclusive(root_start, root_len);
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

            long long qs = s_orient ? (node.second.qry.size - node.second.qry.to)
                                    : node.second.qry.from;
            long long ts =
                target_orient ? (node.second.tgt.size - node.second.tgt.to)
                              : node.second.tgt.from;

            int delta = static_cast<int>(qs - ts);
            int target_start = spos + delta;
            const int target_len = infer_len_for_node(tgt);
            const int target_end = to_end_inclusive(target_start, target_len);

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

template <class RecigMapT>
inline layout_result_t connected_comp_plus_layout_lstsq(
    const RecigMapT &recig) {
    layout_result_t layouts;
    gtl::flat_hash_set<u64> visited;

    auto infer_len_for_node = [&](u64 id) -> int {
        auto out_range = recig.equal_range(id);
        if (out_range.first != out_range.second) {
            const auto &rm = out_range.first->second;
            if (rm.qry.size > 0) return static_cast<int>(rm.qry.size);
        }
        for (const auto &kv : recig) {
            if (kv.second.node != id) continue;
            const auto &rm = kv.second;
            if (rm.tgt.size > 0) return static_cast<int>(rm.tgt.size);
        }
        return 0;
    };

    auto to_end_inclusive = [&](int start, int len) -> int {
        if (len <= 0) return start;
        if (start > std::numeric_limits<int>::max() - (len - 1)) {
            return std::numeric_limits<int>::max();
        }
        if (start < std::numeric_limits<int>::min() + (len - 1)) {
            return std::numeric_limits<int>::min();
        }
        return start + (len - 1);
    };

    for (const auto &nit : recig) {
        if (visited.contains(nit.first)) continue;

        gtl::vector<u64> component;
        gtl::flat_hash_set<u64> comp_visited;
        std::stack<u64> stack;
        stack.push(nit.first);

        while (!stack.empty()) {
            u64 node = stack.top();
            stack.pop();
            if (comp_visited.contains(node)) continue;
            comp_visited.emplace(node);
            visited.emplace(node);
            component.push_back(node);

            auto range = recig.equal_range(node);
            for (auto it = range.first; it != range.second; ++it) {
                if (!comp_visited.contains(it->second.node)) {
                    stack.push(it->second.node);
                }
            }
        }
        if (component.empty()) continue;

        gtl::flat_hash_map<u64, size_t> idx;
        for (size_t i = 0; i < component.size(); ++i) idx[component[i]] = i;
        const size_t n = component.size();

        struct edge_info {
            u64 src;
            u64 tgt;
            int delta;
        };
        gtl::vector<edge_info> edges;

        gtl::flat_hash_map<u64, bool> orientations;
        std::queue<u64> orient_queue;
        orient_queue.push(component[0]);
        orientations[component[0]] = true;

        while (!orient_queue.empty()) {
            u64 src = orient_queue.front();
            orient_queue.pop();
            bool s_orient = orientations[src];

            auto range = recig.equal_range(src);
            for (auto it = range.first; it != range.second; ++it) {
                u64 tgt = it->second.node;
                if (!idx.contains(tgt)) continue;

                bool tgt_orient;
                if (orientations.contains(tgt)) {
                    tgt_orient = orientations[tgt];
                } else {
                    tgt_orient = !s_orient;
                    orientations[tgt] = tgt_orient;
                    orient_queue.push(tgt);
                }

                long long qs =
                    s_orient ? (it->second.qry.size - it->second.qry.to)
                             : it->second.qry.from;
                long long ts =
                    tgt_orient ? (it->second.tgt.size - it->second.tgt.to)
                               : it->second.tgt.from;
                edges.push_back({src, tgt, static_cast<int>(qs - ts)});
            }
        }

        Eigen::MatrixXd L = Eigen::MatrixXd::Zero(n, n);
        Eigen::VectorXd b = Eigen::VectorXd::Zero(n);
        for (const auto &e : edges) {
            const size_t i = idx[e.src];
            const size_t j = idx[e.tgt];
            L(i, i) += 1.0;
            L(j, j) += 1.0;
            L(i, j) -= 1.0;
            L(j, i) -= 1.0;
            b(j) += e.delta;
            b(i) -= e.delta;
        }
        L.row(0).setZero();
        L(0, 0) = 1.0;
        b(0) = 0.0;

        const Eigen::VectorXd x_vec = L.colPivHouseholderQr().solve(b);
        gtl::flat_hash_map<u64, double> x_positions;
        for (size_t i = 0; i < n; ++i) x_positions[component[i]] = x_vec(i);

        std::sort(component.begin(), component.end(), [&](u64 a, u64 bnode) {
            return x_positions[a] < x_positions[bnode];
        });

        layouts.emplace_back();
        const double min_pos = x_positions[component[0]];
        for (u64 node_id : component) {
            const double pos = x_positions[node_id] - min_pos;
            const int start = static_cast<int>(std::round(pos));
            const int len = infer_len_for_node(node_id);
            const int end = to_end_inclusive(start, len);
            const bool orient =
                orientations.contains(node_id) ? orientations[node_id] : true;
            layouts.back().emplace(node_id, std::make_tuple(start, end, orient));
        }
    }
    return layouts;
}
}  // namespace layout_graph_algo

namespace layout_graph_emit {
struct ultra_layout_node_info {
    u64 id;
    int from;
    int to;
    bool ori;
};

struct ultra_layout_job {
    u64 layout_idx;
    gtl::vector<ultra_layout_node_info> lnodes;
};

// append_seqvec_ascii removed — ASCII conversion deferred to print time

template <class AssemblerT>
inline gtl::vector<SequenceContainerType> build_layout_sequence_from_nodes(
    [[maybe_unused]] AssemblerT &as, dbg_index &dbg,
    const gtl::vector<ultra_layout_node_info> &lnodes,
    [[maybe_unused]] int hash_combine_len, u64 layout_idx = 0) {

    // ── Local structs ──────────────────────────────────────────────────
    struct lane_node {
        u64 id = 0;
        int from = 0;
        int to = 0;
        bool is_forward_lane = true;
        const gtl::vector<SequenceContainerType> *raw_seqs = nullptr;
    };

    if (lnodes.empty()) {
        return {};
    }

    // ── Phase 0: Build lane_nodes ──────────────────────────────────────
    gtl::vector<lane_node> nodes;
    nodes.reserve(lnodes.size());

    for (const auto &n : lnodes) {
        auto it = dbg.dbg.find(n.id);
        if (it == dbg.dbg.end()) continue;
        if (it->second.seq.empty()) continue;

        lane_node ln;
        ln.id = n.id;
        ln.from = n.from;
        ln.to = n.to;
        ln.is_forward_lane = n.ori;
        ln.raw_seqs = &it->second.seq;
        nodes.push_back(std::move(ln));
    }
    if (nodes.empty()) return {};

    // Sort all contigs by layout position
    gtl::vector<size_t> sorted_idx(nodes.size());
    for (size_t i = 0; i < nodes.size(); ++i) sorted_idx[i] = i;
    std::sort(sorted_idx.begin(), sorted_idx.end(),
              [&](size_t ia, size_t ib) {
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
            // Pick the contig with more bases
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

    // ── Per-contig core sequences ─────────────────────────────────────
    gtl::flat_hash_map<size_t, gtl::vector<SequenceContainerType>> contig_core_seqs;
    gtl::flat_hash_map<size_t, lcp_read> contig_lr;

    // Build cores for all contigs in layout order, RC only if reverse
    gtl::vector<size_t> contig_order;
    contig_order.reserve(sorted_idx.size());

    for (size_t ci : sorted_idx) {
        const auto &raw = *nodes[ci].raw_seqs;

        // RC fragments if reverse orientation
        gtl::vector<SequenceContainerType> rev_seqs;
        const gtl::vector<SequenceContainerType> *seqs_ptr = &raw;
        if (!nodes[ci].is_forward_lane) {
            rev_seqs.reserve(raw.size());
            for (auto rit = raw.rbegin(); rit != raw.rend(); ++rit)
                rev_seqs.emplace_back(rit->reverse_complement());
            seqs_ptr = &rev_seqs;
        }

        auto [lr, cseqs] = AssemblerT::build_lr_cores_with_seqs(
            *seqs_ptr, as.syncmer_functor, as.lcp_level);

        if (lr.cores.empty()) {
            // Synthetic single-core fallback — concatenate fragments
            // directly without ASCII round-trip
            const auto &frags = *seqs_ptr;
            size_t total_len = 0;
            for (const auto &f : frags) total_len += f.len;
            if (total_len == 0) continue;

            SequenceContainerType synth(frags.begin(), frags.end());
            gtl::vector<SequenceContainerType> synth_vec;
            synth_vec.push_back(std::move(synth));
            lcp_read synth_lr;
            synth_lr.cores.push_back(0);
            synth_lr.pos.emplace_back(0u,
                static_cast<uint32_t>(total_len));
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
        std::fprintf(stderr,
                     "[layout-dag] layout=%llu SKIPPED: no contigs with cores\n",
                     (unsigned long long)layout_idx);
        return {};
    }

    // ── POA graph ──────────────────────────────────────────────────────
    struct poa_node {
        uint32_t core_label;
        size_t   contig_idx;  // which contig's core_seq to emit
        size_t   core_idx;
        uint16_t weight;      // how many contigs pass through
    };

    gtl::vector<poa_node> poa;
    // Adjacency: node -> [(successor, edge_weight)]
    gtl::vector<gtl::vector<std::pair<uint32_t, uint16_t>>> poa_fwd;
    // Label index: core_label -> [node_ids]
    gtl::flat_hash_map<uint32_t, gtl::vector<uint32_t>> label_nodes;
    // All (contig_idx, core_idx) that pass through each node
    gtl::vector<gtl::vector<std::pair<size_t, size_t>>> node_contigs;

    poa.reserve(contig_order.size() * 64);
    poa_fwd.reserve(contig_order.size() * 64);
    node_contigs.reserve(contig_order.size() * 64);

    // ── Phase 1: Add each contig to the POA graph via DP alignment ─────
    // Standard POA: align each new contig's core-label sequence against
    // the existing DAG using DP (LCS on DAG), then merge matched nodes
    // and insert new nodes for unmatched query positions.
    auto add_edge = [&](int32_t from, uint32_t to) {
        if (from < 0) return;
        for (auto &[t, ew] : poa_fwd[static_cast<size_t>(from)]) {
            if (t == to) { ++ew; return; }
        }
        poa_fwd[static_cast<size_t>(from)].push_back({to, 1});
    };

    for (size_t ci : contig_order) {
        const auto &lr = contig_lr[ci];
        const size_t M = lr.cores.size();
        if (M == 0) continue;

        const size_t gn = poa.size();

        // First contig or empty graph: just create nodes
        if (gn == 0) {
            int32_t last_node = -1;
            for (size_t k = 0; k < M; ++k) {
                uint32_t cur = static_cast<uint32_t>(poa.size());
                poa.push_back({lr.cores[k], ci, k, 1});
                poa_fwd.emplace_back();
                node_contigs.push_back({{ci, k}});
                label_nodes[lr.cores[k]].push_back(cur);
                add_edge(last_node, cur);
                last_node = static_cast<int32_t>(cur);
            }
            continue;
        }

        // ── Seed-chain alignment ────────────────────────────────────
        // 1. Compute topological ranks — component by component so
        //    disconnected components get non-interleaving rank ranges.
        //    Components ordered by min layout_from so earlier-in-layout
        //    components get lower ranks.
        gtl::vector<uint32_t> topo_rank(gn, 0);
        {
            gtl::vector<uint32_t> in_deg(gn, 0);
            for (size_t u = 0; u < gn; ++u)
                for (auto &[v, w] : poa_fwd[u])
                    if (v < gn) ++in_deg[v];

            // Union-find for connected components
            gtl::vector<uint32_t> uf(gn);
            for (uint32_t i = 0; i < gn; ++i) uf[i] = i;
            auto uf_find = [&](uint32_t x) -> uint32_t {
                while (uf[x] != x) { uf[x] = uf[uf[x]]; x = uf[x]; }
                return x;
            };
            for (size_t u = 0; u < gn; ++u)
                for (auto &[v, w] : poa_fwd[u])
                    if (v < gn) uf[uf_find(u)] = uf_find(v);

            // Group sources by component root
            gtl::flat_hash_map<uint32_t, gtl::vector<uint32_t>> comp_srcs;
            for (uint32_t i = 0; i < gn; ++i)
                if (in_deg[i] == 0) comp_srcs[uf_find(i)].push_back(i);

            // Sort component roots by min layout_from
            gtl::vector<uint32_t> roots;
            roots.reserve(comp_srcs.size());
            for (auto &[r, _] : comp_srcs) roots.push_back(r);
            std::sort(roots.begin(), roots.end(),
                      [&](uint32_t a, uint32_t b) {
                          int ma = INT_MAX, mb = INT_MAX;
                          for (uint32_t s : comp_srcs[a])
                              ma = std::min(ma,
                                   nodes[poa[s].contig_idx].from);
                          for (uint32_t s : comp_srcs[b])
                              mb = std::min(mb,
                                   nodes[poa[s].contig_idx].from);
                          return ma < mb;
                      });

            // Topo sort each component sequentially
            uint32_t rank = 0;
            for (uint32_t root : roots) {
                std::queue<uint32_t> q;
                for (uint32_t s : comp_srcs[root]) q.push(s);
                while (!q.empty()) {
                    uint32_t u = q.front(); q.pop();
                    topo_rank[u] = rank++;
                    for (auto &[v, w] : poa_fwd[u])
                        if (v < gn && --in_deg[v] == 0) q.push(v);
                }
            }
        }

        // 2. Collect seeds: (rank, query_pos, graph_node)
        struct poa_seed {
            uint32_t rank; size_t k; uint32_t v;
        };
        gtl::vector<poa_seed> seeds;
        for (size_t k = 0; k < M; ++k) {
            auto it = label_nodes.find(lr.cores[k]);
            if (it == label_nodes.end()) continue;
            for (uint32_t v : it->second) {
                if (v >= gn) continue;
                seeds.push_back({topo_rank[v], k, v});
            }
        }

        // 3. Sort by rank asc, k desc (standard 2D strictly-increasing
        //    LIS trick: ties in rank get k-descending so at most one
        //    per rank survives in any chain).
        std::sort(seeds.begin(), seeds.end(),
                  [](const poa_seed &a, const poa_seed &b) {
                      if (a.rank != b.rank) return a.rank < b.rank;
                      return a.k > b.k;
                  });

        // 4. LIS on k (strictly increasing) with traceback
        const size_t S = seeds.size();
        gtl::vector<size_t> tails;       // tails[i] = k of best len-(i+1) chain
        gtl::vector<size_t> tail_idx;    // seed index for that tail
        gtl::vector<int32_t> lis_prev(S, -1);

        for (size_t i = 0; i < S; ++i) {
            size_t k = seeds[i].k;
            auto pos = std::lower_bound(tails.begin(), tails.end(), k);
            size_t lp = static_cast<size_t>(pos - tails.begin());

            lis_prev[i] = (lp > 0)
                ? static_cast<int32_t>(tail_idx[lp - 1]) : -1;

            if (pos == tails.end()) {
                tails.push_back(k);
                tail_idx.push_back(i);
            } else {
                *pos = k;
                tail_idx[lp] = i;
            }
        }

        // 5. Traceback → query_to_node
        gtl::vector<int32_t> query_to_node(M, -1);
        if (!tails.empty()) {
            for (int32_t idx = static_cast<int32_t>(tail_idx.back());
                 idx >= 0; idx = lis_prev[idx])
                query_to_node[seeds[idx].k] =
                    static_cast<int32_t>(seeds[idx].v);
        }

        // 6. Apply: merge matched, create new for unmatched
        int32_t last_node = -1;
        for (size_t k = 0; k < M; ++k) {
            uint32_t cur_node;
            if (query_to_node[k] >= 0) {
                cur_node = static_cast<uint32_t>(query_to_node[k]);
                poa[cur_node].weight++;
                node_contigs[cur_node].push_back({ci, k});
            } else {
                cur_node = static_cast<uint32_t>(poa.size());
                poa.push_back({lr.cores[k], ci, k, 1});
                poa_fwd.emplace_back();
                node_contigs.push_back({{ci, k}});
                label_nodes[lr.cores[k]].push_back(cur_node);
            }
            add_edge(last_node, cur_node);
            last_node = static_cast<int32_t>(cur_node);
        }
    }

    if (poa.empty()) {
        std::fprintf(stderr,
                     "[layout-dag] layout=%llu SKIPPED: POA graph empty\n",
                     (unsigned long long)layout_idx);
        return {};
    }

    // ── Phase 2: Find heaviest path per connected component ────────────
    const size_t n = poa.size();

    // Final topological sort
    gtl::vector<uint32_t> topo_order;
    topo_order.reserve(n);
    {
        gtl::vector<uint32_t> in_deg(n, 0);
        for (size_t u = 0; u < n; ++u)
            for (auto &[v, w] : poa_fwd[u])
                ++in_deg[v];
        std::queue<uint32_t> q;
        for (uint32_t i = 0; i < n; ++i)
            if (in_deg[i] == 0) q.push(i);
        while (!q.empty()) {
            uint32_t u = q.front(); q.pop();
            topo_order.push_back(u);
            for (auto &[v, w] : poa_fwd[u])
                if (--in_deg[v] == 0) q.push(v);
        }
        // If some nodes weren't reached (residual cycles), append them
        if (topo_order.size() < n) {
            std::fprintf(stderr,
                         "[layout-dag] layout=%llu WARNING: %zu/%zu nodes "
                         "stuck in cycles, appending in index order\n",
                         (unsigned long long)layout_idx,
                         n - topo_order.size(), n);
            gtl::flat_hash_set<uint32_t> visited(
                topo_order.begin(), topo_order.end());
            for (uint32_t i = 0; i < n; ++i)
                if (!visited.contains(i)) topo_order.push_back(i);
        }
    }

    // Build reverse adjacency
    gtl::vector<gtl::vector<std::pair<uint32_t, uint16_t>>> poa_rev(n);
    for (size_t u = 0; u < n; ++u)
        for (auto &[v, w] : poa_fwd[u])
            poa_rev[v].push_back({static_cast<uint32_t>(u), w});

    // DP: heaviest path (maximize sum of node weights)
    gtl::vector<int64_t> dp(n, 0);
    gtl::vector<int32_t> parent(n, -1);
    for (uint32_t v : topo_order) {
        dp[v] = poa[v].weight;
        parent[v] = -1;
        for (auto &[u, w] : poa_rev[v]) {
            int64_t cand = dp[u] + poa[v].weight;
            if (cand > dp[v]) {
                dp[v] = cand;
                parent[v] = static_cast<int32_t>(u);
            }
        }
    }

    // Find connected components (undirected BFS)
    gtl::vector<int32_t> comp_id(n, -1);
    int32_t num_comps = 0;
    for (size_t i = 0; i < n; ++i) {
        if (comp_id[i] >= 0) continue;
        std::queue<uint32_t> bfs;
        bfs.push(static_cast<uint32_t>(i));
        comp_id[i] = num_comps;
        while (!bfs.empty()) {
            uint32_t u = bfs.front(); bfs.pop();
            for (auto &[v, w] : poa_fwd[u]) {
                if (comp_id[v] < 0) { comp_id[v] = num_comps; bfs.push(v); }
            }
            for (auto &[u2, w] : poa_rev[u]) {
                if (comp_id[u2] < 0) { comp_id[u2] = num_comps; bfs.push(u2); }
            }
        }
        ++num_comps;
    }

    // For each component: find heaviest path, determine layout position
    struct comp_path {
        int layout_from;
        int layout_to;
        int64_t total_weight;
        gtl::vector<uint32_t> path;
    };
    gtl::vector<comp_path> comp_paths;

    for (int32_t c = 0; c < num_comps; ++c) {
        // Find best sink (no outgoing edges) in this component
        int32_t best_end = -1;
        int64_t best_dp = -1;
        for (size_t i = 0; i < n; ++i) {
            if (comp_id[i] != c) continue;
            if (poa_fwd[i].empty() && dp[i] > best_dp) {
                best_dp = dp[i];
                best_end = static_cast<int32_t>(i);
            }
        }
        // Fallback: max dp in component if no sinks
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

        // Trace back
        gtl::vector<uint32_t> path;
        for (int32_t cur = best_end; cur >= 0; cur = parent[cur])
            path.push_back(static_cast<uint32_t>(cur));
        std::reverse(path.begin(), path.end());

        // Extend the heaviest path at both ends if there are
        // sinks/sources in this component that cover more layout span.
        {
            gtl::flat_hash_set<uint32_t> on_path(path.begin(), path.end());
            int path_max_to = INT_MIN;
            for (uint32_t nid : path)
                path_max_to = std::max(path_max_to,
                                       nodes[poa[nid].contig_idx].to);
            int path_min_from = INT_MAX;
            for (uint32_t nid : path)
                path_min_from = std::min(path_min_from,
                                         nodes[poa[nid].contig_idx].from);

            // Forward extension: find sink with layout_to > path_max_to
            int32_t ext_sink = -1;
            int ext_to = path_max_to;
            for (size_t i = 0; i < n; ++i) {
                if (comp_id[i] != c) continue;
                if (on_path.contains(static_cast<uint32_t>(i))) continue;
                int to = nodes[poa[i].contig_idx].to;
                if (to > ext_to && poa_fwd[i].empty()) {
                    ext_to = to;
                    ext_sink = static_cast<int32_t>(i);
                }
            }
            if (ext_sink >= 0) {
                // Traceback from ext_sink, find where it meets the main path
                gtl::vector<uint32_t> ext_path;
                for (int32_t cur = ext_sink; cur >= 0; cur = parent[cur])
                    ext_path.push_back(static_cast<uint32_t>(cur));
                std::reverse(ext_path.begin(), ext_path.end());
                // Find last node in ext_path that's on the main path
                size_t branch_pos = 0;
                for (size_t i = 0; i < ext_path.size(); ++i)
                    if (on_path.contains(ext_path[i])) branch_pos = i;
                // Append the diverging suffix
                for (size_t i = branch_pos + 1; i < ext_path.size(); ++i)
                    path.push_back(ext_path[i]);
            }

            // Backward extension: find source with layout_from < path_min_from
            int32_t ext_source = -1;
            int ext_from = path_min_from;
            for (size_t i = 0; i < n; ++i) {
                if (comp_id[i] != c) continue;
                if (on_path.contains(static_cast<uint32_t>(i))) continue;
                int from = nodes[poa[i].contig_idx].from;
                if (from < ext_from && poa_rev[i].empty()) {
                    ext_from = from;
                    ext_source = static_cast<int32_t>(i);
                }
            }
            if (ext_source >= 0) {
                // Follow successors greedily from ext_source toward the
                // main path, guided by dp values.
                gtl::vector<uint32_t> prefix;
                uint32_t cur = static_cast<uint32_t>(ext_source);
                while (!on_path.contains(cur)) {
                    prefix.push_back(cur);
                    uint32_t best_succ = UINT32_MAX;
                    int64_t best_dp_val = -1;
                    for (auto &[v, w] : poa_fwd[cur]) {
                        if (dp[v] > best_dp_val) {
                            best_dp_val = dp[v];
                            best_succ = v;
                        }
                    }
                    if (best_succ == UINT32_MAX) break;
                    cur = best_succ;
                }
                if (!prefix.empty()) {
                    prefix.insert(prefix.end(), path.begin(), path.end());
                    path = std::move(prefix);
                }
            }
        }

        // Layout position range from all nodes' contigs in this component
        int min_pos = INT_MAX, max_pos = INT_MIN;
        for (size_t i = 0; i < n; ++i) {
            if (comp_id[i] != c) continue;
            min_pos = std::min(min_pos, nodes[poa[i].contig_idx].from);
            max_pos = std::max(max_pos, nodes[poa[i].contig_idx].to);
        }

        comp_paths.push_back({min_pos, max_pos, best_dp, std::move(path)});
    }

    // Sort components by layout position
    std::sort(comp_paths.begin(), comp_paths.end(),
              [](const auto &a, const auto &b) {
                  return a.layout_from < b.layout_from;
              });


    // ── Phase 3: Collect output sequences ──────────────────────────────
    // For each node on the path, pick a consistent contig to avoid
    // boundary mismatches at merged nodes. Collect CrumbledSubSeq64
    // objects directly — defer ASCII conversion to print time.
    gtl::vector<SequenceContainerType> merged_seqs;
    merged_seqs.reserve(n);

    for (size_t ci = 0; ci < comp_paths.size(); ++ci) {
        if (ci > 0 && !merged_seqs.empty()) {
            // Sentinel: zero-len entry signals N-gap between components
            merged_seqs.emplace_back(0, 0);
        }

        const auto &path = comp_paths[ci].path;
        size_t emit_contig = SIZE_MAX;

        for (size_t pi = 0; pi < path.size(); ++pi) {
            uint32_t nid = path[pi];
            const auto &nc = node_contigs[nid];

            // Try to stay on current contig
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

            // If current contig doesn't pass through this node, switch.
            // Prefer the contig that also appears in the NEXT node.
            if (use_contig == SIZE_MAX) {
                if (pi + 1 < path.size()) {
                    const auto &nc_next = node_contigs[path[pi + 1]];
                    for (auto &[cid, kid] : nc) {
                        for (auto &[cid2, kid2] : nc_next) {
                            if (cid == cid2) {
                                use_contig = cid;
                                use_core = kid;
                                goto found_contig;
                            }
                        }
                    }
                }
                if (!nc.empty()) {
                    use_contig = nc[0].first;
                    use_core = nc[0].second;
                }
            }
            found_contig:

            if (use_contig == SIZE_MAX) continue;
            emit_contig = use_contig;

            auto cs_it = contig_core_seqs.find(use_contig);
            if (cs_it != contig_core_seqs.end() &&
                use_core < cs_it->second.size()) {
                auto &core_seq = cs_it->second[use_core];
                if (core_seq.len > 0) {
                    merged_seqs.push_back(std::move(core_seq));
                }
            }
        }
    }

    return merged_seqs;
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
                if (!job_queue.wait_dequeue_timed(job_idx, WAIT_QUEUE_TIMEOUT)) {
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

template <class AssemblerT>
inline void emit_ultra_merged_graph(AssemblerT &as, dbg_index &dbg, FILE *out) {
    INIT_TIMERS(LayoutEmit);
    auto [recig, layouts] = as.fold_dbg_graph(dbg);
    TIME_CHECKPOINT(LayoutEmit, stderr, "fold_dbg_graph (%zu layouts)", layouts.size());

    static constexpr uint8_t SIDE_START = 1;
    static constexpr uint8_t SIDE_END = 2;
    auto side_label = [&](uint8_t side) -> const char * {
        if (side == SIDE_START) return "start";
        if (side == SIDE_END) return "end";
        return "na";
    };
    auto ori_char = [](bool ori) -> char { return ori ? '+' : '-'; };

    struct supernode_out {
        std::string name;
        gtl::vector<SequenceContainerType> owned_seqs;   // layout seqs (owned)
        const gtl::vector<SequenceContainerType> *ref_seqs = nullptr; // unitig seqs (borrowed)
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
        int min_start = INT_MAX;
        int max_end = INT_MIN;
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

        const int span = std::max(1, max_end - min_start);
        const int fuzzy = std::max<int>(as.olc_lcmer_l * 2, span / 20);
        for (const auto &n : lnodes) {
            uint8_t mask = 0;
            if (n.from <= min_start + fuzzy) mask |= SIDE_START;
            if (n.to >= max_end - fuzzy) mask |= SIDE_END;
            if (mask == 0) {
                int d_start = n.from - min_start;
                if (d_start < 0) d_start = -d_start;
                int d_end = max_end - n.to;
                if (d_end < 0) d_end = -d_end;
                mask = (d_start <= d_end) ? SIDE_START : SIDE_END;
            }
            node_boundary_mask[n.id] |= mask;
            node_layout_ori[n.id] = n.ori;
            nodes_in_layouts.insert(n.id);
        }

        layout_jobs.push_back({layout_idx, std::move(lnodes)});
        ++layout_idx;
    }


    auto layout_seqs =
        build_layout_sequences_parallel(as, dbg, layout_jobs, hash_combine_len);
    TIME_CHECKPOINT(LayoutEmit, stderr, "Build layout sequences (parallel)");
    for (size_t i = 0; i < layout_jobs.size(); ++i) {
        const auto &job = layout_jobs[i];
        const u64 sid = supernodes.size();
        supernodes.push_back({"M" + std::to_string(job.layout_idx),
                              std::move(layout_seqs[i]),
                              nullptr, true, job.layout_idx});
        for (const auto &n : job.lnodes) {
            node_to_super[n.id] = sid;
        }
    }

    for (u64 node_id : dbg.active_set) {
        if (nodes_in_layouts.contains(node_id)) continue;
        auto it = dbg.dbg.find(node_id);
        if (it == dbg.dbg.end()) continue;
        if (it->second.redirected() || it->second.seq.empty()) continue;

        const u64 sid = supernodes.size();
        supernodes.push_back(
            {"C" + std::to_string(node_id),
             {}, &it->second.seq, false, node_id});
        node_to_super[node_id] = sid;
    }

    TIME_CHECKPOINT(LayoutEmit, stderr, "Supernode creation");
    using edge_key = std::tuple<u64, u64, char, char, uint8_t, uint8_t>;
    std::map<edge_key, edge_agg> merged_edges;
    for (const auto &[to_id, to_node] : dbg.dbg) {
        if (to_node.redirected()) continue;
        auto to_super_it = node_to_super.find(to_id);
        if (to_super_it == node_to_super.end()) continue;

        for (u64 from_id : to_node.prev) {
            auto from_super_it = node_to_super.find(from_id);
            if (from_super_it == node_to_super.end()) continue;
            if (from_super_it->second == to_super_it->second) continue;

            bool accept = true;
            uint8_t src_side = 0, dst_side = 0;
            char src_ori = '+', dst_ori = '+';

            if (nodes_in_layouts.contains(from_id)) {
                bool ori = node_layout_ori[from_id];
                src_ori = ori_char(ori);
                src_side = ori ? SIDE_END : SIDE_START;
                if ((node_boundary_mask[from_id] & src_side) == 0) {
                    accept = false;
                }
            }
            if (nodes_in_layouts.contains(to_id)) {
                bool ori = node_layout_ori[to_id];
                dst_ori = ori_char(ori);
                dst_side = ori ? SIDE_START : SIDE_END;
                if ((node_boundary_mask[to_id] & dst_side) == 0) {
                    accept = false;
                }
            }
            if (!accept) continue;

            edge_key k{from_super_it->second,
                       to_super_it->second,
                       src_ori,
                       dst_ori,
                       src_side,
                       dst_side};
            auto &agg = merged_edges[k];
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
    }

    TIME_CHECKPOINT(LayoutEmit, stderr, "Edge aggregation");
    // Print supernodes — convert CrumbledSubSeq64 to ASCII only now
    constexpr size_t N_GAP_LEN = 64;
    charsi chs;
    for (const auto &sn : supernodes) {
        const char *tp = sn.is_layout ? "layout" : "unitig";
        const char *id_tag = sn.is_layout ? "LY" : "UI";
        std::fprintf(out, "S\t%s\t", sn.name.c_str());
        bool any_printed = false;
        for (const auto &s : sn.seqs()) {
            if (s.len == 0) {
                // N-gap sentinel
                for (size_t i = 0; i < N_GAP_LEN; ++i)
                    std::fputc('N', out);
                any_printed = true;
                continue;
            }
            chs.require(s.len + 32);
            s.convert_to_ascii((char *)chs.data);
            std::fwrite(chs.data, 1, s.len, out);
            any_printed = true;
        }
        if (!any_printed) std::fputc('*', out);
        std::fprintf(out, "\tTP:Z:%s\t%s:i:%llu\n",
                     tp, id_tag, (unsigned long long)sn.aux_id);
    }
    for (const auto &[_, edge] : merged_edges) {
        const auto &src_name = supernodes[edge.src_super].name;
        const auto &dst_name = supernodes[edge.dst_super].name;
        std::fprintf(
            out,
            "L\t%s\t%c\t%s\t%c\t0M\tSO:Z:%c\tTO:Z:%c\tSS:Z:%s\tTS:Z:%s\t"
            "EC:i:%zu\tFN:i:%llu\tTN:i:%llu\n",
            src_name.c_str(), edge.src_ori, dst_name.c_str(), edge.dst_ori,
            edge.src_ori, edge.dst_ori, side_label(edge.src_side),
            side_label(edge.dst_side), edge.count,
            (unsigned long long)edge.example_from,
            (unsigned long long)edge.example_to);
    }
    TIME_CHECKPOINT(LayoutEmit, stderr, "GFA printing");
}

template <class AssemblerT>
inline void emit_ultra_merged_graph(AssemblerT &as, dbg_index &dbg) {
    emit_ultra_merged_graph(as, dbg, stdout);
}

}  // namespace layout_graph_emit

// ── Unused legacy code (commented out) ──────────────────────────────
#if 0
inline void print_merged_layouts(
    FILE *out, layout_graph_algo::layout_result_t layouts,
    dbg_index::dbg_map_t &dbg) {
    charsi buf;
    u64 contig_id = 0;

    for (auto &layout : layouts) {
        if (layout.empty()) continue;

        struct LayoutEntry {
            int start;
            int end;
            u64 node_id;
            bool orient;
            dbg_index::dbg_map_t::iterator dbg_it;
        };
        gtl::vector<LayoutEntry> entries;
        entries.reserve(layout.size());

        for (auto &[node_id, start_end_orient] : layout) {
            auto [pos_core, end_core, orient] = start_end_orient;
            auto it = dbg.find(node_id);
            if (it == dbg.end() || it->second.seq.empty()) continue;
            entries.push_back({pos_core, end_core, node_id, orient, it});
        }
        if (entries.empty()) continue;

        std::sort(entries.begin(), entries.end(),
                  [](const auto &a, const auto &b) {
                      return a.start < b.start;
                  });

        double core_to_base_scale = 0.0;
        {
            int min_start = INT_MAX;
            int max_start = INT_MIN;
            size_t total_bases = 0;
            for (const auto &e : entries) {
                min_start = std::min(min_start, e.start);
                max_start = std::max(max_start, e.start);
                size_t seq_len = 0;
                for (const auto &s : e.dbg_it->second.seq) seq_len += s.len;
                total_bases += seq_len;
            }
            const int core_span =
                (max_start >= min_start) ? (max_start - min_start) : 0;
            if (core_span > 0 && total_bases > 0) {
                core_to_base_scale = static_cast<double>(total_bases) /
                                     static_cast<double>(core_span);
            }
        }

        gtl::vector<LayoutEntry> non_contained;
        non_contained.reserve(entries.size());
        int max_end = INT_MIN;
        for (auto &entry : entries) {
            if (entry.end <= max_end) {
                entry.dbg_it->second.seq.clear();
                entry.dbg_it->second.prev.clear();
                continue;
            }
            non_contained.push_back(std::move(entry));
            max_end = std::max(max_end, non_contained.back().end);
        }
        if (non_contained.empty()) continue;

        std::fprintf(out, "S\tM%llu\t", (unsigned long long)contig_id);
        size_t total_len = 0;
        int current_end = non_contained[0].start;
        gtl::vector<u64> prev_edges = non_contained[0].dbg_it->second.prev;

        for (size_t i = 0; i < non_contained.size(); ++i) {
            auto &curr = non_contained[i];
            int overlap_cores = current_end - curr.start;
            if (overlap_cores < 0) overlap_cores = 0;

            size_t bases_to_skip = 0;
            if (core_to_base_scale > 0.0) {
                const double est =
                    static_cast<double>(overlap_cores) * core_to_base_scale;
                if (est > 0.0) bases_to_skip = static_cast<size_t>(est);
            }

            gtl::vector<SequenceContainerType> curr_seqs;
            if (curr.orient) {
                for (auto &s : curr.dbg_it->second.seq) {
                    curr_seqs.emplace_back(std::move(s));
                }
            } else {
                for (auto it = curr.dbg_it->second.seq.rbegin();
                     it != curr.dbg_it->second.seq.rend(); ++it) {
                    curr_seqs.emplace_back(it->reverse_complement());
                }
            }

            curr.dbg_it->second.seq.clear();
            curr.dbg_it->second.prev.clear();

            size_t skipped = 0;
            for (auto &s : curr_seqs) {
                if (skipped >= bases_to_skip) {
                    buf.require(s.len + 32);
                    const char *buf_at = s.convert_to_ascii((char *)buf.data);
                    std::fprintf(out, "%s", buf_at);
                    total_len += s.len;
                } else if (skipped + s.len > bases_to_skip) {
                    const size_t local_skip = bases_to_skip - skipped;
                    buf.require(s.len + 32);
                    const char *buf_at = s.convert_to_ascii((char *)buf.data);
                    std::fprintf(out, "%s", buf_at + local_skip);
                    total_len += (s.len - local_skip);
                }
                skipped += s.len;
            }
            current_end = std::max(current_end, curr.end);
        }

        std::fprintf(out, "\tLN:i:%zu\n", total_len);
        for (u64 prev_id : prev_edges) {
            if (dbg.contains(prev_id)) {
                std::fprintf(out, "L\tC%llu\t+\tM%llu\t+\t0M\n",
                             (unsigned long long)prev_id,
                             (unsigned long long)contig_id);
            }
        }
        contig_id++;
    }
}

template <class RecigMapT>
inline gtl::vector<u64> path_with_good_ovlps_unweighted_bfs(
    const RecigMapT &recig, u64 source, u64 target) {
    if (source == target) return gtl::vector<u64>{source};

    std::stack<u64> q;
    gtl::flat_hash_set<u64> visited;
    gtl::flat_hash_map<u64, u64> parent;

    q.push(source);
    visited.insert(source);

    bool found = false;
    while (!q.empty() && !found) {
        u64 current = q.top();
        q.pop();

        auto range = recig.equal_range(current);
        gtl::vector<decltype(range.first)> rv;
        for (auto it = range.first; it != range.second; ++it) {
            rv.push_back(it);
        }

        std::sort(rv.begin(), rv.end(), [](const auto &a, const auto &b) {
            const auto &ma = a->second;
            const auto &mb = b->second;
            u64 overlap_a = ma.qry.to - ma.qry.from;
            u64 overlap_b = mb.qry.to - mb.qry.from;
            return overlap_a < overlap_b;
        });

        for (auto it = rv.begin(); it != rv.end(); ++it) {
            u64 neighbor = (*it)->second.node;
            if (visited.find(neighbor) != visited.end()) continue;
            visited.insert(neighbor);
            parent[neighbor] = current;

            if (neighbor == target) {
                found = true;
                break;
            }

            q.push(neighbor);
        }
    }

    if (!found) return {};

    gtl::vector<u64> path;
    for (u64 v = target;; v = parent[v]) {
        path.push_back(v);
        if (v == source) break;
    }
    std::reverse(path.begin(), path.end());
    return path;
}

template <class RecigMapT>
inline gtl::vector<u64> shortest_path_unweighted_bfs(
    const RecigMapT &recig, u64 source, u64 target) {
    if (source == target) return gtl::vector<u64>{source};

    std::queue<u64> q;
    gtl::flat_hash_set<u64> visited;
    gtl::flat_hash_map<u64, u64> parent;

    q.push(source);
    visited.insert(source);

    bool found = false;
    while (!q.empty() && !found) {
        u64 current = q.front();
        q.pop();

        auto range = recig.equal_range(current);
        for (auto it = range.first; it != range.second; ++it) {
            u64 neighbor = it->second.node;
            if (visited.find(neighbor) != visited.end()) continue;

            visited.insert(neighbor);
            parent[neighbor] = current;

            if (neighbor == target) {
                found = true;
                break;
            }

            q.push(neighbor);
        }
    }

    if (!found) return {};

    gtl::vector<u64> path;
    for (u64 v = target;; v = parent[v]) {
        path.push_back(v);
        if (v == source) break;
    }
    std::reverse(path.begin(), path.end());
    return path;
}

using half_lcmer = lcmer<LCMER_L / 2>;

struct ctg_builder {
    gtl::vector<SequenceContainerType> prefix;
    gtl::vector<SequenceContainerType> seqs;
    gtl::vector<u32> cores;
    gtl::flat_hash_map<half_lcmer, size_t> idx;

    bool empty() const noexcept { return seqs.empty(); }
};

template <class AssemblerT>
inline auto find_anchors_and_extend_the_ctg_builder(
    AssemblerT &as, ctg_builder &base,
    const gtl::vector<SequenceContainerType> &extension,
    const gtl::vector<SequenceContainerType> &next,
    const std::optional<half_lcmer> &prev_anchor, bool forward)
    -> half_lcmer {
    auto [extension_cores, reextension] = [&]() {
        if (forward) {
            return AssemblerT::build_lr_cores_with_seqs(extension,
                                                        as.syncmer_functor,
                                                        as.lcp_level);
        }
        return AssemblerT::build_rev_lr_cores_with_seqs(extension,
                                                        as.syncmer_functor,
                                                        as.lcp_level);
    }();

    auto mini_extension_index = [&]() {
        gtl::flat_hash_map<half_lcmer, size_t> idx;
        auto &cores = extension_cores.cores;
        assert(cores.size() > half_lcmer::size);
        half_lcmer l{cores.begin()};
        for (auto it = cores.begin() + half_lcmer::size; it != cores.end();
             ++it) {
            idx[l] = std::distance(cores.begin(), it);
            l.update_lcmer(*it);
        }
        return idx;
    }();

    auto left_anchor = [&]() -> half_lcmer {
        if (prev_anchor.has_value()) {
            return prev_anchor.value();
        }
        assert(base.empty());
        return half_lcmer{extension_cores.cores.begin()};
    }();

    auto next_tension_cores = [&]() {
        if (next.empty()) return lcp_read{};
        if (!forward) {
            return AssemblerT::build_lr_cores(next, as.syncmer_functor,
                                              as.lcp_level);
        }
        return AssemblerT::build_rev_lr_cores(next, as.syncmer_functor,
                                              as.lcp_level);
    }();

    auto right_anchor = [&]() -> std::optional<half_lcmer> {
        auto &cores = next_tension_cores.cores;
        if (cores.size() < half_lcmer::size) return {};

        half_lcmer l{cores.begin()};
        for (auto it = cores.rbegin() + half_lcmer::size; it != cores.rend();
             ++it) {
            if (mini_extension_index.contains(l)) return l;
            l = half_lcmer::right_shift(l);
            l.setN(*it);
        }
        return {};
    }();

    if (left_anchor == right_anchor) return right_anchor.value();

    auto left_extension_pos = [&]() {
        auto &cores = extension_cores.cores;
        if (cores.size() < half_lcmer::size) return 0L;

        half_lcmer l{cores.begin()};
        auto it = cores.begin() + half_lcmer::size;
        for (; it != cores.end(); ++it) {
            if (l == left_anchor) return std::distance(cores.begin(), it);
            l.update_lcmer(*it);
        }
        return 0L;
    }();

    assert(left_extension_pos != 0);
    auto &cores = extension_cores.cores;
    assert(left_extension_pos >= half_lcmer::size);
    assert(cores.size() > left_extension_pos);
    half_lcmer l{cores.begin() + (left_extension_pos - half_lcmer::size)};

    for (size_t i = left_extension_pos; i < cores.size(); ++i) {
        l.update_lcmer(cores[i]);
        base.seqs.emplace_back(std::move(reextension[i]));
        base.cores.push_back(cores[i]);
        base.idx.emplace(l, base.seqs.size() - 1);
        if (right_anchor.has_value() && l == right_anchor.value()) break;
    }

    return right_anchor.has_value() ? right_anchor.value() : half_lcmer{};
}

template <class AssemblerT, class RecigMapT>
inline void merge_layouts_ultra(
    AssemblerT &as, FILE *out, const layout_graph_algo::layout_result_t &layouts,
    const RecigMapT &reciprocal, dbg_index &dbgi) {
    auto identify_component_endpoints = [](const auto &layout) {
        u64 leftmost = 0;
        u64 rightmost = 0;
        int min_start = INT_MAX;
        int max_start = INT_MIN;
        int leftmost_span = INT_MIN;
        int rightmost_span = INT_MIN;
        bool start_orient_forward = true;

        for (const auto &[node_id, start_end_orient] : layout) {
            auto [start_core, end_core, orient] = start_end_orient;
            (void)orient;
            const int s = end_core - start_core;

            if (start_core < min_start ||
                (start_core == min_start && s > leftmost_span)) {
                min_start = start_core;
                leftmost = node_id;
                start_orient_forward = orient;
                leftmost_span = s;
            }

            if (start_core > max_start ||
                (start_core == max_start && s > rightmost_span)) {
                max_start = start_core;
                rightmost = node_id;
                rightmost_span = s;
            }
        }
        return std::make_tuple(leftmost, rightmost, start_orient_forward);
    };

    auto &dbg = dbgi.dbg;
    u64 layout_no = 0;
    charsi chs;
    auto single_node_print = [&out, &layout_no, &dbg,
                              &chs](const auto &layout) {
        std::fprintf(out, ">M%llu %zu\n", (unsigned long long)layout_no,
                     layout.size());
        u64 utig = layout.begin()->first;
        auto it = dbg.find(utig);
        if (it != dbg.end()) {
            for (const auto &sc : it->second.seq) {
                chs.require(sc.len);
                sc.convert_to_ascii((char *)chs.data);
                std::fprintf(out, "%s", chs.data);
            }
        }
        std::fprintf(out, "\n");
    };

    for (const auto &layout : layouts) {
        if (layout.size() <= 2) {
            single_node_print(layout);
            ++layout_no;
            continue;
        }

        auto [source, target, oriori] = identify_component_endpoints(layout);
        auto path = shortest_path_unweighted_bfs(reciprocal, source, target);
        if (path.empty()) {
            ++layout_no;
            continue;
        }
        if (path.size() == 1) {
            single_node_print(layout);
            ++layout_no;
            continue;
        }

        ctg_builder base;
        std::optional<half_lcmer> prev_anchor = std::nullopt;
        for (size_t i = 0; i + 1 < path.size(); ++i) {
            u64 left = path[i];
            u64 right = path[i + 1];

            auto left_dbg_it = dbg.find(left);
            auto right_dbg_it = dbg.find(right);
            if (left_dbg_it == dbg.end() || right_dbg_it == dbg.end()) continue;

            prev_anchor = find_anchors_and_extend_the_ctg_builder(
                as, base, left_dbg_it->second.seq, right_dbg_it->second.seq,
                prev_anchor, oriori);
            oriori = !oriori;
            std::fprintf(stderr, "Extending: %llu with %llu\n",
                         (unsigned long long)layout_no,
                         (unsigned long long)path[i]);
        }

        auto last_it = dbg.find(path.back());
        if (last_it != dbg.end()) {
            prev_anchor = find_anchors_and_extend_the_ctg_builder(
                as, base, last_it->second.seq, {}, prev_anchor,
                (path.size() - 1) % 2 == 0);
        }

        std::fprintf(out, ">L%llu\n", (unsigned long long)layout_no);
        for (const SequenceContainerType &s : base.seqs) {
            chs.require(s.len);
            s.convert_to_ascii((char *)chs.data);
            std::fprintf(out, "%s", chs.data);
        }
        std::fprintf(out, "\n");

        for (u64 node : path) {
            auto it = dbg.find(node);
            if (it == dbg.end()) continue;
            std::fprintf(out, ">Z%llu L%llu\n", (unsigned long long)node,
                         (unsigned long long)layout_no);
            for (const SequenceContainerType &s : it->second.seq) {
                chs.require(s.len);
                s.convert_to_ascii((char *)chs.data);
                std::fprintf(out, "%s", chs.data);
            }
            std::fprintf(out, "\n");
        }
        ++layout_no;
    }
}

template <class AssemblerT, class RecigMapT>
inline void merge_layouts_ldbg(
    AssemblerT &as, FILE *out, const layout_graph_algo::layout_result_t &layouts,
    const RecigMapT &reciprocal, dbg_index &dbgi) {
    auto &dbg = dbgi.dbg;
    auto identify_component_endpoints = [](const auto &layout) {
        u64 leftmost = 0;
        u64 rightmost = 0;
        int min_start = INT_MAX;
        int max_start = INT_MIN;
        int leftmost_span = INT_MIN;
        int rightmost_span = INT_MIN;

        for (const auto &[node_id, start_end_orient] : layout) {
            auto [start_core, end_core, orient] = start_end_orient;
            (void)orient;
            const int s = end_core - start_core;

            if (start_core < min_start ||
                (start_core == min_start && s > leftmost_span)) {
                min_start = start_core;
                leftmost = node_id;
                leftmost_span = s;
            }

            if (start_core > max_start ||
                (start_core == max_start && s > rightmost_span)) {
                max_start = start_core;
                rightmost = node_id;
                rightmost_span = s;
            }
        }
        return std::make_pair(leftmost, rightmost);
    };

    auto process_substr_cores =
        [](const lcp_read &lr, u64 &lcmer_index, int thread_id,
           dbg_index_light &dbgi,
           const gtl::vector<SequenceContainerType> &subseqs) {
            if (lr.cores.size() < LCMER_L) return gtl::vector<u64>{};

            lcmer<LCMER_L> l{lr.cores.begin() + 1};  // Skip first core.
            gtl::vector<uint64_t> ll(1, 0);
            ll.reserve(lr.cores.size() - LCMER_L + 3);
            ll.push_back(dbgi.emplace(l, lcmer_index, thread_id));
            for (size_t i = LCMER_L + 1; i < lr.cores.size() - 1; ++i) {
                uint32_t in_core = lr.cores[i];
                l.update_lcmer(in_core);
                ll.push_back(dbgi.emplace(l, lcmer_index, thread_id));
            }
            ll.push_back(0);

            for (size_t i = 1; i < ll.size() - 1; ++i) {
                if (AssemblerT::insert_lcmer_helper(ll, i, dbgi)) {
                    dbgi.dbg.modify_if(
                        ll[i],
                        [&lr, i, &subseqs](
                            dbg_index_light::dbg_map_t::value_type &val) {
                            val.second.seq.emplace_back(
                                subseqs, lr.pos[i + LCMER_L - 1].second,
                                lr.pos[i + LCMER_L].second);
                        });
                }
            }
            return ll;
        };

    u64 new_index = 0;
    u64 layout_no = 0;
    int tid = 0;
    charsi chs;

    for (const auto &layout : layouts) {
        dbg_index_light smal_dbgi{0.0, 0, 1UL};
        u64 src_lcm = 0;
        u64 tgt_lcm = 0;

        auto [source, target] = identify_component_endpoints(layout);
        auto path = shortest_path_unweighted_bfs(reciprocal, source, target);
        if (path.empty()) {
            ++layout_no;
            continue;
        }

        bool plus = true;
        for (size_t i = 0; i < path.size(); ++i) {
            source = path[i];
            auto source_node = dbg.find(source);
            if (source_node == dbg.end()) {
                std::fprintf(stderr, "Source node %llu not found in dbg\n",
                             (unsigned long long)source);
                continue;
            }

            if (plus) {
                auto lrs = AssemblerT::build_lr_cores(source_node->second.seq,
                                                      as.syncmer_functor,
                                                      as.lcp_level);
                auto ll = process_substr_cores(lrs, new_index, tid, smal_dbgi,
                                               source_node->second.seq);
                if (i == 0) {
                    for (u64 nod : ll) {
                        if (nod != 0) {
                            src_lcm = nod;
                            break;
                        }
                    }
                } else if (i == path.size() - 1) {
                    for (auto it = ll.rbegin(); it != ll.rend(); ++it) {
                        u64 nod = *it;
                        if (nod != 0) {
                            tgt_lcm = nod;
                            break;
                        }
                    }
                }
            } else {
                gtl::vector<SequenceContainerType> rev_seqs;
                rev_seqs.reserve(source_node->second.seq.size());
                for (auto it = source_node->second.seq.rbegin();
                     it != source_node->second.seq.rend(); ++it) {
                    rev_seqs.emplace_back(it->reverse_complement());
                }

                auto lrs = AssemblerT::build_lr_cores(rev_seqs,
                                                      as.syncmer_functor,
                                                      as.lcp_level);
                auto ll = process_substr_cores(lrs, new_index, tid, smal_dbgi,
                                               rev_seqs);
                if (i == path.size() - 1) {
                    for (auto it = ll.rbegin(); it != ll.rend(); ++it) {
                        u64 nod = *it;
                        if (nod != 0) {
                            tgt_lcm = nod;
                            break;
                        }
                    }
                }
            }
            plus = !plus;
        }

        auto sp = smal_dbgi.shortest_path_bfs(tgt_lcm, src_lcm);
        if (!sp.empty()) {
            std::fprintf(out, ">M%llu %zu\n", (unsigned long long)layout_no,
                         layout.size());
            for (auto it = sp.rbegin(); it != sp.rend(); ++it) {
                u64 n = *it;
                SequenceContainerType &seq = smal_dbgi.dbg.at(n).seq[0];
                chs.require(seq.len);
                seq.convert_to_ascii((char *)chs.data);
                std::fprintf(out, "%s", chs.data);
            }
            std::fprintf(out, "\n");
        } else {
            std::fprintf(out, ">M%llu %zu NP\n", (unsigned long long)layout_no,
                         layout.size());
            auto first_it = dbgi.dbg.find(path[0]);
            if (first_it != dbgi.dbg.end()) {
                for (auto &seq : first_it->second.seq) {
                    chs.require(seq.len);
                    seq.convert_to_ascii((char *)chs.data);
                    std::fprintf(out, "%s", chs.data);
                }
            }
            std::fprintf(out, "\n");
        }
        ++layout_no;
    }
}
namespace layout_graph_lago{
    struct lago{
        gtl::vector<u64> forward_lane;
        gtl::vector<u64> reverse_lane;


    };
}

// TODO if needed
struct lsg_cleaning {
    void cut_tips(layout_small_graph &lsg) {}
    void pop_bubbles(layout_small_graph &lsg) {}
    void simplify_simple_paths(layout_small_graph &lsg,
                               gtl::vector<u64> &to_fix) {}
    void fix_edges(layout_small_graph &lsg, gtl::vector<u64> &to_fix) {}
    void run(layout_small_graph &lsg, int iterations, int threads) {
        gtl::vector<u64> to_fix;
        simplify_simple_paths(lsg, to_fix);
        fix_edges(lsg, to_fix);
        for (int i = 0; i < iterations; ++i) {
            cut_tips(lsg);
            pop_bubbles(lsg);
            simplify_simple_paths(lsg, to_fix);
            fix_edges(lsg, to_fix);
        }
    }
};
#endif // #if 0

#endif // _LAYOUT_GRAPH_H
