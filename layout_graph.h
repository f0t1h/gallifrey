#ifndef _LAYOUT_GRAPH_H
#define _LAYOUT_GRAPH_H

#include <algorithm>
#include <cstring>
#include <functional>
#include <gtl/phmap_fwd_decl.hpp>
#include <gtl/vector.hpp>
#include <queue>

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
            dups.push_back({id, dup_id, orig.core_count, std::move(seq_copy)});
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

#endif
