#ifndef _LAYOUT_GRAPH_H
#define _LAYOUT_GRAPH_H

#include <gtl/vector.hpp>
#include <gtl/phmap_fwd_decl.hpp>
#include "graph.h"
#include <queue>

using layout_small_graph = gtl::flat_hash_map<uint64_t, dbg_index::bidbg_val>;

// Find the longest path in the DAG by node count.
// Returns pairs of (node_id, sequences) along the path (source → sink order).
inline gtl::vector<std::pair<u64, gtl::vector<SequenceContainerType>>>
longest_path(layout_small_graph &lsg) {
    if (lsg.empty()) return {};

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
        fprintf(stderr, "Warning: layout graph contains cycles (%zu of %zu nodes unreachable)\n",
                lsg.size() - processed, lsg.size());
    }

    // Reconstruct path and collect node_id + sequences (source → sink).
    gtl::vector<u64> path;
    for (u64 cur = best_node;; ) {
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


//TODO if needed
struct lsg_cleaning{
    void cut_tips(layout_small_graph &lsg){

    }
    void pop_bubbles(layout_small_graph &lsg){

    }
    void simplify_simple_paths(layout_small_graph &lsg, gtl::vector<u64> &to_fix){

    }
    void fix_edges(layout_small_graph &lsg, gtl::vector<u64> &to_fix){

    }
    void run(layout_small_graph &lsg, int iterations, int threads){

        gtl::vector<u64> to_fix;
        simplify_simple_paths(lsg, to_fix);
        fix_edges(lsg, to_fix);
        for(int i = 0; i < iterations; ++i){
            cut_tips(lsg);
            pop_bubbles(lsg);
            simplify_simple_paths(lsg, to_fix);
            fix_edges(lsg, to_fix);
        }
    }
};

#endif
