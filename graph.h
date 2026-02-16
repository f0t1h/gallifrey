#ifndef _GRAPH_H_
#define _GRAPH_H_

#define PREV_WITH_VEC
#define HT_SHARD_BC 8

#include <atomic>
#include <blockingconcurrentqueue.h>
#include <cstdint>
#include <cstdio>
#include <gtl/phmap.hpp>
#include <gtl/vector.hpp>
#include <mutex>
#include <string>
#include <thread>
#include <queue>
#include <vector>

#include "asm.h"
#include "bloomfilter.h"
#include "config.h"
#include "lcmer.h"

template <class map, class F>
static void for_each_map(const map &mp, const F &fCallback, int threads) {
    std::vector<std::thread> workers;
    for (int tid = 0; tid < threads; ++tid) {
        workers.emplace_back(
            [&](const map &mp, int tid) {
                for (size_t i = tid; i < mp.subcnt(); i += threads) {
                    mp.with_submap(i, [&](const map::EmbeddedSet &set) {
                        for (const auto &p : set) {
                            fCallback(tid, p);
                        }
                    });
                }
            },
            std::cref(mp), tid);
    }
    for (auto &worker : workers) {
        worker.join();
    }
}

template <std::integral int_type>
consteval int_type nextpowerof2(int_type v) {
    v--;
    for (size_t i = 1; i < sizeof(int_type) * 8; i <<= 1) {
        v |= v >> i;
    }
    v++;
    return v;
}

using u64 = uint64_t;
using u32 = uint32_t;
using DefaultBF =
    bloom_interface<bloom_ext::BlockedLockedBloomFilter<0, std::mutex>>;
template <class BFT = DefaultBF>
struct dbg_index_base {
    using BloomFilterType = BFT;
    static constexpr u64 lcindex_increment =
        nextpowerof2((u64)MAXIMUM_WORK_THREADS);
    using mutex_type = std::mutex;
    struct bidbg_val {
        bidbg_val(unsigned core_count)
            : core_count(core_count), forward_edge_count(0), seq{}, prev{} {}
        bidbg_val(unsigned core_count, u64 prev)
            : core_count(core_count),
              forward_edge_count(0),
              seq{},
              prev{prev} {}

        unsigned core_count : 31l;
        unsigned visited : 1;
        unsigned forward_edge_count;

        gtl::vector<SequenceContainerType> seq;

        constexpr bool redirected() const { return seq.empty(); }
#ifdef PREV_WITH_SET
        gtl::flat_hash_set<u64> prev;
        void emplace_prev(u64 val) { prev.emplace(val); }
#else
#ifdef PREV_WITH_VEC

        gtl::vector<u64> prev;
        void emplace_prev(u64 val) {
            auto it = std::find(prev.begin(), prev.end(), val);
            if (it == prev.end()) {
                prev.push_back(val);
            }
        }
#endif
#endif
        constexpr double avg_core_count() const {
            return (double)core_count / seq.size();
        }
    };
    using dbg_map_t = gtl::parallel_flat_hash_map<
        u64, bidbg_val, std::hash<u64>, gtl::priv::hash_default_eq<u64>,
        std::allocator<std::pair<const u64, bidbg_val>>, HT_SHARD_BC,
        mutex_type>;

    struct lci_type {
        u64 id;
        // u64 moved : 1;
        static constexpr size_t max_id = std::numeric_limits<u64>::max();
        static consteval lci_type empty() { return {0}; }
        //    u64 pos; // TODO remove this
        // lci_type() : id(0), pos(-1) {}
        // lci_type(u64 _id, u64 _pos) : id(_id), pos(_pos) {}
        // lci_type(u64 _id) : id(_id), pos(-1) {}
    };
    using lcmer_index_t = gtl::parallel_flat_hash_map<
        lcmer<LCMER_L>, lci_type, std::hash<lcmer<LCMER_L>>,
        gtl::priv::hash_default_eq<lcmer<LCMER_L>>,
        std::allocator<std::pair<const lcmer<LCMER_L>, lci_type>>, HT_SHARD_BC,
        mutex_type>;

    dbg_map_t dbg;
    lcmer_index_t lci;
    gtl::flat_hash_set<u64> active_set;

    // merge_info_map_t merge_map;
    BloomFilterType bf;

    dbg_index_base(double false_positive_rate, uint64_t projected_count,
                   uint64_t seed)
        : bf{false_positive_rate, projected_count, seed} {}

    // dbg_index(const dbg_index &other)
    lci_type get(const lcmer<LCMER_L> &l) const {
        auto it = lci.find(l);
        return it == lci.end() ? lci_type::empty() : it->second;
        // if (lci.contains(l)) {
        //     return lci.at(l);
        // }
        // // lci.if_contains(l,[&ret](const lcmer_index_t::value_type &v){ ret
        // = v.second;}); return 0;
    }
    void count_forward_edges() {
        for (auto &node : dbg) {
            for (auto &edge : node.second.prev) {
                auto it = dbg.find(edge);
                if (it == dbg.end()) {
                    continue;
                }
                it->second.forward_edge_count++;
            }
        }
    }
    // void simplify_graph(int max_iterations = 1, size_t max_tip_length = 5,
    //                     size_t bubble_diff_eps = 2, size_t min_core_depth =
    //                     3, double bubble_depth_epsilon = 2) {
    //     gtl::flat_hash_map<u64, u64> redirect_deleted;
    //     auto try_simplify_simple_path = [&](u64 id, auto current) -> bool {
    //         // Find a node we cannot simplify
    //         // Iterate backwards and merge until we cannot simplify
    //         // o->o->o->o->o->o -> o->o

    //         if (current == dbg.end()) {
    //             return false;
    //         }
    //         unsigned first_fec = current->second.forward_edge_count;
    //         bool simplified = false;
    //         auto prev_it = dbg.find(current->second.prev[0]);
    //         gtl::vector<bidbg_val *> merged;
    //         unsigned cc_sum = 0;
    //         while (prev_it != dbg.end() && current->second.prev.size() == 1
    //         &&
    //                prev_it->second.forward_edge_count == 1) {
    //             if (prev_it != dbg.end()) {
    //                 // for (auto it = current->second.seq.begin();
    //                 //      it != current->second.seq.end(); ++it) {
    //                 merged.emplace_back(&current->second);
    //                 // }
    //                 cc_sum += current->second.core_count;
    //                 // dbg.erase(current->first);
    //                 simplified = true;
    //             } else {  // Prev is erased
    //             }
    //             current = prev_it;
    //             if (prev_it->second.prev.size() > 0) {
    //                 prev_it = dbg.find(current->second.prev[0]);
    //             } else {
    //                 break;
    //             }
    //         }
    //         if (current->first != id) {
    //             redirect_deleted[id] = current->first;
    //         }

    //         current->second.forward_edge_count = first_fec;
    //         current->second.core_count += cc_sum;
    //         for (auto it = merged.rbegin(); it != merged.rend(); ++it) {
    //             for (auto &&node : (*it)->seq) {
    //                 current->second.seq.emplace_back(std::move(node));
    //             }
    //             (*it)->seq.clear();
    //         }

    //         return simplified;
    //     };
    //     for (int iter = 0; iter < max_iterations; ++iter) {
    //         // Remove simple paths

    //         for (auto &dn : lci) {
    //             auto current = dbg.find(dn.second.id);
    //             if (current == dbg.end()) {
    //                 dn.second.id = lci_type::max_id;
    //                 // dn.second.id = redirect_deleted[dn.second.id];
    //                 // dn.second.moved = 1;
    //                 continue;
    //             }
    //             if (current->second.redirected()) {
    //                 dn.second.moved = true;
    //                 // dn.second.id = lci_type::max_id;
    //             }
    //             // Simple paths
    //             if ((current->second.forward_edge_count != 1 &&
    //                  current->second.prev.size() == 1)) {
    //                 try_simplify_simple_path(dn.second.id, current);
    //             } else if (current->second.prev.size() > 1) {
    //                 for (auto &prev_id : current->second.prev) {
    //                     auto prev_it = dbg.find(prev_id);
    //                     if (prev_it != dbg.end() &&
    //                         prev_it->second.prev.size() == 1) {
    //                         try_simplify_simple_path(prev_id, prev_it);
    //                     }
    //                 }
    //             }

    //             // auto prev_it = dbg.find(node.prev[0]);
    //             // if (prev_it != dbg.end()) {
    //             //     auto &prev_node = prev_it->second;
    //             //     dn.second.pos = prev_node.merged.size();
    //             //     prev_node.merged.push_back(dn.second.id);
    //             //     prev_node.cc = node.cc;
    //             //     dbg.erase(dn.second.id);
    //             // } else {  // Prev is erased
    //             // }
    //             // merge_map.emplace(dn.second, node.prev[0], prev_node)
    //         }
    //         for (auto &node : dbg) {
    //             size_t idx = 0;
    //             int self_idx = -1;
    //             for (auto &p : node.second.prev) {
    //                 auto it = redirect_deleted.find(p);
    //                 if (it != redirect_deleted.end()) {
    //                     p = it->second;
    //                 } else if (p == node.first) {  // Self edge
    //                     node.second.forward_edge_count--;
    //                     self_idx = idx;
    //                 }
    //                 ++idx;
    //             }
    //             if (self_idx != -1) {
    //                 node.second.prev.erase(node.second.prev.begin() +
    //                 self_idx);
    //             }
    //         }

    //         // Tips
    //         for (auto &dn : lci) {
    //             auto current = dbg.find(dn.second.id);
    //             if (current == dbg.end()) {
    //                 dn.second.id = lci_type::max_id;
    //                 continue;
    //             }
    //             if (current->second.avg_core_count() < min_core_depth) {
    //                 // Forlorn
    //                 if (current->second.forward_edge_count == 0 &&
    //                     current->second.prev.size() == 0) {
    //                     dbg.erase(current);
    //                 } else if (current->second.forward_edge_count == 0 &&
    //                            current->second.prev.size() == 1) {
    //                     auto prev_it = dbg.find(current->second.prev[0]);
    //                     if (prev_it != dbg.end()) {
    //                         prev_it->second.forward_edge_count--;
    //                     }
    //                     dbg.erase(current);
    //                 }
    //                 // else if (current->second.fec == 1 &&
    //                 //            current->second.prev.size() == 0) {
    //                 //     dbg.erase(current);
    //                 // }
    //             } else {
    //                 for (u64 &pid : current->second.prev) {
    //                     auto prev_it = dbg.find(pid);
    //                     if (prev_it != dbg.end()) {
    //                         if (prev_it->second.avg_core_count() <
    //                                 min_core_depth &&
    //                             prev_it->second.prev.size() == 0 &&
    //                             prev_it->second.forward_edge_count == 1) {
    //                             dbg.erase(prev_it);
    //                             pid = std::numeric_limits<u64>::max();
    //                         }
    //                     }
    //                 }
    //                 auto &pvec = current->second.prev;
    //                 pvec.erase(std::remove_if(
    //                                pvec.begin(), pvec.end(),
    //                                [](u64 val) {
    //                                    return val ==
    //                                           std::numeric_limits<u64>::max();
    //                                }),
    //                            pvec.end());
    //                 for (u64 &pid : current->second.prev) {
    //                     auto prev_it = dbg.find(pid);
    //                     if (prev_it != dbg.end()) {
    //                         if (prev_it->second.avg_core_count() <
    //                                 current->second.avg_core_count() &&
    //                             prev_it->second.seq.size() < max_tip_length
    //                             && current->second.seq.size() >
    //                                 2 * max_tip_length &&
    //                             prev_it->second.prev.size() == 0 &&
    //                             prev_it->second.forward_edge_count == 1) {
    //                             dbg.erase(prev_it);
    //                             pid = std::numeric_limits<u64>::max();
    //                         }
    //                     }
    //                 }
    //                 pvec.erase(std::remove_if(
    //                                pvec.begin(), pvec.end(),
    //                                [](u64 val) {
    //                                    return val ==
    //                                           std::numeric_limits<u64>::max();
    //                                }),
    //                            pvec.end());
    //             }
    //             // We have non-zero forward edges, and more than one prev
    //             // edges
    //             // if (current->second.fec > 0 &&
    //             //     current->second.prev.size() > 0) {
    //             //     // Check if one of the prev edges is not terminating
    //             //     bool non_terminating_exists = false;
    //             //     size_t min_edge_len = INT_MAX;
    //             //     size_t max_edge_len = INT_MIN;
    //             //     for (auto &prev_id : current->second.prev) {
    //             //         auto prev_it = dbg.find(prev_id);
    //             //         if (prev_it == dbg.end()) {
    //             //             continue;
    //             //         }
    //             //         if (prev_it->second.prev.size() > 0) {
    //             //             non_terminating_exists = true;
    //             //         }
    //             //         min_edge_len =
    //             //             std::min(min_edge_len,
    //             //             prev_it->second.seq.size());
    //             //         max_edge_len =
    //             //             std::max(max_edge_len,
    //             //             prev_it->second.seq.size());
    //             //     }
    //             //     // Remove every tip if non_terminating_exists

    //             //     if (non_terminating_exists) {
    //             //         for (auto &prev_id : current->second.prev) {
    //             //             auto prev_it = dbg.find(prev_id);
    //             //             if (prev_it == dbg.end()) {
    //             //                 // delete from prev
    //             //                 current->second.prev.erase(
    //             // std::remove(current->second.prev.begin(),
    //             // current->second.prev.end(),
    //             //                                 prev_id),
    //             //                     current->second.prev.end());
    //             //                 continue;
    //             //             }
    //             //             if (prev_it->second.prev.size() == 0 &&
    //             //                 prev_it->second.fec == 1 &&
    //             //                 prev_it->second.seq.size() <
    //             max_tip_length)
    //             //                 {
    //             //                 // Erase the node from the
    //             //                 current->second.prev
    //             //                 current->second.prev.erase(
    //             // std::remove(current->second.prev.begin(),
    //             // current->second.prev.end(),
    //             //                                 prev_id),
    //             //                     current->second.prev.end());
    //             //                 dbg.erase(prev_it);
    //             //             }
    //             //         }
    //             //     } else if (max_edge_len - min_edge_len >
    //             max_tip_length)
    //             //     {
    //             //         for (auto &prev_id : current->second.prev) {
    //             //             auto prev_it = dbg.find(prev_id);
    //             //             if (prev_it == dbg.end()) {
    //             //                 continue;
    //             //             }
    //             //             if (prev_it->second.prev.size() == 0 &&
    //             //                 prev_it->second.fec == 1 &&
    //             //                 prev_it->second.seq.size() <
    //             //                     max_tip_length - max_tip_length) {
    //             //                 // Erase the node from the
    //             //                 current->second.prev

    //             //                 current->second.prev.erase(
    //             // std::remove(current->second.prev.begin(),
    //             // current->second.prev.end(),
    //             //                                 prev_id),
    //             //                     current->second.prev.end());
    //             //                 dbg.erase(prev_it);
    //             //             }
    //             //         }
    //             //     }
    //             // }
    //         }

    //         // Buubles
    //         for (auto &dn : lci) {
    //             auto current = dbg.find(dn.second.id);
    //             if (current == dbg.end()) {
    //                 dn.second.id = -1;
    //                 continue;
    //             }

    //             if (current->second.prev.size() >=
    //                 2) {  // Simple case with 2 should expand later to many
    //                 case auto p1 = dbg.find(current->second.prev[0]); auto p2
    //                 = dbg.find(current->second.prev[1]); if (p1 != dbg.end()
    //                 && p2 != dbg.end()) {
    //                     if (p1->second.prev.size() == 1 &&
    //                         p2->second.prev.size() == 1) {
    //                         if (p1->second.prev[0] == p2->second.prev[0]) {
    //                             // BUBBLE FOUND
    //                             if (p2->second.avg_core_count() /
    //                                     p1->second.avg_core_count() <
    //                                 bubble_depth_epsilon) {
    //                                 // Remove p2
    //                                 dbg.erase(p2);
    //                                 current->second.prev.erase(std::next(
    //                                     current->second.prev.begin()));
    //                                 auto nn = dbg.find(p1->second.prev[0]);
    //                                 if (nn != dbg.end()) {
    //                                     nn->second.forward_edge_count--;
    //                                 }
    //                             } else if (p1->second.avg_core_count() /
    //                                            p2->second.avg_core_count() <
    //                                        bubble_depth_epsilon) {
    //                                 // Remove p1
    //                                 dbg.erase(p1);
    //                                 current->second.prev[0] =
    //                                     current->second.prev[1];
    //                                 current->second.prev.erase(std::next(
    //                                     current->second.prev.begin()));

    //                                 auto nn = dbg.find(p2->second.prev[0]);
    //                                 if (nn != dbg.end()) {
    //                                     nn->second.forward_edge_count--;
    //                                 }
    //                             }
    //                         }
    //                     }
    //                 }
    //             }
    //         }
    //         // Death by thousand cuts
    //         for (auto &dn : lci) {
    //             auto current = dbg.find(dn.second.id);
    //             if (current == dbg.end()) {
    //                 continue;
    //             }
    //             if (current->second.avg_core_count() < 2 &&
    //                 current->second.prev.size() == 1 &&
    //                 current->second.forward_edge_count == 1) {
    //                 auto prev_it = dbg.find(current->second.prev[0]);
    //                 if (prev_it != dbg.end())
    //                     prev_it->second.forward_edge_count--;
    //                 dbg.erase(current);
    //             }
    //         }
    //         // Simple Loops

    //         for (auto &dn : lci) {
    //             auto current = dbg.find(dn.second.id);
    //             if (current == dbg.end()) {
    //                 continue;
    //             }
    //             if (current->second.forward_edge_count == 1 &&
    //                 current->second.prev.size() == 1) {
    //                 auto prev_it = dbg.find(current->second.prev[0]);
    //                 if (prev_it != dbg.end()) {
    //                     auto &prepre = prev_it->second.prev;
    //                     for (auto pp = prepre.begin(); pp != prepre.end();
    //                          ++pp) {
    //                         if (*pp == dn.second.id) {  // Loop here
    //                             *pp = std::numeric_limits<u64>::max();
    //                             if (current->second.core_count >=
    //                                 2 * prev_it->second
    //                                         .core_count) {  // TODO: this is
    //                                         not
    //                                                         // right.
    //                                 // Need to duplicate the
    //                                 // sequence
    //                                 auto &tomerge = current->second.seq;
    //                                 for (auto it = tomerge.begin();
    //                                      it != tomerge.end(); ++it) {
    //                                     prev_it->second.seq.emplace_back(
    //                                         std::move(*it));
    //                                 }
    //                             }

    //                             prev_it->second.forward_edge_count--;
    //                             dbg.erase(current);
    //                         }
    //                     }
    //                 }
    //             }
    //         }
    //         // Fix edges
    //         for (auto &node : dbg) {
    //             auto &pvec = node.second.prev;
    //             for (auto it = pvec.begin(); it != pvec.end(); ++it) {
    //                 if (dbg.find(*it) == dbg.end()) {
    //                     *it = std::numeric_limits<u64>::max();
    //                 }
    //             }
    //             pvec.erase(
    //                 std::remove_if(pvec.begin(), pvec.end(),
    //                                [](u64 val) {
    //                                    return val ==
    //                                           std::numeric_limits<u64>::max();
    //                                }),
    //                 pvec.end());
    //         }
    //         redirect_deleted.clear();
    //     }

    //     erase_if(lci, [](const auto &node) { // MARK
    //         return node.second.id == lci_type::max_id;
    //     });
    //     // // Remove bubbles and tips
    //     // for (auto &node : dbg) {
    //     //     if (node.second.prev.size() == 2) {  // can be simplified.
    //     //         u64 p1 = node.second.prev[0];
    //     //         u64 p2 = node.second.prev[1];
    //     //         if (dbg[p1].cc == 1 && dbg[p2].cc == 1 &&
    //     //             dbg[p1].prev.size() == 1 && dbg[p2].prev.size() == 1
    //     &&
    //     //             dbg[p1].prev[0] == dbg[p2].prev[0]) {
    //     //             ///     o
    //     //             ///   o< >o -> o->o->o
    //     //             ///     o
    //     //         } else if (node.second.prev.size() == 0) {
    //     //             // Remove tips
    //     //             /// ..o->o->o..       o->o->o
    //     //             ///       \     ->
    //     //             ///        o
    //     //         }
    //     //     }
    //     // }
    // }
    u64 emplace_nobloomfilter(const lcmer<LCMER_L> &l, u64 &idx, int tid) {
        auto vv = lci.try_emplace(l, idx | tid);
        if (vv.second) {  // Means new item on the table
            idx += lcindex_increment;
            // threadcount.
        }
        return vv.first->second.id;
    }
    u64 emplace(const lcmer<LCMER_L> &l, u64 &idx, int tid) {
        auto new_hash = lci.hash(l);
        if (!bf.contains(new_hash)) {
            bf.insert(new_hash);
            return 0;
        }
        auto vv = lci.try_emplace_with_hash(new_hash, l, idx | tid);
        if (vv.second) {  // Means new item on the table
            idx += lcindex_increment;
            // threadcount.
        }
        return vv.first->second.id;
    }
    gtl::vector<u64> shortest_path_bfs(u64 source,
                                           u64 target) const {
        if (source == target) return gtl::vector<u64>{source};
        if (source == 0 || target == 0) return {};


        std::queue<u64> q;
        gtl::flat_hash_map<u64, u64> parent;
        gtl::flat_hash_map<u64, uint8_t> visited;

        visited[source] = 1;
        q.push(source);

        while (!q.empty()) {
            const u64 cur = q.front();
            q.pop();


            auto it = dbg.find(cur);
            if(it == dbg.end()){
                continue;
            }
            for (u64 nb : it->second.prev) {
                if (visited.find(nb) != visited.end()) continue;
                auto neit = dbg.find(nb);
                if (neit == dbg.end()) continue;

                visited.emplace(nb, 1);
                parent.emplace(nb, cur);

                if (nb == target) {
                    gtl::vector<u64> path;
                    u64 v = target;
                    path.push_back(v);
                    while (v != source) {
                        auto pit = parent.find(v);
                        if (pit == parent.end()) return {};
                        v = pit->second;
                        path.push_back(v);
                    }
                    std::reverse(path.begin(), path.end());
                    return path;
                }

                q.push(nb);
            }
        }
        return {};
    }

};

using dbg_index = dbg_index_base<DefaultBF>;
using dbg_index_light = dbg_index_base<MockBF>;

struct FlatDBG {
    struct Node {
        gtl::vector<u64> prev;
        gtl::vector<u64> next;
        gtl::vector<SequenceContainerType> seq;
    };
    gtl::vector<Node> nodes;
};
struct GraphCleaning {
    using bidbg_val = dbg_index::bidbg_val;
    using lci_type = dbg_index::lci_type;
    using dbg_map_t = dbg_index::dbg_map_t;
    using lcmer_index_t = dbg_index::lcmer_index_t;
    using redir_map = gtl::parallel_flat_hash_map_d<u64, u64, 4, std::mutex>;
    static auto attempt_simplify_one_path(dbg_map_t &dbg,
                                          redir_map &redirect_deleted, u64 id,
                                          auto starting_from, auto &to_delete)
        -> bool {
        // Find a node we cannot simplify
        // Iterate backwards and merge until we cannot simplify
        // o->o->o->o->o->o -> o->o

        auto current = starting_from;
        if (current == dbg.end()) {
            return false;
        }
        unsigned first_fec = current->second.forward_edge_count;
        bool simplified = false;
        auto prev_it = dbg.find(current->second.prev[0]);
        gtl::vector<dbg_map_t::iterator> merged;
        unsigned cc_sum = 0;
        while (prev_it != dbg.end() && current->second.prev.size() == 1 &&
               prev_it->second.forward_edge_count == 1) {
            if (prev_it != dbg.end()) {
                cc_sum += current->second.core_count;
                merged.push_back(current);

                simplified = true;
            } else {  // Prev is erased
            }
            current = prev_it;
            if (prev_it->second.prev.size() > 0) {
                prev_it = dbg.find(current->second.prev[0]);
            } else {
                break;
            }
        }
        if (current->first != id) {
            redirect_deleted[id] = current->first;
        }

        current->second.forward_edge_count = first_fec;
        current->second.core_count += cc_sum;
        for (auto it = merged.rbegin(); it != merged.rend(); ++it) {
            // for (auto &seq : (*it)->second.seq)
            gtl::vector<SequenceContainerType> &seqs = (*it)->second.seq;
            for (auto sit = seqs.begin(); sit != seqs.end(); ++sit) {
                to_delete.emplace(sit->origin_node, current->first);

                current->second.seq.emplace_back(std::move(*sit));
            }
            seqs.clear();

            to_delete.emplace((*it)->first, current->first);
        }

        return simplified;
    };

    static auto attempt_simplify_one_path_lockless(dbg_map_t &dbg,
                                                   redir_map &redirect_deleted,
                                                   u64 id, auto starting_from,
                                                   auto &to_delete) -> bool {
        // Find a node we cannot simplify
        // Iterate backwards and merge until we cannot simplify
        // o->o->o->o->o->o -> o->o
        gtl::LockableBaseImpl<gtl::NullMutex>::DoNothing no_lock;

        auto current = starting_from;
        if (current == dbg.end()) {
            return false;
        }
        unsigned first_fec = current->second.forward_edge_count;
        bool simplified = false;
        auto prev_hash = dbg.hash(current->second.prev[0]);
        auto prev_it = dbg.find(current->second.prev[0], prev_hash, no_lock);
        gtl::vector<dbg_map_t::iterator> merged;
        unsigned cc_sum = 0;
        while (prev_it != dbg.end() && current->second.prev.size() == 1 &&
               prev_it->second.forward_edge_count == 1) {
            if (prev_it != dbg.end()) {
                cc_sum += current->second.core_count;
                merged.push_back(current);

                simplified = true;
            } else {  // Prev is erased
            }
            current = prev_it;
            if (prev_it->second.prev.size() > 0) {
                prev_hash = dbg.hash(current->second.prev[0]);
                prev_it = dbg.find(current->second.prev[0], prev_hash, no_lock);
            } else {
                break;
            }
        }
        if (current->first != id) {
            redirect_deleted[id] = current->first;
        }

        current->second.forward_edge_count = first_fec;
        current->second.core_count += cc_sum;
        for (auto it = merged.rbegin(); it != merged.rend(); ++it) {
            // for (auto &seq : (*it)->second.seq)
            gtl::vector<SequenceContainerType> &seqs = (*it)->second.seq;
            for (auto sit = seqs.begin(); sit != seqs.end(); ++sit) {
                to_delete.emplace(sit->origin_node, current->first);

                current->second.seq.emplace_back(std::move(*sit));
            }
            seqs.clear();

            to_delete.emplace((*it)->first, current->first);
        }

        return simplified;
    };
    // static void dfs_simplify(dbg_index::dbg_map_t &dbg,
    //                          dbg_index::lcmer_index_t &lci,
    //                          gtl::flat_hash_map<u64, u64> &redirect_deleted)
    //                          {
    //     auto dfs_visit = [](dbg_index::dbg_map_t &dm,
    //                         dbg_index::dbg_map_t::value_type &val) -> void {
    //         val.second.visited = true;
    //     };
    //     auto pop = [](gtl::vector<u64> &vec) {
    //         u64 val = vec.back();
    //         vec.pop_back();
    //         return val;
    //     };
    //     auto modify_foo = [](dbg_map_t::value_type &val) {
    //         val.second.visited = true;
    //     };
    //     gtl::vector<u64> to_visit;
    //     to_visit.reserve(1024);
    //     for (auto &dbn : dbg) {
    //         while (!to_visit.empty()) {
    //             auto val = pop(to_visit);
    //             auto it = dbg.find(val);
    //             if (it != dbg.end()) {
    //                 if (it->second.visited) {
    //                     continue;
    //                 }
    //                 dfs_visit(dbg, *it);
    //                 for (u64 neighbour : it->second.prev) {
    //                     to_visit.push_back(neighbour);
    //                 }
    //             }
    //         }
    //         if (!dbn.second.visited) {
    //             to_visit.push_back(dbn.first);
    //         }
    //     }
    // }

    static auto simplify_simple_parallel_mark(dbg_index::dbg_map_t &dbg,
                                              int threads) {
        gtl::vector<gtl::vector<std::pair<u64, dbg_index::dbg_map_t::iterator>>>
            to_simplify_from(threads);
        gtl::LockableBaseImpl<gtl::NullMutex>::DoNothing no_lock;

        for_each_map(
            dbg,
            [&to_simplify_from, &no_lock, &dbg](
                int tid, const dbg_index::dbg_map_t::value_type &item) {
                // printf("TTTT:%d\n", tid);
                if (item.second.forward_edge_count != 1 &&
                    item.second.prev.size() == 1) {
                    auto ih = dbg.hash(item.first);
                    auto iit = dbg.find(item.first, ih, no_lock);
                    to_simplify_from[tid].emplace_back(item.first, iit);
                } else if (item.second.prev.size() > 1) {
                    for (auto &prev_id : item.second.prev) {
                        auto prev_hash = dbg.hash(prev_id);

                        auto prev_it = dbg.find(prev_id, prev_hash, no_lock);
                        if (prev_it != dbg.end() &&
                            prev_it->second.prev.size() == 1) {
                            to_simplify_from[tid].emplace_back(prev_id,
                                                               prev_it);
                        }
                    }
                }
            },
            threads);
        // for (size_t i = 0; i < dbg.subcnt(); ++i) {
        //     dbg.with_submap(
        //         i, [&to_simplify_from, &no_lock,
        //             &dbg](const dbg_index::dbg_map_t::EmbeddedSet &set) {
        //             for (auto it = set.begin(); it != set.end(); ++it) {
        //                 // for (const auto &p : set) {
        //                 if (it->second.forward_edge_count != 1 &&
        //                     it->second.prev.size() == 1) {
        //                     auto ih = dbg.hash(it->first);
        //                     auto iit = dbg.find(it->first, ih, no_lock);
        //                     to_simplify_from.emplace_back(it->first, iit);
        //                 } else if (it->second.prev.size() > 1) {
        //                     for (auto &prev_id : it->second.prev) {
        //                         auto prev_hash = dbg.hash(prev_id);

        //                         auto prev_it =
        //                             dbg.find(prev_id, prev_hash, no_lock);
        //                         if (prev_it != dbg.end() &&
        //                             prev_it->second.prev.size() == 1) {
        //                             to_simplify_from.emplace_back(prev_id,
        //                                                           prev_it);
        //                         }
        //                     }
        //                 }
        //             }
        //         });
        // }
        return to_simplify_from;
    }
    static void simplify_simple_parallel_sweep(
        dbg_index::dbg_map_t &dbg, dbg_index::lcmer_index_t &lci,
        const gtl::vector<
            gtl::vector<std::pair<u64, dbg_index::dbg_map_t::iterator>>>
            &to_simplify,
        redir_map &redirect_deleted, auto &active_set) {
        gtl::LockableBaseImpl<gtl::NullMutex>::DoNothing no_lock;

        // int tid=0;
        gtl::vector<std::thread> threads;
        std::mutex del_mutex;
        for (size_t i = 0; i < to_simplify.size(); ++i) {
            threads.emplace_back(
                [&](int i) {
                    gtl::flat_hash_map<u64, u64> to_delete;
                    // printf("tid: %d\n", tid++);
                    for (const auto &[node, it] : to_simplify[i]) {
                        // printf("%llu\n", node);
                        // auto hv = dbg.hash(node);
                        // auto it = dbg.find(node, hv, no_lock);
                        // attempt_simplify_one_path(dbg,lci, redirect_deleted,
                        // it.first, auto starting_from, auto &to_delete);
                        attempt_simplify_one_path_lockless(
                            dbg, redirect_deleted, node, it, to_delete);
                    }
                    for (auto &node : lci) {
                        auto it = to_delete.find(node.second.id);
                        if (it != to_delete.end()) {
                            if (it->second != 0) {
                                node.second.id =
                                    it->second;  // Redirect to new node
                            }
                        }
                    }
                    std::lock_guard lock(del_mutex);
                    for (const auto &[del, to] : to_delete) {
                        active_set.erase(del);
                    }
                },
                i);
        }

        for (auto &thread : threads) {
            thread.join();
        }

        for (auto &node : dbg) {
            int self_idx = -1;
            for (auto &p : node.second.prev) {
                auto it = redirect_deleted.find(p);
                if (it != redirect_deleted.end()) {
                    p = it->second;
                } else if (p == node.first) {  // Self edge // TODO: reconsider
                                               // how this is done
                    // node.second.forward_edge_count--;
                    // self_idx = idx;
                }
            }
            if (self_idx != -1) {
                node.second.prev.erase(node.second.prev.begin() + self_idx);
            }
        }
    }
    static void simplify_simple_paths(dbg_index::dbg_map_t &dbg,
                                      dbg_index::lcmer_index_t &lci,
                                      redir_map &redirect_deleted,
                                      auto &active_set) {
        // TODO Parallel idea
        // First run a parallel discovery step which does not modify the graph
        // and finds the pairs of nodes we want to simplify (from, to) and
        // extracts the sequences. On a single thread we fix the edges and add
        // sequences to the target. Then in parallel we redirect marked nodes to
        // simplified () Should reduce the time needed 4X (20s to 5s on HG) with
        // 8 threads

        // for (auto &dn : lci) {
        gtl::flat_hash_map<u64, u64> to_delete;
        for (u64 node : active_set) {
            auto current = dbg.find(node);
            if (current == dbg.end()) {
                to_delete.emplace(node, 0);
                continue;
            }

            // Simple paths
            if ((current->second.forward_edge_count != 1 &&
                 current->second.prev.size() == 1)) {
                attempt_simplify_one_path_lockless(dbg, redirect_deleted, node,
                                                   current, to_delete);
            } else if (current->second.prev.size() > 1) {
                for (auto &prev_id : current->second.prev) {
                    auto prev_it = dbg.find(prev_id);
                    if (prev_it != dbg.end() &&
                        prev_it->second.prev.size() == 1) {
                        attempt_simplify_one_path_lockless(
                            dbg, redirect_deleted, prev_id, prev_it, to_delete);
                    }
                }
            }
        }

        for (auto &node : lci) {
            auto it = to_delete.find(node.second.id);
            if (it != to_delete.end()) {
                if (it->second != 0) {
                    node.second.id = it->second;  // Redirect to new node
                }
            }
        }
        for (auto [del, to] : to_delete) {
            active_set.erase(del);
        }

        for (auto &node : dbg) {
            int self_idx = -1;
            for (auto &p : node.second.prev) {
                auto it = redirect_deleted.find(p);
                if (it != redirect_deleted.end()) {
                    p = it->second;
                } else if (p == node.first) {  // Self edge // TODO: reconsider
                                               // how this is done
                    // node.second.forward_edge_count--;
                    // self_idx = idx;
                }
            }
            if (self_idx != -1) {
                node.second.prev.erase(node.second.prev.begin() + self_idx);
            }
        }
    }
    static void cut_tips(dbg_map_t &dbg, [[maybe_unused]] lcmer_index_t &lci,
                         int min_core_depth, size_t max_tip_length,
                         auto &active_set) {
        // for (auto &dn : lci) {
        for (u64 node : active_set) {
            auto current = dbg.find(node);
            if (current == dbg.end()) {
                continue;
            }
            if (current->second.redirected()) {
                continue;
            }
            if (current->second.avg_core_count() < min_core_depth) {
                // Forlorn
                if (current->second.forward_edge_count == 0 &&
                    current->second.prev.size() == 0) {
                    dbg.erase(current);  // MARK
                } else if (current->second.forward_edge_count == 0 &&
                           current->second.prev.size() == 1) {
                    auto prev_it = dbg.find(current->second.prev[0]);
                    if (prev_it != dbg.end() || prev_it->second.redirected()) {
                        prev_it->second.forward_edge_count--;
                    }
                    dbg.erase(current);  // MARK
                }
                // else if (current->second.fec == 1 &&
                //            current->second.prev.size() == 0) {
                //     dbg.erase(current);
                // }
            } else {
                for (u64 &pid : current->second.prev) {
                    auto prev_it = dbg.find(pid);
                    if (prev_it != dbg.end() || prev_it->second.redirected()) {
                        if (prev_it->second.avg_core_count() < min_core_depth &&
                            prev_it->second.prev.size() == 0 &&
                            prev_it->second.forward_edge_count == 1) {
                            dbg.erase(prev_it);  // MARK
                            pid = std::numeric_limits<u64>::max();
                        }
                    }
                }
                auto &pvec = current->second.prev;
                pvec.erase(
                    std::remove_if(pvec.begin(), pvec.end(),
                                   [](u64 val) {
                                       return val ==
                                              std::numeric_limits<u64>::max();
                                   }),
                    pvec.end());
                for (u64 &pid : current->second.prev) {
                    auto prev_it = dbg.find(pid);
                    if (prev_it != dbg.end() || prev_it->second.redirected()) {
                        if (prev_it->second.avg_core_count() <
                                current->second.avg_core_count() &&
                            prev_it->second.seq.size() < max_tip_length &&
                            current->second.seq.size() > 2 * max_tip_length &&
                            prev_it->second.prev.size() == 0 &&
                            prev_it->second.forward_edge_count == 1) {
                            dbg.erase(prev_it);  // MARK
                            pid = std::numeric_limits<u64>::max();
                        }
                    }
                }
                pvec.erase(
                    std::remove_if(pvec.begin(), pvec.end(),
                                   [](u64 val) {
                                       return val ==
                                              std::numeric_limits<u64>::max();
                                   }),
                    pvec.end());
            }
        }
    }

    // Bubble depth delta calculation
    ///
    static void pop_bubbles(dbg_map_t &dbg, [[maybe_unused]] lcmer_index_t &lci,
                            double bubble_depth_epsilon, auto &active_set) {
        for (u64 node : active_set) {
            auto current = dbg.find(node);
            if (current == dbg.end() || current->second.redirected()) {
                continue;
            }

            if (current->second.prev.size() >=
                2) {  // Simple case with 2 should expand later to many case
                auto p1 = dbg.find(current->second.prev[0]);
                auto p2 = dbg.find(current->second.prev[1]);
                if ((p1 != dbg.end() && !p1->second.redirected()) &&
                    (p2 != dbg.end() && !p2->second.redirected())) {
                    if (p1->second.prev.size() == 1 &&
                        p2->second.prev.size() == 1) {
                        if (p1->second.prev[0] == p2->second.prev[0]) {
                            // BUBBLE FOUND
                            if (p2->second.avg_core_count() /
                                    p1->second.avg_core_count() <
                                bubble_depth_epsilon) {
                                // Remove p2
                                dbg.erase(p2);  // MARK
                                current->second.prev.erase(
                                    std::next(current->second.prev.begin()));
                                auto nn = dbg.find(p1->second.prev[0]);
                                if (nn != dbg.end() &&
                                    !nn->second.redirected()) {
                                    nn->second.forward_edge_count--;
                                }
                            } else if (p1->second.avg_core_count() /
                                           p2->second.avg_core_count() <
                                       bubble_depth_epsilon) {
                                // Remove p1
                                dbg.erase(p1);  // MARK
                                current->second.prev[0] =
                                    current->second.prev[1];
                                current->second.prev.erase(
                                    std::next(current->second.prev.begin()));

                                auto nn = dbg.find(p2->second.prev[0]);
                                if (nn != dbg.end() &&
                                    !nn->second.redirected()) {
                                    nn->second.forward_edge_count--;
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    static void cut_low_support_edges(dbg_map_t &dbg,
                                      [[maybe_unused]] lcmer_index_t &lci,
                                      auto &active_set) {
        // Death by thousand cuts
        for (auto &node : active_set) {
            auto current = dbg.find(node);
            if (current == dbg.end()) {
                continue;
            }
            if (current->second.avg_core_count() < 2 &&
                current->second.prev.size() == 1 &&
                current->second.forward_edge_count == 1) {
                auto prev_it = dbg.find(current->second.prev[0]);
                if (prev_it != dbg.end()) prev_it->second.forward_edge_count--;
                dbg.erase(current);
            }
        }
    }
    static void open_simple_loops(dbg_map_t &dbg,
                                  [[maybe_unused]] lcmer_index_t &lci,
                                  auto &active_set) {
        for (auto &node : active_set) {
            auto current = dbg.find(node);
            if (current == dbg.end() || current->second.redirected()) {
                continue;
            }
            if (current->second.forward_edge_count == 1 &&
                current->second.prev.size() == 1) {
                auto prev_it = dbg.find(current->second.prev[0]);
                if (prev_it != dbg.end() && !prev_it->second.redirected()) {
                    auto &prepre = prev_it->second.prev;
                    for (auto pp = prepre.begin(); pp != prepre.end(); ++pp) {
                        if (*pp == node) {  // Loop here
                            *pp = std::numeric_limits<u64>::max();
                            if (current->second.core_count >=
                                2 * prev_it->second
                                        .core_count) {  // TODO: this is
                                                        // not right.
                                // Need to duplicate the
                                // sequence
                                auto &tomerge = current->second.seq;
                                for (auto it = tomerge.begin();
                                     it != tomerge.end(); ++it) {
                                    prev_it->second.seq.emplace_back(
                                        std::move(*it));
                                }
                            }

                            prev_it->second.forward_edge_count--;
                            dbg.erase(current);
                        }
                    }
                }
            }
        }
    }
    static void fix_edges(dbg_map_t &dbg, [[maybe_unused]] lcmer_index_t &lci,
                          [[maybe_unused]] const auto &active_set) {
        // Fix edges
        for (auto &node : dbg) {
            if (node.second.redirected()) continue;
            auto &pvec = node.second.prev;
            for (auto it = pvec.begin(); it != pvec.end(); ++it) {
                if (dbg.find(*it) == dbg.end() ||
                    dbg.find(*it)->second.redirected()) {
                    *it = std::numeric_limits<u64>::max();
                }
            }
            pvec.erase(std::remove_if(
                           pvec.begin(), pvec.end(),
                           [](u64 val) {
                               return val == std::numeric_limits<u64>::max();
                           }),
                       pvec.end());
        }
        // gtl::vector<decltype(lci.begin())> to_purge;
        // for (auto it = lci.begin(); it != lci.end(); ++it) {
        //     auto dbg_it = dbg.find(it->second.id);
        //     if (dbg_it == dbg.end()) {
        //         to_purge.emplace_back(it);
        //     }
        // }

        // for (auto it : to_purge) {
        //     lci.erase(it);
        // }
    }
    // static void parallel_run(dbg_index &dbgi,
    //                          decltype(asm_config::graph_cleaning) cfg,
    //                          int threads) {
    //     gtl::flat_hash_map<u64, u64> redirect_deleted;
    //     simplify_simple_paths(dbgi.dbg, dbgi.lci, redirect_deleted);

    //     // Find Connected Components
    //     // For each connected component simplify
    //     for (int i = 0; i < cfg.iterations; ++i) {
    //         cut_tips(dbgi.dbg, dbgi.lci, 3, cfg.max_tip_length);
    //         pop_bubbles(dbgi.dbg, dbgi.lci,
    //                     cfg.min_bubble_depth_difference_percentage);
    //         // cut_low_support_edges(dbgi.dbg, dbgi.lci);
    //         open_simple_loops(dbgi.dbg, dbgi.lci);
    //         fix_edges(dbgi.dbg, dbgi.lci);
    //         redirect_deleted.clear();
    //         simplify_simple_paths(dbgi.dbg, dbgi.lci, redirect_deleted);
    //     }
    // }
    static void run(dbg_index &dbgi, const asm_config &acfg,
                    bool once = false) {
        auto cfg = acfg.graph_cleaning;
        INIT_TIMERS(GraphCleaning);
        redir_map redirect_deleted;
        decltype(dbgi.active_set) &active_set = dbgi.active_set;
        active_set.reserve(dbgi.dbg.size() + 1);
        for (const auto &node : dbgi.dbg) {
            active_set.insert(node.first);
        }
        // active_set = gtl::flat_hash_set<u64> {dbgi.dbg.begin(),
        // dbgi.dbg.end()}; dfs_simplify(dbgi.dbg, dbgi.lci, redirect_deleted);
        // IF_TIME_CHECKPOINT(!acfg.silent, GraphCleaning, stderr,
        //                    "DFS Simplification");
        if (cfg.iterations > -1) {
            auto to_simplify =
                simplify_simple_parallel_mark(dbgi.dbg, acfg.threads.work);
            simplify_simple_parallel_sweep(dbgi.dbg, dbgi.lci, to_simplify,
                                           redirect_deleted, active_set);
            // simplify_simple_paths(dbgi.dbg, dbgi.lci, redirect_deleted,
            //                       active_set);
            IF_TIME_CHECKPOINT(!acfg.silent, GraphCleaning, stderr,
                               "Initial Simplification");
            for (int i = 1; i < (once ? 2 : cfg.iterations); ++i) {
                fix_edges(dbgi.dbg, dbgi.lci, active_set);
                cut_tips(dbgi.dbg, dbgi.lci, 3, cfg.max_tip_length, active_set);
                // TODO FIX pop_bubbles sequence problem fix
                // pop_bubbles(dbgi.dbg, dbgi.lci,
                //             cfg.min_bubble_depth_difference_percentage,
                //             active_set);
                cut_low_support_edges(dbgi.dbg, dbgi.lci, active_set);
                // open_simple_loops(dbgi.dbg, dbgi.lci, active_set);

                redirect_deleted.clear();
                // erase_if(dbgi.lci, [](const auto &node) {  // MARK
                //     return node.second.id == lci_type::max_id;
                // });
                simplify_simple_paths(dbgi.dbg, dbgi.lci, redirect_deleted,
                                      active_set);
                IF_TIME_CHECKPOINT(!acfg.silent, GraphCleaning, stderr,
                                   "Iteration %d done.", i);
            }
            fix_edges(dbgi.dbg, dbgi.lci, active_set);
        }
    }
};

inline void print_raw_dbg(FILE *out, const dbg_index &dbgi) {
    for (const auto &[lval, edges] : dbgi.dbg) {
        fprintf(out, "S\tC%" PRIu64 "\t%s\tLN:i:%d\tRC:i:%d", lval, "*", 100,
                edges.core_count);
    }

    for (const auto &[lval, edges] : dbgi.dbg) {
        for (const auto &p : edges.prev) {
            fprintf(out, "L\tC%" PRIu64 "\t+\tC%" PRIu64 "\t+\t%dM\n", p, lval,
                    1);
        }
    }
}

__attribute__((format(printf, 2, 3)))
static inline void vec_appendf(std::vector<char> &v, const char *fmt, ...) {
    va_list args, args_copy;
    va_start(args, fmt);
    va_copy(args_copy, args);
    int n = vsnprintf(nullptr, 0, fmt, args);
    va_end(args);
    if (n > 0) {
        size_t old = v.size();
        v.resize(old + static_cast<size_t>(n) + 1);
        vsnprintf(v.data() + old, static_cast<size_t>(n) + 1, fmt, args_copy);
        v.pop_back();  // drop the null terminator
    }
    va_end(args_copy);
}

inline void print_dbg_seqs(FILE *out, const dbg_index &dbgi,
                           int nthreads = 0) {
    if (nthreads <= 0)
        nthreads = std::max(1u, std::thread::hardware_concurrency()) - 1;
    if (nthreads < 1) nthreads = 1;

    // Work item: a pointer into the dbg map (stable while we only read).
    struct WorkItem {
        uint64_t lval;
        const typename dbg_index::bidbg_val *edges;
    };

    using work_q = moodycamel::BlockingConcurrentQueue<WorkItem>;
    work_q work_queue;
    std::atomic<bool> work_done{false};
    std::atomic<int> active_workers{nthreads};
    // ── Worker thread pool ──────────────────────────────────────────
    // Each worker dequeues nodes, converts sequences to ASCII, formats
    // the full S-line + L-lines into a vector<char>, then enqueues it.
    // Format one node (S-line + L-lines) into line_buf.
    auto format_node = [&](std::vector<char> &line_buf, const WorkItem &item) {
        const auto &edges = *item.edges;
        const uint64_t lval = item.lval;
        line_buf.clear();

        // Pre-compute total sequence length.
        size_t sofar = 0;
        for (const auto &s : edges.seq) {
            sofar += s.len;
        }

        // S-line header.
        vec_appendf(line_buf, "S\tC%" PRIu64 "\t", lval);

        // Convert all sequence chunks directly into line_buf.
        {
            size_t seq_start = line_buf.size();
            line_buf.resize(seq_start + sofar + 1);  // +1 for trailing \0
            size_t offset = 0;
            for (const auto &s : edges.seq) {
                s.convert_to_ascii(line_buf.data() + seq_start + offset);
                offset += s.len;
            }
            line_buf.resize(seq_start + sofar);  // trim null terminator
        }

        // S-line trailer.
        vec_appendf(line_buf,
                    "\tLN:i:%zu\tRC:i:%d\tCC:i:%u\tBE:i:%zu\tFE:i:%u\n",
                    sofar, (int)(edges.avg_core_count() * sofar),
                    edges.core_count, edges.prev.size(),
                    edges.forward_edge_count);

        // L-lines for incoming edges.
        for (uint64_t nl : edges.prev) {
            if (dbgi.dbg.contains(nl))
                vec_appendf(line_buf,
                            "L\tC%" PRIu64 "\t+\tC%" PRIu64 "\t+\t%dM\n",
                            nl, lval, 0);
        }
    };


    std::vector<std::thread> workers;
    workers.reserve(nthreads);
    for (int i = 0; i < nthreads; ++i) {
        workers.emplace_back([&] {
            WorkItem item;
            std::vector<char> line_buf;
            line_buf.reserve(1 << 16);

            while (true) {
                if (!work_queue.wait_dequeue_timed(item,
                                                   WAIT_QUEUE_TIMEOUT)) {
                    if (work_done.load()) break;
                    continue;
                }
                format_node(line_buf, item);
                fwrite(line_buf.data(), line_buf.size(), 1, out);
            }

            // Drain remaining work after done flag is set.
            while (work_queue.try_dequeue(item)) {
                const auto &edges = *item.edges;
                const uint64_t lval = item.lval;
                line_buf.clear();

                size_t sofar = 0;
                for (const auto &s : edges.seq) {
                    sofar += s.len;
                }

                vec_appendf(line_buf, "S\tC%" PRIu64 "\t", lval);

                {
                    size_t seq_start = line_buf.size();
                    line_buf.resize(seq_start + sofar + 1);
                    size_t offset = 0;
                    for (const auto &s : edges.seq) {
                        s.convert_to_ascii(line_buf.data() + seq_start +
                                           offset);
                        offset += s.len;
                    }
                    line_buf.resize(seq_start + sofar);
                }

                vec_appendf(
                    line_buf,
                    "\tLN:i:%zu\tRC:i:%d\tCC:i:%u\tBE:i:%zu\tFE:i:%u\n",
                    sofar, (int)(edges.avg_core_count() * sofar),
                    edges.core_count, edges.prev.size(),
                    edges.forward_edge_count);

                for (uint64_t nl : edges.prev) {
                    if (dbgi.dbg.contains(nl))
                        vec_appendf(
                            line_buf,
                            "L\tC%" PRIu64 "\t+\tC%" PRIu64 "\t+\t%dM\n",
                            nl, lval, 0);
                }
                fwrite(line_buf.data(), line_buf.size(), 1, out);

                line_buf.clear();
            }
            active_workers.fetch_sub(1);
        });
    }

    // ── Producer: main thread enqueues work ─────────────────────────
    for (const auto &[lval, edges] : dbgi.dbg) {
        if (edges.core_count == 0 && edges.prev.size() == 0) continue;
        if (edges.seq.size() == 0) continue;

        while (!work_queue.try_enqueue(WorkItem{lval, &edges})) {
            std::this_thread::sleep_for(std::chrono::nanoseconds(10));
        }
    }
    work_done.store(true);

    for (auto &w : workers) w.join();
    // All workers have exited and decremented active_workers to 0;
    // the writer thread will notice and drain the output queue.

}

inline void print_dbg_dot(const dbg_index &dbgi, int i) {
    std::string filename = std::string{"graph"} + std::to_string(i) + ".dot";
    FILE *dot_file = fopen(filename.c_str(), "w");
    if (dot_file) {
        // print_raw_dbg(dot_file, dbgi);
        fprintf(dot_file, "digraph G {\n");
        // Legend
        fprintf(dot_file, "\tlegend [shape=none, label=<\n");
        fprintf(
            dot_file,
            "\t\t<table border=\"0\" cellborder=\"1\" cellspacing=\"0\">\n");
        fprintf(dot_file,
                "\t\t\t<tr><td bgcolor=\"lightgrey\">Node ID</td><td "
                "bgcolor=\"lightgrey\">Forward</td><td "
                "bgcolor=\"lightgrey\">Prev</td><td "
                "bgcolor=\"lightgrey\">Merged count</td><td "
                "bgcolor=\"lightgrey\">Core Count </td></tr>\n");
        fprintf(dot_file, "\t\t</table>\n");
        fprintf(dot_file, "\t>];\n");
        for (auto &node : dbgi.dbg) {
            fprintf(dot_file,
                    "\t%" PRIu64 " [label=\"%" PRIu64 "->%u,%zu,%zu,%d\"];\n",
                    node.first, node.first, node.second.forward_edge_count,
                    node.second.prev.size(), node.second.seq.size(),
                    node.second.core_count);
            for (auto &prev : node.second.prev) {
                fprintf(dot_file, "\t%" PRIu64 " -> %" PRIu64 ";\n", prev,
                        node.first);
            }
        }
        fprintf(dot_file, "}\n");
        fclose(dot_file);
    }
}

#endif  // _GRAPH_H_
