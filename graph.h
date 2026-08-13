#ifndef _GRAPH_H_
#define _GRAPH_H_

#include "gtl/btree.hpp"
#define PREV_WITH_VEC
#define HT_SHARD_BC 8

#include <blockingconcurrentqueue.h>

#include <atomic>
#include <cstdint>
#include <cstdio>
#include <cassert>
#include <cmath>
#include <cstdlib>
#include <gtl/phmap.hpp>
#include <gtl/vector.hpp>
#include <limits>
#include <mutex>
#include <queue>
#include <string>
#include <thread>
#include <vector>
#include <cinttypes>

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
using i64 = int64_t;
using u32 = uint32_t;
using DefaultBF = bloom_interface<bloom_ext::BlockedAtomicBloomFilter<>>;
template <class BFT = DefaultBF>
struct dbg_index_base {
    using BloomFilterType = BFT;
    static constexpr u64 lcindex_increment =
        nextpowerof2((u64)MAXIMUM_WORK_THREADS);
    using mutex_type = std::mutex;
    struct bidbg_val {
        struct counts_t {
            unsigned core_count;
        };

        bidbg_val(unsigned core_count)
            : redirected_flag(false), forward_edge_count(0), counts{core_count},
              seq{}, prev{} {}
        bidbg_val(unsigned core_count, u64 prev)
            : redirected_flag(false), forward_edge_count(0), counts{core_count},
              seq{}, prev{prev} {}

        bool redirected_flag;
        // forward_edge_count is deliberately OUTSIDE the redirect_to union
        // (see RC_REVIEW.md RC1). Adjacent simple-path chains share their
        // boundary node: one chain reads the boundary's forward_edge_count
        // while the neighbouring chain (which owns that node as its seed)
        // redirects it. Keeping fec out of the union means mark_redirected's
        // redirect_to write never aliases the fec a concurrent walker reads,
        // and fec itself is only ever written on representatives — never on a
        // fec!=1 boundary node. It occupies the padding that already followed
        // redirected_flag, so this costs zero extra bytes per node. core_count
        // stays unioned with redirect_to because it is only ever touched by a
        // node's own (single) owning chain, so it cannot race.
        unsigned forward_edge_count;
        union {
            counts_t counts;
            u64 redirect_to;
        };

        gtl::vector<SequenceContainerType> seq;

        bool redirected() const {
            return redirected_flag || (prev.size() == 1 && prev[0] == 0);
        }
        void mark_redirected(u64 to) {
            redirected_flag = true;
            redirect_to = to;
            gtl::vector<SequenceContainerType>().swap(seq);
            gtl::vector<u64>().swap(prev);
            prev.push_back(0);  // 0 is an illegal node id; marks redirected.
        }

        gtl::vector<u64> prev;
        void emplace_prev(u64 val) {
            if (val == 0) return;
            auto it = std::find(prev.begin(), prev.end(), val);
            if (it == prev.end()) {
                prev.push_back(val);
            }
        }

        constexpr double avg_core_count() const {
            // Precondition: live (non-consumed) node. seq.empty() ⟺ the node
            // was absorbed by another during simplification (see
            // mark_redirected) and must be skipped via redirected() before
            // calling this. The assert catches a broken invariant in debug/
            // tests; the release path degrades to 0.0 rather than letting inf
            // propagate into depth comparisons.
            assert(!seq.empty());
            return seq.empty() ? 0.0 : (double)counts.core_count / seq.size();
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
    };

    dbg_map_t dbg;
    gtl::btree_set<u64> active_set;

    // merge_info_map_t merge_map;
    BloomFilterType bf;

    dbg_index_base(double false_positive_rate, uint64_t projected_count,
                   uint64_t seed)
        : bf{false_positive_rate, projected_count, seed} {}

    // dbg_index(const dbg_index &other)
    static constexpr u64 lcmer_id(const LCmer &l) { return l.value(); }

    lci_type get(const LCmer &l) const {
        u64 id = lcmer_id(l);
        gtl::LockableBaseImpl<gtl::NullMutex>::DoNothing no_lock;
        auto &dbg_nc = const_cast<dbg_map_t &>(dbg);
        auto h = dbg_nc.hash(id);
        auto it = dbg_nc.find(id, h, no_lock);
        if (it == dbg_nc.end()) {
            return lci_type::empty();
        }
        for (int hop = 0;
             hop < 128 && it != dbg_nc.end() && it->second.redirected();
             ++hop) {
            const u64 next = it->second.redirect_to;
            if (next == 0 || next == id) {
                return lci_type::empty();
            }
            id = next;
            h = dbg_nc.hash(id);
            it = dbg_nc.find(id, h, no_lock);
        }
        return it == dbg_nc.end() ? lci_type::empty() : lci_type{id};
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
    u64 emplace_nobloomfilter(const LCmer &l, u64 &idx, int tid) {
        (void)idx;
        (void)tid;
        return lcmer_id(l);
    }
    u64 emplace(const LCmer &l, u64 &idx, int tid) {
        (void)idx;
        (void)tid;
        const auto new_hash = std::hash<LCmer>{}(l);
        if (!bf.contains(new_hash)) {
            bf.insert(new_hash);
            return 0;
        }
        return lcmer_id(l);
    }
    gtl::vector<u64> shortest_path_bfs(u64 source, u64 target) const {
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
            if (it == dbg.end()) {
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
    using redir_map = gtl::parallel_flat_hash_map_d<u64, u64, 4, std::mutex>;
    struct avg_core_count_stats {
        double baseline = 0.0;
    };
    static double median_of_sorted(const gtl::vector<double> &vals) {
        if (vals.empty()) return 0.0;
        const size_t mid = vals.size() / 2;
        if (vals.size() % 2 == 1) return vals[mid];
        return (vals[mid - 1] + vals[mid]) / 2.0;
    }

    static avg_core_count_stats collect_active_avg_core_count_stats(
        const dbg_map_t &dbg, const gtl::btree_set<u64> &active_set) {
        avg_core_count_stats stats;
        constexpr size_t kMinDepthSample = 16;
        struct depth_len_sample {
            double depth;
            size_t len_proxy;
            bool topology_stable;
        };
        gtl::vector<depth_len_sample> samples;
        samples.reserve(active_set.size());
        gtl::vector<size_t> lengths;
        lengths.reserve(active_set.size());

        for (u64 node : active_set) {
            auto it = dbg.find(node);
            if (it == dbg.end() || it->second.redirected()) continue;
            if (it->second.seq.empty()) continue;
            const double d = it->second.avg_core_count();
            if (d <= 0.0) continue;

            const size_t len_proxy = it->second.seq.size();
            const bool topology_stable =
                (it->second.prev.size() <= 1 &&
                 it->second.forward_edge_count <= 1);
            samples.push_back({d, len_proxy, topology_stable});
            lengths.push_back(len_proxy);
        }
        if (samples.empty()) return stats;

        std::sort(lengths.begin(), lengths.end());
        const size_t len_q75 = lengths[(lengths.size() - 1) * 3 / 4];

        auto collect_depths = [&](auto pred) {
            gtl::vector<double> out;
            out.reserve(samples.size());
            for (const auto &s : samples) {
                if (pred(s)) out.push_back(s.depth);
            }
            return out;
        };

        gtl::vector<double> depths = collect_depths([&](const auto &s) {
            return s.topology_stable && s.len_proxy >= len_q75;
        });
        if (depths.size() < kMinDepthSample) {
            depths = collect_depths(
                [&](const auto &s) { return s.len_proxy >= len_q75; });
        }
        if (depths.size() < kMinDepthSample) {
            depths = collect_depths([](const auto &) { return true; });
        }

        std::sort(depths.begin(), depths.end());
        const auto qidx = [n = depths.size()](double q) -> size_t {
            const size_t idx =
                static_cast<size_t>(q * static_cast<double>(n - 1));
            return std::min(idx, n - 1);
        };
        const double raw_q1 = depths[qidx(0.25)];
        const double raw_q3 = depths[qidx(0.75)];
        const double iqr = raw_q3 - raw_q1;
        const double lower = raw_q1 - 1.5 * iqr;
        const double upper = raw_q3 + 1.5 * iqr;

        gtl::vector<double> filtered;
        filtered.reserve(depths.size());
        for (double d : depths) {
            if (d >= lower && d <= upper) filtered.push_back(d);
        }
        const auto &pool = filtered.empty() ? depths : filtered;
        stats.baseline = median_of_sorted(pool);
        return stats;
    }

    static auto attempt_simplify_one_path_lockless(dbg_map_t &dbg,
                                                   redir_map &redirect_deleted,
                                                   u64 id,
                                                   auto starting_from) -> bool {
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
               prev_it->first != current->first &&  // never step onto a self-
                                                    // loop (would merge the
                                                    // node into itself forever)
               prev_it->second.forward_edge_count == 1) {
            if (prev_it != dbg.end()) {
                cc_sum += current->second.counts.core_count;
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
        if (!simplified) {
            // Degenerate: nothing was merged (merged empty, cc_sum == 0,
            // current == seed). Mutating the seed here would be a
            // value-identical self-store that races neighbor chains'
            // stop-condition reads of forward_edge_count (CONFIRMATION_C1.md)
            // — the only write that could cross chain boundaries. Skip it.
            return false;
        }
        // simplified implies the walk moved, i.e. current != seed.
        redirect_deleted[id] = current->first;

        current->second.forward_edge_count = first_fec;
        current->second.counts.core_count += cc_sum;
        for (auto it = merged.rbegin(); it != merged.rend(); ++it) {
            // for (auto &seq : (*it)->second.seq)
            gtl::vector<SequenceContainerType> &seqs = (*it)->second.seq;
            for (auto sit = seqs.begin(); sit != seqs.end(); ++sit) {
                current->second.seq.emplace_back(std::move(*sit));
            }
            seqs.clear();
            (*it)->second.mark_redirected(current->first);
        }

        return true;
    };

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
                        if (prev_id == item.first) {
                            continue;  // never seed a self-loop
                        }
                        auto prev_hash = dbg.hash(prev_id);

                        auto prev_it = dbg.find(prev_id, prev_hash, no_lock);
                        // Seed only fec==1 predecessors here. A fec!=1
                        // predecessor is already a type-1 seed (above), so
                        // seeding it again would have two threads walk the
                        // same chain -> RC1 race + double merge. With this
                        // guard, type-1 (fec!=1) and type-2 (fec==1) seeds are
                        // disjoint and every chain is owned by one thread.
                        // (See RC_REVIEW.md RC1.)
                        if (prev_it != dbg.end() &&
                            prev_it->second.prev.size() == 1 &&
                            prev_it->second.forward_edge_count == 1) {
                            to_simplify_from[tid].emplace_back(prev_id,
                                                               prev_it);
                        }
                    }
                }
            },
            threads);
        return to_simplify_from;
    }
    static void simplify_simple_parallel_sweep(
        dbg_index::dbg_map_t &dbg,
        const gtl::vector<
            gtl::vector<std::pair<u64, dbg_index::dbg_map_t::iterator>>>
            &to_simplify,
        redir_map &redirect_deleted, auto &active_set) {
        // Chains are node-disjoint (see simplify_simple_parallel_mark's seed
        // guard), so each thread mutates only the nodes of its own chains with
        // no synchronization. forward_edge_count in particular is written only
        // on representatives strictly inside the writer's own chain: the
        // degenerate-path gate in attempt_simplify_one_path_lockless (the C1
        // fix) removed the one cross-chain self-store, so the sweep is
        // data-race free under the disjointness invariant (CONFIRMATION_C1.md).
        // The shared redirect_deleted is a mutexed parallel map; active_set is
        // pruned after the join (below) rather than via a per-node id list.
        // (See RC_REVIEW.md RC1.)
        gtl::vector<std::thread> threads;
        for (size_t i = 0; i < to_simplify.size(); ++i) {
            threads.emplace_back(
                [&](int i) {
                    for (const auto &[node, it] : to_simplify[i]) {
                        attempt_simplify_one_path_lockless(
                            dbg, redirect_deleted, node, it);
                    }
                },
                i);
        }

        for (auto &thread : threads) {
            thread.join();
        }

        // Prune absorbed nodes from active_set by state, not by a stored id
        // list: every node merged into a representative was marked
        // redirected() during the sweep, and nodes no longer in dbg are gone.
        for (auto it = active_set.begin(); it != active_set.end();) {
            auto nit = dbg.find(*it);
            if (nit == dbg.end() || nit->second.redirected()) {
                it = active_set.erase(it);
            } else {
                ++it;
            }
        }

        for (auto &node : dbg) {
            for (size_t idx = 0; idx < node.second.prev.size(); ++idx) {
                auto &p = node.second.prev[idx];
                auto it = redirect_deleted.find(p);
                if (it != redirect_deleted.end()) {
                    p = it->second;
                }
                if (p == node.first) {
                    // Cycle collapse: the redirect rewrite turned this prev
                    // entry into a self-reference. The self-edge is spurious
                    // bookkeeping (the collapsed cycle edge) — drop it and
                    // keep in/out bookkeeping consistent. (At most one such
                    // entry can exist: emplace_prev dedups.)
                    if (node.second.forward_edge_count > 0) {
                        node.second.forward_edge_count--;
                    }
                    node.second.prev.erase(node.second.prev.begin() + idx);
                    --idx;
                }
            }
        }
    }
    // Unified simple-path simplification: discover node-disjoint chains in
    // parallel (mark), then merge them lock-free (sweep). Used for both the
    // initial pass and every cleaning iteration; threads==1 just means a single
    // partition. Replaces the old divergent single-threaded simplify_simple_paths
    // and the dead attempt_simplify_one_path. (See RC_REVIEW.md RC1.)
    static void simplify_simple_paths(dbg_index::dbg_map_t &dbg,
                                      redir_map &redirect_deleted,
                                      auto &active_set, int threads) {
        auto to_simplify = simplify_simple_parallel_mark(dbg, threads);
        simplify_simple_parallel_sweep(dbg, to_simplify, redirect_deleted,
                                       active_set);
    }

    static void cut_tips(dbg_map_t &dbg, int min_core_depth,
                         size_t max_tip_length, auto &active_set) {
        const auto prune_prev = [&dbg](auto current) {
            auto &pvec = current->second.prev;
            pvec.erase(std::remove_if(pvec.begin(), pvec.end(),
                                      [&dbg](u64 pid) {
                                          auto it = dbg.find(pid);
                                          return it == dbg.end() ||
                                                 it->second.redirected();
                                      }),
                       pvec.end());
            std::sort(pvec.begin(), pvec.end());
        };

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
                    if (prev_it != dbg.end() && !prev_it->second.redirected()) {
                        prev_it->second.forward_edge_count--;
                    }
                    dbg.erase(current);  // MARK
                }
            } else {
                gtl::vector<u64> prev_ids(current->second.prev.begin(),
                                          current->second.prev.end());
                std::sort(prev_ids.begin(), prev_ids.end());
                for (u64 pid : prev_ids) {
                    auto prev_it = dbg.find(pid);
                    if (prev_it != dbg.end() && !prev_it->second.redirected()) {
                        if (prev_it->second.avg_core_count() < min_core_depth &&
                            prev_it->second.prev.size() == 0 &&
                            prev_it->second.forward_edge_count == 1) {
                            dbg.erase(prev_it);  // MARK
                        }
                    }
                }

                current = dbg.find(node);
                if (current == dbg.end() || current->second.redirected()) {
                    continue;
                }
                prune_prev(current);

                prev_ids.assign(current->second.prev.begin(),
                                current->second.prev.end());
                for (u64 pid : prev_ids) {
                    auto prev_it = dbg.find(pid);
                    if (prev_it != dbg.end() && !prev_it->second.redirected()) {
                        if (prev_it->second.avg_core_count() <
                                current->second.avg_core_count() &&
                            prev_it->second.seq.size() < max_tip_length &&
                            current->second.seq.size() > 2 * max_tip_length &&
                            prev_it->second.prev.size() == 0 &&
                            prev_it->second.forward_edge_count == 1) {
                            dbg.erase(prev_it);  // MARK
                        }
                    }
                }

                current = dbg.find(node);
                if (current == dbg.end() || current->second.redirected()) {
                    continue;
                }
                prune_prev(current);
            }
        }
    }

    // Live forward adjacency: successor lists derived from prev edges,
    // sorted per node for deterministic arm walking. Dead/redirected nodes
    // are excluded.
    static gtl::flat_hash_map<u64, gtl::vector<u64>> build_live_fwd(
        const dbg_map_t &dbg) {
        gtl::flat_hash_map<u64, gtl::vector<u64>> fwd;
        for (const auto &[id, node] : dbg) {
            if (node.redirected()) continue;
            for (u64 p : node.prev) {
                auto pit = dbg.find(p);
                if (pit != dbg.end() && !pit->second.redirected()) {
                    fwd[p].push_back(id);
                }
            }
        }
        for (auto &[_, succs] : fwd) {
            std::sort(succs.begin(), succs.end());
        }
        return fwd;
    }

    // Walk a simple forward path from `start` until a reconvergence node
    // (>= 2 predecessors) is reached; return it. Abort (0) on non-simple
    // nodes, divergence, dead/stale edges, cycles, or the node cap. The path
    // is collected excluding the sink.
    static u64 walk_arm_path(
        const dbg_map_t &dbg,
        const gtl::flat_hash_map<u64, gtl::vector<u64>> &fwd, u64 start,
        gtl::vector<u64> &path) {
        path.clear();
        u64 cur = start;
        constexpr size_t kMaxArmNodes = 64;
        for (size_t step = 0; step < kMaxArmNodes; ++step) {
            auto cit = dbg.find(cur);
            if (cit == dbg.end() || cit->second.redirected()) return 0;
            if (cit->second.prev.size() != 1) return 0;
            path.push_back(cur);
            auto nxt_it = fwd.find(cur);
            if (nxt_it == fwd.end() || nxt_it->second.empty()) return 0;
            u64 next = 0;
            for (u64 cand : nxt_it->second) {
                auto nit = dbg.find(cand);
                if (nit == dbg.end() || nit->second.redirected()) {
                    return 0;  // stale edge from an earlier mutation
                }
                if (next != 0) return 0;  // divergence
                next = cand;
            }
            if (next == start) return 0;  // cycle back onto the arm start
            auto nit = dbg.find(next);
            if (nit->second.prev.size() >= 2) return next;  // sink
            cur = next;
        }
        return 0;  // cap exceeded
    }

    // Arm stats: core count and minimum avg_core_count along a path. Dead
    // nodes anywhere on the path abort with {0, 0.0}.
    static std::pair<size_t, double> arm_depth_len(
        const dbg_map_t &dbg, const gtl::vector<u64> &path) {
        size_t len = 0;
        double min_depth = std::numeric_limits<double>::infinity();
        for (u64 nid : path) {
            auto nit = dbg.find(nid);
            if (nit == dbg.end() || nit->second.redirected()) {
                return {0, 0.0};
            }
            len += nit->second.seq.size();
            min_depth = std::min(min_depth, nit->second.avg_core_count());
        }
        return {len, min_depth};
    }

    // Bubble popping (strict fork-2 -> merge-2 bubbles).
    // A bubble is two node-disjoint simple paths S -> ... -> T between a
    // 2-successor fork S and a 2-predecessor merge T (S != T). Post-simplify
    // arms are single nodes, but arms are walked explicitly so nested-branch
    // candidates are rejected instead of misread. The shallower arm (minimum
    // avg_core_count along the arm) is popped when it is at least
    // bubble_depth_epsilon times shallower than the other and the two arm
    // lengths (core counts) are within bubble_length_bound of each other.
    // Strict 2/2 topology is required: among 3+ parallel copies (SD/rDNA
    // tangles) nothing is popped -- deleting one of many real copies is a
    static size_t pop_bubbles(dbg_map_t &dbg, double bubble_depth_epsilon,
                              double bubble_length_bound,
                              size_t max_popped_arm_len, bool dump_pops,
                              bool dump_loops, auto &active_set) {
        auto fwd = build_live_fwd(dbg);

        size_t popped = 0;
        for (u64 sid : active_set) {
            auto sit = dbg.find(sid);
            if (sit == dbg.end() || sit->second.redirected()) continue;
            if (sit->second.forward_edge_count < 2) continue;
            auto ffit = fwd.find(sid);
            if (ffit == fwd.end() || ffit->second.size() != 2) continue;

            gtl::vector<u64> path0, path1;
            const u64 t0 = walk_arm_path(dbg, fwd, ffit->second[0], path0);
            if (t0 == 0) continue;
            const u64 t1 = walk_arm_path(dbg, fwd, ffit->second[1], path1);
            if (t1 == 0) continue;

            // Loop bubble: an arm closes back onto the fork S itself (its
            // walk sinks at S). Fork-join candidates never do (T != S).
            // Loops are NOT popped: depth-ratio popping was measured and
            // rejected (75 pops -> 28 misassemblies on chr18, 19 on .fa;
            // shallow loop arms are real low-coverage repeat units, not
            // sequencing error). Loops are handled by unroll_loops
            // (sequence duplication by predicted copy number) instead.
            if (t0 == sid || t1 == sid) {
                if (dump_loops) {
                    const auto &loop_path = (t0 == sid) ? path0 : path1;
                    const auto &other_path = (t0 == sid) ? path1 : path0;
                    auto [lenL, dL] = arm_depth_len(dbg, loop_path);
                    auto [lenO, dO] = arm_depth_len(dbg, other_path);
                    std::fprintf(stderr,
                                 "[loop] S=%llu loop_arm=%d lenL=%zu "
                                 "lenO=%zu dL=%.3f dO=%.3f S_depth=%.3f "
                                 "S_prev=%zu t0=%llu t1=%llu\n",
                                 (unsigned long long)sid, (t0 == sid) ? 0 : 1,
                                 lenL, lenO, dL, dO,
                                 sit->second.avg_core_count(),
                                 sit->second.prev.size(),
                                 (unsigned long long)t0,
                                 (unsigned long long)t1);
                }
                continue;
            }
            if (t1 != t0) continue;  // fork-join: arms must reconverge

            auto tit = dbg.find(t0);
            if (tit == dbg.end() || tit->second.redirected() ||
                tit->second.prev.size() != 2) {
                continue;  // strict 2-merge only
            }

            auto [len0, d0] = arm_depth_len(dbg, path0);
            auto [len1, d1] = arm_depth_len(dbg, path1);
            if (len0 == 0 || len1 == 0) continue;
            const size_t lmax = std::max(len0, len1);
            const size_t lmin = std::min(len0, len1);
            if ((double)lmax > (double)lmin * bubble_length_bound) {
                continue;  // lengths too divergent
            }
            const double dhi = std::max(d0, d1);
            const double dlo = std::min(d0, d1);
            if (dhi <= 0.0 || dlo / dhi >= bubble_depth_epsilon) continue;

            const size_t dead_len = (d0 <= d1) ? len0 : len1;
            if (max_popped_arm_len > 0 &&
                dead_len > max_popped_arm_len) {
                continue;  // arm too long to be a sequencing-error bubble
            }
            if (dump_pops) {
                std::fprintf(stderr,
                             "[pop] S=%llu T=%llu len0=%zu len1=%zu "
                             "d0=%.3f d1=%.3f\n",
                             (unsigned long long)sid, (unsigned long long)t0,
                             len0, len1, d0, d1);
            }

            // Pop the shallower arm; on a tie, path0 goes (deterministic).
            const gtl::vector<u64> &dead = (d0 <= d1) ? path0 : path1;
            for (u64 nid : dead) {
                dbg.erase(nid);
                active_set.erase(nid);
            }
            auto &tprev = tit->second.prev;
            const u64 sink_in = dead.back();
            auto pit = std::find(tprev.begin(), tprev.end(), sink_in);
            if (pit != tprev.end()) tprev.erase(pit);
            sit->second.forward_edge_count--;
            ++popped;
        }
        return popped;
    }
    // Loop handling is its own mechanism, distinct from fork-join bubble
    // popping: a loop arm is never deleted. Its sequence is duplicated into
    // the fork node as many times as the predicted copy number implies,
    // unrolling the cycle in place. Two shapes are handled:
    //   A -> B (self-loop) -> C   : B -> B edge; B.seq duplicated.
    //   A -> (S -> arm -> S) -> D : a walked simple arm closing back at its
    //                               fork S; (arm_seq + S.seq) duplicated.
    // Copy numbers come from depth ratios against a haploid reference (the
    // exit path's depth or the genome-wide baseline), and unrolling only
    // happens in haploid context (reference close to baseline), so SD/rDNA
    // pileup loops -- real divergent copies -- are left alone.
    static size_t unroll_loops(dbg_map_t &dbg, gtl::btree_set<u64> &active_set,
                               const avg_core_count_stats &depth_stats,
                               size_t max_loop_arm_len, bool dump_unroll) {
        const double baseline = depth_stats.baseline;
        if (baseline <= 0.0) return 0;
        constexpr double kHaploidRefMin = 0.5;
        constexpr double kHaploidRefMax = 1.5;
        constexpr int kMaxTurns = 8;
        auto fwd = build_live_fwd(dbg);
        gtl::vector<u64> path0, path1;
        gtl::vector<u64> to_rm;
        size_t unrolled = 0;

        for (u64 sid : active_set) {
            auto sit = dbg.find(sid);
            if (sit == dbg.end() || sit->second.redirected()) continue;
            if (sit->second.forward_edge_count < 2) continue;
            auto ffit = fwd.find(sid);

            // Shape 1: direct self-loop S -> S.
            bool self_loop = false;
            double ref = baseline;
            for (u64 c : ffit->second) {
                auto cit = dbg.find(c);
                if (cit == dbg.end() || cit->second.redirected()) continue;
                if (c == sid) {
                    self_loop = true;
                    continue;
                }
                const double d = cit->second.avg_core_count();
                if (d > 0.0) {
                    ref = d;
                    break;
                }
            }
            if (self_loop) {
                const double s_depth = sit->second.avg_core_count();
                const double ref_ratio = ref / baseline;
                if (ref_ratio >= kHaploidRefMin &&
                    ref_ratio <= kHaploidRefMax && s_depth > 0.0) {
                    const int turns = std::clamp(
                        static_cast<int>(std::llround(s_depth / ref)), 1,
                        kMaxTurns);
                    if (turns >= 2) {
                        if (dump_unroll) {
                            std::fprintf(stderr,
                                         "[unroll-self] S=%llu turns=%d "
                                         "depth=%.3f ref=%.3f\n",
                                         (unsigned long long)sid, turns,
                                         s_depth, ref);
                        }
                        gtl::vector<SequenceContainerType> orig;
                        for (const auto &s : sit->second.seq) {
                            orig.push_back(SequenceContainerType::copy_of(s));
                        }
                        for (int k = 1; k < turns; ++k) {
                            for (const auto &s : orig) {
                                sit->second.seq.push_back(
                                    SequenceContainerType::copy_of(s));
                            }
                        }
                        sit->second.forward_edge_count--;
                        auto &pv = sit->second.prev;
                        auto pit = std::find(pv.begin(), pv.end(), sid);
                        if (pit != pv.end()) pv.erase(pit);
                        ++unrolled;
                    }
                }
                continue;
            }

            // Shape 2: strict fork; exactly one arm closes back at S.
            if (ffit->second.size() != 2) continue;
            const u64 t0 = walk_arm_path(dbg, fwd, ffit->second[0], path0);
            if (t0 == 0) continue;
            const u64 t1 = walk_arm_path(dbg, fwd, ffit->second[1], path1);
            if (t1 == 0) continue;
            const bool loop0 = (t0 == sid);
            const bool loop1 = (t1 == sid);
            if (loop0 == loop1) continue;  // no loop, or both arms loop
            const gtl::vector<u64> &loop_path = loop0 ? path0 : path1;
            const gtl::vector<u64> &exit_path = loop0 ? path1 : path0;
            auto [lenL, dL] = arm_depth_len(dbg, loop_path);
            auto [lenE, dE] = arm_depth_len(dbg, exit_path);
            if (lenL == 0 || lenE == 0) continue;
            if (max_loop_arm_len > 0 && lenL > max_loop_arm_len) {
                continue;  // loop unit too large: repeat architecture
            }
            const double ref2 = dE > 0.0 ? dE : baseline;
            const double ref_ratio2 = ref2 / baseline;
            if (ref_ratio2 < kHaploidRefMin || ref_ratio2 > kHaploidRefMax) {
                continue;  // not haploid context: pileup copy divergence
            }
            // Copy-number prediction: the loop arm must be a clean integer
            // multiple of the haploid reference (its own traversal count),
            // and S must match one extra traversal. Depth averages are noisy;
            // a sloppy fit means the "loop" is a collapsed multi-copy repeat,
            // not a predictable cycle -- leave it alone.
            const double s_depth = sit->second.avg_core_count();
            const int t_int = static_cast<int>(std::llround(dL / ref2));
            if (t_int < 2 || t_int > kMaxTurns) continue;
            const double loop_err = std::abs(dL - t_int * ref2) / ref2;
            const double curr_err =
                std::abs(s_depth - (t_int + 1) * ref2) / ref2;
            if (loop_err > 0.25 || curr_err > 0.25) continue;
            const int turns = t_int;
            if (dump_unroll) {
                std::fprintf(stderr,
                             "[unroll] S=%llu turns=%d arm_len=%zu dL=%.3f "
                             "ref=%.3f S_depth=%.3f\n",
                             (unsigned long long)sid, turns, lenL, dL, ref2,
                             sit->second.avg_core_count());
            }
            // Unroll: S.seq += (arm_seq + S.seq) x turns.
            gtl::vector<SequenceContainerType> s_orig;
            for (const auto &s : sit->second.seq) {
                s_orig.push_back(SequenceContainerType::copy_of(s));
            }
            gtl::vector<SequenceContainerType> arm_seq;
            for (u64 nid : loop_path) {
                auto nit = dbg.find(nid);
                for (const auto &s : nit->second.seq) {
                    arm_seq.push_back(SequenceContainerType::copy_of(s));
                }
            }
            for (int k = 0; k < turns; ++k) {
                for (const auto &s : arm_seq) {
                    sit->second.seq.push_back(
                        SequenceContainerType::copy_of(s));
                }
                for (const auto &s : s_orig) {
                    sit->second.seq.push_back(
                        SequenceContainerType::copy_of(s));
                }
            }
            // Remove the loop arm; rewire S. dbg.erase is safe here (we
            // iterate active_set, not dbg), but active_set erases are
            // deferred: erasing mid-iteration corrupts the btree walk.
            for (u64 nid : loop_path) {
                dbg.erase(nid);
                to_rm.push_back(nid);
            }
            auto &spv = sit->second.prev;
            auto pit = std::find(spv.begin(), spv.end(), loop_path.back());
            sit->second.forward_edge_count--;
            ++unrolled;
        }
        for (u64 nid : to_rm) {
            active_set.erase(nid);
        }
        return unrolled;
    }

    static void cut_low_support_edges(dbg_map_t &dbg, auto &active_set) {
        // Death by thousand cuts
        for (auto &node : active_set) {
            auto current = dbg.find(node);
            if (current == dbg.end() || current->second.redirected()) {
                // Skip consumed husks: a redirected node's seq is empty and
                // its active union member is redirect_to, not counts.
                continue;
            }
            if (current->second.avg_core_count() < 2 &&
                current->second.prev.size() == 1 &&
                current->second.forward_edge_count == 1) {
                auto prev_it = dbg.find(current->second.prev[0]);
                if (prev_it != dbg.end())
                    prev_it->second.forward_edge_count--;
                dbg.erase(current);
            }
        }
    }
    static void fix_edges(dbg_map_t &dbg) {
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
    }
    static void run(dbg_index &dbgi, const asm_config &acfg) {
        auto cfg = acfg.graph_cleaning;
        INIT_TIMERS(GraphCleaning);
        redir_map redirect_deleted;
        decltype(dbgi.active_set) &active_set = dbgi.active_set;

        for (const auto &node : dbgi.dbg) {
            active_set.insert(node.first);
        }
        if (cfg.iterations > -1) {
            simplify_simple_paths(dbgi.dbg, redirect_deleted, active_set,
                                  acfg.threads.work);
            IF_TIME_CHECKPOINT(!acfg.silent, GraphCleaning, stderr,
                               "Initial Simplification");
            cut_low_support_edges(dbgi.dbg, active_set);
            const auto initial_depth_stats =
                collect_active_avg_core_count_stats(dbgi.dbg, active_set);
            const bool dump_pops =
                std::getenv("LASM_DUMP_POPS") != nullptr;
            const bool dump_loops =
                std::getenv("LASM_DUMP_LOOPS") != nullptr;
            const bool dump_unroll =
                std::getenv("LASM_DUMP_UNROLL") != nullptr;
            for (int i = 1; i < cfg.iterations; ++i) {
                fix_edges(dbgi.dbg);
                cut_tips(dbgi.dbg, 3, cfg.max_tip_length, active_set);
                // Strict fork-2/merge-2 bubble popping (see pop_bubbles).
                const size_t popped = pop_bubbles(
                    dbgi.dbg, cfg.min_bubble_depth_difference_percentage,
                    cfg.min_bubble_length_difference,
                    cfg.max_bubble_arm_length, dump_pops, dump_loops,
                    active_set);
                if (popped > 0)
                    std::fprintf(stderr,
                                 "[graph-clean] popped %zu bubbles\n",
                                 popped);
                if (cfg.enable_loop_opening) {
                    const size_t unrolled = unroll_loops(
                        dbgi.dbg, active_set, initial_depth_stats,
                        cfg.max_bubble_arm_length, dump_unroll);
                    if (unrolled > 0)
                        std::fprintf(stderr,
                                     "[graph-clean] unrolled %zu loops\n",
                                     unrolled);
                }
                redirect_deleted.clear();
                // erase_if(dbgi.lci, [](const auto &node) {  // MARK
                //     return node.second.id == lci_type::max_id;
                // });
                // Iterations run the simplify serially (single partition) --
                // the committed behavior. Parallelizing them is a future perf
                // win once RC3 (a residual cross-chain fec race in the parallel
                // sweep -- see RC_REVIEW.md) is resolved. The initial pass above
                // keeps the committed parallelism (threads.work).
                simplify_simple_paths(dbgi.dbg, redirect_deleted, active_set, 1);
                IF_TIME_CHECKPOINT(!acfg.silent, GraphCleaning, stderr,
                                   "Iteration %d done.", i);
            }
            fix_edges(dbgi.dbg);
        }
    }
};

inline void print_raw_dbg(FILE *out, const dbg_index &dbgi) {
    for (const auto &[lval, edges] : dbgi.dbg) {
        fprintf(out, "S\tC%" PRIu64 "\t%s\tLN:i:%d\tRC:i:%d", lval, "*", 100,
                edges.counts.core_count);
    }

    for (const auto &[lval, edges] : dbgi.dbg) {
        for (const auto &p : edges.prev) {
            fprintf(out, "L\tC%" PRIu64 "\t+\tC%" PRIu64 "\t+\t%dM\n", p, lval,
                    1);
        }
    }
}

__attribute__((format(printf, 2, 3))) static inline void vec_appendf(
    std::vector<char> &v, const char *fmt, ...) {
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

inline void print_dbg_seqs(FILE *out, const dbg_index &dbgi, int nthreads = 0) {
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
                    "\tLN:i:%zu\tRC:i:%d\tCC:i:%u\tBE:i:%zu\tFE:i:%u\n", sofar,
                    (int)(edges.avg_core_count() * sofar),
                    edges.counts.core_count, edges.prev.size(),
                    edges.forward_edge_count);

        // L-lines for incoming edges.
        for (uint64_t nl : edges.prev) {
            if (dbgi.dbg.contains(nl))
                vec_appendf(line_buf,
                            "L\tC%" PRIu64 "\t+\tC%" PRIu64 "\t+\t%dM\n", nl,
                            lval, 0);
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
                if (!work_queue.wait_dequeue_timed(item, WAIT_QUEUE_TIMEOUT)) {
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

                vec_appendf(line_buf,
                            "\tLN:i:%zu\tRC:i:%d\tCC:i:%u\tBE:i:%zu\tFE:i:%u\n",
                            sofar, (int)(edges.avg_core_count() * sofar),
                            edges.counts.core_count, edges.prev.size(),
                            edges.forward_edge_count);

                for (uint64_t nl : edges.prev) {
                    if (dbgi.dbg.contains(nl))
                        vec_appendf(line_buf,
                                    "L\tC%" PRIu64 "\t+\tC%" PRIu64
                                    "\t+\t%dM\n",
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
        if (edges.counts.core_count == 0 && edges.prev.size() == 0) continue;
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
                    node.first, node.first,
                    node.second.forward_edge_count,
                    node.second.prev.size(), node.second.seq.size(),
                    node.second.counts.core_count);
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
