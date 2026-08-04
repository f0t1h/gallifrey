#include <blockingconcurrentqueue.h>
#include <fcntl.h>
#include <htslib/bgzf.h>
#include <htslib/hfile.h>
#include <htslib/hts.h>
#include <htslib/sam.h>
#include <htslib/thread_pool.h>
#include <stdio.h>

#include <algorithm>
#include <chrono>
#include <cstdarg>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <filesystem>
#include <gtl/btree.hpp>
#include <gtl/phmap.hpp>
#include <gtl/phmap_fwd_decl.hpp>
#include <gtl/vector.hpp>
#include <iostream>
#include <map>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

#include "asm.h"
#include "config.h"
#include "ext/doctest.h"
#include "graph.h"
#include "layout_graph.h"
#include "lcmer.h"
#include "no_destructor.h"
#include "lcptools/lcptools_ho.h"
#include "syncmer.h"

using std::string;
enum class input_type { SAM, BAM, FASTX, FASTXGZ, FASTXBZ2, UNKNOWN };

// OLC end-skip policy: cores skipped at each end of a DBG node before
// overlap-anchoring (see find_lcmer_ranges). Written as function-like macros of
// the OLC lcmer length `l` so they scale with --olc-lcmer-l and can be retuned
// at compile time without editing source, e.g.:
//   cmake -DCMAKE_CXX_FLAGS='-DOLC_SKIP_FIRST(l)=0 -DOLC_SKIP_LAST(l)=((l))'
// The old (l, 2l) required >=4*l cores and forced ~29% of good nodes (chr22)
// into singletons; the (1, l/2) default drops singleton contigs 87.0%->84.5%
// and lifts N50 +13% with no extra misassemblies. A node needs
// (l + OLC_SKIP_FIRST(l) + OLC_SKIP_LAST(l)) cores to anchor any overlap.
#ifndef OLC_SKIP_FIRST
#define OLC_SKIP_FIRST(l) (1)
#endif
#ifndef OLC_SKIP_LAST
#define OLC_SKIP_LAST(l) ((l) / 2)
#endif

// Prototype: also keep one-directional (non-reciprocal) overlaps. When A finds
// B but B never finds A, the reverse overlap normally gets dropped. Enabling
// this synthesizes an approximate reverse span so those extra nodes can still
// join a layout. Off by default; build with -DOLC_NONRECIPROCAL=1 to measure.
#ifndef OLC_NONRECIPROCAL
#define OLC_NONRECIPROCAL 0
#endif

// Saturating atomic counting bloom for --thread-bridge read-walk dedup: caps how
// many reads are retained per bridge key so kept walks are bounded by
// O(distinct bridges x cap), independent of sequencing depth. Collisions only
// over-count (retain slightly fewer) which is safe.
struct CountingBloom {
    std::vector<std::atomic<uint8_t>> cells;
    uint64_t mask;
    explicit CountingBloom(unsigned log2 = 22)
        : cells(size_t(1) << log2), mask((uint64_t(1) << log2) - 1) {}
    // Increment key's counter (saturating); return the count BEFORE increment.
    uint8_t bump(uint64_t key) {
        key ^= key >> 33;
        key *= 0xff51afd7ed558ccdULL;
        key ^= key >> 33;
        std::atomic<uint8_t> &c = cells[key & mask];
        uint8_t v = c.load(std::memory_order_relaxed);
        while (v < 255 &&
               !c.compare_exchange_weak(v, uint8_t(v + 1),
                                        std::memory_order_relaxed)) {
        }
        return v;
    }
};

// Stable 64-bit identity of a read qname (FNV-1a). Used for distinct-read bridge
// voting so retained walks never store the qname string — retained-walk memory is
// the dominant cost at whole-genome scale (~100GB reads), so per-walk heap strings
// are not affordable.
static inline u64 hash_qname(const char *s) {
    u64 h = 1469598103934665603ULL;
    for (; *s; ++s) {
        h ^= static_cast<unsigned char>(*s);
        h *= 1099511628211ULL;
    }
    return h;
}

// A retained branch-crossing read walk, stored structure-of-arrays for density:
// `nodes[k]` is the raw DBG id of emitted fragment k (0 = dropped node), and the
// fragment spans the half-open query interval [query_pos[k], query_pos[k+1]).
// Fragments are contiguous, so one boundary per step (+1 trailing) replaces the
// old start/end pair. No qname string is kept — only its hash.
struct retained_read {
    u64 qname_hash = 0;
    gtl::vector<u64> nodes;
    gtl::vector<uint32_t> query_pos;  // size == nodes.size() + 1
};

struct thread_bridge_read_observation {
    u64 qname = 0;  // qname hash (read identity)
    u64 u = 0;
    u64 v = 0;
    int64_t gap = 0;
};

struct thread_anchor_run {
    u64 anchor = 0;
    uint32_t query_start = 0;
    uint32_t query_end = 0;
};

template <class ResolveAnchor>
static void append_thread_bridge_observations(
    const retained_read &read, ResolveAnchor resolve_anchor,
    gtl::vector<thread_bridge_read_observation> &observations) {
    gtl::vector<thread_anchor_run> runs;
    bool previous_was_anchor = false;
    for (size_t k = 0; k < read.nodes.size(); ++k) {
        const u64 anchor = resolve_anchor(read.nodes[k]);
        if (anchor == 0) {
            previous_was_anchor = false;
            continue;
        }
        const uint32_t query_start = read.query_pos[k];
        const uint32_t query_end = read.query_pos[k + 1];
        if (previous_was_anchor && !runs.empty() &&
            runs.back().anchor == anchor) {
            runs.back().query_end = query_end;
        } else {
            runs.push_back({anchor, query_start, query_end});
        }
        previous_was_anchor = true;
    }
    for (size_t i = 1; i < runs.size(); ++i) {
        const auto &source = runs[i - 1];
        const auto &destination = runs[i];
        if (source.anchor == destination.anchor) continue;
        // Fragment intervals are half-open. GL is this anchor-run separation:
        // the destination run's first unrepresented base minus the source
        // run's last represented base. It is not corrected by a guessed
        // supernode-boundary offset; FASTA uses it only after conservative
        // uncertainty gates, so exact boundary-offset correction remains a
        // limitation. Zero-id spans still contribute to this coordinate.
        const int64_t gap = static_cast<int64_t>(destination.query_start) -
                            static_cast<int64_t>(source.query_end);
        observations.push_back(
            {read.qname_hash, source.anchor, destination.anchor, gap});
    }
}

// Collapse repeated occurrences of one directed bridge on one qname to one
// read vote. Revisited anchors can produce several separations; their lower
// empirical median is the single observation retained for that read/pair.
static gtl::vector<thread_bridge_read_observation>
deduplicate_thread_bridge_observations(
    gtl::vector<thread_bridge_read_observation> observations) {
    // Within a (qname,u,v) group all fields but gap are equal, so retaining the
    // lower-median element is the same as rebuilding it with the median gap.
    return layout_graph_emit::collapse_runs_lower_median(
        std::move(observations),
        [](const auto &a, const auto &b) {
            return std::tie(a.qname, a.u, a.v, a.gap) <
                   std::tie(b.qname, b.u, b.v, b.gap);
        },
        [](const auto &a, const auto &b) {
            return a.qname == b.qname && a.u == b.u && a.v == b.v;
        });
}

struct assembler {
    using SYNCMER_HASHER = syncmer::MurMur;
    using SyncmerType = syncmer::GenericSubmerRule<uint32_t, SYNCMER_HASHER>;

    using dqt = moodycamel::BlockingConcurrentQueue<bam1_t *>;
    gtl::vector<string> inputs;
    string output;
    decltype(asm_config::threads) threads;
    int lcp_level;
    double density;
    decltype(asm_config::bloom) bloom;
    uint64_t seed;
    decltype(asm_config::syncmer) syncmer_params;
    SyncmerType syncmer_functor;
    int dbg_lcmer_l;
    int olc_lcmer_l;
    bool parallel_file_reading;
    std::string dump_insert_types_path;
    decltype(asm_config::thread_bridge) thread_bridge;
    // --thread-bridge: compact retained reads. Every raw lcmer (including id 0)
    // carries its emitted-fragment query interval; no read sequence is retained.
    gtl::vector<retained_read> kept_walks;
    // Confirmed read-threaded bridges (u->v DBG unitig reps + support). Filled by
    // run_thread_bridge; consumed at emit to add scaffold edges to the GFA.
    gtl::vector<layout_graph_emit::bridge_edge> confirmed_bridges;

    assembler(const asm_config &cfg)
        : inputs{},
          output{cfg.output},
          threads{cfg.threads},
          lcp_level{cfg.lcp_level},
          bloom{cfg.bloom},
          seed{cfg.seed},
          syncmer_params{cfg.syncmer},
          // The syncmer hash seed is 32-bit; cfg.seed is size_t only so the
          // CLI can accept large values. Deliberate truncation — different
          // seeds may alias mod 2^32, which is acceptable for a hash seed.
          syncmer_functor{syncmer_params, static_cast<unsigned int>(seed)},
          dbg_lcmer_l{cfg.dbg_lcmer_l},
          olc_lcmer_l{cfg.olc_lcmer_l},
          parallel_file_reading{cfg.parallel_file_reading},
          dump_insert_types_path{cfg.dump_insert_types_path},
          thread_bridge{cfg.thread_bridge} {
        inputs.assign(cfg.inputs.begin(), cfg.inputs.end());
        // Overlap detection (fold_dbg_graph → find_lcmer_ranges) queries the
        // DBG, which is keyed by length-`dbg_lcmer_l` lcmers. The litecmer hash
        // is window-length dependent, so querying with a different length misses
        // every key. The two lengths must therefore match for overlaps to be
        // found; warn loudly rather than silently produce an empty graph.
        if (olc_lcmer_l != dbg_lcmer_l) {
            fprintf(stderr,
                    "[config] WARNING: --olc-lcmer-l (%d) != --dbg-lcmer-l (%d). "
                    "Overlap detection queries the DBG built at --dbg-lcmer-l, "
                    "so mismatched lengths lose overlaps and collapse the "
                    "layout. Set them equal unless you are deliberately "
                    "probing this behavior.\n",
                    olc_lcmer_l, dbg_lcmer_l);
        }
    }

    template <class BFT>
    static inline auto insert_lcmer_helper(const gtl::vector<uint64_t> &ll,
                                           const uint32_t &i,
                                           dbg_index_base<BFT> &dbgi) -> bool {
        using index_type = dbg_index_base<BFT>;
        if (ll[i] == 0) {
            return false;
        }
        if (ll[i - 1] == 0) {
            if (ll[i + 1] == 0) {
                return false;
            } else {
                return dbgi.dbg.lazy_emplace_l(
                    ll[i],
                    [&](index_type::dbg_map_t::value_type &val) {
                        val.second.counts.core_count++;
                    },
                    [&](const index_type::dbg_map_t::constructor &ctr) {
                        ctr(ll[i], typename index_type::bidbg_val{1});
                    });
            }
        } else {
            if (ll[i + 1] == 0) {
                return dbgi.dbg.lazy_emplace_l(
                    ll[i],
                    [&](index_type::dbg_map_t::value_type &val) {
                        val.second.emplace_prev(
                            ll[i - 1]);  // TODO verify this CC logic holds
                        // val.second.cc++;
                    },
                    [&](const index_type::dbg_map_t::constructor &ctr) {
                        ctr(ll[i],
                            typename index_type::bidbg_val{0U, ll[i - 1]});
                    });
            } else {
                return dbgi.dbg.lazy_emplace_l(
                    ll[i],
                    [&](index_type::dbg_map_t::value_type &val) {
                        val.second.emplace_prev(ll[i - 1]);
                        val.second.counts.core_count++;
                    },
                    [&](const index_type::dbg_map_t::constructor &ctr) {
                        ctr(ll[i],
                            typename index_type::bidbg_val{1U, ll[i - 1]});
                    });
            }
        }
    }
    // Per-read tally of DBG insertion events (see insert_lcmer_helper_classify).
    // novel/coverage are mutually exclusive per step; branch_open/branch_join
    // are independent flags of the same new-edge event and can co-occur (with
    // each other and with novel), so they are counted separately.
    struct insert_tally {
        uint32_t steps = 0;        // classified insertion positions
        uint32_t novel = 0;        // node newly constructed
        uint32_t coverage = 0;     // existing node, edge already present
        uint32_t branch_open = 0;  // predecessor gained a 2nd+ distinct successor
        uint32_t branch_join = 0;  // current node gained a 2nd+ distinct predecessor
        // The read's coverage-gated node-id walk (0 = lcmer dropped by the
        // coverage>=2 gate). Consecutive non-zero ids are DBG edges; the whole
        // adjacency + node degrees are reconstructable from all reads' walks,
        // so this alone supports the path-level (component-spanning) analysis.
        gtl::vector<uint64_t> walk;
    };

    // Classifying twin of insert_lcmer_helper: performs the IDENTICAL graph
    // mutation (so the built DBG is unchanged) while recording what each
    // insertion did. Only invoked on the --dump-insert-types path. `succ_track`
    // gives the out-edge view needed for branch-open (the DBG stores only prev).
    template <class BFT>
    static inline bool insert_lcmer_helper_classify(
        const gtl::vector<uint64_t> &ll, const uint32_t &i,
        dbg_index_base<BFT> &dbgi, insert_tally &tally) {
        using index_type = dbg_index_base<BFT>;
        const u64 cur = ll[i];
        if (cur == 0) return false;
        const u64 prevn = ll[i - 1];
        const u64 nextn = ll[i + 1];
        const bool has_prev = (prevn != 0);
        const bool has_next = (nextn != 0);
        if (!has_prev && !has_next) return false;  // isolated (matches base)

        tally.steps++;
        bool constructed = false;
        bool edge_added = false;
        size_t prev_size_before = 0;
        dbgi.dbg.lazy_emplace_l(
            cur,
            [&](typename index_type::dbg_map_t::value_type &val) {
                if (has_prev) {
                    prev_size_before = val.second.prev.size();
                    val.second.emplace_prev(prevn);  // no-op if already present
                    edge_added = (val.second.prev.size() != prev_size_before);
                }
                if (has_next) val.second.counts.core_count++;
            },
            [&](const typename index_type::dbg_map_t::constructor &ctr) {
                constructed = true;
                const unsigned cc = has_next ? 1U : 0U;
                if (has_prev) {
                    ctr(cur, typename index_type::bidbg_val{cc, prevn});
                    edge_added = true;  // the constructed node's prev is a new edge
                } else {
                    ctr(cur, typename index_type::bidbg_val{cc});
                }
            });

        if (constructed) {
            tally.novel++;
        } else if (!edge_added) {
            // Existing node, no new in-edge (edge already present, or no prev).
            tally.coverage++;
        }
        if (edge_added && !constructed && prev_size_before >= 1) {
            tally.branch_join++;  // in-degree grew to >= 2
        }
        // Branch-open: did prevn gain a NEW successor, and did it already have one?
        if (edge_added && has_prev) {
            size_t succ_before = 0;
            bool succ_added = false;
            dbgi.succ_track.lazy_emplace_l(
                prevn,
                [&](auto &val) {
                    succ_before = val.second.size();
                    if (std::find(val.second.begin(), val.second.end(), cur) ==
                        val.second.end()) {
                        val.second.push_back(cur);
                        succ_added = true;
                    }
                },
                [&](const auto &ctr) {
                    ctr(prevn, gtl::vector<u64>{cur});
                    succ_added = true;
                });
            if (succ_added && succ_before >= 1) tally.branch_open++;
        }
        return constructed;  // drives sequence attachment (as in the base helper)
    }

    static gtl::vector<uint64_t> find_lcmers(const lcp_read &lr,
                                             const dbg_index &dbgi,
                                             int lcmer_l) {
        LCmer l{lr.cores.begin() + 1, lcmer_l};  // Skip first core;

        gtl::vector<uint64_t> ll(1, 0);
        ll.reserve(lr.cores.size() - lcmer_l + 3);
        dbg_index::lci_type lval = dbgi.get(l);
        ll.push_back(lval.id);
        for (size_t i = static_cast<size_t>(lcmer_l + 1);
             i < lr.cores.size() - 1; ++i) {
            uint32_t in_core = lr.cores[i];
            l.update_lcmer(in_core, lr.cores[i - lcmer_l], lcmer_l);
            lval = dbgi.get(l);
            ll.push_back(lval.id);
        }
        ll.push_back(0);
        return ll;
    }
    static auto find_lcmer_ranges(const lcp_read &lr, const dbg_index &dbgi,
                                  int lcmer_l) {
        // Cores skipped at each node end before overlap-anchoring (see the
        // OLC_SKIP_* macros). Node ends are boundary-noisy, so a small skip
        // avoids false anchors.
        long skip_first = OLC_SKIP_FIRST(lcmer_l);
        long skip_last = OLC_SKIP_LAST(lcmer_l);
        // Need enough cores to form at least one valid lcmer after skipping
        if (lr.cores.size() <
            static_cast<size_t>(lcmer_l + skip_first + skip_last)) {
            return gtl::vector<std::tuple<u64, u64, u64, LCmer, LCmer>>{};
        }

        LCmer l{lr.cores.begin() + skip_first, lcmer_l};  // Skip first core;
        gtl::vector<std::tuple<u64, u64, u64, LCmer, LCmer>> llr;

        dbg_index::lci_type lval = dbgi.get(l);
        llr.emplace_back(lval.id, skip_first, skip_first, l, l);
        // Fix: ensure we don't go out of bounds when accessing i +
        // runtime_lcmer_l
        for (long long i = skip_first + 1;
             i < ((long long)lr.cores.size() - lcmer_l - skip_last); ++i) {
            uint32_t in_core = lr.cores[i + lcmer_l - 1];
            l.update_lcmer(in_core, lr.cores[i - 1], lcmer_l);
            lval = dbgi.get(l);
            if (lval.id == std::get<0>(llr.back())) {  // Extend the previous
                std::get<2>(llr.back()) = i;
                std::get<4>(llr.back()) = l;
            } else if (lval.id != 0) {
                llr.emplace_back(lval.id, i, i, l, l);
            }
        }
        return llr;
    }
    template <class BFT>
    static void process_cores(const lcp_read &lr, dbg_index_base<BFT> &dbgi,
                              int thread_id, uint64_t &lcmer_index,
                              int lcmer_l) {
        LCmer l{lr.cores.begin() + 1, lcmer_l};  // Skip first core;
        gtl::vector<std::tuple<u64, u64, u64>> llr;

        gtl::vector<uint64_t> ll(1, 0);
        ll.reserve(lr.cores.size() - lcmer_l + 3);
        ll.push_back(dbgi.emplace(l, lcmer_index, thread_id));
        for (size_t i = static_cast<size_t>(lcmer_l + 1);
             i < lr.cores.size() - 1; ++i) {
            uint32_t in_core = lr.cores[i];
            l.update_lcmer(in_core, lr.cores[i - lcmer_l], lcmer_l);
            ll.push_back(dbgi.emplace(l, lcmer_index, thread_id));
        }
        ll.push_back(0);

        for (size_t i = 1; i < ll.size() - 1; ++i) {
            if (insert_lcmer_helper(ll, i, dbgi)) {
                // Lets try atomic
                dbgi.dbg.modify_if(
                    ll[i],
                    [&lr, i,
                     lcmer_l](dbg_index_base<BFT>::dbg_map_t::value_type &val) {
                        val.second.seq.emplace_back(
                            lr.data, lr.pos[i + lcmer_l - 1].second,
                            lr.pos[i + lcmer_l].second, val.first);
                    });
            }
        }
    }
    // Identical to process_cores (same emplace/insert/seq mutations) but
    // classifies each insertion into `tally`. Used only on --dump-insert-types.
    template <class BFT>
    static void process_cores_classify(const lcp_read &lr,
                                       dbg_index_base<BFT> &dbgi, int thread_id,
                                       uint64_t &lcmer_index, int lcmer_l,
                                       insert_tally &tally) {
        LCmer l{lr.cores.begin() + 1, lcmer_l};  // Skip first core;
        gtl::vector<uint64_t> ll(1, 0);
        ll.reserve(lr.cores.size() - lcmer_l + 3);
        ll.push_back(dbgi.emplace(l, lcmer_index, thread_id));
        for (size_t i = static_cast<size_t>(lcmer_l + 1);
             i < lr.cores.size() - 1; ++i) {
            uint32_t in_core = lr.cores[i];
            l.update_lcmer(in_core, lr.cores[i - lcmer_l], lcmer_l);
            ll.push_back(dbgi.emplace(l, lcmer_index, thread_id));
        }
        ll.push_back(0);
        tally.walk.assign(ll.begin(), ll.end());

        for (size_t i = 1; i < ll.size() - 1; ++i) {
            if (insert_lcmer_helper_classify(ll, i, dbgi, tally)) {
                dbgi.dbg.modify_if(
                    ll[i],
                    [&lr, i,
                     lcmer_l](dbg_index_base<BFT>::dbg_map_t::value_type &val) {
                        val.second.seq.emplace_back(
                            lr.data, lr.pos[i + lcmer_l - 1].second,
                            lr.pos[i + lcmer_l].second, val.first);
                    });
            }
        }
    }
    // Is `node` currently a branch (in-deg>=2 from prev, or out-deg>=2 from the
    // classify-path successor tracker)? Reflects reads processed so far.
    template <class BFT>
    static bool tb_is_branch(dbg_index_base<BFT> &dbgi, u64 node) {
        bool is_branch = false;
        // The vectors are mutated by other construction workers. Read their
        // sizes only inside phmap's shard-locked callback; a locked find alone
        // would not keep the returned iterator/value protected.
        dbgi.dbg.modify_if(node, [&](auto &entry) {
            is_branch = entry.second.prev.size() >= 2;
        });
        if (is_branch) return true;
        dbgi.succ_track.modify_if(node, [&](auto &entry) {
            is_branch = entry.second.size() >= 2;
        });
        return is_branch;
    }

    // --thread-bridge walk capture: identical DBG mutation as process_cores
    // (uses the classify helper only to maintain the successor tracker + branch
    // info), then STREAMING retention: keep the walk only if it crosses a branch
    // (spans a repeat), and cap retained reads per bridge via the counting bloom.
    template <class BFT>
    static void process_cores_thread(const lcp_read &lr, dbg_index_base<BFT> &dbgi,
                                     int thread_id, uint64_t &lcmer_index,
                                     int lcmer_l,
                                     const char *qname,
                                     gtl::vector<retained_read> &walk_sink,
                                     CountingBloom *bloom, int cap) {
        LCmer l{lr.cores.begin() + 1, lcmer_l};
        gtl::vector<uint64_t> ll(1, 0);
        ll.reserve(lr.cores.size() - lcmer_l + 3);
        ll.push_back(dbgi.emplace(l, lcmer_index, thread_id));
        for (size_t i = static_cast<size_t>(lcmer_l + 1);
             i < lr.cores.size() - 1; ++i) {
            l.update_lcmer(lr.cores[i], lr.cores[i - lcmer_l], lcmer_l);
            ll.push_back(dbgi.emplace(l, lcmer_index, thread_id));
        }
        ll.push_back(0);
        insert_tally tally;  // drives succ_track maintenance; counts unused here
        for (size_t i = 1; i < ll.size() - 1; ++i) {
            if (insert_lcmer_helper_classify(ll, i, dbgi, tally)) {
                dbgi.dbg.modify_if(
                    ll[i],
                    [&lr, i,
                     lcmer_l](dbg_index_base<BFT>::dbg_map_t::value_type &val) {
                        val.second.seq.emplace_back(
                            lr.data, lr.pos[i + lcmer_l - 1].second,
                            lr.pos[i + lcmer_l].second, val.first);
                    });
            }
        }

        // Streaming keep decision: does the walk cross a branch (a repeat)?
        long first_br = -1, last_br = -1;
        for (size_t i = 1; i + 1 < ll.size(); ++i) {
            if (ll[i] && tb_is_branch(dbgi, ll[i])) {
                if (first_br < 0) first_br = static_cast<long>(i);
                last_br = static_cast<long>(i);
            }
        }
        if (first_br < 0) return;  // no branch crossed -> not a bridge candidate

        // Dedup key = the (entry-branch, exit-branch) node pair the walk crosses.
        // These are DBG node ids shared by reads spanning the same tangle, so the
        // cap binds per bridge -> retained reads ~ O(distinct bridges x cap),
        // independent of sequencing depth (unlike flank anchors, which shift with
        // each tiling read).
        if (bloom) {
            const u64 x = ll[first_br], y = ll[last_br];
            const u64 a = std::min(x, y), b = std::max(x, y);
            const u64 key = a * 0x9E3779B97F4A7C15ULL ^ (b + 0x632BE59BD9B4E019ULL);
            if (bloom->bump(key) >= static_cast<uint8_t>(cap)) return;  // over cap
        }
        retained_read retained;
        retained.qname_hash = hash_qname(qname);
        retained.nodes.reserve(ll.size() - 2);
        retained.query_pos.reserve(ll.size() - 1);
        for (size_t i = 1; i + 1 < ll.size(); ++i) {
            // The sequence represented by ll[i] is the post-lcmer fragment
            // [pos[i+l-1].second, pos[i+l].second), matching the fragment
            // inserted above. Retain it even when ll[i] is zero so missing DBG
            // nodes still contribute their exact base span to a bridge gap.
            // Fragments are contiguous (each end == the next start), so we store
            // only the start here and append the final end after the loop.
            retained.nodes.push_back(ll[i]);
            retained.query_pos.push_back(
                lr.pos[i + static_cast<size_t>(lcmer_l) - 1].second);
        }
        // Trailing boundary = query_end of the last retained fragment.
        retained.query_pos.push_back(
            lr.pos[(ll.size() - 2) + static_cast<size_t>(lcmer_l)].second);
        walk_sink.push_back(std::move(retained));
    }

    // Chase a raw (pre-cleaning) node id to its post-cleaning unitig
    // representative (mirrors dbg_index::get but from an id). 0 = dropped.
    static u64 rep_of(const dbg_index &dbg, u64 id) {
        for (int hop = 0; hop < 128 && id != 0; ++hop) {
            auto it = dbg.dbg.find(id);
            if (it == dbg.dbg.end()) return 0;
            if (!it->second.redirected()) return id;
            const u64 nxt = it->second.redirect_to;
            if (nxt == 0 || nxt == id) return 0;
            id = nxt;
        }
        return 0;
    }

    // Post-cleaning read-threading bridge pass: contract each retained read to
    // anchor runs, vote directed anchor adjacencies once per qname, and confirm
    // the unambiguous, well-supported ones.
    std::expected<void, std::string> run_thread_bridge(
        dbg_index &dbg, const std::string &out_path) {
        using layout_graph_emit::bridge_edge;
        using layout_graph_emit::bridge_realization;
        using layout_graph_emit::bridge_realization_label;
        using layout_graph_emit::classify_bridge_realization;
        using layout_graph_emit::signed_gap_summary;
        using layout_graph_emit::summarize_signed_gaps;
        const int min_anchor = thread_bridge.min_anchor;
        // Single-copy depth baseline: median avg_core_count of long, non-branching
        // unitigs. Bridge anchors must be UNIQUE (in/out deg <=1) AND single-copy
        // (avg_core_count <= 1.5x baseline) so bridges cannot route through
        // collapsed-repeat unitigs (the false-join source).
        gtl::vector<double> depths;
        for (auto &kv : dbg.dbg) {
            auto &n = kv.second;
            if (n.redirected() || n.seq.empty()) continue;
            if (n.seq.size() >= 100 && n.prev.size() <= 1 &&
                n.forward_edge_count <= 1) {
                depths.push_back(n.avg_core_count());
            }
        }
        std::sort(depths.begin(), depths.end());
        const double baseline =
            depths.empty() ? 0.0 : depths[depths.size() / 2];
        const double cn_cap = baseline > 0 ? baseline * 1.5 : 1e18;
        auto substantial = [&](u64 rep) -> bool {
            auto it = dbg.dbg.find(rep);
            if (it == dbg.dbg.end() || it->second.redirected() ||
                it->second.seq.size() < static_cast<size_t>(min_anchor))
                return false;
            const auto &n = it->second;
            return n.prev.size() <= 1 && n.forward_edge_count <= 1 &&
                   n.avg_core_count() <= cn_cap;  // unique & single-copy
        };
        gtl::vector<thread_bridge_read_observation> observations;
        for (const auto &read : kept_walks) {
            append_thread_bridge_observations(
                read,
                [&](u64 raw_node) {
                    const u64 rep = rep_of(dbg, raw_node);
                    return rep != 0 && substantial(rep) ? rep : u64{0};
                },
                observations);
        }
        auto distinct =
            deduplicate_thread_bridge_observations(std::move(observations));

        using pair_key = std::pair<u64, u64>;
        gtl::btree_map<pair_key, gtl::vector<thread_bridge_read_observation>>
            by_pair;
        gtl::btree_map<u64, gtl::btree_map<u64, uint32_t>> outgoing, incoming;
        gtl::btree_map<u64, uint32_t> out_tot, in_tot;
        for (const auto &observation : distinct) {
            by_pair[{observation.u, observation.v}].push_back(observation);
            ++outgoing[observation.u][observation.v];
            ++incoming[observation.v][observation.u];
            ++out_tot[observation.u];
            ++in_tot[observation.v];
        }

        const std::string reads_path = out_path + ".reads.tsv";
        const std::string diagnostics_path = out_path + ".diagnostics.tsv";
        FILE *tf = fopen(out_path.c_str(), "w");
        if (!tf)
            return std::unexpected("Cannot open thread-bridge output: " +
                                   out_path);
        FILE *rf = fopen(reads_path.c_str(), "w");
        if (!rf) {
            fclose(tf);
            return std::unexpected("Cannot open thread-bridge reads sidecar: " +
                                   reads_path);
        }
        FILE *df = fopen(diagnostics_path.c_str(), "w");
        if (!df) {
            fclose(rf);
            fclose(tf);
            return std::unexpected(
                "Cannot open thread-bridge diagnostics sidecar: " +
                diagnostics_path);
        }

        auto runner_up = [](const gtl::btree_map<u64, uint32_t> &votes,
                            u64 winner) {
            uint32_t second = 0;
            for (const auto &[candidate, count] : votes)
                if (candidate != winner) second = std::max(second, count);
            return second;
        };
        struct anchor_diagnostics {
            size_t length_cores = 0;
            size_t length_bases = 0;
            double depth = 0;
            size_t in_degree = 0;
            size_t out_degree = 0;
            double copy_number = 0;
        };
        auto anchor_info = [&](u64 id) {
            anchor_diagnostics info;
            auto it = dbg.dbg.find(id);
            if (it == dbg.dbg.end() || it->second.redirected() ||
                it->second.seq.empty())
                return info;
            info.length_cores = it->second.seq.size();
            for (const auto &sequence : it->second.seq)
                info.length_bases += sequence.len;
            info.depth = it->second.avg_core_count();
            info.in_degree = it->second.prev.size();
            info.out_degree = it->second.forward_edge_count;
            info.copy_number = baseline > 0 ? info.depth / baseline : 0.0;
            return info;
        };

        fprintf(tf,
                "u\tv\tsupport\tgap\tcount\tp10\tp90\tmad\trealization\n");
        fprintf(rf, "u\tv\tqname_hash\tgap\n");
        fprintf(df,
                "u\tv\tu_length_cores\tu_length_bases\tu_depth\tu_in_degree\t"
                "u_out_degree\tu_cn\tv_length_cores\tv_length_bases\tv_depth\t"
                "v_in_degree\tv_out_degree\tv_cn\t"
                "support\tu_out_total\tu_out_runner_up\tv_in_total\t"
                "v_in_runner_up\tgap_count\tgap_median\tgap_p10\tgap_p90\t"
                "gap_mad\tconfirmed\trealization\n");

        confirmed_bridges.clear();
        for (auto &[pair, pair_observations] : by_pair) {
            const auto [a, b] = pair;
            gtl::vector<int64_t> gaps;
            gaps.reserve(pair_observations.size());
            for (const auto &observation : pair_observations)
                gaps.push_back(observation.gap);
            const signed_gap_summary gap_summary = summarize_signed_gaps(gaps);
            const uint32_t support = gap_summary.count;
            const uint32_t out_second = runner_up(outgoing[a], b);
            const uint32_t in_second = runner_up(incoming[b], a);
            const bool confirmed =
                support >= static_cast<uint32_t>(thread_bridge.min_support) &&
                support == out_tot[a] && support == in_tot[b];
            const bridge_realization realization =
                classify_bridge_realization(gap_summary);
            if (confirmed) {
                fprintf(tf,
                        "%" PRIu64 "\t%" PRIu64 "\t%u\t%" PRId64
                        "\t%u\t%" PRId64 "\t%" PRId64 "\t%" PRId64
                        "\t%s\n",
                        a, b, support, gap_summary.median, gap_summary.count,
                        gap_summary.p10, gap_summary.p90, gap_summary.mad,
                        bridge_realization_label(realization));
                gtl::vector<layout_graph_emit::bridge_read_gap_observation>
                    read_gaps;
                read_gaps.reserve(pair_observations.size());
                for (const auto &observation : pair_observations)
                    read_gaps.push_back(
                        {observation.qname, observation.gap});
                confirmed_bridges.push_back(
                    {a, b, support, gap_summary.median, gap_summary,
                     realization, std::move(read_gaps)});
                for (const auto &observation : pair_observations)
                    fprintf(rf,
                            "%" PRIu64 "\t%" PRIu64 "\t%016" PRIx64 "\t%" PRId64
                            "\n",
                            a, b, observation.qname, observation.gap);
            }
            const anchor_diagnostics ai = anchor_info(a);
            const anchor_diagnostics bi = anchor_info(b);
            const char *diagnostic_status =
                confirmed
                    ? bridge_realization_label(realization)
                    : (support < static_cast<uint32_t>(thread_bridge.min_support)
                           ? "not_confirmed_support"
                           : "not_confirmed_dominance");
            fprintf(df,
                    "%" PRIu64 "\t%" PRIu64
                    "\t%zu\t%zu\t%.6f\t%zu\t%zu\t%.6f"
                    "\t%zu\t%zu\t%.6f\t%zu\t%zu\t%.6f"
                    "\t%u\t%u\t%u\t%u\t%u"
                    "\t%u\t%" PRId64 "\t%" PRId64 "\t%" PRId64
                    "\t%" PRId64 "\t%s\t%s\n",
                    a, b, ai.length_cores, ai.length_bases, ai.depth,
                    ai.in_degree, ai.out_degree, ai.copy_number,
                    bi.length_cores, bi.length_bases, bi.depth, bi.in_degree,
                    bi.out_degree, bi.copy_number, support, out_tot[a],
                    out_second, in_tot[b], in_second, gap_summary.count,
                    gap_summary.median, gap_summary.p10, gap_summary.p90,
                    gap_summary.mad, confirmed ? "yes" : "no",
                    diagnostic_status);
        }
        fclose(df);
        fclose(rf);
        fclose(tf);

        capture_representative_gap_paths(dbg);

        fprintf(stderr,
                "thread-bridge: %zu retained reads -> %zu distinct read votes "
                "-> %zu confirmed bridges (support>=%d, full dominance) -> %s\n",
                kept_walks.size(), distinct.size(), confirmed_bridges.size(),
                thread_bridge.min_support, out_path.c_str());
        return {};
    }

    // Tier-2 capture: for each confirmed bridge, pick the supporting walk whose
    // tangle gap is closest to the confirmed median gap (ties broken
    // deterministically) and record the cleaned-unitig reps it crossed, with the
    // raw entry/exit node of each. Emit later reconstructs the real gap sequence
    // from these; walks are only re-scanned for confirmed bridges, so no
    // per-observation path storage is needed.
    void capture_representative_gap_paths(dbg_index &dbg) {
        using layout_graph_emit::tangle_segment;
        if (confirmed_bridges.empty()) return;
        const int min_anchor = thread_bridge.min_anchor;
        gtl::vector<double> depths;
        for (auto &kv : dbg.dbg) {
            auto &n = kv.second;
            if (n.redirected() || n.seq.empty()) continue;
            if (n.seq.size() >= 100 && n.prev.size() <= 1 &&
                n.forward_edge_count <= 1)
                depths.push_back(n.avg_core_count());
        }
        std::sort(depths.begin(), depths.end());
        const double cn_cap =
            depths.empty() ? 1e18 : depths[depths.size() / 2] * 1.5;
        auto substantial = [&](u64 rep) -> bool {
            auto it = dbg.dbg.find(rep);
            if (it == dbg.dbg.end() || it->second.redirected() ||
                it->second.seq.size() < static_cast<size_t>(min_anchor))
                return false;
            const auto &n = it->second;
            return n.prev.size() <= 1 && n.forward_edge_count <= 1 &&
                   n.avg_core_count() <= cn_cap;
        };
        auto resolve = [&](u64 raw) -> u64 {
            const u64 rep = rep_of(dbg, raw);
            return rep != 0 && substantial(rep) ? rep : u64{0};
        };

        using pair_key = std::pair<u64, u64>;
        gtl::btree_map<pair_key, size_t> confirmed_index;
        for (size_t i = 0; i < confirmed_bridges.size(); ++i)
            confirmed_index[{confirmed_bridges[i].u, confirmed_bridges[i].v}] = i;
        struct sel_key {
            int64_t dist;
            u64 entry;
            u64 qname;
        };
        gtl::btree_map<pair_key, sel_key> best;

        struct anchor_run {
            u64 anchor;
            size_t first_idx;
            size_t last_idx;
            uint32_t query_start;
            uint32_t query_end;
        };
        gtl::vector<anchor_run> runs;
        for (const auto &read : kept_walks) {
            runs.clear();
            bool previous_was_anchor = false;
            for (size_t k = 0; k < read.nodes.size(); ++k) {
                const u64 anchor = resolve(read.nodes[k]);
                if (anchor == 0) {
                    previous_was_anchor = false;
                    continue;
                }
                const uint32_t query_start = read.query_pos[k];
                const uint32_t query_end = read.query_pos[k + 1];
                if (previous_was_anchor && !runs.empty() &&
                    runs.back().anchor == anchor) {
                    runs.back().last_idx = k;
                    runs.back().query_end = query_end;
                } else {
                    runs.push_back({anchor, k, k, query_start, query_end});
                }
                previous_was_anchor = true;
            }
            for (size_t i = 1; i < runs.size(); ++i) {
                const u64 a = runs[i - 1].anchor, b = runs[i].anchor;
                if (a == b) continue;
                auto ci = confirmed_index.find({a, b});
                if (ci == confirmed_index.end()) continue;
                const int64_t gap = static_cast<int64_t>(runs[i].query_start) -
                                    static_cast<int64_t>(runs[i - 1].query_end);
                const int64_t median = confirmed_bridges[ci->second].gap;
                const int64_t dist = gap >= median ? gap - median : median - gap;
                // Reps strictly between the two anchor runs, consecutive
                // same-rep steps collapsed into one segment.
                gtl::vector<tangle_segment> segments;
                for (size_t k = runs[i - 1].last_idx + 1;
                     k < runs[i].first_idx; ++k) {
                    const u64 rep = rep_of(dbg, read.nodes[k]);
                    if (rep == 0) continue;
                    if (!segments.empty() && segments.back().rep == rep)
                        segments.back().exit_raw = read.nodes[k];
                    else
                        segments.push_back({rep, read.nodes[k], read.nodes[k]});
                }
                const u64 first_entry =
                    segments.empty() ? 0 : segments.front().entry_raw;
                const sel_key key{dist, first_entry, read.qname_hash};
                auto bit = best.find({a, b});
                if (bit == best.end() ||
                    std::tie(key.dist, key.entry, key.qname) <
                        std::tie(bit->second.dist, bit->second.entry,
                                 bit->second.qname)) {
                    best[{a, b}] = key;
                    confirmed_bridges[ci->second].gap_path = std::move(segments);
                }
            }
        }

        size_t with_path = 0, empty_path = 0, single_rep = 0, seg_total = 0,
               seg_max = 0;
        for (const auto &b : confirmed_bridges) {
            if (b.gap_path.empty()) {
                ++empty_path;
                continue;
            }
            ++with_path;
            if (b.gap_path.size() == 1) ++single_rep;
            seg_total += b.gap_path.size();
            seg_max = std::max(seg_max, b.gap_path.size());
        }
        fprintf(stderr,
                "thread-bridge: gap-path captured for %zu/%zu confirmed "
                "(%zu empty, %zu single-rep; avg %.2f reps, max %zu)\n",
                with_path, confirmed_bridges.size(), empty_path, single_rep,
                with_path ? static_cast<double>(seg_total) / with_path : 0.0,
                seg_max);
    }

    static void process_cores_for_basespace(
        const lcp_read &lr, dbg_index &dbgi, [[maybe_unused]] int thread_id,
        [[maybe_unused]] uint64_t &lcmer_index, int lcmer_l) {
        // Keep the DBG node sequences in BASE space (actual substrings), but do
        // NOT convert overlap/layout coordinates to base space during
        // layouting/merging.
        //
        // Layouting will be done in core-index space. Only at the final
        // emission step will we scale core-space layout positions into base
        // offsets.
        LCmer l{lr.cores.begin() + 1, lcmer_l};  // Skip first core;

        uint64_t val = dbgi.get(l).id;

        for (size_t i = static_cast<size_t>(lcmer_l + 1);
             i < lr.cores.size() - 1; ++i) {
            uint64_t from = lr.pos[i - 1].second;
            uint64_t to = lr.pos[i].second;

            dbgi.dbg.modify_if(val, [&](dbg_index::dbg_map_t::value_type &val) {
                if (val.second.seq.size() == 0) {
                    if (to - from > 0) {
                        val.second.seq.emplace_back(lr.data, from, to,
                                                    val.first);
                        // assert((to & 3) + 4 < lr.data_length);
                    }
                }
            });
            uint32_t in_core = lr.cores[i];
            l.update_lcmer(in_core, lr.cores[i - lcmer_l], lcmer_l);
            val = dbgi.get(l).id;
        }
    }

    static std::pair<int, int> get_thread_dist(int threads) {
        if (threads == 0) {
            return {0, 0};
        }
        if (threads < 16) {
            int sam_threads = threads / 2;
            return {sam_threads, threads - sam_threads};
        }
        return {8, threads - 8};
    }

    static auto process_data(
        int thread_id, int lcp_level, dqt &data_queue, dqt &recycle_queue,
        dbg_index &dbgi, std::atomic<bool> &done, int &processed_read_count,
        int &all_reads, auto Process_Func, SyncmerType syncmer_functor,
        int lcmer_l,
        gtl::vector<std::pair<std::string, insert_tally>> *insert_sink = nullptr,
        gtl::vector<retained_read> *walk_sink = nullptr,
        CountingBloom *tb_bloom = nullptr, int tb_cap = 0) {
        uint64_t lcmer_index = 0;
        bam1_t *data_p = NULL;
        lps str = {0, 0, 0, 0};
        charsi chs;
        pos_scratch pos_buf;
        lcp_read lr;
        while (true) {
            if (!data_queue.wait_dequeue_timed(data_p, WAIT_QUEUE_TIMEOUT)) {
                // if (!data_queue.try_dequeue(data_p)) {
                if (done.load()) {
                    break;
                }
                // std::this_thread::sleep_for(std::chrono::nanoseconds(100));
                continue;
            }
            all_reads++;
            lr.clear();
            parse_cores_from_bam_record(bam_get_seq(data_p),
                                        data_p->core.l_qseq, lcp_level, &str,
                                        syncmer_functor, chs, lr, pos_buf);
            // Reads with too few cores are skipped, but the record is always
            // recycled afterwards.
            if (lr.cores.size() >= static_cast<size_t>(lcmer_l + 2)) {
                ++processed_read_count;
                if (walk_sink) {
                    process_cores_thread<DefaultBF>(lr, dbgi, thread_id,
                                                    lcmer_index, lcmer_l,
                                                    bam_get_qname(data_p),
                                                    *walk_sink, tb_bloom, tb_cap);
                } else if (insert_sink) {
                    insert_tally tally;
                    process_cores_classify<DefaultBF>(lr, dbgi, thread_id,
                                                      lcmer_index, lcmer_l,
                                                      tally);
                    insert_sink->emplace_back(bam_get_qname(data_p), tally);
                } else {
                    Process_Func(lr, dbgi, thread_id, lcmer_index, lcmer_l);
                }
            }
            if (!recycle_queue.try_enqueue(data_p)) {
                bam_destroy1(data_p);
                data_p = NULL;
            }
        }
        free_lps(&str);
    }

    [[nodiscard]] auto generate_data(const std::string &input, dqt &data_queue,
                                     dqt &recycle_queue, int threads) {
        samFile *fp = sam_open(input.c_str(), "r");

        if (!fp) {
            std::cerr << "Error: Could not open " << input << std::endl;
            return false;
        }
        hts_set_threads(fp, threads);

        bam_hdr_t *hdr = sam_hdr_read(fp);
        bam1_t *rec = nullptr;
        while (true) {
            // Try to get a recycled record object
            if (rec == nullptr && !recycle_queue.try_dequeue(rec)) {
                // If the recycle queue is empty, allocate a new one as a
                // fallback
                rec = bam_init1();
            }

            if (sam_read1(fp, hdr, rec) >= 0) {
                while (!data_queue.try_enqueue(rec)) {
                    std::this_thread::sleep_for(std::chrono::nanoseconds(10));
                }
                rec = nullptr;
            } else {
                // End of file or error. Destroy the unused record and exit.
                bam_destroy1(rec);
                break;
            }
        }
        bam_hdr_destroy(hdr);
        sam_close(fp);
        // if (pool.pool) {
        //     hts_tpool_destroy(pool.pool);
        // }
        return true;
    }

    // Populates the caller-owned dbgi. The dbg_index lives in the caller's
    // storage (see run_assembly) so it is not heap-allocated here.
    std::expected<void, std::string> construct_dbg(dbg_index &dbgi) {
        dqt data_queue{static_cast<size_t>(threads.work) * 1024, 1, 1};
        dqt recycle_queue{static_cast<size_t>(threads.work) * 1024, 1, 1};
        std::atomic<bool> done{false};

        gtl::vector<std::thread> workers;

        gtl::vector<int> valid_read_counts(threads.work, 0);
        gtl::vector<int> all_read_counts(threads.work, 0);

        // --dump-insert-types: one sink per worker (avoids cross-thread
        // contention); merged and written after join. nullptr => normal path.
        const bool classify_inserts = !dump_insert_types_path.empty();
        std::vector<gtl::vector<std::pair<std::string, insert_tally>>>
            insert_sinks(classify_inserts ? threads.work : 0);
        // --thread-bridge: one walk sink per worker, merged into kept_walks; a
        // shared counting bloom caps retained reads per bridge.
        gtl::vector<gtl::vector<retained_read>>
            walk_sinks(thread_bridge.enable ? threads.work : 0);
        std::unique_ptr<CountingBloom> tb_bloom;
        if (thread_bridge.enable) tb_bloom = std::make_unique<CountingBloom>(22);

        for (int tid = 1; tid <= threads.work; ++tid) {
            workers.emplace_back([&, tid]() {
                process_data(tid, lcp_level, data_queue, recycle_queue, dbgi,
                             done, valid_read_counts[tid - 1],
                             all_read_counts[tid - 1], process_cores<DefaultBF>,
                             syncmer_functor, dbg_lcmer_l,
                             classify_inserts ? &insert_sinks[tid - 1] : nullptr,
                             thread_bridge.enable ? &walk_sinks[tid - 1] : nullptr,
                             tb_bloom.get(), thread_bridge.cap);
            });
        }

        std::atomic<bool> read_failed{false};
        if (parallel_file_reading) {
            gtl::vector<std::thread> readers;
            readers.reserve(inputs.size());
            for (const auto &input : inputs) {
                readers.emplace_back([&, input]() {
                    if (!generate_data(input, std::ref(data_queue),
                                       std::ref(recycle_queue), threads.io)) {
                        fprintf(stderr, "Failed to read from input %s\n", input.c_str());
                        read_failed.store(true);
                    }
                });
            }
            for (auto &reader : readers) {
                reader.join();
            }
        } else {
            for (const auto &input : inputs) {
                if (!generate_data(input, std::ref(data_queue),
                                   std::ref(recycle_queue), threads.io)) {
                    fprintf(stderr, "Failed to read from input %s\n", input.c_str());
                    read_failed.store(true);
                    break;
                }
            }
        }
        if (read_failed.load()) {
            done.store(true);
            for (auto &worker : workers) {
                worker.join();
            }
            return std::unexpected(
                "Error in reading data for DBG construction.");
        }

        done.store(true);
        for (auto &worker : workers) {
            worker.join();
        }

        u64 valid_sum = 0;
        u64 all_sum = 0;
        for (size_t i = 0; i < valid_read_counts.size(); ++i) {
            valid_sum += valid_read_counts[i];
            all_sum += all_read_counts[i];
        }
        fprintf(stderr, "%" PRIu64 " processed\n%" PRIu64 " total reads\n",
                valid_sum, all_sum);

        if (thread_bridge.enable) {
            size_t total = 0;
            for (const auto &s : walk_sinks) total += s.size();
            kept_walks.reserve(total);
            for (auto &s : walk_sinks)
                for (auto &w : s) kept_walks.emplace_back(std::move(w));
            fprintf(stderr, "thread-bridge: retained %zu read walks\n",
                    kept_walks.size());
        }

        if (classify_inserts) {
            FILE *ito = fopen(dump_insert_types_path.c_str(), "w");
            if (!ito) {
                return std::unexpected("Cannot open --dump-insert-types file: " +
                                       dump_insert_types_path);
            }
            fprintf(ito,
                    "read_name\tsteps\tnovel\tcoverage\tbranch_open\tbranch_join\n");
            u64 rows = 0;
            for (const auto &sink : insert_sinks) {
                for (const auto &[name, t] : sink) {
                    fprintf(ito, "%s\t%u\t%u\t%u\t%u\t%u\n", name.c_str(),
                            t.steps, t.novel, t.coverage, t.branch_open,
                            t.branch_join);
                    ++rows;
                }
            }
            fclose(ito);
            fprintf(stderr, "insert-type classification written for %" PRIu64
                            " reads to %s\n",
                    rows, dump_insert_types_path.c_str());

            // Sidecar: per-read node-id walk (0 = coverage-gate drop). Supports
            // the path-level component-spanning analysis in Python.
            const std::string walks_path = dump_insert_types_path + ".walks";
            FILE *wto = fopen(walks_path.c_str(), "w");
            if (!wto) {
                return std::unexpected("Cannot open read-walks file: " +
                                       walks_path);
            }
            for (const auto &sink : insert_sinks) {
                for (const auto &[name, t] : sink) {
                    fputs(name.c_str(), wto);
                    for (const u64 id : t.walk) fprintf(wto, "\t%" PRIu64, id);
                    fputc('\n', wto);
                }
            }
            fclose(wto);
            fprintf(stderr, "read walks written to %s\n", walks_path.c_str());
        }
        return {};
    }

    auto convert2basespace(dbg_index &dbgi) -> int {
        dqt data_queue{static_cast<size_t>(threads.work) * 1024, 1, 1};
        dqt recycle_queue{static_cast<size_t>(threads.work) * 1024, 1, 1};
        std::atomic<bool> done{false};

        gtl::vector<std::thread> workers;
        gtl::vector<int> valid_read_counts(threads.work, 0);
        gtl::vector<int> all_read_counts(threads.work, 0);

        for (int tid = 1; tid <= threads.work; ++tid) {
            workers.emplace_back([&, tid]() {
                process_data(
                    tid, lcp_level, data_queue, recycle_queue, dbgi, done,
                    valid_read_counts[tid - 1], all_read_counts[tid - 1],
                    process_cores_for_basespace, syncmer_functor, dbg_lcmer_l);
            });
        }
        std::atomic<bool> read_failed{false};
        if (parallel_file_reading) {
            gtl::vector<std::thread> readers;
            readers.reserve(inputs.size());
            for (const auto &input : inputs) {
                readers.emplace_back([&, input]() {
                    if (!generate_data(input, std::ref(data_queue),
                                       std::ref(recycle_queue), threads.io)) {
                        read_failed.store(true);
                    }
                });
            }
            for (auto &reader : readers) {
                reader.join();
            }
        } else {
            for (const auto &input : inputs) {
                if (!generate_data(input, std::ref(data_queue),
                                   std::ref(recycle_queue), threads.io)) {
                    read_failed.store(true);
                    break;
                }
            }
        }
        if (read_failed.load()) {
            done.store(true);
            for (auto &worker : workers) {
                worker.join();
            }
            fprintf(stderr, "Error in reading data for basespace conversion\n");
            return -1;
        }

        done.store(true);
        for (auto &worker : workers) {
            worker.join();
        }

        int valid_sum = 0;
        int all_sum = 0;
        for (size_t i = 0; i < valid_read_counts.size(); ++i) {
            valid_sum += valid_read_counts[i];
            all_sum += all_read_counts[i];
        }

        fprintf(stderr, "%d processed\n%d total reads\n", valid_sum, all_sum);
        return 0;
    }

    size_t count_seqv_len(const gtl::vector<SequenceContainerType> &seqv) {
        size_t sum = 0;
        for (const SequenceContainerType &s : seqv) {
            sum += s.len;
        }
        return sum;
    }
    // One documented coordinate space for overlaps: inclusive core indices on a
    // node's (reverse-complemented) core sequence. `len` is the node's core-span
    // length in *steps* (cores - 1), carried so downstream never re-derives it.
    struct CoreSpan {
        u64 from;
        u64 to;
        u64 len;
    };

    // Phase-A record: on some key node, an overlap with `other` covers `self`
    // (coordinates on the key node itself). Only the key node's own span is
    // known at this stage; the partner's span is joined in phase B.
    struct SelfMatch {
        u64 other;
        CoreSpan self;
    };

    // Phase-B record (a reciprocal overlap), keyed by the query node: `target`
    // is the partner node id; `query`/`target_span` are each node's core span in
    // the shared overlap. No reinterpreted fields — names match meaning.
    struct Overlap {
        u64 target;
        CoreSpan query;
        CoreSpan target_span;
    };

    static auto build_lr_cores(const gtl::vector<SequenceContainerType> &seqs,
                               const SyncmerType &syncmer_functor,
                               int lcp_level, lps &str) -> lcp_read {
        lcp_read lr;
        lr.clear();

        init_lps_crumbled(&str, seqs, syncmer_functor);
        lps_deepen(&str, lcp_level);

        for (int i = 0; i < str.size; ++i) {
            lr.cores.emplace_back(str.cores[i].label);
            lr.pos.emplace_back(str.cores[i].start, str.cores[i].end);
        }
        return lr;
    }

    static auto build_lr_cores_with_seqs(
        const gtl::vector<SequenceContainerType> &seqs,
        const SyncmerType &syncmer_functor, int lcp_level, lps &str)
        -> std::pair<lcp_read, gtl::vector<SequenceContainerType>> {
        lcp_read lr;
        lr.clear();
        gtl::vector<SequenceContainerType> core_seqs;

        init_lps_crumbled(&str, seqs, syncmer_functor);
        lps_deepen(&str, lcp_level);

        if (str.size == 0 || seqs.empty()) {
            return {std::move(lr), std::move(core_seqs)};
        }

        lr.cores.reserve(str.size);
        lr.pos.reserve(str.size);
        core_seqs.reserve(str.size);

        // Co-iterate: track current fragment and its cumulative position
        size_t frag_idx = 0;
        size_t frag_cumulative = 0;  // start position of current fragment

        size_t prev_end = 0;  // End position of previous core

        for (int i = 0; i < str.size; ++i) {
            size_t cur_start = str.cores[i].start;
            size_t cur_end = str.cores[i].end;
            lr.cores.emplace_back(str.cores[i].label);
            lr.pos.emplace_back(cur_start, cur_end);

            // Extract non-overlapping portion: [prev_end, cur_end)
            // For the first core, this is [0, cur_end)
            size_t from = prev_end;
            size_t to = cur_end;
            size_t sub_len = to - from;

            core_seqs.emplace_back(sub_len, 0);
            auto &result = core_seqs.back();
            size_t dest_pos = 0;

            // Advance frag_idx until we reach the fragment containing 'from'
            while (frag_idx < seqs.size() &&
                   frag_cumulative + seqs[frag_idx].len <= from) {
                frag_cumulative += seqs[frag_idx].len;
                ++frag_idx;
            }

            // Extract bases from fragments, starting from current position
            size_t cur_frag = frag_idx;
            size_t cur_cumulative = frag_cumulative;

            while (cur_frag < seqs.size() && dest_pos < sub_len) {
                const auto &frag = seqs[cur_frag];
                size_t frag_start = cur_cumulative;
                size_t frag_end = cur_cumulative + frag.len;

                // Compute overlap of [frag_start, frag_end) with [from, to)
                size_t overlap_start = std::max(frag_start, from);
                size_t overlap_end = std::min(frag_end, to);

                if (overlap_start < overlap_end) {
                    size_t local_from = overlap_start - frag_start;
                    size_t copy_len = overlap_end - overlap_start;

                    // Copy bases from frag to result
                    for (size_t b = 0; b < copy_len; ++b) {
                        size_t src_bit =
                            (frag.skip_first_n + local_from + b) * 2;
                        uint64_t code =
                            (frag.data[src_bit >> 6] >> (src_bit & 63)) & 0x3;
                        size_t dst_bit = dest_pos * 2;
                        result.data[dst_bit >> 6] |= code << (dst_bit & 63);
                        ++dest_pos;
                    }
                    result.origin_node = frag.origin_node;
                }
                cur_cumulative += frag.len;
                ++cur_frag;
            }

            prev_end = cur_end;
        }

        return {std::move(lr), std::move(core_seqs)};
    }

    static auto build_rev_lr_cores(
        const gtl::vector<SequenceContainerType> &seqs,
        const SyncmerType &syncmer_functor, int lcp_level, lps &str)
        -> lcp_read {
        lcp_read lr;
        lr.clear();

        gtl::vector<SequenceContainerType> rev_seqs;
        for (auto it = seqs.rbegin(); it != seqs.rend(); ++it) {
            rev_seqs.emplace_back(it->reverse_complement());
        }
        init_lps_crumbled(&str, rev_seqs, syncmer_functor);
        lps_deepen(&str, lcp_level);

        for (int i = 0; i < str.size; ++i) {
            lr.cores.emplace_back(str.cores[i].label);
            lr.pos.emplace_back(str.cores[i].start, str.cores[i].end);
        }
        return lr;
    }

    static auto build_rev_lr_cores_with_seqs(
        const gtl::vector<SequenceContainerType> &seqs,
        const SyncmerType &syncmer_functor, int lcp_level, lps &str)
        -> std::pair<lcp_read, gtl::vector<SequenceContainerType>> {
        lcp_read lr;
        lr.clear();
        gtl::vector<SequenceContainerType> core_seqs;

        gtl::vector<SequenceContainerType> rev_seqs;
        rev_seqs.reserve(seqs.size());
        for (auto it = seqs.rbegin(); it != seqs.rend(); ++it) {
            rev_seqs.emplace_back(it->reverse_complement());
        }

        init_lps_crumbled(&str, rev_seqs, syncmer_functor);
        lps_deepen(&str, lcp_level);

        if (str.size == 0 || rev_seqs.empty()) {
            return {std::move(lr), std::move(core_seqs)};
        }

        lr.cores.reserve(str.size);
        lr.pos.reserve(str.size);
        core_seqs.reserve(str.size);

        // Co-iterate: track current fragment and its cumulative position
        size_t frag_idx = 0;
        size_t frag_cumulative = 0;  // start position of current fragment

        size_t prev_end = 0;  // End position of previous core

        for (int i = 0; i < str.size; ++i) {
            size_t cur_start = str.cores[i].start;
            size_t cur_end = str.cores[i].end;
            lr.cores.emplace_back(str.cores[i].label);
            lr.pos.emplace_back(cur_start, cur_end);

            // Extract non-overlapping portion: [prev_end, cur_end)
            // For the first core, this is [0, cur_end)
            size_t from = prev_end;
            size_t to = cur_end;
            size_t sub_len = to - from;

            core_seqs.emplace_back(sub_len, 0);
            auto &result = core_seqs.back();
            size_t dest_pos = 0;

            // Advance frag_idx until we reach the fragment containing 'from'
            while (frag_idx < rev_seqs.size() &&
                   frag_cumulative + rev_seqs[frag_idx].len <= from) {
                frag_cumulative += rev_seqs[frag_idx].len;
                ++frag_idx;
            }

            // Extract bases from fragments, starting from current position
            size_t cur_frag = frag_idx;
            size_t cur_cumulative = frag_cumulative;

            while (cur_frag < rev_seqs.size() && dest_pos < sub_len) {
                const auto &frag = rev_seqs[cur_frag];
                size_t frag_start = cur_cumulative;
                size_t frag_end = cur_cumulative + frag.len;

                // Compute overlap of [frag_start, frag_end) with [from, to)
                size_t overlap_start = std::max(frag_start, from);
                size_t overlap_end = std::min(frag_end, to);

                if (overlap_start < overlap_end) {
                    size_t local_from = overlap_start - frag_start;
                    size_t copy_len = overlap_end - overlap_start;

                    // Copy bases from frag to result
                    for (size_t b = 0; b < copy_len; ++b) {
                        size_t src_bit =
                            (frag.skip_first_n + local_from + b) * 2;
                        uint64_t code =
                            (frag.data[src_bit >> 6] >> (src_bit & 63)) & 0x3;
                        size_t dst_bit = dest_pos * 2;
                        result.data[dst_bit >> 6] |= code << (dst_bit & 63);
                        ++dest_pos;
                    }
                    result.origin_node = frag.origin_node;
                }
                cur_cumulative += frag.len;
                ++cur_frag;
            }

            prev_end = cur_end;
        }

        return {std::move(lr), std::move(core_seqs)};
    }

    auto fold_dbg_graph(dbg_index &dbgi) {
        fprintf(stderr, "dbg_size: %zu\n", dbgi.dbg.size());

        // Discard a reciprocal match that spans less than this fraction (1/N) of
        // the query node's cores — too short to be a reliable overlap.
        constexpr size_t MIN_OVERLAP_CORE_FRACTION_DIV = 20;

        std::atomic<int> all_n{0};
        std::atomic<int> hits{0};
        std::atomic<int> bad{0};
        // Phase-0 matching diagnostics: why do good nodes fail to match?
        std::atomic<int> too_short{0};   // < 4*olc_lcmer_l cores → no ranges (#4)
        std::atomic<int> no_ranges{0};   // long enough but 0 lcmer hits
        std::atomic<int> frac_drop{0};   // overlap < 5% of cores, filtered (#6)
        gtl::parallel_flat_hash_map_d<u64, gtl::vector<SelfMatch>, HT_SHARD_BC,
                                      std::mutex>
            matched;
        for_each_map(
            dbgi.dbg,
            [&](int tid, const auto &node) {
                (void)tid;
                all_n.fetch_add(1);
                if (dbgi.active_set.find(node.first) == dbgi.active_set.end()) {
                    bad.fetch_add(1);
                    return;
                }

                lps str = {0, 0, 0, 0};
                lcp_read lr;
                // dbgi.dbg.for_each([&](auto &node) {
                if (node.second.redirected()) {
                    bad.fetch_add(1);
                    return;
                }
                // good_n.fetch_add(1);
                lr.clear();
                gtl::vector<SequenceContainerType> rev_seqs;
                for (auto it = node.second.seq.rbegin();
                     it != node.second.seq.rend(); ++it) {
                    rev_seqs.emplace_back(it->reverse_complement());
                }
                init_lps_crumbled(&str, rev_seqs, syncmer_functor);
                lps_deepen(&str, lcp_level);

                for (int i = 0; i < str.size; ++i) {
                    lr.cores.emplace_back(str.cores[i].label);
                    lr.pos.emplace_back(str.cores[i].start, str.cores[i].end);
                }
                if (lr.cores.size() < static_cast<size_t>(olc_lcmer_l + 1)) {
                    bad.fetch_add(1);
                    // good_n.fetch_sub(1);
                    free_lps(&str);
                    return;
                }
                // Keep everything in CORE-INDEX space for layouting.
                // `find_lcmer_ranges()` already returns [from_core,to_core] in
                // core indices (inclusive), which is the correct coordinate
                // system for connected_comp_plus_layout().
                // Matching diagnostic: bucket why a good node yields no match.
                // Mirrors find_lcmer_ranges' gate via the same OLC_SKIP_* macros.
                if (lr.cores.size() <
                    static_cast<size_t>(olc_lcmer_l + OLC_SKIP_FIRST(olc_lcmer_l) +
                                        OLC_SKIP_LAST(olc_lcmer_l))) {
                    too_short.fetch_add(1);  // fails find_lcmer_ranges gate (#4)
                }
                auto found_ranges = find_lcmer_ranges(lr, dbgi, olc_lcmer_l);
                bool any_hit = false;
                for (const auto &[id, from_core, to_core, from_lcmer,
                                  to_lcmer] : found_ranges) {
                    if (id != 0) {
                        if (to_core - from_core <
                            lr.pos.size() / MIN_OVERLAP_CORE_FRACTION_DIV) {
                            frac_drop.fetch_add(1);
                            continue;
                        }
                        any_hit = true;
                        hits.fetch_add(1);

                        // Store the overlap on THIS node in CORE indices
                        // (inclusive; see CoreSpan). Base-space scaling happens
                        // only at emission time. Core-space length is the number
                        // of steps between endpoints: (cores.size() - 1).
                        const u64 core_len_steps =
                            (lr.cores.size() > 0)
                                ? static_cast<u64>(lr.cores.size() - 1)
                                : 0;

                        matched.lazy_emplace_l(
                            node.first,
                            [&](auto &val) {
                                val.second.emplace_back(
                                    id,
                                    CoreSpan{static_cast<u64>(from_core),
                                             static_cast<u64>(to_core),
                                             core_len_steps});
                            },
                            [&](const auto &ctr) {
                                ctr(node.first,
                                    gtl::vector<SelfMatch>{SelfMatch{
                                        id,
                                        CoreSpan{static_cast<u64>(from_core),
                                                 static_cast<u64>(to_core),
                                                 core_len_steps}}});
                            });
                    }
                }
                if (!any_hit &&
                    lr.cores.size() >=
                        static_cast<size_t>(olc_lcmer_l +
                                            OLC_SKIP_FIRST(olc_lcmer_l) +
                                            OLC_SKIP_LAST(olc_lcmer_l)))
                    no_ranges.fetch_add(1);
                free_lps(&str);
            },
            threads.work);

        gtl::btree_multimap<u64, Overlap> reciprocal;
        int non_recip = 0, recip_kept = 0;
        for (const auto &[k, v] : matched) {
            for (const auto &match : v) {
                // Keep the overlap only if the partner records the reverse one
                // too (a reciprocal match). `match.self` is k's core span in the
                // overlap; `ff->self` is the partner's. Both are core-index
                // spans (see CoreSpan); base-space scaling happens at emission.
                auto it = matched.find(match.other);
                auto ff = (it != matched.end())
                              ? std::find_if(it->second.begin(),
                                             it->second.end(),
                                             [k](const auto &val) {
                                                 return val.other == k;
                                             })
                              : gtl::vector<SelfMatch>::const_iterator{};
                if (it != matched.end() && ff != it->second.end()) {
                    reciprocal.emplace(k,
                                       Overlap{match.other, match.self, ff->self});
                    ++recip_kept;
                } else {
                    ++non_recip;  // one-directional match (#7)
#if OLC_NONRECIPROCAL
                    // Prototype: keep it anyway. We only know k's span
                    // (match.self); approximate the partner's span as the same
                    // shape (self-symmetric) so connected_comp_plus_layout can
                    // still position the partner. Coarse — it measures whether
                    // wiring these extra nodes in helps before we compute the
                    // partner's real overlap interval.
                    reciprocal.emplace(
                        k, Overlap{match.other, match.self, match.self});
#endif
                }
            }
        }
        fprintf(stderr, "all:%d\ngood: %d nodes\nhits: %d\n", all_n.load(),
                all_n.load() - bad.load(), hits.load());
        fprintf(stderr,
                "[match-diag] too_short:%d no_ranges:%d frac_drop:%d "
                "recip_kept:%d non_recip:%d\n",
                too_short.load(), no_ranges.load(), frac_drop.load(), recip_kept,
                non_recip);
        auto layouts =
            layout_graph_algo::connected_comp_plus_layout(reciprocal);

        return std::make_pair(reciprocal, layouts);
    }
};

// TODO do not use extensions to identify input type.
input_type identify_input(const std::string &input) {
    std::string_view sv(input);
    if (sv.ends_with(".bam")) {
        return input_type::BAM;
    } else if (sv.ends_with(".sam")) {
        return input_type::SAM;
    } else if (sv.ends_with(".fastq") || sv.ends_with(".fq") ||
               sv.ends_with(".fasta") || sv.ends_with(".fa") ||
               sv.ends_with(".fna")) {
        return input_type::FASTX;
    } else if (sv.ends_with(".fastq.gz") || sv.ends_with(".fq.gz") ||
               sv.ends_with(".fasta.gz") || sv.ends_with(".fa.gz") ||
               sv.ends_with(".fna.gz")) {
        return input_type::FASTXGZ;
    } else if (sv.ends_with(".fastq.bz2") || sv.ends_with(".fq.bz2") ||
               sv.ends_with(".fasta.bz2") || sv.ends_with(".fa.bz2") ||
               sv.ends_with(".fna.bz2")) {
        return input_type::FASTXBZ2;
    }
    return input_type::UNKNOWN;
}

static bool validate_path(const std::string &path) {
    auto slash = path.rfind('/');
    if (slash == std::string::npos) return true;
    std::string dir = path.substr(0, slash);
    struct stat st;
    if (stat(dir.c_str(), &st) != 0 || !S_ISDIR(st.st_mode)) {
        fprintf(stderr, "error: output directory does not exist: %s\n",
                dir.c_str());
        return false;
    }
    return true;
}

int run_assembly(const asm_config &cfg) {
    if (!validate_path(cfg.output)) return 1;

    // Validate the syncmer configuration once, up front, before any syncmer
    // functor is constructed (the functor ctor assumes a valid mask/mode).
    if (auto v = cfg.syncmer.validate(); !v) {
        fprintf(stderr, "Error: %s\n", v.error().c_str());
        return 1;
    }

    if (!cfg.silent) cfg.print_cfg(stderr);
    // free(malloc(1024ULL * 1024 * 1024 * 20)); INTIMIDATE the OS for the
    // ram

    INIT_TIMERS(HIFI_ASSEMBLY);

    assembler as(cfg);
    // The DBG index lives for the whole program; it is held in stack-resident
    // inline storage and its (expensive) destructor is intentionally skipped at
    // exit. See no_destructor.h.
    NoDestructor<dbg_index> dbg_storage(cfg.bloom.false_positive_rate,
                                        cfg.bloom.projected_count, cfg.seed);
    dbg_index &dbg = *dbg_storage;
    if (auto r = as.construct_dbg(dbg); !r) {
        fprintf(stderr, "%s Terminating.\n", r.error().c_str());
        return 1;
    }
    IF_TIME_CHECKPOINT(!cfg.silent, HIFI_ASSEMBLY, stderr, "DBG construction");

    dbg.count_forward_edges();
    IF_TIME_CHECKPOINT(!cfg.silent, HIFI_ASSEMBLY, stderr,
                       "Count forward edges");

    if (cfg.print_precleanup_graph) {
        const std::string raw_path = cfg.output + ".raw.gfa";
        FILE *out = fopen(raw_path.c_str(), "w");
        if (!out) {
            fprintf(stderr, "Cannot open output: %s\n", raw_path.c_str());
            return 1;
        }
        print_dbg_seqs(out, dbg);
        fclose(out);

        IF_TIME_CHECKPOINT(!cfg.silent, HIFI_ASSEMBLY, stderr, "DBG printing");
    }

    GraphCleaning::run(dbg, cfg);

    IF_TIME_CHECKPOINT(!cfg.silent, HIFI_ASSEMBLY, stderr, "Graph cleaning");

    if (cfg.thread_bridge.enable) {
        if (auto r = as.run_thread_bridge(dbg, cfg.output + ".threads.tsv");
            !r) {
            fprintf(stderr, "%s\n", r.error().c_str());
            return 1;
        }
        IF_TIME_CHECKPOINT(!cfg.silent, HIFI_ASSEMBLY, stderr, "Thread-bridge");
    }


    // Dump pre-layout (bidirected) unitig graph, gated on --print-bi-dbg.
    if (cfg.print_prelayout_graph) {
        const std::string pre_path = cfg.output + ".bi.gfa";
        FILE *pre = fopen(pre_path.c_str(), "w");
        if (!pre) {
            fprintf(stderr, "Cannot open output: %s\n", pre_path.c_str());
            return 1;
        }
        print_dbg_seqs(pre, dbg);
        fclose(pre);
    }

    {
        FILE *out = fopen(cfg.output.c_str(), "w");
        if (!out) {
            fprintf(stderr, "Cannot open output: %s\n", cfg.output.c_str());
            return 1;
        }
        // print_dbg_seqs(out, dbg);
        auto emit_result =
            layout_graph_emit::emit_ultra_merged_graph(as, dbg, out);
        fclose(out);
        if (!emit_result) {
            fprintf(stderr, "%s\n", emit_result.error().c_str());
            return 1;
        }
    }

    IF_TIME_CHECKPOINT(!cfg.silent, HIFI_ASSEMBLY, stderr, "Merged/Printed");

    return 0;
}

void tests_mock_fill_seq(dbg_index &dbgi) {
    for (auto &node : dbgi.dbg) {
        node.second.seq.emplace_back();
    }
}

TEST_CASE("build_lr_cores_with_seqs extracts correct sequences") {
    using SYNCMER_HASHER = syncmer::MurMur;
    using SyncmerType = syncmer::GenericSubmerRule<uint32_t, SYNCMER_HASHER>;
    decltype(asm_config::syncmer) spara;
    spara.k = 8;
    spara.density = 0.05;
    SyncmerType syncmer_functor(spara, 42);
    int lcp_level = 3;

    // Create test sequences - must be long enough and varied enough to produce
    // cores
    const char *seq1 =
        "ACGTACGTAAGGCCTTAACCGGTTAAGGCCTTAACCGGTTACGTACGTAAGGCCTT"
        "TGCATGCAGGAATTCCGGAATTCCGGAATTCCTGCATGCAGGAATTCCGGAATTCC"
        "AACCGGTTAACCGGTTAACCGGTTAACCGGTTAACCGGTTAACCGGTTAACCGGTT"
        "GGCCAATTGGCCAATTGGCCAATTGGCCAATTGGCCAATTGGCCAATTGGCCAATT";
    const char *seq2 =
        "TGCATGCATGCAACGTACGTACGTTGCATGCATGCAACGTACGTACGTTGCATGCA"
        "GGAATTCCGGAATTCCAACCGGTTGGAATTCCGGAATTCCAACCGGTTGGAATTCC"
        "AACCGGTTAACCGGTTGGCCAATTAACCGGTTAACCGGTTGGCCAATTAACCGGTT"
        "TTGGCCAATTGGCCAACCGGTTAATTGGCCAATTGGCCAACCGGTTAATTGGCCAA";
    const char *seq3 =
        "GGCCGGCCAATTAATTGGCCGGCCAATTAATTGGCCGGCCAATTAATTGGCCGGCC"
        "ACGTACGTAACCGGTTACGTACGTAACCGGTTACGTACGTAACCGGTTACGTACGT"
        "TGCATGCAGGAATTCCTGCATGCAGGAATTCCTGCATGCAGGAATTCCTGCATGCA"
        "AATTGGCCAATTGGCCAATTGGCCAATTGGCCAATTGGCCAATTGGCCAATTGGCC";

    gtl::vector<SequenceContainerType> seqs;
    seqs.emplace_back(seq1, 0, strlen(seq1), 0,
                      SequenceContainerType::from_ascii{});
    seqs.emplace_back(seq2, 0, strlen(seq2), 0,
                      SequenceContainerType::from_ascii{});
    seqs.emplace_back(seq3, 0, strlen(seq3), 0,
                      SequenceContainerType::from_ascii{});

    // Build concatenated reference string
    std::string concat_ref =
        std::string(seq1) + std::string(seq2) + std::string(seq3);

    SUBCASE("build_lr_cores_with_seqs") {
        lps lps_buf = {0, 0, 0, 0};
        auto [lr, core_seqs] = assembler::build_lr_cores_with_seqs(
            seqs, syncmer_functor, lcp_level, lps_buf);
        auto lr_original =
            assembler::build_lr_cores(seqs, syncmer_functor, lcp_level, lps_buf);
        free_lps(&lps_buf);

        // Ensure we actually have cores to test
        REQUIRE(lr.cores.size() > 0);

        // Check that cores match the original function
        REQUIRE(lr.cores.size() == lr_original.cores.size());
        REQUIRE(lr.pos.size() == lr_original.pos.size());
        for (size_t i = 0; i < lr.cores.size(); ++i) {
            CHECK_EQ(lr.cores[i], lr_original.cores[i]);
            CHECK_EQ(lr.pos[i].first, lr_original.pos[i].first);
            CHECK_EQ(lr.pos[i].second, lr_original.pos[i].second);
        }

        REQUIRE(lr.cores.size() == core_seqs.size());
        REQUIRE(lr.pos.size() == core_seqs.size());

        // Verify non-overlapping extraction: core_seqs[i] should be [prev_end,
        // cur_end)
        char buffer[2048];
        size_t prev_end = 0;
        for (size_t i = 0; i < core_seqs.size(); ++i) {
            size_t cur_end = lr.pos[i].second;
            // Expected is the non-overlapping portion: [prev_end, cur_end)
            std::string expected =
                concat_ref.substr(prev_end, cur_end - prev_end);

            core_seqs[i].convert_to_ascii(buffer);
            std::string actual(buffer);

            INFO("Core ", i, " prev_end=", prev_end, " cur_end=", cur_end);
            CHECK_EQ(expected, actual);

            prev_end = cur_end;
        }
    }

    SUBCASE("build_rev_lr_cores_with_seqs") {
        lps lps_buf = {0, 0, 0, 0};
        auto [lr, core_seqs] = assembler::build_rev_lr_cores_with_seqs(
            seqs, syncmer_functor, lcp_level, lps_buf);
        auto lr_original =
            assembler::build_rev_lr_cores(seqs, syncmer_functor, lcp_level, lps_buf);
        free_lps(&lps_buf);

        // Ensure we actually have cores to test
        REQUIRE(lr.cores.size() > 0);

        // Check that cores match the original function
        REQUIRE(lr.cores.size() == lr_original.cores.size());
        REQUIRE(lr.pos.size() == lr_original.pos.size());
        for (size_t i = 0; i < lr.cores.size(); ++i) {
            CHECK_EQ(lr.cores[i], lr_original.cores[i]);
            CHECK_EQ(lr.pos[i].first, lr_original.pos[i].first);
            CHECK_EQ(lr.pos[i].second, lr_original.pos[i].second);
        }

        // Build reverse-complemented concatenated reference
        gtl::vector<SequenceContainerType> rev_seqs;
        for (auto it = seqs.rbegin(); it != seqs.rend(); ++it) {
            rev_seqs.emplace_back(it->reverse_complement());
        }
        std::string rev_concat_ref;
        char temp_buf[2048];
        for (const auto &rs : rev_seqs) {
            rs.convert_to_ascii(temp_buf);
            rev_concat_ref += std::string(temp_buf);
        }

        REQUIRE(lr.cores.size() == core_seqs.size());
        REQUIRE(lr.pos.size() == core_seqs.size());

        // Verify non-overlapping extraction: core_seqs[i] should be [prev_end,
        // cur_end)
        char buffer[2048];
        size_t prev_end = 0;
        for (size_t i = 0; i < core_seqs.size(); ++i) {
            size_t cur_end = lr.pos[i].second;
            // Expected is the non-overlapping portion: [prev_end, cur_end)
            std::string expected =
                rev_concat_ref.substr(prev_end, cur_end - prev_end);

            core_seqs[i].convert_to_ascii(buffer);
            std::string actual(buffer);

            INFO("Core ", i, " prev_end=", prev_end, " cur_end=", cur_end);
            CHECK_EQ(expected, actual);

            prev_end = cur_end;
        }
    }

    SUBCASE("core sequences span fragment boundaries correctly") {
        // Build cores from vector of 3 sequences
        lps lps_buf = {0, 0, 0, 0};
        auto [lr, core_seqs] = assembler::build_lr_cores_with_seqs(
            seqs, syncmer_functor, lcp_level, lps_buf);
        free_lps(&lps_buf);

        // Ensure we actually have cores to test
        REQUIRE(lr.cores.size() > 0);

        // Get fragment boundaries
        size_t frag1_end = strlen(seq1);
        size_t frag2_end = frag1_end + strlen(seq2);

        // Find non-overlapping portions that span fragment boundaries
        [[maybe_unused]] bool found_spanning_seq = false;
        char buffer[2048];
        size_t prev_end = 0;
        for (size_t i = 0; i < core_seqs.size(); ++i) {
            size_t cur_end = lr.pos[i].second;

            // Check if this non-overlapping portion [prev_end, cur_end) spans a
            // fragment boundary
            bool spans_boundary =
                (prev_end < frag1_end && cur_end > frag1_end) ||
                (prev_end < frag2_end && cur_end > frag2_end);
            if (spans_boundary) {
                found_spanning_seq = true;

                // Verify the extracted sequence matches the non-overlapping
                // portion
                std::string expected =
                    concat_ref.substr(prev_end, cur_end - prev_end);
                core_seqs[i].convert_to_ascii(buffer);
                std::string actual(buffer);

                INFO("Spanning seq ", i, " prev_end=", prev_end,
                     " cur_end=", cur_end, " frag1_end=", frag1_end,
                     " frag2_end=", frag2_end);
                CHECK_EQ(expected, actual);
            }

            prev_end = cur_end;
        }

        // We expect at least one sequence to span fragment boundaries with
        // these test sequences
        INFO("Total cores: ", lr.cores.size());
        CHECK(found_spanning_seq);
    }

    SUBCASE("rev: core sequences span fragment boundaries correctly") {
        // Build cores from vector of 3 sequences (reverse complemented)
        lps lps_buf = {0, 0, 0, 0};
        auto [lr, core_seqs] = assembler::build_rev_lr_cores_with_seqs(
            seqs, syncmer_functor, lcp_level, lps_buf);
        free_lps(&lps_buf);

        // Build reverse-complemented concatenated reference
        gtl::vector<SequenceContainerType> rev_seqs;
        for (auto it = seqs.rbegin(); it != seqs.rend(); ++it) {
            rev_seqs.emplace_back(it->reverse_complement());
        }
        std::string rev_concat_ref;
        char temp_buf[2048];
        for (const auto &rs : rev_seqs) {
            rs.convert_to_ascii(temp_buf);
            rev_concat_ref += std::string(temp_buf);
        }

        // Ensure we actually have cores to test
        REQUIRE(lr.cores.size() > 0);

        // Get fragment boundaries (in rev_seqs order: seq3_rc, seq2_rc,
        // seq1_rc)
        size_t frag1_end = strlen(seq3);  // seq3 is first in reverse order
        size_t frag2_end = frag1_end + strlen(seq2);

        // Find non-overlapping portions that span fragment boundaries
        [[maybe_unused]] bool found_spanning_seq = false;
        char buffer[2048];
        size_t prev_end = 0;
        for (size_t i = 0; i < core_seqs.size(); ++i) {
            size_t cur_end = lr.pos[i].second;

            // Check if this non-overlapping portion [prev_end, cur_end) spans a
            // fragment boundary
            bool spans_boundary =
                (prev_end < frag1_end && cur_end > frag1_end) ||
                (prev_end < frag2_end && cur_end > frag2_end);
            if (spans_boundary) {
                found_spanning_seq = true;

                // Verify the extracted sequence matches the non-overlapping
                // portion
                std::string expected =
                    rev_concat_ref.substr(prev_end, cur_end - prev_end);
                core_seqs[i].convert_to_ascii(buffer);
                std::string actual(buffer);

                INFO("Spanning seq ", i, " prev_end=", prev_end,
                     " cur_end=", cur_end, " frag1_end=", frag1_end,
                     " frag2_end=", frag2_end);
                CHECK_EQ(expected, actual);
            }

            prev_end = cur_end;
        }

        // We expect at least one sequence to span fragment boundaries with
        // these test sequences
        INFO("Total cores: ", lr.cores.size());
        CHECK(found_spanning_seq);
    }
}

TEST_CASE(
    "CrumbledSubSeq64 iterator-range constructor concatenates represented "
    "sequence") {
    const char *src = "TTTTACGTACGTGGGGCCAAATTTCCGG";
    gtl::vector<SequenceContainerType> frags;
    frags.emplace_back(src, 3, 9, 101, SequenceContainerType::from_ascii{});
    frags.emplace_back(src, 10, 18, 102, SequenceContainerType::from_ascii{});
    frags.emplace_back(src, 20, 27, 103, SequenceContainerType::from_ascii{});

    auto to_ascii = [](const SequenceContainerType &s) {
        std::string out;
        out.resize(s.len);
        if (!out.empty()) {
            s.convert_to_ascii(out.data());
        }
        return out;
    };

    SUBCASE("full range") {
        SequenceContainerType merged(frags.begin(), frags.end());
        const std::string expected = std::string(src + 3, 9 - 3) +
                                     std::string(src + 10, 18 - 10) +
                                     std::string(src + 20, 27 - 20);
        CHECK_EQ(merged.len, expected.size());
        CHECK_EQ(to_ascii(merged), expected);
        CHECK_EQ(merged.skip_first_n, 0);
        CHECK_EQ(merged.origin_node, 103);
    }

    SUBCASE("sub range") {
        SequenceContainerType merged(frags.begin() + 1, frags.end());
        const std::string expected =
            std::string(src + 10, 18 - 10) + std::string(src + 20, 27 - 20);
        CHECK_EQ(merged.len, expected.size());
        CHECK_EQ(to_ascii(merged), expected);
        CHECK_EQ(merged.origin_node, 103);
    }

    SUBCASE("empty range") {
        SequenceContainerType merged(frags.begin(), frags.begin());
        CHECK_EQ(merged.len, 0);
        CHECK_EQ(merged.skip_first_n, 0);
        CHECK_EQ(merged.origin_node, 0);
    }
}

// ── Characterization tests: layout_graph_emit core algorithms ──────────
// These pin down the observable behavior of the seed-chain DAG pipeline
// (build_scd_graph → find_heaviest_paths → collect_path_sequences) so the
// planned complexity/quality refactors (docs/layout-graph-plan.md) can be
// proven behavior-preserving. They exercise the pure struct-based functions
// in isolation, without the assembler.
TEST_CASE("layout_graph_emit::build_scd_graph merges shared cores") {
    using namespace layout_graph_emit;

    // Two contigs sharing cores 20 and 30 (contig0: 10,20,30; contig1: 20,30,40).
    // contig ids double as indices into `nodes`.
    gtl::vector<lane_node> nodes(2);
    nodes[0] = {/*id*/ 0, /*from*/ 0, /*to*/ 100, /*fwd*/ true, nullptr};
    nodes[1] = {/*id*/ 1, /*from*/ 10, /*to*/ 110, /*fwd*/ true, nullptr};

    gtl::flat_hash_map<size_t, lcp_read> contig_lr;
    lcp_read lr0;
    lr0.cores = {10u, 20u, 30u};
    lr0.pos = {{0u, 1u}, {1u, 2u}, {2u, 3u}};
    lcp_read lr1;
    lr1.cores = {20u, 30u, 40u};
    lr1.pos = {{0u, 1u}, {1u, 2u}, {2u, 3u}};
    contig_lr[0] = std::move(lr0);
    contig_lr[1] = std::move(lr1);

    gtl::vector<size_t> contig_order = {0, 1};

    gtl::vector<scd_node> scd;
    gtl::vector<scd_adj_t> scd_fwd;
    gtl::flat_hash_map<uint32_t, gtl::vector<uint32_t>> label_nodes;
    gtl::vector<gtl::vector<std::pair<size_t, size_t>>> node_contigs;

    build_scd_graph(contig_order, contig_lr, nodes, scd, scd_fwd, label_nodes,
                    node_contigs);

    // contig0 creates 3 nodes; contig1 merges 20 & 30, adds only 40 → 4 total.
    CHECK_EQ(scd.size(), 4u);
    CHECK_EQ(scd_fwd.size(), scd.size());
    CHECK_EQ(node_contigs.size(), scd.size());

    // Weights: shared cores seen twice, private cores once.
    auto weight_of_label = [&](uint32_t label) -> int {
        for (const auto &n : scd)
            if (n.core_label == label) return n.weight;
        return -1;
    };
    CHECK_EQ(weight_of_label(10u), 1);
    CHECK_EQ(weight_of_label(20u), 2);
    CHECK_EQ(weight_of_label(30u), 2);
    CHECK_EQ(weight_of_label(40u), 1);

    // Merged nodes carry both contigs' provenance.
    for (size_t idx = 0; idx < scd.size(); ++idx) {
        if (scd[idx].core_label == 20u || scd[idx].core_label == 30u)
            CHECK_EQ(node_contigs[idx].size(), 2u);
    }
}

TEST_CASE("layout_graph_emit::find_heaviest_paths follows the heaviest chain") {
    using namespace layout_graph_emit;

    // Graph: 0→1→2 (main chain) plus a light branch 0→3.
    // All weights 1 → main-chain sink (2, dp=3) beats branch sink (3, dp=2).
    gtl::vector<lane_node> nodes(1);
    nodes[0] = {/*id*/ 0, /*from*/ 0, /*to*/ 100, /*fwd*/ true, nullptr};

    gtl::vector<scd_node> scd(4);
    for (auto &n : scd) { n.contig_idx = 0; n.core_idx = 0; n.weight = 1; }
    scd[0].core_label = 0;
    scd[1].core_label = 1;
    scd[2].core_label = 2;
    scd[3].core_label = 3;

    gtl::vector<scd_adj_t> scd_fwd(4);
    scd_fwd[0].push_back({1u, 1});
    scd_fwd[1].push_back({2u, 1});
    scd_fwd[0].push_back({3u, 1});

    auto paths = find_heaviest_paths(scd, scd_fwd, nodes, /*layout_idx*/ 0);

    REQUIRE_EQ(paths.size(), 1u);         // single connected component
    REQUIRE_EQ(paths[0].path.size(), 3u);
    CHECK_EQ(paths[0].path[0], 0u);
    CHECK_EQ(paths[0].path[1], 1u);
    CHECK_EQ(paths[0].path[2], 2u);       // heaviest sink, not the branch
    CHECK_EQ(paths[0].total_weight, 3);
}

TEST_CASE("layout_graph_emit::collect_path_sequences concatenates with gaps") {
    using namespace layout_graph_emit;

    auto make_seq = [](const char *s) {
        return SequenceContainerType(s, 0, std::strlen(s), 0,
                                     SequenceContainerType::from_ascii{});
    };

    // node_contigs: node0→(contig0,core0), node1→(contig0,core1),
    //               node2→(contig1,core0)
    gtl::vector<gtl::vector<std::pair<size_t, size_t>>> node_contigs(3);
    node_contigs[0] = {{0, 0}};
    node_contigs[1] = {{0, 1}};
    node_contigs[2] = {{1, 0}};

    gtl::flat_hash_map<size_t, gtl::vector<SequenceContainerType>>
        contig_core_seqs;
    contig_core_seqs[0].push_back(make_seq("ACGT"));      // len 4
    contig_core_seqs[0].push_back(make_seq("TTTT"));      // len 4
    contig_core_seqs[1].push_back(make_seq("GGGGGG"));    // len 6

    SUBCASE("single component, path order preserved") {
        gtl::vector<comp_path> paths;
        paths.push_back({0, 100, 2, gtl::vector<uint32_t>{0u, 1u}});

        auto comps = collect_path_sequences(paths, node_contigs,
                                            contig_core_seqs, 3);
        REQUIRE_EQ(comps.size(), 1u);
        REQUIRE_EQ(comps[0].size(), 2u);
        CHECK_EQ(comps[0][0].len, 4u);
        CHECK_EQ(comps[0][1].len, 4u);
    }

    SUBCASE("two components stay separate — never gap-joined (R7)") {
        gtl::vector<comp_path> paths;
        paths.push_back({0, 100, 2, gtl::vector<uint32_t>{0u, 1u}});
        paths.push_back({200, 300, 1, gtl::vector<uint32_t>{2u}});

        auto comps = collect_path_sequences(paths, node_contigs,
                                            contig_core_seqs, 3);
        // Two independent groups, emitted as two contigs. No sentinel, no
        // fabricated adjacency between disconnected components.
        REQUIRE_EQ(comps.size(), 2u);
        REQUIRE_EQ(comps[0].size(), 2u);
        CHECK_EQ(comps[0][0].len, 4u);
        CHECK_EQ(comps[0][1].len, 4u);
        REQUIRE_EQ(comps[1].size(), 1u);
        CHECK_EQ(comps[1][0].len, 6u);
        for (const auto &comp : comps)
            for ([[maybe_unused]] const auto &s : comp)
                CHECK(s.len > 0);  // no zero-length gap sentinels remain
    }

    SUBCASE("fallback picks longest representative, not arbitrary first (C12)") {
        // A single-node path whose node has two candidate contigs and no
        // continuity hint. The fallback must pick the longer core sequence
        // (contig 3, len 6), not the first-listed one (contig 5, len 2).
        gtl::vector<gtl::vector<std::pair<size_t, size_t>>> nc(1);
        nc[0] = {{5, 0}, {3, 0}};  // first-listed (5) is the SHORT one

        gtl::flat_hash_map<size_t, gtl::vector<SequenceContainerType>> ccs;
        ccs[5].push_back(make_seq("AC"));      // len 2 (arbitrary nc[0] choice)
        ccs[3].push_back(make_seq("ACGTAC"));  // len 6 (longest representative)

        gtl::vector<comp_path> paths;
        paths.push_back({0, 100, 1, gtl::vector<uint32_t>{0u}});

        auto comps = collect_path_sequences(paths, nc, ccs, 1);
        REQUIRE_EQ(comps.size(), 1u);
        REQUIRE_EQ(comps[0].size(), 1u);
        CHECK_EQ(comps[0][0].len, 6u);  // longest, not nc[0]'s len 2
    }
}

TEST_CASE("layout_graph_emit::resolve_overlapping_components drops duplicates") {
    using namespace layout_graph_emit;
    auto P = [](i64 from, i64 to, int64_t w) {
        return comp_path{from, to, w, gtl::vector<uint32_t>{}};
    };

    SUBCASE("nested lighter component dropped, heavier kept (M79 class)") {
        // Real M79 spans: comp0 [0,86] w153 fully contains comp1 [1,15] w17.
        gtl::vector<comp_path> cps{P(0, 86, 153), P(1, 15, 17)};
        resolve_overlapping_components(cps);
        REQUIRE_EQ(cps.size(), 1u);
        CHECK_EQ(cps[0].layout_from, 0);
        CHECK_EQ(cps[0].total_weight, 153);
    }

    SUBCASE("disjoint components both kept (no-op on the pre-sorted input)") {
        // find_heaviest_paths sorts by layout_from before calling this; the
        // no-drop path preserves that order.
        gtl::vector<comp_path> cps{P(0, 100, 2), P(200, 300, 1)};
        resolve_overlapping_components(cps);
        REQUIRE_EQ(cps.size(), 2u);
        CHECK_EQ(cps[0].layout_from, 0);
        CHECK_EQ(cps[1].layout_from, 200);
    }

    SUBCASE("small overlap (<=50% of shorter) keeps both") {
        // shorter interval len 100; overlap [90,100] = 10 <= 50 ⇒ both kept.
        gtl::vector<comp_path> cps{P(0, 100, 5), P(90, 190, 4)};
        resolve_overlapping_components(cps);
        CHECK_EQ(cps.size(), 2u);
    }

    SUBCASE("major overlap (>50% of shorter) keeps the heavier one") {
        // shorter interval len 100; overlap [10,100] = 90 > 50 ⇒ drop lighter.
        gtl::vector<comp_path> cps{P(0, 100, 4), P(10, 200, 9)};
        resolve_overlapping_components(cps);
        REQUIRE_EQ(cps.size(), 1u);
        CHECK_EQ(cps[0].total_weight, 9);   // heavier survives regardless of order
    }

    SUBCASE("no-op when nothing overlaps or fewer than two") {
        gtl::vector<comp_path> one{P(0, 100, 1)};
        resolve_overlapping_components(one);
        CHECK_EQ(one.size(), 1u);
        gtl::vector<comp_path> none;
        resolve_overlapping_components(none);
        CHECK_EQ(none.size(), 0u);
    }
}

TEST_CASE("layout_graph_emit::PostLayoutCleaning::cut_tips peels short tips") {
    using namespace layout_graph_emit;

    auto add_node = [](supernode_vec &sn, const char *name, bool layout,
                       size_t len) {
        supernode_out s;
        s.name = name;
        s.is_layout = layout;
        s.ref_seqs = nullptr;
        s.aux_id = 0;
        std::string seq(len, 'A');
        s.owned_seqs.emplace_back(seq.c_str(), 0, len, 0,
                                  SequenceContainerType::from_ascii{});
        sn.push_back(std::move(s));
    };
    auto add_edge = [](edge_map &em, u64 s, u64 d) {
        edge_key k{s, d, '+', '+', uint8_t{0}, uint8_t{0}};
        edge_agg a;
        a.src_super = s;
        a.dst_super = d;
        a.count = 1;
        em[k] = a;
    };

    // Topology (M = layout, C = unitig; short len 4, long len 50):
    //   M0 → C1(short, tip)                 → removed
    //   M0 → C2(short) → C3(short, dead-end) → both peeled (C3 then C2)
    //   M0 → C4(long, tip)                  → kept (too long)
    supernode_vec sn;
    add_node(sn, "M0", true, 100);
    add_node(sn, "C1", false, 4);
    add_node(sn, "C2", false, 4);
    add_node(sn, "C3", false, 4);
    add_node(sn, "C4", false, 50);

    edge_map edges;
    add_edge(edges, 0, 1);
    add_edge(edges, 0, 2);
    add_edge(edges, 2, 3);
    add_edge(edges, 0, 4);

    PostLayoutCleaning::run(sn, edges, /*max_tip_len*/ 10);

    CHECK_FALSE(sn[0].name.empty());  // M-node never removed
    CHECK(sn[1].name.empty());        // short tip removed
    CHECK(sn[2].name.empty());        // peeled after C3
    CHECK(sn[3].name.empty());        // short dead-end removed
    CHECK_FALSE(sn[4].name.empty());  // long tip kept

    // Only the M0 → C4 edge survives.
    REQUIRE_EQ(edges.size(), 1u);
    CHECK_EQ(edges.begin()->second.src_super, 0u);
    CHECK_EQ(edges.begin()->second.dst_super, 4u);
}

TEST_CASE(
    "layout_graph_emit::inject_thread_bridge_edges adds 1:1 scaffold edges") {
    using namespace layout_graph_emit;
    // Two layout supernodes; DBG node 100 sits at M0's END boundary, node 200 at
    // M1's START boundary. Bridge 100->200 must add exactly one TP:Z:bridge
    // superedge M0(end)->M1(start) carrying the support.
    supernode_vec sn(2);
    sn[0].name = "M0"; sn[0].is_layout = true; sn[0].aux_id = 0;
    sn[1].name = "M1"; sn[1].is_layout = true; sn[1].aux_id = 0;
    gtl::flat_hash_map<u64, u64> n2s{{100, 0}, {200, 1}};
    gtl::flat_hash_set<u64> in_layouts{100, 200};
    gtl::flat_hash_map<u64, bool> ori{{100, true}, {200, true}};
    // ori=true => src side END, dst side START; masks must include those sides.
    gtl::flat_hash_map<u64, uint8_t> mask{{100, SIDE_END}, {200, SIDE_START}};
    gtl::vector<bridge_edge> bridges{{100, 200, 7}};

    edge_map edges;
    CHECK_EQ(inject_thread_bridge_edges(bridges, sn, n2s, in_layouts, ori, mask,
                                        edges).added,
             1u);
    REQUIRE_EQ(edges.size(), 1u);
    [[maybe_unused]] const auto &e = edges.begin()->second;
    CHECK(e.thread_bridge);
    CHECK_EQ(e.tb_support, 7u);
    CHECK_EQ(e.src_super, 0u);
    CHECK_EQ(e.dst_super, 1u);
    CHECK_EQ(e.src_side, SIDE_END);
    CHECK_EQ(e.dst_side, SIDE_START);
    CHECK_FALSE(thread_bridge_is_fasta_safe(e));

    // Anchor not on the needed boundary side -> dropped (no valid contig end).
    gtl::flat_hash_map<u64, uint8_t> mask_mid{{100, SIDE_START},
                                              {200, SIDE_START}};
    edge_map edges2;
    CHECK_EQ(inject_thread_bridge_edges(bridges, sn, n2s, in_layouts, ori,
                                        mask_mid, edges2).added,
             0u);

    // Ambiguous: two bridges share the M0-END endpoint -> neither injected.
    supernode_vec sn3(3);
    sn3[0].name = "M0"; sn3[0].is_layout = true; sn3[0].aux_id = 0;
    sn3[1].name = "M1"; sn3[1].is_layout = true; sn3[1].aux_id = 0;
    sn3[2].name = "M2"; sn3[2].is_layout = true; sn3[2].aux_id = 0;
    gtl::flat_hash_map<u64, u64> n2s3{{100, 0}, {200, 1}, {300, 2}};
    gtl::flat_hash_set<u64> inl3{100, 200, 300};
    gtl::flat_hash_map<u64, bool> ori3{{100, true}, {200, true}, {300, true}};
    gtl::flat_hash_map<u64, uint8_t> mask3{
        {100, SIDE_END}, {200, SIDE_START}, {300, SIDE_START}};
    gtl::vector<bridge_edge> bridges3{{100, 200, 7}, {100, 300, 7}};
    edge_map edges3;
    CHECK_EQ(inject_thread_bridge_edges(bridges3, sn3, n2s3, inl3, ori3, mask3,
                                        edges3).added,
             0u);

    // A layout edge already connecting the pair -> bridge not duplicated.
    edge_map edges4;
    edge_key existing{0, 1, '+', '+', SIDE_END, SIDE_START};
    edges4[existing].src_super = 0;
    edges4[existing].dst_super = 1;
    edges4[existing].count = 1;
    CHECK_EQ(inject_thread_bridge_edges(bridges, sn, n2s, in_layouts, ori, mask,
                                        edges4).added,
             0u);
}

TEST_CASE("thread-bridge branch queries read vectors through locked callbacks") {
    dbg_index_light dbgi(0.001, 16, 42);
    constexpr u64 in_branch = 10;
    constexpr u64 out_branch = 20;
    constexpr u64 linear = 30;
    dbgi.dbg.emplace(in_branch,
                     dbg_index_light::bidbg_val{1u, 1u});
    dbgi.dbg.emplace(out_branch, dbg_index_light::bidbg_val{1u});
    dbgi.dbg.emplace(linear, dbg_index_light::bidbg_val{1u, 1u});
    dbgi.succ_track.emplace(out_branch, gtl::vector<u64>{2u, 3u});
    dbgi.succ_track.emplace(linear, gtl::vector<u64>{2u});
    dbgi.dbg.modify_if(in_branch, [](auto &entry) {
        entry.second.emplace_prev(2u);
    });

    CHECK(assembler::tb_is_branch(dbgi, in_branch));
    CHECK(assembler::tb_is_branch(dbgi, out_branch));
    CHECK_FALSE(assembler::tb_is_branch(dbgi, linear));
    CHECK_FALSE(assembler::tb_is_branch(dbgi, 999u));
}

TEST_CASE("thread-bridge observations deduplicate read votes and summarize signed gaps") {
    using namespace layout_graph_emit;
    retained_read read{
        /*qname_hash=*/42u,
        /*nodes=*/{1u, 2u, 0u, 3u},
        /*query_pos=*/{0u, 20u, 30u, 90u, 130u}};
    gtl::vector<thread_bridge_read_observation> coordinate_observations;
    append_thread_bridge_observations(
        read,
        [](u64 raw) {
            if (raw == 1 || raw == 2) return u64{10};
            return raw == 3 ? u64{20} : u64{0};
        },
        coordinate_observations);
    REQUIRE_EQ(coordinate_observations.size(), 1u);
    CHECK_EQ(coordinate_observations[0].gap, 60);

    constexpr u64 read1 = 1u, read2 = 2u;
    gtl::vector<thread_bridge_read_observation> raw{
        {read2, 10, 20, 90}, {read1, 10, 20, 50},
        {read1, 10, 20, -10}, {read1, 10, 20, 30},
        {read1, 20, 10, 7}};
    auto distinct = deduplicate_thread_bridge_observations(std::move(raw));
    REQUIRE_EQ(distinct.size(), 3u);
    auto vote = std::find_if(distinct.begin(), distinct.end(), [](const auto &o) {
        return o.qname == read1 && o.u == 10 && o.v == 20;
    });
    REQUIRE(vote != distinct.end());
    CHECK_EQ(vote->gap, 30);

    gtl::vector<int64_t> gaps{-20, 100, 110, 120, 130};
    [[maybe_unused]] const signed_gap_summary summary =
        summarize_signed_gaps(gaps);
    CHECK_EQ(summary.count, 5u);
    CHECK_EQ(summary.median, 110);
    CHECK_EQ(summary.p10, -20);
    CHECK_EQ(summary.p90, 120);
    CHECK_EQ(summary.mad, 10);
    CHECK_EQ(classify_bridge_realization(summary), bridge_realization::gap_safe);

    gtl::vector<int64_t> conflicting_tail{
        127, 128, 128, 128, 128, 129, 129, 129, 129, 1392, 1463, 2423, 4961};
    CHECK_EQ(classify_bridge_realization(summarize_signed_gaps(conflicting_tail)),
             bridge_realization::gap_uncertain);

    gtl::vector<int64_t> overlap{-50, -20, -10, 20, 30};
    CHECK_EQ(classify_bridge_realization(summarize_signed_gaps(overlap)),
             bridge_realization::overlap_unresolved);
    gtl::vector<int64_t> insufficient{100, 110, 120, 130};
    CHECK_EQ(classify_bridge_realization(summarize_signed_gaps(insufficient)),
             bridge_realization::gap_uncertain);

    edge_agg edge;
    edge.thread_bridge = true;
    edge.tb_gap = 110;
    edge.tb_realization = bridge_realization::gap_safe;
    CHECK(thread_bridge_is_fasta_safe(edge));
    edge.tb_realization = bridge_realization::gap_uncertain;
    CHECK_FALSE(thread_bridge_is_fasta_safe(edge));
    edge.tb_realization = bridge_realization::overlap_unresolved;
    CHECK_FALSE(thread_bridge_is_fasta_safe(edge));
}

TEST_CASE("thread-bridge injector aggregates reciprocal projections and rejects conflicts") {
    using namespace layout_graph_emit;
    supernode_vec sn(5);
    for (size_t i = 0; i < sn.size(); ++i) {
        sn[i].name = (i < 3 ? "M" : "C") + std::to_string(i);
        sn[i].is_layout = i < 3;
        sn[i].aux_id = i;
    }
    // 100(+ source) and 101(- destination) both project to M0-END;
    // 200(+ destination) and 201(- source) both project to M1-START.
    gtl::flat_hash_map<u64, u64> n2s{{100, 0}, {101, 0}, {200, 1},
                                     {201, 1}, {300, 2}, {400, 3},
                                     {500, 4}};
    gtl::flat_hash_set<u64> in_layouts{100, 101, 200, 201, 300};
    gtl::flat_hash_map<u64, bool> ori{{100, true},  {101, false},
                                      {200, true},  {201, false},
                                      {300, true}};
    gtl::flat_hash_map<u64, uint8_t> mask{{100, SIDE_END}, {101, SIDE_END},
                                          {200, SIDE_START}, {201, SIDE_START},
                                          {300, SIDE_END}};
    auto with_reads = [](u64 u, u64 v,
                         std::initializer_list<bridge_read_gap_observation>
                             observations) {
        bridge_edge bridge{u, v, static_cast<uint32_t>(observations.size())};
        bridge.read_gap_observations.assign(observations.begin(),
                                            observations.end());
        return bridge;
    };

    SUBCASE("equivalent reciprocal records aggregate") {
        gtl::vector<bridge_edge> bridges{{100, 200, 7, 100},
                                         {201, 101, 5, 120}};
        edge_map edges;
        [[maybe_unused]] const auto result = inject_thread_bridge_edges(
            bridges, sn, n2s, in_layouts, ori, mask, edges);
        CHECK_EQ(result.added, 1u);
        CHECK_EQ(result.aggregated_records, 1u);
        REQUIRE_EQ(edges.size(), 1u);
        CHECK_EQ(edges.begin()->second.tb_support, 12u);
    }

    SUBCASE("projected support deduplicates qname across anchor pairs") {
        gtl::vector<bridge_edge> bridges;
        constexpr u64 read_shared = 77u;
        bridges.push_back(
            with_reads(100, 200, {{read_shared, 100}}));
        bridges.push_back(
            with_reads(201, 101, {{read_shared, 120}}));
        edge_map edges;
        [[maybe_unused]] const auto result = inject_thread_bridge_edges(
            bridges, sn, n2s, in_layouts, ori, mask, edges);
        CHECK_EQ(result.added, 1u);
        REQUIRE_EQ(edges.size(), 1u);
        [[maybe_unused]] const auto &edge = edges.begin()->second;
        CHECK_EQ(edge.tb_support, 1u);
        CHECK_EQ(edge.count, 1u);
        CHECK_EQ(edge.tb_gap_summary.count, 1u);
        // The duplicate qname contributes the lower median min(100, 120).
        CHECK_EQ(edge.tb_gap, 100);
        CHECK_EQ(edge.tb_realization, bridge_realization::gap_uncertain);
    }

    SUBCASE("reciprocal M-C evidence treats C orientation as a wildcard") {
        gtl::vector<bridge_edge> bridges{{100, 400, 5, 100},
                                         {400, 101, 5, 120}};
        edge_map edges;
        [[maybe_unused]] const auto result = inject_thread_bridge_edges(
            bridges, sn, n2s, in_layouts, ori, mask, edges);
        CHECK_EQ(result.added, 1u);
        CHECK_EQ(result.aggregated_records, 1u);
        REQUIRE_EQ(edges.size(), 1u);
        CHECK_EQ(edges.begin()->second.tb_support, 10u);
        CHECK(edges.begin()->second.tb_c_orientation_unresolved);
    }

    SUBCASE("reciprocal C-C evidence aggregates with both orientations unknown") {
        gtl::vector<bridge_edge> bridges{{400, 500, 5, 100},
                                         {500, 400, 5, 120}};
        edge_map edges;
        [[maybe_unused]] const auto result = inject_thread_bridge_edges(
            bridges, sn, n2s, in_layouts, ori, mask, edges);
        CHECK_EQ(result.added, 1u);
        CHECK_EQ(result.aggregated_records, 1u);
        REQUIRE_EQ(edges.size(), 1u);
        CHECK_EQ(edges.begin()->second.tb_support, 10u);
        CHECK(edges.begin()->second.tb_c_orientation_unresolved);
        CHECK_FALSE(thread_bridge_is_fasta_safe(edges.begin()->second));
    }

    SUBCASE("reciprocal side conflict rejects the join") {
        gtl::vector<bridge_edge> bridges{{100, 200, 7, 100},
                                         {200, 100, 7, 100}};
        // Both masks are widened so the conflicting START/END projection gets
        // through endpoint resolution and is rejected as a projection conflict.
        auto both_sides = mask;
        both_sides[100] = SIDE_START | SIDE_END;
        both_sides[200] = SIDE_START | SIDE_END;
        edge_map edges;
        [[maybe_unused]] const auto result = inject_thread_bridge_edges(
            bridges, sn, n2s, in_layouts, ori, both_sides, edges);
        CHECK_EQ(result.added, 0u);
        CHECK_EQ(result.projection_conflict, 1u);
    }

    SUBCASE("one physical side reused as source and destination rejects both") {
        // 100->200 consumes M0-END as source. 300->101 consumes M0-END as
        // destination; direction is deliberately swapped relative to its reuse.
        gtl::vector<bridge_edge> bridges{{100, 200, 7, 100},
                                         {300, 101, 7, 100}};
        edge_map edges;
        [[maybe_unused]] const auto result = inject_thread_bridge_edges(
            bridges, sn, n2s, in_layouts, ori, mask, edges);
        CHECK_EQ(result.added, 0u);
        CHECK_EQ(result.endpoint_reuse, 2u);
        CHECK(edges.empty());
    }

    SUBCASE("standalone C orientation keeps GFA edge but blocks FASTA") {
        gtl::vector<bridge_edge> bridges{{400, 200, 7, 100}};
        edge_map edges;
        [[maybe_unused]] const auto result = inject_thread_bridge_edges(
            bridges, sn, n2s, in_layouts, ori, mask, edges);
        CHECK_EQ(result.added, 1u);
        REQUIRE_EQ(edges.size(), 1u);
        CHECK_EQ(edges.begin()->second.tb_realization,
                 bridge_realization::gap_safe);
        CHECK(edges.begin()->second.tb_c_orientation_unresolved);
        CHECK_FALSE(thread_bridge_is_fasta_safe(edges.begin()->second));
    }

    SUBCASE("unknown C side does not create a false physical-end conflict") {
        gtl::vector<bridge_edge> bridges{{100, 400, 7, 100},
                                         {400, 200, 7, 100}};
        edge_map edges;
        [[maybe_unused]] const auto result = inject_thread_bridge_edges(
            bridges, sn, n2s, in_layouts, ori, mask, edges);
        CHECK_EQ(result.added, 2u);
        CHECK_EQ(result.endpoint_reuse, 0u);
        REQUIRE_EQ(edges.size(), 2u);
        for (const auto &entry : edges) {
            [[maybe_unused]] const auto &edge = entry.second;
            CHECK(edge.tb_c_orientation_unresolved);
            CHECK_FALSE(thread_bridge_is_fasta_safe(edge));
        }
    }

    SUBCASE("equivalent reciprocal existing edge is not duplicated") {
        gtl::vector<bridge_edge> bridges{{100, 200, 7, 100}};
        edge_map edges;
        edge_agg existing;
        existing.src_super = 1;
        existing.dst_super = 0;
        existing.src_ori = '-';
        existing.dst_ori = '-';
        existing.src_side = SIDE_START;
        existing.dst_side = SIDE_END;
        edges[{1, 0, '-', '-', SIDE_START, SIDE_END}] = existing;
        [[maybe_unused]] const auto result = inject_thread_bridge_edges(
            bridges, sn, n2s, in_layouts, ori, mask, edges);
        CHECK_EQ(result.added, 0u);
        CHECK_EQ(result.existing_edge, 1u);
        CHECK_EQ(edges.size(), 1u);
    }
}

TEST_CASE("thread-bridge FASTA requires its paths sidecar") {
    using namespace layout_graph_emit;
    const auto unique = std::to_string(
        std::chrono::steady_clock::now().time_since_epoch().count());
    const std::filesystem::path temp_root =
        std::filesystem::temp_directory_path() /
        ("lasm-thread-bridge-fasta-" + unique);
    const std::filesystem::path fasta_path = temp_root / "bridges.fa";
    std::error_code error;
    REQUIRE(std::filesystem::create_directories(
        std::filesystem::path{fasta_path.string() + ".paths.tsv"}, error));
    REQUIRE_FALSE(error);

    const auto result = write_thread_bridge_fasta(
        supernode_vec{}, edge_map{}, fasta_path.string());
    REQUIRE_FALSE(result.has_value());
    REQUIRE(result.error().find("path sidecar") != std::string::npos);

    std::filesystem::remove_all(temp_root, error);
    CHECK_FALSE(error);
}

TEST_CASE("layout_graph_emit::find_heaviest_paths emits only the connected "
          "heaviest chain, never a disconnected prefix (R1)") {
    using namespace layout_graph_emit;

    // Main heaviest path A(0)→B(1)→C(2). An earlier source D(3) branches to a
    // heavy dead-end E(4) and also to A. Under the old span-extension heuristics
    // (C13) a greedy backward walk from D could try to splice D/E onto the front;
    // R1 removed all extensions, so the result is purely the virtual-source→sink
    // heaviest path — the connected chain [0,1,2]. D(3) and E(4) are never
    // prepended because they are not on that traceback.
    gtl::vector<lane_node> nodes(5);
    nodes[0] = {0, 10, 20, true, nullptr};  // A
    nodes[1] = {1, 20, 30, true, nullptr};  // B
    nodes[2] = {2, 30, 40, true, nullptr};  // C
    nodes[3] = {3, 0, 5, true, nullptr};    // D (earliest source)
    nodes[4] = {4, 1, 6, true, nullptr};    // E

    gtl::vector<scd_node> scd(5);
    for (uint32_t i = 0; i < 5; ++i) {
        scd[i].contig_idx = i;
        scd[i].core_idx = 0;
        scd[i].core_label = i;
    }
    scd[0].weight = 1;  // A
    scd[1].weight = 1;  // B
    scd[2].weight = 1;  // C
    scd[3].weight = 0;  // D (0 → not pulled onto the main path as a parent)
    scd[4].weight = 2;  // E (heavier → greedy walk from D goes here)

    gtl::vector<scd_adj_t> scd_fwd(5);
    scd_fwd[0].push_back({1u, 1});  // A→B
    scd_fwd[1].push_back({2u, 1});  // B→C
    scd_fwd[3].push_back({0u, 1});  // D→A
    scd_fwd[3].push_back({4u, 1});  // D→E (dead-end)

    auto paths = find_heaviest_paths(scd, scd_fwd, nodes, /*layout_idx*/ 0);

    // The result must be exactly the connected chain [0,1,2] — D(3),E(4) are
    // NOT prepended, since with extensions removed the emitted path is just the
    // heaviest edge-connected traceback (0→1 and 1→2 are the only edges here).
    REQUIRE_EQ(paths.size(), 1u);
    REQUIRE_EQ(paths[0].path.size(), 3u);
    CHECK_EQ(paths[0].path[0], 0u);
    CHECK_EQ(paths[0].path[1], 1u);
    CHECK_EQ(paths[0].path[2], 2u);
}

TEST_CASE("layout_graph_emit::find_heaviest_paths does not forward-extend to a "
          "wider-span lighter sink (R1)") {
    using namespace layout_graph_emit;

    // Heaviest chain A(0)→B(1)→C(2) (C is heavy). B also branches to D(3), a
    // sink whose `to` (50) is wider than the whole main path (max_to 30) but
    // whose path weight is lower. The OLD forward extension spliced D on after
    // B — emitting [0,1,2,3] to grab the extra span. R1 removed extensions, so
    // the result is exactly the heaviest edge-connected chain [0,1,2]; the
    // wider-but-lighter D is dropped. This is the behavior R1 deliberately
    // changes (verify with QUAST that dropping such span does not regress).
    gtl::vector<lane_node> nodes(4);
    nodes[0] = {0, 0, 10, true, nullptr};   // A
    nodes[1] = {1, 10, 20, true, nullptr};  // B
    nodes[2] = {2, 20, 30, true, nullptr};  // C (main sink)
    nodes[3] = {3, 20, 50, true, nullptr};  // D (wider `to`, lighter)

    gtl::vector<scd_node> scd(4);
    for (uint32_t i = 0; i < 4; ++i) {
        scd[i].contig_idx = i;
        scd[i].core_idx = 0;
        scd[i].core_label = i;
    }
    scd[0].weight = 1;  // A
    scd[1].weight = 1;  // B
    scd[2].weight = 2;  // C → dp[C]=4 beats dp[D]=3
    scd[3].weight = 1;  // D

    gtl::vector<scd_adj_t> scd_fwd(4);
    scd_fwd[0].push_back({1u, 1});  // A→B
    scd_fwd[1].push_back({2u, 1});  // B→C
    scd_fwd[1].push_back({3u, 1});  // B→D (wider-span lighter branch)

    auto paths = find_heaviest_paths(scd, scd_fwd, nodes, /*layout_idx*/ 0);

    REQUIRE_EQ(paths.size(), 1u);
    REQUIRE_EQ(paths[0].path.size(), 3u);   // NOT 4 — D is not appended
    CHECK_EQ(paths[0].path[0], 0u);
    CHECK_EQ(paths[0].path[1], 1u);
    CHECK_EQ(paths[0].path[2], 2u);
    CHECK_EQ(paths[0].total_weight, 4);
}

TEST_CASE("layout_graph_emit::resolve_layout_orientations policies (C10)") {
    using namespace layout_graph_emit;
    using layout_graph_algo::layout_result_t;

    auto mk = [](i64 a, i64 b, bool o) { return std::make_tuple(a, b, o); };

    SUBCASE("no conflict → both policies agree, conflicts=0") {
        layout_result_t layouts;
        layouts.emplace_back();
        layouts.back().emplace(100u, mk(0, 10, true));
        layouts.back().emplace(200u, mk(0, 10, false));

        size_t conf = 123;
        auto fw = resolve_layout_orientations(
            layouts, orientation_policy::first_wins, &conf);
        CHECK_EQ(conf, 0u);
        CHECK_EQ(fw.at(100u), true);
        CHECK_EQ(fw.at(200u), false);

        auto mv = resolve_layout_orientations(
            layouts, orientation_policy::majority_vote, &conf);
        CHECK_EQ(conf, 0u);
        CHECK_EQ(mv.at(100u), true);
        CHECK_EQ(mv.at(200u), false);
    }

    SUBCASE("conflict → first_wins and majority_vote differ") {
        // node 100: layout0 says '+', layouts 1&2 say '-' (majority '-').
        // node 200: only layout0, '+', no conflict.
        layout_result_t layouts;
        layouts.emplace_back();
        layouts.back().emplace(100u, mk(0, 10, true));
        layouts.back().emplace(200u, mk(0, 10, true));
        layouts.emplace_back();
        layouts.back().emplace(100u, mk(0, 10, false));
        layouts.emplace_back();
        layouts.back().emplace(100u, mk(0, 10, false));

        size_t conf = 0;
        auto fw = resolve_layout_orientations(
            layouts, orientation_policy::first_wins, &conf);
        CHECK_EQ(conf, 1u);            // exactly node 100 conflicts
        CHECK_EQ(fw.at(100u), true);   // layout 0 claimed it first
        CHECK_EQ(fw.at(200u), true);

        auto mv = resolve_layout_orientations(
            layouts, orientation_policy::majority_vote, &conf);
        CHECK_EQ(conf, 1u);
        CHECK_EQ(mv.at(100u), false);  // 2 '-' outvote 1 '+'
        CHECK_EQ(mv.at(200u), true);
    }
}

TEST_CASE("layout_graph_emit::find_contained_cnodes keeps standalone/circular "
          "components (C3)") {
    using namespace layout_graph_emit;

    // Component zoo (expected outcome):
    //   10            isolated node (no edges at all)              -> keep
    //   20 <-> 21     closed circular pair (neighbors internal)    -> keep
    //   30 -> L1      adjacent to layout interior only             -> contained
    //   40 -> L2      adjacent to a layout boundary                -> keep
    //   50 -> 99      adjacent to a non-layout node                -> keep
    dbg_index dbgi(0.001, 1024, 42);
    auto add_node = [&](u64 id, std::initializer_list<u64> prevs) {
        dbgi.dbg.emplace(id, dbg_index::bidbg_val{100u});
        auto it = dbgi.dbg.find(id);
        for (u64 p : prevs) it->second.prev.push_back(p);
        it->second.seq.emplace_back(64, 0);  // non-empty seq => valid candidate
        dbgi.active_set.insert(id);
    };

    constexpr u64 L1 = 1001, L2 = 1002;  // layout nodes
    add_node(10, {});
    add_node(20, {21});
    add_node(21, {20});
    add_node(30, {L1});
    add_node(40, {L2});
    add_node(50, {99});
    add_node(L1, {});
    add_node(L2, {});
    // Node 99: present in the DBG but neither active nor part of a layout.
    dbgi.dbg.emplace(99, dbg_index::bidbg_val{100u});

    gtl::flat_hash_set<u64> nodes_in_layouts{L1, L2};
    gtl::flat_hash_map<u64, uint8_t> node_boundary_mask;
    gtl::flat_hash_map<u64, bool> node_layout_ori;
    node_boundary_mask[L1] = 0;          // interior: no side bits set
    node_layout_ori[L1] = false;
    node_boundary_mask[L2] = SIDE_START;  // boundary on the side a
    node_layout_ori[L2] = false;          // predecessor looks for

    gtl::flat_hash_set<u64> candidate_set;
    auto contained =
        find_contained_cnodes(dbgi, nodes_in_layouts, node_boundary_mask,
                              node_layout_ori, candidate_set);

    CHECK(candidate_set.size() == 6);  // all active non-layout nodes
    CHECK(contained.size() == 1);
    CHECK(contained.contains(30));   // only the interior-adjacent component
    CHECK(!contained.contains(10));  // isolated: kept (was dropped pre-C3-fix)
    CHECK(!contained.contains(20));  // closed circular pair: kept (was dropped)
    CHECK(!contained.contains(21));
    CHECK(!contained.contains(40));  // boundary-adjacent: kept
    CHECK(!contained.contains(50));  // non-layout-adjacent: kept
}

TEST_CASE("GraphCleaning::simplify_simple_paths terminates on self-loops (C1b)"
          * doctest::timeout(10)) {
    // A self-loop node with fec==1 used to send the backward walk into an
    // infinite merge-itself loop (worker spun inside
    // attempt_simplify_one_path_lockless growing `merged`; main parked in
    // join). The walk must now skip self steps, and the post-sweep redirect
    // rewrite must drop self-references.
    dbg_index dbgi(0.001, 1024, 7);
    auto add_node = [&](u64 id, std::initializer_list<u64> prevs,
                        unsigned fec) {
        dbgi.dbg.emplace(id, dbg_index::bidbg_val{100u});
        auto it = dbgi.dbg.find(id);
        for (u64 p : prevs) it->second.prev.push_back(p);
        it->second.forward_edge_count = fec;
        dbgi.active_set.insert(id);
    };
    add_node(100, {100}, 1);       // A: self-loop, fec==1
    add_node(200, {100, 300}, 0);  // A0: multi-prev (type-2 seeds A)
    add_node(300, {}, 1);          // B

    GraphCleaning::redir_map redirect_deleted;
    GraphCleaning::simplify_simple_paths(dbgi.dbg, redirect_deleted,
                                         dbgi.active_set, 4);

    auto a = dbgi.dbg.find(100);
    REQUIRE(a != dbgi.dbg.end());
    for ([[maybe_unused]] u64 p : a->second.prev) {
        CHECK(p != 100);  // no self-reference survives the rewrite
    }
    CHECK(a->second.prev.empty());  // its only edge was the self-edge
    CHECK(a->second.forward_edge_count == 0);  // bookkeeping stays consistent

    [[maybe_unused]] auto a0 = dbgi.dbg.find(200);
    REQUIRE(a0 != dbgi.dbg.end());
    CHECK(a0->second.prev.size() == 2);  // untouched
    CHECK(dbgi.dbg.find(300) != dbgi.dbg.end());
}

TEST_CASE("GraphCleaning::simplify_simple_paths leaves degenerate boundary "
          "seeds untouched (C1)") {
    // Configuration: boundary node S (fec=2, prev={T}, T.fec=2) is a type-1
    // seed whose backward walk is degenerate (merges nothing), while a
    // reader chain stops at S. The degenerate path must be a pure no-op:
    // mutating S would race the reader chain's stop-condition read of S.fec
    // (see CONFIRMATION_C1.md). Under a TSan build this doubles as the race
    // regression test (pre-fix it reports graph.h read/write on S.fec).
    dbg_index dbgi(0.001, 1024, 42);
    auto add_node = [&](u64 id, std::initializer_list<u64> prevs,
                        unsigned fec) {
        dbgi.dbg.emplace(id, dbg_index::bidbg_val{100u});
        auto it = dbgi.dbg.find(id);
        for (u64 p : prevs) it->second.prev.push_back(p);
        it->second.forward_edge_count = fec;
        dbgi.active_set.insert(id);
    };
    constexpr u64 T = 1, D = 2, S = 3;
    add_node(T, {}, 2);
    add_node(D, {T}, 0);
    add_node(S, {T}, 2);
    // Reader chain: X(=109) -> 103 -> 102 -> 101 -> S (backward walk order).
    add_node(101, {S}, 1);
    add_node(102, {101}, 1);
    add_node(103, {102}, 1);
    add_node(109, {103}, 2);
    add_node(110, {109}, 0);
    add_node(111, {109}, 0);

    GraphCleaning::redir_map redirect_deleted;
    GraphCleaning::simplify_simple_paths(dbgi.dbg, redirect_deleted,
                                         dbgi.active_set, 4);

    // Degenerate boundary seed S: completely untouched.
    [[maybe_unused]] auto s = dbgi.dbg.find(S);
    REQUIRE(s != dbgi.dbg.end());
    CHECK(s->second.forward_edge_count == 2);
    CHECK(s->second.counts.core_count == 100);
    CHECK(s->second.prev.size() == 1);
    CHECK(s->second.prev[0] == T);
    CHECK(!s->second.redirected());

    // Reader chain merged into its representative (node 101, next to S).
    [[maybe_unused]] auto rep = dbgi.dbg.find(101);
    REQUIRE(rep != dbgi.dbg.end());
    CHECK(rep->second.forward_edge_count == 2);      // inherited seed's fec
    CHECK(rep->second.counts.core_count == 400);     // 100 + 3 merged nodes
    for (u64 merged : {102, 103, 109}) {
        [[maybe_unused]] auto m = dbgi.dbg.find(merged);
        REQUIRE(m != dbgi.dbg.end());
        CHECK(m->second.redirected());
    }
    CHECK(redirect_deleted.find(109)->second == 101);
}

// ── Fix 9: MinCode syncmer rule arithmetic ──────────────────────────────
// Uses the Identity hash so the "minimum-hash" submer is just the minimum
// integer submer, making expected results hand-computable. k=4, s=2 ⇒ three
// 2-base (4-bit, mask 0xF) windows at bit offsets 0,2,4. For a k-mer `val`:
//   w0 = val & 0xF,  w1 = (val>>2) & 0xF,  w2 = (val>>4) & 0xF
// min_index = argmin(w0,w1,w2) (first on ties). Closed ⇒ min_index in {0, k-s=2};
// Open(t) ⇒ min_index == t.
TEST_CASE("Syncmer rule (closed/open): s-mer window arithmetic") {
    using Rule = syncmer::GenericSubmerRule<uint32_t, syncmer::Identity>;

    // val=0x1B=0b00011011: w0=11, w1=6, w2=1 ⇒ min at index 2
    // val=0x8F=0b10001111: w0=15, w1=3, w2=8 ⇒ min at index 1
    // val=0xF0=0b11110000: w0=0,  w1=12, w2=15 ⇒ min at index 0
    SUBCASE("closed") {
        decltype(asm_config::syncmer) p;
        p.k = 4;
        p.smer.s = 2;  // closed: no t
        Rule sr{p, 0};
        CHECK(sr(0x1B));        // min index 2 == k-s ⇒ closed syncmer
        CHECK_FALSE(sr(0x8F));  // min index 1 (interior) ⇒ not
        CHECK(sr(0xF0));        // min index 0 ⇒ closed syncmer
    }
    SUBCASE("open t=1") {
        decltype(asm_config::syncmer) p;
        p.k = 4;
        p.smer.s = 2;
        p.smer.t = 1;
        Rule sr{p, 0};
        CHECK(sr(0x8F));        // min index 1 == t
        CHECK_FALSE(sr(0x1B));  // min index 2 != 1
        CHECK_FALSE(sr(0xF0));  // min index 0 != 1
    }
    SUBCASE("open t=0 matches window 0 (first base)") {
        decltype(asm_config::syncmer) p;
        p.k = 4;
        p.smer.s = 2;
        p.smer.t = 0;
        Rule sr{p, 0};
        CHECK(sr(0xF0));        // min index 0 == t
        CHECK_FALSE(sr(0x8F));  // min index 1 != 0
    }
    // The last window (index k-s) must actually be examined: a k-mer whose
    // unique minimum submer sits in the final window is a closed syncmer.
    SUBCASE("last window is examined") {
        decltype(asm_config::syncmer) p;
        p.k = 4;
        p.smer.s = 2;
        Rule sr{p, 0};
        CHECK(sr(0x1B));  // min submer (=1) is window 2, the last one
    }
}

// --canonical: applying the rule to min(kmer, reverse-complement) must make a
// k-mer and its reverse complement always agree, and the flag must actually
// change behavior vs the default (directional) rule.
TEST_CASE("Strand-canonical selection co-selects a k-mer and its rev-comp") {
    using Rule = syncmer::GenericSubmerRule<uint32_t, syncmer::MurMur>;
    const int k = 8;
    decltype(asm_config::syncmer) p;
    p.k = k;
    p.smer.s = 3;  // closed syncmer
    p.canonical = true;
    Rule canon{p, 42};
    p.canonical = false;
    Rule direct{p, 42};

    [[maybe_unused]] bool saw_directional_diff = false;
    for (uint32_t x = 0; x < (1u << (2 * k)); ++x) {
        uint32_t rc = reverse_complement_kmer<uint32_t>(x, k);
        CHECK(canon(x) == canon(rc));  // co-selected under --canonical
        if (direct(x) != direct(rc)) saw_directional_diff = true;
    }
    CHECK(saw_directional_diff);  // the flag genuinely changes selection
}

// ── Spaced-seed compaction: '0' positions are non-existent in selection ───
// Mask "11011" (k=5) removes field 2 (bits 4-5; the middle position), so
// selection runs on the dense 4-mer of kept fields:
//   compact(v) = f0 | f1<<2 | f3<<4 | f4<<6      (f_i = (v >> 2i) & 3)
// With Identity hash and s=2 there are k_eff-s+1 = 3 s-windows over the
// compacted value; closed ⇒ min at window 0 or 2.
TEST_CASE("Spaced-seed compaction: masked positions are non-existent") {
    using Rule = syncmer::GenericSubmerRule<uint32_t, syncmer::Identity>;
    decltype(asm_config::syncmer) p;
    p.k = 5;
    p.smer.s = 2;
    p.spaced_seed_mask = "11011";

    SUBCASE("compact() gathers kept fields in order") {
        Rule sr{p, 0};
        CHECK(sr.k_eff == 4);
        CHECK(sr.compact(0x7B) == 0x1B);   // f0=3 f1=2 f2=3 f3=1 f4=0
        CHECK(sr.compact(0x3FF) == 0xFF);  // all fields 3
    }
    SUBCASE("closed syncmer windows live in compacted space") {
        Rule sr{p, 0};
        CHECK(sr(0x7B));         // compact 0x1B: min submer at last window
        CHECK_FALSE(sr(0x20F));  // compact 0x8F: min at interior window 1
        CHECK(sr(0x3C0));        // compact 0xF0: min at window 0
    }
    SUBCASE("k-mers differing only at the masked position co-select") {
        Rule sr{p, 0};
        // Same kept fields as 0x7B, with the masked field 2 set to 0/1/2.
        CHECK(sr(0x4B) == sr(0x7B));
        CHECK(sr(0x5B) == sr(0x7B));
        CHECK(sr(0x6B) == sr(0x7B));
    }
    SUBCASE("open syncmer t indexes compacted windows") {
        decltype(asm_config::syncmer) po = p;
        po.smer.t = 1;
        Rule sr{po, 0};
        CHECK(sr(0x20F));        // compact 0x8F: min at window 1 == t
        CHECK_FALSE(sr(0x3C0));  // compact 0xF0: min at window 0 != t
        po.smer.t = 0;
        Rule sr0{po, 0};
        CHECK(sr0(0x3C0));
        CHECK_FALSE(sr0(0x20F));
    }
    SUBCASE("'2'/'3' keep their position but still bit-mask") {
        decltype(asm_config::syncmer) p2 = p;
        p2.spaced_seed_mask = "11121";  // field 2 keeps only its high bit
        Rule sr{p2, 0};
        CHECK(sr.k_eff == 5);              // nothing removed
        CHECK(sr.compact(0x1AB) == 0x1AB); // no '0' ⇒ compaction is identity
        // field 2 = 0b10 vs 0b11 both bit-mask to 0b10 ⇒ co-selected.
        CHECK(sr(0x61) == sr(0x71));
    }
}

// --canonical with a removing mask: a k-mer and its reverse complement must
// still always agree (canonicalization runs on the raw k-mer, compaction is
// baked into the LuT).
TEST_CASE("Spaced-seed compaction with --canonical co-selects rev-comps") {
    using Rule = syncmer::GenericSubmerRule<uint32_t, syncmer::MurMur>;
    const int k = 5;
    decltype(asm_config::syncmer) p;
    p.k = k;
    p.smer.s = 2;
    p.canonical = true;
    p.spaced_seed_mask = "11011";
    Rule canon{p, 42};
    for (uint32_t x = 0; x < (1u << (2 * k)); ++x) {
        [[maybe_unused]] uint32_t rc = reverse_complement_kmer<uint32_t>(x, k);
        CHECK(canon(x) == canon(rc));
    }
}

// Mincode submer with compaction: the density threshold applies to the hash
// of the compacted value. ScaleHash (v << 26) makes the threshold exact:
// density 0.5 ⇒ threshold 2^31 ⇒ selected iff compact(v) < 32 (no wrap for
// compact values < 64).
TEST_CASE("Spaced-seed compaction: mincode hashes the compacted value") {
    struct ScaleHash {
        unsigned seed;
        ScaleHash(unsigned s) : seed{s} {}
        unsigned operator()(unsigned v) const { return v << 26; }
    };
    using Rule = syncmer::GenericSubmerRule<uint32_t, ScaleHash>;
    decltype(asm_config::syncmer) p;
    p.k = 5;
    p.density = 0.5;
    p.spaced_seed_mask = "11011";
    Rule sr{p, 0};
    CHECK(sr(0x4F));        // compact 31 < 32
    CHECK(sr(0x7F));        // same compact (masked field 2) ⇒ co-selected
    CHECK_FALSE(sr(0x80));  // compact 32
    CHECK_FALSE(sr(0xB0));  // same compact 32
}

// Validation must reject s/t bounds against the EFFECTIVE k (k minus the
// removed positions), and reject masks that remove every position.
TEST_CASE("Syncmer validation uses effective k after compaction") {
    decltype(asm_config::syncmer) p;
    p.k = 5;
    p.smer.s = 2;
    p.spaced_seed_mask = "11011";  // k_eff = 4

    auto ek = p.effective_k_checked();
    REQUIRE(ek.has_value());
    CHECK(*ek == 4);

    SUBCASE("t bounded by effective-k - s + 1") {
        p.smer.t = 2;  // == k_eff - s ⇒ last valid window
        CHECK(p.validate().has_value());
        p.smer.t = 3;  // == k_eff - s + 1 ⇒ out of range
        CHECK_FALSE(p.validate().has_value());
    }
    SUBCASE("s bounded by effective k") {
        p.smer.s = 3;
        CHECK(p.validate().has_value());
        p.smer.s = 4;  // == k_eff
        CHECK_FALSE(p.validate().has_value());
    }
    SUBCASE("mask removing every position is rejected") {
        p.spaced_seed_mask = "00000";
        CHECK_FALSE(p.validate().has_value());
    }
    SUBCASE("short mask leaves leading positions kept") {
        p.spaced_seed_mask = "101";  // removes field 1 only; fields 3,4 kept
        auto ek2 = p.effective_k_checked();
        REQUIRE(ek2.has_value());
        CHECK(*ek2 == 4);
        CHECK(p.validate().has_value());  // s=2 < 4
    }
    SUBCASE("mask without '0' keeps full k") {
        p.spaced_seed_mask = "11121";
        auto ek3 = p.effective_k_checked();
        REQUIRE(ek3.has_value());
        CHECK(*ek3 == 5);
        p.smer.s = 4;
        CHECK(p.validate().has_value());  // s=4 < 5 (old check used k)
    }
}

// ── Fix 10: bloom projected_count widening (H13) ────────────────────────
// A human-scale projected count (n > INT32_MAX) previously narrowed to a
// negative int, so ComputeNumBits produced a negative-then-wrapped bit count
// and make_unique threw bad_alloc — the filter could not even be constructed.
// With int64_t n it constructs cleanly and sizes correctly. A high eps keeps
// the test allocation to a few MB while still exercising the >INT32_MAX path.
TEST_CASE("BlockedAtomicBloomFilter handles n > INT32_MAX without narrowing") {
    const int64_t n = 3'000'000'000LL;  // > INT32_MAX (2.147e9)
    bloom_ext::BlockedAtomicBloomFilter<> bf(n, 0.99f);  // eps→small m (~8 MB)
    CHECK(bf.n == n);            // full value retained, not narrowed
    CHECK(bf.m > 10'000'000ULL); // large positive bit count from the full n
    CHECK(bf.num_words > 0);
}

// The clamp keeps a pathological n from requesting a petabyte allocation or
// tripping the double->uint64 cast; ComputeNumBits saturates at the ceiling.
TEST_CASE("BlockedAtomicBloomFilter clamps an absurd projected count") {
    bloom_ext::BlockedAtomicBloomFilter<> bf(1, 0.99f);
    bf.n = 1'000'000'000'000'000LL;  // 1e15
    CHECK(bf.ComputeNumBits() == static_cast<uint64_t>(1e13));  // clamped
}
