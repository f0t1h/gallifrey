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
#include <gtl/btree.hpp>
#include <gtl/phmap.hpp>
#include <gtl/phmap_fwd_decl.hpp>
#include <gtl/vector.hpp>
#include <iostream>
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
// CMake knob bridge: LASM_OLC_SKIP_LAST_VALUE (plain value macro; function-like
// defines cannot cross the compiler command line) -> OLC_SKIP_LAST(l).
#ifdef LASM_OLC_SKIP_LAST_VALUE
#define OLC_SKIP_LAST(l) (LASM_OLC_SKIP_LAST_VALUE)
#endif
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
          parallel_file_reading{cfg.parallel_file_reading} {
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
                                  int lcmer_l, long skip_first = -1,
                                  long skip_last = -1) {
        // Cores skipped at each node end before overlap-anchoring (see the
        // OLC_SKIP_* macros). Node ends are boundary-noisy, so a small skip
        // avoids false anchors. Callers may override (e.g. 0/0 for the relaxed
        // reverse lookup of short partners).
        if (skip_first < 0) skip_first = OLC_SKIP_FIRST(lcmer_l);
        if (skip_last < 0) skip_last = OLC_SKIP_LAST(lcmer_l);
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

    static auto process_data(int thread_id, int lcp_level, dqt &data_queue,
                             dqt &recycle_queue, dbg_index &dbgi,
                             std::atomic<bool> &done, int &processed_read_count,
                             int &all_reads, auto Process_Func,
                             SyncmerType syncmer_functor, int lcmer_l) {
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
                Process_Func(lr, dbgi, thread_id, lcmer_index, lcmer_l);
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

        for (int tid = 1; tid <= threads.work; ++tid) {
            workers.emplace_back([&, tid]() {
                process_data(tid, lcp_level, data_queue, recycle_queue, dbgi,
                             done, valid_read_counts[tid - 1],
                             all_read_counts[tid - 1], process_cores<DefaultBF>,
                             syncmer_functor, dbg_lcmer_l);
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

        // Discard a match whose overlap spans less than this fraction (1/N) of
        // the SHORTER node's cores — too short to be a reliable overlap. The
        // gate is enforced at the reciprocity stage where both nodes' spans and
        // totals are known; the old query-relative form made containment
        // overlaps (short node inside long node) structurally impossible.
        constexpr size_t MIN_OVERLAP_CORE_FRACTION_DIV = 20;

        std::atomic<int> all_n{0};
        std::atomic<int> hits{0};
        std::atomic<int> bad{0};
        // Phase-0 matching diagnostics: why do good nodes fail to match?
        std::atomic<int> too_short{0};   // < 4*olc_lcmer_l cores → no ranges (#4)
        std::atomic<int> no_ranges{0};   // long enough but 0 lcmer hits
        std::atomic<int> frac_drop{0};   // reciprocal but < 5% of shorter node
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
                const int pre_deepen_size = str.size;
                lps_deepen(&str, lcp_level);

                for (int i = 0; i < str.size; ++i) {
                    lr.cores.emplace_back(str.cores[i].label);
                    lr.pos.emplace_back(str.cores[i].start, str.cores[i].end);
                }
                if (lr.cores.size() < static_cast<size_t>(olc_lcmer_l + 1)) {
                    if (getenv("LASM_DUMP_COREADMIT")) {
                        fprintf(stderr, "[coreadmit] %" PRIu64
                                "\tpre=%d\tpost=%zu\n",
                                node.first, pre_deepen_size, lr.cores.size());
                    }
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
        int non_recip = 0, recip_kept = 0, relaxed_kept = 0;
        // Absolute overlap-span floor in cores: half an lcmer window. Below
        // this, an "overlap" rests on less than one anchor's worth of
        // progression and is treated as noise (see the fraction gate).
        const u64 span_floor = std::max<u64>(2, (u64)olc_lcmer_l / 2);
        // rc-cores against the DBG with skips disabled, its total core count).
        // Filled lazily only for partners of one-directional matches.
        struct relaxed_info {
            gtl::vector<std::tuple<u64, u64, u64>> ranges;
            u64 total_cores = 0;
        };
        gtl::flat_hash_map<u64, relaxed_info> relaxed_cache;
        auto relaxed_ranges_for = [&](u64 nid) -> relaxed_info & {
            auto [it, inserted] = relaxed_cache.try_emplace(nid);
            if (!inserted) return it->second;
            auto nit = dbgi.dbg.find(nid);
            if (nit == dbgi.dbg.end() || nit->second.redirected() ||
                nit->second.seq.empty()) {
                return it->second;
            }
            lps rstr = {0, 0, 0, 0};
            lcp_read b_lr =
                build_rev_lr_cores(nit->second.seq, syncmer_functor, lcp_level,
                                   rstr);
            it->second.total_cores = static_cast<u64>(
                (b_lr.cores.size() > 0) ? b_lr.cores.size() - 1 : 0);
            if (b_lr.cores.size() >= static_cast<size_t>(olc_lcmer_l + 1)) {
                auto rr = find_lcmer_ranges(b_lr, dbgi, olc_lcmer_l, 0, 0);
                for (const auto &[rid, fr, to, fl, tl] : rr) {
                    if (rid != 0) it->second.ranges.emplace_back(rid, fr, to);
                }
            }
            free_lps(&rstr);
            return it->second;
        };
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
                    // Span-symmetric fraction gate: the overlap extent must be
                    // >= 1/MIN_OVERLAP_CORE_FRACTION_DIV of the SHORTER node's
                    // total cores, and >= span_floor cores in ABSOLUTE terms.
                    // The absolute floor matters at deeper LCP levels: cores
                    // are sparser, so 5% of a short node can be 1-2 cores —
                    // shorter than one anchor window and pure noise.
                    const u64 span_a = match.self.to - match.self.from;
                    const u64 span_b = ff->self.to - ff->self.from;
                    const u64 shorter_total =
                        std::min(match.self.len, ff->self.len);
                    if (std::min(span_a, span_b) <
                        std::max(shorter_total / MIN_OVERLAP_CORE_FRACTION_DIV,
                                 span_floor)) {
                        frac_drop.fetch_add(1);
                        continue;
                    }
                    // NOTE: a full-full rejection (both spans >= 95% of their
                    // nodes) was measured and REVERTED: it dropped 20 pairs,
                    // misassemblies 41->40 (noise), dup ratio worse
                    // (1.062->1.104). Near-identical SD copies are
                    // indistinguishable at the lcmer-key level; separating
                    // them needs read flank-association (future work).
                    reciprocal.emplace(k,
                                       Overlap{match.other, match.self, ff->self});
                    ++recip_kept;
                } else {
                    // Relaxed reverse lookup: A→B was stored but B→A was not.
                    // B may have been too short for the standard skip gates, or
                    // its end noise hid A. Re-run B's match with skips
                    // disabled; if A appears, keep the overlap with B's REAL
                    // span under the same symmetric fraction gate. (The old
                    // OLC_NONRECIPROCAL prototype faked B's span
                    // self-symmetrically instead.)
                    relaxed_info &rinfo = relaxed_ranges_for(match.other);
                    const auto best = std::max_element(
                        rinfo.ranges.begin(), rinfo.ranges.end(),
                        [](const auto &a, const auto &b) {
                            return std::get<2>(a) - std::get<1>(a) <
                                   std::get<2>(b) - std::get<1>(b);
                        });
                    if (best != rinfo.ranges.end() &&
                        std::get<0>(*best) == k) {
                        const u64 span_a = match.self.to - match.self.from;
                        const u64 span_b = std::get<2>(*best) - std::get<1>(*best);
                        const u64 shorter_total =
                            std::min(match.self.len, rinfo.total_cores);
                        if (std::min(span_a, span_b) >=
                            std::max(shorter_total / MIN_OVERLAP_CORE_FRACTION_DIV,
                                     span_floor)) {
                            reciprocal.emplace(
                                k,
                                Overlap{match.other, match.self,
                                        CoreSpan{std::get<1>(*best),
                                                 std::get<2>(*best),
                                                 rinfo.total_cores}});
                        } else {
                            ++non_recip;
                        }
                    } else {
                        ++non_recip;
                    }
                }
            }
        }
        if (getenv("LASM_ADJ_JOIN")) {
            // DBG adjacency as join evidence for OLC-ineligible crumbs: nodes
            // whose core span is below olc_lcmer_l+1 cannot re-derive enough
            // cores to form lcmers and are invisible to the overlap stage
            // (measured: 1086 of 1524 active nodes on chr18). Real reads
            // traversed the edge, so attach such crumbs to their
            // OLC-eligible neighbors at the boundary. Two OLC-eligible nodes
            // never join via adjacency alone -- the tangle misjoin risk is
            // bounded by crumb-sized spans.
            const u64 min_cores = static_cast<u64>(olc_lcmer_l + 1);
            size_t adj_joins = 0;
            const bool adj_dump = getenv("LASM_ADJ_DUMP") != nullptr;
            const char *adj_fasta = getenv("LASM_ADJ_FASTA");
            // Pileup guard: attaching a crumb to a collapsed repeat hub
            // (partner depth far above haploid) joins it across genomic
            // copies -- the measured chimera source. Median over live nodes.
            gtl::vector<double> depths;
            depths.reserve(dbgi.active_set.size());
            for (u64 nid : dbgi.active_set) {
                auto nit = dbgi.dbg.find(nid);
                if (nit == dbgi.dbg.end() || nit->second.redirected() ||
                    nit->second.seq.empty()) {
                    continue;
                }
                depths.push_back(nit->second.avg_core_count());
            }
            double depth_median = 0.0;
            if (!depths.empty()) {
                std::sort(depths.begin(), depths.end());
                depth_median =
                    0.5 * (depths[(depths.size() - 1) / 2] +
                           depths[depths.size() / 2]);
            }
            const double kPileupMaxMult = 3.0;
            const double pileup_cap = depth_median * kPileupMaxMult;
            for (const auto &[id, node] : dbgi.dbg) {
                if (node.redirected()) continue;
                const u64 span = static_cast<u64>(node.seq.size());
                if (span >= min_cores) continue;  // OLC-eligible: skip
                if (span == 0) continue;
                // Hub guard: a crumb with >=2 big prevs would join two
                // OLC-eligible nodes through a tiny span (chimera risk in
                // tangles). Only single-prev crumbs attach along the prev
                // direction; multi-prev merge points are left to OLC.
                if (node.prev.size() != 1) continue;
                const u64 p = node.prev[0];
                auto pit = dbgi.dbg.find(p);
                if (pit == dbgi.dbg.end() || pit->second.redirected()) {
                    continue;
                }
                const u64 pspan = static_cast<u64>(pit->second.seq.size());
                auto write_seq = [&](FILE *f, const char *tag, u64 nid,
                                     const auto &seqs) {
                    std::fprintf(f, ">%s_%llu\n", tag,
                                 (unsigned long long)nid);
                    for (const auto &s : seqs) {
                        std::vector<char> buf(s.len + 1);
                        s.convert_to_ascii(buf.data());
                        fwrite(buf.data(), 1, s.len, f);
                    }
                    fputc('\n', f);
                };
                auto attach = [&](u64 q, u64 qspan, const char *dir) {
                    auto qit0 = dbgi.dbg.find(q);
                    if (qit0 == dbgi.dbg.end() ||
                        qit0->second.redirected() ||
                        qit0->second.seq.empty()) {
                        return;
                    }
                    const double q_depth = qit0->second.avg_core_count();
                    if (pileup_cap > 0.0 && q_depth > pileup_cap) {
                        return;  // repeat hub: ambiguous across copies
                    }
                    // Zero-width tail placement: the crumb abuts the
                    // partner's end (no overlap). Partner-side span {0,0,q}
                    // puts the crumb at partner.start + qspan in the BFS;
                    // crumb-side target span {q,0,q} places the partner at
                    // crumb.start - qspan in the reverse direction.
                    reciprocal.emplace(
                        q, Overlap{id, CoreSpan{0, 0, qspan},
                                   CoreSpan{0, span, span}});
                    reciprocal.emplace(
                        id, Overlap{q, CoreSpan{0, span, span},
                                    CoreSpan{qspan, 0, qspan}});
                    ++adj_joins;
                    if (adj_dump) {
                        fprintf(stderr,
                                "[adj] crumb=%llu partner=%llu dir=%s "
                                "span=%zu pspan=%zu\n",
                                (unsigned long long)id,
                                (unsigned long long)q, dir, (size_t)span,
                                (size_t)qspan);
                    }
                    if (adj_fasta) {
                        FILE *f = std::fopen(adj_fasta, "a");
                        write_seq(f, "crumb", id, node.seq);
                        auto qit = dbgi.dbg.find(q);
                        if (qit != dbgi.dbg.end()) {
                            write_seq(f, "partner", q, qit->second.seq);
                        }
                        fclose(f);
                    }
                };
                attach(p, pspan, "prev");
                // Forward attachments: the last hop crumb->big lives in the
                // big node's prev list, invisible from the crumb's own prevs.
                const auto fwds = GraphCleaning::build_live_fwd(dbgi.dbg);
                auto fit = fwds.find(id);
                if (fit != fwds.end()) {
                    u64 big_fwds = 0;
                    for (u64 q : fit->second) {
                        auto qit = dbgi.dbg.find(q);
                        if (qit != dbgi.dbg.end() &&
                            !qit->second.redirected() &&
                            static_cast<u64>(qit->second.seq.size()) >=
                                min_cores) {
                            ++big_fwds;
                        }
                    }
                    if (big_fwds < 2) {
                        for (u64 q : fit->second) {
                            auto qit = dbgi.dbg.find(q);
                            if (qit == dbgi.dbg.end() ||
                                qit->second.redirected()) {
                                continue;
                            }
                            attach(
                                q,
                                static_cast<u64>(qit->second.seq.size()),
                                "fwd");
                        }
                    }
                }
            }
            fprintf(stderr, "[adj-join] added %zu crumb joins\n", adj_joins);
        }
        if (getenv("LASM_DUMP_MATCHMULT")) {
            // Distinct reciprocal partners per node: the tangle-vs-unique
            gtl::flat_hash_map<u64, gtl::flat_hash_set<u64>> partners;
            for (const auto &[k, ov] : reciprocal) {
                partners[k].insert(ov.target);
            }
            size_t h1 = 0, h2 = 0, h3 = 0, h5 = 0, h10 = 0, hbig = 0;
            for (const auto &[k, v] : partners) {
                const size_t n = v.size();
                if (n <= 1) ++h1;
                else if (n <= 2) ++h2;
                else if (n <= 3) ++h3;
                else if (n <= 5) ++h5;
                else if (n <= 10) ++h10;
                else ++hbig;
            }
            fprintf(stderr,
                    "[match-mult] nodes by distinct-partners: 1:%zu 2:%zu "
                    "3:%zu 4-5:%zu 6-10:%zu >10:%zu\n",
                    h1, h2, h3, h5, h10, hbig);
        }
        if (getenv("LASM_DUMP_SPANCOV")) {
            // Per-(query, other): run_count, span_sum, max_run, query_cores.
            // Same-locus overlaps are one near-total run; divergent-copy
            // overlaps fragment into many short runs.
            const std::string sc_path = output + ".spancov.tsv";
            FILE *sc = fopen(sc_path.c_str(), "w");
            if (sc) {
                fprintf(sc, "query\tother\truns\tspan_sum\tmax_run\tquery_cores\n");
                for (const auto &[k, v] : matched) {
                    gtl::flat_hash_map<u64, gtl::vector<u64>> per_other;
                    u64 qcores = 0;
                    for (const auto &m : v) {
                        per_other[m.other].push_back(m.self.to - m.self.from);
                        qcores = std::max(qcores, m.self.len);
                    }
                    for (const auto &[o, spans] : per_other) {
                        u64 sum = 0, mx = 0;
                        for (u64 s : spans) {
                            sum += s;
                            mx = std::max(mx, s);
                        }
                        fprintf(sc,
                                "%" PRIu64 "\t%" PRIu64 "\t%zu\t%" PRIu64
                                "\t%" PRIu64 "\t%" PRIu64 "\n",
                                k, o, spans.size(), sum, mx, qcores);
                    }
                }
                fclose(sc);
            }
        }
        fprintf(stderr, "all:%d\ngood: %d nodes\nhits: %d\n", all_n.load(),
                all_n.load() - bad.load(), hits.load());
        fprintf(stderr,
                "[match-diag] too_short:%d no_ranges:%d frac_drop:%d "
                "recip_kept:%d relaxed_kept:%d non_recip:%d\n",
                too_short.load(), no_ranges.load(), frac_drop.load(),
                recip_kept, relaxed_kept, non_recip);
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

    GraphCleaning::unit_mult_fn_t unit_mult_fn;
    if (cfg.graph_cleaning.enable_loop_opening) {
        // Pre-simplify per-lcmer multiplicity snapshot: the oracle used by
        // unroll_loops to check that a loop unit's internal lcmers actually
        // occur at the predicted copy number in the reads (degeneracy-
        // artifact cycles over low-complexity sequence have inflated node
        // depths but ~1x internal lcmer multiplicities).
        gtl::flat_hash_map<u64, unsigned> lcmer_mult;
        gtl::vector<unsigned> mults;
        lcmer_mult.reserve(dbg.dbg.size());
        mults.reserve(dbg.dbg.size());
        for (const auto &[id, node] : dbg.dbg) {
            if (node.redirected()) continue;
            lcmer_mult.emplace(id, node.counts.core_count);
            mults.push_back(node.counts.core_count);
        }
        if (mults.size() >= 2) {
            std::sort(mults.begin(), mults.end());
            const double mult_median =
                0.5 * (mults[(mults.size() - 1) / 2] +
                       mults[mults.size() / 2]);
            if (mult_median > 0.0) {
                unit_mult_fn =
                    [&, lcmer_l = as.dbg_lcmer_l](
                        const gtl::vector<SequenceContainerType> &s_seq,
                        const gtl::vector<SequenceContainerType> &arm_seq)
                    -> double {
                    gtl::vector<SequenceContainerType> unit_seqs;
                    unit_seqs.reserve(s_seq.size() + arm_seq.size());
                    for (const auto &s : s_seq) {
                        unit_seqs.push_back(
                            SequenceContainerType::copy_of(s));
                    }
                    for (const auto &s : arm_seq) {
                        unit_seqs.push_back(
                            SequenceContainerType::copy_of(s));
                    }
                    lps str = {0, 0, 0, 0};
                    lcp_read lr = assembler::build_lr_cores(
                        unit_seqs, as.syncmer_functor, as.lcp_level, str);
                    std::vector<double> vals;
                    if (lr.cores.size() >
                        static_cast<size_t>(lcmer_l + 1)) {
                        LCmer l{lr.cores.begin() + 1, lcmer_l};
                        auto lookup = [&](const LCmer &lc) {
                            auto mit = lcmer_mult.find(lc.value());
                            if (mit != lcmer_mult.end()) {
                                vals.push_back(
                                    static_cast<double>(mit->second));
                            }
                        };
                        lookup(l);
                        for (size_t i = static_cast<size_t>(lcmer_l + 1);
                             i < lr.cores.size() - 1; ++i) {
                            l.update_lcmer(lr.cores[i],
                                           lr.cores[i - lcmer_l], lcmer_l);
                            lookup(l);
                        }
                    }
                    if (vals.size() < 4) return 1.0;  // too few: no evidence
                    std::sort(vals.begin(), vals.end());
                    const double med =
                        0.5 * (vals[(vals.size() - 1) / 2] +
                               vals[vals.size() / 2]);
                    return med / mult_median;
                };
            }
        }
    }
    GraphCleaning::run(dbg, cfg, unit_mult_fn ? &unit_mult_fn : nullptr);

    IF_TIME_CHECKPOINT(!cfg.silent, HIFI_ASSEMBLY, stderr, "Graph cleaning");

    // Diagnostic dump (LASM_DUMP_NODEDEPTH=1): per-node stats for repeat-
    // boundary calibration — avg_core_count, seq bases, fragment count, in/out
    // degree. Fragments-per-base scales with collapsed copy number even though
    // avg_core_count (core_count/frags) stays flat under collapse.
    if (getenv("LASM_DUMP_NODEDEPTH")) {
        const std::string nd_path = cfg.output + ".nodedepth.tsv";
        FILE *nd = fopen(nd_path.c_str(), "w");
        if (nd) {
            fprintf(nd, "id\tavg_core\tseq_bases\tfrags\tprev\tfec\n");
            for (u64 node : dbg.active_set) {
                auto it = dbg.dbg.find(node);
                if (it == dbg.dbg.end() || it->second.redirected() ||
                    it->second.seq.empty()) {
                    continue;
                }
                size_t bases = 0;
                for (const auto &s : it->second.seq) bases += s.len;
                fprintf(nd,
                        "%" PRIu64 "\t%.2f\t%zu\t%zu\t%zu\t%u\n", node,
                        it->second.avg_core_count(), bases,
                        it->second.seq.size(), it->second.prev.size(),
                        it->second.forward_edge_count);
            }
            fclose(nd);
        }
    }

    {
        FILE *out = fopen(cfg.output.c_str(), "w");
        if (!out) {
            fprintf(stderr, "Cannot open output: %s\n", cfg.output.c_str());
            return 1;
        }
        layout_graph_emit::emit_ultra_merged_graph(as, dbg, out);
        fclose(out);
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
