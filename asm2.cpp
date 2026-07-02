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

struct assembler {
    using SYNCMER_HASHER = syncmer::MurMur;
    using SyncmerType = syncmer::GenericSyncmerRule<uint32_t, SYNCMER_HASHER>;

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
          syncmer_functor{syncmer_params, static_cast<unsigned int>(seed)},
          dbg_lcmer_l{cfg.dbg_lcmer_l},
          olc_lcmer_l{cfg.dbg_lcmer_l},
          parallel_file_reading{cfg.parallel_file_reading} {
        inputs.assign(cfg.inputs.begin(), cfg.inputs.end());
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
                                  int lcmer_l) {
        long skip_first = lcmer_l;
        long skip_last = 2 * lcmer_l;
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
    struct rev_match {
        u64 node;
        struct ival {
            u64 from;
            u64 to;
            u64 size;
        };
        ival tgt;
        ival qry;

        // Matching LCMER indices.
        //
        // Convention:
        // - *_query_* refer to the QUERY/SOURCE node (the map key / "from"
        // node)
        // - *_target_* refer to the TARGET node (rev_match::node)
        //
        // These indices are later used to "switch" between SequenceContainer
        // streams at overlap boundaries without base-offset arithmetic.
        LCmer first_matching_lcmer_query;
        LCmer last_matching_lcmer_query;
        LCmer first_matching_lcmer_target;
        LCmer last_matching_lcmer_target;
    };

    auto btg_connected_components(gtl::btree_multimap<u64, rev_match> &recig) {
        gtl::vector<gtl::vector<u64>> components;
        gtl::flat_hash_set<u64> visited;
        for (const auto &[from, to] : recig) {
            gtl::vector<u64> component;
            gtl::vector<u64> stack;
            stack.push_back(from);
            while (!stack.empty()) {
                u64 node = stack.back();
                stack.pop_back();
                if (visited.find(node) != visited.end()) {
                    continue;
                }
                visited.emplace(node);
                component.push_back(node);
                auto range = recig.equal_range(node);
                for (auto it = range.first; it != range.second; ++it) {
                    if (visited.find(it->second.node) == visited.end()) {
                        stack.push_back(it->second.node);
                    }
                }
            }
            if (component.size() > 1) {
                components.push_back(component);
            }
        }
        return components;
    }


    auto add_lcmers(const gtl::vector<SequenceContainerType> &subseqs,
                    auto &node, dbg_index &dbgi, u64 &idx, int tid) {
        lps str = {0, 0, 0, 0};
        lcp_read lr;
        // dbgi.dbg.for_each([&](auto &node) {
        if (node.second.redirected()) {
            return;
        }
        // good_n.fetch_add(1);
        lr.clear();

        // Extract a subsequence [from, to) from the concatenated rev_seqs
        // Returns a single CrumbledSubSeq64 containing the extracted range

        init_lps_crumbled(&str, subseqs, syncmer_functor);
        lps_deepen(&str, lcp_level);

        for (int i = 0; i < str.size; ++i) {
            lr.cores.emplace_back(str.cores[i].label);
            lr.pos.emplace_back(str.cores[i].start, str.cores[i].end);
        }


        if (lr.cores.size() < static_cast<size_t>(olc_lcmer_l + 1)) {
            // good_n.fetch_sub(1);
            free_lps(&str);
            return;
        }
        gtl::vector<u64> ll = {0};
        LCmer l{lr.cores.begin() + 1, olc_lcmer_l};  // Skip first core;
        ll.push_back(dbgi.emplace_nobloomfilter(l, idx, tid));
        for (size_t i = 2; i + olc_lcmer_l + 1 < lr.cores.size(); ++i) {
            uint32_t in_core = lr.cores[i + olc_lcmer_l];
            l.update_lcmer(in_core, lr.cores[i], olc_lcmer_l);
            ll.push_back(dbgi.emplace_nobloomfilter(l, idx, tid));
        }
        ll.push_back(0);

        for (size_t i = 1; i < ll.size() - 1; ++i) {
            if (insert_lcmer_helper(ll, i, dbgi)) {
                dbgi.dbg.modify_if(
                    ll[i], [&lr, i, &subseqs,
                            this](dbg_index::dbg_map_t::value_type &val) {
                        val.second.seq.emplace_back(
                            subseqs, lr.pos[i + olc_lcmer_l - 1].second,
                            lr.pos[i + olc_lcmer_l].second);
                    });
            }
        }
    }

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

    auto find_rev_lcmers(const gtl::vector<SequenceContainerType> &subseqs,
                         dbg_index &dbgi, int tid, u64 &idx)
        -> gtl::vector<SequenceContainerType> {
        lps str = {0, 0, 0, 0};
        lcp_read lr;
        lr.clear();

        gtl::vector<SequenceContainerType> rev_seqs;
        for (auto it = subseqs.rbegin(); it != subseqs.rend(); ++it) {
            rev_seqs.emplace_back(it->reverse_complement());
        }
        init_lps_crumbled(&str, rev_seqs, syncmer_functor);
        lps_deepen(&str, lcp_level);

        for (int i = 0; i < str.size; ++i) {
            lr.cores.emplace_back(str.cores[i].label);
            lr.pos.emplace_back(str.cores[i].start, str.cores[i].end);
        }
        if (lr.cores.size() < static_cast<size_t>(olc_lcmer_l + 2)) {
            free_lps(&str);
            return {};
        }
        LCmer l{lr.cores.begin() + 1, olc_lcmer_l};  // Skip first core;
        gtl::vector<u64> ll;

        u64 lval = dbgi.emplace_nobloomfilter(l, idx, tid);
        ll.emplace_back(lval);
        for (size_t i = 2; i + olc_lcmer_l + 1 < lr.cores.size(); ++i) {
            uint32_t in_core = lr.cores[i + olc_lcmer_l];
            l.update_lcmer(in_core, lr.cores[i], olc_lcmer_l);
            lval = dbgi.emplace_nobloomfilter(l, idx, tid);
            ll.emplace_back(lval);
        }
        gtl::vector<SequenceContainerType> seqs;
        for (size_t i = 1; i < ll.size() - 1; ++i) {
            seqs.emplace_back(rev_seqs, lr.pos[i + olc_lcmer_l - 1].second,
                              lr.pos[i + olc_lcmer_l].second);
            seqs.back().origin_node = ll[i];
        }
        // TODO: index it if needed
        return seqs;
    }
    auto fold_dbg_graph(dbg_index &dbgi) {
        fprintf(stderr, "dbg_size: %zu\n", dbgi.dbg.size());

        std::atomic<int> all_n{0};
        std::atomic<int> hits{0};
        std::atomic<int> bad{0};
        gtl::parallel_flat_hash_map_d<u64, gtl::vector<rev_match>, HT_SHARD_BC,
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
                auto found_ranges = find_lcmer_ranges(lr, dbgi, olc_lcmer_l);
                for (const auto &[id, from_core, to_core, from_lcmer,
                                  to_lcmer] : found_ranges) {
                    if (id != 0) {
                        if (to_core - from_core <
                            lr.pos.size() / 20) {  // TODO: parameterize this
                            continue;
                        }
                        hits.fetch_add(1);

                        // Store match intervals in CORE indices (inclusive).
                        // IMPORTANT: rev_match fields are u64 and we keep them
                        // semantically as core-space indices here. Base-space
                        // scaling happens at emission time.
                        //
                        // We also store "length" fields in core-space to keep
                        // conversions consistent. Core-space length is number
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
                                    rev_match::ival{static_cast<u64>(from_core),
                                                    static_cast<u64>(to_core),
                                                    core_len_steps},
                                    rev_match::ival{0, 0, 0}, from_lcmer,
                                    to_lcmer, LCmer{}, LCmer{});
                            },
                            [&](const auto &ctr) {
                                ctr(node.first,
                                    gtl::vector<rev_match>{rev_match{
                                        id,
                                        rev_match::ival{
                                            static_cast<u64>(from_core),
                                            static_cast<u64>(to_core),
                                            core_len_steps},
                                        rev_match::ival{0, 0, 0}, from_lcmer,
                                        to_lcmer, LCmer{}, LCmer{}}});
                            });
                    }
                }
                free_lps(&str);
            },
            threads.work);

        gtl::btree_multimap<u64, rev_match> reciprocal;
        for (const auto &[k, v] : matched) {
            for (const auto &match : v) {
                // if (k < match.node) {  // smaller node is the key
                auto it = matched.find(
                    match.node);  // and matched node also has overlaps
                if (it != matched.end()) {
                    auto ff = std::find_if(
                        it->second.begin(), it->second.end(),
                        [k](const auto &val) { return val.node == k; });
                    if (ff != it->second.end()) {  // and has the reverse match
                        // Store both query and target coordinates in CORE-INDEX
                        // space. match.from/to are core indices on k's rev_comp
                        // (query) ff->from/to are core indices on match.node's
                        // rev_comp (target)
                        //
                        // NOTE: rev_match uses generic u64 fields; here we
                        // interpret them as core-space units. Base-space
                        // scaling will be applied only when emitting merged
                        // sequences / reporting.
                        rev_match full_match;
                        full_match.node = match.node;
                        full_match.tgt.from =
                            ff->tgt.from;  // ts (target start) in cores
                        full_match.tgt.to =
                            ff->tgt.to;  // te (target end) in cores
                        full_match.tgt.size =
                            ff->tgt.size;  // tlen in cores (steps)
                        full_match.qry.from =
                            match.tgt.from;  // qs (query start) in cores
                        full_match.qry.to =
                            match.tgt.to;  // qe (query end) in cores
                        full_match.qry.size =
                            match.tgt.size;  // qlen in cores (steps)
                        full_match.first_matching_lcmer_query =
                            match.first_matching_lcmer_query;
                        full_match.last_matching_lcmer_query =
                            match.last_matching_lcmer_query;
                        full_match.first_matching_lcmer_target =
                            ff->first_matching_lcmer_query;
                        full_match.last_matching_lcmer_target =
                            ff->last_matching_lcmer_query;
                        reciprocal.emplace(k, full_match);

                    } else {  // Handle after reciprocals are good
                    }
                }
            }
        }
        fprintf(stderr, "all:%d\ngood: %d nodes\nhits: %d\n", all_n.load(),
                all_n.load() - bad.load(), hits.load());
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
        const std::string raw_path = cfg.output + ".raw";
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


    // Dump pre-layout unitig graph
    {
        const std::string pre_path = cfg.output + ".pre_layout";
        FILE *pre = fopen(pre_path.c_str(), "w");
        if (!pre) {
            fprintf(stderr, "Cannot open output: %s\n", pre_path.c_str());
            return 1;
        }
        print_dbg_seqs(pre, dbg);
        fclose(pre);
    }

    // Later stages commented out for profiling graph building
    {
        FILE *out = fopen(cfg.output.c_str(), "w");
        if (!out) {
            fprintf(stderr, "Cannot open output: %s\n", cfg.output.c_str());
            return 1;
        }
        // print_dbg_seqs(out, dbg);
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
    using SyncmerType = syncmer::GenericSyncmerRule<uint32_t, SYNCMER_HASHER>;
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

        auto merged = collect_path_sequences(paths, node_contigs,
                                             contig_core_seqs, 3);
        REQUIRE_EQ(merged.size(), 2u);
        CHECK_EQ(merged[0].len, 4u);
        CHECK_EQ(merged[1].len, 4u);
    }

    SUBCASE("two components separated by a zero-length gap") {
        gtl::vector<comp_path> paths;
        paths.push_back({0, 100, 2, gtl::vector<uint32_t>{0u, 1u}});
        paths.push_back({200, 300, 1, gtl::vector<uint32_t>{2u}});

        auto merged = collect_path_sequences(paths, node_contigs,
                                             contig_core_seqs, 3);
        // seqA, seqB, GAP, seqC
        REQUIRE_EQ(merged.size(), 4u);
        CHECK_EQ(merged[0].len, 4u);
        CHECK_EQ(merged[1].len, 4u);
        CHECK_EQ(merged[2].len, 0u);   // gap sentinel between components
        CHECK_EQ(merged[3].len, 6u);
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

TEST_CASE("layout_graph_emit::find_heaviest_paths drops non-reconnecting "
          "backward extension (C13)") {
    using namespace layout_graph_emit;

    // Main heaviest path A(0)→B(1)→C(2). An earlier source D(3) branches to a
    // heavy dead-end E(4) and also to A. The greedy backward walk from D prefers
    // E (higher dp) and dead-ends there before reaching the path, so nothing
    // should be prepended — the result must stay the connected chain [0,1,2].
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
    // NOT prepended, since the greedy backward walk dead-ends at E before
    // reaching the path (0→1 and 1→2 are the only edges among these nodes).
    REQUIRE_EQ(paths.size(), 1u);
    REQUIRE_EQ(paths[0].path.size(), 3u);
    CHECK_EQ(paths[0].path[0], 0u);
    CHECK_EQ(paths[0].path[1], 1u);
    CHECK_EQ(paths[0].path[2], 2u);
}
