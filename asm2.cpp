#include <blockingconcurrentqueue.h>
#include <fcntl.h>
#include <htslib/bgzf.h>
#include <htslib/hfile.h>
#include <htslib/hts.h>
#include <htslib/sam.h>
#include <htslib/thread_pool.h>
#include <stdio.h>

#include <Eigen/Dense>
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
#include <limits>
#include <optional>
#include <queue>
#include <stack>
#include <string>
#include <tuple>
#include <utility>

#include "asm.h"
#include "bgzf.h"
#include "config.h"
#include "ext/doctest.h"
#include "graph.h"
#include "lcmer.h"
#include "lcptools/lcptools_ho.h"
#include "syncmer.h"

using std::string;
enum class input_type { SAM, BAM, FASTX, FASTXGZ, FASTXBZ2, UNKNOWN };

template <class T>
struct mptr {
    static constexpr size_t vmask = 3;
    static constexpr size_t pmask = ~vmask;
    T *ptr;
    ~mptr() {
        delete unwrap();
        ptr = 0;
    }
    mptr(T val) : ptr{new T{val}} {}
    mptr(T val, int message) : ptr{new T{val}} { write_message(message); }
    mptr(T *ptr) : ptr{ptr} {}
    mptr(T *ptr, int message) : ptr{ptr} { write_message(message); }
    mptr(mptr &&other) {
        ptr = other.ptr;
        other.ptr = 0;
    }
    mptr(const mptr &other) = delete;
    mptr &operator=(const mptr &other) = delete;
    mptr &&operator=(mptr &&other) {
        this->ptr = other.ptr;
        other.ptr = 0;
    }
    constexpr const T *unwrap() const { return (T *)((size_t)(ptr)&pmask); }
    constexpr T *unwrap() { return (T *)((size_t)(ptr)&pmask); }
    constexpr T operator*() const {
        if (ptr) return *unwrap();
        throw 42;
    }
    constexpr int get_message() const { return (size_t)(ptr)&vmask; }
    void write_message(int message) {
        ptr = (T *)((size_t)(ptr) | (message & vmask));
    }
};

struct assembler {
    using SYNCMER_HASHER = syncmer::MurMur;
    using SyncmerType = syncmer::GenericSyncmerRule<uint32_t, SYNCMER_HASHER>;

    using cqt = moodycamel::ConcurrentQueue<Chunk>;
    using dqt = moodycamel::BlockingConcurrentQueue<bam1_t *>;
    string input;
    string output;
    decltype(asm_config::threads) threads;
    int lcp_level;
    double density;
    decltype(asm_config::bloom) bloom;
    uint64_t seed;
    decltype(asm_config::syncmer) syncmer_params;
    SyncmerType syncmer_functor;

    assembler(const asm_config &cfg)
        : input{cfg.input},
          output{cfg.output},
          threads{cfg.threads},
          lcp_level{cfg.lcp_level},
          bloom{cfg.bloom},
          seed{cfg.seed},
          syncmer_params{cfg.syncmer},
          syncmer_functor{syncmer_params, static_cast<unsigned int>(seed)} {}

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
                        val.second.core_count++;
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
                        val.second.core_count++;
                    },
                    [&](const index_type::dbg_map_t::constructor &ctr) {
                        ctr(ll[i],
                            typename index_type::bidbg_val{1U, ll[i - 1]});
                    });
            }
        }
    }
    static gtl::vector<uint64_t> find_lcmers(const lcp_read &lr,
                                             const dbg_index &dbgi) {
        lcmer<LCMER_L> l{lr.cores.begin() + 1};  // Skip first core;

        gtl::vector<uint64_t> ll(1, 0);
        ll.reserve(lr.cores.size() - LCMER_L + 3);
        dbg_index::lci_type lval = dbgi.get(l);
        ll.push_back(lval.id);
        for (size_t i = LCMER_L + 1; i < lr.cores.size() - 1; ++i) {
            uint32_t in_core = lr.cores[i];
            l.update_lcmer(in_core);
            lval = dbgi.get(l);
            ll.push_back(lval.id);
        }
        ll.push_back(0);
        return ll;
    }
    static auto find_lcmer_ranges(const lcp_read &lr, const dbg_index &dbgi) {
        long skip_first = LCMER_L;
        long skip_last = 2 * LCMER_L;
        // Need enough cores to form at least one valid lcmer after skipping
        if (lr.cores.size() <
            static_cast<size_t>(LCMER_L + skip_first + skip_last)) {
            return gtl::vector<
                std::tuple<u64, u64, u64, lcmer<LCMER_L>, lcmer<LCMER_L>>>{};
        }

        lcmer<LCMER_L> l{lr.cores.begin() + skip_first};  // Skip first core;
        gtl::vector<std::tuple<u64, u64, u64, lcmer<LCMER_L>, lcmer<LCMER_L>>>
            llr;

        dbg_index::lci_type lval = dbgi.get(l);
        llr.emplace_back(lval.id, skip_first, skip_first, l, l);
        // Fix: ensure we don't go out of bounds when accessing i + LCMER_L
        for (long long i = skip_first + 1;
             i < ((long long)lr.cores.size() - LCMER_L - skip_last); ++i) {
            uint32_t in_core = lr.cores[i + LCMER_L - 1];
            l.update_lcmer(in_core);
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
                              int thread_id, uint64_t &lcmer_index) {
        lcmer<LCMER_L> l{lr.cores.begin() + 1};  // Skip first core;
        gtl::vector<std::tuple<u64, u64, u64>> llr;

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
            if (insert_lcmer_helper(ll, i, dbgi)) {
                // Lets try atomic
                dbgi.dbg.modify_if(
                    ll[i],
                    [&lr, i](dbg_index_base<BFT>::dbg_map_t::value_type &val) {
                        val.second.seq.emplace_back(
                            lr.data, lr.pos[i + LCMER_L - 1].second,
                            lr.pos[i + LCMER_L].second, val.first);
                    });
            }
        }
    }
    static void process_cores_for_basespace(
        const lcp_read &lr, dbg_index &dbgi, [[maybe_unused]] int thread_id,
        [[maybe_unused]] uint64_t &lcmer_index) {
        // Keep the DBG node sequences in BASE space (actual substrings), but do
        // NOT convert overlap/layout coordinates to base space during
        // layouting/merging.
        //
        // Layouting will be done in core-index space. Only at the final
        // emission step will we scale core-space layout positions into base
        // offsets.
        lcmer<LCMER_L> l{lr.cores.begin() + 1};  // Skip first core;

        uint64_t val = dbgi.get(l).id;

        for (size_t i = LCMER_L + 1; i < lr.cores.size() - 1; ++i) {
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
            l.update_lcmer(in_core);
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

    // static auto process_data_from_chunks(int thread_id, int lcp_level,
    //                                      string filename, cqt &chunk_queue,
    //                                      dbg_index &dbgi,
    //                                      std::atomic<bool> &done,
    //                                      int &processed_read_count,
    //                                      int &all_reads) {
    //     uint64_t lcmer_index = 0;
    // }
    static auto process_data(
        int thread_id, int lcp_level, dqt &data_queue, dqt &recycle_queue,
        dbg_index &dbgi, std::atomic<bool> &done, int &processed_read_count,
        int &all_reads,
        decltype(assembler::process_cores<DefaultBF>) Process_Func,
        SyncmerType syncmer_functor) {
        uint64_t lcmer_index = 0;
        bam1_t *data_p = NULL;
        lps str = {0, 0, 0, 0};
        charsi chs;
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
                                        syncmer_functor, chs, lr);
            if (lr.cores.size() < LCMER_L + 2) {
                // printf("Skipping read %lld, cores.size()=%zu\n", here,
                // lr.cores.size());
                goto pd_cleanup;
            }
            ++processed_read_count;
            Process_Func(lr, dbgi, thread_id, lcmer_index);
        pd_cleanup:
            if (!recycle_queue.try_dequeue(data_p)) {
                bam_destroy1(data_p);
                data_p = NULL;
            }
        }
        free_lps(&str);
    }

    auto generate_chunks(std::string input, cqt &chunk_queue) {
        BlockIterator iter(input);
        Block blk;
        // fetch chunks until no more blocks
        while (iter.next(blk)) {
            // first block in chunk
            int64_t first_vstart = blk.vstart();
            int64_t last_vend = blk.vend();

            size_t count = 1;
            while (count < BAMREAD_CHUNK_SIZE && iter.next(blk)) {
                last_vend = blk.vend();
                ++count;
            }
            while (!chunk_queue.try_enqueue({first_vstart, last_vend}));
        }
    }

    [[nodiscard]] auto generate_data(const std::string &input, dqt &data_queue,
                                     dqt &recycle_queue, int threads) {
        samFile *fp = sam_open(input.c_str(), "r");

        if (!fp) {
            std::cerr << "Error: Could not open " << input << std::endl;
            return false;
        }
        hts_set_threads(fp, threads);
        // htsThreadPool pool = {NULL, 0};
        // if (threads > 1) {
        //     pool.pool = hts_tpool_init(threads);
        //     // hts_set_thread_pool(fp, &pool);
        //     hts_set_cache_size(fp, 1024 * 1024 * 1024);
        //     hts_set_opt(fp, HTS_OPT_THREAD_POOL, &pool);
        // }

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

    auto construct_dbg() {
        dqt data_queue{static_cast<size_t>(threads.work) * 1024, 1, 1};
        dqt recycle_queue{static_cast<size_t>(threads.work) * 1024, 1, 1};
        std::atomic<bool> done{false};

        dbg_index *dbg_pt = new dbg_index{bloom.false_positive_rate,
                                          bloom.projected_count, seed};
        dbg_index &dbgi = *dbg_pt;

        gtl::vector<std::thread> workers;

        gtl::vector<int> valid_read_counts(threads.work, 0);
        gtl::vector<int> all_read_counts(threads.work, 0);

        for (int tid = 1; tid <= threads.work; ++tid) {
            workers.emplace_back(
                process_data, tid, lcp_level, std::ref(data_queue),
                std::ref(recycle_queue), std::ref(dbgi), std::ref(done),
                std::ref(valid_read_counts[tid - 1]),
                std::ref(all_read_counts[tid - 1]), process_cores<DefaultBF>,
                std::cref(syncmer_functor));
        }

        if (!generate_data(input, std::ref(data_queue), std::ref(recycle_queue),
                           threads.io)) {
            //Should print error in generate_data
            exit(-1);
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
        fprintf(stderr, "%llu processed\n%llu total reads\n", valid_sum,
                all_sum);
        return dbg_pt;
    }

    auto convert2basespace(dbg_index &dbgi) -> int {
        dqt data_queue{static_cast<size_t>(threads.work) * 1024, 1, 1};
        dqt recycle_queue{static_cast<size_t>(threads.work) * 1024, 1, 1};
        std::atomic<bool> done{false};

        gtl::vector<std::thread> workers;
        gtl::vector<int> valid_read_counts(threads.work, 0);
        gtl::vector<int> all_read_counts(threads.work, 0);

        for (int tid = 1; tid <= threads.work; ++tid) {
            workers.emplace_back(
                process_data, tid, lcp_level, std::ref(data_queue),
                std::ref(recycle_queue), std::ref(dbgi), std::ref(done),
                std::ref(valid_read_counts[tid - 1]),
                std::ref(all_read_counts[tid - 1]), process_cores_for_basespace,
                std::cref(syncmer_functor));
        }
        if (!generate_data(input, std::ref(data_queue), std::ref(recycle_queue),
                           threads.io)) {
            fprintf(stderr, "Error in reading data for basespace conversion\n");
            exit(-1);
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
    // I need to recompute the LCP, it is likely more expensive to find relavant
    // nodes.
    // void iterative_re_dbg(dbg_index &dbgi) {
    //     gtl::vector<std::pair<u64, size_t>> node_id_len;
    //     for (const auto &[k, v] : dbgi.dbg) {
    //         node_id_len.emplace_back(k, v.seq.size());
    //     }
    //     std::sort(
    //         node_id_len.begin(), node_id_len.end(),
    //         [](const auto &a, const auto &b) { return a.second < b.second;
    //         });
    //     // Check this again TODO.
    //     dbg_index new_index(dbgi.bf.implementation.m,
    //     dbgi.bf.implementation.n,
    //                         dbgi.bf.random_seed);

    //     auto it = node_id_len.begin();
    //     assert(it != node_id_len.end());
    //     lps str;

    //     auto move_forward_unitig = [&](u64 id) {
    //         auto node_it = dbgi.dbg.find(id);
    //         if (node_it != dbgi.dbg.end()) {
    //             int len = count_seqv_len(node_it->second.seq);
    //             prepare_lps(&str, len);
    //             init_lps_crumbled(&str, node_it->second.seq,
    //             syncmer_functor); lps_deepen(&str, lcp_level);
    //         }
    //     };
    //     auto rev_contig_present_in_the_index =
    //         [&](const dbg_index::dbg_map_t::value_type &node, auto ctg) {
    //             return true;
    //         };
    //     move_forward_unitig(node_id_len.front().first);
    //     ++it;
    //     while (it != node_id_len.end()) {
    //         auto rev_ctg = calculate_reverse_utig(*it);
    //         if (rev_contig_present_in_the_index(*it, rev_ctg)) {
    //         } else {
    //             move_forward_unitig(it->first);
    //         }
    //         ++it;
    //     }
    // }
    struct rev_match {
        u64 node;
        struct ival {
            u64 from;
            u64 to;
            u64 size;
        };
        ival tgt;
        ival qry;
        // u64 from;  // target start in bases (ts)
        // u64 to;    // target end in bases (te)
        // u64 size;  // target length in bases (tlen)
        // u64 qs;    // query start in bases
        // u64 qe;    // query end in bases
        // u64 qlen;  // query length in bases (steps)

        // Matching LCMER indices.
        //
        // Convention:
        // - *_query_* refer to the QUERY/SOURCE node (the map key / "from"
        // node)
        // - *_target_* refer to the TARGET node (rev_match::node)
        //
        // These indices are later used to "switch" between SequenceContainer
        // streams at overlap boundaries without base-offset arithmetic.
        lcmer<LCMER_L> first_matching_lcmer_query;
        lcmer<LCMER_L> last_matching_lcmer_query;
        lcmer<LCMER_L> first_matching_lcmer_target;
        lcmer<LCMER_L> last_matching_lcmer_target;
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

    std::tuple<int, int> get_start_positions(
        u64 current, u64 neighbor, gtl::btree_multimap<u64, rev_match> rmatch,
        int cgroup) {
        int qs = -9999, ts = -9999;
        int qe = -9999, te = -9999;
        int ql = -9999, tl = -9999;

        auto range_c = rmatch.equal_range(current);
        for (auto it = range_c.first; it != range_c.second; ++it) {
            if (it->second.node == neighbor) {
                qe = it->second.tgt.to;
                qs = it->second.tgt.from;
                ql = it->second.tgt.size;
                break;
            }
        }
        // Find the specific reciprocal entry for (neighbor -> current)
        auto range_t = rmatch.equal_range(neighbor);
        for (auto it = range_t.first; it != range_t.second; ++it) {
            if (it->second.node == current) {
                ts = it->second.tgt.from;
                te = it->second.tgt.to;
                tl = it->second.tgt.size;
                break;
            }
        }
        int ngroup = 1 - cgroup;
        if (cgroup == 1) {
            qs = ql - qe;
        }
        if (ngroup == 1) {
            ts = tl - te;
        }
        return std::make_tuple(qs, ts);
    }

    auto connected_comp_plus_layout(
        const gtl::btree_multimap<u64, rev_match> &recig) {
        // Layout entries are now (start, end, orient) where start/end are
        // inclusive. We infer per-node length from incident edges (qlen for
        // outgoing, size for incoming).
        gtl::vector<gtl::flat_hash_map<u64, std::tuple<int, int, bool>>>
            layouts;
        std::queue<std::pair<u64, rev_match>> q;
        gtl::flat_hash_set<u64> visited;

        auto infer_len_for_node = [&](u64 id) -> int {
            // Prefer outgoing query length (qlen).
            auto out_range = recig.equal_range(id);
            if (out_range.first != out_range.second) {
                const auto &rm = out_range.first->second;
                if (rm.qry.size > 0) return static_cast<int>(rm.qry.size);
            }
            // Fallback: any incoming edge's target length (size).
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
            if (visited.contains(nit.first)) {
                continue;
            }
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
                std::pair<u64, rev_match> node = q.front();
                q.pop();

                u64 src = node.first;
                u64 tgt = node.second.node;

                if (visited.contains(tgt)) {
                    continue;
                }
                visited.emplace(tgt);

                auto [spos, epos, s_orient] = layouts.back().at(src);
                (void)epos;
                bool target_orient = !s_orient;

                // Coordinates are on reverse complement, convert based on
                // layout orientation
                long long qs, ts;
                if (s_orient) {
                    // Source is forward in layout: convert rev_comp coords to
                    // forward
                    qs = node.second.qry.size - node.second.qry.to;
                } else {
                    // Source is reverse in layout: rev_comp coords are correct
                    qs = node.second.qry.from;
                }
                if (target_orient) {
                    // Target is forward in layout: convert rev_comp coords to
                    // forward
                    ts = node.second.tgt.size - node.second.tgt.to;
                } else {
                    // Target is reverse in layout: rev_comp coords are correct
                    ts = node.second.tgt.from;
                }

                int delta = static_cast<int>(qs - ts);
                int target_start = spos + delta;

                const int target_len = infer_len_for_node(tgt);
                const int target_end =
                    to_end_inclusive(target_start, target_len);

                layouts.back().emplace(
                    tgt,
                    std::make_tuple(target_start, target_end, target_orient));

                auto range = recig.equal_range(tgt);
                for (auto it = range.first; it != range.second; ++it) {
                    if (!visited.contains(it->second.node)) {
                        q.push(*it);
                    }
                }
            }
        }
        return layouts;
    }

    /// Least-squares layout computation (based on order_prototype.py
    /// order_lstsq). Instead of BFS traversal, this method considers ALL edges
    /// simultaneously and finds globally optimal positions that minimize
    /// squared error.
    auto connected_comp_plus_layout_lstsq(
        const gtl::btree_multimap<u64, rev_match> &recig) {
        gtl::vector<gtl::flat_hash_map<u64, std::tuple<int, int, bool>>>
            layouts;
        gtl::flat_hash_set<u64> visited;

        auto infer_len_for_node = [&](u64 id) -> int {
            // Prefer outgoing query length (qlen).
            auto out_range = recig.equal_range(id);
            if (out_range.first != out_range.second) {
                const auto &rm = out_range.first->second;
                if (rm.qry.size > 0) return static_cast<int>(rm.qry.size);
            }
            // Fallback: any incoming edge's target length (size).
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

        // First, find all connected components using DFS
        for (const auto &nit : recig) {
            if (visited.contains(nit.first)) {
                continue;
            }

            // Collect all nodes in this connected component
            gtl::vector<u64> component;
            gtl::flat_hash_set<u64> comp_visited;
            std::stack<u64> stack;
            stack.push(nit.first);

            while (!stack.empty()) {
                u64 node = stack.top();
                stack.pop();
                if (comp_visited.contains(node)) {
                    continue;
                }
                comp_visited.emplace(node);
                visited.emplace(node);
                component.push_back(node);

                // Add neighbors (outgoing edges)
                auto range = recig.equal_range(node);
                for (auto it = range.first; it != range.second; ++it) {
                    if (!comp_visited.contains(it->second.node)) {
                        stack.push(it->second.node);
                    }
                }
            }

            if (component.empty()) {
                continue;
            }

            // Build index mapping: node_id -> index
            gtl::flat_hash_map<u64, size_t> idx;
            for (size_t i = 0; i < component.size(); ++i) {
                idx[component[i]] = i;
            }
            size_t n = component.size();

            // Collect edges with deltas and track orientation
            // edges: (src, tgt, delta, src_orient, tgt_orient)
            // We need to track orientation as we traverse
            struct edge_info {
                u64 src;
                u64 tgt;
                int delta;
                bool src_orient;
                bool tgt_orient;
            };
            gtl::vector<edge_info> edges;

            // Use BFS to determine orientations for all nodes
            gtl::flat_hash_map<u64, bool> orientations;
            std::queue<u64> orient_queue;
            orient_queue.push(component[0]);
            orientations[component[0]] = true;  // First node is forward

            while (!orient_queue.empty()) {
                u64 src = orient_queue.front();
                orient_queue.pop();
                bool s_orient = orientations[src];

                auto range = recig.equal_range(src);
                for (auto it = range.first; it != range.second; ++it) {
                    u64 tgt = it->second.node;
                    if (!idx.contains(tgt)) continue;  // Not in this component

                    bool tgt_orient;
                    if (orientations.contains(tgt)) {
                        tgt_orient = orientations[tgt];
                    } else {
                        tgt_orient = !s_orient;  // Alternate orientation
                        orientations[tgt] = tgt_orient;
                        orient_queue.push(tgt);
                    }

                    // Compute delta based on orientations
                    long long qs, ts;
                    if (s_orient) {
                        qs = it->second.qry.size - it->second.qry.to;
                    } else {
                        qs = it->second.qry.from;
                    }
                    if (tgt_orient) {
                        ts = it->second.tgt.size - it->second.tgt.to;
                    } else {
                        ts = it->second.tgt.from;
                    }

                    int delta = static_cast<int>(qs - ts);
                    edges.push_back({src, tgt, delta, s_orient, tgt_orient});
                }
            }

            // Build Laplacian matrix L and vector b
            // L is the graph Laplacian: for each edge (q, t, d):
            //   L[i,i] += 1, L[j,j] += 1, L[i,j] -= 1, L[j,i] -= 1
            //   b[j] += d, b[i] -= d
            Eigen::MatrixXd L = Eigen::MatrixXd::Zero(n, n);
            Eigen::VectorXd b = Eigen::VectorXd::Zero(n);

            for (const auto &e : edges) {
                size_t i = idx[e.src];
                size_t j = idx[e.tgt];
                L(i, i) += 1.0;
                L(j, j) += 1.0;
                L(i, j) -= 1.0;
                L(j, i) -= 1.0;
                b(j) += e.delta;
                b(i) -= e.delta;
            }

            // Pin the first node at position 0
            L.row(0).setZero();
            L(0, 0) = 1.0;
            b(0) = 0.0;

            // Solve L * x = b
            Eigen::VectorXd x_vec = L.colPivHouseholderQr().solve(b);

            // Extract positions
            gtl::flat_hash_map<u64, double> x_positions;
            for (size_t i = 0; i < n; ++i) {
                x_positions[component[i]] = x_vec(i);
            }

            // Sort component by computed positions
            std::sort(component.begin(), component.end(),
                      [&x_positions](u64 a, u64 b) {
                          return x_positions[a] < x_positions[b];
                      });

            // Build the layout map
            layouts.emplace_back();
            double min_pos = x_positions[component[0]];

            for (u64 node_id : component) {
                double pos = x_positions[node_id] - min_pos;
                int start = static_cast<int>(std::round(pos));
                int len = infer_len_for_node(node_id);
                int end = to_end_inclusive(start, len);
                bool orient = orientations.contains(node_id)
                                  ? orientations[node_id]
                                  : true;
                layouts.back().emplace(node_id,
                                       std::make_tuple(start, end, orient));
            }
        }
        return layouts;
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
        if (lr.cores.size() < LCMER_L + 1) {
            // good_n.fetch_sub(1);
            free_lps(&str);
            return;
        }
        gtl::vector<u64> ll = {0};
        lcmer<LCMER_L> l{lr.cores.begin() + 1};  // Skip first core;
        ll.push_back(dbgi.emplace_nobloomfilter(l, idx, tid));
        for (size_t i = 2; i < lr.cores.size() - LCMER_L - 1; ++i) {
            uint32_t in_core = lr.cores[i + LCMER_L];
            l.update_lcmer(in_core);
            ll.push_back(dbgi.emplace_nobloomfilter(l, idx, tid));
        }
        ll.push_back(0);

        for (size_t i = 1; i < ll.size() - 1; ++i) {
            if (insert_lcmer_helper(ll, i, dbgi)) {
                dbgi.dbg.modify_if(
                    ll[i],
                    [&lr, i, &subseqs](dbg_index::dbg_map_t::value_type &val) {
                        val.second.seq.emplace_back(
                            subseqs, lr.pos[i + LCMER_L - 1].second,
                            lr.pos[i + LCMER_L].second);
                    });
            }
        }
    }

    static auto build_lr_cores(const gtl::vector<SequenceContainerType> &seqs,
                               const SyncmerType &syncmer_functor,
                               int lcp_level) -> lcp_read {
        lps str = {0, 0, 0, 0};
        lcp_read lr;
        lr.clear();

        init_lps_crumbled(&str, seqs, syncmer_functor);
        lps_deepen(&str, lcp_level);

        for (int i = 0; i < str.size; ++i) {
            lr.cores.emplace_back(str.cores[i].label);
            lr.pos.emplace_back(str.cores[i].start, str.cores[i].end);
        }
        // if (lr.cores.size() < LCMER_L + 2) {
        //     // free_lps(&str);
        //     return {};
        // }
        // free_lps(&str);
        return lr;
    }

    static auto build_lr_cores_with_seqs(
        const gtl::vector<SequenceContainerType> &seqs,
        const SyncmerType &syncmer_functor, int lcp_level)
        -> std::pair<lcp_read, gtl::vector<SequenceContainerType>> {
        lps str = {0, 0, 0, 0};
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
        const SyncmerType &syncmer_functor, int lcp_level) -> lcp_read {
        lps str = {0, 0, 0, 0};
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
        // if (lr.cores.size() < LCMER_L + 2) {
        //     // free_lps(&str);
        //     return {};
        // }
        // free_lps(&str);
        return lr;
    }

    static auto build_rev_lr_cores_with_seqs(
        const gtl::vector<SequenceContainerType> &seqs,
        const SyncmerType &syncmer_functor, int lcp_level)
        -> std::pair<lcp_read, gtl::vector<SequenceContainerType>> {
        lps str = {0, 0, 0, 0};
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
        if (lr.cores.size() < LCMER_L + 2) {
            free_lps(&str);
            return {};
        }
        lcmer<LCMER_L> l{lr.cores.begin() + 1};  // Skip first core;
        gtl::vector<u64> ll;

        u64 lval = dbgi.emplace_nobloomfilter(l, idx, tid);
        ll.emplace_back(lval);
        for (size_t i = 2; i < lr.cores.size() - LCMER_L - 1; ++i) {
            uint32_t in_core = lr.cores[i + LCMER_L];
            l.update_lcmer(in_core);
            lval = dbgi.emplace_nobloomfilter(l, idx, tid);
            ll.emplace_back(lval);
        }
        gtl::vector<SequenceContainerType> seqs;
        for (size_t i = 1; i < ll.size() - 1; ++i) {
            seqs.emplace_back(rev_seqs, lr.pos[i + LCMER_L - 1].second,
                              lr.pos[i + LCMER_L].second);
            seqs.back().origin_node = ll[i];
        }
        // TODO: index it if needed
        return seqs;
    }
    auto fold_dbg_graph(dbg_index &dbgi) {
        fprintf(stderr, "dbg_size: %zu\n", dbgi.dbg.size());
        fprintf(stderr, "lci_size: %zu\n", dbgi.lci.size());

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
                if (lr.cores.size() < LCMER_L + 1) {
                    bad.fetch_add(1);
                    // good_n.fetch_sub(1);
                    free_lps(&str);
                    return;
                }
                // Keep everything in CORE-INDEX space for layouting.
                // `find_lcmer_ranges()` already returns [from_core,to_core] in
                // core indices (inclusive), which is the correct coordinate
                // system for connected_comp_plus_layout().
                auto found_ranges = find_lcmer_ranges(lr, dbgi);
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

                        matched[node.first].emplace_back(
                            id,
                            rev_match::ival{static_cast<u64>(from_core),
                                            static_cast<u64>(to_core),
                                            core_len_steps},
                            rev_match::ival{0, 0, 0}, from_lcmer, to_lcmer);
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

                        // printf("%llu\t%llu\t%llu\t%llu\t%llu\t%llu\n", k,
                        //        match.node, match.from, match.to,
                        //        ff->from, ff->to);
                        // printf(
                        //     "C%llu\t%llu\t%llu\t%llu\t+\tC%llu\t%llu\t%"
                        //     "llu\t%llu\t%llu\t%llu\t%d\n",
                        //     k, match.size + LCMER_L, match.from,
                        //     match.to + LCMER_L, match.node, ff->size +
                        //     LCMER_L, ff->from, ff->to + LCMER_L, ff->to -
                        //     ff->from, ff->to - ff->from, 250);
                    } else {  // Handle after reciprocals are good
                        // printf(
                        //     "C%llu\t%llu\t%llu\t%llu\t+\tC%llu\t%d\t%d\t%d\t%"
                        //     "d\t%d\t%d\n",
                        //     k, match.size + LCMER_L, match.from, match.to,
                        //     match.node, -1, -1, -1, -1, -1, 250);
                    }
                    // }
                }
                // printf("%llu\t%llu\t%llu\t%llu\n", k, match.node, match.from,
                //        match.to);
            }
        }
        fprintf(stderr, "all:%d\ngood: %d nodes\nhits: %d\n", all_n.load(),
                all_n.load() - bad.load(), hits.load());
        // auto layouts = connected_comp_plus_layout(reciprocal);
        auto layouts = connected_comp_plus_layout_lstsq(reciprocal);
#if 1
        int cc = 0;
        for (auto &l : layouts) {
            fprintf(stdout, "# Layout %d:\n", cc++);
            for (auto &[k, v] : l) {
                auto it = dbgi.dbg.find(k);
                if (it != dbgi.dbg.end()) {
                    auto [start_core, end_core, orient] = v;
                    (void)end_core;
                    fprintf(stdout, "C%llu\t%d\t%d\t%c\n", k, start_core,
                            start_core + end_core, (orient ? '+' : '-'));
                } else {
                    auto [start_core, end_core, orient] = v;
                    (void)end_core;
                    (void)orient;
                    fprintf(stdout, "C%llu\t%d\t(not found)\n", k, start_core);
                }
            }
        }
#endif
        return std::make_pair(reciprocal, layouts);
    }

    struct LayoutEntry {
        u64 id;
        int start;
        int end;
        dbg_index::dbg_map_t::iterator dbg_it;
        bool orient;
    };

    static gtl::vector<u64> path_with_good_ovlps_unweighted_bfs(
        const gtl::btree_multimap<u64, rev_match> &recig, u64 source,
        u64 target) {
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
            std::sort(rv.begin(), rv.end(), [&](const auto &a, const auto &b) {
                const rev_match &ma = a->second;
                const rev_match &mb = b->second;
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
    static gtl::vector<u64> shortest_path_unweighted_bfs(
        const gtl::btree_multimap<u64, rev_match> &recig, u64 source,
        u64 target) {
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
    // static auto layout_shortest_path(
    //     const gtl::btree_multimap<u64, rev_match> &reciprocal,
    //     gtl::vector<LayoutEntry> &entries){ std::sort(entries.begin(),
    //     entries.end(),
    //               [](const auto &a, const auto &b) { return a.start <
    //               b.start; });
    //     u64 source = entries.front().id;
    //     u64 target = entries.back().id;
    //     auto find_shortest_path = [&reciprocal](u64 source, u64 target) {
    //         // Unweighted shortest path (fewest edges) via BFS.
    //         // Track parents to avoid copying full paths for every queued
    //         node. std::queue<u64> q; gtl::flat_hash_set<u64> visited;
    //         gtl::flat_hash_map<u64, u64> parent;

    //         q.push(source);
    //         visited.insert(source);

    //         bool found = false;
    //         while (!q.empty() && !found) {
    //             u64 current = q.front();
    //             q.pop();

    //             if (current == target) {
    //                 found = true;
    //                 break;
    //             }

    //             auto range = reciprocal.equal_range(current);
    //             for (auto it = range.first; it != range.second; ++it) {
    //                 u64 neighbor = it->second.node;
    //                 if (visited.find(neighbor) != visited.end()) continue;

    //                 visited.insert(neighbor);
    //                 parent[neighbor] = current;

    //                 if (neighbor == target) {
    //                     found = true;
    //                     break;
    //                 }

    //                 q.push(neighbor);
    //             }
    //         }

    //         if (!found) return gtl::vector<u64>{};

    //         gtl::vector<u64> path;
    //         for (u64 v = target;; v = parent[v]) {
    //             path.push_back(v);
    //             if (v == source) break;
    //         }
    //         std::reverse(path.begin(), path.end());
    //         return path;
    //     };
    //     return 0;
    // }

    using half_lcmer = lcmer<LCMER_L / 2>;
    struct ctg_builder {
        gtl::vector<SequenceContainerType> prefix;
        gtl::vector<SequenceContainerType> seqs;
        gtl::vector<u32> cores;

        gtl::flat_hash_map<half_lcmer, size_t> idx;
        bool empty() const noexcept { return seqs.empty(); }
    };
    auto find_anchors_and_extend_the_ctg_builder(
        ctg_builder &base, const gtl::vector<SequenceContainerType> &extension,
        const gtl::vector<SequenceContainerType> &next,
        const std::optional<half_lcmer> &prev_anchor, bool forward)
        -> half_lcmer {
        auto [extension_cores, reextension] = [&]() {
            if (forward) {
                return build_lr_cores_with_seqs(extension, syncmer_functor,
                                                lcp_level);
            } else {
                return build_rev_lr_cores_with_seqs(extension, syncmer_functor,
                                                    lcp_level);
            }
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
            } else {
                assert(base.empty());
                return half_lcmer{extension_cores.cores.begin()};
            }
        }();

        auto next_tension_cores = [&]() {
            if (next.empty()) {
                return lcp_read{};
            }
            if (!forward) {
                return build_lr_cores(next, syncmer_functor, lcp_level);
            } else {
                return build_rev_lr_cores(next, syncmer_functor, lcp_level);
            }
        }();
        auto right_anchor = [&]() -> std::optional<half_lcmer> {
            auto &cores = next_tension_cores.cores;
            if (cores.size() < half_lcmer::size) {
                return {};
            }
            half_lcmer l{cores.begin()};
            for (auto it = cores.rbegin() + half_lcmer::size;
                 it != cores.rend(); ++it) {
                if (mini_extension_index.contains(l)) {
                    return l;
                }
                // l.update_lcmer(*it);
                l = half_lcmer::right_shift(l);
                l.setN(*it);
            }
            // Couldn't find.
            return {};
        }();

        if (left_anchor == right_anchor) {
            return right_anchor.value();
        }

        auto left_extension_pos = [&]() {
            auto &cores = extension_cores.cores;
            if (cores.size() < half_lcmer::size) {
                return 0L;
            }
            half_lcmer l{cores.begin()};
            auto it = cores.begin() + half_lcmer::size;
            for (; it != cores.end(); ++it) {
                if (l == left_anchor) {
                    return std::distance(cores.begin(), it);
                }
                l.update_lcmer(*it);
            }
            return 0L;
        }();
        assert(left_extension_pos != 0);

        // Extend until right_anchor
        // Following process_cores pattern: store non-overlapping sequence
        // from end of previous core to end of current core
        auto &cores = extension_cores.cores;

        // Initialize lcmer to the left_anchor window
        assert(left_extension_pos >= half_lcmer::size);
        assert(cores.size() > left_extension_pos);
        half_lcmer l{cores.begin() + (left_extension_pos - half_lcmer::size)};

        // Start adding cores from left_extension_pos
        // Following process_cores pattern: for core at position i, the
        // transition fragment spans from pos[i-1].second to pos[i].second (end
        // of prev core to end of current core). This creates non-overlapping
        // fragments.
        size_t i = left_extension_pos;
        for (; i < cores.size(); ++i) {
            // Update lcmer to include the new core
            l.update_lcmer(cores[i]);

            // Extract the transition fragment using position-based range,
            // matching process_cores semantics: from end of previous core to
            // end of current core
            base.seqs.emplace_back(std::move(reextension[i]));
            base.cores.push_back(cores[i]);
            base.idx.emplace(l, base.seqs.size() - 1);

            if (right_anchor.has_value() && l == right_anchor.value()) {
                break;
            }
        }
        // assert(i < cores.size() && "right_anchor not found in extension");
        return (right_anchor.has_value() ? right_anchor.value() : half_lcmer{});
        // Algorithm
        // Given base, extension and next
        // Find anchor of extension on base, find anchor of next on extension,
        // extend base until 2nd anchor. Accept previous anchors.
    }

    void merge_layouts_ultra(
        FILE *out,
        const gtl::vector<gtl::flat_hash_map<u64, std::tuple<int, int, bool>>>
            &layouts,
        const gtl::btree_multimap<u64, rev_match> &reciprocal,
        dbg_index &dbgi) {
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
        u64 lcmer_index = 0;
        auto single_node_print = [&out, &layout_no, &dbg,
                                  &chs](decltype(*layouts.begin()) &layout) {
            fprintf(out, ">M%llu %zu\n", layout_no, layout.size());
            u64 utig = layout.begin()->first;
            auto it = dbg.find(utig);
            if (it != dbg.end()) {
                for (const auto &sc : it->second.seq) {
                    chs.require(sc.len);
                    sc.convert_to_ascii((char *)chs.data);
                    fprintf(out, "%s", chs.data);
                }
            }
            fprintf(out, "\n");
        };

        auto iterate_until_lcmer_is_found =
            [&](const lcp_read &lr, const lcmer<LCMER_L> &needle) -> long long {
            // Search entire sequence thoroughly - no skipping
            if (lr.cores.size() < LCMER_L) {
                return -2;
            }
            lcmer<LCMER_L> l{lr.cores.begin()};
            if (needle == l) {
                return 0;
            }
            // Search all positions from 1 to the last valid lcmer position
            for (long long i = 1;
                 i <= static_cast<long long>(lr.cores.size() - LCMER_L); ++i) {
                uint32_t in_core = lr.cores[i + LCMER_L - 1];
                l.update_lcmer(in_core);
                if (needle == l) {
                    return i;
                }
            }
            return -1;
        };

        auto find_rev_match =
            [](const gtl::btree_multimap<u64, rev_match> &reciprocal, u64 left,
               u64 right) -> rev_match {
            auto range = reciprocal.equal_range(left);
            for (auto it = range.first; it != range.second; ++it) {
                if (it->second.node == right) {
                    return it->second;
                }
            }
            throw std::runtime_error("No rev match found");
        };

        for (const auto &layout : layouts) {
            if (layout.size() <= 2) {
                // Trivial layout. For now lets just print them as is. We
                // can try to recover some edge sequence later.
                single_node_print(layout);
            } else {
                // printf("LLLL: %llu\n", layout_no);
                // Find leftmost start and rightmost end on the layout.
                gtl::vector<gtl::vector<SequenceContainerType>>
                    merged_seqs;  // per entry we get hte sequences. Even
                                  // forward

                auto [source, target, oriori] =
                    identify_component_endpoints(layout);
                auto path =
                    shortest_path_unweighted_bfs(reciprocal, source, target);
                assert(path.size() > 0);
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

                    rev_match left_right_match =
                        find_rev_match(reciprocal, left, right);

                    auto left_dbg_it = dbg.find(left);
                    auto right_dbg_it = dbg.find(right);
                    assert(left_dbg_it != dbg.end());
                    prev_anchor = find_anchors_and_extend_the_ctg_builder(
                        base, left_dbg_it->second.seq, right_dbg_it->second.seq,
                        prev_anchor, oriori);
                    oriori = !oriori;
                    fprintf(stderr, "Extending: %llu with %llu\n", layout_no,
                            path[i]);
                    // auto l_lr = build_lr_cores(left_dbg_it->second.seq,
                    //                            syncmer_functor, lcp_level);
                    // auto r_lr = build_lr_cores(dbg.find(right)->second.seq,
                    //                            syncmer_functor, lcp_level);

                    // auto l_rev_lr = build_rev_lr_cores(
                    //     dbg.find(left)->second.seq, syncmer_functor,
                    //     lcp_level);

                    // auto r_rev_lr =
                    //     build_rev_lr_cores(dbg.find(right)->second.seq,
                    //                        syncmer_functor, lcp_level);

                    // // Helper to search all sequences as fallback
                    // auto search_with_fallback =
                    //     [&](const lcmer<LCMER_L> &needle,
                    //         const lcp_read &primary,
                    //         const char *primary_name) -> long long {
                    //     long long pos =
                    //         iterate_until_lcmer_is_found(primary, needle);
                    //     if (pos >= 0) {
                    //         printf("F: %s at %lld\n", primary_name, pos);
                    //         return pos;
                    //     }
                    //     printf("NF %s: ", primary_name);
                    //     for (u32 core : needle.data) {
                    //         printf("%u ", core);
                    //     }
                    //     printf("\n");
                    //     // Not found in primary, try all other sequences
                    //     // (silently)

                    //     if (&primary != &l_lr) {
                    //         pos = iterate_until_lcmer_is_found(l_lr, needle);
                    //         if (pos >= 0) {
                    //             printf("F: l_lr at %lld (expected %s)\n",
                    //             pos,
                    //                    primary_name);
                    //         }
                    //     }

                    //     if (&primary != &r_lr) {
                    //         pos = iterate_until_lcmer_is_found(r_lr, needle);
                    //         if (pos >= 0) {
                    //             printf("F: r_lr at %lld (expected %s)\n",
                    //             pos,
                    //                    primary_name);
                    //         }
                    //     }

                    //     if (&primary != &l_rev_lr) {
                    //         pos =
                    //             iterate_until_lcmer_is_found(l_rev_lr,
                    //             needle);
                    //         if (pos >= 0) {
                    //             printf("F: l_rev_lr at %lld (expected %s)\n",
                    //                    pos, primary_name);
                    //         }
                    //     }

                    //     if (&primary != &r_rev_lr) {
                    //         pos =
                    //             iterate_until_lcmer_is_found(r_rev_lr,
                    //             needle);
                    //         if (pos >= 0) {
                    //             printf("F: r_rev_lr at %lld (expected %s)\n",
                    //                    pos, primary_name);
                    //         }
                    //     }

                    //     // printf("NF: %s - Needle: ", primary_name);
                    //     // for (size_t i = 0; i < LCMER_L; ++i) {
                    //     //     printf("%u ", needle.data[i]);
                    //     // }
                    //     // printf("\n");
                    //     return -1;
                    // };

                    // long long needle_posf = search_with_fallback(
                    //     left_right_match.first_matching_lcmer_target, l_lr,
                    //     "l_lr");
                    // long long needle_pos = search_with_fallback(
                    //     left_right_match.last_matching_lcmer_target, l_lr,
                    //     "l_lr");
                    // long long rneedle_posf = search_with_fallback(
                    //     left_right_match.first_matching_lcmer_query, r_lr,
                    //     "r_lr");
                    // long long rneedle_pos = search_with_fallback(
                    //     left_right_match.last_matching_lcmer_query, r_lr,
                    //     "r_lr");

                    // long long rev_needle_posf = search_with_fallback(
                    //     left_right_match.first_matching_lcmer_query,
                    //     l_rev_lr, "l_rev_lr");
                    // long long rev_needle_pos = search_with_fallback(
                    //     left_right_match.last_matching_lcmer_query, l_rev_lr,
                    //     "l_rev_lr");
                    // long long rrev_needle_posf = search_with_fallback(
                    //     left_right_match.first_matching_lcmer_target,
                    //     r_rev_lr, "r_rev_lr");
                    // long long rrev_needle_pos = search_with_fallback(
                    //     left_right_match.last_matching_lcmer_target,
                    //     r_rev_lr, "r_rev_lr");

                    // // If we see any -1, lets print first LCMER_L cores
                    // if (rev_needle_pos == -1) {
                    //     for (u32 c : l_rev_lr.cores) {
                    //         printf("%u ", c);
                    //     }
                    //     printf("\n");
                    //     for (u32 c : r_lr.cores) {
                    //         printf("%u ", c);
                    //     }
                    //     printf("\n");
                    //     for (u32 c : r_rev_lr.cores) {
                    //         printf("%u ", c);
                    //     }
                    //     printf("\n");
                    //     for (u32 c : l_lr.cores) {
                    //         printf("%u ", c);
                    //     }
                    //     printf("\n");
                    //     // for(int i = 0; i < 2*LCMER_L; ++i){
                    //     //     printf("%u ", l_lr.cores[i]);

                    //     // }
                    //     // printf("\n");
                    //     // for(int i = 0; i < 2*LCMER_L; ++i){
                    //     //     printf("%u ", r_rev_lr.cores[i]);

                    //     // }
                    //     // printf("\n!!!\n");
                    // }

                    // // print rev_match
                    // printf(
                    //     "#RM %llu -> %llu : qry[%llu,%llu] tgt[%llu,%llu] "
                    //     "qlen:%llu tlen:%llu ",
                    //     left, right, left_right_match.qry.from,
                    //     left_right_match.qry.to, left_right_match.tgt.from,
                    //     left_right_match.tgt.to, left_right_match.qry.size,
                    //     left_right_match.tgt.size);
                    // printf("\t\t%lld:%lld %lld:%lld -- %lld:%lld
                    // %lld:%lld\n",
                    //        needle_posf, needle_pos, rneedle_posf,
                    //        rneedle_pos, rev_needle_posf, rev_needle_pos,
                    //        rrev_needle_posf, rrev_needle_pos);
                }
                prev_anchor = find_anchors_and_extend_the_ctg_builder(
                    base, dbg.find(path.back())->second.seq, {}, prev_anchor,
                    (path.size() - 1) % 2 == 0);
                fprintf(out, ">L%llu\n", layout_no);
                for (const SequenceContainerType &s : base.seqs) {
                    chs.require(s.len);
                    s.convert_to_ascii((char *)chs.data);
                    fprintf(out, "%s", chs.data);
                }
                fprintf(out, "\n");
                for (u64 node : path) {
                    auto it = dbg.find(node);
                    fprintf(out, ">Z%llu L%llu\n", node, layout_no);
                    for (const SequenceContainerType &s : it->second.seq) {
                        chs.require(s.len);
                        s.convert_to_ascii((char *)chs.data);
                        fprintf(out, "%s", chs.data);
                    }
                    fprintf(out, "\n");
                }
            }
            ++layout_no;
        }
    }
    void merge_layouts_ldbg(
        FILE *out,
        const gtl::vector<gtl::flat_hash_map<u64, std::tuple<int, int, bool>>>
            &layouts,
        const gtl::btree_multimap<u64, rev_match> &reciprocal,
        dbg_index &dbgi) {
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
        // auto find_rev = [&reciprocal](u64 src, u64 step) {
        //     auto range = reciprocal.equal_range(src);
        //     for (auto it = range.first; it != range.second; ++it) {
        //         if (it->second.node == step) {
        //             return it;
        //         }
        //     }
        //     return reciprocal.end();
        // };
        auto process_substr_cores =
            [](const lcp_read &lr, u64 &lcmer_index, int thread_id,
               dbg_index_light &dbgi,
               const gtl::vector<SequenceContainerType> &subseqs) {
                if (lr.cores.size() < LCMER_L) {
                    return gtl::vector<u64>{};
                }
                lcmer<LCMER_L> l{lr.cores.begin() + 1};  // Skip first core;
                gtl::vector<std::tuple<u64, u64, u64>> llr;

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
                    if (insert_lcmer_helper(ll, i, dbgi)) {
                        // Lets try atomic
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
        // u64 prev_lm_query = 0;
        int tid = 0;

        u64 src_lcm = 0, tgt_lcm = 0;

        charsi chs;
        for (const auto &layout : layouts) {
            dbg_index_light smal_dbgi{0.0, 0, 1UL};
            // Find leftmost start and rightmost end on the layout.
            gtl::vector<gtl::vector<SequenceContainerType>>
                merged_seqs;  // per entry we get hte sequences. Even
                              // forward

            auto [source, target] = identify_component_endpoints(layout);
            auto path =
                shortest_path_unweighted_bfs(reciprocal, source, target);
            bool plus = true;
            source = path[0];
            for (size_t i = 0; i < path.size(); ++i) {
                source = path[i];
                // fprintf(out, ">M%llu%c %llu\n", source, (plus ? '+' :
                // '-'),
                //         layout_no);

                auto source_node = dbg.find(source);
                if (source_node == dbg.end()) {
                    fprintf(stderr, "Source node %llu not found in dbg\n",
                            source);
                    continue;
                }
                if (plus) {
                    auto lrs = build_lr_cores(source_node->second.seq,
                                              syncmer_functor, lcp_level);
                    // process_cores(lrs, smal_dbgi, tid, new_index);
                    auto ll =
                        process_substr_cores(lrs, new_index, tid, smal_dbgi,
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
                            u64 &nod = *it;
                            if (nod != 0) {
                                tgt_lcm = nod;
                                break;
                            }
                        }
                    }
                    // for (const SequenceContainerType &seq :
                    //      source_node->second.seq) {
                    //     chs.require(seq.len);
                    //     seq.convert_to_ascii((char *)chs.data);
                    //     fprintf(out, "%s",
                    //             (const char *)chs.data +
                    //             seq.skip_first_n);
                    // }
                    // fprintf(out, "\n");
                } else {
                    gtl::vector<SequenceContainerType> rev_seqs;
                    for (auto it = source_node->second.seq.rbegin();
                         it != source_node->second.seq.rend(); ++it) {
                        rev_seqs.emplace_back(it->reverse_complement());
                    }
                    auto lrs =
                        build_lr_cores(rev_seqs, syncmer_functor, lcp_level);
                    auto ll = process_substr_cores(lrs, new_index, tid,
                                                   smal_dbgi, rev_seqs);
                    if (i == path.size() - 1) {
                        for (auto it = ll.rbegin(); it != ll.rend(); ++it) {
                            u64 &nod = *it;
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
            if (sp.size() > 0) {
                fprintf(out, ">M%llu %zu\n", layout_no, layout.size());
                for (auto it = sp.rbegin(); it != sp.rend(); ++it) {
                    u64 n = *it;
                    SequenceContainerType &seq = smal_dbgi.dbg.at(n).seq[0];
                    chs.require(seq.len);
                    seq.convert_to_ascii((char *)chs.data);
                    fprintf(out, "%s", chs.data);
                }
                fprintf(out, "\n");
            } else {
                fprintf(out, ">M%llu %zu NP\n", layout_no, layout.size());
                auto &seqs = dbgi.dbg.find(path[0])->second.seq;
                for (auto it = seqs.begin(); it != seqs.end(); ++it) {
                    auto &seq = *it;
                    chs.require(seq.len);
                    seq.convert_to_ascii((char *)chs.data);
                    fprintf(out, "%s", chs.data);
                }
                fprintf(out, "\n");
            }
            // for (size_t i = 1; i < path.size(); ++i) {
            //     fprintf(out, ">M%llu%c %llu\n", source, (plus?'+':
            //     '-'),layout_no); u64 step = path[i]; auto revit =
            //     find_rev(source, step); assert(revit !=
            //     reciprocal.end()); const rev_match &rv = revit->second;

            //     auto source_node = dbg.find(source);
            //     if (source_node == dbg.end()) {
            //         fprintf(stderr, "Source node %llu not found in
            //         dbg\n",
            //                 source);
            //         continue;
            //     }

            //     if (plus) {
            //         for (const SequenceContainerType &seq :
            //              source_node->second.seq) {
            //             chs.require(seq.len);
            //             seq.convert_to_ascii((char *)chs.data);
            //             fprintf(out, "%s",
            //                     (const char *)chs.data +
            //                     seq.skip_first_n);
            //         }
            //         fprintf(out, "\n");
            //     } else {
            //         auto lrs = find_rev_lcmers(source_node->second.seq,
            //         dbgi,
            //                                    13, new_index);
            //         auto it = lrs.begin();

            //         while (it != lrs.end()) {
            //             SequenceContainerType &seq = *it;
            //             chs.require(seq.len);
            //             seq.convert_to_ascii((char *)chs.data);
            //             fprintf(out, "%s",
            //                     (const char *)chs.data +
            //                     seq.skip_first_n);
            //             ++it;
            //         }
            //         fprintf(out, "\n");
            //     }
            //     prev_lm_query = rv.first_matching_lcmer_target;
            //     plus = !plus;
            //     source = step;
            // }

            ++layout_no;
        }
    }

    static void print_merged_layouts(
        FILE *out,
        gtl::vector<gtl::flat_hash_map<u64, std::tuple<int, int, bool>>>
            layouts,
        dbg_index::dbg_map_t &dbg) {
        charsi buf;
        u64 contig_id = 0;

        // Diagnostics: `pos` in `layouts` is produced by layouting
        // reciprocal matches which are typically in "core-space" /
        // index-space, while `seq_len` below is in bases. If these units
        // are mixed, overlap trimming can delete/duplicate large sequence,
        // potentially leading to misassemblies. We don't change behavior
        // here, but we warn.

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

            // In core-space mode:
            // - `pos` is in core units (layout space)
            // - sequence lengths are in bases
            //
            // We keep start/end in core-space here and only convert
            // core->base when actually skipping bases during concatenation
            // (at output time).
            for (auto &[node_id, start_end_orient] : layout) {
                auto [pos_core, end_core, orient] = start_end_orient;
                auto it = dbg.find(node_id);
                if (it == dbg.end() || it->second.seq.empty()) continue;

                entries.push_back({pos_core, end_core, node_id, orient, it});
            }

            if (entries.empty()) continue;

            std::sort(
                entries.begin(), entries.end(),
                [](const auto &a, const auto &b) { return a.start < b.start; });

            // Compute a core->base scaling factor for this component, used
            // only at emission time. We estimate: bases_per_core_step ≈
            // (total_bases_across_entries / core_span_across_entries). This
            // avoids mixing units during layouting/ordering.
            double core_to_base_scale = 0.0;
            {
                int min_start = INT_MAX;
                int max_start = INT_MIN;
                size_t total_bases = 0;
                for (const auto &e : entries) {
                    min_start = std::min(min_start, e.start);
                    max_start = std::max(max_start, e.start);
                    // base length of this node
                    size_t seq_len = 0;
                    for (const auto &s : e.dbg_it->second.seq) {
                        seq_len += s.len;
                    }
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

            fprintf(out, "S\tM%" PRIu64 "\t", contig_id);

            size_t total_len = 0;
            int current_end = non_contained[0].start;

            gtl::vector<u64> prev_edges = non_contained[0].dbg_it->second.prev;

            for (size_t i = 0; i < non_contained.size(); ++i) {
                auto &curr = non_contained[i];

                // Overlap is computed in CORE units.
                int overlap_cores = current_end - curr.start;
                if (overlap_cores < 0) overlap_cores = 0;

                // Convert core overlap into an approximate base skip at
                // emission time. We use a per-component scaling based on
                // total base length / total core span. This keeps all
                // layout reasoning in core space until the final emission.
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
                        const char *buf_at =
                            s.convert_to_ascii((char *)buf.data);
                        fprintf(out, "%s", buf_at);
                        total_len += s.len;
                    } else if (skipped + s.len > bases_to_skip) {
                        // Partial skip: print only the suffix beyond
                        // bases_to_skip. NOTE: The previous code printed
                        // the entire segment, duplicating overlaps and
                        // inflating assemblies. This fix assumes ASCII
                        // encoding with 1 char/base.
                        const size_t local_skip = bases_to_skip - skipped;
                        buf.require(s.len + 32);
                        const char *buf_at =
                            s.convert_to_ascii((char *)buf.data);
                        fprintf(out, "%s", buf_at + local_skip);
                        total_len += (s.len - local_skip);
                    }
                    skipped += s.len;
                }

                current_end = std::max(current_end, curr.end);
            }

            fprintf(out, "\tLN:i:%zu\n", total_len);

            // Print incoming edges from first contig
            for (u64 prev_id : prev_edges) {
                if (dbg.contains(prev_id)) {
                    fprintf(out, "L\tC%" PRIu64 "\t+\tM%" PRIu64 "\t+\t0M\n",
                            prev_id, contig_id);
                }
            }

            contig_id++;
        }
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

int run_assembly(const asm_config &cfg) {
    if (!cfg.silent) cfg.print_cfg(stderr);
    // free(malloc(1024ULL * 1024 * 1024 * 20)); INTIMIDATE the OS for the
    // ram

    INIT_TIMERS(HIFI_ASSEMBLY);

    // auto inty = identify_input(cfg.input);

    // if (inty == input_type::UNKNOWN) {
    //     fprintf(stderr, "Unknown input type %s\n", cfg.input.c_str());
    //     return -1;
    // }

    assembler as(cfg);
    dbg_index *dbgp = as.construct_dbg();
    dbg_index &dbg = *dbgp;
    IF_TIME_CHECKPOINT(!cfg.silent, HIFI_ASSEMBLY, stderr, "DBG construction");

    dbg.count_forward_edges();
    IF_TIME_CHECKPOINT(!cfg.silent, HIFI_ASSEMBLY, stderr,
                       "Count forward edges");

    if (cfg.print_precleanup_graph) {
        FILE *out = fopen((cfg.output + ".raw").c_str(), "w");
        print_dbg_seqs(out, dbg);
        fclose(out);

        IF_TIME_CHECKPOINT(!cfg.silent, HIFI_ASSEMBLY, stderr, "DBG printing");
    }

    GraphCleaning::run(dbg, cfg);

    IF_TIME_CHECKPOINT(!cfg.silent, HIFI_ASSEMBLY, stderr, "Graph cleaning");

    auto [recig, layouts] = as.fold_dbg_graph(dbg);
    IF_TIME_CHECKPOINT(!cfg.silent, HIFI_ASSEMBLY, stderr, "Fold graph");

    // GraphCleaning::run(dbg, cfg, );

    FILE *out = fopen(cfg.output.c_str(), "w");

    // assembler::print_merged_layouts(out, layouts, dbg.dbg);
    as.merge_layouts_ultra(out, layouts, recig, dbg);
    // assembler::print_merged_layouts_weave_majority_orientation(out,
    // recig, std::move(layouts), dbg.dbg);
    //  Print remaining unmerged nodes (their sequences weren't cleared)
    // print_dbg_seqs(out, dbg);
    fclose(out);
    IF_TIME_CHECKPOINT(!cfg.silent, HIFI_ASSEMBLY, stderr, "DBG printing");

    return 0;
}

void tests_mock_fill_seq(dbg_index &dbgi) {
    for (auto &node : dbgi.dbg) {
        node.second.seq.emplace_back();
    }
}

// TEST_CASE("Graph simplification simple") {
//     dbg_index dbgi{0.95, 10, 0};
//     gtl::vector<u32> mock_cores;
//     int MOCK_CORE_COUNT = 40;
//     for (int i = 0; i < MOCK_CORE_COUNT; i++) {
//         mock_cores.push_back(MurmurHash3_32(&i, sizeof(int), 156));
//     }
//     lcmer<LCMER_L> l(mock_cores.begin());
//     u64 idx = 0;

//     gtl::vector<u64> ll(1, 0);
//     ll.reserve(MOCK_CORE_COUNT - LCMER_L + 3);
//     ll.push_back(dbgi.emplace_nobloomfilter(l, idx, 0));
//     for (int i = LCMER_L + 1; i < MOCK_CORE_COUNT - 1; ++i) {
//         uint32_t in_core = mock_cores[i];
//         l.update_lcmer(in_core);
//         ll.push_back(dbgi.emplace_nobloomfilter(l, idx, 0));
//     }
//     ll.push_back(0);

//     for (size_t i = 1; i <= ll.size() - 1; i++) {
//         assembler::insert_lcmer_helper(ll, i, dbgi);
//     }

//     SUBCASE("Check ll vec") {
//         for (size_t i = 1; i < ll.size() - 1; ++i) {
//             CHECK(ll[i] == (i - 1) * 128);
//         }
//     }
//     SUBCASE("Pre simplify dbg") {
//         for (auto &node : dbgi.dbg) {
//             INFO(node.first, "\t.\t.\t", node.second.seq.size());
//             CHECK(node.second.core_count < 2);
//             CHECK(node.second.prev.size() < 2);
//         }
//         CAPTURE(dbgi.dbg.size());
//         CHECK(dbgi.lci.size() == MOCK_CORE_COUNT - 2 - LCMER_L + 1);
//     }

//     tests_mock_fill_seq(dbgi);
//     for (auto &n : dbgi.dbg) {  // Main Path has 10X depth
//         n.second.core_count *= 10;
//     }
//     dbgi.count_forward_edges();

//     GraphCleaning::run(dbgi, {.graph_cleaning = {1, 10, 10, 2}});
//     SUBCASE("Post simplify dbg") {
//         for (auto &node : dbgi.dbg) {
//             INFO(node.first, "\t.\t.\t", node.second.seq.size());
//         }
//         CAPTURE(dbgi.dbg.size());
//         CHECK(dbgi.dbg.size() == 1);
//         CHECK(dbgi.dbg.begin()->second.seq.size() == 30);
//     }

//     // // Write Dot file
//     // FILE *dot_file = fopen("graph.dot", "w");

//     // if (dot_file) {
//     //     // print_raw_dbg(dot_file, dbgi);
//     //     fprintf(dot_file, "digraph G {\n");
//     //     for (auto &node : dbgi.dbg) {
//     //         fprintf(dot_file, "\t%" PRIu64 " [label=\"%" PRIu64
//     "\"];\n",
//     //         node.first,
//     //                 node.first);
//     //         for (auto &prev : node.second.prev) {
//     //             fprintf(dot_file, "\t%" PRIu64 " -> %" PRIu64 ";\n",
//     prev,
//     //             node.first);
//     //         }
//     //     }
//     //     fprintf(dot_file, "}\n");
//     //     fclose(dot_file);
//     // }
// }

// TEST_CASE("Graph simplification Bubble") {
//     dbg_index dbgi{0.95, 10, 0};
//     gtl::vector<u32> mock_cores;
//     int MOCK_CORE_COUNT = 40;
//     for (int i = 0; i < MOCK_CORE_COUNT; i++) {
//         mock_cores.push_back(MurmurHash3_32(&i, sizeof(int), 156));
//     }
//     lcmer<LCMER_L> l(mock_cores.begin());
//     u64 idx = 0;

//     gtl::vector<u64> ll(1, 0);
//     ll.reserve(MOCK_CORE_COUNT - LCMER_L + 3);
//     ll.push_back(dbgi.emplace_nobloomfilter(l, idx, 0));
//     for (int i = LCMER_L; i < MOCK_CORE_COUNT - 1; ++i) {
//         uint32_t in_core = mock_cores[i];
//         l.update_lcmer(in_core);
//         ll.push_back(dbgi.emplace_nobloomfilter(l, idx, 0));
//     }
//     ll.push_back(0);
//     for (size_t i = 1; i <= ll.size() - 1; i++) {
//         assembler::insert_lcmer_helper(ll, i, dbgi);
//     }

//     ll.clear();
//     mock_cores[MOCK_CORE_COUNT / 2] = 123123123;
//     l = lcmer<LCMER_L>(mock_cores.begin());

//     ll.push_back(dbgi.emplace_nobloomfilter(l, idx, 0));
//     for (int i = LCMER_L; i < MOCK_CORE_COUNT - 1; ++i) {
//         uint32_t in_core = mock_cores[i];
//         l.update_lcmer(in_core);
//         ll.push_back(dbgi.emplace_nobloomfilter(l, idx, 0));
//     }
//     ll.push_back(0);
//     for (size_t i = 1; i <= ll.size() - 1; i++) {
//         assembler::insert_lcmer_helper(ll, i, dbgi);
//     }

//     dbgi.count_forward_edges();
//     for (auto &n : dbgi.dbg) {  // Main Path has 10X depth
//         n.second.core_count *= 10;
//     }
//     tests_mock_fill_seq(dbgi);
//     print_dbg_dot(dbgi, 0);

//     GraphCleaning::run(dbgi, {.graph_cleaning = {1, 10, 10, 2}});
//     print_dbg_dot(dbgi, 1);

//     SUBCASE("Post simplify dbg iter 2") {
//         for (auto &node : dbgi.dbg) {
//             INFO(node.first, "\t.\t.\t", node.second.seq.size());
//         }
//         CAPTURE(dbgi.dbg.size());
//         CHECK(dbgi.dbg.size() == 1);
//     }

//     print_dbg_dot(dbgi, 2);
// }

// TEST_CASE("Graph simplification Tips (Outgoing)") {
//     dbg_index dbgi{0.95, 10, 0};
//     gtl::vector<u32> mock_cores;
//     int MOCK_CORE_COUNT = 30;
//     for (int i = 0; i < MOCK_CORE_COUNT; i++) {
//         mock_cores.push_back(MurmurHash3_32(&i, sizeof(int), 156));
//     }
//     lcmer<LCMER_L> l(mock_cores.begin());
//     u64 idx = 0;

//     gtl::vector<u64> ll(1, 0);
//     ll.reserve(MOCK_CORE_COUNT - LCMER_L + 3);
//     ll.push_back(dbgi.emplace_nobloomfilter(l, idx, 0));
//     for (int i = LCMER_L; i < MOCK_CORE_COUNT - 1; ++i) {
//         uint32_t in_core = mock_cores[i];
//         l.update_lcmer(in_core);
//         ll.push_back(dbgi.emplace_nobloomfilter(l, idx, 0));
//     }
//     ll.push_back(0);
//     for (size_t i = 1; i <= ll.size() - 1; i++) {
//         assembler::insert_lcmer_helper(ll, i, dbgi);
//     }
//     for (auto &n : dbgi.dbg) {  // Main Path has 10X depth
//         n.second.core_count *= 10;
//     }
//     ll.clear();
//     mock_cores[12] = 123123123;
//     l = lcmer<LCMER_L>(mock_cores.begin());

//     ll.push_back(dbgi.emplace_nobloomfilter(l, idx, 0));
//     for (int i = LCMER_L; i < MOCK_CORE_COUNT - 1; ++i) {
//         uint32_t in_core = mock_cores[i];
//         l.update_lcmer(in_core);
//         ll.push_back(dbgi.emplace_nobloomfilter(l, idx, 0));
//     }
//     ll.push_back(0);
//     for (size_t i = 1; i <= ll.size() / 2 - 1; i++) {
//         assembler::insert_lcmer_helper(ll, i, dbgi);
//     }

//     dbgi.count_forward_edges();

//     tests_mock_fill_seq(dbgi);
//     print_dbg_dot(dbgi, 3);

//     GraphCleaning::run(dbgi, {.graph_cleaning = {1, 10, 10, 2}});
//     SUBCASE("Post simplify dbg") {
//         for (auto &node : dbgi.dbg) {
//             INFO(node.first, "\t.\t.\t", node.second.seq.size());
//         }
//         CAPTURE(dbgi.dbg.size());
//         CHECK(dbgi.dbg.size() == 1);
//     }
//     print_dbg_dot(dbgi, 4);
// }
// TEST_CASE("Graph simplification Tips (Incoming)") {
//     dbg_index dbgi{0.95, 10, 0};
//     gtl::vector<u32> mock_cores;
//     lcmer<LCMER_L> l;
//     gtl::vector<u64> ll;

//     u64 idx = 0;
//     int MOCK_CORE_COUNT = 30;
//     for (int i = 0; i < MOCK_CORE_COUNT; i++) {
//         mock_cores.push_back(MurmurHash3_32(&i, sizeof(int), 156));
//     }
//     for (int DEPTH = 0; DEPTH < 10; ++DEPTH) {
//         l = lcmer<LCMER_L>(mock_cores.begin());

//         ll.push_back(0);

//         ll.reserve(MOCK_CORE_COUNT - LCMER_L + 3);
//         ll.push_back(dbgi.emplace_nobloomfilter(l, idx, 0));
//         for (int i = LCMER_L; i < MOCK_CORE_COUNT - 1; ++i) {
//             uint32_t in_core = mock_cores[i];
//             l.update_lcmer(in_core);
//             ll.push_back(dbgi.emplace_nobloomfilter(l, idx, 0));
//         }
//         ll.push_back(0);
//         for (size_t i = 1; i <= ll.size() - 1; i++) {
//             assembler::insert_lcmer_helper(ll, i, dbgi);
//         }

//         ll.clear();
//     }
//     ll.push_back(0);
//     mock_cores[8] = 123123123;
//     l = lcmer<LCMER_L>(mock_cores.begin() + 7);

//     ll.push_back(dbgi.emplace_nobloomfilter(l, idx, 0));
//     for (int i = LCMER_L + 7; i < MOCK_CORE_COUNT - 1; ++i) {
//         uint32_t in_core = mock_cores[i];
//         l.update_lcmer(in_core);
//         ll.push_back(dbgi.emplace_nobloomfilter(l, idx, 0));
//     }
//     ll.push_back(0);
//     for (size_t i = 1; i <= ll.size() / 2 - 1; i++) {
//         assembler::insert_lcmer_helper(ll, i, dbgi);
//     }

//     dbgi.count_forward_edges();

//     tests_mock_fill_seq(dbgi);
//     print_dbg_dot(dbgi, 6);
//     GraphCleaning::run(dbgi, {.graph_cleaning = {1, 10, 10, 2}});

//     SUBCASE("Post simplify dbg") {
//         for (auto &node : dbgi.dbg) {
//             INFO(node.first, "\t.\t.\t", node.second.seq.size());
//         }
//         CAPTURE(dbgi.dbg.size());
//         CHECK(dbgi.dbg.size() == 1);
//     }
//     print_dbg_dot(dbgi, 7);
// }

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
        auto [lr, core_seqs] = assembler::build_lr_cores_with_seqs(
            seqs, syncmer_functor, lcp_level);
        auto lr_original =
            assembler::build_lr_cores(seqs, syncmer_functor, lcp_level);

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
        auto [lr, core_seqs] = assembler::build_rev_lr_cores_with_seqs(
            seqs, syncmer_functor, lcp_level);
        auto lr_original =
            assembler::build_rev_lr_cores(seqs, syncmer_functor, lcp_level);

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
        auto [lr, core_seqs] = assembler::build_lr_cores_with_seqs(
            seqs, syncmer_functor, lcp_level);

        // Ensure we actually have cores to test
        REQUIRE(lr.cores.size() > 0);

        // Get fragment boundaries
        size_t frag1_end = strlen(seq1);
        size_t frag2_end = frag1_end + strlen(seq2);

        // Find non-overlapping portions that span fragment boundaries
        bool found_spanning_seq = false;
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
        auto [lr, core_seqs] = assembler::build_rev_lr_cores_with_seqs(
            seqs, syncmer_functor, lcp_level);

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
        bool found_spanning_seq = false;
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
