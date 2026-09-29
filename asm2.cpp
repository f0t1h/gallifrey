#include <blockingconcurrentqueue.h>
#include <fcntl.h>
#include <htslib/bgzf.h>
#include <htslib/hfile.h>
#include <htslib/hts.h>
#include <htslib/sam.h>
#include <htslib/thread_pool.h>
#include <stdio.h>
#include <unistd.h>
#if defined(__GLIBC__)
#include <malloc.h>
#endif

#include <algorithm>
#include <chrono>
#include <cstdarg>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <gtl/btree.hpp>
#include <gtl/phmap.hpp>
#include <gtl/phmap_fwd_decl.hpp>
#include <gtl/vector.hpp>
#include <iostream>
#include <memory>
#include <random>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

#include "asm.h"
#include "config.h"
#include "fastx.h"
#include "graph.h"
#include "layout_graph.h"
#include "lcmer.h"
#include "lcptools/lcptools_ho.h"
#include "no_destructor.h"
#include "syncmer.h"

using std::string;
// OLC end-skip policy: cores skipped at each end of a DBG node before
// overlap-anchoring (see find_lcmer_ranges). Function-like macros of the OLC
// lcmer length `l`. A node needs (l + OLC_SKIP_FIRST(l) + OLC_SKIP_LAST(l))
// cores to anchor any overlap.
#ifndef OLC_SKIP_FIRST
#define OLC_SKIP_FIRST(l) (1)
#endif
#ifndef OLC_SKIP_LAST
#define OLC_SKIP_LAST(l) ((l) / 2)
#endif

struct assembler {
    using SYNCMER_HASHER = syncmer::MurMur;
    using SyncmerType = syncmer::GenericSubmerRule<uint32_t, SYNCMER_HASHER>;

    template <class Rec>
    using queue_of = moodycamel::BlockingConcurrentQueue<Rec *>;
    // Records from the readers to the workers; a full queue blocks readers.
    template <class Rec>
    using data_queue_of = bounded_queue<Rec *>;
    using hts_queue = queue_of<bam1_t>;
    using fastx_queue = queue_of<seq_record>;
    // (path, is_compressed) -- the flag is only meaningful for FASTX inputs.
    using stage_input = std::pair<string, bool>;
    gtl::vector<string> inputs;
    string output;
    decltype(asm_config::threads) threads;
    int lcp_level;
    decltype(asm_config::bloom) bloom;
    uint32_t seed;
    decltype(asm_config::syncmer) syncmer_params;
    SyncmerType syncmer_functor;
    int lcmer_l;
    bool parallel_file_reading;
    bool drop_crumbs;
    bool prime_reverse;
    bool strict_layout_labels;
    bool end_joins;

    assembler(const asm_config &cfg)
        : inputs{},
          output{cfg.output},
          threads{cfg.threads},
          lcp_level{cfg.lcp_level},
          bloom{cfg.bloom},
          seed{cfg.seed},
          syncmer_params{cfg.syncmer},
          syncmer_functor{syncmer_params, seed},
          lcmer_l{cfg.lcmer_l},
          parallel_file_reading{cfg.parallel_file_reading},
          drop_crumbs{cfg.drop_crumbs},
          prime_reverse{cfg.prime_reverse},
          strict_layout_labels{cfg.strict_layout_labels},
          end_joins{cfg.end_joins} {
        inputs.assign(cfg.inputs.begin(), cfg.inputs.end());
    }

    // Seeded (--multipass segment windows): the window adds no count and marks
    // its node seeded; edges are made as for reads. A new node gets its first
    // fragment, fragment(key), while the insert holds the submap lock.
    template <class BFT, bool Seeded = false, class Fragment>
    static inline auto insert_lcmer_helper(const gtl::vector<uint64_t> &ll,
                                           const uint32_t &i,
                                           dbg_index_base<BFT> &dbgi,
                                           Fragment &&fragment) -> bool {
        using index_type = dbg_index_base<BFT>;
        using val_type = typename index_type::bidbg_val;
        constexpr unsigned kCount = [] {
            if constexpr (Seeded) {
                return 0U;
            } else {
                return 1U;
            }
        }();
        const auto make = [&](auto... args) {
            val_type v{args...};
            if constexpr (Seeded) v.seeded = true;
            v.seq.push_back(fragment(ll[i]));
            return v;
        };
        const auto mark = [](val_type &v) {
            if constexpr (Seeded) v.seeded = true;
        };
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
                        val.second.counts.core_count += kCount;
                        mark(val.second);
                    },
                    [&](const index_type::dbg_map_t::constructor &ctr) {
                        ctr(ll[i], make(kCount));
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
                        mark(val.second);
                    },
                    [&](const index_type::dbg_map_t::constructor &ctr) {
                        ctr(ll[i], make(0U, ll[i - 1]));
                    });
            } else {
                return dbgi.dbg.lazy_emplace_l(
                    ll[i],
                    [&](index_type::dbg_map_t::value_type &val) {
                        val.second.emplace_prev(ll[i - 1]);
                        val.second.counts.core_count += kCount;
                        mark(val.second);
                    },
                    [&](const index_type::dbg_map_t::constructor &ctr) {
                        ctr(ll[i], make(kCount, ll[i - 1]));
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
        // Cores skipped at each node end before overlap-anchoring (OLC_SKIP_*
        // macros).
        const long skip_first = OLC_SKIP_FIRST(lcmer_l);
        const long skip_last = OLC_SKIP_LAST(lcmer_l);
        // Need enough cores to form at least one valid lcmer after skipping
        if (lr.cores.size() <
            static_cast<size_t>(lcmer_l + skip_first + skip_last)) {
            return gtl::vector<std::tuple<u64, u64, u64, LCmer, LCmer>>{};
        }

        LCmer l{lr.cores.begin() + skip_first, lcmer_l};  // Skip first core;
        gtl::vector<std::tuple<u64, u64, u64, LCmer, LCmer>> llr;

        dbg_index::lci_type lval = dbgi.get(l);
        llr.emplace_back(lval.id, skip_first, skip_first, l, l);
        // Loop bound keeps i + lcmer_l in range.
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

        // Node ll[i] is the window of cores i..i+L-1. It stores the bases its
        // own last core adds, [end(i+L-2), end(i+L-1)): they lie inside the
        // key window, so a successor's fragment starts where this one ends on
        // every path, and spelling a path concatenates fragments exactly.
        for (size_t i = 1; i < ll.size() - 1; ++i) {
            insert_lcmer_helper(ll, i, dbgi, [&](u64 key) {
                return SequenceContainerType(
                    lr.data, lr.pos[i + lcmer_l - 2].second,
                    lr.pos[i + lcmer_l - 1].second, key);
            });
        }
    }

    // --prime-reverse: parses the read's reverse complement from its nibbles
    // and inserts those lcmer keys into the Bloom filter only (no node, count
    // or edge), so the first sighting of the sequence on the other strand is
    // admitted. The windows are the ones process_cores would take. `rc` and
    // `lr` are caller scratch.
    template <class BFT>
    static void prime_reverse_complement(const uint8_t *nibbles, size_t len,
                                         int lcp_level, core_buf *str,
                                         const SyncmerType &syncmer_functor,
                                         charsi &chs, pos_scratch &pos_buf,
                                         gtl::vector<uint8_t> &rc, lcp_read &lr,
                                         dbg_index_base<BFT> &dbgi,
                                         int lcmer_l) {
        revcomp_nibbles(nibbles, static_cast<uint32_t>(len), rc);
        lr.clear();
        parse_cores_from_nibbles(rc.data(), len, lcp_level, str,
                                 syncmer_functor, chs, lr, pos_buf);
        if (lr.cores.size() < static_cast<size_t>(lcmer_l + 2)) {
            return;
        }
        uint64_t unused_index = 0;
        LCmer l{lr.cores.begin() + 1, lcmer_l};  // Skip first core;
        dbgi.emplace(l, unused_index, 0);
        for (size_t i = static_cast<size_t>(lcmer_l + 1);
             i < lr.cores.size() - 1; ++i) {
            l.update_lcmer(lr.cores[i], lr.cores[i - lcmer_l], lcmer_l);
            dbgi.emplace(l, unused_index, 0);
        }
    }

    // ── --multipass: seeding a pass with the previous pass's segments ─────
    // Base j of seq[from, to), or of its reverse complement.
    static uint8_t range_base(const SequenceContainerType &seq, size_t from,
                              size_t to, bool rc, size_t j) {
        using view = KmerViewOnCrumbs64<uint64_t>;
        const size_t skip = static_cast<size_t>(seq.skip_first_n);
        if (rc) {
            return static_cast<uint8_t>(
                3 ^ view::base_at(seq.data, skip + to - 1 - j));
        }
        return view::base_at(seq.data, skip + from + j);
    }

    // Cores of seq[from, to) (or of its reverse complement) at lcp_level,
    // positions in the range's own coordinates. The run-length stream is
    // read straight from the 2-bit words, so the range is never copied.
    static void parse_range(const SequenceContainerType &seq, size_t from,
                            size_t to, bool rc, int lcp_level, core_buf *str,
                            const SyncmerType &syncmer_functor, charsi &rle,
                            pos_scratch &pos, lcp_read &lr) {
        lr.clear();
        const size_t len = to - from;
        rle.require(len + kOnehotRlePad);
        pos.require(len + kOnehotRlePad);
        int64_t rle_len = 0;
        uint8_t prev = 255;
        for (size_t j = 0; j < len; ++j) {
            const uint8_t b = range_base(seq, from, to, rc, j);
            if (b != prev) {
                rle.data[rle_len] = b;
                pos.data[rle_len] = static_cast<uint32_t>(j);
                ++rle_len;
                prev = b;
            }
        }
        parse1_rle_syncmers(rle.data, pos.data, rle_len,
                            static_cast<int64_t>(len), str, syncmer_functor);
        str->deepen(syncmer_functor, lcp_level);
        lr.reserve(str->size);
        str->append_to(lr);
    }

    // Bases [a, b) of the range seq[from, to) in the orientation it was
    // parsed in, as a node fragment.
    static SequenceContainerType range_fragment(
        const SequenceContainerType &seq, size_t from, size_t to, bool rc,
        size_t a, size_t b, u64 origin_node) {
        SequenceContainerType out{b - a, 0};
        out.origin_node = origin_node;
        for (size_t j = a; j < b; ++j) {
            const size_t bit = (j - a) * 2;
            out.data[bit >> 6] |=
                static_cast<uint64_t>(range_base(seq, from, to, rc, j))
                << (bit & 63);
        }
        return out;
    }

    // The windows process_cores takes from a read, seeded: keyed without the
    // Bloom filter, no count, node marked seeded. A new node's fragment comes
    // from the parsed range. Returns the number of windows. `ll` is scratch.
    template <class BFT>
    static size_t process_seeded_cores(const lcp_read &lr,
                                       const SequenceContainerType &seq,
                                       size_t from, size_t to, bool rc,
                                       dbg_index_base<BFT> &dbgi, int lcmer_l,
                                       gtl::vector<uint64_t> &ll) {
        if (lr.cores.size() < static_cast<size_t>(lcmer_l + 2)) return 0;
        uint64_t unused_index = 0;
        LCmer l{lr.cores.begin() + 1, lcmer_l};  // Skip first core;
        ll.assign(1, 0);
        ll.push_back(dbgi.emplace_nobloomfilter(l, unused_index, 0));
        for (size_t i = static_cast<size_t>(lcmer_l + 1);
             i < lr.cores.size() - 1; ++i) {
            l.update_lcmer(lr.cores[i], lr.cores[i - lcmer_l], lcmer_l);
            ll.push_back(dbgi.emplace_nobloomfilter(l, unused_index, 0));
        }
        ll.push_back(0);
        for (size_t i = 1; i < ll.size() - 1; ++i) {
            insert_lcmer_helper<BFT, true>(ll, i, dbgi, [&](u64 key) {
                return range_fragment(seq, from, to, rc,
                                      lr.pos[i + lcmer_l - 2].second,
                                      lr.pos[i + lcmer_l - 1].second, key);
            });
        }
        return ll.size() - 2;
    }

    // Segments are seeded in chunks, so parse buffers stay chunk-sized; the
    // overlap is far above one window plus the cores a parse loses at each
    // end, so every window of a segment lies inside some chunk.
    static constexpr size_t kSeedChunkBases = size_t{1} << 20;
    static constexpr size_t kSeedChunkOverlap = size_t{1} << 16;

    // Seeds dbgi with every window of every segment, in both orientations
    // (the LCP parse is not strand-symmetric). Run after all reads, so read
    // nodes keep their read fragments. Returns {chunks, windows}.
    std::pair<size_t, size_t> insert_segments(
        dbg_index &dbgi, const gtl::vector<SequenceContainerType> &segments,
        size_t chunk_bases = kSeedChunkBases,
        size_t chunk_overlap = kSeedChunkOverlap) const {
        struct chunk {
            size_t seg, from, to;
        };
        gtl::vector<chunk> chunks;
        for (size_t s = 0; s < segments.size(); ++s) {
            const size_t len = segments[s].len;
            for (size_t from = 0; from < len;
                 from += chunk_bases - chunk_overlap) {
                const size_t to = std::min(len, from + chunk_bases);
                chunks.push_back({s, from, to});
                if (to == len) break;
            }
        }
        std::atomic<size_t> next{0};
        std::atomic<size_t> windows{0};
        gtl::vector<std::thread> workers;
        for (int t = 0; t < threads.work; ++t) {
            workers.emplace_back([&]() {
                core_buf str;
                charsi rle;
                pos_scratch pos;
                lcp_read lr;
                gtl::vector<uint64_t> ll;
                size_t n = 0;
                for (size_t c = next++; c < chunks.size(); c = next++) {
                    const auto &[s, from, to] = chunks[c];
                    for (bool rc : {false, true}) {
                        parse_range(segments[s], from, to, rc, lcp_level, &str,
                                    syncmer_functor, rle, pos, lr);
                        n += process_seeded_cores(lr, segments[s], from, to, rc,
                                                  dbgi, lcmer_l, ll);
                    }
                }
                windows += n;
            });
        }
        for (auto &w : workers) w.join();
        return {chunks.size(), windows.load()};
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

    // Record traits. The worker loop only needs the packed nt16 sequence and
    // its base count, so each input path keeps its own record type and the
    // htslib (BAM/CRAM/SAM) path stays copy-free.
    struct hts_record_ops {
        using record = bam1_t;
        static record *create() { return bam_init1(); }
        static void destroy(record *r) { bam_destroy1(r); }
        static const uint8_t *seq(const record *r) { return bam_get_seq(r); }
        static size_t len(const record *r) {
            return static_cast<size_t>(r->core.l_qseq);
        }
    };
    struct fastx_record_ops {
        using record = seq_record;
        static record *create() { return seq_record_init(); }
        static void destroy(record *r) { seq_record_destroy(r); }
        static const uint8_t *seq(const record *r) { return r->nibbles; }
        static size_t len(const record *r) { return r->len; }
    };

    template <class Ops>
    static auto process_data(int thread_id, int lcp_level,
                             data_queue_of<typename Ops::record> &data_queue,
                             queue_of<typename Ops::record> &recycle_queue,
                             dbg_index &dbgi, std::atomic<bool> &done,
                             int &processed_read_count, int &all_reads,
                             auto Process_Func, SyncmerType syncmer_functor,
                             int lcmer_l, bool prime_reverse) {
        uint64_t lcmer_index = 0;
        typename Ops::record *data_p = NULL;
        core_buf str;
        charsi chs;
        pos_scratch pos_buf;
        lcp_read lr;
        gtl::vector<uint8_t> rc_nibbles;
        while (true) {
            if (!data_queue.wait_dequeue_timed(data_p, WAIT_QUEUE_TIMEOUT)) {
                if (done.load()) {
                    break;
                }
                continue;
            }
            all_reads++;
            lr.clear();
            parse_cores_from_nibbles(Ops::seq(data_p), Ops::len(data_p),
                                     lcp_level, &str, syncmer_functor, chs, lr,
                                     pos_buf);
            // Skip reads with too few cores; the record is recycled regardless.
            if (lr.cores.size() >= static_cast<size_t>(lcmer_l + 2)) {
                ++processed_read_count;
                Process_Func(lr, dbgi, thread_id, lcmer_index, lcmer_l);
                if (prime_reverse) {
                    prime_reverse_complement(Ops::seq(data_p), Ops::len(data_p),
                                             lcp_level, &str, syncmer_functor,
                                             chs, pos_buf, rc_nibbles, lr, dbgi,
                                             lcmer_l);
                }
            }
            if (!recycle_queue.try_enqueue(data_p)) {
                Ops::destroy(data_p);
                data_p = NULL;
            }
        }
    }

    // Producer loop shared by every reader: takes a recycled record (or
    // allocates one), lets `fill` populate it from the source, and enqueues it
    // for the workers. Returns when `fill` reports end of input; whether that
    // end was clean is the caller's to determine from its source.
    template <class Ops, class Fill>
    static void enqueue_records(data_queue_of<typename Ops::record> &data_queue,
                                queue_of<typename Ops::record> &recycle_queue,
                                Fill &&fill) {
        typename Ops::record *rec = nullptr;
        while (true) {
            if (rec == nullptr && !recycle_queue.try_dequeue(rec)) {
                rec = Ops::create();
            }
            if (!fill(rec)) {
                Ops::destroy(rec);
                return;
            }
            data_queue.enqueue(rec);
            rec = nullptr;
        }
    }

    // htslib path: BAM, CRAM, SAM.
    [[nodiscard]] static bool read_hts(const stage_input &input,
                                       data_queue_of<bam1_t> &data_queue,
                                       hts_queue &recycle_queue, int threads) {
        samFile *fp = sam_open(input.first.c_str(), "r");
        if (!fp) {
            std::cerr << "Error: Could not open " << input.first << std::endl;
            return false;
        }
        hts_set_threads(fp, threads);
        bam_hdr_t *hdr = sam_hdr_read(fp);
        if (!hdr) {
            std::cerr << "Error: Could not read header of " << input.first
                      << std::endl;
            sam_close(fp);
            return false;
        }
        // sam_read1: >= 0 record, -1 clean EOF, < -1 truncated/corrupt input.
        // A truncated file must fail the run, not silently assemble a subset.
        int rc = 0;
        enqueue_records<hts_record_ops>(data_queue, recycle_queue,
                                        [&](bam1_t *rec) {
                                            rc = sam_read1(fp, hdr, rec);
                                            return rc >= 0;
                                        });
        bam_hdr_destroy(hdr);
        sam_close(fp);
        if (rc < -1) {
            std::cerr << "Error: truncated or corrupt input " << input.first
                      << " (sam_read1 rc=" << rc << ")" << std::endl;
            return false;
        }
        return true;
    }

    // Native FASTA/FASTQ path: no bam1_t rebuild, no quality array, bases are
    // packed into the record's nibble buffer during the block scan.
    [[nodiscard]] static bool read_fastx(const stage_input &input,
                                         data_queue_of<seq_record> &data_queue,
                                         fastx_queue &recycle_queue,
                                         int threads) {
        fastx_reader reader;
        if (!reader.open(input.first, input.second, threads)) {
            std::cerr << "Error: Could not open " << input.first << std::endl;
            return false;
        }
        enqueue_records<fastx_record_ops>(
            data_queue, recycle_queue,
            [&](seq_record *rec) { return reader.next(*rec); });
        if (reader.failed()) {
            std::cerr << "Error: malformed or unreadable FASTA/FASTQ input "
                      << input.first << std::endl;
            return false;
        }
        return true;
    }

    template <class Ops>
    static void drain_and_destroy(auto &q) {
        typename Ops::record *rec = nullptr;
        while (q.try_dequeue(rec)) {
            Ops::destroy(rec);
        }
    }

    // One worker pool over one group of same-format inputs. Returns false if
    // any reader failed. Read counts are accumulated into the caller's sums.
    template <class Ops, class Reader>
    bool run_stage(const gtl::vector<stage_input> &stage_inputs, Reader reader,
                   dbg_index &dbgi, auto Process_Func, u64 &valid_sum,
                   u64 &all_sum) {
        // Records in flight per worker before the readers block. Each holds a
        // whole read, so the bound keeps the buffered reads small.
        constexpr size_t kRecordsPerWorker = 64;
        const size_t capacity =
            kRecordsPerWorker * static_cast<size_t>(threads.work);
        // One reader thread per input in parallel reading mode, else one.
        const size_t producers = [&]() -> size_t {
            if (parallel_file_reading) return stage_inputs.size();
            return 1;
        }();
        data_queue_of<typename Ops::record> data_queue{capacity, producers};
        queue_of<typename Ops::record> recycle_queue{capacity, 1, 1};
        std::atomic<bool> done{false};

        gtl::vector<std::thread> workers;
        gtl::vector<int> valid_read_counts(threads.work, 0);
        gtl::vector<int> all_read_counts(threads.work, 0);

        for (int tid = 1; tid <= threads.work; ++tid) {
            workers.emplace_back([&, tid]() {
                process_data<Ops>(tid, lcp_level, data_queue, recycle_queue,
                                  dbgi, done, valid_read_counts[tid - 1],
                                  all_read_counts[tid - 1], Process_Func,
                                  syncmer_functor, lcmer_l, prime_reverse);
            });
        }

        std::atomic<bool> read_failed{false};
        if (parallel_file_reading) {
            gtl::vector<std::thread> readers;
            readers.reserve(stage_inputs.size());
            for (const auto &input : stage_inputs) {
                readers.emplace_back([&, input]() {
                    if (!reader(input, data_queue, recycle_queue, threads.io)) {
                        fprintf(stderr, "Failed to read from input %s\n",
                                input.first.c_str());
                        read_failed.store(true);
                    }
                });
            }
            for (auto &r : readers) {
                r.join();
            }
        } else {
            for (const auto &input : stage_inputs) {
                if (!reader(input, data_queue, recycle_queue, threads.io)) {
                    fprintf(stderr, "Failed to read from input %s\n",
                            input.first.c_str());
                    read_failed.store(true);
                    break;
                }
            }
        }

        done.store(true);
        for (auto &worker : workers) {
            worker.join();
        }
        drain_and_destroy<Ops>(data_queue);
        drain_and_destroy<Ops>(recycle_queue);

        for (size_t i = 0; i < valid_read_counts.size(); ++i) {
            valid_sum += static_cast<u64>(valid_read_counts[i]);
            all_sum += static_cast<u64>(all_read_counts[i]);
        }
        return !read_failed.load();
    }

    // Runs `Process_Func` over every input, routing each file to the htslib
    // reader or the native FASTX reader according to its detected format.
    bool run_over_inputs(dbg_index &dbgi, auto Process_Func, u64 &valid_sum,
                         u64 &all_sum) {
        gtl::vector<stage_input> hts_inputs;
        gtl::vector<stage_input> fastx_inputs;
        for (const auto &input : inputs) {
            switch (detect_input_format(input)) {
                case input_format::htslib:
                    hts_inputs.emplace_back(input, false);
                    break;
                case input_format::fastx_plain:
                    fastx_inputs.emplace_back(input, false);
                    break;
                case input_format::fastx_compressed:
                    fastx_inputs.emplace_back(input, true);
                    break;
                case input_format::invalid:
                    fprintf(stderr, "Error: cannot open input %s\n",
                            input.c_str());
                    return false;
            }
        }
        bool ok = true;
        if (!hts_inputs.empty()) {
            ok = run_stage<hts_record_ops>(hts_inputs, read_hts, dbgi,
                                           Process_Func, valid_sum, all_sum);
        }
        if (ok && !fastx_inputs.empty()) {
            ok = run_stage<fastx_record_ops>(fastx_inputs, read_fastx, dbgi,
                                             Process_Func, valid_sum, all_sum);
        }
        return ok;
    }

    // Populates the caller-owned dbgi (storage lives in run_assembly).
    std::expected<void, std::string> construct_dbg(dbg_index &dbgi) {
        u64 valid_sum = 0;
        u64 all_sum = 0;
        if (!run_over_inputs(dbgi, process_cores<DefaultBF>, valid_sum,
                             all_sum)) {
            return std::unexpected(
                "Error in reading data for DBG construction.");
        }
        fprintf(stderr, "%" PRIu64 " processed\n%" PRIu64 " total reads\n",
                valid_sum, all_sum);
        const auto &bf = dbgi.bf.implementation;
        const double fill = bf.fill_fraction();
        // Exact for the blocked layout, unlike the unblocked fill^(k - 1).
        const double est_fpr = bf.fpr_for_fill(fill);
        const double keys = bf.distinct_keys_for_fill(fill);
        fprintf(stderr,
                "Bloom filter: %.4f of bits set, estimated false-positive "
                "rate %.2e (about %.4g distinct keys)\n",
                fill, est_fpr, keys);
        // The filter meets --false-positive-probability while it holds at
        // most --projected-element-count distinct keys. Far above that a
        // false positive admits a first-sighting lcmer as a node, so the
        // "seen twice" gate weakens. k <= 1 sets no bits: the filter is
        // deliberately off (eps near 1).
        // Warn an order of magnitude above the configured rate, not at
        // run-to-run fill noise.
        constexpr double kWarnFactor = 10;
        if (bf.k > 1 && est_fpr > kWarnFactor * bloom.false_positive_rate) {
            fprintf(stderr,
                    "Warning: Bloom filter estimated false-positive rate "
                    "%.2e exceeds 10x the configured %.2e; it holds about "
                    "%.3g distinct keys, so raise --projected-element-count "
                    "(now %" PRIu64 ") to at least that\n",
                    est_fpr, bloom.false_positive_rate, keys,
                    bloom.projected_count);
        }
        return {};
    }

    size_t count_seqv_len(const gtl::vector<SequenceContainerType> &seqv) {
        size_t sum = 0;
        for (const SequenceContainerType &s : seqv) {
            sum += s.len;
        }
        return sum;
    }
    // Overlap coordinate space: inclusive core indices on a node's
    // (reverse-complemented) core sequence. `len` is the core-span length in
    // steps (cores - 1).
    struct CoreSpan {
        u64 from;
        u64 to;
        u64 len;
    };

    // Phase-A record: on the key node, an overlap with `other` covers `self`
    // (coordinates on the key node). The partner's span is joined in phase B.
    struct SelfMatch {
        u64 other;
        CoreSpan self;
    };

    // Phase-B record (reciprocal overlap), keyed by the query node: `target`
    // is the partner id; `query`/`target_span` are each node's core span in
    // the shared overlap.
    struct Overlap {
        u64 target;
        CoreSpan query;
        CoreSpan target_span;
    };

    static auto build_lr_cores(const gtl::vector<SequenceContainerType> &seqs,
                               const SyncmerType &syncmer_functor,
                               int lcp_level, core_buf &str) -> lcp_read {
        lcp_read lr;
        lr.clear();

        init_lps_crumbled(&str, seqs, syncmer_functor);
        str.deepen(syncmer_functor, lcp_level);

        str.append_to(lr);
        return lr;
    }

    static auto build_lr_cores_with_seqs(
        const gtl::vector<SequenceContainerType> &seqs,
        const SyncmerType &syncmer_functor, int lcp_level, core_buf &str)
        -> std::pair<lcp_read, gtl::vector<SequenceContainerType>> {
        lcp_read lr;
        lr.clear();
        gtl::vector<SequenceContainerType> core_seqs;

        init_lps_crumbled(&str, seqs, syncmer_functor);
        str.deepen(syncmer_functor, lcp_level);

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

        for (int64_t i = 0; i < str.size; ++i) {
            size_t cur_start = str.start(i);
            size_t cur_end = str.end(i);
            lr.cores.emplace_back(str.label[i]);
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
        const SyncmerType &syncmer_functor, int lcp_level, core_buf &str)
        -> lcp_read {
        lcp_read lr;
        lr.clear();

        gtl::vector<SequenceContainerType> rev_seqs;
        for (auto it = seqs.rbegin(); it != seqs.rend(); ++it) {
            rev_seqs.emplace_back(it->reverse_complement());
        }
        init_lps_crumbled(&str, rev_seqs, syncmer_functor);
        str.deepen(syncmer_functor, lcp_level);

        str.append_to(lr);
        return lr;
    }

    static auto build_rev_lr_cores_with_seqs(
        const gtl::vector<SequenceContainerType> &seqs,
        const SyncmerType &syncmer_functor, int lcp_level, core_buf &str)
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
        str.deepen(syncmer_functor, lcp_level);

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

        for (int64_t i = 0; i < str.size; ++i) {
            size_t cur_start = str.start(i);
            size_t cur_end = str.end(i);
            lr.cores.emplace_back(str.label[i]);
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

        std::atomic<int> all_n{0};
        std::atomic<int> hits{0};
        std::atomic<int> bad{0};
        // Phase-0 matching diagnostics.
        std::atomic<int> too_short{0};   // < 4*lcmer_l cores → no ranges (#4)
        std::atomic<int> no_ranges{0};   // long enough but 0 lcmer hits
        std::atomic<int> short_drop{0};  // reciprocal but below span_floor
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

                core_buf str;
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
                str.deepen(syncmer_functor, lcp_level);

                str.append_to(lr);
                if (lr.cores.size() < static_cast<size_t>(lcmer_l + 1)) {
                    bad.fetch_add(1);
                    // good_n.fetch_sub(1);
                    return;
                }
                // Layout works in core-index space; find_lcmer_ranges returns
                // inclusive [from_core,to_core]. Diagnostic below mirrors its
                // gate via OLC_SKIP_*.
                if (lr.cores.size() <
                    static_cast<size_t>(lcmer_l + OLC_SKIP_FIRST(lcmer_l) +
                                        OLC_SKIP_LAST(lcmer_l))) {
                    too_short.fetch_add(1);  // fails find_lcmer_ranges gate (#4)
                }
                auto found_ranges = find_lcmer_ranges(lr, dbgi, lcmer_l);
                bool any_hit = false;
                for (const auto &[id, from_core, to_core, from_lcmer,
                                  to_lcmer] : found_ranges) {
                    if (id != 0) {
                        any_hit = true;
                        hits.fetch_add(1);

                        // Store the overlap on THIS node in CORE indices
                        // (inclusive; see CoreSpan). Core-space length is
                        // (cores.size() - 1) steps.
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
                        static_cast<size_t>(lcmer_l + OLC_SKIP_FIRST(lcmer_l) +
                                            OLC_SKIP_LAST(lcmer_l)))
                    no_ranges.fetch_add(1);
            },
            threads.work);

        gtl::btree_multimap<u64, Overlap> reciprocal;
        int non_recip = 0, recip_kept = 0;
        // Absolute overlap-span floor in cores: half an lcmer window; below
        // this an overlap is treated as noise.
        const u64 span_floor = std::max<u64>(2, (u64)lcmer_l / 2);
        for (const auto &[k, v] : matched) {
            for (const auto &match : v) {
                // Keep the overlap only if the partner records the reverse one
                // too. `match.self` is k's core span; `ff->self` is the
                // partner's.
                auto it = matched.find(match.other);
                auto ff = (it != matched.end())
                              ? std::find_if(it->second.begin(),
                                             it->second.end(),
                                             [k](const auto &val) {
                                                 return val.other == k;
                                             })
                              : gtl::vector<SelfMatch>::const_iterator{};
                if (it != matched.end() && ff != it->second.end()) {
                    // Overlap extent must be >= span_floor cores on both
                    // nodes. No length-relative gate: a fraction of the
                    // shorter node rejected the few-kb end overlaps that join
                    // Mb-scale nodes.
                    const u64 span_a = match.self.to - match.self.from;
                    const u64 span_b = ff->self.to - ff->self.from;
                    if (std::min(span_a, span_b) < span_floor) {
                        short_drop.fetch_add(1);
                        continue;
                    }
                    reciprocal.emplace(k,
                                       Overlap{match.other, match.self, ff->self});
                    ++recip_kept;
                } else {
                    ++non_recip;
                }
            }
        }

        fprintf(stderr, "all:%d\ngood: %d nodes\nhits: %d\n", all_n.load(),
                all_n.load() - bad.load(), hits.load());
        fprintf(stderr,
                "[match-diag] too_short:%d no_ranges:%d short_drop:%d "
                "recip_kept:%d non_recip:%d\n",
                too_short.load(), no_ranges.load(), short_drop.load(),
                recip_kept, non_recip);
        auto layouts =
            layout_graph_algo::connected_comp_plus_layout(reciprocal);
        gtl::flat_hash_set<u64> dropped;
        if (drop_crumbs) dropped = find_droppable_crumbs(dbgi);
        return std::make_tuple(std::move(reciprocal), std::move(layouts),
                               std::move(dropped));
    }

    // Crumbs: nodes with fewer than lcmer_l+1 cores, too short to form
    // an lcmer and therefore invisible to overlap detection. A crumb with a
    // single live predecessor that continues into a longer node is not
    // emitted as its own C segment. Crumbs that sit between two long nodes
    // (a hub) or whose only continuation is a collapsed repeat (depth above
    // 3x the live median) are kept.
    gtl::flat_hash_set<u64> find_droppable_crumbs(const dbg_index &dbgi) const {
        const u64 min_cores = static_cast<u64>(lcmer_l + 1);
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
            depth_median = 0.5 * (depths[(depths.size() - 1) / 2] +
                                  depths[depths.size() / 2]);
        }
        constexpr double kPileupMaxMult = 3.0;
        const double pileup_cap = depth_median * kPileupMaxMult;
        const auto fwds =
            GraphCleaning::build_live_fwd(dbgi.dbg, dbgi.active_set);
        auto live = [&](u64 nid) {
            auto it = dbgi.dbg.find(nid);
            return (it == dbgi.dbg.end() || it->second.redirected())
                       ? dbgi.dbg.end()
                       : it;
        };

        gtl::flat_hash_set<u64> dropped;
        for (u64 id : dbgi.active_set) {
            auto it = live(id);
            if (it == dbgi.dbg.end()) continue;
            const auto &node = it->second;
            const u64 span = static_cast<u64>(node.seq.size());
            if (span == 0 || span >= min_cores) continue;
            if (node.prev.size() != 1 || live(node.prev[0]) == dbgi.dbg.end()) {
                continue;
            }
            auto fit = fwds.find(id);
            if (fit == fwds.end()) continue;
            u64 long_succs = 0;
            bool continues = false;
            for (u64 q : fit->second) {
                auto qit = live(q);
                if (qit == dbgi.dbg.end() || qit->second.seq.empty()) continue;
                if (static_cast<u64>(qit->second.seq.size()) >= min_cores) {
                    ++long_succs;
                }
                if (pileup_cap <= 0.0 ||
                    qit->second.avg_core_count() <= pileup_cap) {
                    continues = true;
                }
            }
            if (long_succs >= 2 || !continues) continue;
            dropped.insert(id);
        }
        fprintf(stderr, "[crumbs] dropped %zu crumb nodes\n", dropped.size());
        return dropped;
    }
};

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

// Checks shared by both drivers, before any work; prints the configuration.
static bool prepare_run(const asm_config &cfg) {
    if (!validate_path(cfg.output)) return false;

    // Validate the syncmer configuration once, up front, before any syncmer
    // functor is constructed (the functor ctor assumes a valid rule).
    if (auto v = cfg.syncmer.validate(); !v) {
        fprintf(stderr, "Error: %s\n", v.error().c_str());
        return false;
    }

    if (!cfg.silent) cfg.print_cfg(stderr);
    return true;
}

// Debug dump of the DBG's unitigs (--print-raw-dbg, --print-bi-dbg).
static bool dump_dbg(const std::string &path, const dbg_index &dbg) {
    FILE *out = fopen(path.c_str(), "w");
    if (!out) {
        fprintf(stderr, "Cannot open output: %s\n", path.c_str());
        return false;
    }
    print_dbg_seqs(out, dbg);
    fclose(out);
    return true;
}

// Layout of the cleaned DBG, written as the output GFA.
static bool write_output_gfa(const asm_config &cfg, assembler &as,
                             dbg_index &dbg) {
    FILE *out = fopen(cfg.output.c_str(), "w");
    if (!out) {
        fprintf(stderr, "Cannot open output: %s\n", cfg.output.c_str());
        return false;
    }
    layout_graph_emit::emit_ultra_merged_graph(as, dbg, out);
    fclose(out);
    return true;
}

static int run_passes(const asm_config &cfg);

int run_assembly(const asm_config &cfg) {
    if (!cfg.multipass.empty()) return run_passes(cfg);
    if (!prepare_run(cfg)) return 1;
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

    dbg.count_forward_edges(cfg.threads.work);
    IF_TIME_CHECKPOINT(!cfg.silent, HIFI_ASSEMBLY, stderr,
                       "Count forward edges");

    if (cfg.print_precleanup_graph) {
        if (!dump_dbg(cfg.output + ".raw.gfa", dbg)) return 1;
        IF_TIME_CHECKPOINT(!cfg.silent, HIFI_ASSEMBLY, stderr, "DBG printing");
    }

    GraphCleaning::run(dbg, cfg);

    IF_TIME_CHECKPOINT(!cfg.silent, HIFI_ASSEMBLY, stderr, "Graph cleaning");

    // Dump pre-layout (bidirected) unitig graph, gated on --print-bi-dbg.
    if (cfg.print_prelayout_graph && !dump_dbg(cfg.output + ".bi.gfa", dbg)) {
        return 1;
    }

    if (!write_output_gfa(cfg, as, dbg)) return 1;

    IF_TIME_CHECKPOINT(!cfg.silent, HIFI_ASSEMBLY, stderr, "Merged/Printed");

    return 0;
}

static size_t total_bases(const gtl::vector<SequenceContainerType> &seqs) {
    size_t sum = 0;
    for (const auto &s : seqs) sum += s.len;
    return sum;
}

// Sequences of the segments a layout of the cleaned DBG would write: the
// input of the next pass.
static gtl::vector<SequenceContainerType> collect_segments(assembler &as,
                                                           dbg_index &dbg) {
    gtl::vector<SequenceContainerType> segments;
    layout_graph_emit::build_merged_graph(
        as, dbg, "Segment collection",
        [&](const layout_graph_emit::supernode_vec &supernodes,
            const layout_graph_emit::edge_map &) {
            segments = layout_graph_emit::segment_sequences(supernodes,
                                                            as.threads.work);
        });
    return segments;
}

// --multipass driver. Pass p uses the stages of a single pass at its own LCP
// level and lcmer length: reads into a fresh graph (construct_dbg), then the
// previous pass's segments as seeded windows (insert_segments), edge counts,
// cleaning and layout. Earlier passes keep only their segments and free
// their graph; the last pass writes the outputs.
static int run_passes(const asm_config &cfg) {
    if (!prepare_run(cfg)) return 1;

    INIT_TIMERS(HIFI_ASSEMBLY);
    INIT_TIMERS(MULTIPASS);

    assembler as(cfg);
    const std::vector<pass_params> &passes = cfg.multipass;
    const size_t n = passes.size();
    // The previous pass's segments: the only data kept between passes.
    gtl::vector<SequenceContainerType> segments;
    for (size_t p = 0; p < n; ++p) {
        const bool last = p + 1 == n;
        as.lcp_level = passes[p].lcp_level;
        as.lcmer_l = passes[p].lcmer_l;
        if (!cfg.silent)
            fprintf(stderr,
                    "[multipass] pass %zu/%zu: level %d, lcmer-l %d, seeding "
                    "%zu segments (%zu bp, both strands)\n",
                    p + 1, n, as.lcp_level, as.lcmer_l, segments.size(),
                    total_bases(segments));
        // A fresh index (map and Bloom filter) per pass.
        auto graph = std::make_unique<dbg_index>(
            cfg.bloom.false_positive_rate, cfg.bloom.projected_count, cfg.seed);
        dbg_index &dbg = *graph;
        if (auto r = as.construct_dbg(dbg); !r) {
            fprintf(stderr, "%s Terminating.\n", r.error().c_str());
            return 1;
        }
        IF_TIME_CHECKPOINT(!cfg.silent, HIFI_ASSEMBLY, stderr,
                           "DBG construction");

        if (!segments.empty()) {
            const auto [chunks, windows] = as.insert_segments(dbg, segments);
            gtl::vector<SequenceContainerType>().swap(segments);
            IF_TIME_CHECKPOINT(!cfg.silent, HIFI_ASSEMBLY, stderr,
                               "Segment seeding (%zu chunks, %zu windows)",
                               chunks, windows);
        }

        dbg.count_forward_edges(cfg.threads.work);
        IF_TIME_CHECKPOINT(!cfg.silent, HIFI_ASSEMBLY, stderr,
                           "Count forward edges");

        if (last && cfg.print_precleanup_graph) {
            if (!dump_dbg(cfg.output + ".raw.gfa", dbg)) return 1;
            IF_TIME_CHECKPOINT(!cfg.silent, HIFI_ASSEMBLY, stderr,
                               "DBG printing");
        }

        GraphCleaning::run(dbg, cfg);
        IF_TIME_CHECKPOINT(!cfg.silent, HIFI_ASSEMBLY, stderr,
                           "Graph cleaning");

        if (!last) {
            segments = collect_segments(as, dbg);
            dbg.clear_parallel(cfg.threads.work);
            graph.reset();
#if defined(__GLIBC__)
            // Return the freed pass graph to the OS; glibc otherwise keeps
            // most of it in its arenas while the next pass runs.
            malloc_trim(0);
#endif
            IF_TIME_CHECKPOINT(!cfg.silent, MULTIPASS, stderr,
                               "pass %zu/%zu (%d:%d) done", p + 1, n,
                               as.lcp_level, as.lcmer_l);
            if (!cfg.silent)
                fprintf(stderr,
                        "[multipass] pass %zu/%zu: %zu segments (%zu bp) kept "
                        "for the next pass\n",
                        p + 1, n, segments.size(), total_bases(segments));
            continue;
        }

        if (cfg.print_prelayout_graph &&
            !dump_dbg(cfg.output + ".bi.gfa", dbg)) {
            return 1;
        }
        if (!write_output_gfa(cfg, as, dbg)) return 1;
        IF_TIME_CHECKPOINT(!cfg.silent, HIFI_ASSEMBLY, stderr,
                           "Merged/Printed");
        IF_TIME_CHECKPOINT(!cfg.silent, MULTIPASS, stderr,
                           "pass %zu/%zu (%d:%d) done", p + 1, n, as.lcp_level,
                           as.lcmer_l);
        if (!cfg.silent)
            fprintf(stderr, "[multipass] pass %zu/%zu: output written\n", p + 1,
                    n);
        // The last graph is never torn down: the process exits once the GFA
        // is written, so teardown is pure overhead (as NoDestructor).
        (void)graph.release();
    }
    return 0;
}

void tests_mock_fill_seq(dbg_index &dbgi) {
    for (auto &node : dbgi.dbg) {
        node.second.seq.emplace_back();
    }
}

#ifdef LASM_UNIT_TESTS
#include "tests/asm2.test.inc"
#endif
