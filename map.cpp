#include <htslib/hts.h>

#include <algorithm>
#include <chrono>
#include <cstdint>
#include <deque>
#include <iostream>
#include <thread>

#include "asm.h"
#include "blockingconcurrentqueue.h"
#include "config.h"
#include "graph.h"
#include "gtl/phmap.hpp"
#include "gtl/vector.hpp"
#include "htslib/sam.h"
#include "lcmer.h"
#include "lcptools/lcptools_ho.h"
#include "syncmer.h"

struct SWResult {
    int score;
    int start_i, start_j;
    int end_i, end_j;
    gtl::vector<uint32_t> alignedA;
    gtl::vector<uint32_t> alignedB;
};

static inline int sw_match(uint32_t a, uint32_t b) { return (a == b ? 2 : -1); }

//
// Forward SW DP, returns max score and coordinates (end of local alignment)
// Uses O(m) memory
//
void sw_forward(const gtl::vector<uint32_t>& A, const gtl::vector<uint32_t>& B,
                int gap, gtl::vector<int>& prev, gtl::vector<int>& curr,
                int& best_score, int& best_i, int& best_j) {
    int n = A.size();
    int m = B.size();

    prev.assign(m + 1, 0);
    curr.assign(m + 1, 0);

    best_score = 0;
    best_i = best_j = 0;

    for (int i = 1; i <= n; i++) {
        curr[0] = 0;
        for (int j = 1; j <= m; j++) {
            int diag = prev[j - 1] + sw_match(A[i - 1], B[j - 1]);
            int up = prev[j] + gap;
            int left = curr[j - 1] + gap;
            int s = std::max({0, diag, up, left});
            curr[j] = s;

            if (s > best_score) {
                best_score = s;
                best_i = i;
                best_j = j;
            }
        }
        std::swap(prev, curr);
    }
}

//
// Reverse DP: run SW on reversed sequences to find alignment start
// Also O(m) memory
//
void sw_reverse(const gtl::vector<uint32_t>& A, const gtl::vector<uint32_t>& B,
                int gap, int end_i, int end_j, [[maybe_unused]] int best_score,
                gtl::vector<int>& prev, gtl::vector<int>& curr, int& start_i,
                int& start_j) {
    gtl::vector<uint32_t> Ar(A.begin(), A.begin() + end_i);
    gtl::vector<uint32_t> Br(B.begin(), B.begin() + end_j);

    std::reverse(Ar.begin(), Ar.end());
    std::reverse(Br.begin(), Br.end());

    int bi, bj;
    int score;
    sw_forward(Ar, Br, gap, prev, curr, score, bi, bj);

    start_i = end_i - bi;
    start_j = end_j - bj;
}

SWResult smith_waterman_linear(const gtl::vector<uint32_t>& A,
                               const gtl::vector<uint32_t>& B, int gap = -1) {
    gtl::vector<int> prev, curr;

    // 1. Forward DP gives alignment endpoint
    int best_score, end_i, end_j;
    sw_forward(A, B, gap, prev, curr, best_score, end_i, end_j);

    // 2. Reverse DP gives alignment start
    int start_i, start_j;
    sw_reverse(A, B, gap, end_i, end_j, best_score, prev, curr, start_i,
               start_j);

    // 3. Extract subsequences to align exactly
    gtl::vector<uint32_t> subA(A.begin() + start_i, A.begin() + end_i);
    gtl::vector<uint32_t> subB(B.begin() + start_j, B.begin() + end_j);

    // 4. Do a full DP on sub-block (small)
    int na = subA.size();
    int nb = subB.size();
    gtl::vector<gtl::vector<int>> H(na + 1, gtl::vector<int>(nb + 1));

    int max_i = 0, max_j = 0;
    for (int i = 1; i <= na; i++) {
        for (int j = 1; j <= nb; j++) {
            int diag = H[i - 1][j - 1] + sw_match(subA[i - 1], subB[j - 1]);
            int up = H[i - 1][j] + gap;
            int left = H[i][j - 1] + gap;
            H[i][j] = std::max({0, diag, up, left});
            if (H[i][j] > H[max_i][max_j]) {
                max_i = i;
                max_j = j;
            }
        }
    }

    // 5. Traceback on the tiny block
    int i = max_i, j = max_j;
    gtl::vector<uint32_t> aliA, aliB;

    while (i > 0 && j > 0 && H[i][j] > 0) {
        if (H[i][j] == H[i - 1][j - 1] + sw_match(subA[i - 1], subB[j - 1])) {
            aliA.push_back(subA[i - 1]);
            aliB.push_back(subB[j - 1]);
            i--;
            j--;
        } else if (H[i][j] == H[i - 1][j] + gap) {
            aliA.push_back(subA[i - 1]);
            aliB.push_back(UINT32_MAX);
            i--;
        } else {
            aliA.push_back(UINT32_MAX);
            aliB.push_back(subB[j - 1]);
            j--;
        }
    }

    std::reverse(aliA.begin(), aliA.end());
    std::reverse(aliB.begin(), aliB.end());

    return {best_score, start_i, start_j, end_i, end_j, aliA, aliB};
}
struct Minimizer {
    uint64_t hash;
    uint32_t pos;
};

static const uint64_t P = 0x9e3779b97f4a7c15ULL;

inline uint64_t mix_uint(uint32_t x) {
    uint64_t h = x;
    h ^= h >> 33;
    h *= 0xff51afd7ed558ccdULL;
    h ^= h >> 33;
    h *= 0xc4ceb9fe1a85ec53ULL;
    h ^= h >> 33;
    return h;
}

uint64_t hash_kmer(const gtl::vector<uint32_t>& seq, int i, int K) {
    uint64_t h = 0;
    for (int j = 0; j < K; j++) h = h * P ^ mix_uint(seq[i + j]);
    return h;
}

gtl::vector<Minimizer> compute_minimizers(const gtl::vector<uint32_t>& seq,
                                          int K, int W) {
    int N = seq.size();
    if (N < K) return {};
    int num_kmers = N - K + 1;
    if (num_kmers < W) return {};

    gtl::vector<uint64_t> hashes(num_kmers);

    for (int i = 0; i < num_kmers; i++) hashes[i] = hash_kmer(seq, i, K);

    std::deque<int> dq;
    gtl::vector<Minimizer> mins;
    mins.reserve(num_kmers / W + 10);

    for (int i = 0; i < num_kmers; i++) {
        while (!dq.empty() && hashes[dq.back()] >= hashes[i]) dq.pop_back();

        dq.push_back(i);

        int left = i - W + 1;
        if (dq.front() < left) dq.pop_front();

        if (i >= W - 1) {
            int pos = dq.front();
            if (!mins.empty()) {
                if (mins.back().pos == (uint32_t)pos) {
                    continue;
                }
            }
            mins.push_back({hashes[pos], (uint32_t)pos});
        }
    }

    return mins;
}
struct SeedLoc {
    uint32_t contig;
    uint32_t pos;
};

using MinimizerIndex = gtl::parallel_flat_hash_map<
    uint64_t, gtl::vector<SeedLoc>, std::hash<uint64_t>,
    std::equal_to<uint64_t>,
    std::allocator<std::pair<const uint64_t, gtl::vector<SeedLoc>>>, 8,
    std::mutex>;

MinimizerIndex build_minimizer_index(
    const gtl::vector<std::pair<std::string, lcp_read>>& lcp_ctgs, int K, int W,
    int tcount) {
    MinimizerIndex index;

    // Pre-allocate based on heuristic:
    // assume about (#bases / W) minimizers total
    size_t total_bases = 0;
    for (const auto& c : lcp_ctgs) total_bases += c.second.cores.size();
    index.reserve(total_bases / W + 16);

    auto work = [tcount, &lcp_ctgs, &index, K, W](int tid) {
        for (uint32_t cid = tid; cid < lcp_ctgs.size(); cid += tcount) {
            const auto& ref = lcp_ctgs[cid].second.cores;

            gtl::vector<Minimizer> mins = compute_minimizers(ref, K, W);

            for (const auto& m : mins) {
                index.lazy_emplace_l(
                    m.hash,
                    [cid, &m](auto& val) {
                        val.second.emplace_back(cid, m.pos);
                    },
                    [cid, &m](MinimizerIndex::constructor ctr) {
                        ctr(m.hash, gtl::vector<SeedLoc>{{cid, m.pos}});
                    });
                // auto it = index.find(m.hash);

                // if (it == index.end()) {
                //     gtl::vector<SeedLoc> v;
                //     v.push_back({cid, m.pos});
                //     index.emplace(m.hash, std::move(v));
                // } else {
                //     it->second.push_back({cid, m.pos});
                // }
            }
        }
    };

    gtl::vector<std::thread> threads;
    for (int i = 0; i < tcount; ++i) {
        threads.emplace_back(work, i);
    }
    for (auto& t : threads) {
        t.join();
    }
    return index;
}
struct mapper {
    using SYNCMER_HASHER = syncmer::MurMur;
    using SyncmerType = syncmer::GenericSyncmerRule<uint32_t, SYNCMER_HASHER>;
    const asm_config cfg;
    SyncmerType syncmer_functor;

    mapper(const asm_config& cfg)
        : cfg{cfg}, syncmer_functor(cfg.syncmer, cfg.seed) {}
    auto read_hts(const std::string& path, bool reversed) {
        std::mutex gm;
        gtl::vector<std::pair<std::string, lcp_read>> seqs;

        auto atomic_emplace = [&](auto& vec, const std::string& name,
                                  lcp_read&& element) {
            std::lock_guard<std::mutex> lock(gm);
            return vec.emplace_back(name, element);
        };
        moodycamel::BlockingConcurrentQueue<bam1_t*> queue{1024, 1, 1};
        std::atomic<bool> done(false);
        auto work = [&queue, &done, &seqs, this, reversed, atomic_emplace]() {
            charsi chs;
            pos_scratch pos_buf;
            bam1_t* rec = NULL;
            lps str = {0, 0, 0, 0};

            while (true) {
                if (queue.wait_dequeue_timed(rec,
                                             std::chrono::milliseconds(1))) {
                    lcp_read read;
                    parse_cores_from_bam_record(
                        bam_get_seq(rec), rec->core.l_qseq, cfg.lcp_level, &str,
                        syncmer_functor, chs, read, pos_buf, reversed);
                    atomic_emplace(seqs, bam_get_qname(rec), std::move(read));
                    bam_destroy1(rec);
                    lps_clear(&str);
                    rec = NULL;
                } else {
                    if (done.load()) {
                        return;
                    }
                    std::this_thread::sleep_for(std::chrono::milliseconds(1));
                }
            }
        };

        htsFile* fp = hts_open(path.c_str(), "r");
        if (fp == NULL) {
            std::cerr << "Failed to open file: " << path << std::endl;
            exit(-1);
        }
        hts_set_threads(fp, cfg.threads.io);
        bam_hdr_t* hdr = NULL;
        bam1_t* rec = nullptr;
        gtl::vector<std::thread> threads;
        for (int i = 0; i < cfg.threads.work; ++i) {
            threads.emplace_back(work);
        }
        while (true) {
            // Try to get a recycled record object
            if (rec == nullptr) {
                rec = bam_init1();
            }
            if (sam_read1(fp, hdr, rec) >= 0) {
                while (!queue.try_enqueue(rec)) {
                    std::this_thread::sleep_for(std::chrono::milliseconds(1));
                }

                rec = nullptr;
            } else {
                // End of file or error. Destroy the unused record and exit.
                bam_destroy1(rec);
                rec = nullptr;
                break;
            }
        }
        done.store(true);
        for (auto& thread : threads) {
            thread.join();
        }
        bam_hdr_destroy(hdr);
        sam_close(fp);
        return seqs;
    }
};
struct SeedHit {
    uint32_t qpos;
    uint32_t contig;
    uint32_t rpos;
};
gtl::vector<SeedHit> find_query_seeds(const gtl::vector<uint32_t>& query,
                                      const MinimizerIndex& index, int K,
                                      int W) {
    gtl::vector<SeedHit> hits;

    gtl::vector<Minimizer> qmins = compute_minimizers(query, K, W);

    // Reserve approx (#minimizers * average bucket size)
    hits.reserve(qmins.size() * 2);

    for (const auto& qm : qmins) {
        auto it = index.find(qm.hash);
        if (it == index.end()) continue;  // no match
        const gtl::vector<SeedLoc>& bucket = it->second;
        for (const SeedLoc& loc : bucket) {
            hits.push_back({qm.pos, loc.contig, loc.pos});
        }
    }

    return hits;
}

struct Chain {
    gtl::vector<int> indices;  // seed indices in original hits vector
    int score;
};

// Input: hits for one contig (unsorted)
// Output: multiple best chains

struct Interval {
    u32 ctg;
    u32 s;
    u32 e;
};
struct Segment {
    Interval qry;
    Interval ref;
    int elements;
};
gtl::vector<Segment> chain_multiple(u32 me, gtl::vector<SeedHit>& hits) {
    // Simplify hit streaks
    gtl::vector<Segment> segments;
    if (hits.empty()) {
        return {};
    }
    if (hits.size() == 1) {
        return {
            {Segment{Interval{me, hits[0].qpos, hits[0].qpos},
                     Interval{hits[0].contig, hits[0].rpos, hits[0].rpos}, 1}}};
    }

    std::sort(hits.begin(), hits.end(), [](const SeedHit& a, const SeedHit& b) {
        if (a.contig != b.contig) return a.contig < b.contig;
        if (a.qpos != b.qpos) return a.qpos < b.qpos;
        return a.rpos < b.rpos;
    });
    auto left = hits.begin();
    auto it = std::next(left);

    while (it != hits.end() && std::next(it) != hits.end()) {
        if (it->contig == left->contig) {
            ++it;
        } else {
            segments.emplace_back(Interval{me, left->qpos, it->qpos},
                                  Interval{left->contig, left->rpos, it->rpos},
                                  (int)(it - left));
            ++it;
            left = it;
            ++it;
        }
    }
    if (left + 1 != it)
        segments.emplace_back(Interval{me, left->qpos, it->qpos},
                              Interval{left->contig, left->rpos, it->rpos},
                              (int)(it - left));

    return segments;
}

int run_mapping(const asm_config& param) {
    INIT_TIMERS(Mapping)
    mapper map(param);

    if (param.inputs.empty()) {
        fprintf(stderr, "Error: at least one input file is required.\n");
        return -1;
    }
    auto ref = map.read_hts(param.inputs.front(), false);
    auto qry = map.read_hts(param.output, false);
    TIME_CHECKPOINT(Mapping, stderr, "Read input files %s and %s",
                    param.inputs.front().c_str(), param.output.c_str());
    // auto ref = map.read_hts(
    //     param.output);  // I am reusing output here for passing str, I
    //     dont want
    //                     // to rewiret or modify the asm_confgi

    int K = 3;
    int W = 2;
    // gtl::vector<gtl::flat_hash_set<uint32_t>> rs;
    // gtl::vector<gtl::flat_hash_set<uint32_t>> qs;

    auto mindex = build_minimizer_index(ref, K, W, param.threads.work);
    TIME_CHECKPOINT(Mapping, stderr, "Built index");

    std::atomic<bool> done{false};

    moodycamel::BlockingConcurrentQueue<size_t> jobqueue{qry.size() * 2};
    moodycamel::BlockingConcurrentQueue<
        std::tuple<size_t, gtl::vector<SeedHit>, gtl::vector<Segment>>>
        resultqueue{qry.size() * 2};

    auto worker_fnc = [&done, &jobqueue, &qry, &mindex, K, W, &resultqueue]() {
        size_t job = 0;
        while (true) {
            if (jobqueue.wait_dequeue_timed(job,
                                            std::chrono::milliseconds(1))) {
                auto hits =
                    find_query_seeds(qry.at(job).second.cores, mindex, K, W);
                auto chains = chain_multiple(job, hits);
                while (!resultqueue.try_enqueue(
                    std::make_tuple(job, std::move(hits), std::move(chains)))) {
                    std::this_thread::sleep_for(std::chrono::milliseconds(1));
                }
            } else {
                if (done.load() && jobqueue.size_approx() == 0) {
                    break;
                }
            }
        }
    };
    gtl::vector<std::thread> threads;
    for (size_t i = 0; i < qry.size(); ++i) {
        jobqueue.enqueue(i);
    }
    for (int i = 0; i < param.threads.work; ++i) {
        threads.emplace_back(worker_fnc);
    }

    size_t pcount = 0;
    std::tuple<size_t, gtl::vector<SeedHit>, gtl::vector<Segment>> qitem;
    while (pcount < qry.size()) {
        if (resultqueue.wait_dequeue_timed(qitem,
                                           std::chrono::milliseconds(1))) {
            printf("Chains: %zu[\n", std::get<2>(qitem).size());

            int i = 0;
            for (const auto& chain : std::get<2>(qitem)) {
                printf("Chain %d:(%s:%d-%d) (%s:%d-%d) #%d\n", i,
                       qry[chain.qry.ctg].first.c_str(), chain.qry.s,
                       chain.qry.e, ref[chain.ref.ctg].first.c_str(),
                       chain.ref.s, chain.ref.e, chain.elements);
                ++i;
            }

            printf("]\nHits: %zu\n", std::get<1>(qitem).size());
            printf("%zu: [", std::get<0>(qitem));
            for (const auto& hit : std::get<1>(qitem)) {
                printf("(%s, %d, %d) ", ref[hit.contig].first.c_str(), hit.qpos,
                       hit.rpos);
            }
            printf("]\n");
            ++pcount;
        }
    }
    done.store(true);
    for (std::thread& t : threads) {
        t.join();
    }
    TIME_CHECKPOINT(Mapping, stderr, "Finished mapping");
    // for (const lcp_read& lr : ref) {
    //     rs.emplace_back(lr.cores.begin(), lr.cores.end());
    // }
    // for (const lcp_read& lr : qry) {
    //     qs.emplace_back(lr.cores.begin(), lr.cores.end());
    // }

    // for (size_t i = 0; i < rs.size(); ++i) {
    //     for (size_t j = 0; j < qs.size(); ++j) {
    //         gtl::vector<uint32_t> out;
    //         std::ranges::set_intersection(rs[i], qs[j],
    //                                       std::back_inserter(out));
    //         printf("%zu ", out.size());
    //     }
    //     printf("\n");
    // }
    // printf("\n\n-------------\n\n");

    // for (const lcp_read& lr: ref) {
    //     if(lr.cores.size() < 4){
    //         ++i;
    //         continue;
    //     }
    //     gtl::vector<std::pair<size_t, int>> nice_matches;
    //     int j = 0;
    //     for (const lcp_read& lq : qry) {
    //         if(lq.cores.size() < 4){
    //             ++j;
    //             continue;
    //         }
    //         auto results =  smith_waterman_linear(lr.cores, lq.cores);
    //         if(results.score > 0){
    //             nice_matches.emplace_back(j, results.score);
    //         }
    //         ++j;
    //     }
    //     // Print the matches in python list(tuple) format
    //     printf("%d: [", i);
    //     for(auto& match : nice_matches){
    //         printf("(%zu, %d)", match.first, match.second);
    //     }
    //     printf("]\n");
    //     ++i;
    // }
    return 0;
}
