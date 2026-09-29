#ifndef _BAMDOCTOR_H_
#define _BAMDOCTOR_H_

#include <string>
#include <vector>
#include <functional>
#include <htslib/sam.h>
#include <htslib/hts.h>
#include <htslib/thread_pool.h>
#include <stdio.h>
#include "config.h"
#include <thread>
#include <atomic>
#include <concurrentqueue.h>
#include <gtl/vector.hpp>
#include "asm.h"
#include "lcmer.h"

#include "syncmer.h"

enum class SamReadTestingMode {
    BASIC,
    POOL,
};

inline int run_sam_read_testing(const std::string &input, int threads,
                                int lcp_level, SamReadTestingMode mode,
                                const std::string &job_type,
                                const syncmer_params &sr_params) {
    INIT_TIMERS(BAM_READING);
    printf("Input: %s\n", input.c_str());
    printf("Threads: %d\n", threads);
    printf("LCP level: %d\n", lcp_level);
    printf("Mode: %s\n", mode == SamReadTestingMode::BASIC ? "BASIC" : "POOL");


    // Read bam file with htslib
    samFile *fp = sam_open(input.c_str(), "r");
    if(!fp){
        fprintf(stderr, "Error: Could not open %s\n", input.c_str());
        return 1;
    }
    sam_hdr_t *hdr = sam_hdr_read(fp);
    if(!hdr){
        fprintf(stderr, "Error: Could not read header from file %s\n", input.c_str());
        sam_close(fp);
        return 1;
    }
    bam1_t *bamdata = bam_init1();
    if(!bamdata){
        fprintf(stderr, "Error: Could not initialize bamdata\n");
        sam_hdr_destroy(hdr);
        sam_close(fp);
        return 1;
    }
    htsThreadPool pool = {NULL, 0};
    if(mode == SamReadTestingMode::BASIC){
        hts_set_threads(fp, threads);
    }else if(mode == SamReadTestingMode::POOL){

        if(threads > 1){
            pool.pool = hts_tpool_init(threads);
            if(!pool.pool){
                fprintf(stderr, "Error: Could not initialize thread pool\n");
                bam_destroy1(bamdata);
                sam_hdr_destroy(hdr);
                sam_close(fp);
                return 1;
            }
            if(hts_set_opt(fp, HTS_OPT_THREAD_POOL, &pool) < 0){
                fprintf(stderr, "Error: Could not set thread pool\n");
                hts_tpool_destroy(pool.pool);
                bam_destroy1(bamdata);
                sam_hdr_destroy(hdr);
                sam_close(fp);
                return 1;
            }
        }
    }
    moodycamel::ConcurrentQueue<bam1_t *> data_queue {1024, 1, static_cast<size_t>(threads)};
    moodycamel::ConcurrentQueue<bam1_t *> recycle_queue {1024, 1, static_cast<size_t>(threads)};
    gtl::vector<std::thread> workers;
    std::atomic<bool> done{false};
    // One slot per worker, sized by the thread count.
    std::vector<int> counts(threads, 0);
    // Same rule and seed as the assembly pipeline at its default
    // --random-seed (GenericSubmerRule over MurMur).
    syncmer::GenericSubmerRule<uint32_t, syncmer::MurMur> sr{
        sr_params, kDefaultRandomSeed};
    // Per-worker scratch (charsi/pos_scratch) is passed as arguments, not
    // captured, so the shared job closure stays race-free.
    auto job = [&]()
        -> std::function<int(bam1_t *, core_buf *, charsi &, pos_scratch &)> {
        if (job_type == "read_count"){
            return [](bam1_t *, core_buf *, charsi &, pos_scratch &) -> int {
                return 1;
            };
        }
        else if (job_type == "core_count"){
            return [&, lcp_level](bam1_t *bamdata, core_buf *str, charsi &chs,
                                  pos_scratch &pos_buf) -> int {
                lcp_read lr;
                lr.clear();
                parse_cores_from_nibbles(bam_get_seq(bamdata), bamdata->core.l_qseq, lcp_level, str, sr, chs, lr, pos_buf);
                // if (lr.cores.size() < LCMER_L + 2) {
                //     return 0;
                // }

                return (int)lr.cores.size();
            };
        }
        return [](bam1_t *, core_buf *, charsi &, pos_scratch &) -> int {
            return 0;
        };
    }();

    for(int i = 0; i < threads; ++i){
        workers.emplace_back([&, i](){
            bam1_t *bamdata = NULL;
            core_buf str;
            charsi chs;
            pos_scratch pos_buf;
            while(!done.load() || data_queue.size_approx() > 0){
                if(!data_queue.try_dequeue(bamdata)){
                    std::this_thread::sleep_for(std::chrono::nanoseconds(10));
                    continue;
                }
                counts[i] += job(bamdata, &str, chs, pos_buf);
                // Hand the buffer back for the producer to refill.
                if(!recycle_queue.try_enqueue(bamdata)){
                    bam_destroy1(bamdata);
                }
            }
        });
    }

    // Fill bamdata, hand it to a worker, then take the next buffer from the
    // recycle queue (or allocate one).
    while(sam_read1(fp, hdr, bamdata) >= 0){
        while(!data_queue.try_enqueue(bamdata)){
            std::this_thread::sleep_for(std::chrono::nanoseconds(100));
        }
        if(!recycle_queue.try_dequeue(bamdata)){
            bamdata = bam_init1();
        }
    }
    done.store(true);
    for(auto &worker : workers){
        worker.join();
    }
    size_t count = 0;
    for(int i = 0; i < threads; ++i){
        count += counts[i];
    }
    // Reclaim every buffer after the joins: drain data_queue and
    // recycle_queue, then destroy the spare bamdata read buffer.
    bam1_t *item;
    while(data_queue.try_dequeue(item)){
        bam_destroy1(item);
    }
    while(recycle_queue.try_dequeue(item)){
        bam_destroy1(item);
    }
    bam_destroy1(bamdata);
    TIME_CHECKPOINT(BAM_READING, stderr, "Read %zu reads", count);

    sam_hdr_destroy(hdr);
    sam_close(fp);
    if(pool.pool){
        hts_tpool_destroy(pool.pool);
    }

    return 0;
}


#endif
