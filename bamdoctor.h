#ifndef _BAMDOCTOR_H_
#define _BAMDOCTOR_H_

#include <string>
#include <htslib/sam.h>
#include <htslib/hts.h>
#include <htslib/thread_pool.h>
#include <stdio.h>
#include "config.h"
#include <thread>
#include <atomic>
#include <concurrentqueue.h>
#include <gtl/vector.hpp>
#include "lcmer.h"

#include "syncmer.h"

enum class SamReadTestingMode {
    BASIC,
    POOL,
};

inline int run_sam_read_testing(const std::string &input, int threads, int lcp_level, SamReadTestingMode mode, const std::string &job_type, double density){
    INIT_TIMERS(BAM_READING);
    //Explain the program nicely formatted brief
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
        return 1;
    }
    bam1_t *bamdata = bam_init1();
    if(!bamdata){
        fprintf(stderr, "Error: Could not initialize bamdata\n");
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
                return 1;
            }
            if(hts_set_opt(fp, HTS_OPT_THREAD_POOL, &pool) < 0){
                fprintf(stderr, "Error: Could not set thread pool\n");
                return 1;
            }
        }
    }
    moodycamel::ConcurrentQueue<bam1_t *> data_queue {1024, 1, static_cast<size_t>(threads)};
    moodycamel::ConcurrentQueue<bam1_t *> recycle_queue {1024, 1, static_cast<size_t>(threads)};
    gtl::vector<std::thread> workers;
    std::atomic<bool> done{false};
    int counts[1024] = {0};
    syncmer::DensitySyncmerRule<uint32_t, syncmer::MurMur> sr{8, density, 423, static_cast<unsigned int>(-1)};
    charsi chs;
    auto job = [&]() -> std::function<int(bam1_t *, lps *)> {
        if (job_type == "read_count"){
            return [](bam1_t *, lps *) -> int {
                return 1;
            };
        }
        else if (job_type == "core_count"){
            return [&, lcp_level](bam1_t *bamdata, lps *str) -> int {
                lcp_read lr;
                lr.clear();
                syncmer_lcp::parse_cores_from_bam_record(bam_get_seq(bamdata), bamdata->core.l_qseq, lcp_level, str, sr, chs, lr);
                // if (lr.cores.size() < LCMER_L + 2) {
                //     return 0;
                // }

                return (int)lr.cores.size();
            };
        }
        return [](bam1_t *, lps *) -> int {
            return 0;
        };
    }();

    for(int i = 0; i < threads; ++i){
        workers.emplace_back([&, i](){
            bam1_t *bamdata = NULL;
            lps str = {0, 0, 0, nullptr};
            while(!done.load() || data_queue.size_approx() > 0){
                if(!data_queue.try_dequeue(bamdata)){
                    std::this_thread::sleep_for(std::chrono::nanoseconds(10));
                    continue;
                }
                counts[i] += job(bamdata, &str);
                // ++counts[i];
                if(!recycle_queue.try_enqueue(bamdata)){
                    bam_destroy1(bamdata);
                }
            }
            if(str.cores != nullptr){
                free_lps(&str);
            }
        });
    }

    while(sam_read1(fp, hdr, bamdata) >= 0){
        while(!data_queue.try_enqueue(bamdata)){
            std::this_thread::sleep_for(std::chrono::nanoseconds(100));
        }
        if(!recycle_queue.try_enqueue(bamdata)){
            // bam_destroy1(bamdata);
            bamdata = bam_init1();
        }
    }
    done.store(true);
    for(auto &worker : workers){
        worker.join();
    }
    bam1_t *item;
    size_t count =  0;
    while(data_queue.try_dequeue(item)){
        ++count;
        // bam_destroy1(item);
    }
    for(int i = 0; i < threads; ++i){
        count += counts[i];
    }
    // while(recycle_queue.try_dequeue(item)){
    //     bam_destroy1(item);
    // }
    TIME_CHECKPOINT(BAM_READING, stderr, "Read %zu reads", count);

    // bam_destroy1(bamdata);
    sam_hdr_destroy(hdr);
    sam_close(fp);

    return 0;
}


#endif
