#ifndef _LASM_CONFIG_H_
#define _LASM_CONFIG_H_

#ifndef BAMREAD_CHUNK_SIZE
#define BAMREAD_CHUNK_SIZE 256ULL
#endif

#ifndef WAIT_QUEUE_TIMEOUT_MS
#define WAIT_QUEUE_TIMEOUT_MS 2
#endif

#ifndef WAIT_QUEUE_TIMEOUT
#define WAIT_QUEUE_TIMEOUT std::chrono::milliseconds(WAIT_QUEUE_TIMEOUT_MS)
#endif

#ifndef MAXIMUM_WORK_THREADS
#define MAXIMUM_WORK_THREADS 127
#endif

#define INIT_TIMERS(name)                                                   \
    std::chrono::high_resolution_clock::time_point name##_chrono_start =    \
        std::chrono::high_resolution_clock::now();                          \
    std::chrono::high_resolution_clock::time_point name##_last_checkpoint = \
        name##_chrono_start;                                                \
    clock_t name##_clock_start = std::clock();                              \
    clock_t name##_last_clock = name##_clock_start;

#define IF_TIME_CHECKPOINT(cond, name, stream, fmt, ...)                       \
    do {                                                                       \
        if (cond) {                                                            \
            auto name##_now = std::chrono::high_resolution_clock::now();       \
            clock_t name##_clock_now = std::clock();                           \
                                                                               \
            long long name##_total_elapsed =                                   \
                std::chrono::duration_cast<std::chrono::milliseconds>(         \
                    name##_now - name##_chrono_start)                          \
                    .count();                                                  \
            long long name##_interval_elapsed =                                \
                std::chrono::duration_cast<std::chrono::milliseconds>(         \
                    name##_now - name##_last_checkpoint)                       \
                    .count();                                                  \
                                                                               \
            double name##_total_cpu =                                          \
                1000.0 * (name##_clock_now - name##_clock_start) /             \
                CLOCKS_PER_SEC;                                                \
            double name##_interval_cpu =                                       \
                1000.0 * (name##_clock_now - name##_last_clock) /              \
                CLOCKS_PER_SEC;                                                \
                                                                               \
            double name##_total_util =                                         \
                (name##_total_elapsed > 0)                                     \
                    ? (name##_total_cpu / name##_total_elapsed) * 100.0        \
                    : 0.0;                                                     \
            double name##_interval_util =                                      \
                (name##_interval_elapsed > 0)                                  \
                    ? (name##_interval_cpu / name##_interval_elapsed) * 100.0  \
                    : 0.0;                                                     \
                                                                               \
            fprintf(stream, "[Checkpoint: %s] ", #name);                       \
            fprintf(stream, fmt __VA_OPT__(, ) __VA_ARGS__);                   \
            fprintf(stream, "\n");                                             \
            fprintf(stream,                                                    \
                    "  Total:   Wall = %lld ms, CPU = %.2f ms, Utilization = " \
                    "%.2f%%\n",                                                \
                    name##_total_elapsed, name##_total_cpu,                    \
                    name##_total_util);                                        \
            fprintf(stream,                                                    \
                    "  Current: Wall = %lld ms, CPU = %.2f ms, Utilization = " \
                    "%.2f%%\n\n",                                              \
                    name##_interval_elapsed, name##_interval_cpu,              \
                    name##_interval_util);                                     \
                                                                               \
            name##_last_checkpoint = name##_now;                               \
            name##_last_clock = name##_clock_now;                              \
        }                                                                      \
    } while (0)
#define TIME_CHECKPOINT(...) IF_TIME_CHECKPOINT(true, __VA_ARGS__)

#endif
