#include "bgzf.h"
#include "ext/termbox2.h"
#include "gtl/phmap.hpp"
#include "histogram.h"
#include "lcmer.h"
#include "syncmer.h"
template <class T>
struct stat_accu {
    T count;
    T mean;
    T dS;
    T min;
    T max;

    stat_accu() : count(0), mean(0), dS(0), min(0), max(0) {}

    void update(T val) {
        count++;
        auto mean_dif = (val - mean) / count;
        auto new_mean = mean + mean_dif;

        auto dS_dif = (val - new_mean) * (val - mean);
        auto new_dS = dS + dS_dif;

        mean = new_mean;
        dS = new_dS;
    }

    constexpr T variance() const { return dS / count; }

    constexpr T stddev() const { return sqrt(variance()); }
};

int run_stats(const std::string &input, int update_frequency, int lcmer_level, double density) {
    BAMReader reader(input);
    bam1_t *bamdata = reader.init_data();
    // uint8_t *data;
    // uint64_t read_count;
    uint64_t base_count = 0;

    gtl::flat_hash_map<std::string, stat_accu<double>> stats;
    uint64_t read_count = 0;

    Screen screen(2, 2);
    int read_len_hist = screen.addHistogram<int>(0, 0, "Read Lengths");
    int core_count_hist = screen.addHistogram<int>(0, 1, "Core counts");
    int core_len_hist = screen.addHistogram<int>(1, 0, "Core lengths");
    int text = screen.addText<uint64_t, uint64_t, double, double, double,
                              double, double, double>(
        1, 1,
        "read count: {}\nbase-count: {}\nlength(µ,σ): {:.2f} "
        "{:.2f}\ncores(µ,σ): {:.2f} {:.2f}\ncore-lengths(µ,σ): {:.2f} "
        "{:.2f}");

    struct tb_event ev;
    lps str = {0, 0, 0, 0};
    syncmer::DensitySyncmerRule<uint32_t, syncmer::MurMur> syncmer_functor{
        8,
        density, 423U, static_cast<unsigned int>(-1)};
    charsi chs;
    lcp_read lcps;
    while (reader.read(bamdata)) {
        lcps.clear();
         syncmer_lcp::parse_cores_from_bam_record(
            bam_get_seq(bamdata), bamdata->core.l_qseq, lcmer_level, &str,
            syncmer_functor, chs, lcps);
        stats["length"].update(bamdata->core.l_qseq);
        stats["core-count"].update(lcps.cores.size());
        for (const auto &len : lcps.pos) {
            screen.pushHistogram<int>(core_len_hist, len.second - len.first);
            stats["core-length"].update(len.second - len.first);
        }
        base_count += bamdata->core.l_qseq;
        screen.pushHistogram<int>(read_len_hist, bamdata->core.l_qseq);
        screen.pushHistogram<int>(core_count_hist, lcps.cores.size());
        screen.setTextValues<uint64_t, uint64_t, double, double, double, double,
                             double, double>(
            text, read_count, base_count, stats["length"].mean,
            stats["length"].stddev(), stats["core-count"].mean,
            stats["core-count"].stddev(), stats["core-length"].mean,
            stats["core-length"].stddev());
        if (read_count % update_frequency == 0) {
            screen.update();
        }
        ++read_count;
        // screen.update();
        while (tb_peek_event(&ev, 0) != TB_ERR_NO_EVENT) {
            if (ev.type == TB_EVENT_KEY &&
                (ev.key == TB_KEY_CTRL_C || ev.key == TB_KEY_ESC ||
                 ev.key == TB_KEY_ENTER))
                goto enough;
        }
    }
enough:
    while (tb_peek_event(&ev, 100) == TB_ERR_NO_EVENT) {
        if (ev.type == TB_EVENT_KEY) goto enough;
    }
    lps_clear(&str);
    return 0;
}
