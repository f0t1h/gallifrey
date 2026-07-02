#include "lcmer.h"

#include <bit>
#include <cstdint>

#include "asm.h"
#include "ext/doctest.h"
#include "ext/lcptools/lcptools_ho.h"
#include "gtl/phmap.hpp"
#include "htslib/hts.h"
#include "htslib/sam.h"
#include "nibble2base.hpp"
#include "syncmer.h"

void init_core1bam(struct core *cr, const uint8_t *data, uint64_t distance,
                   uint64_t si, uint64_t ei) {
    cr->start = si;
    cr->end = ei;
    cr->label = 0;
    cr->label |= ((distance - 2) << 12);
    cr->label |= (bam_seqi(data, si) << 8);
    cr->label |= (bam_seqi(data, si + distance - 2) << 4);
    cr->label |= (bam_seqi(data, si + distance - 1));
    cr->bit_rep = lcp_bit_type{0x8000000000000000, cr->label};
    cr->bit_size = 2 * distance;
}


uint64_t parse1_bam_lmax(const uint8_t *data, int64_t len, struct core *cores) {
    int64_t it1 = 0;
    int64_t it2 = len;

    int core_index = 0;
    // int64_t last_invalid_char_index = -1;

    // find lcp cores
    for (; it1 + 2 < len; it1++) {
        if (bam_seqi(data, it1) == bam_seqi(data, it1 + 1)) {
            continue;
        }
        if (0 == it1) {
            continue;
        }

        // check for LMAX
        if (it1 + 3 < len && bam_seqi(data, it1) < bam_seqi(data, (it1 + 1)) &&
            bam_seqi(data, (it1 + 1)) > bam_seqi(data, (it1 + 2)) &&
            bam_seqi(data, (it1 - 1)) <= bam_seqi(data, (it1)) &&
            bam_seqi(data, (it1 + 2)) >= bam_seqi(data, (it1 + 3))) {

            // create LMAX core
            it2 = it1 + 3;
            init_core1bam(&(cores[core_index]), data, it2 - it1, it1, it2);
            core_index++;

            continue;
        }
        // cout << "None\n";
    }
    // cout << "CORESSS" << it1 << "\t" << it2 << "\n";
    return core_index;
}
int onehotto2bitDNA2(const uint8_t *data, char *ret) {
    static char table[] = {'A', 'C', 'G', 'T'};
    for (size_t i = 0;; ++i) {
        int clz = std::countl_zero(static_cast<unsigned>(data[i]));
        int ctz = std::countr_zero(static_cast<unsigned>(data[i]));

        if (clz < 0 || ctz > 3) {
            ret[i] = 0;
            return i;
        }
        ret[2 * i] = table[27 - clz];
        ret[2 * i + 1] = table[ctz];
    }
    return -1;
}
uint8_t *onehotto2bitDNA(uint8_t *data, size_t len) {
    uint8_t *ret = (uint8_t *)malloc(len);
    for (size_t i = 0; i < len / 2; ++i) {
        ret[2 * i] = 27 - std::countl_zero(static_cast<unsigned>(data[i]));
        ret[2 * i + 1] = std::countr_zero(static_cast<unsigned>(data[i]));
    }
    return ret;
}
void onehotto2bitDNA_bivec_chs(const uint8_t *data, size_t len, charsi &chs) {
    chs.require(len);
#if defined(__ARM_NEON) || defined(__ARM_NEON__)
    static const uint8_t decode4_tbl_arr[16] = {
        0,0,1,0, 2,0,0,0, 3,0,0,0, 0,0,0,0};
    const uint8x16_t tbl = vld1q_u8(decode4_tbl_arr);
    const uint8x16_t mask = vdupq_n_u8(0x0F);
    size_t i = 0;
    for (; i + 16 <= len / 2; i += 16) {
        uint8x16_t v = vld1q_u8(data + i);
        uint8x16_t hi = vqtbl1q_u8(tbl, vshrq_n_u8(v, 4));
        uint8x16_t lo = vqtbl1q_u8(tbl, vandq_u8(v, mask));
        uint8x16x2_t interleaved = vzipq_u8(hi, lo);
        vst1q_u8(chs.data + 2 * i,      interleaved.val[0]);
        vst1q_u8(chs.data + 2 * i + 16, interleaved.val[1]);
    }
    for (; i < len / 2; ++i) {
        chs.data[2 * i]     = 27 - std::countl_zero(static_cast<unsigned>(data[i]));
        chs.data[2 * i + 1] = std::countr_zero(static_cast<unsigned>(data[i]));
    }
#elif defined(__SSSE3__)
    static const uint8_t decode4_tbl_arr[16] = {
        0,0,1,0, 2,0,0,0, 3,0,0,0, 0,0,0,0};
    const __m128i tbl = _mm_loadu_si128(
        reinterpret_cast<const __m128i *>(decode4_tbl_arr));
    const __m128i mask = _mm_set1_epi8(0x0F);
    size_t i = 0;
    for (; i + 16 <= len / 2; i += 16) {
        __m128i v = _mm_loadu_si128(
            reinterpret_cast<const __m128i *>(data + i));
        __m128i hi = _mm_shuffle_epi8(tbl,
            _mm_and_si128(_mm_srli_epi16(v, 4), mask));
        __m128i lo = _mm_shuffle_epi8(tbl, _mm_and_si128(v, mask));
        __m128i interleaved_lo = _mm_unpacklo_epi8(hi, lo);
        __m128i interleaved_hi = _mm_unpackhi_epi8(hi, lo);
        _mm_storeu_si128(
            reinterpret_cast<__m128i *>(chs.data + 2 * i), interleaved_lo);
        _mm_storeu_si128(
            reinterpret_cast<__m128i *>(chs.data + 2 * i + 16), interleaved_hi);
    }
    for (; i < len / 2; ++i) {
        chs.data[2 * i]     = 27 - std::countl_zero(static_cast<unsigned>(data[i]));
        chs.data[2 * i + 1] = std::countr_zero(static_cast<unsigned>(data[i]));
    }
#else
    for (size_t i = 0; i < len / 2; ++i) {
        chs.data[2 * i]     = 27 - std::countl_zero(static_cast<unsigned>(data[i]));
        chs.data[2 * i + 1] = std::countr_zero(static_cast<unsigned>(data[i]));
    }
#endif
    if (len % 2 == 1) {
        chs.data[len - 1] = 27 - std::countl_zero(static_cast<unsigned>(data[len / 2]));
    }
}
void onehotto2bitDNA_bivec_chs_reverse(const uint8_t *data, size_t len,
                                       charsi &chs) {
    chs.require(len);
#if defined(__ARM_NEON) || defined(__ARM_NEON__)
    static const uint8_t rc_tbl_arr[16] = {
        0,3,2,0, 1,0,0,0, 0,0,0,0, 0,0,0,0};
    const uint8x16_t tbl = vld1q_u8(rc_tbl_arr);
    const uint8x16_t mask = vdupq_n_u8(0x0F);
    size_t half = len / 2;
    size_t i = 0;
    for (; i + 16 <= half; i += 16) {
        uint8x16_t v = vld1q_u8(data + i);
        uint8x16_t lo = vqtbl1q_u8(tbl, vandq_u8(v, mask));
        uint8x16_t hi = vqtbl1q_u8(tbl, vshrq_n_u8(v, 4));
        uint8x16x2_t interleaved = vzipq_u8(lo, hi);
        uint8x16_t r0 = vrev64q_u8(interleaved.val[1]);
        r0 = vcombine_u8(vget_high_u8(r0), vget_low_u8(r0));
        uint8x16_t r1 = vrev64q_u8(interleaved.val[0]);
        r1 = vcombine_u8(vget_high_u8(r1), vget_low_u8(r1));
        size_t out_pos = len - 2 * (i + 16);
        vst1q_u8(chs.data + out_pos,      r0);
        vst1q_u8(chs.data + out_pos + 16, r1);
    }
    size_t mid = half - 1;
    for (; i < half; ++i) {
        chs.data[2 * (mid - i)]     = std::countl_zero(static_cast<unsigned>(data[i])) - 24;
        chs.data[2 * (mid - i) + 1] = 3 - std::countr_zero(static_cast<unsigned>(data[i]));
    }
#elif defined(__SSSE3__)
    static const uint8_t rc_tbl_arr[16] = {
        0,3,2,0, 1,0,0,0, 0,0,0,0, 0,0,0,0};
    const __m128i tbl = _mm_loadu_si128(
        reinterpret_cast<const __m128i *>(rc_tbl_arr));
    const __m128i mask = _mm_set1_epi8(0x0F);
    // Byte-reverse pattern for _mm_shuffle_epi8: reverses 16 bytes
    const __m128i rev_pattern = _mm_set_epi8(
        0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15);
    size_t half = len / 2;
    size_t i = 0;
    for (; i + 16 <= half; i += 16) {
        __m128i v = _mm_loadu_si128(
            reinterpret_cast<const __m128i *>(data + i));
        __m128i lo = _mm_shuffle_epi8(tbl, _mm_and_si128(v, mask));
        __m128i hi = _mm_shuffle_epi8(tbl,
            _mm_and_si128(_mm_srli_epi16(v, 4), mask));
        // Interleave: lo nibble (complement) first, then hi nibble (complement)
        // Note: reverse order (lo, hi) because we're complementing
        __m128i interleaved_lo = _mm_unpacklo_epi8(lo, hi);
        __m128i interleaved_hi = _mm_unpackhi_epi8(lo, hi);
        // Reverse both 16-byte vectors
        __m128i r0 = _mm_shuffle_epi8(interleaved_hi, rev_pattern);
        __m128i r1 = _mm_shuffle_epi8(interleaved_lo, rev_pattern);
        size_t out_pos = len - 2 * (i + 16);
        _mm_storeu_si128(
            reinterpret_cast<__m128i *>(chs.data + out_pos), r0);
        _mm_storeu_si128(
            reinterpret_cast<__m128i *>(chs.data + out_pos + 16), r1);
    }
    size_t mid = half - 1;
    for (; i < half; ++i) {
        chs.data[2 * (mid - i)]     = std::countl_zero(static_cast<unsigned>(data[i])) - 24;
        chs.data[2 * (mid - i) + 1] = 3 - std::countr_zero(static_cast<unsigned>(data[i]));
    }
#else
    size_t mid = len / 2 - 1;
    for (size_t i = 0; i < len / 2; ++i) {
        chs.data[2 * (mid - i)]     = std::countl_zero(static_cast<unsigned>(data[i])) - 24;
        chs.data[2 * (mid - i) + 1] = 3 - std::countr_zero(static_cast<unsigned>(data[i]));
    }
#endif
    if (len % 2 == 1) {
        chs.data[0] = std::countl_zero(static_cast<unsigned>(data[len / 2])) - 24;
    }
}


// TESTS //

TEST_CASE("Nibbles to Crumbs") {
    gtl::vector<uint8_t> bytes;
    gtl::flat_hash_map<uint8_t, std::pair<uint8_t, uint8_t>> expected;
    static constexpr int tab[] = {4, 0, 1, 4, 2, 4, 4, 4,
                                  3, 4, 4, 4, 4, 4, 4, 4};

    for (int high : {1, 2, 4, 8}) {
        for (int low : {1, 2, 4, 8}) {
            uint8_t byte = (high << 4) | low;
            bytes.push_back(byte);
            expected[byte] = {tab[high], tab[low]};
        }
    }
    SUBCASE("Even reverse") {
        charsi c;
        c.require(2 * bytes.size());
        onehotto2bitDNA_bivec_chs_reverse(bytes.data(), 2 * bytes.size(), c);
        [[maybe_unused]] size_t bs = bytes.size() - 1;
        for (unsigned i = 0; i < bytes.size(); ++i) {
            CHECK(3 - expected[bytes[bs - i]].first == c.data[2 * i]);
            CHECK(3 - expected[bytes[bs - i]].second == c.data[2 * i + 1]);
        }
    }
    bytes.push_back((1 << 4));
    expected[bytes.back()] = {tab[1], 0};
}

// }

TEST_CASE("K-mer view on Crumbs") {
    gtl::vector<uint8_t> bytes{{0, 1, 2, 3, 0, 1, 2, 3}};
    SUBCASE("K==4") {
        static constexpr int K = 4;
        KmerViewOnInts<uint64_t, 2> kview(K);
        size_t i = kview.fill(bytes.data());
        [[maybe_unused]] uint32_t expected =
            bytes[0] << 6 | bytes[1] << 4 | bytes[2] << 2 | bytes[3];
        for (; i < bytes.size(); ++i) {
            CHECK(expected == kview.get());
            kview.move_get(bytes.data(), i);
            expected <<= 2;
            expected |= bytes[i];
            expected &= 0xff;
        }
    }
    SUBCASE("K==5") {
        static constexpr int K = 5;
        KmerViewOnInts<uint64_t, 2> kview(K);
        size_t i = kview.fill(bytes.data());
        [[maybe_unused]] uint32_t expected = bytes[0] << 8 | bytes[1] << 6 | bytes[2] << 4 |
                            bytes[3] << 2 | bytes[4];
        for (; i < bytes.size(); ++i) {
            CHECK(expected == kview.get());
            kview.move_get(bytes.data(), i);
            expected <<= 2;
            expected |= bytes[i];
            expected &= 0x3ff;
        }
    }
}

TEST_CASE("Character table on kmerview on bytes") {
    using kv = KmerViewOnBytes<16, uint32_t, 2>;
    (void) kv::k;
    CHECK_EQ(kv::tab['A'], 0);
    CHECK_EQ(kv::tab['G'], 2);
    CHECK_EQ(kv::tab['C'], 1);
    CHECK_EQ(kv::tab['T'], 3);

    CHECK_EQ(kv::tab['a'], 0);
    CHECK_EQ(kv::tab['g'], 2);
    CHECK_EQ(kv::tab['c'], 1);
    CHECK_EQ(kv::tab['t'], 3);

    CHECK_EQ(kv::tab['u'], 3);
    CHECK_EQ(kv::tab['U'], 3);

    CHECK_EQ(kv::tab['N'], 4);
    CHECK_EQ(kv::tab['Z'], 4);
}

TEST_CASE("BAM Read Parsing") {
    std::string input = "../tests/test_single.bam";
    samFile *fp = sam_open(input.c_str(), "r");
    if (fp == NULL) {
        FAIL("Failed to open BAM file:", input);
    }

    bam_hdr_t *hdr = sam_hdr_read(fp);
    bam1_t *rec = bam_init1();
    if (sam_read1(fp, hdr, rec) <= 0) {
        FAIL("Failed to read BAM record");
    }

    char *seq_gt = (char *)malloc(rec->core.l_qseq + 1);
    for (int i = 0; i < rec->core.l_qseq; ++i) {
        seq_gt[i] = seq_nt16_str[bam_seqi(bam_get_seq(rec), i)];
    }
    seq_gt[rec->core.l_qseq] = '\0';
}

TEST_CASE("BAM Read Parsing With Crumbs") {
    std::string input = "../tests/test_single.bam";
    samFile *fp = sam_open(input.c_str(), "r");
    if (fp == NULL) {
        FAIL("Failed to open BAM file:", input);
    }

    bam_hdr_t *hdr = sam_hdr_read(fp);
    bam1_t *rec = bam_init1();
    if (sam_read1(fp, hdr, rec) <= 0) {
        FAIL("Failed to read BAM record");
    }

    char *seq_gt = (char *)malloc(rec->core.l_qseq + 1);
    for (int i = 0; i < rec->core.l_qseq; ++i) {
        seq_gt[i] = seq_nt16_str[bam_seqi(bam_get_seq(rec), i)];
    }
    seq_gt[rec->core.l_qseq] = '\0';
    SUBCASE(
        "Converting Nibble-encoded data to crumbs and decoding it to bases") {
        int seqlen = strlen(seq_gt);
        srand(time(NULL));
        int N_RANDOM_CASES = 50;

        gtl::vector<CrumbledSubSeq64> copied;

        gtl::vector<size_t> froms;
        gtl::vector<size_t> tos;
        for (int casei = 0; casei < N_RANDOM_CASES; ++casei) {
            int to = rand() % (seqlen - 100) + 100;
            int back = 100 + (rand() % 12);
            int from = (to > back) ? (to - back) : 0;
            copied.emplace_back(bam_get_seq(rec), from, to, 0,
                                CrumbledSubSeq64::from_nibbles{});
            froms.push_back(from);
            tos.push_back(to);
        }
        char buffer[1024 * 1024];
        for (int casei = 0; casei < N_RANDOM_CASES; ++casei) {
            auto &seq = copied[casei];
            size_t from = froms[casei];
            size_t to = tos[casei];
            [[maybe_unused]] char *starting_from = seq.convert_to_ascii(buffer);
            INFO("from: ", from, " to: ", to, "seq.len: ", seq.len);
            CHECK_EQ(std::string{seq_gt}.substr(from, to - from),
                     std::string{starting_from});
            CrumbledSubSeq64 fr_chr_cr{seq_gt, from, to, 0,
                                       CrumbledSubSeq64::from_ascii{}};
            starting_from = fr_chr_cr.convert_to_ascii(buffer);
            CHECK_EQ(std::string{seq_gt}.substr(from, to - from),
                     std::string{starting_from});
        }
    }

    SUBCASE(
        "Compare two windowing functions to see if they produce the same "
        "lcmers") {
        lps str1 = {0, 0, 0, nullptr};
        lps str2 = {0, 0, 0, nullptr};
        using SYNCMER_HASHER = syncmer::MurMur;
        using SyncmerType =
            syncmer::GenericSyncmerRule<uint32_t, SYNCMER_HASHER>;
        decltype(asm_config::syncmer) spara;
        spara.k = 8;
        spara.density = 0.05;
        SyncmerType fnctr(spara, 11);
        charsi chs;
        pos_scratch pos_buf;
        lcp_read lr;
        parse_cores_from_bam_record(bam_get_seq(rec), rec->core.l_qseq, 3,
                                    &str1, fnctr, chs, lr, pos_buf);
        gtl::vector<CrumbledSubSeq64> ss;
        ss.emplace_back(bam_get_seq(rec), 0, rec->core.l_qseq, 0,
                        CrumbledSubSeq64::from_nibbles{});

        lcp_read lr2;
        init_lps_crumbled(&str2, ss, fnctr);
        // printf("\n");
        lps_deepen(&str2, 2);
        for (int i = 0; i < str2.size; ++i) {
            lr2.cores.emplace_back(str2.cores[i].label);
            lr2.pos.emplace_back(str2.cores[i].start, str2.cores[i].end);
        }
        CHECK(lr.cores.size() == lr2.cores.size());
    }

    SUBCASE("LcpSubsample re-LCP consistency (parse vs crumbled)") {
        using SYNCMER_HASHER = syncmer::MurMur;
        using SyncmerType =
            syncmer::GenericSyncmerRule<uint32_t, SYNCMER_HASHER>;
        decltype(asm_config::syncmer) spara;
        spara.k = 8;
        spara.lcp_low_level = 2;
        spara.density = 0.3;
        SyncmerType fnctr(spara, 11);
        const int final_level = 3;

        charsi chs;
        pos_scratch pos_buf;
        lcp_read lr;
        lps str1 = {0, 0, 0, nullptr};
        parse_cores_from_bam_record(bam_get_seq(rec), rec->core.l_qseq,
                                    final_level, &str1, fnctr, chs, lr, pos_buf);

        gtl::vector<CrumbledSubSeq64> ss;
        ss.emplace_back(bam_get_seq(rec), 0, rec->core.l_qseq, 0,
                        CrumbledSubSeq64::from_nibbles{});
        lps str2 = {0, 0, 0, nullptr};
        lcp_read lr2;
        init_lps_crumbled(&str2, ss, fnctr);
        lps_deepen(&str2, final_level);
        for (int i = 0; i < str2.size; ++i) {
            lr2.cores.emplace_back(str2.cores[i].label);
            lr2.pos.emplace_back(str2.cores[i].start, str2.cores[i].end);
        }
        // Same sequence + same params must yield identical cores (labels and
        // positions), else node re-parses cannot match the DBG index -> few hits.
        REQUIRE(lr.cores.size() == lr2.cores.size());
        for (size_t i = 0; i < lr.cores.size(); ++i) {
            CHECK(lr.cores[i] == lr2.cores[i]);
            CHECK(lr.pos[i].first == lr2.pos[i].first);
            CHECK(lr.pos[i].second == lr2.pos[i].second);
        }
        free_lps(&str1);
        free_lps(&str2);
    }

    bam_hdr_destroy(hdr);
    sam_close(fp);
    free(seq_gt);
}

// Reference scalar-only two-pass: decode then RLE.
// Used as ground truth for testing onehot_to_rle.
[[maybe_unused]] static int64_t reference_onehot_to_rle(const uint8_t *data, int64_t len,
                                        uint8_t *rle_data, uint32_t *rle_pos) {
    if (len <= 0) return 0;
    // Scalar decode table: BAM one-hot nibble → 2-bit
    static constexpr uint8_t tbl[16] = {
        0,0,1,0, 2,0,0,0, 3,0,0,0, 0,0,0,0};
    // Decode all bases
    gtl::vector<uint8_t> decoded(len);
    for (int64_t i = 0; i < len / 2; ++i) {
        decoded[2 * i]     = tbl[data[i] >> 4];
        decoded[2 * i + 1] = tbl[data[i] & 0x0F];
    }
    if (len % 2 == 1) {
        decoded[len - 1] = tbl[data[len / 2] >> 4];
    }
    // Scalar RLE compaction
    rle_data[0] = decoded[0];
    rle_pos[0] = 0;
    int64_t out = 1;
    for (int64_t i = 1; i < len; ++i) {
        if (decoded[i] != decoded[i - 1]) {
            rle_data[out] = decoded[i];
            rle_pos[out] = static_cast<uint32_t>(i);
            ++out;
        }
    }
    return out;
}

TEST_CASE("onehot_to_rle correctness") {
    // Encode a synthetic BAM nibble sequence from known bases
    auto encode_bam = [](const gtl::vector<uint8_t> &bases) {
        // bases: 0=A, 1=C, 2=G, 3=T → BAM nibbles: 1,2,4,8
        static constexpr uint8_t to_bam[4] = {1, 2, 4, 8};
        int64_t len = bases.size();
        gtl::vector<uint8_t> bam_data((len + 1) / 2);
        for (int64_t i = 0; i < len / 2; ++i) {
            bam_data[i] = (to_bam[bases[2*i]] << 4) | to_bam[bases[2*i+1]];
        }
        if (len % 2 == 1) {
            bam_data[len / 2] = to_bam[bases[len - 1]] << 4;
        }
        return bam_data;
    };

    auto run_test = [&](const gtl::vector<uint8_t> &bases, [[maybe_unused]] const char *label) {
        INFO(label, " len=", bases.size());
        auto bam_data = encode_bam(bases);
        int64_t len = bases.size();

        gtl::vector<uint8_t> rle_ref(len), rle_test(len);
        gtl::vector<uint32_t> pos_ref(len), pos_test(len);

        int64_t n_ref  = reference_onehot_to_rle(bam_data.data(), len,
                                                  rle_ref.data(), pos_ref.data());
        int64_t n_test = onehot_to_rle(bam_data.data(), len,
                                       rle_test.data(), pos_test.data());

        CHECK(n_ref == n_test);
        for (int64_t i = 0; i < std::min(n_ref, n_test); ++i) {
            INFO("  index ", i);
            CHECK(rle_ref[i] == rle_test[i]);
            CHECK(pos_ref[i] == pos_test[i]);
        }
    };

    SUBCASE("Short sequences") {
        run_test({0}, "single A");
        run_test({0, 1}, "AC");
        run_test({0, 0}, "AA");
        run_test({0, 1, 2, 3}, "ACGT");
        run_test({3, 3, 3, 3}, "TTTT");
        run_test({0, 0, 1, 1, 2, 2, 3, 3}, "AACCGGTT");
    }

    SUBCASE("Odd length") {
        run_test({0, 1, 2}, "ACG (odd)");
        run_test({0, 1, 2, 3, 0}, "ACGTA (odd)");
        gtl::vector<uint8_t> odd33(33, 0);
        odd33[16] = 1; odd33[32] = 2;
        run_test(odd33, "33 bases (odd, spans SIMD boundary)");
    }

    SUBCASE("Homopolymer run (all same, exercises fast-path skip)") {
        run_test(gtl::vector<uint8_t>(64, 2), "64x G");
        run_test(gtl::vector<uint8_t>(100, 0), "100x A");
    }

    SUBCASE("Alternating (max changes)") {
        gtl::vector<uint8_t> alt(128);
        for (size_t i = 0; i < alt.size(); ++i) alt[i] = i % 4;
        run_test(alt, "128 alternating ACGT");
    }

    SUBCASE("Exactly 32 bases (one SIMD iteration)") {
        gtl::vector<uint8_t> v(32);
        for (size_t i = 0; i < 32; ++i) v[i] = (i / 8) % 4;
        run_test(v, "32 bases, 4 blocks of 8");
    }

    SUBCASE("Boundary: 31, 32, 33 bases") {
        for (int n : {31, 32, 33}) {
            gtl::vector<uint8_t> v(n);
            for (int i = 0; i < n; ++i) v[i] = (i * 7) % 4;
            run_test(v, (std::string("len=") + std::to_string(n)).c_str());
        }
    }

    SUBCASE("Large random sequence") {
        srand(42);
        for (int trial = 0; trial < 10; ++trial) {
            int len = 500 + rand() % 15000;
            gtl::vector<uint8_t> bases(len);
            for (int i = 0; i < len; ++i) bases[i] = rand() % 4;
            run_test(bases, (std::string("random trial ") +
                             std::to_string(trial)).c_str());
        }
    }

    SUBCASE("Real BAM record") {
        // Try multiple paths for test data
        const char *paths[] = {
            "../tests/test_single.bam",
            "../../data/remth.bam",
            "../../../data/remth.bam",
            "/Users/fatihkaraoglanoglu/Workspace/data/remth.bam",
            nullptr
        };
        samFile *fp = nullptr;
        std::string input;
        for (int p = 0; paths[p]; ++p) {
            input = paths[p];
            fp = sam_open(input.c_str(), "r");
            if (fp) break;
        }
        if (fp == NULL) {
            MESSAGE("Skipping real BAM test: no test BAM found");
            return;
        }
        bam_hdr_t *hdr = sam_hdr_read(fp);
        bam1_t *rec = bam_init1();
        int count = 0;
        while (sam_read1(fp, hdr, rec) > 0 && count < 20) {
            const uint8_t *data = bam_get_seq(rec);
            int64_t len = rec->core.l_qseq;
            INFO("BAM record ", count, " len=", len);

            gtl::vector<uint8_t> rle_ref(len), rle_test(len);
            gtl::vector<uint32_t> pos_ref(len), pos_test(len);

            int64_t n_ref  = reference_onehot_to_rle(data, len,
                                rle_ref.data(), pos_ref.data());
            int64_t n_test = onehot_to_rle(data, len,
                                rle_test.data(), pos_test.data());

            CHECK(n_ref == n_test);
            for (int64_t i = 0; i < std::min(n_ref, n_test); ++i) {
                if (rle_ref[i] != rle_test[i] || pos_ref[i] != pos_test[i]) {
                    INFO("  mismatch at rle index ", i,
                         ": ref base=", (int)rle_ref[i], " pos=", pos_ref[i],
                         " vs test base=", (int)rle_test[i], " pos=", pos_test[i]);
                    CHECK(rle_ref[i] == rle_test[i]);
                    CHECK(pos_ref[i] == pos_test[i]);
                    break;  // stop at first mismatch for readable output
                }
            }
            ++count;
        }
        bam_destroy1(rec);
        bam_hdr_destroy(hdr);
        sam_close(fp);
    }
}
