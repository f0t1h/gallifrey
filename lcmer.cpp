#include "lcmer.h"

#include <bit>
#include <cstdint>

#include "asm.h"
#include "ext/doctest.h"
#include "ext/lcptools/lcptools_ho.h"
#include "fastx.h"
#include "gtl/phmap.hpp"
#include "htslib/hts.h"
#include "htslib/sam.h"
#include "nibble2base.hpp"
#include "syncmer.h"

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
        // Interleave lo nibble first, then hi nibble (reversed order for the
        // complement).
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

#ifndef LASM_TEST_DATA_DIR
#define LASM_TEST_DATA_DIR "tests"  // set by CMake for lasm_tests
#endif
// Single unmapped HiFi read (17.4 kb), no reference sequences in the header.
[[maybe_unused]] static const char *const kTestSingleBam =
    LASM_TEST_DATA_DIR "/test_single.bam";

TEST_CASE("BAM Read Parsing") {
    samFile *fp = sam_open(kTestSingleBam, "r");
    REQUIRE_MESSAGE(fp != NULL, "Failed to open BAM file: ", kTestSingleBam);

    bam_hdr_t *hdr = sam_hdr_read(fp);
    REQUIRE(hdr != NULL);
    bam1_t *rec = bam_init1();
    REQUIRE(sam_read1(fp, hdr, rec) > 0);
    CHECK((rec->core.flag & BAM_FUNMAP) != 0);
    CHECK_EQ(rec->core.l_qseq, 17430);

    std::string seq_gt(rec->core.l_qseq, '\0');
    for (int i = 0; i < rec->core.l_qseq; ++i) {
        seq_gt[i] = seq_nt16_str[bam_seqi(bam_get_seq(rec), i)];
    }
    CHECK_EQ(seq_gt.find_first_not_of("ACGT"), std::string::npos);
    CHECK_EQ(sam_read1(fp, hdr, rec), -1);  // exactly one record

    bam_destroy1(rec);
    bam_hdr_destroy(hdr);
    sam_close(fp);
}

TEST_CASE("BAM Read Parsing With Crumbs") {
    samFile *fp = sam_open(kTestSingleBam, "r");
    REQUIRE_MESSAGE(fp != NULL, "Failed to open BAM file: ", kTestSingleBam);

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
        if (seqlen <= 112) {
            return;  // record too short for the random windows below
        }
        srand(42);  // deterministic: this test must be reproducible
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
        core_buf str1;
        core_buf str2;
        using SYNCMER_HASHER = syncmer::MurMur;
        using SyncmerType =
            syncmer::GenericSubmerRule<uint32_t, SYNCMER_HASHER>;
        decltype(asm_config::syncmer) spara;
        spara.k = 8;
        spara.smer.s = 3;
        spara.smer.t = 2;
        SyncmerType fnctr(spara, 11);
        charsi chs;
        pos_scratch pos_buf;
        lcp_read lr;
        constexpr int LCP_LEVEL = 3;
        parse_cores_from_nibbles(bam_get_seq(rec), rec->core.l_qseq, LCP_LEVEL,
                                 &str1, fnctr, chs, lr, pos_buf);
        gtl::vector<CrumbledSubSeq64> ss;
        ss.emplace_back(bam_get_seq(rec), 0, rec->core.l_qseq, 0,
                        CrumbledSubSeq64::from_nibbles{});

        lcp_read lr2;
        init_lps_crumbled(&str2, ss, fnctr);
        // Deepen to the same level as the read path (parse_cores_from_nibbles
        // deepens internally); comparing different levels is meaningless.
        str2.deepen(fnctr, LCP_LEVEL);
        str2.append_to(lr2);
        REQUIRE(lr.cores.size() == lr2.cores.size());
        for (size_t i = 0; i < lr.cores.size(); ++i) {
            INFO("core index ", i);
            CHECK(lr.cores[i] == lr2.cores[i]);
            CHECK(lr.pos[i].first == lr2.pos[i].first);
            CHECK(lr.pos[i].second == lr2.pos[i].second);
        }
    }

    bam_destroy1(rec);
    bam_hdr_destroy(hdr);
    sam_close(fp);
    free(seq_gt);
}

TEST_CASE("BAM and crumbled parsers agree on homopolymer-rich sequence") {
    std::string bases;
    {
        const char *runs[] = {"AAAAAAA",  "CC",     "GGGGGGGGGG", "T",
                              "ACGTACGT", "CCCCCCCC", "GG",        "TTTTTTT",
                              "AACCGGTTAACCGGTT", "AAAA", "CCCC",  "GGGG",
                              "TTTT",     "ACGT"};
        for (int r = 0; r < 12; ++r) {
            for (const char *run : runs) bases += run;
        }
    }
    const int64_t len = static_cast<int64_t>(bases.size());

    // BAM one-hot nibbles: A=1, C=2, G=4, T=8 (high nibble = even index).
    gtl::vector<uint8_t> nib((len + 1) / 2, 0);
    auto code = [](char c) {
        switch (c) {
            case 'A': return 1;
            case 'C': return 2;
            case 'G': return 4;
            default: return 8;
        }
    };
    for (int64_t i = 0; i < len; ++i) {
        if (i % 2 == 0) {
            nib[i / 2] |= code(bases[i]) << 4;
        } else {
            nib[i / 2] |= code(bases[i]);
        }
    }

    using SyncmerType = syncmer::GenericSubmerRule<uint32_t, syncmer::MurMur>;
    decltype(asm_config::syncmer) spara;
    spara.k = 6;       // small k so homopolymer windows are frequent
    spara.smer.s = 2;  // closed syncmer: selects ~2/5 of the windows
    SyncmerType fnctr(spara, 7);

    core_buf str1;
    core_buf str2;
    charsi chs;
    pos_scratch pos_buf;
    lcp_read lr;
    parse_cores_from_nibbles(nib.data(), len, 1, &str1, fnctr, chs, lr,
                             pos_buf);

    gtl::vector<CrumbledSubSeq64> ss;
    ss.emplace_back(bases.c_str(), 0, bases.size(), 0,
                    CrumbledSubSeq64::from_ascii{});
    lcp_read lr2;
    init_lps_crumbled(&str2, ss, fnctr);
    str2.append_to(lr2);

    REQUIRE(!lr.cores.empty());
    REQUIRE(lr.cores.size() == lr2.cores.size());
    for (size_t i = 0; i < lr.cores.size(); ++i) {
        INFO("core ", i);
        CHECK(lr.cores[i] == lr2.cores[i]);
        CHECK(lr.pos[i].first == lr2.pos[i].first);
        CHECK(lr.pos[i].second == lr2.pos[i].second);
    }
    // True-span sanity: every span lies inside the sequence and covers at
    // least k bases (an RLE k-mer spans >= k original bases).
    for (size_t i = 0; i < lr.pos.size(); ++i) {
        INFO("span ", i);
        CHECK(lr.pos[i].second <= static_cast<uint64_t>(len));
        CHECK(lr.pos[i].second - lr.pos[i].first >=
              static_cast<uint64_t>(spara.k));
    }
}

TEST_CASE("core_buf deepening matches lcptools lps_deepen on a read") {
    samFile *fp = sam_open(kTestSingleBam, "r");
    REQUIRE_MESSAGE(fp != NULL, "Failed to open BAM file: ", kTestSingleBam);
    bam_hdr_t *hdr = sam_hdr_read(fp);
    bam1_t *rec = bam_init1();
    REQUIRE(sam_read1(fp, hdr, rec) > 0);
    const uint8_t *seq = bam_get_seq(rec);
    const int64_t len = rec->core.l_qseq;

    using SyncmerType = syncmer::GenericSubmerRule<uint32_t, syncmer::MurMur>;
    for (int k : {6, 11, 16}) {
        decltype(asm_config::syncmer) spara;
        spara.k = k;
        spara.apply_default_rule();
        SyncmerType fnctr(spara, 1453);
        charsi chs;
        pos_scratch pos_buf;
        chs.require(len + kOnehotRlePad);
        pos_buf.require(len + kOnehotRlePad);
        const int64_t rle_len = onehot_to_rle(seq, len, chs.data, pos_buf.data);
        for (int level : {2, 3, 4, 5, 6}) {
            INFO("k ", k, " level ", level);
            core_buf got;
            parse1_rle_syncmers(chs.data, pos_buf.data, rle_len, len, &got,
                                fnctr);
            REQUIRE(got.size > 100);
            // Reference: the same level-1 k-mer cores as struct core.
            struct lps ref;
            ref.level = 1;
            ref.size = static_cast<int>(got.size);
            ref.capacity = ref.size;
            ref.cores = (struct core *)malloc(ref.size * sizeof(struct core));
            for (int i = 0; i < ref.size; ++i)
                init_core4(&ref.cores[i], 2 * k,
                           0x8000000000000000ULL | got.label[i], got.label[i],
                           got.start(i), got.end(i));
            lps_deepen(&ref, level);
            got.deepen(fnctr, level);
            REQUIRE(got.level == ref.level);
            REQUIRE(got.size == ref.size);
            for (int i = 0; i < ref.size; ++i) {
                INFO("core ", i);
                CHECK(got.label[i] == ref.cores[i].label);
                CHECK(got.start(i) == ref.cores[i].start);
                CHECK(got.end(i) == ref.cores[i].end);
                CHECK(got.bit_rep[i] == ref.cores[i].bit_rep);
                CHECK(got.bit_size[i] == ref.cores[i].bit_size);
            }
            free(ref.cores);
        }
    }
    bam_destroy1(rec);
    bam_hdr_destroy(hdr);
    sam_close(fp);
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

        // std::vector: its max_size() check bounds the size, so GCC does not
        // warn about a size_t overflow in gtl::vector's fill
        // (-Wstringop-overflow).
        const size_t cap = bases.size() + kOnehotRlePad;
        gtl::vector<uint8_t> rle_ref(bases.size());
        gtl::vector<uint32_t> pos_ref(bases.size());
        std::vector<uint8_t> rle_test(cap);
        std::vector<uint32_t> pos_test(cap);

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
        samFile *fp = sam_open(kTestSingleBam, "r");
        REQUIRE_MESSAGE(fp != NULL, "Failed to open BAM file: ",
                        kTestSingleBam);
        bam_hdr_t *hdr = sam_hdr_read(fp);
        bam1_t *rec = bam_init1();
        int count = 0;
        while (sam_read1(fp, hdr, rec) > 0 && count < 20) {
            const uint8_t *data = bam_get_seq(rec);
            int64_t len = rec->core.l_qseq;
            INFO("BAM record ", count, " len=", len);

            gtl::vector<uint8_t> rle_ref(len), rle_test(len + kOnehotRlePad);
            gtl::vector<uint32_t> pos_ref(len), pos_test(len + kOnehotRlePad);

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

// The FASTX path exists only because it is cheaper than htslib; it is correct
// only if it produces the exact same packed sequence -- and therefore the exact
// same core stream -- as the BAM path for the same bases.
TEST_CASE("FASTX reader agrees with the htslib path") {
    samFile *fp = sam_open(kTestSingleBam, "r");
    REQUIRE_MESSAGE(fp != NULL, "Failed to open BAM file: ", kTestSingleBam);
    bam_hdr_t *hdr = sam_hdr_read(fp);
    bam1_t *rec = bam_init1();
    REQUIRE(sam_read1(fp, hdr, rec) > 0);
    const int len = rec->core.l_qseq;
    const uint8_t *bam_seq = bam_get_seq(rec);

    std::string bases(static_cast<size_t>(len), '\0');
    for (int i = 0; i < len; ++i) {
        bases[static_cast<size_t>(i)] = seq_nt16_str[bam_seqi(bam_seq, i)];
    }

    auto wrap = [&](const std::string &seq, size_t cols, const char *eol) {
        std::string out;
        for (size_t i = 0; i < seq.size(); i += cols) {
            out += seq.substr(i, cols);
            out += eol;
        }
        return out;
    };
    auto write_file = [](const std::string &path, const std::string &content) {
        FILE *f = fopen(path.c_str(), "w");
        REQUIRE(f != nullptr);
        [[maybe_unused]] const size_t wrote =
            fwrite(content.data(), 1, content.size(), f);
        fclose(f);
        REQUIRE(wrote == content.size());
    };

    // Private scratch directory so concurrent test runs cannot collide.
    char dir_template[] = "/tmp/lasm_fastx_test.XXXXXX";
    [[maybe_unused]] const char *made = mkdtemp(dir_template);
    REQUIRE(made != nullptr);
    const std::string dir = std::string(dir_template) + "/";
    const std::string fa = dir + "lasm_fastx_test.fa";
    const std::string fa_crlf = dir + "lasm_fastx_test_crlf.fa";
    const std::string fq = dir + "lasm_fastx_test.fq";

    // Two records per file: the reader must not fuse or drop them.
    write_file(fa, ">r1 comment\n" + wrap(bases, 60, "\n") + ">r2\n" +
                       wrap(bases, 60, "\n"));
    // CRLF + lowercase bases must decode identically.
    std::string lower = bases;
    for (char &c : lower) c = static_cast<char>(tolower(c));
    write_file(fa_crlf, ">r1\r\n" + wrap(lower, 70, "\r\n") + ">r2\r\n" +
                            wrap(lower, 70, "\r\n"));
    // Quality line made entirely of '@' -- the classic FASTQ parser trap.
    const std::string qual(bases.size(), '@');
    write_file(fq, "@r1\n" + bases + "\n+\n" + qual + "\n@r2\n" + bases +
                       "\n+r2\n" + qual + "\n");
    // Wrapped FASTQ with CRLF: every quality line starts with '@', and the
    // quality span crosses many line breaks and block boundaries.
    const std::string fq_crlf = dir + "lasm_fastx_test_crlf.fq";
    write_file(fq_crlf, "@r1\r\n" + wrap(bases, 80, "\r\n") + "+\r\n" +
                            wrap(qual, 80, "\r\n") + "@r2\r\n" +
                            wrap(bases, 80, "\r\n") + "+\r\n" +
                            wrap(qual, 80, "\r\n"));

    auto check_file = [&](const std::string &path, size_t block_size) {
        INFO("path=", path, " block=", block_size);
        fastx_reader reader(block_size);
        [[maybe_unused]] const bool opened = reader.open(path, false, 1);
        REQUIRE(opened);
        seq_record sr;
        for (int r = 0; r < 2; ++r) {
            INFO("record ", r);
            REQUIRE(reader.next(sr));
            CHECK_EQ(sr.len, static_cast<uint32_t>(len));
            CHECK_EQ(memcmp(sr.nibbles, bam_seq, static_cast<size_t>(len) / 2),
                     0);
        }
        CHECK_FALSE(reader.next(sr));
        CHECK_FALSE(reader.failed());
        free(sr.nibbles);
    };

    // Block sizes far below a record length exercise every state transition
    // across a buffer boundary.
    for (size_t block :
         {size_t(7), size_t(64), size_t(4096), size_t(1) << 20}) {
        check_file(fa, block);
        check_file(fa_crlf, block);
        check_file(fq, block);
        check_file(fq_crlf, block);
    }

    SUBCASE("non-ACGT bases and mid-window whitespace decode as htslib does") {
        // Exercises the scalar fallback inside the vectorized packer: N/IUPAC
        // codes and odd-length lines that push the packer off even offsets.
        const std::string dirty =
            "ACGTacgtNNNNRYKM"
            "ACGTACGTACGTACGTA"
            "nnnACGT"
            "W";
        const std::string path = dir + "lasm_fastx_dirty.fa";
        write_file(path, ">d\n" + dirty.substr(0, 16) + "\n" +
                             dirty.substr(16, 17) + "\n" + dirty.substr(33) +
                             "\n");
        fastx_reader reader(13);
        [[maybe_unused]] const bool opened = reader.open(path, false, 1);
        REQUIRE(opened);
        seq_record sr;
        REQUIRE(reader.next(sr));
        REQUIRE_EQ(sr.len, static_cast<uint32_t>(dirty.size()));
        for (size_t i = 0; i < dirty.size(); ++i) {
            INFO("base ", i, " '", dirty[i], "'");
            CHECK_EQ(bam_seqi(sr.nibbles, i),
                     seq_nt16_table[static_cast<uint8_t>(dirty[i])]);
        }
        CHECK_FALSE(reader.next(sr));
        free(sr.nibbles);
        remove(path.c_str());
    }

    SUBCASE("core streams are identical") {
        using SYNCMER_HASHER = syncmer::MurMur;
        using SyncmerType =
            syncmer::GenericSubmerRule<uint32_t, SYNCMER_HASHER>;
        syncmer_params spara;
        spara.k = 15;
        spara.smer.s = 5;
        spara.smer.t = 5;
        SyncmerType fnctr(spara, 11);

        fastx_reader reader;
        [[maybe_unused]] const bool opened = reader.open(fq, false, 1);
        REQUIRE(opened);
        seq_record sr;
        REQUIRE(reader.next(sr));

        core_buf str1;
        core_buf str2;
        charsi chs1, chs2;
        pos_scratch pos1, pos2;
        lcp_read lr_bam, lr_fastx;
        constexpr int LCP_LEVEL = 3;
        parse_cores_from_nibbles(bam_seq, static_cast<size_t>(len), LCP_LEVEL,
                                 &str1, fnctr, chs1, lr_bam, pos1);
        parse_cores_from_nibbles(sr.nibbles, sr.len, LCP_LEVEL, &str2, fnctr,
                                 chs2, lr_fastx, pos2);
        REQUIRE(lr_bam.cores.size() > 0);
        REQUIRE_EQ(lr_fastx.cores.size(), lr_bam.cores.size());
        CHECK_EQ(memcmp(lr_fastx.cores.data(), lr_bam.cores.data(),
                        lr_bam.cores.size() * sizeof(uint32_t)),
                 0);
        CHECK(lr_fastx.pos == lr_bam.pos);
        free(sr.nibbles);
    }

    remove(fa.c_str());
    remove(fa_crlf.c_str());
    remove(fq.c_str());
    remove(fq_crlf.c_str());
    rmdir(dir_template);
    bam_destroy1(rec);
    bam_hdr_destroy(hdr);
    sam_close(fp);
}
