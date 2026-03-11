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
    cr->bit_rep = CL(uint128_t){0x8000000000000000, cr->label};
    cr->bit_size = 2 * distance;
}

// syncmer::syncmer_precomputed<
//     syncmer::SYNCMER_RULE<SYNCMER_K, SYNCMER_S, SYNCMER_T>>
//     table;

// uint64_t parse1_bam_syncmers_old(const uint8_t *data, int64_t len,
//                                  struct core *cores) {
//     int core_index = 0;
//     KmerViewOnNibbles<SYNCMER_K, uint32_t, 2> kmer_view;
//     int64_t it1 = kmer_view.fill(data);
//     // find lcp cores
//     for (; it1 + SYNCMER_K < len; it1++) {
//         // printf("%d\t%d\n", kmer_view.get(), table[kmer_view.get()]);
//         if (table[kmer_view.get()]) {
//             init_core_kmer<uint32_t, 2>(&(cores[core_index]),
//             kmer_view.get(),
//                                         it1, it1 + SYNCMER_K);
//             core_index++;
//         }
//         kmer_view.move_get(data, it1);
//         // cout << "None\n";
//     }
//     // cout << "CORESSS" << it1 << "\t" << it2 << "\n";
//     return core_index;
// }
// uint64_t parse1_kseq_syncmers(const char *data, int64_t len,
//                               struct core *cores) {
//     int core_index = 0;
//     KmerViewOnBytes<SYNCMER_K, uint32_t, 2> kmer_view;
//     int64_t it1 = kmer_view.fill(data);
//     // find lcp cores
//     for (; it1 + SYNCMER_K < len; it1++) {
//         // printf("%d\t%d\n", kmer_view.get(), table[kmer_view.get()]);
//         if (table[kmer_view.get()]) {
//             init_core_kmer<uint32_t, 2>(&(cores[core_index]),
//             kmer_view.get(),
//                                         it1, it1 + SYNCMER_K);
//             core_index++;
//         }
//         kmer_view.move_get(data, it1);
//         // cout << "None\n";
//     }
//     // cout << "CORESSS" << it1 << "\t" << it2 << "\n";
//     return core_index;
// }

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
            // // check if there is any SSEQ cores left behind
            // if (it2 < it1 && last_invalid_char_index < it2 - 1) {
            //   init_core1bam(&(cores[core_index]), data, it1 - it2 + 2, it2 -
            //   1,
            //                 it1 + 1);
            //   core_index++;
            // }

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
// void init_lps_bam_old(lps *lps_ptr, const uint8_t *data, size_t len) {
//     lps_ptr->level = 1;
//     lps_ptr->size = 0;
//     lps_ptr->capacity = (len / CONSTANT_FACTOR);
//     lps_ptr->cores =
//         (struct core *)malloc((len / CONSTANT_FACTOR) * sizeof(struct core));
//     lps_ptr->size = parse1_bam_syncmers_old(data, len, lps_ptr->cores);
// }

// void init_lps_kseq(lps *lps_ptr, const char *data, size_t len) {
//     lps_ptr->level = 1;
//     lps_ptr->size = 0;
//     lps_ptr->cores =
//         (struct core *)malloc((len / CONSTANT_FACTOR) * sizeof(struct core));
//     lps_ptr->size = parse1_kseq_syncmers(data, len, lps_ptr->cores);
// }
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
// // uint8_t *onehotto2bitDNA_bivec(const uint8_t *data, size_t len) {
// //     uint8_t *ret = (uint8_t *)malloc(len);
// // #pragma clang loop vectorize(enable) interleave(enable)
// //     for (size_t i = 0; i < len / 2; ++i) {
// //         ret[2 * i] = 27 - std::countl_zero(static_cast<unsigned>(data[i]));
// //         ret[2 * i + 1] = std::countr_zero(static_cast<unsigned>(data[i]));
// //     }
// //     if (len % 2 == 1) {
// //         ret[len - 1] = 27 - std::countl_zero(static_cast<unsigned>(data[len / 2]));
// //     }
// //     return ret;
// }

// lcp_read parse_cores_from_kseq(const char *s, size_t len, int lcp_level) {
//     // uint64_t estimated_core_size = (int)(len / std::pow(1.5, lcp_level));

//     lps str;
//     init_lps_kseq(&str, s, len);
//     lps_deepen(&str, lcp_level);
//     lcp_read lr;
//     lr.reserve(str.size);
//     for (int i = 0; i < str.size; ++i) {
//         lr.cores.emplace_back(str.cores[i].label);
//         lr.pos.emplace_back(str.cores[i].start, str.cores[i].end);
//     }
//     free_lps(&str);
//     return lr;
// }

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
    // SUBCASE("Even") {
    //     uint8_t *conversion =
    //         onehotto2bitDNA_bivec(bytes.data(), 2 * bytes.size());
    //     for (unsigned i = 0; i < bytes.size(); ++i) {
    //         CHECK(expected[bytes[i]].first == conversion[2 * i]);
    //         CHECK(expected[bytes[i]].second == conversion[2 * i + 1]);
    //     }
    //     free(conversion);
    // }
    SUBCASE("Even reverse") {
        charsi c;
        c.require(2 * bytes.size());
        onehotto2bitDNA_bivec_chs_reverse(bytes.data(), 2 * bytes.size(), c);
        size_t bs = bytes.size() - 1;
        for (unsigned i = 0; i < bytes.size(); ++i) {
            CHECK(3 - expected[bytes[bs - i]].first == c.data[2 * i]);
            CHECK(3 - expected[bytes[bs - i]].second == c.data[2 * i + 1]);
        }
    }
    bytes.push_back((1 << 4));
    expected[bytes.back()] = {tab[1], 0};
    // SUBCASE("Odd") {
    //     uint8_t *conversion =
    //         onehotto2bitDNA_bivec(bytes.data(), 2 * bytes.size() - 1);
    //     for (unsigned i = 0; i < bytes.size(); ++i) {
    //         CHECK(expected[bytes[i]].first == conversion[2 * i]);
    //         CHECK(expected[bytes[i]].second == conversion[2 * i + 1]);
    //     }
    //     free(conversion);
    // }
}

// TEST_CASE("CrumbledSubSeq reverse complement") {
//     auto expected_rc = [](std::string_view seq) {
//         auto complement = [](char base) {
//             switch (base) {
//                 case 'A':
//                     return 'T';
//                 case 'C':
//                     return 'G';
//                 case 'G':
//                     return 'C';
//                 case 'T':
//                     return 'A';
//                 default:
//                     return 'N';
//             }
//         };
//         std::string rc;
//         rc.reserve(seq.size());
//         for (auto it = seq.rbegin(); it != seq.rend(); ++it) {
//             rc.push_back(complement(*it));
//         }
//         return rc;
//     };

//     SUBCASE("ASCII constructor maintains order") {
//         const std::string sequence =
//             "CTCGCGACCTGGCTGCCCAGATGGGCCTGTTCCCCAGCCCCGAGCTCCTAGGCTATGACAGCTTAC"
//             "CAAGGTTCCGGCCACTTGTCAGATAGCAGAGTCGGCCACGGGGCTCAAGGCTGACTTGTACTTCCC"
//             "TCTGGCTGCTCGTAGGACGCCACCCTGAGAATGGACCGCAGTACCAGCCTGCTGT";
//         CrumbledSubSeq64 crumb{sequence.c_str(), 0, sequence.size(),
//                                CrumbledSubSeq64::from_ascii{}};
//         std::string forward(sequence.size(), '\0');
//         crumb.convert_to_ascii(forward.data());
//         CHECK_EQ(sequence, forward);

//         auto rc = crumb.reverse_complement();
//         std::string rc_ascii(sequence.size(), '\0');
//         char *starting_from = rc.convert_to_ascii(rc_ascii.data());
//         CHECK_EQ(expected_rc(sequence), std::string{starting_from});
//     }

//     SUBCASE("Nibble constructor respects skip offset") {
//         const std::string raw =
//             "CGGGGTGAGCACCTCAGGAAAAGCAGGCTGGGGGTCGGGGGCAAGCTCCGCTTTGGCCGGGCTGAG"
//             "GTCAAGGCACCTGTGAGTGAGTTGGAGTGAGGGGAGGAGCGGCAGGGCCCTGGAGTAGAGTCCAGG"
//             "GGAGGTCCCAAGCCCAGGGATAAGCTCTGTCATGGGGTGAGAGACAGGGCAGAGACAGACTGGAGG"
//             "CTGGGAGGCCAGGAGAAGGGCTGGAGTCAGGGACAGGGGTCTATACCCACCCTCCCCCCGCCACCC"
//             "CACTGGTGTCTGAGAAAGGACTCCCACAGGGTCCCCAGAACTTGGTCTGAGTTCCTCAGTGTGCAC"
//             "TGGTCACTACTAGCATGGACTGGAAAGCCGGGCGTGGAG";
//         auto encode_bam = [](std::string_view seq) {
//             gtl::vector<uint8_t> packed((seq.size() + 1) / 2, 0);
//             for (size_t i = 0; i < seq.size(); ++i) {
//                 uint8_t nibble = seq_nt16_table[(uint8_t)seq[i]];
//                 if ((i & 1U) == 0) {
//                     packed[i / 2] = nibble << 4;
//                 } else {
//                     packed[i / 2] |= nibble & 0x0F;
//                 }
//             }
//             return packed;
//         };
//         auto bam_bytes = encode_bam(raw);
//         for (size_t i = 0; i < raw.size() - 8; ++i) {
//             size_t from = 3;  // force skip_first_n != 0
//             size_t to = 3 + 5 + i;
//             CrumbledSubSeq64 crumb{bam_bytes.data(), from, to,
//                                    CrumbledSubSeq64::from_nibbles{}};
//             std::string view(to - from, '\0');
//             char *starting_from = crumb.convert_to_ascii(view.data());
//             CHECK_EQ(raw.substr(from, to - from),
//             std::string{starting_from});

//             auto rc = crumb.reverse_complement();
//             std::string rc_ascii(view.size(), '\0');
//             rc.convert_to_ascii(rc_ascii.data());
//             CHECK_EQ(expected_rc(raw.substr(from, to - from)), rc_ascii);

//             CrumbledSubSeq64 crumb64{bam_bytes.data(), from, to,
//                                      CrumbledSubSeq64::from_nibbles{}};
//             auto rc64_old = crumb64.reverse_complement_old();
//             std::string rc_old_ascii(view.size(), '\0');
//             char *starting_from_old =
//                 rc64_old.convert_to_ascii(rc_old_ascii.data());
//             CHECK_EQ(expected_rc(raw.substr(from, to - from)),
//                      starting_from_old);
//             INFO(crumb64.skip_first_n);

//             auto rc64_fast = crumb64.reverse_complement();
//             std::string rc_fast_ascii(view.size(), '\0');
//             starting_from = rc64_fast.convert_to_ascii(rc_fast_ascii.data());
//             INFO(rc64_fast.skip_first_n);
//             INFO((to - from));
//             CHECK_EQ(expected_rc(raw.substr(from, to - from)),
//                      std::string{starting_from});
//             CHECK_EQ(std::string{starting_from_old},
//                      std::string{starting_from});
//         }
//     }
// }

TEST_CASE("K-mer view on Crumbs") {
    gtl::vector<uint8_t> bytes{{0, 1, 2, 3, 0, 1, 2, 3}};
    SUBCASE("K==4") {
        static constexpr int K = 4;
        KmerViewOnInts<uint64_t, 2> kview(K);
        size_t i = kview.fill(bytes.data());
        uint32_t expected =
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
        uint32_t expected = bytes[0] << 8 | bytes[1] << 6 | bytes[2] << 4 |
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
    //     SUBCASE("Simd Nibble2base") {
    //         char *seq_simd = (char *)malloc(rec->core.l_qseq + 1);
    //         nibble2base(bam_get_seq(rec), seq_simd, rec->core.l_qseq);

    //         CHECK_EQ(std::string{seq_gt}, std::string{seq_simd});
    //         free(seq_simd);
    //     }

    //     SUBCASE("Translating intervals from nibbles to bytes") {
    //         SUBCASE("EVEN") {
    //             int from = 168;
    //             int to = 266;
    //             char *seq = (char *)malloc(to - from + 1);
    //             nibble2base(bam_get_seq(rec) + from / 2, seq, to - from);

    //             CHECK_EQ(std::string{seq_gt}.substr(from, to - from),
    //                      std::string{seq});
    //             free(seq);
    //         }
    //         SUBCASE("ODD_EVEN") {
    //             int from = 167;
    //             int to = 266;
    //             char *seq = (char *)malloc(to - from + 2);
    //             nibble2base(bam_get_seq(rec) + from / 2, seq, to - from + 1);

    //             CHECK_EQ(std::string{seq_gt}.substr(from, to - from),
    //                      std::string{seq}.substr(1));
    //             free(seq);
    //         }
    //         SUBCASE("ODD_ODD") {
    //             int from = 167;
    //             int to = 267;
    //             char *seq = (char *)malloc(to - from + 2);
    //             nibble2base(bam_get_seq(rec) + from / 2, seq, to - from + 1);

    //             CHECK_EQ(std::string{seq_gt}.substr(from, to - from),
    //                      std::string{seq}.substr(1));
    //             free(seq);
    //         }
    //         SUBCASE("EVEN_ODD") {
    //             int from = 168;
    //             int to = 267;
    //             char *seq = (char *)malloc(to - from + 2);
    //             nibble2base(bam_get_seq(rec) + from / 2, seq, to - from);

    //             CHECK_EQ(std::string{seq_gt}.substr(from, to - from),
    //                      std::string{seq});
    //             free(seq);
    //         }
    //         SUBCASE("ANY") {
    //             int seqlen = strlen(seq_gt);
    //             srand(time(NULL));
    //             for (int casei = 0; casei < 50; ++casei) {
    //                 int to = rand() % (seqlen - 100) + 100;
    //                 int from = rand() % (to - 100);
    //                 char *seq = (char *)malloc(to - from + 2);

    //                 if (from % 2 == 0) {
    //                     nibble2base(bam_get_seq(rec) + from / 2, seq, to -
    //                     from); seq[to - from] = '\0';
    //                     CHECK_EQ(std::string{seq_gt}.substr(from, to - from),
    //                              std::string{seq});
    //                 } else {
    //                     nibble2base(bam_get_seq(rec) + from / 2, seq,
    //                                 to - from + 1);
    //                     seq[to - from + 1] = '\0';
    //                     CHECK_EQ(std::string{seq_gt}.substr(from, to - from),
    //                              std::string{seq}.substr(1));
    //                 }
    //                 free(seq);
    //             }
    //         }
    //     }
    //     SUBCASE("Copying Nibble-encoded data and decoding it to bases") {
    //         int seqlen = strlen(seq_gt);
    //         srand(time(NULL));
    //         int N_RANDOM_CASES = 50;
    //         gtl::vector<bam_nibble_data> copied;

    //         gtl::vector<size_t> froms;
    //         gtl::vector<size_t> tos;
    //         for (int casei = 0; casei < N_RANDOM_CASES; ++casei) {
    //             int to = rand() % (seqlen - 100) + 100;
    //             int from = rand() % (to - 100);
    //             copied.emplace_back(bam_get_seq(rec), from, to);
    //             froms.push_back(from);
    //             tos.push_back(to);
    //         }
    //         char buffer[1024 * 1024];
    //         for (int casei = 0; casei < N_RANDOM_CASES; ++casei) {
    //             auto &seq = copied[casei];
    //             size_t from = froms[casei];
    //             size_t to = tos[casei];
    //             seq.convert_to_ascii(buffer);
    //             INFO("from: ", from, " to: ", to, "seq.len: ", seq.len);
    //             CHECK_EQ(std::string{seq_gt}.substr(from, to - from),
    //                      std::string{buffer}.substr(seq.skip_first ? 1 : 0));
    //         }
    //     }
    //     bam_hdr_destroy(hdr);
    //     sam_close(fp);
    //     free(seq_gt);
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
            char *starting_from = seq.convert_to_ascii(buffer);
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
        // printf(">R%llu\n", node.first);
        // charsi chs;
        // for(const auto &seq : rev_seqs) {
        //     chs.require(seq.len+1);
        //     seq.convert_to_ascii((char* ) chs.data);
        //     printf("%s", chs.data + seq.skip_first_n);
        // }
        // printf("\n");
        lps_deepen(&str2, 2);
        for (int i = 0; i < str2.size; ++i) {
            lr2.cores.emplace_back(str2.cores[i].label);
            lr2.pos.emplace_back(str2.cores[i].start, str2.cores[i].end);
        }
        CHECK(lr.cores.size() == lr2.cores.size());
    }

    bam_hdr_destroy(hdr);
    sam_close(fp);
    free(seq_gt);
}

// Reference scalar-only two-pass: decode then RLE.
// Used as ground truth for testing onehot_to_rle.
static int64_t reference_onehot_to_rle(const uint8_t *data, int64_t len,
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

    auto run_test = [&](const gtl::vector<uint8_t> &bases, const char *label) {
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
