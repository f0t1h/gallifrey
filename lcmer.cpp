#include "lcmer.h"

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
    cr->bit_rep = 0x8000000000000000 | cr->label;
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
        int clz = __builtin_clz(data[i]);
        int ctz = __builtin_ctz(data[i]);

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
        ret[2 * i] = 27 - __builtin_clz(data[i]);
        ret[2 * i + 1] = __builtin_ctz(data[i]);
    }
    return ret;
}
void onehotto2bitDNA_bivec_chs(const uint8_t *data, size_t len, charsi &chs) {
    chs.require(len);

    for (size_t i = 0; i < len / 2; ++i) {
        chs.data[2 * i] = 27 - __builtin_clz(data[i]);
        chs.data[2 * i + 1] = __builtin_ctz(data[i]);
    }
    if (len % 2 == 1) {
        chs.data[len - 1] = 27 - __builtin_clz(data[len / 2]);
    }
}
void onehotto2bitDNA_bivec_chs_reverse(const uint8_t *data, size_t len,
                                       charsi &chs) {
    chs.require(len);
    size_t mid = len / 2 - 1;
    for (size_t i = 0; i < len / 2; ++i) {
        chs.data[2 * (mid - i)] = __builtin_clz(data[i]) - 24;
        chs.data[2 * (mid - i) + 1] = 3 - __builtin_ctz(data[i]);
    }
    if (len % 2 == 1) {
        chs.data[0] = __builtin_clz(data[len / 2]) - 24;
    }
}
uint8_t *onehotto2bitDNA_bivec(const uint8_t *data, size_t len) {
    uint8_t *ret = (uint8_t *)malloc(len);
#pragma clang loop vectorize(enable) interleave(enable)
    for (size_t i = 0; i < len / 2; ++i) {
        ret[2 * i] = 27 - __builtin_clz(data[i]);
        ret[2 * i + 1] = __builtin_ctz(data[i]);
    }
    if (len % 2 == 1) {
        ret[len - 1] = 27 - __builtin_clz(data[len / 2]);
    }
    return ret;
}

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
    SUBCASE("Even") {
        uint8_t *conversion =
            onehotto2bitDNA_bivec(bytes.data(), 2 * bytes.size());
        for (unsigned i = 0; i < bytes.size(); ++i) {
            CHECK(expected[bytes[i]].first == conversion[2 * i]);
            CHECK(expected[bytes[i]].second == conversion[2 * i + 1]);
        }
        free(conversion);
    }
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
    SUBCASE("Odd") {
        uint8_t *conversion =
            onehotto2bitDNA_bivec(bytes.data(), 2 * bytes.size() - 1);
        for (unsigned i = 0; i < bytes.size(); ++i) {
            CHECK(expected[bytes[i]].first == conversion[2 * i]);
            CHECK(expected[bytes[i]].second == conversion[2 * i + 1]);
        }
        free(conversion);
    }
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
            int from = to - 100 - (rand() % 12);
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
        lps str1 = {0,0,0,0};
        lps str2 = {0,0,0,0};
        using SYNCMER_HASHER = syncmer::MurMur;
        using SyncmerType =
            syncmer::GenericSyncmerRule<uint32_t, SYNCMER_HASHER>;
        decltype(asm_config::syncmer) spara;
        spara.k = 8;
        spara.density = 0.05;
        SyncmerType fnctr(spara, 11);
        charsi chs;
        lcp_read lr;
        parse_cores_from_bam_record(bam_get_seq(rec), rec->core.l_qseq, 3,
                                    &str1, fnctr, chs, lr);
        gtl::vector<CrumbledSubSeq64> ss;
        ss.emplace_back(bam_get_seq(rec), 0, rec->core.l_qseq, 0,
                        CrumbledSubSeq64::from_nibbles{});

        lcp_read lr2;
        syncmer_lcp::init_lps_crumbled(&str2, ss, fnctr);
        // printf(">R%llu\n", node.first);
        // charsi chs;
        // for(const auto &seq : rev_seqs) {
        //     chs.require(seq.len+1);
        //     seq.convert_to_ascii((char* ) chs.data);
        //     printf("%s", chs.data + seq.skip_first_n);
        // }
        // printf("\n");
        lps_deepen(&str2, 3);
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

// ---------------------------------------------------------------------------
// Helpers for base_lcp parse1 tests
// ---------------------------------------------------------------------------

/// Pack ASCII bases (A/C/G/T) into BAM 4-bit nibble format (2 bases per byte).
/// `out` must be at least (len+1)/2 bytes.
static void ascii_to_bam_nibble(const char *ascii, size_t len, uint8_t *out) {
    static const uint8_t lut[256] = {
        ['A'] = 1, ['a'] = 1, ['C'] = 2, ['c'] = 2,
        ['G'] = 4, ['g'] = 4, ['T'] = 8, ['t'] = 8,
    };
    memset(out, 0, (len + 1) / 2);
    for (size_t i = 0; i < len; ++i) {
        out[i >> 1] |= lut[(unsigned char)ascii[i]] << (((~i) & 1) << 2);
    }
}

/// Expected 2-bit code for a base character (A=0 C=1 G=2 T=3).
static int base_code(char c) {
    switch (c) {
        case 'A': case 'a': return 0;
        case 'C': case 'c': return 1;
        case 'G': case 'g': return 2;
        case 'T': case 't': return 3;
        default: return -1;
    }
}

/// Run base_lcp::parse1_frombam on an ASCII sequence.
static uint64_t run_parse1_bam(const char *ascii, size_t len,
                               struct core *cores, uint64_t offset = 0) {
    uint8_t *nibble = new uint8_t[(len + 1) / 2]();
    ascii_to_bam_nibble(ascii, len, nibble);
    uint64_t n = base_lcp::parse1_frombam(nibble, (int64_t)len, cores, offset);
    delete[] nibble;
    return n;
}

/// Run base_lcp::parse1_from_SCT_vec with a SINGLE CrumbledSubSeq64 from ASCII.
static uint64_t run_parse1_sct_single(const char *ascii, size_t len,
                                      struct core *cores,
                                      uint64_t offset = 0) {
    gtl::vector<CrumbledSubSeq64> sseqs;
    sseqs.emplace_back(ascii, (size_t)0, len, /*origin_node=*/0,
                       CrumbledSubSeq64::from_ascii{});
    return base_lcp::parse1_from_SCT_vec(sseqs, cores, offset);
}

/// Run base_lcp::parse1_from_SCT_vec with the ASCII string split into
/// `nfrags` consecutive fragments.
static uint64_t run_parse1_sct_multi(const char *ascii, size_t len,
                                     size_t nfrags, struct core *cores,
                                     uint64_t offset = 0) {
    gtl::vector<CrumbledSubSeq64> sseqs;
    size_t chunk = len / nfrags;
    for (size_t i = 0; i < nfrags; ++i) {
        size_t from = i * chunk;
        size_t to = (i + 1 == nfrags) ? len : from + chunk;
        sseqs.emplace_back(ascii, from, to, /*origin_node=*/0,
                           CrumbledSubSeq64::from_ascii{});
    }
    return base_lcp::parse1_from_SCT_vec(sseqs, cores, offset);
}

/// Compare two core arrays element by element using doctest.
static void assert_cores_equal(const struct core *a, uint64_t na,
                               const struct core *b, uint64_t nb,
                               const char *label = "") {
    INFO(label);
    REQUIRE(na == nb);
    for (uint64_t i = 0; i < na; ++i) {
        INFO("core index ", i);
        CHECK(a[i].start    == b[i].start);
        CHECK(a[i].end      == b[i].end);
        CHECK(a[i].label    == b[i].label);
        CHECK(a[i].bit_rep  == b[i].bit_rep);
        CHECK(a[i].bit_size == b[i].bit_size);
    }
}

/// Compute the label that init_core1hot would produce for a 3-base core
/// (distance=3) at ASCII positions [p, p+3).
static uint64_t expected_label3(const char *ascii, size_t p) {
    return ((uint64_t)(3 - 2) << 6)
         | ((uint64_t)base_code(ascii[p])     << 4)
         | ((uint64_t)base_code(ascii[p + 1]) << 2)
         | ((uint64_t)base_code(ascii[p + 2]));
}

// ---------------------------------------------------------------------------
// parse1_frombam — individual core-type tests
// ---------------------------------------------------------------------------

TEST_CASE("parse1_frombam: LMIN cores") {
    // "CAG": C(1)>A(0)<G(2) → single LMIN
    const char *seq = "CAG";
    struct core cores[16];
    uint64_t n = run_parse1_bam(seq, strlen(seq), cores);

    REQUIRE(n == 1);
    CHECK(cores[0].start == 0);
    CHECK(cores[0].end   == 3);
    CHECK(cores[0].label == expected_label3(seq, 0));
}

TEST_CASE("parse1_frombam: RINT cores") {
    // "ACTTTTGA": RINT at [1,7) with middle_count=4 → distance 6
    const char *seq = "ACTTTTGA";
    struct core cores[16];
    uint64_t n = run_parse1_bam(seq, strlen(seq), cores);

    bool found = false;
    for (uint64_t i = 0; i < n; ++i) {
        uint64_t dist = ((cores[i].label >> 6) & 0x3F) + 2;
        if (dist == 6 && cores[i].start == 1 && cores[i].end == 7) {
            found = true;
            CHECK(((cores[i].label >> 4) & 3) == 1);  // first = C
            CHECK(((cores[i].label >> 2) & 3) == 3);  // penult = T
            CHECK(((cores[i].label)      & 3) == 2);  // last = G
        }
    }
    CHECK(found);
}

TEST_CASE("parse1_frombam: LMAX cores") {
    // "CAGTCAG" → LMIN(0,3), LMAX(2,5), LMIN(4,7)
    const char *seq = "CAGTCAG";
    struct core cores[16];
    uint64_t n = run_parse1_bam(seq, strlen(seq), cores);

    REQUIRE(n == 3);
    CHECK(cores[0].start == 0);  CHECK(cores[0].end == 3);
    CHECK(cores[0].label == expected_label3(seq, 0));
    CHECK(cores[1].start == 2);  CHECK(cores[1].end == 5);
    CHECK(cores[1].label == expected_label3(seq, 2));
    CHECK(cores[2].start == 4);  CHECK(cores[2].end == 7);
    CHECK(cores[2].label == expected_label3(seq, 4));
}

TEST_CASE("parse1_frombam: SSEQ cores") {
    // "GACGTACGA" → LMIN(0,3), SSEQ(2,5), LMIN(4,7)
    const char *seq = "GACGTACGA";
    struct core cores[32];
    uint64_t n = run_parse1_bam(seq, strlen(seq), cores);

    REQUIRE(n == 3);
    CHECK(cores[0].start == 0);  CHECK(cores[0].end == 3);
    CHECK(cores[1].start == 2);  CHECK(cores[1].end == 5);
    CHECK(cores[2].start == 4);  CHECK(cores[2].end == 7);
}

TEST_CASE("parse1_frombam: minimum length") {
    struct core cores[8];
    CHECK(run_parse1_bam("CAG", 3, cores) == 1);
    CHECK(run_parse1_bam("CA",  2, cores) == 0);
    CHECK(run_parse1_bam("A",   1, cores) == 0);
}

// ---------------------------------------------------------------------------
// parse1_from_SCT_vec — mirror of the above
// ---------------------------------------------------------------------------

TEST_CASE("parse1_from_SCT_vec: LMIN cores") {
    const char *seq = "CAG";
    struct core cores[16];
    uint64_t n = run_parse1_sct_single(seq, strlen(seq), cores);

    REQUIRE(n == 1);
    CHECK(cores[0].start == 0);
    CHECK(cores[0].end   == 3);
    CHECK(cores[0].label == expected_label3(seq, 0));
}

TEST_CASE("parse1_from_SCT_vec: RINT cores") {
    const char *seq = "ACTTTTGA";
    struct core cores[16];
    uint64_t n = run_parse1_sct_single(seq, strlen(seq), cores);

    bool found = false;
    for (uint64_t i = 0; i < n; ++i) {
        uint64_t dist = ((cores[i].label >> 6) & 0x3F) + 2;
        if (dist == 6 && cores[i].start == 1 && cores[i].end == 7) {
            found = true;
            CHECK(((cores[i].label >> 4) & 3) == 1);
            CHECK(((cores[i].label >> 2) & 3) == 3);
            CHECK(((cores[i].label)      & 3) == 2);
        }
    }
    CHECK(found);
}

TEST_CASE("parse1_from_SCT_vec: LMAX cores") {
    const char *seq = "CAGTCAG";
    struct core cores[16];
    uint64_t n = run_parse1_sct_single(seq, strlen(seq), cores);

    REQUIRE(n == 3);
    CHECK(cores[0].start == 0);  CHECK(cores[0].end == 3);
    CHECK(cores[1].start == 2);  CHECK(cores[1].end == 5);
    CHECK(cores[2].start == 4);  CHECK(cores[2].end == 7);
}

TEST_CASE("parse1_from_SCT_vec: SSEQ cores") {
    const char *seq = "GACGTACGA";
    struct core cores[32];
    uint64_t n = run_parse1_sct_single(seq, strlen(seq), cores);

    REQUIRE(n == 3);
    CHECK(cores[0].start == 0);  CHECK(cores[0].end == 3);
    CHECK(cores[1].start == 2);  CHECK(cores[1].end == 5);
    CHECK(cores[2].start == 4);  CHECK(cores[2].end == 7);
}

TEST_CASE("parse1_from_SCT_vec: minimum / empty") {
    struct core cores[8];
    CHECK(run_parse1_sct_single("CAG", 3, cores) == 1);
    CHECK(run_parse1_sct_single("CA",  2, cores) == 0);
    CHECK(run_parse1_sct_single("A",   1, cores) == 0);

    gtl::vector<CrumbledSubSeq64> empty;
    CHECK(base_lcp::parse1_from_SCT_vec(empty, cores, 0) == 0);
}

// ---------------------------------------------------------------------------
// Equivalence: prove parse1_frombam == parse1_from_SCT_vec
// ---------------------------------------------------------------------------

TEST_CASE("parse1 equivalence: BAM vs single-fragment SCT") {
    const char *sequences[] = {
        "CAG",
        "ACTTTTGA",
        "CAGTCAG",
        "GACGTACGA",
        "ACGTACGTACGTACGT",
        "TTTACGTTTTTGCAGTCAGTCATTTG",
        "GCATTTTTTTTTTTTTTTTCGATCGATCGA",
        "ACGACGACGACG",
        "TGCATGCATGCATGCATGCA",
        "AACCGGTTAACCGGTT",
        "GATTACAGATTAGATTACAGATTAGATTACAGATTA",
    };

    for (const char *seq : sequences) {
        size_t len = strlen(seq);
        struct core bam_cores[256];
        struct core sct_cores[256];
        uint64_t nb = run_parse1_bam(seq, len, bam_cores);
        uint64_t ns = run_parse1_sct_single(seq, len, sct_cores);
        assert_cores_equal(bam_cores, nb, sct_cores, ns, seq);
    }
}

TEST_CASE("parse1 equivalence: multi-fragment SCT matches single-fragment") {
    const char *sequences[] = {
        "CAGTCAG",
        "ACTTTTGA",
        "GACGTACGA",
        "ACGTACGTACGTACGT",
        "TTTACGTTTTTGCAGTCAGTCATTTG",
        "GATTACAGATTAGATTACAGATTAGATTACAGATTA",
    };

    for (const char *seq : sequences) {
        size_t len = strlen(seq);
        struct core single_cores[256];
        uint64_t n_single = run_parse1_sct_single(seq, len, single_cores);

        SUBCASE("2 fragments") {
            struct core multi_cores[256];
            uint64_t n_multi = run_parse1_sct_multi(seq, len, 2, multi_cores);
            assert_cores_equal(single_cores, n_single,
                               multi_cores, n_multi, seq);
        }
        SUBCASE("3 fragments") {
            struct core multi_cores[256];
            uint64_t n_multi = run_parse1_sct_multi(seq, len, 3, multi_cores);
            assert_cores_equal(single_cores, n_single,
                               multi_cores, n_multi, seq);
        }
        if (len >= 5) {
            SUBCASE("5 fragments") {
                struct core multi_cores[256];
                uint64_t n_multi =
                    run_parse1_sct_multi(seq, len, 5, multi_cores);
                assert_cores_equal(single_cores, n_single,
                                   multi_cores, n_multi, seq);
            }
        }
    }
}

TEST_CASE("parse1 equivalence: offset propagation") {
    const char *seq = "CAGTCAG";
    size_t len = strlen(seq);
    uint64_t offset = 1000;

    struct core bam_cores[16], sct_cores[16];
    uint64_t nb = run_parse1_bam(seq, len, bam_cores, offset);
    uint64_t ns = run_parse1_sct_single(seq, len, sct_cores, offset);
    assert_cores_equal(bam_cores, nb, sct_cores, ns, "offset propagation");

    for (uint64_t i = 0; i < nb; ++i) {
        CHECK(bam_cores[i].start >= offset);
        CHECK(bam_cores[i].end   >= offset);
    }
}

TEST_CASE("parse1 equivalence: long RINT run") {
    // "C" + 50x"T" + "G"  →  single RINT spanning [0,52)
    std::string seq(52, 'T');
    seq[0]  = 'C';
    seq[51] = 'G';

    struct core bam_cores[64], sct_cores[64];
    uint64_t nb = run_parse1_bam(seq.c_str(), seq.size(), bam_cores);
    uint64_t ns = run_parse1_sct_single(seq.c_str(), seq.size(), sct_cores);
    assert_cores_equal(bam_cores, nb, sct_cores, ns, "long RINT");

    REQUIRE(nb == 1);
    CHECK(bam_cores[0].start == 0);
    CHECK(bam_cores[0].end   == 52);
    uint64_t dist = ((bam_cores[0].label >> 6) & 0x3F) + 2;
    CHECK(dist == 52);
    CHECK(((bam_cores[0].label >> 4) & 3) == 1);  // C
    CHECK(((bam_cores[0].label >> 2) & 3) == 3);  // T
    CHECK(((bam_cores[0].label)      & 3) == 2);  // G
}

TEST_CASE("parse1 equivalence: long RINT split across many fragments") {
    std::string seq(52, 'T');
    seq[0]  = 'C';
    seq[51] = 'G';

    struct core single_cores[64], multi_cores[64];
    uint64_t n_single =
        run_parse1_sct_single(seq.c_str(), seq.size(), single_cores);
    uint64_t n_multi =
        run_parse1_sct_multi(seq.c_str(), seq.size(), 26, multi_cores);
    assert_cores_equal(single_cores, n_single, multi_cores, n_multi,
                       "long RINT 26 frags");
}

TEST_CASE("init_lps + lps_deepen equivalence: BAM vs crumbled path") {
    const char *sequences[] = {
        "CAGTCAGTCAGTCAGTCAGT",
        "ACGTACGTACGTACGTACGT",
        "GATTACAGATTAGATTACAGATTAGATTACAGATTA",
        "TTTACGTTTTTGCAGTCAGTCATTTGCATTTTTTT",
    };

    for (const char *seq : sequences) {
        size_t len = strlen(seq);
        INFO("sequence: ", seq);

        // --- BAM path ---
        uint8_t *nibble = new uint8_t[(len + 1) / 2]();
        ascii_to_bam_nibble(seq, len, nibble);
        lps str_bam = {0, 0, 0, nullptr};
        base_lcp::init_lps_bam(&str_bam, nibble, (int)len);
        lps_deepen(&str_bam, 2);

        // --- Crumbled path ---
        gtl::vector<CrumbledSubSeq64> sseqs;
        sseqs.emplace_back(seq, (size_t)0, len, 0,
                           CrumbledSubSeq64::from_ascii{});
        lps str_cr = {0, 0, 0, nullptr};
        // syncmer_functor is unused by base_lcp, but the API requires it
        using SyncmerType =
            syncmer::GenericSyncmerRule<uint32_t, syncmer::MurMur>;
        syncmer_params spara;
        spara.k = 8;
        spara.density = 0.9;
        SyncmerType fnctr(spara, 11);
        base_lcp::init_lps_crumbled(&str_cr, sseqs, fnctr);
        lps_deepen(&str_cr, 2);

        REQUIRE(str_bam.size == str_cr.size);
        for (int i = 0; i < str_bam.size; ++i) {
            INFO("deepened core ", i);
            CHECK(str_bam.cores[i].start   == str_cr.cores[i].start);
            CHECK(str_bam.cores[i].end     == str_cr.cores[i].end);
            CHECK(str_bam.cores[i].label   == str_cr.cores[i].label);
            CHECK(str_bam.cores[i].bit_rep == str_cr.cores[i].bit_rep);
        }

        free_lps(&str_bam);
        free_lps(&str_cr);
        delete[] nibble;
    }
}

TEST_CASE("parse1 equivalence: real BAM read") {
    std::string input = "../tests/test_single.bam";
    samFile *fp = sam_open(input.c_str(), "r");
    if (fp == NULL) {
        MESSAGE("Skipping — could not open ", input);
        return;
    }
    bam_hdr_t *hdr = sam_hdr_read(fp);
    bam1_t *rec = bam_init1();
    if (sam_read1(fp, hdr, rec) <= 0) {
        bam_hdr_destroy(hdr);
        sam_close(fp);
        FAIL("Failed to read BAM record");
    }

    int seqlen = rec->core.l_qseq;
    const uint8_t *bam_data = bam_get_seq(rec);

    // ASCII for CrumbledSubSeq64
    char *ascii = (char *)malloc(seqlen + 1);
    for (int i = 0; i < seqlen; ++i)
        ascii[i] = seq_nt16_str[bam_seqi(bam_data, i)];
    ascii[seqlen] = '\0';

    size_t cap = (size_t)(seqlen / 1.2);
    struct core *bam_cores = (struct core *)malloc(cap * sizeof(struct core));
    uint64_t nb = base_lcp::parse1_frombam(bam_data, seqlen, bam_cores, 0);

    gtl::vector<CrumbledSubSeq64> sseqs;
    sseqs.emplace_back(bam_data, (size_t)0, (size_t)seqlen, 0,
                       CrumbledSubSeq64::from_nibbles{});
    struct core *sct_cores = (struct core *)malloc(cap * sizeof(struct core));
    uint64_t ns = base_lcp::parse1_from_SCT_vec(sseqs, sct_cores, 0);

    INFO("read length=", seqlen, "  bam_cores=", nb, "  sct_cores=", ns);
    REQUIRE(nb == ns);
    for (uint64_t i = 0; i < nb; ++i) {
        INFO("core ", i);
        CHECK(bam_cores[i].start    == sct_cores[i].start);
        CHECK(bam_cores[i].end      == sct_cores[i].end);
        CHECK(bam_cores[i].label    == sct_cores[i].label);
        CHECK(bam_cores[i].bit_rep  == sct_cores[i].bit_rep);
        CHECK(bam_cores[i].bit_size == sct_cores[i].bit_size);
    }

    SUBCASE("10-fragment split of real read") {
        gtl::vector<CrumbledSubSeq64> frags;
        size_t chunk = seqlen / 10;
        for (int i = 0; i < 10; ++i) {
            size_t from = i * chunk;
            size_t to = (i == 9) ? (size_t)seqlen : from + chunk;
            frags.emplace_back(ascii, from, to, 0,
                               CrumbledSubSeq64::from_ascii{});
        }
        struct core *frag_cores =
            (struct core *)malloc(cap * sizeof(struct core));
        uint64_t nf = base_lcp::parse1_from_SCT_vec(frags, frag_cores, 0);

        INFO("10 frags: ", nf, " vs single: ", ns);
        REQUIRE(nf == ns);
        for (uint64_t i = 0; i < nf; ++i) {
            INFO("core ", i);
            CHECK(frag_cores[i].start   == sct_cores[i].start);
            CHECK(frag_cores[i].end     == sct_cores[i].end);
            CHECK(frag_cores[i].label   == sct_cores[i].label);
            CHECK(frag_cores[i].bit_rep == sct_cores[i].bit_rep);
        }
        free(frag_cores);
    }

    free(bam_cores);
    free(sct_cores);
    free(ascii);
    bam_destroy1(rec);
    bam_hdr_destroy(hdr);
    sam_close(fp);
}

// ---------------------------------------------------------------------------
// base_lcp::parse_cores_from_bam_record (pure_lcp path) — forward
// ---------------------------------------------------------------------------

TEST_CASE("base_lcp::parse_cores_from_bam_record forward matches parse1_frombam") {
    std::string input = "../tests/test_single.bam";
    samFile *fp = sam_open(input.c_str(), "r");
    if (fp == NULL) {
        MESSAGE("Skipping — could not open ", input);
        return;
    }
    bam_hdr_t *hdr = sam_hdr_read(fp);
    bam1_t *rec = bam_init1();
    if (sam_read1(fp, hdr, rec) <= 0) {
        bam_hdr_destroy(hdr);
        sam_close(fp);
        FAIL("Failed to read BAM record");
    }

    int seqlen = rec->core.l_qseq;
    const uint8_t *bam_data = bam_get_seq(rec);

    // --- Direct parse1_frombam on original nibble data ---
    size_t cap = (size_t)(seqlen / 1.2);
    struct core *direct_cores =
        (struct core *)malloc(cap * sizeof(struct core));
    uint64_t nd = base_lcp::parse1_frombam(bam_data, seqlen, direct_cores, 0);

    // --- Through parse_cores_from_bam_record (forward, pure_lcp path) ---
    lps str = {0, 0, 0, nullptr};
    charsi chs;
    lcp_read lr;
    base_lcp::parse_cores_from_bam_record(bam_data, seqlen, 1, &str,
                                          0 /*unused*/, chs, lr, false);

    // Level-1 cores from parse_cores_from_bam_record should match
    // parse1_frombam exactly (before deepening, both are level 1).
    // parse_cores_from_bam_record deepens to lcp_level=1 which is a no-op on
    // level-1 data, so str.cores should equal direct_cores.
    INFO("seqlen=", seqlen, " direct=", nd, " via_record=", str.size);
    REQUIRE(nd == (uint64_t)str.size);
    for (uint64_t i = 0; i < nd; ++i) {
        INFO("core ", i);
        CHECK(direct_cores[i].start    == str.cores[i].start);
        CHECK(direct_cores[i].end      == str.cores[i].end);
        CHECK(direct_cores[i].label    == str.cores[i].label);
        CHECK(direct_cores[i].bit_rep  == str.cores[i].bit_rep);
        CHECK(direct_cores[i].bit_size == str.cores[i].bit_size);
    }

    free(direct_cores);
    free_lps(&str);
    bam_destroy1(rec);
    bam_hdr_destroy(hdr);
    sam_close(fp);
}

TEST_CASE("base_lcp::parse_cores_from_bam_record forward on synthetic data") {
    // Verify the fix without needing a BAM file — use synthetic sequences.
    const char *sequences[] = {
        "CAGTCAGTCAGTCAGTCAGT",
        "ACGTACGTACGTACGTACGT",
        "GATTACAGATTAGATTACAGATTAGATTACAGATTA",
        "TTTACGTTTTTGCAGTCAGTCATTTGCATTTTTTT",
        "GCATTTTTTTTTTTTTTTTCGATCGATCGA",
    };

    for (const char *seq : sequences) {
        size_t len = strlen(seq);
        INFO("sequence: ", seq);

        // Pack to nibble format
        uint8_t *nibble = new uint8_t[(len + 1) / 2]();
        ascii_to_bam_nibble(seq, len, nibble);

        // Direct parse1_frombam
        struct core direct_cores[256];
        uint64_t nd = base_lcp::parse1_frombam(nibble, (int64_t)len,
                                               direct_cores, 0);

        // Through parse_cores_from_bam_record (forward)
        lps str = {0, 0, 0, nullptr};
        charsi chs;
        lcp_read lr;
        base_lcp::parse_cores_from_bam_record(nibble, len, 1, &str,
                                              0 /*unused*/, chs, lr, false);

        REQUIRE(nd == (uint64_t)str.size);
        for (uint64_t i = 0; i < nd; ++i) {
            INFO("core ", i);
            CHECK(direct_cores[i].start    == str.cores[i].start);
            CHECK(direct_cores[i].end      == str.cores[i].end);
            CHECK(direct_cores[i].label    == str.cores[i].label);
            CHECK(direct_cores[i].bit_rep  == str.cores[i].bit_rep);
            CHECK(direct_cores[i].bit_size == str.cores[i].bit_size);
        }

        free_lps(&str);
        delete[] nibble;
    }
}

TEST_CASE("base_lcp::init_lps_bam does not leak on repeated calls") {
    // Call init_lps_bam twice on the same lps — second call should reuse
    // the buffer via prepare_lps instead of leaking.
    const char *seq = "CAGTCAGTCAGTCAGTCAGT";
    size_t len = strlen(seq);
    uint8_t *nibble = new uint8_t[(len + 1) / 2]();
    ascii_to_bam_nibble(seq, len, nibble);

    lps str = {0, 0, 0, nullptr};
    base_lcp::init_lps_bam(&str, nibble, (int)len);
    uint64_t first_size = str.size;
    CHECK(first_size > 0);

    // Second call — should not crash or leak (prepare_lps reuses buffer)
    base_lcp::init_lps_bam(&str, nibble, (int)len);
    CHECK(str.size == (int)first_size);

    free_lps(&str);
    delete[] nibble;
}
