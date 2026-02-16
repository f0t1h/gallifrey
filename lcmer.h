#pragma once
#include "htslib/hts.h"
#ifndef LCMER_H
#define LCMER_H
#include <htslib/sam.h>

#include <algorithm>
#include <array>
#include <cstdint>
#include <cstring>
#include <gtl/vector.hpp>
#include <iterator>
#include <string>

#include "ext/lcptools/lcptools_ho.h"
#include "nibble2base.hpp"
using SequenceContainerType = CrumbledSubSeq64;

// Makes sure we have enough scratch space
template <class T, double GF = 2, size_t ALIGN_BLOCK_BITS = 5>
struct charsi_base {
    static constexpr size_t align_block_bits = ALIGN_BLOCK_BITS;
    static constexpr size_t align_block_mask = (1 << align_block_bits) - 1;
    T *data;
    size_t capacity;
    charsi_base() : data{NULL}, capacity{0} {}
    // Rule of 3 just delete
    charsi_base(const charsi_base &) = delete;
    charsi_base &operator=(const charsi_base &) = delete;
    // Rule of 5
    charsi_base(charsi_base &&other) noexcept
        : data{other.data}, capacity{other.capacity} {
        other.data = nullptr;
        other.capacity = 0;
    }
    charsi_base &operator=(charsi_base &&other) noexcept {
        if (this != &other) {
            free(data);
            data = other.data;
            capacity = other.capacity;
            other.data = nullptr;
            other.capacity = 0;
        }
        return *this;
    }

    ~charsi_base() { free(data); }
    int require(size_t cap) {
        if (cap > capacity) {
            capacity = cap * GF;
            free(data);
            size_t target_cap = capacity * sizeof(T);
            target_cap = (target_cap | align_block_mask) + 1;
            data = (decltype(data))std::aligned_alloc(32, target_cap);
        }
        return data == NULL;
    }
};
using charsi = charsi_base<uint8_t, 2.0, 5>;
int onehotto2bitDNA2(const uint8_t *data, char *ret);
struct sv {
    char *from;
    char *to;
    sv() : from(0), to(0) {}
    sv(char *from, char *to) : from(from), to(to) {}

    void expand_until(char c) {
        while (*to != 0 && *to != c) {
            to++;
        }
    }
    void skip_next(char c) {
        expand_until(c);
        if (*to) ++to;
        from = to;
    }
    template <int N, char DEL>
    void get() {
        for (int i = 0; i < N; ++i) {
            skip_next(DEL);
        }
        expand_until(DEL);
    }
    // operator std::string() { return {from, to}; }
};

struct lcp_read {
    gtl::vector<uint32_t> cores;
    gtl::vector<std::pair<uint32_t, uint32_t>> pos;
    const uint8_t *data;
    void reserve(size_t size) {
        cores.reserve(size);
        pos.reserve(size);
    }
    void clear() {
        cores.clear();
        pos.clear();
        data = NULL;
    }
};

template <int N>
struct lcmer {
    // static uint32_t MaskedMurmur3(const std::array<uint32_t, N> &key,
    //                               uint32_t seed) {
    //     //        const uint8_t *data = (const uint8_t *)key;

    //     uint32_t h1 = seed;

    //     const uint32_t c1 = 0xcc9e2d51;
    //     const uint32_t c2 = 0x1b873593;

    //     // Body: Process blocks of 4 bytes at a time
    //     const uint32_t *blocks = key.data() + N;

    //     for (int i = -N; i; i++) {
    //         uint32_t k1 = lcmer<N>::mask[i + N] ? blocks[i] : 0;

    //         k1 *= c1;
    //         k1 = (k1 << 15) | (k1 >> (32 - 15));
    //         k1 *= c2;

    //         h1 ^= k1;
    //         h1 = (h1 << 15) | (h1 >> (32 - 15));
    //         h1 = h1 * 5 + 0xe6546b64;
    //     }
    //     // Finalization: Mix the hash to ensure the last few bits are fully
    //     // mixed
    //     // h1 ^= (N*4);

    //     /* fmix32 */
    //     h1 ^= h1 >> 16;
    //     h1 *= 0x85ebca6b;
    //     h1 ^= h1 >> 13;
    //     h1 *= 0xc2b2ae35;
    //     h1 ^= h1 >> 16;
    //     return h1;
    // }
    static constexpr size_t size = N;
    std::array<uint32_t, N> data;

    lcmer() {}
    lcmer(const lcmer<N> &other) {
        memcpy(this->data.data(), other.data.data(), N * sizeof(uint32_t));
    }
    lcmer(std::forward_iterator auto it) {
        std::copy(it, it + N, std::begin(data));
    }

    lcmer<N> operator=(const lcmer<N> &other) {
        memcpy(this->data.data(), other.data.data(), N * sizeof(uint32_t));
        return *this;
    }
    bool operator==(const lcmer &other) const {
        return memcmp(data.data(), other.data.data(), N * sizeof(uint32_t)) ==
               0;
        // for (int i = 0; i < N; ++i) {
        //     if (data[i] != other.data[i]) return false;
        // }
        // return true;
    }
    bool operator<(const lcmer &other) const {
        for (int i = 0; i < N; ++i) {
            if (data[i] != other.data[i]) {
                if (data[i] < other.data[i]) {
                    return true;
                } else {
                    return false;
                }
            }
        }
        return false;
    }
    void print_name(FILE *out) const {
        fprintf(out, "%0x", data[0]);
        for (int z = 1; z < N; ++z) {
            fprintf(out, "%0x", data[z]);
        }
    }

    void set0(uint32_t new_core) { data[N - 1] = new_core; }
    void setN(uint32_t new_core) { data[0] = new_core; }
    static lcmer<N> left_shift(lcmer<N> l) {
        for (int i = 1; i < N; ++i) {
            l.data[i - 1] = l.data[i];
        }
        return l;
    }

    static lcmer<N> right_shift(lcmer<N> l) {
        for (int i = N - 1; i > 0; --i) {
            l.data[i] = l.data[i - 1];
        }
        return l;
    }

    uint32_t update_lcmer(uint32_t new_core) {
        uint32_t old_core = data[0];
        for (int i = 1; i < N; ++i) {
            data[i - 1] = data[i];
        }
        data[N - 1] = new_core;
        return old_core;
    }
};

template <int N>
struct std::hash<lcmer<N>> {
    std::size_t operator()(const lcmer<N> &l) const noexcept {
        return MurmurHash3_32(&l.data, N * sizeof(uint32_t), 12312342);
    }
};

template <typename INT_TYPE = uint32_t, INT_TYPE BIT_WIDTH = 2>
struct KmerViewOnInts {
    using int_type = INT_TYPE;
    int_type k;
    int_type KMASK;
    int_type value;
    uint8_t prev;

    static constexpr INT_TYPE BIT_MASK = -1;
    // static constexpr INT_TYPE BIT_MASK = ((INT_TYPE)1 << BIT_WIDTH) - 1;

    KmerViewOnInts(int_type k)
        : k{k},
          KMASK{((int_type)1 << (BIT_WIDTH * k)) - 1},
          value{0},
          prev{0xff} {}

    int fill(const uint8_t *data_address) {
        int_type i;
        for (i = 0; i < k; i++) {
            move_get(data_address, i);
        }
        return i;
    }

    int fill_rle(const uint8_t *data_address, int_type cap) {
        int_type i;
        prev = 0xff;
        int_type it = move_get_rle(data_address, 0, cap);
        for (i = 0; i < k && it < cap; i++) {
            it = move_get_rle(data_address, it, cap);
        }
        return it;
    }

    int_type move_get(const uint8_t *data_address, int_type index) {
        value <<= BIT_WIDTH;
        value &= KMASK;
        value |= (data_address[index] & BIT_MASK);
        return index + 1;
    }

    int_type move_get_rle2(uint8_t prev, const uint8_t *data_address,
                           int_type index, int_type cap) {
        while (index < cap && prev == data_address[index]) {
            index++;
        }
        value <<= BIT_WIDTH;
        value &= KMASK;
        value |= (data_address[index] & BIT_MASK);
        return index + 1;
    }

    int_type move_get_rle(const uint8_t *data_address, int_type index,
                          int_type cap) {
        const uint8_t *p = data_address + index;
        const uint8_t *const end = data_address + cap;

        while (p < end && *p == prev) {
            ++p;
        }
        if (p == end) {
            return cap;
        }
        prev = *p;
        value = ((value << BIT_WIDTH) & KMASK) | (*p & BIT_MASK);
        return static_cast<int_type>((p - data_address) + 1);
    }
    constexpr int_type get() const noexcept { return value; }
};

template <class INT_TYPE = uint32_t>
struct KmerViewOnCrumbs {
    using int_type = INT_TYPE;
    int_type k;
    int_type kmask;
    int_type value;
    uint8_t prev;
    static constexpr uint8_t mask = 3U;
    KmerViewOnCrumbs(int_type k)
        : k{k}, kmask{((int_type)1 << (2 * k)) - 1}, value{0}, prev{255} {}

    int fill_rle(const uint8_t *data_address, int_type cap) {
        int_type i;
        int_type it = move_get_rle(data_address, 0, cap);
        for (i = 0; i < k && it < cap; i++) {
            it = move_get_rle(data_address, it, cap);
        }
        return it;
    }

    int_type move_get_rle(const uint8_t *data_address, int_type index,
                          int_type cap) {
        size_t byte_idx = index / 4;
        size_t crum_idx = index % 4;
        uint8_t at = (data_address[byte_idx] >> (crum_idx * 2)) & mask;
        while (index < cap && at == prev) {
            ++index;
            byte_idx = index / 4;
            crum_idx = index % 4;
            at = (data_address[byte_idx] >> (crum_idx * 2)) & mask;
        }
        if (index == cap) {
            return index;
        }
        value = ((value << 2) & kmask) | at;
        prev = at;
        return index + 1;
    }
    constexpr int_type get() const noexcept { return value; }
};

// Crumb view for 64-bit packed sequences (2 bits/base, 32 bases/word).
template <class INT_TYPE = uint64_t>
struct KmerViewOnCrumbs64 {
    using int_type = INT_TYPE;
    int_type k;
    int_type kmask;
    int_type value;
    uint8_t prev;
    static constexpr uint8_t mask = 3U;

    KmerViewOnCrumbs64(int_type k)
        : k{k}, kmask{((int_type)1 << (2 * k)) - 1}, value{0}, prev{255} {}

    static uint8_t base_at(const uint64_t *data_words, size_t index) {
        const size_t bit_index = index * 2;
        return static_cast<uint8_t>(
            (data_words[bit_index >> 6] >> (bit_index & 63)) & mask);
    }

    int fill_rle(const uint64_t *data_words, int_type cap,
                 int_type skip_first_n) {
        int_type i;
        int_type it = move_get_rle(data_words, skip_first_n, cap);
        for (i = 0; i < k && it < cap; i++) {
            it = move_get_rle(data_words, it, cap);
        }
        return it;
    }

    int_type move_get_rle(const uint64_t *data_words, int_type index,
                          int_type cap) {
        uint8_t at = base_at(data_words, index);
        while (index < cap && at == prev) {
            ++index;
            if (index == cap) {
                break;
            }
            at = base_at(data_words, index);
        }
        if (index == cap) {
            return index;
        }
        value = ((value << 2) & kmask) | at;
        prev = at;
        return index + 1;
    }
    constexpr int_type get() const noexcept { return value; }
    char *to_string() const {
        static const char *lut = "ACGT";
        char *buffer = malloc(k + 1);
        for (int i = 0; i < k; i++) {
            buffer[i] = lut[(value >> (2 * (k - i - 1))) & 3];
        }
        buffer[k] = '\0';
        return buffer;
    }
};

template <int _K, typename INT_TYPE = uint32_t, int BIT_WIDTH = 2>
struct KmerViewOnNibbles {
    static constexpr int tab[] = {4, 0, 1, 4, 2, 4, 4, 4,
                                  3, 4, 4, 4, 4, 4, 4, 4};

    using int_type = INT_TYPE;
    static constexpr int k = _K;
    int_type value;
    static constexpr int_type KMASK = ((int_type)1 << (BIT_WIDTH * k)) - 1;
    KmerViewOnNibbles() : value{0} {}

    int fill(const uint8_t *data_address) {
        int i;
        for (i = 0; i < k; i++) {
            move_get(data_address, i);
        }
        return i;
    }
    int_type move_get(const uint8_t *data_address, int index) {
        value <<= BIT_WIDTH;
        if constexpr (KMASK != -1) value &= KMASK;
        if constexpr (BIT_WIDTH == 2) {
            value |= tab[bam_seqi(data_address, index)];
        } else {
            value |= bam_seqi(data_address, index);
        }
        return value;
    }

    constexpr int_type get() const noexcept { return value; }
};

template <int _K, typename INT_TYPE = uint32_t, int BIT_WIDTH = 2>
struct KmerViewOnBytes {
    static constexpr const int tab[256] = {
        4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
        4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
        4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 0, 4, 1, 4,
        4, 4, 2, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 3, 3, 4, 4, 4, 4, 4, 4,
        4, 4, 4, 4, 4, 0, 4, 1, 4, 4, 4, 2, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
        4, 3, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
        4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
        4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
        4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
        4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
        4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4};

    using int_type = INT_TYPE;
    static constexpr int k = _K;
    int_type value;
    static constexpr int_type KMASK = ((int_type)1 << (BIT_WIDTH * k)) - 1;
    KmerViewOnBytes() : value{0} {}

    int fill(const char *data) {
        int i;
        for (i = 0; i < k; i++) {
            move_get(data, i);
        }
        return i;
    }
    int_type move_get(const char *data, int index) {
        value <<= BIT_WIDTH;
        if constexpr (KMASK != -1) value &= KMASK;
        if constexpr (BIT_WIDTH == 2) {
            value |= tab[static_cast<int>(data[index])];
        } else {
            value |= data[index];
        }
        return value;
    }

    constexpr int_type get() const noexcept { return value; }
};

void init_core1bam(struct core *cr, const uint8_t *data, uint64_t distance,
                   uint64_t si, uint64_t ei);

template <typename _T, int WIDTH = 2>
void init_core_kmer(struct core *cr, _T Kval, uint64_t si, uint64_t ei, int k) {
    cr->start = si;
    cr->end = ei;
    cr->label = 0;
    cr->label |= Kval;
    cr->bit_rep = 0x8000000000000000 | cr->label;
    cr->bit_size = WIDTH * k;
}

template <typename _T>
void init_core_kmer_xor(struct core *cr, _T Kval, uint64_t si, uint64_t ei,
                        int k) {
    cr->start = si;
    cr->end = ei;
    cr->label = 0;
    cr->label |= (Kval ^ reverse_complement_kmer(Kval, k));
    cr->bit_rep = 0x8000000000000000 | cr->label;
    cr->bit_size = 2 * k;
}

uint8_t *onehotto2bitDNA(const uint8_t *data, size_t len);
uint8_t *onehotto2bitDNA_bivec(const uint8_t *data, size_t len);

void onehotto2bitDNA_bivec_chs(const uint8_t *data, size_t len, charsi &c);
void onehotto2bitDNA_bivec_chs_reverse(const uint8_t *data, size_t len,
                                       charsi &c);

static void init_core1hot(struct core *cr, const uint8_t *data, uint64_t from,
                          uint64_t distance, uint64_t start_index,
                          uint64_t end_index) {
    cr->start = start_index;
    cr->end = end_index;
    cr->label = 0;
    cr->label |= ((distance - 2) << 6);
    cr->label |= (seq_nt16_int[bam_seqi(data, from)] << 4);
    cr->label |= (seq_nt16_int[bam_seqi(data, from + distance - 2)] << 2);
    cr->label |= (seq_nt16_int[bam_seqi(data, from + distance - 1)]);
    cr->bit_rep = 0x8000000000000000 | cr->label;
    cr->bit_size = 2 * distance;
}

inline void prepare_lps(lps *lps_ptr, size_t len) {
    lps_ptr->level = 1;
    if (lps_ptr->cores == nullptr) {
        lps_ptr->size = 0;
        lps_ptr->capacity = (len / CONSTANT_FACTOR);
        lps_ptr->cores = (struct core *)malloc((len / CONSTANT_FACTOR) *
                                               sizeof(struct core));
    } else {
        if (lps_ptr->capacity < len / CONSTANT_FACTOR) {
            lps_ptr->capacity = len / CONSTANT_FACTOR;
            lps_ptr->cores = (struct core *)realloc(
                lps_ptr->cores,
                (len / CONSTANT_FACTOR) * sizeof(struct core));
        }
    }
}
struct base_lcp {
    static uint64_t parse1_frombam(const uint8_t *data, int64_t len,
                                   struct core *cores, uint64_t offset) {
        int64_t begin = 0;
        int64_t end = len;
        int64_t it1 = begin;
        int64_t it2 = end;
        int64_t core_index = 0;
        int64_t last_invalid_char_index = -1;

        // find lcp cores
        for (; it1 + 2 < end; it1++) {
            // skip invalid character
            int base = seq_nt16_int[bam_seqi(data, it1)];
            int next_base = seq_nt16_int[bam_seqi(data, it1 + 1)];
            int next_next_base = seq_nt16_int[bam_seqi(data, it1 + 2)];

            if (base == -1) {
                last_invalid_char_index = it1 - begin;
                continue;
            }

            if (base == next_base) {
                continue;
            }

            // check for RINT core
            if (next_base == next_next_base) {
                // count middle characters
                uint32_t middle_count = 1;
                int64_t temp = it1 + 2;
                while (temp < end && seq_nt16_int[bam_seqi(data, temp - 1)] ==
                                         seq_nt16_int[bam_seqi(data, temp)]) {
                    temp++;
                    middle_count++;
                }
                if (temp != end) {
                    // check if there is any SSEQ cores left behind
                    if (it2 < it1 &&
                        last_invalid_char_index < it2 - begin - 1) {
                        init_core1hot(&(cores[core_index]), data, it2 - 1,
                                      it1 - it2 + 2, it2 - begin - 1 + offset,
                                      it1 - begin + 1 + offset);
                        core_index++;
                    }

                    // create RINT core
                    it2 = it1 + 2 + middle_count;
                    init_core1hot(&(cores[core_index]), data, it1,
                                  2 + middle_count, it1 - begin + offset,
                                  it2 - begin + offset);
                    core_index++;

                    continue;
                }
            }

            if (base > next_base && next_base < next_next_base) {
                // check if there is any SSEQ cores left behind
                if (it2 < it1 && last_invalid_char_index < it2 - begin - 1) {
                    init_core1hot(&(cores[core_index]), data, it2 - 1,
                                  it1 - it2 + 2, it2 - begin - 1 + offset,
                                  it1 - begin + 1 + offset);
                    core_index++;
                }

                // create LMIN core
                it2 = it1 + 3;
                init_core1hot(&(cores[core_index]), data, it1, 3,
                              it1 - begin + offset, it2 - begin + offset);
                core_index++;

                continue;
            }

            if (begin == it1) {
                continue;
            }

            // check for LMAX
            if (it1 + 3 < end && base < next_base &&
                next_base > next_next_base &&
                seq_nt16_int[bam_seqi(data, it1 - 1)] <= base &&
                next_next_base >= seq_nt16_int[bam_seqi(data, it1 + 3)]) {
                // check if there is any SSEQ cores left behind
                if (it2 < it1 && last_invalid_char_index < it2 - begin - 1) {
                    init_core1hot(&(cores[core_index]), data, it2 - 1,
                                  it1 - it2 + 2, it2 - begin - 1 + offset,
                                  it1 - begin + 1 + offset);
                    core_index++;
                }

                // create LMAX core
                it2 = it1 + 3;
                init_core1hot(&(cores[core_index]), data, it1, 3,
                              it1 - begin + offset, it2 - begin + offset);
                core_index++;

                continue;
            }
        }

        return core_index;
    }
    static uint64_t parse1_from_SCT_vec(const gtl::vector<SequenceContainerType> &sseqs, struct core *cores, uint64_t offset){
        if (sseqs.empty()) return 0;

        int64_t total_len = 0;
        for (const auto &s : sseqs) total_len += static_cast<int64_t>(s.len);
        if (total_len < 3) return 0;

        // Sequential cursor over the concatenated fragments
        size_t cur_frag = 0, cur_local = 0;

        auto cursor_valid = [&]() -> bool {
            return cur_frag < sseqs.size();
        };

        auto read_at_cursor = [&]() -> int {
            const auto &f = sseqs[cur_frag];
            size_t bi = (static_cast<size_t>(f.skip_first_n) + cur_local) * 2;
            return static_cast<int>((f.data[bi >> 6] >> (bi & 63)) & 0x3);
        };

        auto advance_cursor = [&]() {
            if (++cur_local >= sseqs[cur_frag].len) {
                cur_frag++;
                cur_local = 0;
            }
        };

        auto read_and_advance = [&]() -> int {
            int b = read_at_cursor();
            advance_cursor();
            return b;
        };

        // Init core from three explicit base values (mirrors init_core1hot)
        auto init_core_from_bases = [](struct core *cr, int first_base,
                                       int penult_base, int last_base,
                                       uint64_t distance, uint64_t start_index,
                                       uint64_t end_index) {
            cr->start = start_index;
            cr->end = end_index;
            cr->label = 0;
            cr->label |= ((distance - 2) << 6);
            cr->label |= (static_cast<uint64_t>(first_base) << 4);
            cr->label |= (static_cast<uint64_t>(penult_base) << 2);
            cr->label |= static_cast<uint64_t>(last_base);
            cr->bit_rep = 0x8000000000000000 | cr->label;
            cr->bit_size = 2 * distance;
        };

        // Sliding window: buf[0]=it1-1, buf[1]=it1, buf[2]=it1+1,
        //                  buf[3]=it1+2, buf[4]=it1+3
        int buf[5];
        buf[0] = -1; // sentinel, no base before position 0
        buf[1] = read_and_advance();
        buf[2] = read_and_advance();
        buf[3] = read_and_advance();
        buf[4] = cursor_valid() ? read_and_advance() : -1;

        int64_t begin = 0;
        int64_t end = total_len;
        int64_t it1 = begin;
        int64_t it2 = end;
        int64_t core_index = 0;
        int cached_it2_base = -1; // base at position it2-1, set when it2 changes

        // find lcp cores
        for (; it1 + 2 < end; it1++) {
            // shift the window (skip on first iteration — already filled)
            if (it1 > begin) {
                buf[0] = buf[1];
                buf[1] = buf[2];
                buf[2] = buf[3];
                buf[3] = buf[4];
                buf[4] = cursor_valid() ? read_and_advance() : -1;
            }

            int base = buf[1];
            int next_base = buf[2];
            int next_next_base = buf[3];

            if (base == next_base) {
                continue;
            }

            // check for RINT core
            if (next_base == next_next_base) {
                // count middle characters using buffer then cursor fork
                uint32_t middle_count = 1;
                int64_t temp = it1 + 2;
                int rint_last_base = -1;

                // save cursor so the forward scan doesn't consume it
                size_t save_frag = cur_frag, save_local = cur_local;

                int prev_scan = buf[2]; // base[it1+1]
                while (temp < end) {
                    int curr_scan;
                    int64_t buf_idx = temp - it1 + 1;
                    if (buf_idx <= 4) {
                        curr_scan = buf[buf_idx];
                    } else {
                        if (!cursor_valid()) break;
                        curr_scan = read_and_advance();
                    }
                    if (prev_scan != curr_scan) {
                        rint_last_base = curr_scan;
                        break;
                    }
                    temp++;
                    middle_count++;
                    prev_scan = curr_scan;
                }

                // restore cursor
                cur_frag = save_frag;
                cur_local = save_local;

                if (temp != end) {
                    // check if there is any SSEQ cores left behind
                    if (it2 < it1) {
                        init_core_from_bases(&(cores[core_index]),
                                             cached_it2_base, buf[0], buf[1],
                                             it1 - it2 + 2,
                                             it2 - begin - 1 + offset,
                                             it1 - begin + 1 + offset);
                        core_index++;
                    }

                    // create RINT core
                    it2 = it1 + 2 + middle_count;
                    cached_it2_base = rint_last_base;
                    init_core_from_bases(&(cores[core_index]),
                                         buf[1], next_base, rint_last_base,
                                         2 + middle_count,
                                         it1 - begin + offset,
                                         it2 - begin + offset);
                    core_index++;

                    continue;
                }
            }

            if (base > next_base && next_base < next_next_base) {
                // check if there is any SSEQ cores left behind
                if (it2 < it1) {
                    init_core_from_bases(&(cores[core_index]),
                                         cached_it2_base, buf[0], buf[1],
                                         it1 - it2 + 2,
                                         it2 - begin - 1 + offset,
                                         it1 - begin + 1 + offset);
                    core_index++;
                }

                // create LMIN core
                it2 = it1 + 3;
                cached_it2_base = buf[3];
                init_core_from_bases(&(cores[core_index]),
                                     buf[1], buf[2], buf[3],
                                     3, it1 - begin + offset,
                                     it2 - begin + offset);
                core_index++;

                continue;
            }

            if (begin == it1) {
                continue;
            }

            // check for LMAX
            if (it1 + 3 < end && base < next_base &&
                next_base > next_next_base &&
                buf[0] <= base &&
                next_next_base >= buf[4]) {
                // check if there is any SSEQ cores left behind
                if (it2 < it1) {
                    init_core_from_bases(&(cores[core_index]),
                                         cached_it2_base, buf[0], buf[1],
                                         it1 - it2 + 2,
                                         it2 - begin - 1 + offset,
                                         it1 - begin + 1 + offset);
                    core_index++;
                }

                // create LMAX core
                it2 = it1 + 3;
                cached_it2_base = buf[3];
                init_core_from_bases(&(cores[core_index]),
                                     buf[1], buf[2], buf[3],
                                     3, it1 - begin + offset,
                                     it2 - begin + offset);
                core_index++;

                continue;
            }
        }

        return core_index;
    }
    static void init_lps_crumbled(
        lps *lps_ptr, const gtl::vector<SequenceContainerType> &sseqs,
        const auto &syncmer_functor) {
        size_t len = 0;

        for (auto &sseq : sseqs) {
            len += sseq.len;
        }
        prepare_lps(lps_ptr, len);
        lps_ptr->size =
            parse1_from_SCT_vec(sseqs, lps_ptr->cores, 0);
        assert(lps_ptr->size < lps_ptr->capacity);
    }
    static void init_lps_bam(struct lps *lps_ptr, const uint8_t *str, int len) {
        prepare_lps(lps_ptr, len);
        lps_ptr->size = parse1_frombam(str, len, lps_ptr->cores, 0);
    }

    static void parse_cores_from_bam_record(const uint8_t *data, size_t len,
                                            int lcp_level, lps *str,
                                            const auto &syncmer_functor,
                                            charsi &chs, lcp_read &lr,
                                            bool reversed = false) {
        // parse1_frombam (and init_core1hot) read bases via bam_seqi() +
        // seq_nt16_int[], so they require nibble-packed data — pass directly.
        (void)syncmer_functor;
        (void)chs;
        (void)reversed; // TODO: reverse-complement not yet supported for pure_lcp BAM path
        init_lps_bam(str, data, len);
        lps_deepen(str, lcp_level);

        for (int i = 0; i < str->size; ++i) {
            lr.cores.emplace_back(str->cores[i].label);
            lr.pos.emplace_back(str->cores[i].start, str->cores[i].end);
        }
        lr.data = data;
    }
};

struct syncmer_lcp {
    static uint64_t parse1_bam_syncmers(const uint8_t *data, int64_t len,
                                        struct lps *str,
                                        const auto &syncmer_functor) {
        int core_index = 0;
        const auto k = syncmer_functor.k();
        KmerViewOnInts<uint64_t, 2> kmer_view{k};
        int64_t it1 = kmer_view.fill_rle(data, len);
        if (str->capacity < len * CONSTANT_FACTOR) {
            str->capacity = 2 * len * CONSTANT_FACTOR;
            str->cores = (struct core *)realloc(
                str->cores, (str->capacity) * sizeof(struct core));
        }
        const int64_t limit = len - k;
        auto *cores = str->cores;
        for (; it1 < limit;) {
            const auto kmer = kmer_view.get();
            if (syncmer_functor(kmer)) {
                init_core_kmer<uint32_t>(&(cores[core_index]), kmer, it1,
                                         it1 + k, k);
                core_index++;
            }
            it1 = kmer_view.move_get_rle(data, it1, len);
        }
        return core_index;
    }
    static uint64_t parse1_CrumbledSubSeqs_syncmers(
        const gtl::vector<SequenceContainerType> &sseqs, struct lps *str,
        const auto &syncmer_functor) {
        int core_index = 0;
        if (sseqs.empty()) {
            fprintf(stderr, "Error: Empty sequence list\n");
            return core_index;
        }
        KmerViewOnCrumbs64<uint64_t> kmer_view{syncmer_functor.k()};
        size_t it1 = kmer_view.fill_rle(sseqs[0].data,
                                        sseqs[0].len + sseqs[0].skip_first_n,
                                        sseqs[0].skip_first_n);
        size_t l_so_far = 0;
        for (size_t i = 0; i < sseqs.size(); ++i) {
            for (; it1 < sseqs[i].len + sseqs[i].skip_first_n;) {
                if (syncmer_functor(kmer_view.get())) {
                    init_core_kmer<uint32_t>(
                        &(str->cores[core_index]), kmer_view.get(),
                        l_so_far + it1 - sseqs[i].skip_first_n,
                        l_so_far + it1 + syncmer_functor.k() -
                            sseqs[i].skip_first_n,
                        syncmer_functor.k());
                    core_index++;
                }
                it1 = kmer_view.move_get_rle(
                    sseqs[i].data, it1, sseqs[i].len + sseqs[i].skip_first_n);
            }
            if (sseqs.size() > i + 1) {
                it1 = sseqs[i + 1].skip_first_n;
            }
            l_so_far += sseqs[i].len;
        }
        return core_index;
    }

    uint64_t parse1_bam_lmax(const uint8_t *data, int64_t len,
                             struct core *cores);

    /// @brief Initializes an lps struct with a given sequence in the bam
    /// compressed
    // void init_lps_bam(lps *lps_ptr, const uint8_t *data, size_t len);
    //

    static void init_lps_bam(lps *lps_ptr, const uint8_t *data, size_t len,
                             const auto &syncmer_functor) {
        prepare_lps(lps_ptr, len);
        lps_ptr->size =
            parse1_bam_syncmers(data, len, lps_ptr, syncmer_functor);
    }

    static void init_lps_crumbled(
        lps *lps_ptr, const gtl::vector<SequenceContainerType> &sseqs,
        const auto &syncmer_functor) {
        size_t len = 0;

        for (auto &sseq : sseqs) {
            len += sseq.len;
        }
        prepare_lps(lps_ptr, len);
        lps_ptr->size =
            parse1_CrumbledSubSeqs_syncmers(sseqs, lps_ptr, syncmer_functor);
        assert(lps_ptr->size < lps_ptr->capacity);
    }

    static void parse_cores_from_bam_record(const uint8_t *data, size_t len,
                                            int lcp_level, lps *str,
                                            const auto &syncmer_functor,
                                            charsi &chs, lcp_read &lr,
                                            bool reversed = false) {
        chs.require(len);
        if (reversed) {
            onehotto2bitDNA_bivec_chs_reverse(data, len, chs);
        } else {
            onehotto2bitDNA_bivec_chs(data, len, chs);
        }

        init_lps_bam(str, chs.data, len, syncmer_functor);
        lps_deepen(str, lcp_level);

        for (int i = 0; i < str->size; ++i) {
            lr.cores.emplace_back(str->cores[i].label);
            lr.pos.emplace_back(str->cores[i].start, str->cores[i].end);
        }
        lr.data = data;
    }
};
inline void init_lps_crumbled(
    lps *lps_ptr, const gtl::vector<SequenceContainerType> &sseqs,
    const auto &syncmer_functor) {
    if(syncmer_functor.pure_lcp){
        base_lcp::init_lps_crumbled(lps_ptr, sseqs, syncmer_functor);
    }
    else{
        syncmer_lcp::init_lps_crumbled(lps_ptr, sseqs, syncmer_functor);
    }
}
inline void parse_cores_from_bam_record(const uint8_t *data, size_t len,
                                        int lcp_level, lps *str,
                                        const auto &syncmer_functor,
                                        charsi &chs, lcp_read &lr,
                                        bool reversed = false) {
    if (syncmer_functor.pure_lcp) {
        base_lcp::parse_cores_from_bam_record(
            data, len, lcp_level, str, syncmer_functor, chs, lr, reversed);
    } else {
        syncmer_lcp::parse_cores_from_bam_record(
            data, len, lcp_level, str, syncmer_functor, chs, lr, reversed);
    }
}

template <int lanes = 4>
inline std::array<std::pair<uint8_t *, uint8_t *>, lanes> scalar_rlc(charsi &c1,
                                                                     charsi &c2,
                                                                     int len) {
    int lane_len = len / 4;

    uint8_t *rptrs[lanes + 1] = {
        // TODO alignments
        c1.data,
        c1.data + lane_len,
        c1.data + 2 * lane_len,
        c1.data + 3 * lane_len,
        c1.data + 4 * lane_len,
    };
    uint8_t *wptrs[lanes + 1] = {
        c2.data,
        c2.data + lane_len,
        c2.data + 2 * lane_len,
        c2.data + 3 * lane_len,
        c2.data + 4 * lane_len,
    };
    uint8_t *endings[lanes] = {0, 0, 0, 0};
    for (int l = 0; l < lanes; ++l) {
        uint8_t *_from = rptrs[l];
        uint8_t *_to = rptrs[l + 1];
        uint8_t *writ = wptrs[l];
        *writ = *_from;
        for (uint8_t *it = _from + 1; it < _to; ++it) {
            if (*it != *writ) {
                ++writ;
                *writ = *it;
            }
        }
        endings[l] = writ;
    }
    return {
        {wptrs[0], endings[0]},
        {wptrs[1], endings[1]},
        {wptrs[2], endings[2]},
        {wptrs[3], endings[3]},
    };
}
template <int lanes = 4>
inline uint64_t parse1_bam_syncmer_with_lanes(
    std::array<std::pair<uint8_t *, uint8_t *>, lanes> data, int64_t len,
    struct lps *str, const auto &syncmer_functor) {
    const auto k = syncmer_functor.k();

    if (str->capacity < len * CONSTANT_FACTOR) {
        str->capacity = 2 * len * CONSTANT_FACTOR;
        str->cores = (struct core *)realloc(
            str->cores, (str->capacity) * sizeof(struct core));
    }
    auto *cores = str->cores;
    for (int l = 0; l < lanes; ++l) {
        int core_index = 0;
        KmerViewOnInts<uint64_t, 2> kmer_view{k};
        int64_t it1 = [&data, len, l, k, &kmer_view]() {
            if (l == 0) {
                return kmer_view.fill(data[l].first,
                                      len);  // Fill from beginning
            } else {
                return kmer_view.fill(data[l - 1].second - k + 1,
                                      len);  // Fill from previous ending-K
            }
        }();
        const uint8_t *ending = data[l].second;
        for (uint8_t *it = data[l].first + it1; it != ending; ++it) {
            const auto kmer = kmer_view.get();
            if (syncmer_functor(kmer)) {
                init_core_kmer<uint32_t>(&(cores[core_index]), kmer, it1,
                                         it1 + k,
                                         k);  // TODO fix si, ei
                core_index++;
            }
            it1 = kmer_view.move_get(data, it1, len);
        }
    }
    return 0;
}
inline void parse_cores_from_bam_record2(const uint8_t *data, size_t len,
                                         int lcp_level, lps *str,
                                         const auto &syncmer_functor,
                                         charsi &chs, charsi &chs2,
                                         lcp_read &lr) {
    chs.require(len);
    chs2.require(len);
    onehotto2bitDNA_bivec_chs(data, len, chs);
    std::array<std::pair<uint8_t *, uint8_t *>, 4> lane_boundary =
        scalar_rlc<4>(chs, chs2, len);
    prepare_lps(str, len);
    parse1_bam_syncmer_with_lanes(lane_boundary, len, str, syncmer_functor);
    lps_deepen(str, lcp_level);

    for (int i = 0; i < str->size; ++i) {
        lr.cores.emplace_back(str->cores[i].label);
        lr.pos.emplace_back(str->cores[i].start, str->cores[i].end);
    }
    lr.data = data;
}
lcp_read parse_cores_from_kseq(const char *data, size_t len, int lcp_level);
#endif
