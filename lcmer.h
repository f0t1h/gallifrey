#pragma once
#include <cstdint>

#include "nibble2base.hpp"
#ifndef LCMER_H
#define LCMER_H
#include <htslib/sam.h>

#include <algorithm>
#include <array>
#include <bit>
#include <cstring>
#include <gtl/vector.hpp>
#include <iterator>
#include <string>
#include <type_traits>

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
using pos_scratch = charsi_base<uint32_t, 2.0, 5>;
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

    uint32_t update_lcmer(uint32_t new_core,
                          [[maybe_unused]] uint32_t outcore) {
        uint32_t old_core = data[0];
        for (int i = 1; i < N; ++i) {
            data[i - 1] = data[i];
        }
        data[N - 1] = new_core;
        return old_core;
    }
};

struct litecmer {
    struct words_t {
        uint32_t hi;
        uint32_t lo;
    };
    words_t words;
    static constexpr uint32_t kPolyBase = 0x9e3779b1u;
    static_assert(sizeof(words_t) == sizeof(uint64_t));
    static_assert(std::is_trivially_copyable_v<words_t>);

    litecmer() : words{0, 0} {}
    litecmer(const litecmer &other) : words(other.words) {}
    litecmer(std::forward_iterator auto it, int n) : words{0, 0} {
        if (n > 1) {
            uint32_t h = 0;
            for (int i = 0; i < n - 1; ++i) {
                h = h * kPolyBase + *(it++);
            }
            words.hi = h;
            words.lo = *it;
        } else {
            words.lo = *it;
        }
    }

    litecmer &operator=(const litecmer &other) {
        words = other.words;
        return *this;
    }

    constexpr uint32_t hi() const { return words.hi; }
    constexpr uint32_t lo() const { return words.lo; }
    constexpr uint64_t native_bitcast_value() const {
        return std::bit_cast<uint64_t>(words);
    }
    constexpr uint64_t value() const {
        return (static_cast<uint64_t>(words.hi) << 32) |
               static_cast<uint64_t>(words.lo);
    }

    static uint32_t poly_pow(uint32_t base, int exp) {
        uint32_t out = 1;
        for (int i = 0; i < exp; ++i) {
            out *= base;
        }
        return out;
    }

    static litecmer make_poly(std::forward_iterator auto it, int n) {
        return litecmer(it, n);
    }

    static litecmer make_xor(std::forward_iterator auto it, int n) {
        litecmer l;
        for (int i = 0; i < n - 1; ++i) {
            l.words.hi ^= *(it++);
        }
        l.words.lo = *it;
        return l;
    }

    void update_lcmer_xor(uint32_t new_core, uint32_t out_core) {
        words.hi ^= words.lo;
        words.hi ^= out_core;
        words.lo = new_core;
    }

    // Order-sensitive rolling update for first N-1 cores.
    // Requires state initialized by make_poly().
    void update_lcmer_poly(uint32_t new_core, uint32_t out_core, int n) {
        if (n > 1) {
            const uint32_t kOutMul = poly_pow(kPolyBase, n - 2);
            uint32_t h = words.hi;
            h -= out_core * kOutMul;
            h = h * kPolyBase + words.lo;
            words.hi = h;
        }
        words.lo = new_core;
    }

    void update_lcmer(uint32_t new_core, uint32_t out_core, int n) {
        update_lcmer_poly(new_core, out_core, n);
    }

    constexpr auto operator<=>(const litecmer &other) const {
        return value() <=> other.value();
    };
    constexpr bool operator==(const litecmer &other) const {
        return value() == other.value();
    }
};
using LCmer = litecmer;

template <int N>
struct combined_lcmer {
    lcmer<N> full;
    litecmer lite;

    combined_lcmer(std::forward_iterator auto it) : full(it), lite(it, N) {}

    void update_lcmer(uint32_t new_core, uint32_t out_core) {
        full.update_lcmer(new_core, out_core);
        lite.update_lcmer(new_core, out_core, N);
    }
};

template <>
struct std::hash<litecmer> {
    std::size_t operator()(const litecmer &l) const noexcept {
        return l.value();
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

    static inline uint8_t base_at(const uint64_t *data_words, size_t index) {
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
    cr->bit_size = WIDTH * k;
    cr->bit_rep = uint128_t{0x8000000000000000, cr->label};
}

template <typename _T>
void init_core_kmer_xor(struct core *cr, _T Kval, uint64_t si, uint64_t ei,
                        int k) {
    cr->start = si;
    cr->end = ei;
    cr->label = 0;
    cr->label |= (Kval ^ reverse_complement_kmer(Kval, k));
    cr->bit_rep = uint128_t{0x8000000000000000, cr->label};
    cr->bit_size = 2 * k;
}

uint8_t *onehotto2bitDNA(const uint8_t *data, size_t len);
uint8_t *onehotto2bitDNA_bivec(const uint8_t *data, size_t len);

void onehotto2bitDNA_bivec_chs(const uint8_t *data, size_t len, charsi &c);
void onehotto2bitDNA_bivec_chs_reverse(const uint8_t *data, size_t len,
                                       charsi &c);

inline int64_t rle_compact(const uint8_t *data, int64_t len, uint8_t *rle_data,
                           uint32_t *rle_pos) {
    if (len <= 0) return 0;
    rle_data[0] = data[0];
    rle_pos[0] = 0;
    int64_t out = 1;
#if defined(__ARM_NEON) || defined(__ARM_NEON__)
    int64_t i = 1;
    for (; i + 16 <= len; i += 16) {
        uint8x16_t curr = vld1q_u8(data + i);
        uint8x16_t prev = vld1q_u8(data + i - 1);
        uint8x16_t neq = vmvnq_u8(vceqq_u8(curr, prev));
        uint8x8_t narrowed = vshrn_n_u16(vreinterpretq_u16_u8(neq), 4);
        uint64_t mask64 = vget_lane_u64(vreinterpret_u64_u8(narrowed), 0);
        if (mask64 == 0) continue;
        while (mask64) {
            int nibble_idx = __builtin_ctzll(mask64) >> 2;
            rle_data[out] = data[i + nibble_idx];
            rle_pos[out] = static_cast<uint32_t>(i + nibble_idx);
            ++out;
            mask64 &= ~(0x0FULL << (nibble_idx * 4));
        }
    }
    for (; i < len; ++i) {
        if (data[i] != data[i - 1]) {
            rle_data[out] = data[i];
            rle_pos[out] = static_cast<uint32_t>(i);
            ++out;
        }
    }
#else
    for (int64_t i = 1; i < len; ++i) {
        if (data[i] != data[i - 1]) {
            rle_data[out] = data[i];
            rle_pos[out] = static_cast<uint32_t>(i);
            ++out;
        }
    }
#endif
    return out;
}

// Fused BAM one-hot decode + RLE compaction in a single pass.
// Input: raw BAM 4-bit encoded data (2 bases per byte), length in bases.
// Output: rle_data/rle_pos filled with RLE-compacted 2-bit bases.

#if defined(__ARM_NEON) || defined(__ARM_NEON__)
// Binary reduction: compact 16 one-hot bytes into contiguous packed output.
//
// Each input byte is either 0 (dead/unchanged) or (lane_idx << 4) |
// one_hot_value. One-hot values {1,2,4,8} fit in the low nibble; lane index
// 0-15 in the high nibble.
//
// The reduction merges adjacent pairs in 3 rounds (16→8→4→2), packing kept
// bytes contiguously using widen + variable-shift + OR. Because one-hot values
// occupy distinct bit positions within their nibble, OR never causes collision.
//
// After 3 rounds: two uint64 values, each holding up to 8 packed bytes from the
// LSB. Unpack: extract index (>>4) and value (&0xF), decode one-hot→2-bit,
// store via NEON.
static inline int compact16_emit(
    uint8x16_t data_vec,             // 16 one-hot bytes (BAM nibble values)
    uint8x16_t neq_vec,              // 0xFF where changed, 0x00 where same
    const uint8x16_t &idx_template,  // precomputed (i<<4) for i=0..15
    const uint8x16_t &decode_tbl,    // one-hot→2-bit decode table (16 entries)
    uint8_t *rle_data,               // output: decoded 2-bit values
    uint32_t *rle_pos,               // output: base positions
    uint32_t pos_base  // position offset for this half (chunk_start or +16)
) {
    // Build combined byte per lane: (lane_idx << 4) | value for changed, 0 for
    // unchanged.
    uint8x16_t live = vandq_u8(data_vec, neq_vec);
    uint8x16_t combined = vorrq_u8(live, vandq_u8(neq_vec, idx_template));

    // --- Round 1: 16 bytes → 8 uint16 (merge adjacent pairs) ---
    // Deinterleave: even[k] = combined[2k], odd[k] = combined[2k+1]
    uint8x8x2_t uzp = vuzp_u8(vget_low_u8(combined), vget_high_u8(combined));
    uint8x8_t r1e = uzp.val[0];  // first of each pair
    uint8x8_t r1o = uzp.val[1];  // second of each pair
    uint8x8_t r1e_alive = vcgt_u8(r1e, vdup_n_u8(0));
    // Widen to uint16 so two bytes can coexist in one value
    uint16x8_t r1l = vmovl_u8(r1e);
    uint16x8_t r1r = vmovl_u8(r1o);
    // If left is alive, shift right element up by 8 bits (1 byte); else keep at
    // position 0
    int16x8_t r1s = vreinterpretq_s16_u16(
        vshlq_n_u16(vmovl_u8(vshr_n_u8(r1e_alive, 7)), 3));
    uint16x8_t r1m = vorrq_u16(r1l, vshlq_u16(r1r, r1s));
    // Count kept elements per pair (0, 1, or 2)
    uint8x8_t r1o_alive = vcgt_u8(r1o, vdup_n_u8(0));
    uint8x8_t r1cnt = vadd_u8(vshr_n_u8(r1e_alive, 7), vshr_n_u8(r1o_alive, 7));

    // --- Round 2: 8 uint16 → 4 uint32 (merge adjacent pairs of pairs) ---
    uint16x4x2_t uzp16 = vuzp_u16(vget_low_u16(r1m), vget_high_u16(r1m));
    uint32x4_t r2l = vmovl_u16(uzp16.val[0]);  // quads' left halves
    uint32x4_t r2r = vmovl_u16(uzp16.val[1]);  // quads' right halves
    // Left-pair counts: r1cnt at even positions
    uint8x8x2_t cnt_uzp = vuzp_u8(r1cnt, r1cnt);
    uint8x8_t r2lc = cnt_uzp.val[0];
    // Shift right half by (left_count * 8) bits
    uint16x4_t r2lc16 = vget_low_u16(vmovl_u8(r2lc));
    int32x4_t r2s = vreinterpretq_s32_u32(vmulq_n_u32(vmovl_u16(r2lc16), 8));
    uint32x4_t r2m = vorrq_u32(r2l, vshlq_u32(r2r, r2s));
    // Quad counts
    uint8x8_t r2rc = cnt_uzp.val[1];
    uint8x8_t r2cnt = vadd_u8(r2lc, r2rc);

    // --- Round 3: 4 uint32 → 2 uint64 (merge adjacent quads into octets) ---
    uint32x2x2_t uzp32 = vuzp_u32(vget_low_u32(r2m), vget_high_u32(r2m));
    uint64x2_t r3l = vmovl_u32(uzp32.val[0]);
    uint64x2_t r3r = vmovl_u32(uzp32.val[1]);
    int lc_a = vget_lane_u8(r2cnt, 0);  // left quad count of octet 0
    int lc_b = vget_lane_u8(r2cnt, 2);  // left quad count of octet 1
    int rc_a = vget_lane_u8(r2cnt, 1);  // right quad count of octet 0
    int rc_b = vget_lane_u8(r2cnt, 3);  // right quad count of octet 1
    int64x2_t r3s = {(int64_t)lc_a * 8, (int64_t)lc_b * 8};
    uint64x2_t r3m = vorrq_u64(r3l, vshlq_u64(r3r, r3s));

    int n_lo = lc_a + rc_a;  // packed byte count in r3m[0]
    int n_hi = lc_b + rc_b;  // packed byte count in r3m[1]
    int n_total = n_lo + n_hi;
    if (n_total == 0) return 0;

    // --- Unpack: extract values + indices from packed bytes, store via NEON
    // ---
    uint8x8_t blo = vreinterpret_u8_u64(vget_low_u64(r3m));
    uint8x8_t bhi = vreinterpret_u8_u64(vget_high_u64(r3m));

    // Values: low nibble is one-hot → decode to 2-bit via table lookup
    vst1_u8(rle_data, vqtbl1_u8(decode_tbl, vand_u8(blo, vdup_n_u8(0x0F))));
    if (n_hi > 0)
        vst1_u8(rle_data + n_lo,
                vqtbl1_u8(decode_tbl, vand_u8(bhi, vdup_n_u8(0x0F))));

    // Positions: high nibble is lane index → widen to uint32, add pos_base
    uint32x4_t cs = vdupq_n_u32(pos_base);
    uint16x8_t ilo16 = vmovl_u8(vshr_n_u8(blo, 4));
    vst1q_u32(rle_pos, vaddq_u32(vmovl_u16(vget_low_u16(ilo16)), cs));
    vst1q_u32(rle_pos + 4, vaddq_u32(vmovl_u16(vget_high_u16(ilo16)), cs));
    if (n_hi > 0) {
        uint16x8_t ihi16 = vmovl_u8(vshr_n_u8(bhi, 4));
        vst1q_u32(rle_pos + n_lo,
                  vaddq_u32(vmovl_u16(vget_low_u16(ihi16)), cs));
        vst1q_u32(rle_pos + n_lo + 4,
                  vaddq_u32(vmovl_u16(vget_high_u16(ihi16)), cs));
    }

    return n_total;
}
#endif

inline int64_t onehot_to_rle(const uint8_t *data, int64_t len,
                             uint8_t *rle_data, uint32_t *rle_pos) {
    if (len <= 0) return 0;

    static const uint8_t decode4_tbl_arr[16] = {0, 0, 1, 0, 2, 0, 0, 0,
                                                3, 0, 0, 0, 0, 0, 0, 0};

    // Decode and emit first base
    rle_data[0] = decode4_tbl_arr[data[0] >> 4];
    rle_pos[0] = 0;
    int64_t out = 1;

#if defined(__ARM_NEON) || defined(__ARM_NEON__)
    const uint8x16_t lo_mask = vdupq_n_u8(0x0F);
    const uint8x16_t decode_tbl = vld1q_u8(decode4_tbl_arr);

    // Lane index template: high nibble = lane index, used to tag each base with
    // its position
    static const uint8_t idx_arr[16] = {0x00, 0x10, 0x20, 0x30, 0x40, 0x50,
                                        0x60, 0x70, 0x80, 0x90, 0xA0, 0xB0,
                                        0xC0, 0xD0, 0xE0, 0xF0};
    const uint8x16_t idx_template = vld1q_u8(idx_arr);

    // Track last one-hot value for cross-chunk comparison.
    // One-hot format: A=1, C=2, G=4, T=8.
    uint8x16_t prev_last = vdupq_n_u8(data[0] >> 4);
    int64_t byte_idx = 0;

    // Each iteration processes 16 input bytes → 32 bases.
    for (; byte_idx + 16 <= len / 2; byte_idx += 16) {
        // --- Phase 1: Unpack BAM nibbles to one-hot bytes (no decode) ---
        uint8x16_t v = vld1q_u8(data + byte_idx);
        uint8x16_t hi = vshrq_n_u8(v, 4);      // high nibbles as bytes
        uint8x16_t lo = vandq_u8(v, lo_mask);  // low nibbles as bytes
        uint8x16x2_t z = vzipq_u8(hi, lo);
        uint8x16_t d0 = z.val[0];  // one-hot bases [0..15]
        uint8x16_t d1 = z.val[1];  // one-hot bases [16..31]

        // --- Phase 2: Build previous-base vectors for RLE comparison ---
        uint8x16_t p0 = vextq_u8(prev_last, d0, 15);
        uint8x16_t p1 = vextq_u8(d0, d1, 15);
        prev_last = d1;

        // --- Phase 3: Compare each base with its predecessor ---
        uint8x16_t neq0 = vmvnq_u8(vceqq_u8(d0, p0));
        uint8x16_t neq1 = vmvnq_u8(vceqq_u8(d1, p1));

        // Base 0 was already emitted; mark as unchanged.
        if (byte_idx == 0) neq0 = vsetq_lane_u8(0, neq0, 0);

        // Fast path: skip if all 32 bases match their predecessors.
        if (vmaxvq_u8(vorrq_u8(neq0, neq1)) == 0) continue;

        int64_t chunk_start = byte_idx * 2;

        // --- Phase 4: Binary reduction emit (no scalar ctz loop) ---
        if (vmaxvq_u8(neq0) != 0) {
            out += compact16_emit(d0, neq0, idx_template, decode_tbl,
                                  rle_data + out, rle_pos + out,
                                  static_cast<uint32_t>(chunk_start));
        }
        if (vmaxvq_u8(neq1) != 0) {
            out += compact16_emit(d1, neq1, idx_template, decode_tbl,
                                  rle_data + out, rle_pos + out,
                                  static_cast<uint32_t>(chunk_start + 16));
        }
    }

    // Scalar tail: remaining input bytes
    for (; byte_idx < len / 2; ++byte_idx) {
        uint8_t hi_base = 27 - std::countl_zero(static_cast<unsigned>(data[byte_idx]));
        uint8_t lo_base = std::countr_zero(static_cast<unsigned>(data[byte_idx]));
        int64_t hi_pos = byte_idx * 2;
        int64_t lo_pos = byte_idx * 2 + 1;
        if (hi_pos > 0) {
            if (hi_base != rle_data[out - 1]) {
                rle_data[out] = hi_base;
                rle_pos[out] = static_cast<uint32_t>(hi_pos);
                ++out;
            }
        }
        if (lo_pos < len) {
            if (lo_base != rle_data[out - 1]) {
                rle_data[out] = lo_base;
                rle_pos[out] = static_cast<uint32_t>(lo_pos);
                ++out;
            }
        }
    }
    if (len % 2 == 1) {
        uint8_t last_base = 27 - std::countl_zero(static_cast<unsigned>(data[len / 2]));
        if (last_base != rle_data[out - 1]) {
            rle_data[out] = last_base;
            rle_pos[out] = static_cast<uint32_t>(len - 1);
            ++out;
        }
    }
#else
    // Scalar fallback
    for (int64_t byte_idx = 0; byte_idx < len / 2; ++byte_idx) {
        uint8_t hi_base = 27 - std::countl_zero(static_cast<unsigned>(data[byte_idx]));
        uint8_t lo_base = std::countr_zero(static_cast<unsigned>(data[byte_idx]));
        int64_t hi_pos = byte_idx * 2;
        int64_t lo_pos = byte_idx * 2 + 1;
        if (hi_pos > 0 || byte_idx > 0) {  // skip base 0 (already emitted)
            if (hi_base != rle_data[out - 1]) {
                rle_data[out] = hi_base;
                rle_pos[out] = static_cast<uint32_t>(hi_pos);
                ++out;
            }
        }
        if (lo_pos < len) {
            if (lo_base != rle_data[out - 1]) {
                rle_data[out] = lo_base;
                rle_pos[out] = static_cast<uint32_t>(lo_pos);
                ++out;
            }
        }
    }
    if (len % 2 == 1) {
        uint8_t last_base = 27 - std::countl_zero(static_cast<unsigned>(data[len / 2]));
        if (last_base != rle_data[out - 1]) {
            rle_data[out] = last_base;
            rle_pos[out] = static_cast<uint32_t>(len - 1);
            ++out;
        }
    }
#endif
    return out;
}

inline uint64_t parse1_bam_syncmers(const uint8_t *rle_data,
                                    const uint32_t *rle_pos, int64_t rle_len,
                                    struct lps *str,
                                    const auto &syncmer_functor) {
    const auto k = syncmer_functor.k();
    if (rle_len < k) return 0;

    if (str->capacity < rle_len * CONSTANT_FACTOR) {
        str->capacity = 2 * rle_len * CONSTANT_FACTOR;
        str->cores = (struct core *)realloc(
            str->cores, (str->capacity) * sizeof(struct core));
    }

    uint64_t KMASK = ((uint64_t)1 << (2 * k)) - 1;
    uint64_t value = 0;
    for (int64_t i = 0; i < k; ++i) {
        value = ((value << 2) & KMASK) | rle_data[i];
    }

    int core_index = 0;
    auto *cores = str->cores;
    for (int64_t i = k; i <= rle_len; ++i) {
        if (syncmer_functor(value)) {
            uint32_t pos = rle_pos[i - 1] + 1;
            init_core_kmer<uint32_t>(&(cores[core_index]), value, pos, pos + k,
                                     k);
            core_index++;
        }
        if (i < rle_len) {
            value = ((value << 2) & KMASK) | rle_data[i];
        }
    }
    return core_index;
}

inline uint64_t parse1_bam_syncmers(const uint8_t *data, int64_t len,
                                    struct lps *str,
                                    const auto &syncmer_functor,
                                    charsi &rle_buf, pos_scratch &pos_buf) {
    rle_buf.require(len);
    pos_buf.require(len);
    int64_t rle_len = rle_compact(data, len, rle_buf.data, pos_buf.data);
    return parse1_bam_syncmers(rle_buf.data, pos_buf.data, rle_len, str,
                               syncmer_functor);
}
inline uint64_t parse1_CrumbledSubSeqs_syncmers(
    const gtl::vector<SequenceContainerType> &sseqs, struct lps *str,
    const auto &syncmer_functor) {
    int core_index = 0;
    if (sseqs.empty()) {
        fprintf(stderr, "Error: Empty sequence list\n");
        return core_index;
    }
    KmerViewOnCrumbs64<uint64_t> kmer_view{syncmer_functor.k()};
    size_t it1 =
        kmer_view.fill_rle(sseqs[0].data, sseqs[0].len + sseqs[0].skip_first_n,
                           sseqs[0].skip_first_n);
    size_t l_so_far = 0;
    for (size_t i = 0; i < sseqs.size(); ++i) {
        for (; it1 < sseqs[i].len + sseqs[i].skip_first_n;) {
            if (syncmer_functor(kmer_view.get())) {
                init_core_kmer<uint32_t>(&(str->cores[core_index]),
                                         kmer_view.get(),
                                         l_so_far + it1 - sseqs[i].skip_first_n,
                                         l_so_far + it1 + syncmer_functor.k() -
                                             sseqs[i].skip_first_n,
                                         syncmer_functor.k());
                core_index++;
            }
            it1 = kmer_view.move_get_rle(sseqs[i].data, it1,
                                         sseqs[i].len + sseqs[i].skip_first_n);
        }
        if (sseqs.size() > i + 1) {
            it1 = sseqs[i + 1].skip_first_n;
        }
        l_so_far += sseqs[i].len;
    }
    return core_index;
}

uint64_t parse1_bam_lmax(const uint8_t *data, int64_t len, struct core *cores);

/// @brief Initializes an lps struct with a given sequence in the bam compressed
// void init_lps_bam(lps *lps_ptr, const uint8_t *data, size_t len);
//
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
                lps_ptr->cores, (len / CONSTANT_FACTOR) * sizeof(struct core));
        }
    }
}
inline void init_lps_bam(lps *lps_ptr, const uint8_t *data, size_t len,
                         const auto &syncmer_functor, charsi &rle_buf,
                         pos_scratch &pos_buf) {
    prepare_lps(lps_ptr, len);
    lps_ptr->size = parse1_bam_syncmers(data, len, lps_ptr, syncmer_functor,
                                        rle_buf, pos_buf);
}

inline void init_lps_crumbled(lps *lps_ptr,
                              const gtl::vector<SequenceContainerType> &sseqs,
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

inline void parse_cores_from_bam_record(const uint8_t *data, size_t len,
                                        int lcp_level, lps *str,
                                        const auto &syncmer_functor,
                                        charsi &chs, lcp_read &lr,
                                        pos_scratch &pos_buf,
                                        bool reversed = false) {
    chs.require(len);
    pos_buf.require(len);
    prepare_lps(str, len);
    if (reversed) {
        onehotto2bitDNA_bivec_chs_reverse(data, len, chs);
        str->size = parse1_bam_syncmers(chs.data, len, str, syncmer_functor,
                                        chs, pos_buf);
    } else {
        int64_t rle_len = onehot_to_rle(data, len, chs.data, pos_buf.data);
        str->size = parse1_bam_syncmers(chs.data, pos_buf.data, rle_len, str,
                                        syncmer_functor);
    }
    lps_deepen(str, lcp_level);

    for (int i = 0; i < str->size; ++i) {
        lr.cores.emplace_back(str->cores[i].label);
        lr.pos.emplace_back(str->cores[i].start, str->cores[i].end);
    }
    lr.data = data;
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
