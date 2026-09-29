#pragma once
#include <cstdint>
#if defined(__SSSE3__)
#include <tmmintrin.h>
#endif

#include "nibble2base.hpp"
#ifndef LCMER_H
#define LCMER_H
#include <htslib/sam.h>
#include <stdint.h>

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
#include "syncmer.h"
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
            data = (decltype(data))aligned_alloc(32, target_cap);
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
    const uint8_t *data = nullptr;
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

// Cores of one LCP level (lcptools lps_soa), owning its arrays. Level-1
// cores are the selected k-mers: label = k-mer code, span = start | end << 32
// in base coordinates.
struct core_buf : lps_soa {
    core_buf() { lps_soa_init(this); }
    ~core_buf() { lps_soa_free(this); }
    core_buf(const core_buf &) = delete;
    core_buf &operator=(const core_buf &) = delete;
    uint32_t start(int64_t i) const { return static_cast<uint32_t>(span[i]); }
    uint32_t end(int64_t i) const {
        return static_cast<uint32_t>(span[i] >> 32);
    }
    // Room for n cores; grows by at least 1.5x so per-read calls stay rare.
    void reserve(int64_t n) {
        if (n <= capacity) return;
        if (lps_soa_reserve(this, std::max(n, capacity + capacity / 2)) != 0) {
            fprintf(stderr, "core_buf: out of memory\n");
            abort();
        }
    }
    // Deepens level-1 k-mer cores of `syncmer_functor` to lcp_level.
    void deepen(const auto &syncmer_functor, int lcp_level) {
        lps_soa_deepen(this, 2 * static_cast<int>(syncmer_functor.k()),
                       lcp_level);
    }
    void append_to(lcp_read &lr) const {
        lr.cores.insert(lr.cores.end(), label, label + size);
        for (int64_t i = 0; i < size; ++i)
            lr.pos.emplace_back(start(i), end(i));
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

template <>
struct std::hash<litecmer> {
    std::size_t operator()(const litecmer &l) const noexcept {
        return l.value();
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
        // k-1 more insertions (k total incl. the initial one), so the first
        // window checked by the caller is runs[0..k-1] — the same window
        // parse1_rle_syncmers checks first.
        for (i = 0; i < k - 1 && it < cap; i++) {
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

uint8_t *onehotto2bitDNA(const uint8_t *data, size_t len);
// uint8_t *onehotto2bitDNA_bivec(const uint8_t *data, size_t len);

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
#elif defined(__SSSE3__)
    int64_t i = 1;
    for (; i + 16 <= len; i += 16) {
        __m128i curr =
            _mm_loadu_si128(reinterpret_cast<const __m128i *>(data + i));
        __m128i prev =
            _mm_loadu_si128(reinterpret_cast<const __m128i *>(data + i - 1));
        int mask = ~_mm_movemask_epi8(_mm_cmpeq_epi8(curr, prev)) & 0xFFFF;
        if (mask == 0) continue;
        while (mask) {
            int idx = __builtin_ctz(mask);
            rle_data[out] = data[i + idx];
            rle_pos[out] = static_cast<uint32_t>(i + idx);
            ++out;
            mask &= mask - 1;
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

// Fused nt16 one-hot decode + RLE compaction in a single pass.
// Input: packed nt16 nibbles (2 bases per byte, as bam_get_seq() or the FASTX
// reader produce), length in bases.
// Output: rle_data/rle_pos filled with RLE-compacted 2-bit bases. Runs are
// split on the one-hot nibble and decoded through kOnehotDecode (anything
// but A/C/G/T, e.g. N, decodes as A) identically on every path.
// The SIMD paths store whole vectors past the last written entry: callers
// must provide len + kOnehotRlePad writable entries in rle_data and rle_pos.
inline constexpr int64_t kOnehotRlePad = 32;
inline constexpr uint8_t kOnehotDecode[16] = {0, 0, 1, 0, 2, 0, 0, 0,
                                              3, 0, 0, 0, 0, 0, 0, 0};

#if defined(__ARM_NEON) || defined(__ARM_NEON__)
// Binary reduction: compact 16 one-hot bytes into contiguous packed output.
// Each byte is 0 (unchanged) or (lane_idx << 4) | one_hot_value, one-hot
// {1,2,4,8} in the low nibble, lane 0-15 in the high nibble. Merge adjacent
// pairs in 3 rounds (16→8→4→2) via widen + variable-shift + OR, yielding two
// uint64 of up to 8 packed bytes each; unpack index (>>4) and value (&0xF),
// decode one-hot→2-bit, store via NEON.
static inline int compact16_emit(
    uint8x16_t data_vec,             // 16 one-hot bytes (nt16 nibble values)
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
#elif defined(__SSSE3__)
// Left-pack table: for an 8-bit lane mask m, idx[m] lists the set lanes in
// ascending order, padded with 0x80 (pshufb writes 0 for those lanes).
struct onehot_left_pack_table {
    uint8_t idx[256][8];
    constexpr onehot_left_pack_table() : idx{} {
        for (int m = 0; m < 256; ++m) {
            int n = 0;
            for (int lane = 0; lane < 8; ++lane) {
                if (m & (1 << lane)) idx[m][n++] = static_cast<uint8_t>(lane);
            }
            while (n < 8) idx[m][n++] = 0x80;
        }
    }
};
inline constexpr onehot_left_pack_table kOnehotLeftPack{};

// Emits the run starts among 16 bases without a per-run loop. d: one-hot
// bytes; mask bit i set = lane i starts a run. Each 8-lane half is packed
// by one pshufb whose control bytes are also the lanes' indices, so the
// same control, widened, gives the positions. Stores are unconditional
// (up to 16 entries past the returned count).
static inline int onehot_emit16(__m128i d, unsigned mask,
                                const __m128i &decode_tbl, uint8_t *rle_data,
                                uint32_t *rle_pos, uint32_t pos_base) {
    const unsigned m_lo = mask & 0xFF, m_hi = mask >> 8;
    const int c_lo = std::popcount(m_lo);
    uint64_t lo_idx, hi_idx;
    std::memcpy(&lo_idx, kOnehotLeftPack.idx[m_lo], 8);
    std::memcpy(&hi_idx, kOnehotLeftPack.idx[m_hi], 8);
    // High half lanes are 8..15; 0x80 padding stays >= 0x80 after the add.
    const __m128i ctrl =
        _mm_set_epi64x(static_cast<long long>(hi_idx + 0x0808080808080808ULL),
                       static_cast<long long>(lo_idx));
    const __m128i codes =
        _mm_shuffle_epi8(decode_tbl, _mm_and_si128(d, _mm_set1_epi8(0x0F)));
    const __m128i packed = _mm_shuffle_epi8(codes, ctrl);
    _mm_storel_epi64(reinterpret_cast<__m128i *>(rle_data), packed);
    _mm_storel_epi64(reinterpret_cast<__m128i *>(rle_data + c_lo),
                     _mm_srli_si128(packed, 8));
#if defined(__AVX2__)
    const __m256i base = _mm256_set1_epi32(static_cast<int>(pos_base));
    _mm256_storeu_si256(reinterpret_cast<__m256i *>(rle_pos),
                        _mm256_add_epi32(_mm256_cvtepu8_epi32(ctrl), base));
    _mm256_storeu_si256(
        reinterpret_cast<__m256i *>(rle_pos + c_lo),
        _mm256_add_epi32(_mm256_cvtepu8_epi32(_mm_srli_si128(ctrl, 8)), base));
#else
    const __m128i zero = _mm_setzero_si128();
    const __m128i base = _mm_set1_epi32(static_cast<int>(pos_base));
    const __m128i w_lo = _mm_unpacklo_epi8(ctrl, zero);  // 8 x u16
    const __m128i w_hi = _mm_unpackhi_epi8(ctrl, zero);
    auto store4 = [&](uint32_t *dst, __m128i v) {
        _mm_storeu_si128(reinterpret_cast<__m128i *>(dst),
                         _mm_add_epi32(v, base));
    };
    store4(rle_pos, _mm_unpacklo_epi16(w_lo, zero));
    store4(rle_pos + 4, _mm_unpackhi_epi16(w_lo, zero));
    store4(rle_pos + c_lo, _mm_unpacklo_epi16(w_hi, zero));
    store4(rle_pos + c_lo + 4, _mm_unpackhi_epi16(w_hi, zero));
#endif
    return c_lo + std::popcount(m_hi);
}
#endif

// Scalar RLE of bases [from, len): `prev` is the one-hot nibble of base
// from - 1. Shared by every path for the tail and as the portable fallback.
static inline int64_t onehot_to_rle_scalar(const uint8_t *data, int64_t len,
                                           int64_t from, uint8_t prev,
                                           uint8_t *rle_data, uint32_t *rle_pos,
                                           int64_t out) {
    for (int64_t i = from; i < len; ++i) {
        const uint8_t nib =
            (i & 1) ? (data[i >> 1] & 0x0F) : (data[i >> 1] >> 4);
        if (nib != prev) {
            rle_data[out] = kOnehotDecode[nib];
            rle_pos[out] = static_cast<uint32_t>(i);
            ++out;
        }
        prev = nib;
    }
    return out;
}

inline int64_t onehot_to_rle(const uint8_t *data, int64_t len,
                             uint8_t *rle_data, uint32_t *rle_pos) {
    if (len <= 0) return 0;

    // Base 0 always starts a run.
    rle_data[0] = kOnehotDecode[data[0] >> 4];
    rle_pos[0] = 0;
    int64_t out = 1;
    int64_t byte_idx = 0;

#if defined(__ARM_NEON) || defined(__ARM_NEON__)
    const uint8x16_t lo_mask = vdupq_n_u8(0x0F);
    const uint8x16_t decode_tbl = vld1q_u8(kOnehotDecode);

    // Lane index template: high nibble = lane index, used to tag each base with
    // its position
    static const uint8_t idx_arr[16] = {0x00, 0x10, 0x20, 0x30, 0x40, 0x50,
                                        0x60, 0x70, 0x80, 0x90, 0xA0, 0xB0,
                                        0xC0, 0xD0, 0xE0, 0xF0};
    const uint8x16_t idx_template = vld1q_u8(idx_arr);

    // Track last one-hot value for cross-chunk comparison.
    // One-hot format: A=1, C=2, G=4, T=8.
    uint8x16_t prev_last = vdupq_n_u8(data[0] >> 4);

    // Each iteration processes 16 input bytes → 32 bases.
    for (; byte_idx + 16 <= len / 2; byte_idx += 16) {
        // --- Phase 1: Unpack nt16 nibbles to one-hot bytes (no decode) ---
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
#elif defined(__SSSE3__)
    const __m128i lo_mask = _mm_set1_epi8(0x0F);
    const __m128i decode_tbl =
        _mm_loadu_si128(reinterpret_cast<const __m128i *>(kOnehotDecode));
    // Last one-hot base of the previous block, in lane 15.
    __m128i prev_last = _mm_set1_epi8(static_cast<char>(data[0] >> 4));

    // Each iteration processes 16 input bytes → 32 bases.
    for (; byte_idx + 16 <= len / 2; byte_idx += 16) {
        // Unpack nt16 nibbles to one-hot bytes, bases [0..15] and [16..31].
        const __m128i v =
            _mm_loadu_si128(reinterpret_cast<const __m128i *>(data + byte_idx));
        const __m128i hi = _mm_and_si128(_mm_srli_epi16(v, 4), lo_mask);
        const __m128i lo = _mm_and_si128(v, lo_mask);
        const __m128i d0 = _mm_unpacklo_epi8(hi, lo);
        const __m128i d1 = _mm_unpackhi_epi8(hi, lo);
        // Predecessor of each base; a lane starts a run if it differs.
        const __m128i p0 = _mm_alignr_epi8(d0, prev_last, 15);
        const __m128i p1 = _mm_alignr_epi8(d1, d0, 15);
        prev_last = d1;
        unsigned m0 = ~_mm_movemask_epi8(_mm_cmpeq_epi8(d0, p0)) & 0xFFFFu;
        const unsigned m1 =
            ~_mm_movemask_epi8(_mm_cmpeq_epi8(d1, p1)) & 0xFFFFu;
        if (byte_idx == 0) m0 &= ~1u;  // base 0 was already emitted

        const uint32_t chunk_start = static_cast<uint32_t>(byte_idx * 2);
        out += onehot_emit16(d0, m0, decode_tbl, rle_data + out, rle_pos + out,
                             chunk_start);
        out += onehot_emit16(d1, m1, decode_tbl, rle_data + out, rle_pos + out,
                             chunk_start + 16);
    }
#endif
    // Remaining bases (all of them without SIMD). The predecessor of base
    // 2 * byte_idx is the low nibble of the previous byte.
    const int64_t from = std::max<int64_t>(1, 2 * byte_idx);
    const uint8_t prev =
        (from & 1) ? (data[0] >> 4) : (data[byte_idx - 1] & 0x0F);
    return onehot_to_rle_scalar(data, len, from, prev, rle_data, rle_pos, out);
}

// Packs 64 selection bytes (each 0 or 1) into a 64-bit mask, bit i = s[i].
inline uint64_t sel_mask64(const uint8_t *s) {
#if defined(__ARM_NEON) || defined(__ARM_NEON__)
    // Shift byte i left by i % 8, then three pairwise-add levels sum each run
    // of 8 bytes into one byte: byte j of the result holds bits 8j..8j+7.
    const uint8x16x4_t v = vld1q_u8_x4(s);
    static constexpr int8_t kShift[16] = {0, 1, 2, 3, 4, 5, 6, 7,
                                          0, 1, 2, 3, 4, 5, 6, 7};
    const int8x16_t sh = vld1q_s8(kShift);
    const uint8x16_t p =
        vpaddq_u8(vshlq_u8(v.val[0], sh), vshlq_u8(v.val[1], sh));
    const uint8x16_t q =
        vpaddq_u8(vshlq_u8(v.val[2], sh), vshlq_u8(v.val[3], sh));
    const uint8x16_t r = vpaddq_u8(vpaddq_u8(p, q), vdupq_n_u8(0));
    return vgetq_lane_u64(vreinterpretq_u64_u8(r), 0);
#elif defined(__AVX2__)
    // Shifting each 16-bit lane left by 7 moves bit 0 of every byte to its
    // bit 7, which movemask collects.
    auto half = [](const uint8_t *p) {
        const __m256i v =
            _mm256_loadu_si256(reinterpret_cast<const __m256i *>(p));
        return static_cast<uint32_t>(
            _mm256_movemask_epi8(_mm256_slli_epi16(v, 7)));
    };
    return half(s) | (uint64_t{half(s + 32)} << 32);
#elif defined(__SSE2__)
    uint64_t mask = 0;
    for (int b = 0; b < 4; ++b) {
        const __m128i v =
            _mm_loadu_si128(reinterpret_cast<const __m128i *>(s + 16 * b));
        mask |= uint64_t{static_cast<uint32_t>(
                    _mm_movemask_epi8(_mm_slli_epi16(v, 7)))}
                << (16 * b);
    }
    return mask;
#else
    // 8 bytes of 0/1 → 8 bits via one multiply per word.
    uint64_t mask = 0;
    for (int b = 0; b < 8; ++b) {
        uint64_t w;
        std::memcpy(&w, s + 8 * b, sizeof w);
        mask |= ((w * 0x0102040810204080ULL) >> 56) << (8 * b);
    }
    return mask;
#endif
}

// Level-1 cores of an RLE stream into str (resets it): the selected k-mers.
inline uint64_t parse1_rle_syncmers(const uint8_t *rle_data,
                                    const uint32_t *rle_pos, int64_t rle_len,
                                    int64_t seq_len, core_buf *str,
                                    const auto &syncmer_functor) {
    str->level = 1;
    str->size = 0;
    const auto k = syncmer_functor.k();
    if (rle_len < k) return 0;

    // Selection for every window at once, then emit only the selected ones:
    // iterate set bits of a 64-window mask (sel_mask64); each k-mer code is 4
    // loads of 4-symbol codes. A core spans its RLE-space k-mer (runs
    // [e-k+1, e]) in original base coordinates: start of the first run up to
    // the start of the next run.
    thread_local syncmer::live_scratch scratch;
    constexpr int64_t kBlock = 64;
    scratch.require(static_cast<size_t>(rle_len + kBlock));
    uint8_t *sel = scratch.sel.data();
    std::memset(sel, 0, k - 1);
    std::memset(sel + rle_len, 0, kBlock);
    syncmer_functor.select_live(rle_data, rle_len, scratch);
    const uint8_t *q = scratch.codes();
    const uint32_t kmask = static_cast<uint32_t>((1ULL << (2 * k)) - 1);
    // Room for exactly the selected k-mers (a fraction of the windows).
    int64_t selected = 0;
    for (int64_t base = 0; base < rle_len; base += kBlock)
        selected += std::popcount(sel_mask64(sel + base));
    str->reserve(selected);
    int64_t core_index = 0;
    uint32_t *label = str->label;
    uint64_t *span = str->span;
    for (int64_t base = 0; base < rle_len; base += kBlock) {
        uint64_t mask = sel_mask64(sel + base);
        while (mask) {
            const int64_t e = base + std::countr_zero(mask);
            mask &= mask - 1;
            const uint32_t si = rle_pos[e - k + 1];
            const uint32_t ei = (e + 1 < rle_len)
                                    ? rle_pos[e + 1]
                                    : static_cast<uint32_t>(seq_len);
            label[core_index] = syncmer::live_code(q, e, kmask);
            span[core_index] = si | static_cast<uint64_t>(ei) << 32;
            core_index++;
        }
    }
    str->size = core_index;
    return core_index;
}

inline uint64_t parse1_2bit_syncmers(const uint8_t *data, int64_t len,
                                     core_buf *str, const auto &syncmer_functor,
                                     charsi &rle_buf, pos_scratch &pos_buf) {
    rle_buf.require(len);
    pos_buf.require(len);
    int64_t rle_len = rle_compact(data, len, rle_buf.data, pos_buf.data);
    return parse1_rle_syncmers(rle_buf.data, pos_buf.data, rle_len, len, str,
                               syncmer_functor);
}
inline uint64_t parse1_CrumbledSubSeqs_syncmers(
    const gtl::vector<SequenceContainerType> &sseqs, core_buf *str,
    const auto &syncmer_functor) {
    if (sseqs.empty()) {
        fprintf(stderr, "Error: Empty sequence list\n");
        str->level = 1;
        str->size = 0;
        return 0;
    }
    // Materialize the RLE stream of the concatenated fragments in original
    // (l_so_far) coordinates, then window it with the same code as the read
    // path. A run is NOT restarted when the next fragment begins with the
    // run's base.
    size_t total_len = 0;
    for (const auto &s : sseqs) total_len += s.len;
    charsi rle_data;
    pos_scratch rle_pos;
    rle_data.require(total_len);
    rle_pos.require(total_len);
    int64_t rle_len = 0;
    size_t l_so_far = 0;
    uint8_t prev_base = 255;
    for (const auto &s : sseqs) {
        const size_t cap = s.len + s.skip_first_n;
        for (size_t idx = s.skip_first_n; idx < cap; ++idx) {
            const uint8_t b =
                KmerViewOnCrumbs64<uint64_t>::base_at(s.data, idx);
            if (b != prev_base) {
                rle_data.data[rle_len] = b;
                rle_pos.data[rle_len] =
                    static_cast<uint32_t>(l_so_far + idx - s.skip_first_n);
                ++rle_len;
                prev_base = b;
            }
        }
        l_so_far += s.len;
    }
    return parse1_rle_syncmers(rle_data.data, rle_pos.data, rle_len,
                               static_cast<int64_t>(total_len), str,
                               syncmer_functor);
}

inline void init_lps_crumbled(core_buf *str,
                              const gtl::vector<SequenceContainerType> &sseqs,
                              const auto &syncmer_functor) {
    parse1_CrumbledSubSeqs_syncmers(sseqs, str, syncmer_functor);
}

// Complement of a BAM nibble (A<->T, C<->G, anything else N).
constexpr uint8_t nibble_complement(uint8_t n) {
    constexpr uint8_t comp[16] = {15, 8,  4,  15, 2,  15, 15, 15,
                                  1,  15, 15, 15, 15, 15, 15, 15};
    return comp[n];
}

// Reverse complement of a byte of two packed nibbles: the nibbles swapped and
// each complemented.
inline constexpr std::array<uint8_t, 256> kRevcompNibblePair = [] {
    std::array<uint8_t, 256> table{};
    for (int b = 0; b < 256; ++b) {
        table[b] = static_cast<uint8_t>((nibble_complement(b & 0xF) << 4) |
                                        nibble_complement(b >> 4));
    }
    return table;
}();

// Reverse complement of len packed nibbles into dst, a byte at a time. For an
// odd len the source's padding nibble would come first, so every output byte
// takes its two nibbles from neighbouring reversed bytes, and the last low
// nibble is 0.
inline void revcomp_nibbles(const uint8_t *src, uint32_t len,
                            gtl::vector<uint8_t> &dst) {
    const uint32_t bytes = (len + 1) / 2;
    dst.resize(bytes);
    if (len % 2 == 0) {
        for (uint32_t j = 0; j < bytes; ++j) {
            dst[j] = kRevcompNibblePair[src[bytes - 1 - j]];
        }
        return;
    }
    for (uint32_t j = 0; j + 1 < bytes; ++j) {
        dst[j] =
            static_cast<uint8_t>((kRevcompNibblePair[src[bytes - 1 - j]] << 4) |
                                 (kRevcompNibblePair[src[bytes - 2 - j]] >> 4));
    }
    dst[bytes - 1] = static_cast<uint8_t>(kRevcompNibblePair[src[0]] << 4);
}

inline void parse_cores_from_nibbles(const uint8_t *data, size_t len,
                                     int lcp_level, core_buf *str,
                                     const auto &syncmer_functor, charsi &chs,
                                     lcp_read &lr, pos_scratch &pos_buf,
                                     bool reversed = false) {
    chs.require(len + kOnehotRlePad);
    pos_buf.require(len + kOnehotRlePad);
    if (reversed) {
        onehotto2bitDNA_bivec_chs_reverse(data, len, chs);
        parse1_2bit_syncmers(chs.data, len, str, syncmer_functor, chs, pos_buf);
    } else {
        int64_t rle_len = onehot_to_rle(data, len, chs.data, pos_buf.data);
        parse1_rle_syncmers(chs.data, pos_buf.data, rle_len, len, str,
                            syncmer_functor);
    }
    str->deepen(syncmer_functor, lcp_level);
    str->append_to(lr);
    lr.data = data;
}
#endif
