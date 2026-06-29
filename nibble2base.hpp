#ifndef NIBBLE2BASE_HPP
#define NIBBLE2BASE_HPP

#include <algorithm>
#include <cassert>
#include <cstdint>
#include <cstring>
#include <gtl/vector.hpp>
#include <iterator>


#include "htslib/hts.h"

#include <htslib/sam.h>



static const uint8_t decode4_tbl[16] = {
    /* 0 */ 0,
    /* 1 */ 0,
    /* 2 */ 1,
    /* 3 */ 0,
    /* 4 */ 2,
    /* 5 */ 0,
    /* 6 */ 0,
    /* 7 */ 0,
    /* 8 */ 3,
    /* 9 */ 0,
    /* A */ 0,
    /* B */ 0,
    /* C */ 0,
    /* D */ 0,
    /* E */ 0,
    /* F */ 0};

#if (defined(__ARM_NEON) || defined(__ARM_NEON__))

#include <arm_neon.h>
#include <stddef.h>
#include <stdint.h>

inline void onehotto2bitDNA_neon(
    const uint8_t *__restrict data, size_t len, uint8_t *__restrict low,
    uint8_t *__restrict high) {
    const uint8x16_t tbl = vld1q_u8(decode4_tbl);
    const uint8x16_t mask = vdupq_n_u8(0x0F);
    size_t i = 0;
    for (; i + 16 <= len; i += 16) {
        uint8x16_t v = vld1q_u8(data + i);
        uint8x16_t lo = vandq_u8(v, mask);
        uint8x16_t hi = vshrq_n_u8(v, 4);
        uint8x16_t dlo = vqtbl1q_u8(tbl, lo);
        uint8x16_t dhi = vqtbl1q_u8(tbl, hi);
        vst1q_u8(low + i, dlo);
        vst1q_u8(high + i, dhi);
    }
    for (; i < len; ++i) {
        uint8_t b = data[i];
        low[i] = decode4_tbl[b & 0x0F];
        high[i] = decode4_tbl[b >> 4];
    }
}

#elif defined(__SSSE3__)

#include <tmmintrin.h>
#include <stddef.h>
#include <stdint.h>

inline void onehotto2bitDNA_sse(
    const uint8_t *__restrict data, size_t len, uint8_t *__restrict low,
    uint8_t *__restrict high) {
    const __m128i tbl = _mm_loadu_si128(reinterpret_cast<const __m128i *>(decode4_tbl));
    const __m128i mask = _mm_set1_epi8(0x0F);
    size_t i = 0;
    for (; i + 16 <= len; i += 16) {
        __m128i v = _mm_loadu_si128(reinterpret_cast<const __m128i *>(data + i));
        __m128i lo = _mm_and_si128(v, mask);
        __m128i hi = _mm_and_si128(_mm_srli_epi16(v, 4), mask);
        __m128i dlo = _mm_shuffle_epi8(tbl, lo);
        __m128i dhi = _mm_shuffle_epi8(tbl, hi);
        _mm_storeu_si128(reinterpret_cast<__m128i *>(low + i), dlo);
        _mm_storeu_si128(reinterpret_cast<__m128i *>(high + i), dhi);
    }
    for (; i < len; ++i) {
        uint8_t b = data[i];
        low[i] = decode4_tbl[b & 0x0F];
        high[i] = decode4_tbl[b >> 4];
    }
}

#endif

inline void onehotto2bitDNA_vector(const uint8_t *__restrict data, size_t len,
                                   uint8_t *__restrict low,
                                   uint8_t *__restrict high) {
#if (defined(__ARM_NEON) || defined(__ARM_NEON__))
    return onehotto2bitDNA_neon(data, len, low, high);
#elif defined(__SSSE3__)
    return onehotto2bitDNA_sse(data, len, low, high);
#else
    for (size_t i = 0; i < len; ++i) {
        uint8_t b = data[i];
        low[i] = decode4_tbl[b & 0x0F];
        high[i] = decode4_tbl[b >> 4];
    }
#endif
}
// Uses naive conversions for now can be optimized later if necessary
struct CrumbledSubSeq {
    size_t len : 62;
    uint8_t skip_first_n : 2;
    uint8_t *data;
    struct from_ascii {};
    struct from_nibbles {};
    static constexpr size_t storage_size(size_t bases) noexcept {
        return bases / 4 + 4;
    }
    uint8_t *nibble2crumb(uint8_t *dest, const uint8_t *src, size_t from,
                          size_t to) {
        size_t _from = from & (~3);
        size_t _to = to & (~3);
        if (_to != to) {
            _to += 4;
        }
        for (size_t i = _from; i < _to; i += 4) {
            size_t write_idx = (i - _from) / 4;
            dest[write_idx] = seq_nt16_int[bam_seqi(src, i)] |
                              seq_nt16_int[bam_seqi(src, i + 1)] << 2 |
                              seq_nt16_int[bam_seqi(src, i + 2)] << 4 |
                              seq_nt16_int[bam_seqi(src, i + 3)] << 6;
            // dest[write_idx] = 0xff;
        }
        skip_first_n = from - _from;

        return dest;
    }
    uint8_t *chars2crumb(uint8_t *dest, const char *src, size_t from,
                         size_t to) {
        for (size_t i = from; i < to; i++) {
            size_t write_index = (i - from) / 4;
            size_t write_crumb = (i - from) % 4;
            if (write_crumb == 0) {
                dest[write_index] = 0;
            }
            switch (src[i]) {
                case 'a':
                case 'A':
                    break;
                case 'c':
                case 'C':
                    dest[write_index] |= 1 << (2 * write_crumb);
                    break;
                case 'g':
                case 'G':
                    dest[write_index] |= 2 << (2 * write_crumb);
                    break;
                case 'U':
                case 'u':
                case 'T':
                case 't':
                    dest[write_index] |= 3 << (2 * write_crumb);
                    break;
                default:
                    break;
            }
        }
        return dest;
    }
    char *crumb2base(char *dest, const uint8_t *src, size_t len) const {
        static constexpr const char *lut = "ACGT";
        for (size_t i = skip_first_n; i < len + skip_first_n; ++i) {
            dest[i - skip_first_n] = lut[(src[i / 4] >> (2 * (i % 4))) & 3];
        }
        return dest;
    }
    CrumbledSubSeq(const char *nibbles, size_t from, size_t to,
                   from_ascii) noexcept
        : len((to - from)),
          skip_first_n(0),
          data(new uint8_t[len / 4 + 4]) {  //(skip_first?1:0)]) {
        chars2crumb(data, nibbles, from, to);
    }
    CrumbledSubSeq(const uint8_t *nibbles, size_t from, size_t to,
                   from_nibbles) noexcept
        : len((to - from)),
          skip_first_n(from & 3),
          data(new uint8_t[len / 4 + 4]) {  //(skip_first?1:0)]) {
        nibble2crumb(data, nibbles, from, to);
    }
    CrumbledSubSeq(const uint8_t *nibbles, size_t from, size_t to) noexcept
        : len((to - from)),
          skip_first_n(from & 3),
          data(new uint8_t[len / 4 + 4]) {  //(skip_first?1:0)]) {
        nibble2crumb(data, nibbles, from, to);
    }
    CrumbledSubSeq(size_t length, uint8_t skip) noexcept
        : len(length),
          skip_first_n(skip),
          data(new uint8_t[storage_size(length)]()) {}
    CrumbledSubSeq() noexcept : CrumbledSubSeq(0, 0) {}
    CrumbledSubSeq(CrumbledSubSeq &&other) noexcept
        : len(other.len), skip_first_n(other.skip_first_n), data(other.data) {
        other.data = nullptr;
    }
    CrumbledSubSeq(const CrumbledSubSeq &other) = delete;
    CrumbledSubSeq &operator=(const CrumbledSubSeq &other) = delete;
    ~CrumbledSubSeq() { delete[] data; }
    CrumbledSubSeq &operator=(CrumbledSubSeq &&other) noexcept {
        if (this != &other) {
            delete[] data;
            data = other.data;
            skip_first_n = other.skip_first_n;
            len = other.len;
            other.data = nullptr;
        }
        return *this;
    }
    char *convert_to_ascii(char *buffer) const {
        crumb2base(buffer, data, len);
        buffer[len] = '\0';
        return buffer;
    }
    [[nodiscard]] CrumbledSubSeq reverse_complement() const {
        CrumbledSubSeq rc{len, 0};
        for (size_t i = 0; i < len; ++i) {
            const size_t src_pos = len - 1 - i + skip_first_n;
            const size_t src_byte = src_pos / 4;
            const size_t src_shift = (src_pos % 4) * 2;
            const uint8_t code = (data[src_byte] >> src_shift) & 0x3;
            const uint8_t complement = static_cast<uint8_t>(code ^ 0x3);
            const size_t dst_byte = i / 4;
            const size_t dst_shift = (i % 4) * 2;
            rc.data[dst_byte] |= static_cast<uint8_t>(complement << dst_shift);
        }
        return rc;
    }
};

// https://github.com/ruanjue/wtdbg2/blob/b77c5657c8095412317e4a20fe3668f5bde6b1ac/dna.h#L108
template <class T>
T reverse_complement_kmer(T kmer, int K) {
    uint32_t shift_count = (sizeof(T) * 8) - 2 * K;
    kmer = ~kmer;
    kmer = ((kmer & 0x3333333333333333LLU) << 2) |
           ((kmer & 0xCCCCCCCCCCCCCCCCLLU) >> 2);
    kmer = ((kmer & 0x0F0F0F0F0F0F0F0FLLU) << 4) |
           ((kmer & 0xF0F0F0F0F0F0F0F0LLU) >> 4);
    if constexpr (sizeof(T) == 8) {
        kmer = __builtin_bswap64(kmer);
    } else if constexpr (sizeof(T) == 4) {
        kmer = __builtin_bswap32(kmer);
    } else if constexpr (sizeof(T) == 2) {
        kmer = __builtin_bswap16(kmer);
    }
    kmer = kmer >> shift_count;
    return kmer;
}

// Small-string-optimization (SSO) toggle. When LASM_CRUMBLED_SSO is defined,
// CrumbledSubSeq64 stores its packed words inline (in a union over the data
// pointer slot) whenever they fit in LASM_CRUMBLED_SSO_WORDS words, avoiding a
// heap allocation and the dependent pointer-chase on access. Otherwise it keeps
// the original pointer-to-heap behaviour. All sequence-data access goes through
// words(); skip_first_n semantics are identical in both modes (the inline
// buffer holds the exact same bit layout a heap buffer would).

#ifdef LASM_CRUMBLED_SSO
#ifndef LASM_CRUMBLED_SSO_WORDS
#define LASM_CRUMBLED_SSO_WORDS 14  // 14*8 = 112 B inline -> 128 B struct
#endif

#endif
// 64-bit packed variant: 32 bases per word, 2 bits per base.
struct CrumbledSubSeq64 {
#ifdef LASM_CRUMBLED_SSO
    size_t len : 57;
    bool is_inline : 1;
    int skip_first_n : 6;  // offset within first word in bases (0-31)
    uint64_t origin_node;
    union {
        uint64_t *heap_;
        uint64_t inline_[LASM_CRUMBLED_SSO_WORDS];
    } store_;
#else
    size_t len : 58;
    int skip_first_n : 6;  // offset within first word in bases (0-31)
    uint64_t *data;
    uint64_t origin_node;
#endif
    struct from_ascii {};
    struct from_nibbles {};

    // Single point of access to the packed words (inline or heap).
    const uint64_t *words() const {
#ifdef LASM_CRUMBLED_SSO
        return is_inline ? store_.inline_ : store_.heap_;
#else
        return data;
#endif
    }
    uint64_t *words() {
#ifdef LASM_CRUMBLED_SSO
        return is_inline ? store_.inline_ : store_.heap_;
#else
        return data;
#endif
    }
    // Set up zeroed storage for `nwords` words (inline if it fits, else heap).
    void alloc_words_(size_t nwords) {
#ifdef LASM_CRUMBLED_SSO
        if (nwords <= LASM_CRUMBLED_SSO_WORDS) {
            is_inline = true;
            std::memset(store_.inline_, 0, sizeof(store_.inline_));
        } else {
            is_inline = false;
            store_.heap_ = new uint64_t[nwords]();
        }
#else
        data = new uint64_t[nwords]();
#endif
    }
    void free_words_() {
#ifdef LASM_CRUMBLED_SSO
        if (!is_inline) delete[] store_.heap_;
#else
        delete[] data;
#endif
    }
    // Never split (heap-backed): each logical subseq is exactly one record, so
    // a core is never continued. Mirrors CrumbledFixedSeq::continued() so the
    // per-core consumers use one code path in both modes.
    bool continued() const { return false; }
    void set_continued(bool) {}  // no-op: never split

    // Uniform append API (single-record here; CrumbledFixedSeq splits).
    template <class Vec>
    static void append_nibbles(Vec &v, const uint8_t *nibbles, size_t from,
                               size_t to, uint64_t origin) {
        if (to > from) v.emplace_back(nibbles, from, to, origin, from_nibbles{});
    }
    template <class Vec>
    static void append_ascii(Vec &v, const char *bases, size_t from, size_t to,
                             uint64_t origin) {
        if (to > from) v.emplace_back(bases, from, to, origin, from_ascii{});
    }
    template <class Vec>
    static void append_subseqs(Vec &v, const gtl::vector<CrumbledSubSeq64> &s,
                               size_t from, size_t to) {
        v.emplace_back(s, from, to);
    }
    template <class Vec>
    static void append_subseqs(Vec &v, const gtl::vector<CrumbledSubSeq64> &s,
                               size_t start_idx, size_t from, size_t to) {
        v.emplace_back(s, start_idx, from, to);
    }
    template <class Vec, class It>
    static void append_concat(Vec &v, It first, It last) {
        v.emplace_back(first, last);
    }

    static constexpr size_t bases_per_word = 32;
    static constexpr size_t storage_size(size_t bases) noexcept {
        return bases / bases_per_word + 1;  // extra word for padding
    }
    template <std::forward_iterator It>
    static size_t combined_length(It first, It last) {
        size_t total = 0;
        for (auto it = first; it != last; ++it) {
            total += static_cast<size_t>((*it).len);
        }
        return total;
    }

    static inline uint8_t base_to_code(char c) {
        switch (c) {
            case 'a':
            case 'A':
                return 0;
            case 'c':
            case 'C':
                return 1;
            case 'g':
            case 'G':
                return 2;
            case 'u':
            case 'U':
            case 't':
            case 'T':
                return 3;
            default:
                return 0;
        }
    }

    template <class FetchCode>
    void pack_range(uint64_t *dest, size_t aligned_from, size_t to,
                    FetchCode &&fetch) {
        const size_t total_bases = to - aligned_from;
        for (size_t i = 0; i < total_bases; ++i) {
            const uint64_t code = fetch(aligned_from + i) & 0x3;
            const size_t bit_index = i * 2;
            dest[bit_index >> 6] |= code << (bit_index & 63);
        }
    }

    CrumbledSubSeq64(const uint8_t *nibbles, size_t from, size_t to,
                     uint64_t origin_node, from_nibbles) noexcept
        : len(to - from),
#ifdef LASM_CRUMBLED_SSO
          is_inline(true),
#endif
          skip_first_n(static_cast<uint8_t>(from & 31)),
          origin_node(origin_node) {
        const size_t aligned_from = from & ~static_cast<size_t>(31);
        const size_t total_bases = (to - aligned_from);
        alloc_words_(storage_size(total_bases));
        pack_range(words(), aligned_from, to, [&](size_t idx) {
            return static_cast<uint64_t>(seq_nt16_int[bam_seqi(nibbles, idx)]);
        });
    }

    CrumbledSubSeq64(const uint8_t *nibbles, size_t from, size_t to,
                     uint64_t origin_node) noexcept
        : CrumbledSubSeq64(nibbles, from, to, origin_node, from_nibbles{}) {}

    CrumbledSubSeq64(const char *bases, size_t from, size_t to,
                     uint64_t origin_node, from_ascii) noexcept
        : len(to - from),
#ifdef LASM_CRUMBLED_SSO
          is_inline(true),
#endif
          skip_first_n(static_cast<uint8_t>(from & 31)),
          origin_node(origin_node) {
        const size_t aligned_from = from & ~static_cast<size_t>(31);
        const size_t total_bases = (to - aligned_from);
        alloc_words_(storage_size(total_bases));
        pack_range(words(), aligned_from, to, [&](size_t idx) {
            return static_cast<uint64_t>(base_to_code(bases[idx]));
        });
    }

    CrumbledSubSeq64(size_t length, uint8_t skip_bases) noexcept
        : len(length),
#ifdef LASM_CRUMBLED_SSO
          is_inline(true),
#endif
          skip_first_n(static_cast<uint8_t>(skip_bases & 63)),
          origin_node(0) {
        alloc_words_(storage_size(length + skip_bases));
    }
    CrumbledSubSeq64() noexcept : CrumbledSubSeq64(0, 0) {}
    CrumbledSubSeq64(CrumbledSubSeq64 &&other) noexcept
        : len(other.len),
#ifdef LASM_CRUMBLED_SSO
          is_inline(other.is_inline),
#endif
          skip_first_n(other.skip_first_n),
          origin_node(other.origin_node) {
#ifdef LASM_CRUMBLED_SSO
        store_ = other.store_;       // bitwise: copies inline data or pointer
        other.is_inline = true;      // neutralize: source dtor frees nothing
        other.len = 0;
#else
        data = other.data;
        other.data = nullptr;
#endif
        other.origin_node = 0;
    }

    template <std::forward_iterator It>
    CrumbledSubSeq64(It first, It last)
        : len(combined_length(first, last)),
#ifdef LASM_CRUMBLED_SSO
          is_inline(true),
#endif
          skip_first_n(0),
          origin_node(0) {
        alloc_words_(storage_size(len + 32));
        size_t dest_pos = 0;
        for (auto it = first; it != last; ++it) {
            const auto &frag = *it;
            for (size_t b = 0; b < frag.len; ++b) {
                const size_t src_bit =
                    (static_cast<size_t>(frag.skip_first_n) + b) * 2;
                const uint64_t code =
                    (frag.words()[src_bit >> 6] >> (src_bit & 63)) & 0x3;
                const size_t dst_bit = dest_pos * 2;
                words()[dst_bit >> 6] |= code << (dst_bit & 63);
                ++dest_pos;
            }
            origin_node = frag.origin_node;
        }
    }

    CrumbledSubSeq64(const gtl::vector<CrumbledSubSeq64> &subseqs, size_t start_idx, size_t from,
                     size_t to)
        : len{to - from},
#ifdef LASM_CRUMBLED_SSO
          is_inline{true},
#endif
          skip_first_n{0}, origin_node{0} {
        alloc_words_(storage_size(to - from + 32));
        size_t sub_len = (from < to) ? (to - from) : 0;
        if (sub_len == 0 || subseqs.empty()) {
            return;
        }

        size_t cumulative = 0;
        size_t dest_pos = 0;

        for (size_t i = start_idx; i < subseqs.size() && dest_pos < sub_len; ++i) {
            const auto &frag = subseqs[i];
            size_t frag_start = cumulative;
            size_t frag_end = cumulative + frag.len;

            // Compute the overlap of [frag_start, frag_end) with [from, to)
            size_t overlap_start = std::max(frag_start, from);
            size_t overlap_end = std::min(frag_end, to);

            if (overlap_start < overlap_end) {
                // Offset within this fragment
                size_t local_from = overlap_start - frag_start;
                size_t copy_len = overlap_end - overlap_start;

                // Copy bases from frag to result
                for (size_t b = 0; b < copy_len; ++b) {
                    size_t src_bit = (frag.skip_first_n + local_from + b) * 2;
                    uint64_t code =
                        (frag.words()[src_bit >> 6] >> (src_bit & 63)) & 0x3;
                    size_t dst_bit = dest_pos * 2;
                    words()[dst_bit >> 6] |= code << (dst_bit & 63);
                    ++dest_pos;
                }
                origin_node = frag.origin_node;
            }
            cumulative += frag.len;
        }
    };
    CrumbledSubSeq64(const gtl::vector<CrumbledSubSeq64> &subseqs, size_t from,
                     size_t to)
        : len{to - from},
#ifdef LASM_CRUMBLED_SSO
          is_inline{true},
#endif
          skip_first_n{0}, origin_node{0} {
        alloc_words_(storage_size(to - from + 32));
        size_t sub_len = (from < to) ? (to - from) : 0;
        if (sub_len == 0 || subseqs.empty()) {
            return;
        }

        size_t cumulative = 0;
        size_t dest_pos = 0;

        for (size_t i = 0; i < subseqs.size() && dest_pos < sub_len; ++i) {
            const auto &frag = subseqs[i];
            size_t frag_start = cumulative;
            size_t frag_end = cumulative + frag.len;

            // Compute the overlap of [frag_start, frag_end) with [from, to)
            size_t overlap_start = std::max(frag_start, from);
            size_t overlap_end = std::min(frag_end, to);

            if (overlap_start < overlap_end) {
                // Offset within this fragment
                size_t local_from = overlap_start - frag_start;
                size_t copy_len = overlap_end - overlap_start;

                // Copy bases from frag to result
                for (size_t b = 0; b < copy_len; ++b) {
                    size_t src_bit = (frag.skip_first_n + local_from + b) * 2;
                    uint64_t code =
                        (frag.words()[src_bit >> 6] >> (src_bit & 63)) & 0x3;
                    size_t dst_bit = dest_pos * 2;
                    words()[dst_bit >> 6] |= code << (dst_bit & 63);
                    ++dest_pos;
                }
                origin_node = frag.origin_node;
            }
            cumulative += frag.len;
        }
    };
    static CrumbledSubSeq64 copy_of(const CrumbledSubSeq64 &other) {
        CrumbledSubSeq64 r(other.len,
                           static_cast<uint8_t>(other.skip_first_n));
        const size_t words = storage_size(other.len + other.skip_first_n);
        std::memcpy(r.words(), other.words(), words * sizeof(uint64_t));
        r.origin_node = other.origin_node;
        return r;
    }
    CrumbledSubSeq64(const CrumbledSubSeq64 &other) = delete;
    CrumbledSubSeq64 &operator=(const CrumbledSubSeq64 &other) = delete;
    ~CrumbledSubSeq64() { free_words_(); }
    CrumbledSubSeq64 &operator=(CrumbledSubSeq64 &&other) noexcept {
        if (this != &other) {
            free_words_();
            len = other.len;
            skip_first_n = other.skip_first_n;
            origin_node = other.origin_node;
#ifdef LASM_CRUMBLED_SSO
            is_inline = other.is_inline;
            store_ = other.store_;
            other.is_inline = true;
            other.len = 0;
#else
            data = other.data;
            other.data = nullptr;
#endif
            other.origin_node = 0;
        }
        return *this;
    }

    char *convert_to_ascii(char *buffer) const {
        static constexpr const char *lut = "ACGT";
        size_t base_index = 0;

        const uint64_t *w = words();
        auto scalar_emit = [&](size_t count) {
            const size_t bit_base = static_cast<size_t>(skip_first_n) * 2;
            for (size_t i = 0; i < count; ++i) {
                const size_t bit_index = bit_base + (base_index + i) * 2;
                const uint64_t word = w[bit_index >> 6];
                buffer[base_index + i] = lut[(word >> (bit_index & 63)) & 0x3];
            }
            base_index += count;
        };

        scalar_emit(len);
        buffer[len] = '\0';
        return buffer;
    }

    [[nodiscard]] CrumbledSubSeq64 reverse_complement() const {
        CrumbledSubSeq64 rc{len, 0};
        const uint64_t *src = words();
        uint64_t *dst = rc.words();
        for (size_t i = 0; i < len; ++i) {
            const size_t src_bit = (skip_first_n + (len - 1 - i)) * 2;
            const uint64_t code = (src[src_bit >> 6] >> (src_bit & 63)) & 0x3;
            const uint64_t complement = code ^ 0x3;
            const size_t dst_bit = i * 2;
            dst[dst_bit >> 6] |= complement << (dst_bit & 63);
        }
        return rc;
    }
    [[nodiscard]] CrumbledSubSeq64 reverse_complement_old() const {
        assert(skip_first_n >= 0);
        CrumbledSubSeq64 rc{len, 31};
        const uint64_t *src = words();
        uint64_t *dst = rc.words();
        if (len < 32) {
            // Fix the skip_first
            uint64_t word = (src[0] >> skip_first_n * 2);
            if (storage_size(len + skip_first_n) > 1) {
                word |= (src[1] << (64 - skip_first_n * 2));
            }

            dst[0] = reverse_complement_kmer(word, len);
            rc.skip_first_n = 0;
            return rc;
        }
        size_t num_blocks = storage_size(len + skip_first_n);
        for (size_t i = 0; i < num_blocks - 1; ++i) {
            dst[num_blocks - i - 1] = reverse_complement_kmer(src[i], 32);
        }
        if ((len + skip_first_n) % 32 != 0)
            dst[0] = reverse_complement_kmer(src[num_blocks - 1],
                                             (len + skip_first_n) % 32);
        rc.skip_first_n = (num_blocks * 32 - len) - skip_first_n;
        return rc;
    }
};

#ifdef LASM_CRUMBLED_FIXED
#ifndef LASM_CRUMBLED_FIXED_WORDS
#define LASM_CRUMBLED_FIXED_WORDS 8
#endif

// Fixed-capacity, never-heap variant. Each record stores up to N*32 bases
// inline (2 bits/base). A logical subsequence longer than the capacity is
// stored as several consecutive records in the same vector; since every
// CrumbledSubSeq64 vector is consumed as a concatenation, this is transparent.
// The struct is trivially copyable: no pointer, no heap, no destructor work.
// skip_first_n is always 0 (records pack from bit 0); the field is retained so
// the concatenation/k-mer consumers that read it keep working unchanged.
template <int N>
struct CrumbledFixedSeq {
    static constexpr size_t bases_per_word = 32;
    static constexpr size_t capacity_bases = static_cast<size_t>(N) * bases_per_word;

    size_t len : 57;       // bases stored in THIS record (0 .. capacity_bases)
    bool cont_ : 1;        // true => the NEXT record continues the same core
    int skip_first_n : 6;  // always 0
    uint64_t origin_node;
    uint64_t words_[N];

    struct from_ascii {};
    struct from_nibbles {};

    const uint64_t *words() const { return words_; }
    uint64_t *words() { return words_; }
    // True when the following record in the vector is part of the same logical
    // core/subseq (set by the splitting append_* helpers). Concatenation
    // consumers ignore it; per-core consumers use it to regroup split cores.
    bool continued() const { return cont_; }
    void set_continued(bool b) { cont_ = b; }

    static constexpr size_t storage_size(size_t bases) noexcept {
        return bases / bases_per_word + 1;
    }
    static inline uint8_t base_to_code(char c) {
        switch (c) {
            case 'c':
            case 'C':
                return 1;
            case 'g':
            case 'G':
                return 2;
            case 'u':
            case 'U':
            case 't':
            case 'T':
                return 3;
            default:
                return 0;
        }
    }
    template <std::forward_iterator It>
    static size_t combined_length(It first, It last) {
        size_t total = 0;
        for (auto it = first; it != last; ++it)
            total += static_cast<size_t>((*it).len);
        return total;
    }

    CrumbledFixedSeq() noexcept
        : len(0), cont_(false), skip_first_n(0), origin_node(0) {
        std::memset(words_, 0, sizeof(words_));
    }
    // Zeroed record of `length` bases (length <= capacity_bases).
    CrumbledFixedSeq(size_t length, uint8_t /*skip*/) noexcept
        : len(length), cont_(false), skip_first_n(0), origin_node(0) {
        assert(length <= capacity_bases);
        std::memset(words_, 0, sizeof(words_));
    }
    // Single-record packers. The range MUST fit (<= capacity_bases); callers
    // that may exceed it use the append_* helpers below, which split.
    CrumbledFixedSeq(const uint8_t *nibbles, size_t from, size_t to,
                     uint64_t origin, from_nibbles) noexcept
        : len(to - from), cont_(false), skip_first_n(0), origin_node(origin) {
        assert(len <= capacity_bases);
        std::memset(words_, 0, sizeof(words_));
        for (size_t i = 0; i < len; ++i) {
            const uint64_t code =
                static_cast<uint64_t>(seq_nt16_int[bam_seqi(nibbles, from + i)]) & 0x3;
            const size_t bit = i * 2;
            words_[bit >> 6] |= code << (bit & 63);
        }
    }
    CrumbledFixedSeq(const uint8_t *nibbles, size_t from, size_t to,
                     uint64_t origin) noexcept
        : CrumbledFixedSeq(nibbles, from, to, origin, from_nibbles{}) {}
    CrumbledFixedSeq(const char *bases, size_t from, size_t to, uint64_t origin,
                     from_ascii) noexcept
        : len(to - from), cont_(false), skip_first_n(0), origin_node(origin) {
        assert(len <= capacity_bases);
        std::memset(words_, 0, sizeof(words_));
        for (size_t i = 0; i < len; ++i) {
            const uint64_t code = base_to_code(bases[from + i]) & 0x3;
            const size_t bit = i * 2;
            words_[bit >> 6] |= code << (bit & 63);
        }
    }
    // Iterator-concat into a SINGLE record (combined length must fit). Used by
    // tests / callers that want a single object; long merges use append_concat.
    template <std::forward_iterator It>
    CrumbledFixedSeq(It first, It last)
        : len(combined_length(first, last)),
          cont_(false),
          skip_first_n(0),
          origin_node(0) {
        assert(len <= capacity_bases);
        std::memset(words_, 0, sizeof(words_));
        size_t dest = 0;
        for (auto it = first; it != last; ++it) {
            const auto &frag = *it;
            for (size_t b = 0; b < frag.len; ++b) {
                const size_t sbit = (frag.skip_first_n + b) * 2;
                const uint64_t code = (frag.words()[sbit >> 6] >> (sbit & 63)) & 0x3;
                const size_t dbit = dest * 2;
                words_[dbit >> 6] |= code << (dbit & 63);
                ++dest;
            }
            origin_node = frag.origin_node;
        }
    }

    char *convert_to_ascii(char *buffer) const {
        static constexpr const char *lut = "ACGT";
        for (size_t i = 0; i < len; ++i) {
            const size_t bit = i * 2;
            buffer[i] = lut[(words_[bit >> 6] >> (bit & 63)) & 0x3];
        }
        buffer[len] = '\0';
        return buffer;
    }
    // Length-preserving: input is <= capacity_bases, so output is one record.
    [[nodiscard]] CrumbledFixedSeq reverse_complement() const {
        CrumbledFixedSeq rc(len, 0);
        for (size_t i = 0; i < len; ++i) {
            const size_t sbit = (len - 1 - i) * 2;
            const uint64_t code = (words_[sbit >> 6] >> (sbit & 63)) & 0x3;
            const size_t dbit = i * 2;
            rc.words_[dbit >> 6] |= (code ^ 0x3) << (dbit & 63);
        }
        return rc;  // origin_node intentionally left 0 (matches old behaviour)
    }
    static CrumbledFixedSeq copy_of(const CrumbledFixedSeq &other) {
        return other;  // trivially copyable
    }

    // ---- splitting append helpers -------------------------------------
    // Universal primitive: append ceil(total_len/capacity_bases) records,
    // pulling each logical base code from fetch(global_index).
    template <class Vec, class Fetch>
    static void append_fetch(Vec &v, size_t total_len, uint64_t origin,
                             Fetch &&fetch) {
        size_t done = 0;
        while (done < total_len) {
            const size_t chunk = std::min(capacity_bases, total_len - done);
            v.emplace_back(chunk, static_cast<uint8_t>(0));
            CrumbledFixedSeq &rec = v.back();
            rec.origin_node = origin;
            for (size_t b = 0; b < chunk; ++b) {
                const uint64_t code = fetch(done + b) & 0x3;
                const size_t bit = b * 2;
                rec.words_[bit >> 6] |= code << (bit & 63);
            }
            done += chunk;
            rec.cont_ = (done < total_len);  // more records follow => continued
        }
    }
    template <class Vec>
    static void append_nibbles(Vec &v, const uint8_t *nibbles, size_t from,
                               size_t to, uint64_t origin) {
        if (to <= from) return;
        append_fetch(v, to - from, origin, [&](size_t i) -> uint64_t {
            return seq_nt16_int[bam_seqi(nibbles, from + i)];
        });
    }
    template <class Vec>
    static void append_ascii(Vec &v, const char *bases, size_t from, size_t to,
                             uint64_t origin) {
        if (to <= from) return;
        append_fetch(v, to - from, origin, [&](size_t i) -> uint64_t {
            return base_to_code(bases[from + i]);
        });
    }
    // Extract concatenated range [from,to) of `subseqs` into split records.
    template <class Vec>
    static void append_subseqs(Vec &v, const gtl::vector<CrumbledFixedSeq> &subseqs,
                               size_t start_idx, size_t from, size_t to) {
        const size_t sub_len = (from < to) ? (to - from) : 0;
        if (sub_len == 0 || subseqs.empty()) return;
        // Pre-create zeroed records so any range not covered by a fragment
        // stays 0 (matches the old buffer-then-scatter-fill behaviour).
        const size_t first_rec = v.size();
        for (size_t done = 0; done < sub_len; done += capacity_bases)
            v.emplace_back(std::min(capacity_bases, sub_len - done),
                           static_cast<uint8_t>(0));
        // Mark all but the last record of THIS subseq as continued.
        for (size_t r = first_rec; r + 1 < v.size(); ++r) v[r].cont_ = true;
        size_t cumulative = 0, dest_pos = 0;
        for (size_t i = start_idx; i < subseqs.size() && dest_pos < sub_len; ++i) {
            const auto &frag = subseqs[i];
            const size_t fs = cumulative, fe = cumulative + frag.len;
            const size_t os = std::max(fs, from), oe = std::min(fe, to);
            if (os < oe) {
                const size_t local_from = os - fs, clen = oe - os;
                for (size_t b = 0; b < clen; ++b) {
                    const size_t sbit = (frag.skip_first_n + local_from + b) * 2;
                    const uint64_t code =
                        (frag.words()[sbit >> 6] >> (sbit & 63)) & 0x3;
                    CrumbledFixedSeq &rec = v[first_rec + dest_pos / capacity_bases];
                    const size_t off = (dest_pos % capacity_bases) * 2;
                    rec.words_[off >> 6] |= code << (off & 63);
                    rec.origin_node = frag.origin_node;
                    ++dest_pos;
                }
            }
            cumulative += frag.len;
        }
    }
    template <class Vec>
    static void append_subseqs(Vec &v, const gtl::vector<CrumbledFixedSeq> &subseqs,
                               size_t from, size_t to) {
        append_subseqs(v, subseqs, 0, from, to);
    }
    // Append the concatenation of [first,last) as one continued run of records.
    template <class Vec, std::forward_iterator It>
    static void append_concat(Vec &v, It first, It last) {
        const size_t total = combined_length(first, last);
        if (total == 0) return;
        uint64_t origin = 0;
        for (auto it = first; it != last; ++it) origin = (*it).origin_node;
        // Single fetch over the concatenation via a running fragment cursor.
        It cur = first;
        size_t base_in_frag = 0;
        append_fetch(v, total, origin, [&](size_t) -> uint64_t {
            while (cur != last && base_in_frag >= (*cur).len) {
                ++cur;
                base_in_frag = 0;
            }
            const auto &frag = *cur;
            const size_t sbit = (frag.skip_first_n + base_in_frag) * 2;
            ++base_in_frag;
            return (frag.words()[sbit >> 6] >> (sbit & 63)) & 0x3;
        });
    }
};
#endif  // LASM_CRUMBLED_FIXED

struct seqveciterator {
    size_t idx = 0;
    size_t bidx = 0;
    size_t bases_before = 0;
    void move_only(const gtl::vector<CrumbledSubSeq64> &sseqs,
                   uint64_t until) {
        size_t new_idx = idx;
        size_t new_bidx = bidx;
        while (new_idx < sseqs.size() &&
               bases_before + (sseqs[new_idx].len - new_bidx) < until) {
            bases_before += sseqs[new_idx].len - new_bidx;
            new_bidx = 0;
            new_idx++;
        }
        new_bidx += until - bases_before;
        idx = new_idx;
        bidx = new_bidx;
        bases_before = until;
    }
    CrumbledSubSeq64 get_and_move(
        const gtl::vector<CrumbledSubSeq64> &sseqs, uint64_t until) {
        // Save bases_before before the loop mutates it, so we can
        // compute the full extraction range for the constructor.
        const size_t saved_bases_before = bases_before;
        size_t new_idx = idx;
        size_t new_bidx = bidx;
        while (new_idx < sseqs.size() &&
               bases_before + (sseqs[new_idx].len - new_bidx) < until) {
            bases_before += sseqs[new_idx].len - new_bidx;
            new_bidx = 0;
            new_idx++;
        }
        new_bidx += until - bases_before;
        // 'from' and 'to' are positions relative to the concatenation
        // starting at start_idx (idx).  We start at bidx within the
        // first fragment and need (until - saved_bases_before) bases.
        CrumbledSubSeq64 seq{sseqs, idx, bidx,
                                  bidx + (until - saved_bases_before)};
        idx = new_idx;
        bidx = new_bidx;
        bases_before = until;
        return seq;
    }
};

#endif
