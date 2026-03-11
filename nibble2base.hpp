#ifndef NIBBLE2BASE_HPP
#define NIBBLE2BASE_HPP
// #include <sys/_types/_u_char.h>
#include <cassert>
#include <concepts>
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


inline void onehotto2bitDNA_neon(  // TODO VERIFY TEST THIS
    const uint8_t *__restrict data, size_t len, uint8_t *__restrict low,
    uint8_t *__restrict high) {
    // Load decode table into a vector register
    const uint8x16_t tbl = vld1q_u8(decode4_tbl);
    const uint8x16_t mask = vdupq_n_u8(0x0F);

    size_t i = 0;

    // Main SIMD loop: 16 bytes at a time
    for (; i + 16 <= len; i += 16) {
        // Load 16 input bytes
        uint8x16_t v = vld1q_u8(data + i);

        // Extract nibbles
        uint8x16_t lo = vandq_u8(v, mask);
        uint8x16_t hi = vshrq_n_u8(v, 4);

        // Table lookup
        uint8x16_t dlo = vqtbl1q_u8(tbl, lo);
        uint8x16_t dhi = vqtbl1q_u8(tbl, hi);

        // Store results
        vst1q_u8(low + i, dlo);
        vst1q_u8(high + i, dhi);
    }

    // Scalar tail
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

// 64-bit packed variant: 32 bases per word, 2 bits per base.
struct CrumbledSubSeq64 {
    size_t len : 58;
    int skip_first_n : 6;  // offset within first word in bases (0-31)
    uint64_t *data;
    uint64_t origin_node;
    struct from_ascii {};
    struct from_nibbles {};

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
          skip_first_n(static_cast<uint8_t>(from & 31)),
          data(nullptr),
          origin_node(origin_node) {
        const size_t aligned_from = from & ~static_cast<size_t>(31);
        const size_t total_bases = (to - aligned_from);
        data = new uint64_t[storage_size(total_bases)]();
        pack_range(data, aligned_from, to, [&](size_t idx) {
            return static_cast<uint64_t>(seq_nt16_int[bam_seqi(nibbles, idx)]);
        });
    }

    CrumbledSubSeq64(const uint8_t *nibbles, size_t from, size_t to,
                     uint64_t origin_node) noexcept
        : CrumbledSubSeq64(nibbles, from, to, origin_node, from_nibbles{}) {}

    CrumbledSubSeq64(const char *bases, size_t from, size_t to,
                     uint64_t origin_node, from_ascii) noexcept
        : len(to - from),
          skip_first_n(static_cast<uint8_t>(from & 31)),
          data(nullptr),
          origin_node(origin_node) {
        const size_t aligned_from = from & ~static_cast<size_t>(31);
        const size_t total_bases = (to - aligned_from);
        data = new uint64_t[storage_size(total_bases)]();
        pack_range(data, aligned_from, to, [&](size_t idx) {
            return static_cast<uint64_t>(base_to_code(bases[idx]));
        });
    }

    CrumbledSubSeq64(size_t length, uint8_t skip_bases) noexcept
        : len(length),
          skip_first_n(static_cast<uint8_t>(skip_bases & 63)),
          data(new uint64_t[storage_size(length + skip_bases)]()),
          origin_node(0) {}
    CrumbledSubSeq64() noexcept : CrumbledSubSeq64(0, 0) {}
    CrumbledSubSeq64(CrumbledSubSeq64 &&other) noexcept
        : len(other.len),
          skip_first_n(other.skip_first_n),
          data(other.data),
          origin_node(other.origin_node) {
        other.data = nullptr;
        other.origin_node = 0;
    }

    template <std::forward_iterator It>
    CrumbledSubSeq64(It first, It last)
        : len(combined_length(first, last)),
          skip_first_n(0),
          data(new uint64_t[storage_size(len + 32)]()),
          origin_node(0) {
        size_t dest_pos = 0;
        for (auto it = first; it != last; ++it) {
            const auto &frag = *it;
            for (size_t b = 0; b < frag.len; ++b) {
                const size_t src_bit =
                    (static_cast<size_t>(frag.skip_first_n) + b) * 2;
                const uint64_t code =
                    (frag.data[src_bit >> 6] >> (src_bit & 63)) & 0x3;
                const size_t dst_bit = dest_pos * 2;
                data[dst_bit >> 6] |= code << (dst_bit & 63);
                ++dest_pos;
            }
            origin_node = frag.origin_node;
        }
    }

    CrumbledSubSeq64(const gtl::vector<CrumbledSubSeq64> &subseqs, size_t start_idx, size_t from,
                     size_t to)
        : len{to - from}, skip_first_n{0}, data{new uint64_t[storage_size(to - from + 32)]()}, origin_node{0} {
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
                        (frag.data[src_bit >> 6] >> (src_bit & 63)) & 0x3;
                    size_t dst_bit = dest_pos * 2;
                    data[dst_bit >> 6] |= code << (dst_bit & 63);
                    ++dest_pos;
                }
                origin_node = frag.origin_node;
            }
            cumulative += frag.len;
        }
    };
    CrumbledSubSeq64(const gtl::vector<CrumbledSubSeq64> &subseqs, size_t from,
                     size_t to)
        : len{to - from}, skip_first_n{0}, data{new uint64_t[storage_size(to - from + 32)]()}, origin_node{0} {
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
                        (frag.data[src_bit >> 6] >> (src_bit & 63)) & 0x3;
                    size_t dst_bit = dest_pos * 2;
                    data[dst_bit >> 6] |= code << (dst_bit & 63);
                    ++dest_pos;
                }
                origin_node = frag.origin_node;
            }
            cumulative += frag.len;
        }
    };
    CrumbledSubSeq64(size_t len, int skip_first_n, uint64_t *data,
                     uint64_t origin_node) noexcept
        : len{len},
          skip_first_n{skip_first_n},
          data(data),
          origin_node(origin_node) {}
    static CrumbledSubSeq64 copy_of(const CrumbledSubSeq64 &other) {
        const size_t words = storage_size(other.len + other.skip_first_n);
        uint64_t *copied_data = new uint64_t[words]();
        std::memcpy(copied_data, other.data, words * sizeof(uint64_t));
        return CrumbledSubSeq64(other.len, other.skip_first_n, copied_data,
                                other.origin_node);
    }
    CrumbledSubSeq64(const CrumbledSubSeq64 &other) = delete;
    CrumbledSubSeq64 &operator=(const CrumbledSubSeq64 &other) = delete;
    ~CrumbledSubSeq64() { delete[] data; }
    CrumbledSubSeq64 &operator=(CrumbledSubSeq64 &&other) noexcept {
        if (this != &other) {
            delete[] data;
            data = other.data;
            skip_first_n = other.skip_first_n;
            len = other.len;
            origin_node = other.origin_node;
            other.data = nullptr;
            other.origin_node = 0;
        }
        return *this;
    }

    char *convert_to_ascii(char *buffer) const {
        static constexpr const char *lut = "ACGT";
        size_t base_index = 0;

        auto scalar_emit = [&](size_t count) {
            const size_t bit_base = static_cast<size_t>(skip_first_n) * 2;
            for (size_t i = 0; i < count; ++i) {
                const size_t bit_index = bit_base + (base_index + i) * 2;
                const uint64_t word = data[bit_index >> 6];
                buffer[base_index + i] = lut[(word >> (bit_index & 63)) & 0x3];
            }
            base_index += count;
        };

#if 0 && (defined(__ARM_NEON) || defined(__ARM_NEON__))
        // Align to word boundary so the NEON path works for any skip_first_n.
        if (len >= bases_per_word) {
            const size_t total = len + skip_first_n;
            const size_t blocks =
                (total + bases_per_word - 1) / bases_per_word;
            const uint64_t *word_ptr = data;
            auto *out = reinterpret_cast<uint8_t *>(buffer);

            static const uint8x8_t ascii_table =
                {'A', 'C', 'G', 'T', 'A', 'C', 'G', 'T'};
            const uint8x8_t mask = vdup_n_u8(0x3);

            for (size_t b = 0; b < blocks; ++b) {
                const uint8x8_t packed =
                    vld1_u8(reinterpret_cast<const uint8_t *>(word_ptr));

                const uint8x8_t c0 = vand_u8(packed, mask);
                const uint8x8_t c1 = vand_u8(vshr_n_u8(packed, 2), mask);
                const uint8x8_t c2 = vand_u8(vshr_n_u8(packed, 4), mask);
                const uint8x8_t c3 = vand_u8(vshr_n_u8(packed, 6), mask);

                uint8x8x4_t expanded;
                expanded.val[0] = vtbl1_u8(ascii_table, c0);
                expanded.val[1] = vtbl1_u8(ascii_table, c1);
                expanded.val[2] = vtbl1_u8(ascii_table, c2);
                expanded.val[3] = vtbl1_u8(ascii_table, c3);

                vst4_u8(out, expanded);

                word_ptr++;
                out += bases_per_word;
            }
            buffer[total] = '\0';
            return buffer + skip_first_n;
        }
#endif
        scalar_emit(len);
        buffer[len] = '\0';
        return buffer;
    }

    [[nodiscard]] CrumbledSubSeq64 reverse_complement() const {
        CrumbledSubSeq64 rc{len, 0};
        for (size_t i = 0; i < len; ++i) {
            const size_t src_bit = (skip_first_n + (len - 1 - i)) * 2;
            const uint64_t code = (data[src_bit >> 6] >> (src_bit & 63)) & 0x3;
            const uint64_t complement = code ^ 0x3;
            const size_t dst_bit = i * 2;
            rc.data[dst_bit >> 6] |= complement << (dst_bit & 63);
        }
        return rc;
    }
    [[nodiscard]] CrumbledSubSeq64 reverse_complement_old() const {
        assert(skip_first_n >= 0);
        CrumbledSubSeq64 rc{len, 31};
        if (len < 32) {
            // Fix the skip_first
            uint64_t word = (data[0] >> skip_first_n * 2);
            if (storage_size(len + skip_first_n) > 1) {
                word |= (data[1] << (64 - skip_first_n * 2));
            }

            rc.data[0] = reverse_complement_kmer(word, len);
            rc.skip_first_n = 0;
            return rc;
        }
        size_t num_blocks = storage_size(len + skip_first_n);
        for (size_t i = 0; i < num_blocks - 1; ++i) {
            rc.data[num_blocks - i - 1] = reverse_complement_kmer(data[i], 32);
        }
        if ((len + skip_first_n) % 32 != 0)
            rc.data[0] = reverse_complement_kmer(data[num_blocks - 1],
                                                 (len + skip_first_n) % 32);
        rc.skip_first_n = (num_blocks * 32 - len) - skip_first_n;
        return rc;
    }
};
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
