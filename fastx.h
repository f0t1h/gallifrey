#ifndef _LASM_FASTX_H_
#define _LASM_FASTX_H_

#include <fcntl.h>
#include <htslib/bgzf.h>
#include <htslib/hfile.h>
#include <htslib/hts.h>
#include <unistd.h>

#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>

#if defined(__ARM_NEON) || defined(__ARM_NEON__)
#include <arm_neon.h>
#elif defined(__SSSE3__)
#include <tmmintrin.h>
#endif

// Native FASTA/FASTQ input path.
//
// htslib can read FASTA/FASTQ, but sam_read1() pays for it: a kstring copy per
// line, a full bam1_t rebuild, and an l_qseq-byte quality array the assembler
// never looks at. The assembler only needs the packed 4-bit sequence and its
// length, so this reader scans blocks with memchr and packs bases straight into
// the record buffer in the *same* nt16 nibble layout bam_get_seq() uses (A=1,
// C=2, G=4, T=8, other=15, high nibble first). Everything downstream --
// onehot_to_rle(), CrumbledSubSeq64's from_nibbles ctor -- is therefore
// byte-identical to the htslib path.

#ifndef FASTX_BLOCK_SIZE
#define FASTX_BLOCK_SIZE (1u << 20)
#endif

// Recycled sequence record: the FASTX counterpart of bam1_t. Buffers grow only.
struct seq_record {
    uint8_t *nibbles = nullptr;
    uint32_t len = 0;  // bases
    uint32_t cap = 0;  // allocated bytes

    void reserve_bases(uint32_t bases) {
        // +8 bytes of slack so nibble writers/readers that round up to a word
        // never touch unallocated memory.
        const uint32_t need = bases / 2 + 8;
        if (need <= cap) return;
        uint32_t ncap = cap ? cap : 8192;
        while (ncap < need) ncap += ncap / 2;
        uint8_t *p = static_cast<uint8_t *>(realloc(nibbles, ncap));
        if (!p) {
            fprintf(stderr, "seq_record: out of memory (%u bytes)\n", ncap);
            abort();
        }
        nibbles = p;
        cap = ncap;
    }
};

inline seq_record *seq_record_init() { return new seq_record(); }

inline void seq_record_destroy(seq_record *r) {
    if (!r) return;
    free(r->nibbles);
    delete r;
}

enum class input_format { htslib, fastx_plain, fastx_compressed, invalid };

// Content sniffing (magic bytes), not file extension: htslib's detector already
// distinguishes bam/cram/sam/fasta/fastq and their compression.
inline input_format detect_input_format(const std::string &path) {
    hFILE *hf = hopen(path.c_str(), "r");
    if (!hf) return input_format::invalid;
    htsFormat fmt;
    const int rc = hts_detect_format(hf, &fmt);
    if (hclose(hf) != 0 || rc < 0) return input_format::invalid;
    if (fmt.format == fasta_format || fmt.format == fastq_format) {
        return fmt.compression == no_compression ? input_format::fastx_plain
                                                 : input_format::fastx_compressed;
    }
    return input_format::htslib;
}

// Block source: raw read(2) for uncompressed files (zero copies beyond the
// kernel), htslib's BGZF for gzip/bgzf (multi-threaded inflate when the file is
// bgzf-blocked).
class fastx_source {
    int fd_ = -1;
    BGZF *bgz_ = nullptr;

public:
    fastx_source() = default;
    fastx_source(const fastx_source &) = delete;
    fastx_source &operator=(const fastx_source &) = delete;
    ~fastx_source() { close(); }

    bool open(const std::string &path, bool compressed, int io_threads) {
        if (compressed) {
            bgz_ = bgzf_open(path.c_str(), "r");
            if (!bgz_) return false;
            if (io_threads > 1) bgzf_mt(bgz_, io_threads, 256);
            return true;
        }
        fd_ = ::open(path.c_str(), O_RDONLY);
        if (fd_ < 0) return false;
#if defined(POSIX_FADV_SEQUENTIAL)
        posix_fadvise(fd_, 0, 0, POSIX_FADV_SEQUENTIAL);
#endif
        return true;
    }

    // Returns bytes read, 0 at EOF, -1 on error.
    int64_t read(uint8_t *buf, size_t n) {
        if (bgz_) return bgzf_read(bgz_, buf, n);
        while (true) {
            const ssize_t got = ::read(fd_, buf, n);
            if (got < 0 && errno == EINTR) continue;
            return got;
        }
    }

    void close() {
        if (bgz_) {
            bgzf_close(bgz_);
            bgz_ = nullptr;
        }
        if (fd_ >= 0) {
            ::close(fd_);
            fd_ = -1;
        }
    }
};

class fastx_reader {
    // kseq record grammar, driven incrementally so a record may straddle any
    // number of blocks: header line -> sequence lines -> ('+' line -> quality
    // bytes) for FASTQ. A record ends at a line starting with '>' or '@', or
    // once `len` quality bytes have been consumed.
    enum class state { start, header, seq, plus, qual };

    fastx_source src_;
    uint8_t *buf_ = nullptr;
    size_t pos_ = 0;
    size_t end_ = 0;
    state st_ = state::start;
    bool line_start_ = true;
    bool eof_ = false;
    bool failed_ = false;
    uint32_t qual_got_ = 0;
    size_t block_size_ = FASTX_BLOCK_SIZE;
    std::string name_;
    bool name_done_ = false;

    // Refills the block buffer. Returns false at EOF or on error.
    bool refill() {
        if (eof_) return false;
        const int64_t got = src_.read(buf_, block_size_);
        if (got <= 0) {
            eof_ = true;
            if (got < 0) failed_ = true;
            return false;
        }
        pos_ = 0;
        end_ = static_cast<size_t>(got);
        return true;
    }

    static inline void pack_one(uint8_t *dst, uint32_t &len, uint8_t c) {
        const uint8_t nib = seq_nt16_table[c];
        if (len & 1) {
            dst[len >> 1] = static_cast<uint8_t>(dst[len >> 1] | nib);
        } else {
            dst[len >> 1] = static_cast<uint8_t>(nib << 4);
        }
        ++len;
    }

    // Packs [p, q) of ASCII bases onto the end of `rec`, dropping any byte at
    // or below ' ' (CR, LF, spaces, tabs). This is the hot loop of the reader:
    // at ~1 byte per base it, not the I/O, bounds throughput, hence the SIMD
    // block for the overwhelmingly common case (16 plain ACGT bases landing on
    // an even base offset).
    static void pack_bases(seq_record &rec, const uint8_t *p, const uint8_t *q) {
        rec.reserve_bases(rec.len + static_cast<uint32_t>(q - p));
        uint8_t *dst = rec.nibbles;
        uint32_t len = rec.len;
#if defined(__ARM_NEON) || defined(__ARM_NEON__)
        // (c >> 1) & 3 maps A/C/T/G (either case) to 0/1/2/3; the table turns
        // that into the nt16 one-hot nibble.
        static const uint8_t code2nib[16] = {1, 2, 8, 4, 0, 0, 0, 0,
                                             0, 0, 0, 0, 0, 0, 0, 0};
        const uint8x16_t lut = vld1q_u8(code2nib);
        while ((len & 1) == 0 && q - p >= 16) {
            const uint8x16_t v = vld1q_u8(p);
            const uint8x16_t lower = vorrq_u8(v, vdupq_n_u8(0x20));
            const uint8x16_t ok = vorrq_u8(
                vorrq_u8(vceqq_u8(lower, vdupq_n_u8('a')),
                         vceqq_u8(lower, vdupq_n_u8('c'))),
                vorrq_u8(vceqq_u8(lower, vdupq_n_u8('g')),
                         vceqq_u8(lower, vdupq_n_u8('t'))));
            if (vminvq_u8(ok) != 0xFF) {
                // Whitespace, N, or IUPAC in this window: pack it one base at
                // a time, then resume the vector loop.
                const uint8_t *stop = p + 16;
                for (; p != stop; ++p) {
                    if (*p > ' ') pack_one(dst, len, *p);
                }
                continue;
            }
            const uint8x16_t code =
                vandq_u8(vshrq_n_u8(v, 1), vdupq_n_u8(3));
            const uint8x16_t nib = vqtbl1q_u8(lut, code);
            const uint8x16_t hi = vuzp1q_u8(nib, nib);  // even bases
            const uint8x16_t lo = vuzp2q_u8(nib, nib);  // odd bases
            vst1_u8(dst + (len >> 1),
                    vget_low_u8(vorrq_u8(vshlq_n_u8(hi, 4), lo)));
            len += 16;
            p += 16;
        }
#elif defined(__SSSE3__)
        static const uint8_t code2nib[16] = {1, 2, 8, 4, 0, 0, 0, 0,
                                             0, 0, 0, 0, 0, 0, 0, 0};
        const __m128i lut =
            _mm_loadu_si128(reinterpret_cast<const __m128i *>(code2nib));
        // Gather the even / odd lanes of the nibble vector into the low 8
        // bytes so a single 64-bit store emits 16 packed bases.
        const __m128i even_idx = _mm_setr_epi8(0, 2, 4, 6, 8, 10, 12, 14, -1,
                                               -1, -1, -1, -1, -1, -1, -1);
        const __m128i odd_idx = _mm_setr_epi8(1, 3, 5, 7, 9, 11, 13, 15, -1,
                                              -1, -1, -1, -1, -1, -1, -1);
        while ((len & 1) == 0 && q - p >= 16) {
            const __m128i v =
                _mm_loadu_si128(reinterpret_cast<const __m128i *>(p));
            const __m128i lower = _mm_or_si128(v, _mm_set1_epi8(0x20));
            const __m128i ok = _mm_or_si128(
                _mm_or_si128(_mm_cmpeq_epi8(lower, _mm_set1_epi8('a')),
                             _mm_cmpeq_epi8(lower, _mm_set1_epi8('c'))),
                _mm_or_si128(_mm_cmpeq_epi8(lower, _mm_set1_epi8('g')),
                             _mm_cmpeq_epi8(lower, _mm_set1_epi8('t'))));
            if (_mm_movemask_epi8(ok) != 0xFFFF) {
                const uint8_t *stop = p + 16;
                for (; p != stop; ++p) {
                    if (*p > ' ') pack_one(dst, len, *p);
                }
                continue;
            }
            const __m128i code =
                _mm_and_si128(_mm_srli_epi16(v, 1), _mm_set1_epi8(3));
            const __m128i nib = _mm_shuffle_epi8(lut, code);
            const __m128i hi = _mm_shuffle_epi8(nib, even_idx);
            const __m128i lo = _mm_shuffle_epi8(nib, odd_idx);
            // Nibbles are <= 8, so the 16-bit shift cannot carry across bytes.
            _mm_storel_epi64(reinterpret_cast<__m128i *>(dst + (len >> 1)),
                             _mm_or_si128(_mm_slli_epi16(hi, 4), lo));
            len += 16;
            p += 16;
        }
#endif
        for (; p != q; ++p) {
            const uint8_t c = *p;
            if (c <= ' ') continue;
            pack_one(dst, len, c);
        }
        rec.len = len;
    }

public:
    // `block_size` is a tuning/testing knob: record boundaries are handled
    // incrementally, so any positive size parses identically.
    explicit fastx_reader(size_t block_size = FASTX_BLOCK_SIZE)
        : block_size_(block_size) {}
    fastx_reader(const fastx_reader &) = delete;
    fastx_reader &operator=(const fastx_reader &) = delete;
    ~fastx_reader() { free(buf_); }

    bool open(const std::string &path, bool compressed, int io_threads) {
        if (!buf_) {
            buf_ = static_cast<uint8_t *>(malloc(block_size_));
            if (!buf_) return false;
        }
        return src_.open(path, compressed, io_threads);
    }

    bool failed() const { return failed_; }

    // Name of the record returned by the last successful next(): the header
    // line truncated at the first whitespace.
    const std::string &name() const { return name_; }

    // Fills `rec` with the next record's packed sequence. Returns false at EOF.
    bool next(seq_record &rec) {
        rec.len = 0;
        while (true) {
            if (pos_ == end_ && !refill()) {
                // End of stream: emit whatever the last record accumulated.
                st_ = state::start;
                if (rec.len == 0) return false;
                if (rec.len & 1) rec.nibbles[rec.len >> 1] &= 0xF0;
                return true;
            }
            const uint8_t *base = buf_;
            switch (st_) {
                case state::start: {
                    // Skip blank lines, then consume the record marker.
                    while (pos_ < end_ && base[pos_] <= ' ') ++pos_;
                    if (pos_ == end_) break;
                    const uint8_t c = base[pos_];
                    if (c != '>' && c != '@') {
                        failed_ = true;
                        eof_ = true;
                        return false;
                    }
                    ++pos_;
                    name_.clear();
                    name_done_ = false;
                    st_ = state::header;
                    break;
                }
                case state::header: {
                    const uint8_t *nl = static_cast<const uint8_t *>(
                        memchr(base + pos_, '\n', end_ - pos_));
                    const size_t stop =
                        nl ? static_cast<size_t>(nl - base) : end_;
                    if (!name_done_) {
                        // Name is the header up to the first whitespace; the
                        // comment after it is discarded.
                        size_t w = pos_;
                        while (w < stop && base[w] > ' ') ++w;
                        name_.append(reinterpret_cast<const char *>(base + pos_),
                                     w - pos_);
                        if (w < stop) name_done_ = true;
                    }
                    if (!nl) {
                        pos_ = end_;
                        break;
                    }
                    pos_ = stop + 1;
                    st_ = state::seq;
                    line_start_ = true;
                    break;
                }
                case state::plus: {
                    const uint8_t *nl = static_cast<const uint8_t *>(
                        memchr(base + pos_, '\n', end_ - pos_));
                    if (!nl) {
                        pos_ = end_;
                        break;
                    }
                    pos_ = static_cast<size_t>(nl - base) + 1;
                    st_ = state::qual;
                    qual_got_ = 0;
                    break;
                }
                case state::seq: {
                    if (line_start_) {
                        const uint8_t c = base[pos_];
                        if (c == '>' || c == '@') {
                            // Next record's header; leave it unconsumed.
                            st_ = state::start;
                            if (rec.len & 1) rec.nibbles[rec.len >> 1] &= 0xF0;
                            return true;
                        }
                        if (c == '+') {
                            st_ = state::plus;
                            break;
                        }
                    }
                    const uint8_t *nl = static_cast<const uint8_t *>(
                        memchr(base + pos_, '\n', end_ - pos_));
                    const size_t stop =
                        nl ? static_cast<size_t>(nl - base) : end_;
                    pack_bases(rec, base + pos_, base + stop);
                    pos_ = nl ? stop + 1 : end_;
                    line_start_ = (nl != nullptr);
                    break;
                }
                case state::qual: {
                    // Quality is never used; only its length matters (it must
                    // equal the sequence length, which is what makes a quality
                    // line starting with '@' unambiguous). Skip it in spans:
                    // one memchr finds the next newline, everything before it
                    // counts. Single-line quality is one memchr over len
                    // bytes; multi-line quality still works. '\r' never counts.
                    while (pos_ < end_ && qual_got_ < rec.len) {
                        const size_t need = rec.len - qual_got_;
                        const size_t span =
                            end_ - pos_ < need ? end_ - pos_ : need;
                        const uint8_t *p = base + pos_;
                        const uint8_t *nl =
                            static_cast<const uint8_t *>(memchr(p, '\n', span));
                        if (!nl) {
                            // A '\r' ending the span is a CRLF split across
                            // blocks (or the need-cut landing on one).
                            qual_got_ += static_cast<uint32_t>(
                                p[span - 1] == '\r' ? span - 1 : span);
                            pos_ += span;
                            continue;
                        }
                        size_t n = static_cast<size_t>(nl - p);
                        if (n > 0 && nl[-1] == '\r') --n;
                        qual_got_ += static_cast<uint32_t>(n);
                        pos_ = static_cast<size_t>(nl - base) + 1;
                    }
                    if (qual_got_ == rec.len) {
                        st_ = state::start;
                        if (rec.len & 1) rec.nibbles[rec.len >> 1] &= 0xF0;
                        return true;
                    }
                    break;
                }
            }
        }
    }
};

#endif
