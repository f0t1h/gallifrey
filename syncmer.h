#ifndef _SYNCMER_H_
#define _SYNCMER_H_

#include <cassert>
#include <concepts>
#include <cstdint>
#include <limits>

#include "asm.h"
#include "gtl/vector.hpp"
#define KOMIHASH_NS_CUSTOM komihash
#include <komihash.h>
namespace syncmer {
template <class Rule>
concept SyncmerRule = requires(Rule k) {
    { k(unsigned{}) } -> std::same_as<bool>;
};
template <class HF, class INT_TYPE>
concept HashFunction = requires(HF h) {
    { HF{INT_TYPE{}} } -> std::same_as<HF>;
    { h(INT_TYPE{}) } -> std::same_as<INT_TYPE>;
};

struct Komi {
    uint64_t seed;
    Komi(uint64_t seed) : seed{seed} {}
    uint64_t operator()(uint64_t val) const {
        return komihash::komihash(&val, sizeof(uint64_t), seed);
    }
};
// uint32_t MurmurHash3_32(const void *key, int len, uint32_t seed);
struct MurMur {
    static uint32_t MurmurHash3_32(const void *key, int len, uint32_t seed) {
        const uint8_t *data = (const uint8_t *)key;
        const int nblocks = len / 4;

        uint32_t h1 = seed;

        const uint32_t c1 = 0xcc9e2d51;
        const uint32_t c2 = 0x1b873593;

        // Body: Process blocks of 4 bytes at a time
        const uint32_t *blocks = (const uint32_t *)(data + nblocks * 4);

        for (int i = -nblocks; i; i++) {
            uint32_t k1 = blocks[i];

            k1 *= c1;
            k1 = (k1 << 15) | (k1 >> (32 - 15));
            k1 *= c2;

            h1 ^= k1;
            h1 = (h1 << 15) | (h1 >> (32 - 15));
            h1 = h1 * 5 + 0xe6546b64;
        }

        // Tail: Process remaining bytes
        const uint8_t *tail = (const uint8_t *)(data + nblocks * 4);

        uint32_t k1 = 0;

        switch (len & 3) {
            case 3:
                k1 ^= tail[2] << 16;

            case 2:
                k1 ^= tail[1] << 8;

            case 1:
                k1 ^= tail[0];
                k1 *= c1;
                k1 = (k1 << 15) | (k1 >> (32 - 15));
                k1 *= c2;
                h1 ^= k1;
        }

        // Finalization: Mix the hash to ensure the last few bits are fully
        // mixed
        h1 ^= len;

        /* fmix32 */
        h1 ^= h1 >> 16;
        h1 *= 0x85ebca6b;
        h1 ^= h1 >> 13;
        h1 *= 0xc2b2ae35;
        h1 ^= h1 >> 16;
        return h1;
    }
    unsigned seed;
    MurMur(unsigned seed) : seed{seed} {}
    unsigned operator()(unsigned value) const {
        return MurmurHash3_32((const char *)&value, sizeof(value), seed);
    }
};
struct Identity {
    unsigned seed;
    Identity(unsigned seed) : seed{seed} {}
    unsigned operator()(unsigned value) const { return value; }
};
template <class INT_TYPE, HashFunction<INT_TYPE> HASH_FUNCTION>
struct DensitySyncmerRule {
    using int_type = INT_TYPE;
    syncmer_params cfg;
    int_type seed;
    HASH_FUNCTION hf;
    int_type mask;
    constexpr int_type k() const { return cfg.k; }
    constexpr double density() const { return cfg.density.value(); }
    constexpr int_type threshold() const {
        return density() * std::numeric_limits<int_type>::max();
    }
    DensitySyncmerRule(const syncmer_params &cfg, int_type seed)
        : cfg{cfg},
          seed{seed},
          hf{seed},
          mask{(int_type)cfg.compute_spaced_mask()} {}
    DensitySyncmerRule(int_type k, double density, INT_TYPE seed, INT_TYPE mask)
        : cfg{.k = k, .density = density}, seed{seed}, hf{seed}, mask{mask} {}
    constexpr bool operator()(int_type val) const {
        return hf(val & mask) < threshold();
    }
};

// Builds a LuT for the given syncmer rule should work fine for K <= 11, 12
template <class INT_TYPE, HashFunction<INT_TYPE> HASH_FUNCTION>
struct GenericSyncmerRule {
    static constexpr double pure_lcp_density_threshold = 0.5;
    using int_type = INT_TYPE;
    syncmer_params cfg;
    int_type seed;
    HASH_FUNCTION hf;
    int_type s_mask;
    int_type k_mask;
    int64_t table_size;
    bool pure_lcp;
    gtl::vector<bool> table;

    constexpr int_type k() const { return cfg.k; }
    constexpr int_type s() const { return cfg.mincode.s.value(); }
    constexpr int_type t() const { return cfg.mincode.t.value(); }
    constexpr double density() const { return cfg.density.value(); }
    constexpr int_type threshold() const {
        return density() * std::numeric_limits<int_type>::max();
    }
    constexpr auto type() const -> syncmer_params::typeenum {
        return cfg.type();
    }
    GenericSyncmerRule(const syncmer_params &cfg, int_type seed)
        : cfg{cfg},
          seed{seed},
          hf{seed},
          k_mask{(INT_TYPE)cfg.compute_spaced_mask()},
          table_size{2 << (2 * k() - 1)},
          pure_lcp(false),
          table(table_size, 0) {
        if (cfg.type() == syncmer_params::typeenum::LcpSubsample) {
            // Core selection is done by LCP+subsample; the lookup table is
            // never queried, so skip building it (and its misleading density
            // log).
            return;
        }
        if (cfg.mincode.s.has_value()) {
            s_mask = (1U << s()) - 1U;
        } else {
            s_mask = -1;
        }
        size_t tc = 0;
        for (int i = 0; i < table_size; ++i) {
            table[i] = is_syncmer(i & k_mask);
            tc += table[i];
        }
        fprintf(stderr, "Real K-mer density: %lf\n", (double)tc / table_size);
        if ((double)tc / table_size > pure_lcp_density_threshold) {
            fprintf(stderr, "Density of %lf > %lf, switching to pure LCP!\n",
                    (double)tc / table_size, pure_lcp_density_threshold);
            pure_lcp = true;
        }
    }
    constexpr bool operator()(int_type val) const { return table[val]; }
    constexpr bool is_syncmer(int_type val) const {
        switch (type()) {
            case syncmer_params::typeenum::MinCodeClosed:  // Intentional
                                                           // fallthrough
            case syncmer_params::typeenum::MinCodeOpen: {
                int_type min_index = -1;
                int_type min_val = std::numeric_limits<int_type>::max();
                for (int_type i = 0; i < k() - s(); ++i) {
                    int_type masked_s = (val >> i) & s_mask;
                    int_type hashed_s = hf(masked_s);
                    if (hashed_s < min_val) {
                        min_index = i;
                        min_val = hashed_s;
                    }
                }
                if (type() == syncmer_params::typeenum::MinCodeClosed) {
                    return min_index == 0 || min_index == k() - s() - 1;
                } else {
                    return t() == min_index;
                }
            }
            case syncmer_params::typeenum::Density: {
                return hf(val & k_mask) < threshold();
            }
            case syncmer_params::typeenum::LcpSubsample:
                // Core selection is done by regular LCP + subsampling, not by
                // this k-mer functor; it is constructed but never queried.
                return false;
            case syncmer_params::typeenum::Invalid:
            default:
                // Unreachable: syncmer_params::validate() rejects Invalid
                // before any functor is constructed.
                assert(false && "syncmer type validated before construction");
                __builtin_unreachable();
        }
    }
};

}  // namespace syncmer

#endif
