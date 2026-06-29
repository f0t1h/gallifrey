#pragma once

#ifndef _BENCH_LAYOUT_H_
#define _BENCH_LAYOUT_H_

#include <cstdint>
#include <string>

namespace CLI {
class App;
}

// "bench-layout" subcommand: measures sequence-data access speed for three
// CrumbledSubSeq64 storage layouts (scattered heap pointer, SSO-128 inline,
// and a contiguous arena) under sequential and shuffled access, on real cores
// extracted from a FASTA. Registers itself with the parent app on construction.
class BenchLayoutProgram {
    std::string input;
    int syncmer_k = 10;
    int submer_s = 4;
    int submer_t = 5;
    int lcp_level = 5;
    uint64_t seed = 0;
    uint64_t max_reads = 30000;
    int reps = 5;

public:
    explicit BenchLayoutProgram(CLI::App &parent);
    BenchLayoutProgram(const BenchLayoutProgram &) = delete;
    BenchLayoutProgram &operator=(const BenchLayoutProgram &) = delete;
};

#endif
