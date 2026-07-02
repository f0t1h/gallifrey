#pragma once

#ifndef _DUMP_CORES_H_
#define _DUMP_CORES_H_

#include <cstdint>
#include <string>

namespace CLI {
class App;
}

// "dump-cores" subcommand: owns its parameters and registers itself with the
// parent app on construction (see lasm.cpp main). Copy/move deleted so the
// callback's captured `this` cannot dangle.
class DumpCoresProgram {
    std::string input;
    std::string output;
    int syncmer_k = 11;
    int submer_s = 4;
    int submer_t = 5;
    int lcp_level = 4;
    uint64_t seed = 0;
    int lcp_low_level = 0;   // >0 selects LCP-subsample mode
    double density = 0.2;    // with lcp_low_level: fraction of cores kept

public:
    explicit DumpCoresProgram(CLI::App &parent);
    DumpCoresProgram(const DumpCoresProgram &) = delete;
    DumpCoresProgram &operator=(const DumpCoresProgram &) = delete;
};

#endif
