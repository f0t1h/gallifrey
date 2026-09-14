# HiFi assembler

## Build

Requirements:

- CMake >= 3.20
- C++23 compiler
- `htslib`
- `zlib`
- POSIX threads

Build (defaults to `Release`):

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j
```

The main executable is:

```bash
./build/lasm
```

### Build types

Select a build type with `-DCMAKE_BUILD_TYPE=<type>`. Each writes to its own
build directory so you can keep them side by side.

| Type | Flags (summary) | 
|------|-----------------|
| `Release` (default) | `-O3 -g -ftree-vectorize -DNDEBUG` |
| `RelWithDebInfo` | `-O2 -g -pg -DNDEBUG` | 
| `Debug` | `-O0 -g -DDEBUG` + UBSan (`-fsanitize=undefined,null,bounds -fsanitize-trap=all`) | 
| `MinSizeRel` | `-Os -DNDEBUG` |

### Build options
- `LASM_NATIVE` — add `-march=native`

```bash
# Default build
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j

# Fastest local build for this machine
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DLASM_NATIVE=ON
cmake --build build -j

# Profiling build with native tuning (produces gmon.out after a run)
cmake -S . -B build-prof -DCMAKE_BUILD_TYPE=RelWithDebInfo -DLASM_NATIVE=ON
cmake --build build-prof -j
```

## Input

`lasm hifi` accepts one or more read files. Supported formats:

- BAM, CRAM, SAM (unmapped reads are fine; CRAM with mapped reads needs the reference available to htslib)
- FASTA and FASTQ, plain or gzipped

The format is detected from the file content, so any file name works. Read
names and quality values are ignored; only sequences are used. Different
formats can be mixed in one invocation.

```bash
lasm hifi <reads1> [reads2 ...] -o output.gfa -l <lcp_level> --density <kmer_density>
```

## Basic Usage

Minimal example:

```bash
build/lasm hifi reads.bam -o assembly.gfa -l 4 --density 0.05
```

Example that uses syncmer:

```bash
build/lasm hifi reads.bam \
  -o assembly.gfa \
  -l 3 \
  -t 8 \
  -@ 4 \
  --dbg-lcmer-l 8 \
  --syncmer-k 11 \
  --submer-s 4 \
  --submer-t 5 \
  --iterations 2 
```

## Main Parameters

Required:

- `input`: one or more BAM/CRAM/SAM or FASTA/FASTQ(.gz) read files
- `-o, --output`: output GFA path
- `-l, --level`: LCP/core depth

Core assembly controls:

- `--dbg-lcmer-l`: rolling `lcmer` window length used for DBG construction
- `--syncmer-k`: syncmer k-mer length
- `--density`: density-based syncmer selection
- `--submer-s`: submer size for mincode/open-syncmer rules
- `--submer-t`: submer position for mincode/open-syncmer rules
- `--spaced-seed-mask`: optional spaced-seed mask
- `--random-seed`: hash/random seed

Parallelism:

- `-t, --threads`: worker thread count
- `-@, --io-threads`: decompression thread count per input file (BAM/CRAM/SAM and gzipped FASTA/FASTQ; no effect on uncompressed FASTA/FASTQ)
- `--file-reading-mode`: `sequential` or `parallel`

Graph cleaning:

- `--iterations`: number of graph-cleaning passes (default 6)
- `--max-tip-length`: maximum tip length in cores
- `--max-bubble-arm-length`: maximum arm length (cores) of a bubble arm to pop; `0` disables bubble popping
- `--min-bubble-depth-difference-ratio`: pop the shallower arm when its depth ratio to the deeper arm is below this (default 2, i.e. depth is not a gate); `--min-bubble-length-difference`: maximum arm length ratio (default 2)

Layout:

- `--drop-crumbs` / `--keep-crumbs`: drop crumb nodes (too short for overlap detection) that have a single predecessor and continue into a longer node, instead of emitting them as `C` segments (default on; halves the contig count at a small genome-fraction cost)

Output and debugging:

- `--print-raw-dbg`: write the pre-cleaning graph with `.raw.gfa` suffix
- `--print-bi-dbg`: write the graph before layout with `.bi.gfa` suffix
- `--silent`: reduce progress logging

## Output Files

Primary output:

- `<output>`: final merged layout/consensus graph in GFA format
- M[0-9]* nodes indicate layout nodes which are strand resolved
- C[0-9]* nodes are not strand resolved (will have duplicates)

Optional auxiliary outputs:

- `<output>.raw.gfa`: raw DBG before cleaning, if `--print-raw-dbg` is enabled
- `<output>.bi.gfa`: cleaned DBG before layout merging if `--print-bi-dbg` is enabled
