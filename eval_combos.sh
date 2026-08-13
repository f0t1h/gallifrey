#!/usr/bin/env bash
# eval_combos.sh — run lasm over a grid of graph-cleaning / join feature
# combinations and score each with QUAST.
#
# Grid axes (all current experimental knobs, all default-off in lasm):
#   LASM_ADJ_JOIN          crumb adjacency joins (off|on)
#   --open-simple-loops    loop unrolling w/ multiplicity oracle (off|on)
#   LASM_UNROLL_TOL        unroll integer-fit tolerance (0.25|0.5; only
#                          meaningful with loops on)
#   --max-bubble-arm-length 0   disables bubble popping (bub on|off)
#
# Usage:
#   ./eval_combos.sh <input.bam|.fa> <reference.fa> <outdir> [-- extra lasm args...]
#
# Examples:
#   ./eval_combos.sh ~/Workspace/lasm-data/unmap18.bam \
#       ~/Workspace/lasm-data/chr18_chm13.fa /tmp/eval18 \
#       -- -l5 --dbg-lcmer-l 6 --olc-lcmer-l 6 --syncmer-k 11 \
#          --submer-s 3 --submer-t 7
#
#   ./eval_combos.sh ~/Workspace/lasm-data/chr22.tksm.hifi.fa \
#       ~/Workspace/lasm-data/chr22_hg002pat.fa /tmp/evalfa \
#       -- -l4 --syncmer-k 11 --submer-s 3 --submer-t 3
#
# Env knobs:
#   COMBO_ONLY=<name>   run a single combo (see names in summary)
#   JOBS=<n>            concurrent lasm runs (default 2; raise on big-RAM
#                       machines — whole-genome runs are memory-heavy)
#   LASM=<path>         lasm binary (default ./build/lasm)
#   THREADS="<flags>"   lasm thread flags (default "-t10 -@8")
#
# The summary includes lasm_s and quast_s wall times per combo, plus a
# host/input/date header line. Runs are resume-safe: re-running skips combos
# with an existing .metrics file, so an interrupted whole-genome grid
# continues where it stopped.
# Final ranked table: <outdir>/summary.txt (sorted by misassemblies asc,
# then NGA50 desc).
set -euo pipefail

INPUTS=${1:?usage: eval_combos.sh '<input1 [input2 ...]>' <ref> <outdir> [-- extra lasm args...]}
INPUTS=${INPUTS//\~\//$HOME\/}   # expand ~/ in every path
REF=${2:?missing reference}
OUT=${3:?missing outdir}
shift 3
EXTRA=()
if [ $# -gt 0 ] && [ "$1" = "--" ]; then
    shift
    EXTRA=("$@")
fi

LASM=${LASM:-./build/lasm}
JOBS=${JOBS:-2}
THREADS=${THREADS:-"-t10 -@8"}
mkdir -p "$OUT"

# Build the combo grid: name|extra_args|env.
COMBOS=()
for adj in off on; do
    for loops in off on; do
        for bub in on off; do
            if [ "$loops" = on ]; then tols="025 05"; else tols="025"; fi
            for tol in $tols; do
                name="adj_${adj}_loops_${loops}_bub_${bub}_tol${tol}"
                args=()
                envs=()
                [ "$adj" = on ] && envs+=(LASM_ADJ_JOIN=1)
                if [ "$loops" = on ]; then
                    args+=(--open-simple-loops)
                    [ "$tol" = 05 ] && envs+=(LASM_UNROLL_TOL=0.5)
                fi
                [ "$bub" = off ] && args+=(--max-bubble-arm-length 0)
                COMBOS+=("${name}|${args[*]-}|${envs[*]-}")
            done
        done
    done
done

if [ -n "${COMBO_ONLY:-}" ]; then
    FOUND=
    for c in "${COMBOS[@]}"; do
        [[ "$c" == "${COMBO_ONLY}|"* ]] && { COMBOS=("$c"); FOUND=1; break; }
    done
    [ -n "$FOUND" ] || { echo "unknown combo: $COMBO_ONLY" >&2; exit 2; }
fi

extract_metrics() {
    local name=$1
    local tsv="$OUT/quast_$name/report.tsv"
    if [ ! -f "$tsv" ]; then
        printf '%s\tFAIL\t0\t0\t0\t0\t0\t0\t0\t0\n' "$name" > "$OUT/$name.metrics"
        return 0
    fi
    awk -F'\t' -v n="$name" '
        $1=="# contigs" {c=$2}
        $1=="N50" {n50=$2}
        $1=="NGA50" {nga=$2}
        $1=="# misassemblies" {m=$2}
        $1=="# local misassemblies" {lm=$2}
        $1=="Genome fraction (%)" {gf=$2}
        $1=="Duplication ratio" {dr=$2}
        $1=="# mismatches per 100 kbp" {mm=$2}
        $1=="Largest contig" {lc=$2}
        END { printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n",
                    n, c, n50, nga, m, lm, gf, dr, mm, lc }
    ' "$tsv" > "$OUT/$name.metrics"
}

# (extract_metrics is also inlined in the xargs runner below; this copy is
# kept for single-shot manual use.)
export OUT REF INPUTS LASM THREADS EXTRA_STR="${EXTRA[*]-}"
printf '%s\n' "${COMBOS[@]}" |
    xargs -P "$JOBS" -n1 bash -c '
        spec=$1
        IFS="|" read -r name args envs <<< "$spec"
        [ -f "$OUT/$name.metrics" ] && { echo "skip $name (done)"; exit 0; }
        echo "run  $name"
        read -ra extra <<< "${EXTRA_STR:-}"
        t0=$(date +%s)
        # shellcheck disable=SC2086  # envs/extra/THREADS are intentionally word-split
        if ! env $envs "$LASM" hifi $INPUTS -o "$OUT/$name.gfa" \
                "${extra[@]}" $args $THREADS > "$OUT/$name.log" 2>&1; then
            printf "%s\tFAIL\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\n" "$name" > "$OUT/$name.metrics"
            exit 0
        fi
        t1=$(date +%s)
        lasm_s=$((t1 - t0))
        if [ ! -s "$OUT/$name.gfa" ]; then
            printf "%s\tEMPTY\t0\t0\t0\t0\t0\t0\t0\t0\t%s\t0\n" "$name" "$lasm_s" > "$OUT/$name.metrics"
            exit 0
        fi
        gfatools gfa2fa "$OUT/$name.gfa" > "$OUT/$name.fa"
        t2=$(date +%s)
        quast.py -r "$REF" "$OUT/$name.fa" -o "$OUT/quast_$name" -t 4 \
            --min-contig 500 --no-plots --no-html --no-icarus \
            > "$OUT/$name.quast.log" 2>&1
        t3=$(date +%s)
        quast_s=$((t3 - t2))
        tsv="$OUT/quast_$name/report.tsv"
        if [ ! -f "$tsv" ]; then
            printf "%s\tFAIL\t0\t0\t0\t0\t0\t0\t0\t0\t%s\t%s\n" "$name" "$lasm_s" "$quast_s" > "$OUT/$name.metrics"
            exit 0
        fi
        awk -F"\t" -v n="$name" -v ls="$lasm_s" -v qs="$quast_s" "
            \$1==\"# contigs\" {c=\$2}
            \$1==\"N50\" {n50=\$2}
            \$1==\"NGA50\" {nga=\$2}
            \$1==\"# misassemblies\" {m=\$2}
            \$1==\"# local misassemblies\" {lm=\$2}
            \$1==\"Genome fraction (%)\" {gf=\$2}
            \$1==\"Duplication ratio\" {dr=\$2}
            \$1==\"# mismatches per 100 kbp\" {mm=\$2}
            \$1==\"Largest contig\" {lc=\$2}
            END { printf \"%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n\",
                        n, c, n50, nga, m, lm, gf, dr, mm, lc, ls, qs }
        " "$tsv" > "$OUT/$name.metrics"
    ' _

{
    n_in=$(printf '%s\n' $INPUTS | wc -l | tr -d ' ')
    if [ "$n_in" = 1 ]; then in_label="$(basename $INPUTS)"; else in_label="$n_in inputs: $(basename ${INPUTS%% *})..."; fi
    printf '# host=%s  input=%s  date=%s\n' "$(hostname)" "$in_label" "$(date -u +%Y-%m-%dT%H:%MZ)"
    printf 'name\tcontigs\tN50\tNGA50\tmisasm\tlocals\tgenfrac\tdup\tmism_100kbp\tlargest\tlasm_s\tquast_s\n'
    cat "$OUT"/*.metrics | sort -t$'\t' -k5,5n -k4,4rn
} | column -t -s$'\t' > "$OUT/summary.txt"

cat "$OUT/summary.txt"
