#!/usr/bin/env bash
# gfa_mc.sh <outdir> — per-combo M/C node and base counts from EXISTING GFA
# files. Backfill companion to eval_combos.sh: runs against a completed or
# in-flight grid without re-running any combo.
#
# Output: <outdir>/mc_report.txt (+ stdout), one row per <name>.gfa found.
set -euo pipefail

OUT=${1:?usage: gfa_mc.sh <outdir>}
{
    printf 'name\tM_nodes\tM_bases\tC_nodes\tC_bases\n'
    for gfa in "$OUT"/*.gfa; do
        [ -e "$gfa" ] || continue
        name=$(basename "$gfa" .gfa)
        awk -F'\t' -v n="$name" '
            $1=="S" && $2 ~ /^M/ { mn++; mb += length($3) }
            $1=="S" && $2 ~ /^C/ { cn++; cb += length($3) }
            END { printf "%s\t%d\t%d\t%d\t%d\n", n, mn+0, mb+0, cn+0, cb+0 }
        ' "$gfa"
    done
} | column -t -s$'\t' | tee "$OUT/mc_report.txt"
