#!/usr/bin/env bash
# Convert one yamet simulation 5-col cell file to an allc-formatted file amet can read.
#
# yamet cell format: chr, pos, m, t, rate  (per row = one CpG)
# allc format:       chr, pos_1based, strand, context, m, t, methylated_flag
#
# The sim_data files use pos_1based starting at 1; the within_cell/between_cell
# files use pos_0based starting at 0. We auto-detect from the first non-empty row
# and emit pos as 1-based so amet's allc parser maps it back to the correct
# 0-based CpG start (cpg_start = allc_pos - 1).
set -euo pipefail

if [[ $# -ne 2 ]]; then
    echo "usage: $0 <in.tsv> <out.allc.tsv>" >&2
    exit 2
fi

in="$1"
out="$2"

first_pos="$(awk 'NF >= 5 {print $2; exit}' "$in")"
if [[ "$first_pos" == "0" ]]; then
    shift_to_1based=1
else
    shift_to_1based=0
fi

awk -v s="$shift_to_1based" 'BEGIN{OFS="\t"} NF >= 5 {
    pos = $2 + s
    flag = ($3 > 0) ? 1 : 0
    print $1, pos, "+", "CG", $3, $4, flag
}' "$in" > "$out"
