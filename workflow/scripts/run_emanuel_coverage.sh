#!/usr/bin/env bash
# Score Emanuel Sonder's coverage simulations with amet, one amet run per ncpg,
# then combine the per-ncpg cell_feature outputs into one all_cells.tsv.gz.
#
# Usage:
#   run_emanuel_coverage.sh <amet> <sim_data_dir> <convert_sim.sh> \
#                           <out_base> <i_max_lag> <min_cpgs> <min_cells> <threads>
#
# The ncpg grid is read from the cpgPositions_<ncpg>.tsv files the generator
# wrote into <sim_data_dir>, so it stays in sync with the parameter grid.
set -euo pipefail

AMET="${1:?amet binary required}"
SIM_DATA="${2:?sim_data dir required}"
CONVERT="${3:?convert_sim.sh required}"
OUT_BASE="${4:?output base dir required}"
I_MAX_LAG="${5:?i-max-lag required}"
MIN_CPGS="${6:?min-cpgs-per-feature required}"
MIN_CELLS="${7:?min-cells-per-group required}"
THREADS="${8:?threads required}"

mkdir -p "$OUT_BASE"

ncpgs=()
for f in "$SIM_DATA"/cpgPositions_*.tsv; do
    [[ -e "$f" ]] || { echo "no cpgPositions_*.tsv in $SIM_DATA" >&2; exit 1; }
    n="$(basename "$f" .tsv)"
    ncpgs+=("${n#cpgPositions_}")
done

for ncpg in "${ncpgs[@]}"; do
    out="$OUT_BASE/$ncpg"
    stage="$out/stage"
    mkdir -p "$stage/cells"

    cell_files=( "$SIM_DATA"/sim_cell_*_*_"${ncpg}"_*.tsv )
    if (( ${#cell_files[@]} == 0 )); then
        echo "no sim_cell files for ncpg=$ncpg" >&2
        continue
    fi
    printf '%s\n' "${cell_files[@]}" | xargs -P "$THREADS" -I{} bash -c '
        src="$1"; stage="$2"; conv="$3"
        base="$(basename "$src" .tsv)"
        "$conv" "$src" "$stage/cells/${base}.allc.tsv"
    ' _ {} "$stage" "$CONVERT"

    # CpG reference: cpgPositions is 1-based, shift to 0-based.
    awk 'BEGIN{OFS="\t"} {print $1, $2 - 1}' \
        "$SIM_DATA/cpgPositions_${ncpg}.tsv" > "$stage/cpg.tsv"

    # One feature spanning every CpG (0-based half-open).
    printf "chrSim\t0\t%d\tall_%d\n" "$ncpg" "$ncpg" > "$stage/features.bed"

    # Manifest. group = transition matrix, the last filename field.
    {
        printf "cell_id\tpath\tgroup\n"
        for f in "$stage/cells"/sim_cell_*.allc.tsv; do
            cid="$(basename "$f" .allc.tsv)"
            printf "%s\t%s\t%s\n" "$cid" "$f" "${cid##*_}"
        done
    } > "$stage/cells.tsv"

    "$AMET" \
        --cells "$stage/cells.tsv" \
        --features "$stage/features.bed" \
        --cpg-reference "$stage/cpg.tsv" \
        --output-prefix "$out/amet" \
        --min-cpgs-per-feature "$MIN_CPGS" \
        --min-cells-per-group "$MIN_CELLS" \
        --i-max-lag "$I_MAX_LAG" \
        --threads "$THREADS"
    echo "scored ncpg=$ncpg"
done

# Combine per-ncpg cell_feature outputs; keep a single header.
combined="$OUT_BASE/all_cells.tsv.gz"
{
    first=1
    for ncpg in "${ncpgs[@]}"; do
        f="$OUT_BASE/$ncpg/amet.cell_feature.tsv.gz"
        if (( first )); then
            zcat "$f"
            first=0
        else
            zcat "$f" | tail -n +2
        fi
    done
} | gzip > "$combined"
echo "wrote $combined"
