#!/usr/bin/env bash
## Pull a prototype subset of cells from barbara to local.
## Read-only on barbara: only ssh + rsync, never writes there.
##
## Usage:  bash workflow/scripts/sync_from_barbara.sh argelaguet
##
## Reads workflow/config/datasets.yaml for paths and prototype subset.
## Writes to results/<dataset>/cells/.

set -euo pipefail

if [[ $# -lt 1 ]]; then
    echo "usage: $0 <argelaguet|ecker|crc>" >&2
    exit 1
fi

dataset="$1"
repo_root="$(cd "$(dirname "$0")/../../.." && pwd)"
config="$repo_root/workflow/config/datasets.yaml"
out_root="$repo_root/results/$dataset"
mkdir -p "$out_root/cells" "$out_root/features"

host=$(awk '/^barbara:/{f=1; next} f && /^  host:/{print $2; exit}' "$config")

case "$dataset" in
    argelaguet)
        meta_remote=$(awk '/^argelaguet:/{f=1; next} f && /^  meta_remote:/{print $2; exit}' "$config")
        cpg_remote=$(awk '/^argelaguet:/{f=1; next} f && /^  cpg_level_remote:/{print $2; exit}' "$config")
        features_remote=$(awk '/^argelaguet:/{f=1; next} f && /^  features_remote:/{print $2; exit}' "$config")
        cells_per_group=$(awk '/^prototype:/{f=1; next} f && /^  cells_per_group:/{print $2; exit}' "$config")

        echo "[sync] pulling argelaguet metadata"
        rsync -av "$host:$meta_remote" "$out_root/sample_metadata.txt"

        ## Prototype lineages must match argelaguet.proto_lineages in datasets.yaml.
        echo "[sync] selecting prototype cells (4 lineages x $cells_per_group)"
        awk -F'\t' -v n="$cells_per_group" '
            NR == 1 {
                for (i = 1; i <= NF; i++) col[$i] = i
                next
            }
            $col["pass_metQC"] == "TRUE" {
                lin = $col["lineage10x"]
                if (lin == "Epiblast" || lin == "Primitive_endoderm" \
                    || lin == "Nascent_mesoderm" || lin == "Rostral_neurectoderm") {
                    if (++count[lin] <= n) print $col["id_met"]
                }
            }
        ' "$out_root/sample_metadata.txt" > "$out_root/cells_to_pull.txt"

        echo "[sync] pulling $(wc -l < "$out_root/cells_to_pull.txt") cell files"
        sed 's,$,.tsv.gz,' "$out_root/cells_to_pull.txt" \
            | rsync -av --files-from=- "$host:$cpg_remote/" "$out_root/cells/"

        echo "[sync] pulling gastrulation feature BEDs"
        rsync -av --include='*.bed' --exclude='*' "$host:$features_remote/" "$out_root/features/"

        echo "[sync] pulling mm10 generic BEDs (matches yamet's _ARGELAGUET_MM10_ANNOTATIONS)"
        mkdir -p "$out_root/mm10"
        rsync -av \
            --include='h3k4me3.bed.gz' --include='h3k9me3.bed.gz' \
            --include='h3k4me1.bed.gz' --include='h3k27me3.bed.gz' \
            --include='h3k27ac.bed.gz' \
            --include='genes.bed.gz' --include='lines.bed.gz' \
            --include='sines.bed.gz' --include='promoters.bed.gz' \
            --exclude='*' \
            "$host:/home/imallona/src/yamet/workflow/mm10/" "$out_root/mm10/"

        echo "[sync] argelaguet done: $(ls "$out_root/cells" | wc -l) cells, $(ls "$out_root/features"/*.bed 2>/dev/null | wc -l) gastro BEDs, $(ls "$out_root/mm10"/*.bed.gz 2>/dev/null | wc -l) mm10 BEDs"
        ;;

    crc)
        raw_remote=$(awk '/^crc:/{f=1; next} f && /^  raw_remote:/{print $2; exit}' "$config")
        proto_patients=$(awk '/^crc:/{f=1; next} f && /^  proto_patients:/{p=1; next}
                              p && /^    -/{gsub(/^    - /, ""); print}' "$config" | tr '\n' '|' | sed 's/|$//')
        proto_locations=$(awk '/^crc:/{f=1; next} f && /^  proto_locations:/{p=1; next}
                               p && /^    -/{gsub(/^    - /, ""); print}' "$config" | tr '\n' '|' | sed 's/|$//')
        cells_per_group=$(awk '/^prototype:/{f=1; next} f && /^  cells_per_group:/{print $2; exit}' "$config")

        echo "[sync] CRC raw_remote=$raw_remote"
        mkdir -p "$out_root/raw"

        ## Build a small file list of singleC files matching prototype
        ## patients/locations, sorted by ascending file size on barbara so
        ## laptop smoke-tests pick the lightest cells first. Take up to N
        ## cells per (patient, location).
        echo "[sync] enumerating prototype-subset cells on barbara (smallest first)"
        ssh "$host" "ls -laS --reverse $raw_remote/" \
          | awk -F' ' -v p="$proto_patients" -v l="$proto_locations" -v n="$cells_per_group" '
              /singleC/ {
                  fname = $NF
                  split(fname, a, "_")
                  pat = a[3]; loc = a[4]; sub(/[0-9]+$/, "", loc)
                  if (pat ~ "^("p")$" && loc ~ "^("l")$") {
                      key = pat "_" loc
                      if (++count[key] <= n) print fname
                  }
              }
            ' > "$out_root/cells_to_pull.txt"

        echo "[sync] pulling $(wc -l < "$out_root/cells_to_pull.txt") singleC files"
        rsync -av --files-from="$out_root/cells_to_pull.txt" "$host:$raw_remote/" "$out_root/raw/"

        echo "[sync] crc done: $(ls "$out_root/raw" | wc -l) files in $out_root/raw"
        ;;

    ecker)
        raw_remote=$(awk '/^ecker:/{f=1; next} f && /^  raw_remote:/{print $2; exit}' "$config")
        meta_remote=$(awk '/^ecker:/{f=1; next} f && /^  meta_remote:/{print $2; exit}' "$config")
        nemo_meta_remote=$(awk '/^ecker:/{f=1; next} f && /^  nemo_meta_remote:/{print $2; exit}' "$config")

        echo "[sync] Ecker raw_remote=$raw_remote"
        mkdir -p "$out_root/raw"

        ## Originals only (paper xlsx + NeMo TSV). yamet's harmonized
        ## meta.tsv.gz is intentionally NOT pulled here.
        echo "[sync] pulling paper supplement xlsx"
        rsync -av "$host:$meta_remote" "$out_root/$(basename "$meta_remote")"
        echo "[sync] pulling original NeMo MOp_Metadata.tsv.gz"
        rsync -av "$host:$nemo_meta_remote" "$out_root/nemo_meta.tsv.gz"

        ## Pull a small fixed slice of tars (alphabetic head). The amet
        ## workflow's manifest rule merges xlsx+NeMo TSV downstream and
        ## prunes anything not in proto_cell_types, so under-coverage of a
        ## given sub_type just yields a smaller manifest -- never an error.
        n_tars=60
        echo "[sync] pulling first $n_tars tars (alphabetic, smoke-test slice)"
        ssh "$host" "ls $raw_remote/*.tar | head -n $n_tars" \
          | xargs -I{} basename {} > "$out_root/cells_to_pull.txt"
        rsync -av --files-from="$out_root/cells_to_pull.txt" "$host:$raw_remote/" "$out_root/raw/"

        echo "[sync] ecker done: $(ls "$out_root/raw"/*.tar 2>/dev/null | wc -l) tars in $out_root/raw"
        ;;

    *)
        echo "unknown dataset: $dataset" >&2
        exit 1
        ;;
esac
