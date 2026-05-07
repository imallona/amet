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
repo_root="$(cd "$(dirname "$0")/../.." && pwd)"
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

        echo "[sync] argelaguet done: $(ls "$out_root/cells" | wc -l) cells, $(ls "$out_root/features"/*.bed 2>/dev/null | wc -l) BEDs"
        ;;

    ecker|crc)
        echo "$dataset rsync not implemented yet" >&2
        exit 2
        ;;

    *)
        echo "unknown dataset: $dataset" >&2
        exit 1
        ;;
esac
