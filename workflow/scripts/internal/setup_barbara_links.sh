#!/usr/bin/env bash
## Set up amet's results/{dataset}/ directories as symlinks to the existing
## yamet data tree on barbara, so we don't re-download anything.
##
## Run this on barbara *before* `make all`. It is idempotent: existing
## symlinks pointing at the right targets are left in place; existing
## non-symlinks are left in place too (we never overwrite real data).
##
## Source paths assume yamet is at $YAMET (default ~/src/yamet/workflow).
## Override via env: YAMET=/some/other/path bash workflow/scripts/setup_barbara_links.sh
##
## After this script, cells/features/raw paths under results/ resolve to the
## yamet directories. amet's snakemake rules (filter_argelaguet_metadata,
## make_manifest_*, ecker_extract_tar, ...) read those files transparently.

set -euo pipefail

YAMET="${YAMET:-$HOME/src/yamet/workflow}"
YAMET_HG19_CURATED="${YAMET_HG19_CURATED:-$HOME/src/yamet/hg19}"

repo_root="$(cd "$(dirname "$0")/../../.." && pwd)"
res="$repo_root/results"

link() {
    local src="$1" dst="$2"
    if [[ ! -e "$src" ]]; then
        echo "[setup] skip: $src does not exist"
        return
    fi
    if [[ -L "$dst" ]]; then
        local cur
        cur="$(readlink -f "$dst")"
        if [[ "$cur" == "$(readlink -f "$src")" ]]; then
            echo "[setup] ok: $dst -> $cur"
            return
        fi
        echo "[setup] replacing stale symlink $dst (was -> $cur)"
        rm "$dst"
    elif [[ -e "$dst" ]]; then
        echo "[setup] skip: $dst exists and is not a symlink (not touching)"
        return
    fi
    mkdir -p "$(dirname "$dst")"
    ln -s "$src" "$dst"
    echo "[setup] linked: $dst -> $src"
}

echo "[setup] yamet root: $YAMET"
echo "[setup] amet results root: $res"

## Argelaguet
link "$YAMET/argelaguet/met/cpg_level"           "$res/argelaguet/cells"
link "$YAMET/argelaguet/sample_metadata.txt"     "$res/argelaguet/sample_metadata.txt"
link "$YAMET/argelaguet/features/genomic_contexts" "$res/argelaguet/features"
link "$YAMET/mm10"                               "$res/argelaguet/mm10"

## CRC. hg19 = yamet/workflow/hg19 (cpgIslandExt, SCNAs, genes/lines/sines as
## .bed.gz). hg19_curated = yamet/hg19 (chromHMM, ChIP, lamin, PMD as .bed).
link "$YAMET/data/crc/raw"                       "$res/crc/raw"
link "$YAMET/hg19"                               "$res/crc/hg19"
link "$YAMET_HG19_CURATED"                       "$res/crc/hg19_curated"

## Ecker. mm10 shares the same dir as argelaguet (yamet's curated mm10 tree
## with ENCODE ChIP, UCSC repeats, and genes/promoters built off Gencode).
link "$YAMET/data/brain/raw"                     "$res/ecker/raw"
link "$YAMET/mm10"                               "$res/ecker/mm10"
link "$YAMET/data/brain/raw/41586_2020_3182_MOESM9_ESM.xlsx" \
     "$res/ecker/41586_2020_3182_MOESM9_ESM.xlsx"
link "$YAMET/data/brain/raw/MOp_Metadata.tsv.gz" "$res/ecker/nemo_meta.tsv.gz"

echo "[setup] done."
