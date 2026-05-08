"""
Liu et al. 2021 mouse brain methylation atlas (CEMBA).

DNA methylation atlas of the mouse brain at single-cell resolution.
Liu et al. https://www.nature.com/articles/s41586-020-03182-8
NeMo: u19_cemba sncell mCseq. Assembly: mm10, Ensembl-style chrom naming
(1, 2, ..., X, Y, L, M).

Per-cell raw files are tar archives containing one methylpy/allc tsv.gz with
the standard 7-column allc layout. amet reads the inner tsv.gz via its 'allc'
parser; the workflow's only preprocessing is a tar extraction (no sequence
content is rewritten).

Two acquisition paths:

  (a) rsync_from_barbara: pull raw .tar files already present on barbara into
      results/ecker/raw/ and extract per-cell tsv.gz into results/ecker/cells/.
  (b) download_from_nemo: cold-start path mirroring yamet's download chain.
      Pulls metadata (NeMo TSV + paper xlsx supplement), derives the per-tar
      URL list, then wget-downloads tars into results/ecker/raw/.
"""

ECKER_DATA = op.join(RESULTS, "ecker")
ECKER_RAW = op.join(ECKER_DATA, "raw")
ECKER_CELLS = op.join(ECKER_DATA, "cells")
ECKER_RUN = op.join(RESULTS, config["ecker"]["run_name"])
ECKER_RUN_NAME = config["ecker"]["run_name"]


rule ecker_download_nemo_metadata:
    """Cold-start: NeMo metadata TSV (cell -> tar URL mapping)."""
    conda:
        op.join("..", "envs", "bedtools.yml")
    output:
        meta = op.join(ECKER_DATA, "nemo_meta.tsv.gz"),
    params:
        url = "https://data.nemoarchive.org/biccn/grant/u19_cemba/cemba/epigenome/sncell/mCseq/mouse/processed/analysis/EckerRen_Mouse_MOp_methylation_ATAC/metadata/mc/MOp_Metadata.tsv.gz",
    log:
        op.join(ECKER_DATA, "logs", "download_nemo_meta.log"),
    shell:
        """
        mkdir -p $(dirname {output.meta})
        curl -sSL {params.url} -o {output.meta} 2> {log}
        """


rule ecker_download_paper_metadata:
    """Cold-start: paper supplementary xlsx with sub_type / cell_class / region."""
    conda:
        op.join("..", "envs", "bedtools.yml")
    output:
        xlsx = op.join(ECKER_DATA, "41586_2020_3182_MOESM9_ESM.xlsx"),
    params:
        url = "https://static-content.springer.com/esm/art%3A10.1038%2Fs41586-020-03182-8/MediaObjects/41586_2020_3182_MOESM9_ESM.xlsx",
    log:
        op.join(ECKER_DATA, "logs", "download_paper_meta.log"),
    shell:
        """
        mkdir -p $(dirname {output.xlsx})
        curl -sSL {params.url} -o {output.xlsx} 2> {log}
        """


rule ecker_compose_metadata:
    """Merge NeMo TSV + paper xlsx into a single cells metadata table.
    Keeps the columns needed for downstream filtering and grouping."""
    conda:
        op.join("..", "envs", "python.yml")
    input:
        nemo = op.join(ECKER_DATA, "nemo_meta.tsv.gz"),
        paper = op.join(ECKER_DATA, "41586_2020_3182_MOESM9_ESM.xlsx"),
    output:
        meta = op.join(ECKER_DATA, "meta.tsv.gz"),
    log:
        op.join(ECKER_DATA, "logs", "compose_metadata.log"),
    shell:
        """
        python {workflow.basedir}/scripts/compose_ecker_metadata.py \
            --nemo {input.nemo} --paper {input.paper} \
            --out {output.meta} &> {log}
        """


rule ecker_derive_tar_urls:
    """Cold-start: turn the NeMo TSV column with allc paths into tar URLs."""
    conda:
        op.join("..", "envs", "bedtools.yml")
    input:
        nemo = op.join(ECKER_DATA, "nemo_meta.tsv.gz"),
    output:
        urls = op.join(ECKER_DATA, "tar_urls.txt"),
    params:
        base = "https://data.nemoarchive.org/biccn/grant/u19_cemba/cemba/epigenome/sncell/mCseq/mouse/processed/counts/",
    log:
        op.join(ECKER_DATA, "logs", "derive_urls.log"),
    shell:
        r"""
        mkdir -p $(dirname {output.urls})
        zcat {input.nemo} \
          | cut -f29 \
          | grep -v AllcPath \
          | awk -F'/allc/' 'BEGIN{{OFS=""}} {{
                gsub(/"/, "", $0);
                n = split($1, a, "/");
                sub(/\.tsv\.gz$/, ".tsv.tar", $NF);
                print "{params.base}" a[n-1] "/" a[n] "/" $NF
            }}' \
          > {output.urls} 2> {log}
        """


checkpoint ecker_download_tars:
    """Cold-start: download per-cell tars from NeMo. Skips files already in
    ECKER_RAW so partial runs can resume."""
    conda:
        op.join("..", "envs", "bedtools.yml")
    input:
        urls = op.join(ECKER_DATA, "tar_urls.txt"),
    output:
        flag = op.join(ECKER_DATA, "download.flag"),
    params:
        raw = ECKER_RAW,
    log:
        op.join(ECKER_DATA, "logs", "download_tars.log"),
    shell:
        r"""
        mkdir -p {params.raw}
        exec &> {log}
        set -euo pipefail
        while read -r url; do
            fname=$(basename "$url")
            if [[ -f {params.raw}/$fname ]]; then
                echo "skip $fname"
            else
                echo "fetch $url"
                wget --quiet --execute=robots=off \
                     --directory-prefix={params.raw} "$url"
            fi
        done < {input.urls}
        touch {output.flag}
        """


checkpoint ecker_make_manifest:
    """Read meta.tsv.gz + the .tar files in ECKER_RAW; derive per-cell allc
    paths under ECKER_CELLS (extracted on demand by ecker_extract_tar).
    Prototype mode keeps only ecker.region_filter (default MOp) and the cell
    types in ecker.proto_cell_types, capped at prototype.cells_per_group."""
    conda:
        op.join("..", "envs", "python.yml")
    input:
        meta = op.join(ECKER_DATA, "meta.tsv.gz"),
    output:
        manifest = op.join(ECKER_DATA, "cells.tsv"),
    params:
        raw_dir = ECKER_RAW,
        cells_dir = ECKER_CELLS,
        region = config["ecker"]["region_filter"],
        proto_cell_types = ",".join(config["ecker"]["proto_cell_types"]),
        cells_per_group = config["prototype"]["cells_per_group"],
        group_col = config["ecker"]["group_column"],
        prototype = "true" if config["prototype"]["enabled"] else "false",
    log:
        op.join(ECKER_DATA, "logs", "manifest.log"),
    shell:
        """
        python {workflow.basedir}/scripts/make_manifest_ecker.py \
            --meta {input.meta} \
            --raw_dir {params.raw_dir} \
            --cells_dir {params.cells_dir} \
            --region {params.region} \
            --proto_cell_types {params.proto_cell_types} \
            --cells_per_group {params.cells_per_group} \
            --group_col {params.group_col} \
            --prototype {params.prototype} \
            --out {output.manifest} &> {log}
        """


rule ecker_extract_tar:
    """Extract a single allc tsv.gz from one tar into ECKER_CELLS/<cell_id>.tsv.gz.
    No content rewrite: we only unpack the inner gzipped tsv as-is so the allc
    parser sees the original methylpy bytes."""
    conda:
        op.join("..", "envs", "bedtools.yml")
    input:
        tar = op.join(ECKER_RAW, "{cell_id}.tsv.tar"),
    output:
        tsv = op.join(ECKER_CELLS, "{cell_id}.tsv.gz"),
    log:
        op.join(ECKER_DATA, "logs", "extract_{cell_id}.log"),
    shell:
        r"""
        mkdir -p $(dirname {output.tsv})
        # The tar contains <cell_id>/<cell_id>.tsv.gz and a sibling .idx;
        # find the .tsv.gz member and stream it as-is to disk.
        member=$(tar -tf {input.tar} | grep -E '\.tsv\.gz$' | head -1)
        if [[ -z "$member" ]]; then
            echo "no .tsv.gz inside {input.tar}" > {log}; exit 1
        fi
        tar -xOf {input.tar} "$member" > {output.tsv} 2> {log}
        """


rule ecker_make_windows_bed:
    """chr19 fixed-size windows on mm10. Output is in Ensembl naming (no chr
    prefix) to match the cells; this also mirrors what argelaguet does."""
    conda:
        op.join("..", "envs", "bedtools.yml")
    input:
        sizes = op.join(REFS, "mm10_ucsc", "chr19.sizes"),
    output:
        bed = op.join(ECKER_RUN, "beds", "windows.bed"),
    params:
        win_size = config["ecker"]["window_size"],
    log:
        op.join(ECKER_RUN, "logs", "make_windows.log"),
    shell:
        r"""
        mkdir -p $(dirname {output.bed})
        bedtools makewindows -g {input.sizes} -w {params.win_size} \
          | awk 'BEGIN{{OFS="\t"}} {{sub(/^chr/, "", $1); print $1, $2, $3, "win_"NR}}' \
          | sort -k1,1 -k2,2n > {output.bed} 2> {log}
        """


rule ecker_make_genes_bed:
    """Whole-gene loci on mm10 chr19, Ensembl chrom naming, capped to
    prototype.features_subset and stamped with feature_id. Pulls a small UCSC
    track to avoid a heavyweight GTF dependency in the prototype."""
    conda:
        op.join("..", "envs", "bedtools.yml")
    output:
        bed = op.join(ECKER_RUN, "beds", "genes.bed"),
    params:
        url = "http://hgdownload.cse.ucsc.edu/goldenpath/mm10/database/refGene.txt.gz",
        n_max = config["prototype"]["features_subset"],
    log:
        op.join(ECKER_RUN, "logs", "make_genes.log"),
    shell:
        r"""
        mkdir -p $(dirname {output.bed})
        curl -sSL {params.url} 2> {log} \
          | gunzip -c \
          | awk 'BEGIN{{OFS="\t"}} $3 == "chr19" {{print "19", $5, $6, $13}}' \
          | sort -k1,1 -k2,2n -u \
          | head -n {params.n_max} \
          | awk 'BEGIN{{OFS="\t"; k=0}} {{k++; print $1, $2, $3, "gene_" k}}' \
          > {output.bed} 2>> {log}
        """


def _ecker_cell_tsvs(wildcards):
    """Resolve per-cell tsv.gz paths from cells.tsv. Triggered through the
    manifest checkpoint so DAG expansion happens after the manifest exists."""
    import csv
    manifest = checkpoints.ecker_make_manifest.get().output.manifest
    with open(manifest) as f:
        reader = csv.DictReader(f, delimiter="\t")
        return [row["path"] for row in reader]


rule ecker_run_amet_on_scope:
    """Run amet on a {scope} BED. scope is 'windows' or 'genes'."""
    wildcard_constraints:
        scope = "windows|genes",
    conda:
        op.join("..", "envs", "bedtools.yml")
    input:
        binary = AMET,
        cells = op.join(ECKER_DATA, "cells.tsv"),
        cell_files = _ecker_cell_tsvs,
        genome = op.join(REFS, "mm10_ensembl", "chr19.fa"),
        bed = op.join(ECKER_RUN, "beds", "{scope}.bed"),
    output:
        cell_feature = op.join(ECKER_RUN, ECKER_RUN_NAME + ".{scope}.cell_feature.tsv.gz"),
        feature = op.join(ECKER_RUN, ECKER_RUN_NAME + ".{scope}.feature.tsv.gz"),
    params:
        prefix = op.join(ECKER_RUN, ECKER_RUN_NAME + ".{scope}"),
        i_max_lag = config["amet"]["i_max_lag"],
        min_cpgs = config["amet"]["min_cpgs_per_feature"],
        min_cells = config["amet"]["min_cells_per_group"],
        thresh = config["amet"]["meth_call_threshold"],
    threads: min(workflow.cores, 8)
    log:
        op.join(ECKER_RUN, "logs", "amet_{scope}.log"),
    shell:
        """
        mkdir -p $(dirname {params.prefix})
        {input.binary} \
            --genome {input.genome} \
            --cells {input.cells} \
            --features {input.bed} \
            --output-prefix {params.prefix} \
            --i-max-lag {params.i_max_lag} \
            --min-cpgs-per-feature {params.min_cpgs} \
            --min-cells-per-group {params.min_cells} \
            --meth-call-threshold {params.thresh} \
            --threads {threads} &> {log}
        """


rule ecker_render_report:
    """Single paper-style report combining feature-level and window-level views."""
    conda:
        op.join("..", "envs", "r-tools.yml")
    input:
        rmd = op.join(REPO_ROOT, "workflow", "Rmd", "ecker.Rmd"),
        feat_cell_feature = op.join(ECKER_RUN, ECKER_RUN_NAME + ".genes.cell_feature.tsv.gz"),
        feat_feature = op.join(ECKER_RUN, ECKER_RUN_NAME + ".genes.feature.tsv.gz"),
        feat_bed = op.join(ECKER_RUN, "beds", "genes.bed"),
        win_cell_feature = op.join(ECKER_RUN, ECKER_RUN_NAME + ".windows.cell_feature.tsv.gz"),
        win_feature = op.join(ECKER_RUN, ECKER_RUN_NAME + ".windows.feature.tsv.gz"),
        win_bed = op.join(ECKER_RUN, "beds", "windows.bed"),
        manifest = op.join(ECKER_DATA, "cells.tsv"),
    output:
        html = op.join(ECKER_RUN, ECKER_RUN_NAME + ".html"),
    params:
        out_dir = ECKER_RUN,
        run_name = ECKER_RUN_NAME,
    log:
        op.join(ECKER_RUN, "logs", "render.log"),
    shell:
        r"""
        mkdir -p {params.out_dir}
        Rscript -e 'rmarkdown::render("{input.rmd}",
            output_file="{params.run_name}.html",
            output_dir="{params.out_dir}",
            params=list(
                feat_cell_feature="{input.feat_cell_feature}",
                feat_feature="{input.feat_feature}",
                feat_bed="{input.feat_bed}",
                win_cell_feature="{input.win_cell_feature}",
                win_feature="{input.win_feature}",
                win_bed="{input.win_bed}",
                manifest="{input.manifest}",
                out_dir="{params.out_dir}"),
            quiet=TRUE)' &> {log}
        """
