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
  (b) download_from_nemo: cold-start path. Pulls metadata (NeMo TSV + paper
      xlsx supplement), derives the per-tar URL list, then wget-downloads
      tars into results/ecker/raw/.
"""

ECKER_DATA = op.join(RESULTS, "ecker")
ECKER_RAW = op.join(ECKER_DATA, "raw")
ECKER_CELLS = op.join(ECKER_DATA, "cells")
ECKER_RUN = op.join(RESULTS, dataset_run_name("ecker"))
ECKER_RUN_NAME = dataset_run_name("ecker")

## Ecker annotations dict. Annotation name == wildcard {annotation}; outer
## key is unused at the wildcard level.
ECKER_ANNOTATIONS = {
    "chip":      ["h3k4me3", "h3k9me3", "h3k27me3", "h3k4me1", "h3k27ac"],
    "genes":     ["genes"],
    "lines":     ["lines"],
    "sines":     ["sines"],
    "promoters": ["promoters"],
}

## All annotation BEDs come from results/ecker/mm10/ (symlinked to a shared
## mm10 tree by setup_barbara_links.sh).
ECKER_MM10 = op.join(ECKER_DATA, "mm10")
_ECKER_ALL_ANN_NAMES = sorted({a for cat in ECKER_ANNOTATIONS
                               for a in ECKER_ANNOTATIONS[cat]})

ECKER_STRATIFY_BY = ["region", "sub_type"]


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
        region_filter = config["ecker"]["region_filter"],
        proto_cell_types = proto_csv("ecker", "proto_cell_types"),
        proto_regions = proto_csv("ecker", "proto_regions"),
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
            --region_filter {params.region_filter} \
            --proto_cell_types "{params.proto_cell_types}" \
            --proto_regions "{params.proto_regions}" \
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
    """Whole-genome windows from the Ensembl sizes so contig naming (1..19,
    X, Y, MT) matches the FASTA amet scores against."""
    conda:
        op.join("..", "envs", "bedtools.yml")
    input:
        sizes = op.join(REFS, "mm10_ensembl", "genome.sizes"),
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
          | awk 'BEGIN{{OFS="\t"}}
                 {{sub(/^chr/, "", $1);
                   if ($1 ~ /^(X|Y|M|MT)$/) next;
                   print $1, $2, $3, "win_"NR}}' \
          | sort -k1,1 -k2,2n > {output.bed} 2> {log}
        """


rule ecker_stage_annotation_bed:
    """Gunzip the source <annotation>.bed.gz from results/ecker/mm10/, strip
    'chr' prefix, drop MT/X/Y, and stamp each interval with
    feature_id = <annotation>_<index>."""
    conda:
        op.join("..", "envs", "bedtools.yml")
    wildcard_constraints:
        annotation = "|".join(_ECKER_ALL_ANN_NAMES),
    input:
        bed = op.join(ECKER_MM10, "{annotation}.bed.gz"),
    output:
        bed = op.join(ECKER_RUN, "beds", "{annotation}.bed"),
    log:
        op.join(ECKER_RUN, "logs", "stage_bed_{annotation}.log"),
    shell:
        r"""
        mkdir -p $(dirname {output.bed})
        echo "[stage_bed] {wildcards.annotation}" > {log}
        zcat {input.bed} \
          | awk -v ann={wildcards.annotation} '
                 BEGIN{{OFS="\t"; k=0}}
                 {{ chr=$1; sub(/^chr/, "", chr);
                    if (chr ~ /^(X|Y|M|MT)$/) next;
                    k++;
                    print chr, $2, $3, ann "_" k }}' \
          | sort -k1,1 -k2,2n > {output.bed} 2>> {log}
        echo "[stage_bed] kept $(wc -l < {output.bed}) intervals" >> {log}
        """


rule ecker_window_annotation_per_annotation:
    """Fraction of each window covered by one annotation's intervals.
    Produces a single-column file with header == {annotation} so columns
    can be paste-merged downstream."""
    wildcard_constraints:
        annotation = "|".join(_ECKER_ALL_ANN_NAMES),
    conda:
        op.join("..", "envs", "bedtools.yml")
    input:
        windows = op.join(ECKER_RUN, "beds", "windows.bed"),
        annotation = op.join(ECKER_RUN, "beds", "{annotation}.bed"),
    output:
        frac = temp(op.join(ECKER_RUN, "beds",
                            "windows_annotation_{annotation}.frac")),
    log:
        op.join(ECKER_RUN, "logs", "window_annotation_{annotation}.log"),
    shell:
        r"""
        mkdir -p $(dirname {output.frac})
        {{
          echo "{wildcards.annotation}"
          bedtools coverage -a {input.windows} -b {input.annotation} | cut -f7
        }} > {output.frac} 2> {log}
        """


rule ecker_combine_window_annotations:
    """Paste the windows BED (chrom/start/end/feature_id) with per-annotation
    fraction columns. Output is a gzipped TSV with one row per window plus a
    single header row."""
    conda:
        op.join("..", "envs", "bedtools.yml")
    input:
        windows = op.join(ECKER_RUN, "beds", "windows.bed"),
        fracs = expand(
            op.join(ECKER_RUN, "beds",
                    "windows_annotation_{annotation}.frac"),
            annotation = _ECKER_ALL_ANN_NAMES),
    output:
        tsv = op.join(ECKER_RUN, "beds", "windows_annotation.tsv.gz"),
    log:
        op.join(ECKER_RUN, "logs", "combine_window_annotations.log"),
    shell:
        r"""
        mkdir -p $(dirname {output.tsv})
        tmp=$(mktemp)
        {{
          echo -e "chrom\tstart\tend\tfeature_id" > "$tmp"
          cat "$tmp" {input.windows} \
            | paste - {input.fracs} \
            | gzip -c > {output.tsv}
        }} 2> {log}
        rm -f "$tmp"
        """


def _ecker_all_cell_tsvs(wildcards):
    """All per-cell tsv.gz paths from the manifest checkpoint."""
    import csv
    manifest = checkpoints.ecker_make_manifest.get().output.manifest
    with open(manifest) as f:
        reader = csv.DictReader(f, delimiter="\t")
        return [row["path"] for row in reader]


rule ecker_per_combo_manifest:
    """Sub-manifest for one (region, sub_type) combo."""
    conda:
        op.join("..", "envs", "python.yml")
    input:
        cells = op.join(ECKER_DATA, "cells.tsv"),
    output:
        manifest = op.join(ECKER_DATA, "manifests",
                           "{region}_{sub_type}.tsv"),
    params:
        max_cells = max_cells_per_combo(),
    log:
        op.join(ECKER_DATA, "logs",
                "manifest_{region}_{sub_type}.log"),
    shell:
        """
        python {workflow.basedir}/scripts/ecker_subset_manifest.py \
            --cells {input.cells} \
            --region {wildcards.region} \
            --sub-type {wildcards.sub_type} \
            --max-cells {params.max_cells} \
            --out {output.manifest} &> {log}
        """


def _ecker_combo_cell_tsvs(wildcards):
    """Per-cell tsv.gz paths for one (region, sub_type) combo. Reads the
    per-combo sub-manifest produced by ecker_per_combo_manifest."""
    import csv
    sub_path = op.join(ECKER_DATA, "manifests",
                       f"{wildcards.region}_{wildcards.sub_type}.tsv")
    if not op.exists(sub_path):
        ## Snakemake hasn't built it yet; cell_files dependency is satisfied
        ## at execution time via the manifest input dependency.
        return []
    with open(sub_path) as f:
        reader = csv.DictReader(f, delimiter="\t")
        return [row["path"] for row in reader]


rule run_amet_on_ecker_features:
    """Run amet on one (annotation, region, sub_type) combo."""
    wildcard_constraints:
        annotation = "|".join(_ECKER_ALL_ANN_NAMES),
    conda:
        op.join("..", "envs", "bedtools.yml")
    input:
        binary = AMET,
        cells = op.join(ECKER_DATA, "manifests",
                        "{region}_{sub_type}.tsv"),
        cell_files = _ecker_combo_cell_tsvs,
        genome = op.join(REFS, "mm10_ensembl", "genome.fa"),
        cpg = op.join(REFS, "mm10_ensembl", "genome.fa.cpg"),
        bed = op.join(ECKER_RUN, "beds", "{annotation}.bed"),
    output:
        cell_feature = op.join(
            ECKER_RUN, "features",
            "{annotation}_{region}_{sub_type}.cell_feature.tsv.gz"),
        feature = op.join(
            ECKER_RUN, "features",
            "{annotation}_{region}_{sub_type}.feature.tsv.gz"),
    params:
        prefix = op.join(
            ECKER_RUN, "features",
            "{annotation}_{region}_{sub_type}"),
        i_max_lag = config["amet"]["i_max_lag"],
        min_cpgs = config["amet"]["min_cpgs_per_feature"],
        min_cells = min_cells_per_group(),
        thresh = config["amet"]["meth_call_threshold"],
    threads: min(workflow.cores, 4)
    log:
        op.join(ECKER_RUN, "logs",
                "amet_{annotation}_{region}_{sub_type}.log"),
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


rule run_amet_on_ecker_windows:
    """Run amet on whole-genome windows over all cells (no per-stratum wildcards)."""
    conda:
        op.join("..", "envs", "bedtools.yml")
    input:
        binary = AMET,
        cells = op.join(ECKER_DATA, "cells.tsv"),
        cell_files = _ecker_all_cell_tsvs,
        genome = op.join(REFS, "mm10_ensembl", "genome.fa"),
        cpg = op.join(REFS, "mm10_ensembl", "genome.fa.cpg"),
        bed = op.join(ECKER_RUN, "beds", "windows.bed"),
    output:
        cell_feature = op.join(
            ECKER_RUN, "windows", "all.cell_feature.tsv.gz"),
        feature = op.join(
            ECKER_RUN, "windows", "all.feature.tsv.gz"),
    params:
        prefix = op.join(ECKER_RUN, "windows", "all"),
        i_max_lag = config["amet"]["i_max_lag"],
        min_cpgs = config["amet"]["min_cpgs_per_feature"],
        min_cells = min_cells_per_group(),
        thresh = config["amet"]["meth_call_threshold"],
    threads: min(workflow.cores, 8)
    log:
        op.join(ECKER_RUN, "logs", "amet_windows.log"),
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


def _ecker_combos():
    """(region, sub_type) pairs from cells.tsv after the manifest checkpoint.
    Sanitizes both fields by replacing space with '-'."""
    import csv
    manifest_path = checkpoints.ecker_make_manifest.get().output.manifest
    pairs = set()
    with open(manifest_path) as f:
        for row in csv.DictReader(f, delimiter="\t"):
            sr = row.get("region")
            st = row.get("sub_type")
            if sr and st:
                pairs.add((str(sr).replace(" ", "-"),
                           str(st).replace(" ", "-")))
    return sorted(pairs)


def list_ecker_features_outputs(wildcards):
    combos = _ecker_combos()
    out = []
    for ann in _ECKER_ALL_ANN_NAMES:
        for sr, st in combos:
            out.append(op.join(ECKER_RUN, "features",
                               f"{ann}_{sr}_{st}.cell_feature.tsv.gz"))
            out.append(op.join(ECKER_RUN, "features",
                               f"{ann}_{sr}_{st}.feature.tsv.gz"))
    return out


def _ecker_render_shell(with_windows_annotation = False):
    helpers = op.join(REPO_ROOT, "workflow", "scripts", "render_logging.R")
    i_max_lag = config["amet"]["i_max_lag"]
    ann_line = (
        '                windows_annotation="{input.windows_annotation}",\n'
        if with_windows_annotation else ""
    )
    return rf"""
        mkdir -p {{params.out_dir}}
        export AMET_RENDER_HELPERS="{helpers}"
        Rscript -e 'rmarkdown::render("{{input.rmd}}",
            output_file="{{params.rmd_name}}.html",
            output_dir="{{params.out_dir}}",
            knit_root_dir="{{params.out_dir}}",
            params=list(
                features_dir="{{params.features_dir}}",
                win_cell_feature="{{input.win_cell_feature}}",
                win_feature="{{input.win_feature}}",
                win_bed="{{input.win_bed}}",
{ann_line}                manifest="{{input.manifest}}",
                out_dir="{{params.out_dir}}",
                log_path="{{log}}",
                threads={{threads}},
                i_max_lag={i_max_lag}),
            quiet=TRUE)' &> {{log}}
        """


## The three analytical Ecker Rmds are independent (no cross-Rmd RDS chain);
## fig_ecker.Rmd consumes their RDS/CSV intermediates. RDS/CSV files are
## declared as snakemake outputs/inputs so the graph captures the wiring.


rule render_ecker:
    conda:
        op.join("..", "envs", "r-tools.yml")
    input:
        rmd = op.join(REPO_ROOT, "workflow", "Rmd", "ecker.Rmd"),
        scripts = RMD_SHARED_SCRIPTS + [EMBEDDING_UTILS_R, DRIVER_UTILS_R],
        features = list_ecker_features_outputs,
        win_cell_feature = op.join(ECKER_RUN, "windows", "all.cell_feature.tsv.gz"),
        win_feature = op.join(ECKER_RUN, "windows", "all.feature.tsv.gz"),
        win_bed = op.join(ECKER_RUN, "beds", "windows.bed"),
        manifest = op.join(ECKER_DATA, "cells.tsv"),
    output:
        html = op.join(ECKER_RUN, "ecker.html"),
        entropy = op.join(ECKER_RUN, "ecker_entropy.rds"),
        groups_meta = op.join(ECKER_RUN, "ecker_groups_meta.rds"),
        cell_matrices = op.join(ECKER_RUN, "ecker_cell_matrices.rds"),
        umap_cell_adjS = op.join(ECKER_RUN, "ecker_umap_cell_i_total.rds"),
        umap_cell_meth = op.join(ECKER_RUN, "ecker_umap_cell_meth.rds"),
        umap_grp_jsd = op.join(ECKER_RUN, "ecker_umap_grp_jsd.rds"),
    params:
        rmd_name = "ecker",
        out_dir = ECKER_RUN,
        features_dir = op.join(ECKER_RUN, "features"),
    log:
        op.join(ECKER_RUN, "logs", "render_ecker.log"),
    threads: 4
    shell:
        _ecker_render_shell()


rule render_ecker_windows:
    conda:
        op.join("..", "envs", "r-tools.yml")
    input:
        rmd = op.join(REPO_ROOT, "workflow", "Rmd", "ecker_windows.Rmd"),
        scripts = RMD_SHARED_SCRIPTS,
        win_cell_feature = op.join(ECKER_RUN, "windows", "all.cell_feature.tsv.gz"),
        win_feature = op.join(ECKER_RUN, "windows", "all.feature.tsv.gz"),
        win_bed = op.join(ECKER_RUN, "beds", "windows.bed"),
        windows_annotation = op.join(ECKER_RUN, "beds", "windows_annotation.tsv.gz"),
        manifest = op.join(ECKER_DATA, "cells.tsv"),
    output:
        html = op.join(ECKER_RUN, "ecker_windows.html"),
        per_cell_summary = op.join(ECKER_RUN, "ecker_windows_per_cell_summary.csv"),
    params:
        rmd_name = "ecker_windows",
        out_dir = ECKER_RUN,
        features_dir = op.join(ECKER_RUN, "features"),
    log:
        op.join(ECKER_RUN, "logs", "render_ecker_windows.log"),
    threads: 4
    shell:
        _ecker_render_shell(with_windows_annotation = True)


rule render_ecker_embeddings:
    conda:
        op.join("..", "envs", "r-tools.yml")
    input:
        rmd = op.join(REPO_ROOT, "workflow", "Rmd", "ecker_embeddings.Rmd"),
        scripts = RMD_SHARED_SCRIPTS + [EMBEDDING_UTILS_R],
        win_cell_feature = op.join(ECKER_RUN, "windows", "all.cell_feature.tsv.gz"),
        win_feature = op.join(ECKER_RUN, "windows", "all.feature.tsv.gz"),
        win_bed = op.join(ECKER_RUN, "beds", "windows.bed"),
        windows_annotation = op.join(ECKER_RUN, "beds", "windows_annotation.tsv.gz"),
        manifest = op.join(ECKER_DATA, "cells.tsv"),
    output:
        html = op.join(ECKER_RUN, "ecker_embeddings.html"),
        umap_windows = op.join(ECKER_RUN, "ecker_umap_windows_i_total.rds"),
        per_cell_summary = op.join(ECKER_RUN, "ecker_embeddings_per_cell_summary.csv"),
        win_varexp = op.join(ECKER_RUN, "ecker_win_varexp.csv"),
        diagnostics = op.join(ECKER_RUN, "ecker_embedding_diagnostics.csv"),
    params:
        rmd_name = "ecker_embeddings",
        out_dir = ECKER_RUN,
        features_dir = op.join(ECKER_RUN, "features"),
    log:
        op.join(ECKER_RUN, "logs", "render_ecker_embeddings.log"),
    threads: 4
    shell:
        _ecker_render_shell(with_windows_annotation = True)


rule render_fig_ecker_rmd:
    """Render fig_ecker.Rmd; consumes RDS/CSV intermediates from the three
    analytical rules above."""
    conda:
        op.join("..", "envs", "r-tools.yml")
    input:
        rmd = op.join(REPO_ROOT, "workflow", "Rmd", "fig_ecker.Rmd"),
        scripts = RMD_SHARED_SCRIPTS + [DRIVER_UTILS_R],
        entropy = op.join(ECKER_RUN, "ecker_entropy.rds"),
        groups_meta = op.join(ECKER_RUN, "ecker_groups_meta.rds"),
        umap_windows = op.join(ECKER_RUN, "ecker_umap_windows_i_total.rds"),
        win_varexp = op.join(ECKER_RUN, "ecker_win_varexp.csv"),
        per_cell_summary = op.join(ECKER_RUN, "ecker_embeddings_per_cell_summary.csv"),
        diagnostics = op.join(ECKER_RUN, "ecker_embedding_diagnostics.csv"),
        win_cell_feature = op.join(ECKER_RUN, "windows", "all.cell_feature.tsv.gz"),
        win_feature = op.join(ECKER_RUN, "windows", "all.feature.tsv.gz"),
        win_bed = op.join(ECKER_RUN, "beds", "windows.bed"),
        windows_annotation = op.join(ECKER_RUN, "beds", "windows_annotation.tsv.gz"),
        manifest = op.join(ECKER_DATA, "cells.tsv"),
    output:
        html = op.join(ECKER_RUN, "fig_ecker.html"),
    params:
        rmd_name = "fig_ecker",
        out_dir = ECKER_RUN,
        features_dir = op.join(ECKER_RUN, "features"),
    log:
        op.join(ECKER_RUN, "logs", "render_fig_ecker.log"),
    threads: 4
    shell:
        _ecker_render_shell(with_windows_annotation = True)
