"""
Argelaguet et al. 2019 mouse gastrulation scNMT-seq: methylation entropy by feature.

Multi-omics profiling of mouse gastrulation at single cell resolution
Ricard Argelaguet, Stephen J Clark, Hisham Mohammed, et al.
https://pmc.ncbi.nlm.nih.gov/articles/PMC6924995
https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE121708

Assembly: GRCm38 (mm10). Reads native scNMT cpg_level files via amet's scnmt
parser; no harmonization step. amet derives the CpG reference internally
from the FASTA (--genome) and resolves chr/no-chr prefix mismatches in its
CpgReference, so no chrom-name munging happens at the workflow level.

Path convention:
  results/argelaguet/        raw inputs pulled from barbara (shared across runs)
  results/<run_name>/        amet outputs + paper-figure report for one run
"""

ARG_DATA = op.join(RESULTS, "argelaguet")
ARG_CELLS = op.join(ARG_DATA, "cells")
ARG_FEATURES_DIR = op.join(ARG_DATA, "features")
ARG_MM10_DIR = op.join(ARG_DATA, "mm10")

ARG_RUN_NAME = dataset_run_name("argelaguet")
ARG_RUN = op.join(RESULTS, ARG_RUN_NAME)

## Annotation set: outer key = category, inner key = annotation name; the
## annotation name is the {annotation} wildcard.
_ARGELAGUET_MM10_ANNOTATIONS = {
    "chip":      ["h3k4me3", "h3k9me3", "h3k27me3", "h3k4me1", "h3k27ac"],
    "genes":     ["genes"],
    "lines":     ["lines"],
    "sines":     ["sines"],
    "promoters": ["promoters"],
}

## Gastrulation-specific source filenames inside the scnmt_gastrulation
## features tarball. Keys are sanitized annotation wildcard values.
_ARGELAGUET_GASTRO_BEDS = {
    "enh-E75-Ect":        "H3K27ac_distal_E7.5_Ect_intersect12.bed",
    "enh-E75-End":        "H3K27ac_distal_E7.5_End_intersect12.bed",
    "enh-E75-Mes":        "H3K27ac_distal_E7.5_Mes_intersect12.bed",
    "enh-E75-union":      "H3K27ac_distal_E7.5_union_intersect12.bed",
    "h3k4me3-E75-Ect":    "H3K4me3_E7.5_Ect.bed",
    "h3k4me3-E75-End":    "H3K4me3_E7.5_End.bed",
    "h3k4me3-E75-Mes":    "H3K4me3_E7.5_Mes.bed",
    "h3k4me3-E75-common": "H3K4me3_E7.5_common.bed",
    "esc-p300":           "ESC_p300.bed",
    "esc-dhs":            "ESC_DHS.bed",
}

_ARGELAGUET_GASTRO_ANNOTATIONS = {
    "enh_gastro":  [k for k in _ARGELAGUET_GASTRO_BEDS if k.startswith("enh-")],
    "h3k4me3_E75": [k for k in _ARGELAGUET_GASTRO_BEDS if k.startswith("h3k4me3-")],
    "esc":         [k for k in _ARGELAGUET_GASTRO_BEDS if k.startswith("esc-")],
}

ARGELAGUET_ANNOTATIONS = {
    **_ARGELAGUET_MM10_ANNOTATIONS,
    **_ARGELAGUET_GASTRO_ANNOTATIONS,
}

_MM10_ANN_NAMES = [a for cat in _ARGELAGUET_MM10_ANNOTATIONS
                   for a in _ARGELAGUET_MM10_ANNOTATIONS[cat]]

_ALL_ARGELAGUET_ANN_NAMES = sorted(
    set(_MM10_ANN_NAMES) | set(_ARGELAGUET_GASTRO_BEDS.keys())
)

ARGELAGUET_STRATIFY_BY = ["stage", "lineage"]
ARGELAGUET_MAX_CELLS = max_cells_per_combo()


rule filter_argelaguet_metadata:
    """Drop QC-failing cells, drop TET TKO plates, keep stage/lineage/embryo."""
    conda:
        op.join("..", "envs", "python.yml")
    input:
        meta = op.join(ARG_DATA, "sample_metadata.txt"),
    output:
        meta = op.join(ARG_DATA, "meta.tsv.gz"),
    log:
        op.join(ARG_DATA, "logs", "filter_metadata.log"),
    shell:
        """
        python {workflow.basedir}/scripts/filter_argelaguet_metadata.py \
            --meta_in {input.meta} \
            --meta_out {output.meta} &> {log}
        """


checkpoint make_argelaguet_manifest:
    """Build amet cells.tsv from the filtered metadata + local cell files.
    Checkpoint so the per-(stage, lineage) DAG fans out from the manifest."""
    conda:
        op.join("..", "envs", "r-tools.yml")
    input:
        meta = op.join(ARG_DATA, "meta.tsv.gz"),
    output:
        manifest = op.join(ARG_DATA, "cells.tsv"),
    params:
        cells_dir = ARG_CELLS,
        group_col = config["argelaguet"]["group_column"],
        proto_stages = proto_csv("argelaguet", "proto_stages"),
        proto_lineages = proto_csv("argelaguet", "proto_lineages"),
        prototype = "true" if config["prototype"]["enabled"] else "false",
    log:
        op.join(ARG_DATA, "logs", "manifest.log"),
    shell:
        """
        Rscript {workflow.basedir}/scripts/make_manifest_argelaguet.R \
            --metadata {input.meta} \
            --cells_dir {params.cells_dir} \
            --group_col {params.group_col} \
            --proto_stages "{params.proto_stages}" \
            --proto_lineages "{params.proto_lineages}" \
            --prototype {params.prototype} \
            --out {output.manifest} &> {log}
        """


rule argelaguet_per_combo_manifest:
    """Subset cells.tsv to one (sanitized stage, sanitized lineage) pair.
    Plate-stratified top-N-by-coverage selection."""
    conda:
        op.join("..", "envs", "python.yml")
    input:
        cells = op.join(ARG_DATA, "cells.tsv"),
    output:
        manifest = op.join(ARG_DATA, "manifests",
                           "{stage}_{lineage}.tsv"),
    params:
        max_cells = ARGELAGUET_MAX_CELLS,
    log:
        op.join(ARG_DATA, "logs", "manifest_{stage}_{lineage}.log"),
    shell:
        """
        python {workflow.basedir}/scripts/argelaguet_subset_manifest.py \
            --cells {input.cells} \
            --stage {wildcards.stage} \
            --lineage {wildcards.lineage} \
            --max-cells {params.max_cells} \
            --out {output.manifest} &> {log}
        """


def _argelaguet_bed_source(wildcards):
    """Resolve the BED source path for an annotation wildcard.

    mm10 annotations come from results/argelaguet/mm10/<ann>.bed.gz;
    gastrulation-specific annotations come from
    results/argelaguet/features/<filename> per the _ARGELAGUET_GASTRO_BEDS
    map.
    """
    if wildcards.annotation in _MM10_ANN_NAMES:
        return op.join(ARG_MM10_DIR, f"{wildcards.annotation}.bed.gz")
    return op.join(ARG_FEATURES_DIR,
                   _ARGELAGUET_GASTRO_BEDS[wildcards.annotation])


rule argelaguet_filter_annotation_bed:
    """Stage one annotation BED: gunzip if needed, strip 'chr' prefix, and
    stamp each peak with feature_id = <annotation>_<index>. Whole-genome."""
    conda:
        op.join("..", "envs", "bedtools.yml")
    wildcard_constraints:
        annotation = "|".join(_ALL_ARGELAGUET_ANN_NAMES),
    input:
        bed = _argelaguet_bed_source,
    output:
        bed = op.join(ARG_RUN, "beds", "{annotation}.bed"),
    log:
        op.join(ARG_RUN, "logs", "filter_bed_{annotation}.log"),
    shell:
        r"""
        mkdir -p $(dirname {output.bed})
        echo "[filter_bed] {wildcards.annotation}: whole-genome" > {log}
        if [[ "{input.bed}" == *.gz ]]; then
            CAT="zcat"
        else
            CAT="cat"
        fi
        $CAT {input.bed} | awk -v ann={wildcards.annotation} '
             BEGIN{{OFS="\t"; k=0}}
             {{ chr=$1; sub(/^chr/, "", chr);
                if (chr ~ /^(X|Y|M|MT)$/) next;
                k++;
                print chr, $2, $3, ann "_" k }}' > {output.bed}
        echo "[filter_bed] kept $(wc -l < {output.bed}) intervals" >> {log}
        """


rule argelaguet_make_windows:
    """Tiled whole-genome windows for the windows/embeddings analysis.
    Strips 'chr' to match Ensembl-named cells."""
    conda:
        op.join("..", "envs", "bedtools.yml")
    input:
        sizes = op.join(REFS, "mm10_ucsc", "genome.sizes"),
    output:
        bed = op.join(ARG_RUN, "beds", "windows.bed"),
    params:
        win_size = config["argelaguet"]["window_size"],
    shell:
        r"""
        bedtools makewindows -g {input.sizes} -w {params.win_size} \
          | awk 'BEGIN{{OFS="\t"}}
                 {{sub(/^chr/,"",$1);
                   if ($1 ~ /^(X|Y|M|MT)$/) next;
                   print $1, $2, $3, "win_"NR}}' \
          | sort -k1,1 -k2,2n > {output.bed}
        """


rule chr19_sizes:
    """Chrom-sizes from the chr19 FASTA. bedtools makewindows needs it; emitted
    here so per-dataset rules can use it without duplicating awk."""
    conda:
        op.join("..", "envs", "bedtools.yml")
    input:
        fa = op.join(REFS, "{source}", "chr19.fa"),
    output:
        sizes = op.join(REFS, "{source}", "chr19.sizes"),
    shell:
        r"""
        awk 'BEGIN{{n=""; len=0}}
             /^>/{{if(n) print n"\t"len; n=substr($1,2); len=0; next}}
             {{len+=length($0)}}
             END{{if(n) print n"\t"len}}' {input.fa} > {output.sizes}
        """


rule run_amet_on_argelaguet_features:
    """Run amet on one (annotation, stage, lineage) combo. Wildcards:
    {annotation, stage, lineage}, where stage and lineage are sanitized
    strings (gsub '[ ._]' '-')."""
    wildcard_constraints:
        annotation = "|".join(_ALL_ARGELAGUET_ANN_NAMES),
    conda:
        op.join("..", "envs", "bedtools.yml")
    input:
        binary = AMET,
        cells = op.join(ARG_DATA, "manifests", "{stage}_{lineage}.tsv"),
        genome = op.join(REFS, "mm10_ucsc", "genome.fa"),
        cpg = op.join(REFS, "mm10_ucsc", "genome.fa.cpg"),
        bed = op.join(ARG_RUN, "beds", "{annotation}.bed"),
    output:
        cell_feature = op.join(
            ARG_RUN, "features",
            "{annotation}_{stage}_{lineage}.cell_feature.tsv.gz"),
        feature = op.join(
            ARG_RUN, "features",
            "{annotation}_{stage}_{lineage}.feature.tsv.gz"),
    params:
        prefix = op.join(
            ARG_RUN, "features",
            "{annotation}_{stage}_{lineage}"),
        i_max_lag = config["amet"]["i_max_lag"],
        min_cpgs = config["amet"]["min_cpgs_per_feature"],
        min_cells = min_cells_per_group(),
        thresh = config["amet"]["meth_call_threshold"],
    threads: min(workflow.cores, 4)
    log:
        op.join(ARG_RUN, "logs",
                "amet_{annotation}_{stage}_{lineage}.log"),
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


rule argelaguet_window_annotation_per_annotation:
    """Per-window fractional coverage by one annotation. bedtools coverage's
    7th column is the fraction of the window covered by features in the
    annotation BED. Header line carries the annotation name so the combine
    step can paste columns by position."""
    conda:
        op.join("..", "envs", "bedtools.yml")
    wildcard_constraints:
        annotation = "|".join(_ALL_ARGELAGUET_ANN_NAMES),
    input:
        windows = op.join(ARG_RUN, "beds", "windows.bed"),
        annotation = op.join(ARG_RUN, "beds", "{annotation}.bed"),
    output:
        frac = temp(op.join(ARG_RUN, "beds", "annotation_cov",
                            "{annotation}.frac")),
    log:
        op.join(ARG_RUN, "logs", "window_annotation_{annotation}.log"),
    shell:
        r"""
        mkdir -p $(dirname {output.frac})
        {{
          echo "{wildcards.annotation}"
          bedtools coverage -a {input.windows} -b {input.annotation} | cut -f7
        }} > {output.frac} 2> {log}
        """


rule argelaguet_combine_window_annotations:
    """Stitch the per-annotation fractional-coverage columns onto the windows
    BED. Output is a header-tagged TSV: chrom, start, end, feature_id, then
    one column per annotation. Drives the per-window annotation matrix used
    by the Rmds that colour windows by genomic context."""
    conda:
        op.join("..", "envs", "bedtools.yml")
    input:
        windows = op.join(ARG_RUN, "beds", "windows.bed"),
        fracs = expand(
            op.join(ARG_RUN, "beds", "annotation_cov", "{annotation}.frac"),
            annotation = _ALL_ARGELAGUET_ANN_NAMES,
        ),
    output:
        annotation = op.join(ARG_RUN, "beds", "windows_annotation.tsv.gz"),
    log:
        op.join(ARG_RUN, "logs", "combine_window_annotations.log"),
    shell:
        r"""
        mkdir -p $(dirname {output.annotation})
        tmp_header=$(mktemp)
        tmp_bed=$(mktemp)
        {{
          printf "chrom\tstart\tend\tfeature_id\n" > $tmp_header
          cat $tmp_header {input.windows} > $tmp_bed
          paste $tmp_bed {input.fracs} | gzip -c > {output.annotation}
          echo "[combine] wrote $(zcat {output.annotation} | wc -l) rows"
        }} > {log} 2>&1
        rm -f $tmp_header $tmp_bed
        """


rule run_amet_on_argelaguet_windows:
    """Run amet over all cells on chr19 windows: one big run, all cells,
    no stratum wildcard."""
    conda:
        op.join("..", "envs", "bedtools.yml")
    input:
        binary = AMET,
        cells = op.join(ARG_DATA, "cells.tsv"),
        genome = op.join(REFS, "mm10_ucsc", "genome.fa"),
        cpg = op.join(REFS, "mm10_ucsc", "genome.fa.cpg"),
        bed = op.join(ARG_RUN, "beds", "windows.bed"),
    output:
        cell_feature = op.join(
            ARG_RUN, "windows", "all.cell_feature.tsv.gz"),
        feature = op.join(
            ARG_RUN, "windows", "all.feature.tsv.gz"),
    params:
        prefix = op.join(ARG_RUN, "windows", "all"),
        i_max_lag = config["amet"]["i_max_lag"],
        min_cpgs = config["amet"]["min_cpgs_per_feature"],
        min_cells = min_cells_per_group(),
        thresh = config["amet"]["meth_call_threshold"],
    threads: min(workflow.cores, 8)
    log:
        op.join(ARG_RUN, "logs", "amet_windows.log"),
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


def _argelaguet_combos():
    """Read cells.tsv after the manifest checkpoint and return sorted unique
    (sanitized stage, sanitized lineage) pairs that have at least one cell."""
    import csv
    import re

    def san(x):
        return re.sub(r"[ ._]", "-", str(x))

    manifest_path = checkpoints.make_argelaguet_manifest.get().output.manifest
    pairs = set()
    with open(manifest_path) as f:
        for row in csv.DictReader(f, delimiter="\t"):
            stage = row.get("stage")
            lineage = row.get("lineage10x")
            if stage and lineage:
                pairs.add((san(stage), san(lineage)))
    return sorted(pairs)


def list_argelaguet_features_outputs(wildcards):
    """All (annotation x stage x lineage) amet output files."""
    combos = _argelaguet_combos()
    out = []
    for ann in _ALL_ARGELAGUET_ANN_NAMES:
        for stage, lineage in combos:
            out.append(op.join(ARG_RUN, "features",
                               f"{ann}_{stage}_{lineage}.cell_feature.tsv.gz"))
            out.append(op.join(ARG_RUN, "features",
                               f"{ann}_{stage}_{lineage}.feature.tsv.gz"))
    return out


def _argelaguet_render_shell(with_windows_annotation = False):
    """Shell template for an Argelaguet Rmd render. windows_annotation is
    optional because the per-feature Rmd does not need a per-window matrix."""
    helpers = op.join(REPO_ROOT, "workflow", "scripts", "render_logging.R")
    i_max_lag = config["amet"]["i_max_lag"]
    extra = ''
    if with_windows_annotation:
        extra = ',\n                windows_annotation="{input.windows_annotation}"'
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
                manifest="{{input.manifest}}",
                out_dir="{{params.out_dir}}",
                log_path="{{log}}",
                threads={{threads}},
                i_max_lag={i_max_lag}{extra}),
            quiet=TRUE)' &> {{log}}
        """


## The three analytical Argelaguet Rmds are independent (no cross-Rmd RDS
## chain); fig_argelaguet.Rmd consumes their RDS/CSV intermediates. RDS/CSV
## files are declared as snakemake outputs/inputs so the graph captures the
## wiring.


rule render_argelaguet:
    conda:
        op.join("..", "envs", "r-tools.yml")
    input:
        rmd = op.join(REPO_ROOT, "workflow", "Rmd", "argelaguet.Rmd"),
        scripts = RMD_SHARED_SCRIPTS + [DRIVER_UTILS_R, EMBEDDING_UTILS_R],
        features = list_argelaguet_features_outputs,
        win_cell_feature = op.join(ARG_RUN, "windows", "all.cell_feature.tsv.gz"),
        win_feature = op.join(ARG_RUN, "windows", "all.feature.tsv.gz"),
        win_bed = op.join(ARG_RUN, "beds", "windows.bed"),
        manifest = op.join(ARG_DATA, "cells.tsv"),
    output:
        html = op.join(ARG_RUN, "argelaguet.html"),
        entropy = op.join(ARG_RUN, "argelaguet_entropy.rds"),
        groups_meta = op.join(ARG_RUN, "argelaguet_groups_meta.rds"),
        cell_matrices = op.join(ARG_RUN, "argelaguet_cell_matrices.rds"),
        umap_cell_adjS = op.join(ARG_RUN, "argelaguet_umap_cell_i_total.rds"),
        umap_cell_meth = op.join(ARG_RUN, "argelaguet_umap_cell_meth.rds"),
        umap_grp_jsd = op.join(ARG_RUN, "argelaguet_umap_grp_jsd.rds"),
    params:
        rmd_name = "argelaguet",
        out_dir = ARG_RUN,
        features_dir = op.join(ARG_RUN, "features"),
    log:
        op.join(ARG_RUN, "logs", "render_argelaguet.log"),
    threads: 4
    shell:
        _argelaguet_render_shell()


rule render_argelaguet_windows:
    conda:
        op.join("..", "envs", "r-tools.yml")
    input:
        rmd = op.join(REPO_ROOT, "workflow", "Rmd", "argelaguet_windows.Rmd"),
        scripts = RMD_SHARED_SCRIPTS,
        win_cell_feature = op.join(ARG_RUN, "windows", "all.cell_feature.tsv.gz"),
        win_feature = op.join(ARG_RUN, "windows", "all.feature.tsv.gz"),
        win_bed = op.join(ARG_RUN, "beds", "windows.bed"),
        windows_annotation = op.join(ARG_RUN, "beds", "windows_annotation.tsv.gz"),
        manifest = op.join(ARG_DATA, "cells.tsv"),
    output:
        html = op.join(ARG_RUN, "argelaguet_windows.html"),
        per_cell_summary = op.join(ARG_RUN, "argelaguet_windows_per_cell_summary.csv"),
    params:
        rmd_name = "argelaguet_windows",
        out_dir = ARG_RUN,
        features_dir = op.join(ARG_RUN, "features"),
    log:
        op.join(ARG_RUN, "logs", "render_argelaguet_windows.log"),
    threads: 4
    shell:
        _argelaguet_render_shell(with_windows_annotation = True)


rule render_argelaguet_embeddings:
    conda:
        op.join("..", "envs", "r-tools.yml")
    input:
        rmd = op.join(REPO_ROOT, "workflow", "Rmd", "argelaguet_embeddings.Rmd"),
        scripts = RMD_SHARED_SCRIPTS + [EMBEDDING_UTILS_R],
        win_cell_feature = op.join(ARG_RUN, "windows", "all.cell_feature.tsv.gz"),
        win_feature = op.join(ARG_RUN, "windows", "all.feature.tsv.gz"),
        win_bed = op.join(ARG_RUN, "beds", "windows.bed"),
        windows_annotation = op.join(ARG_RUN, "beds", "windows_annotation.tsv.gz"),
        manifest = op.join(ARG_DATA, "cells.tsv"),
    output:
        html = op.join(ARG_RUN, "argelaguet_embeddings.html"),
        umap_windows = op.join(ARG_RUN, "argelaguet_umap_windows_i_total.rds"),
        per_cell_summary = op.join(ARG_RUN, "argelaguet_embeddings_per_cell_summary.csv"),
        win_varexp = op.join(ARG_RUN, "argelaguet_win_varexp.csv"),
    params:
        rmd_name = "argelaguet_embeddings",
        out_dir = ARG_RUN,
        features_dir = op.join(ARG_RUN, "features"),
    log:
        op.join(ARG_RUN, "logs", "render_argelaguet_embeddings.log"),
    threads: 4
    shell:
        _argelaguet_render_shell(with_windows_annotation = True)


rule render_fig_argelaguet_rmd:
    """Render fig_argelaguet.Rmd; consumes RDS/CSV intermediates from the
    three analytical rules above."""
    conda:
        op.join("..", "envs", "r-tools.yml")
    input:
        rmd = op.join(REPO_ROOT, "workflow", "Rmd", "fig_argelaguet.Rmd"),
        scripts = RMD_SHARED_SCRIPTS + [DRIVER_UTILS_R],
        entropy = op.join(ARG_RUN, "argelaguet_entropy.rds"),
        groups_meta = op.join(ARG_RUN, "argelaguet_groups_meta.rds"),
        cell_matrices = op.join(ARG_RUN, "argelaguet_cell_matrices.rds"),
        umap_cell_adjS = op.join(ARG_RUN, "argelaguet_umap_cell_i_total.rds"),
        per_cell_summary = op.join(ARG_RUN, "argelaguet_embeddings_per_cell_summary.csv"),
        win_varexp = op.join(ARG_RUN, "argelaguet_win_varexp.csv"),
        win_cell_feature = op.join(ARG_RUN, "windows", "all.cell_feature.tsv.gz"),
        win_feature = op.join(ARG_RUN, "windows", "all.feature.tsv.gz"),
        win_bed = op.join(ARG_RUN, "beds", "windows.bed"),
        windows_annotation = op.join(ARG_RUN, "beds", "windows_annotation.tsv.gz"),
        manifest = op.join(ARG_DATA, "cells.tsv"),
    output:
        html = op.join(ARG_RUN, "fig_argelaguet.html"),
    params:
        rmd_name = "fig_argelaguet",
        out_dir = ARG_RUN,
        features_dir = op.join(ARG_RUN, "features"),
    log:
        op.join(ARG_RUN, "logs", "render_fig_argelaguet.log"),
    threads: 4
    shell:
        _argelaguet_render_shell(with_windows_annotation = True)
