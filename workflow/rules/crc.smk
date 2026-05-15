"""
Bian et al. 2018 colorectal cancer scTrioSeq2 single-cell methylation.

Single-cell multiomics sequencing and analyses of human colorectal cancer.
Bian et al. https://www.science.org/doi/10.1126/science.aao3791
GEO: GSE97693. Assembly: hg19.

Per-cell raw files are Bian-style singleC TSVs (Bismark-derived, 10 columns,
header line). amet reads them via its 'bismark' parser; no harmonization step.

Two acquisition paths:

  (a) rsync_from_barbara: rsync the raw originals already present on barbara
      into results/crc/raw/. This is the normal path while we have access to
      the existing copy.
  (b) download_from_geo: cold-start path. Pulls the GSM list from SRA
      (esearch + efetch) and per-sample directories from GEO. Used when
      starting from scratch with no barbara mirror.

Both paths populate the same results/crc/raw/ directory; the manifest rule
reads whatever is there.
"""

CRC_DATA = op.join(RESULTS, "crc")
CRC_RAW = op.join(CRC_DATA, "raw")
CRC_BEDS = op.join(CRC_DATA, "beds")
CRC_RUN = op.join(RESULTS, dataset_run_name("crc"))
CRC_RUN_NAME = dataset_run_name("crc")

## CRC annotations dict. Keys: cat (outer), values: list of subcats (inner).
## Wildcards: {subcat}_{cat}_{patient}_{location}.
##
## Locally we only have cpgIslandExt (downloaded from UCSC). The rest live on
## barbara as bed files; they get pulled in when running there. The dict keeps
## the full surface so rules referencing wildcard combos remain stable.
CRC_ANNOTATIONS = {
    "pmd":          ["pmds", "hmds"],
    "hmm": [
        "0_Enhancer", "2_Enhancer",
        "11_Promoter", "12_Promoter",
        "1_Transcribed", "4_Transcribed",
        "5_RegPermissive", "7_RegPermissive",
        "6_LowConfidence",
        "3_Quiescent", "8_Quiescent", "10_Quiescent",
        "9_ConstitutiveHet", "13_ConstitutiveHet",
    ],
    "chip":         ["H3K27me3", "H3K9me3", "H3K4me3"],
    "lad":          ["laminb1"],
    "genes":        ["genes"],
    "lines":        ["lines"],
    "sines":        ["sines"],
    "cpgIslandExt": ["cpgIslandExt"],
    "scna":         ["crc01_nc_scna", "crc01_gain_scna", "crc01_lost_scna"],
}

_CRC_LOCAL_ANNOTATIONS = CRC_ANNOTATIONS
_CRC_LOCAL_PAIRS = [(sc, c) for c, subs in _CRC_LOCAL_ANNOTATIONS.items() for sc in subs]
_CRC_ALL_PAIRS = [(sc, c) for c, subs in CRC_ANNOTATIONS.items() for sc in subs]
_CRC_SUBCAT_RE = "|".join(sorted({sc for sc, _ in _CRC_ALL_PAIRS}))
_CRC_CAT_RE = "|".join(sorted({c  for _, c  in _CRC_ALL_PAIRS}))

## hg19 BED trees. Both symlinked in by setup_barbara_links.sh:
##   hg19         (cpgIslandExt, SCNAs, genes/lines/sines.bed.gz)
##   hg19_curated (chromHMM, ChIP, lamin, PMD as plain .bed)
CRC_HG19_WF = op.join(CRC_DATA, "hg19")
CRC_HG19_CURATED = op.join(CRC_DATA, "hg19_curated")

def _crc_yamet_bed_path(subcat, cat):
    """Resolve a (subcat, cat) pair to its source file in the symlinked BED
    trees. Returned paths may be plain BED or gzipped; the staging rule
    handles both."""
    if cat in ("hmm", "chip", "lad", "pmd"):
        return op.join(CRC_HG19_CURATED, f"{subcat}.{cat}.bed")
    if cat == "scna":
        return op.join(CRC_HG19_WF, f"{subcat}.{cat}.bed")
    if cat == "cpgIslandExt":
        return op.join(CRC_HG19_WF, "cpgIslandExt.cpgIslandExt.bed")
    if cat in ("genes", "lines", "sines"):
        stem = {"genes": "genes", "lines": "rmsk.lines", "sines": "rmsk.sines"}[cat]
        return op.join(CRC_HG19_WF, f"{stem}.bed.gz")
    raise ValueError(f"no source mapping for ({subcat}, {cat})")


rule crc_download_accessors:
    """Pull the bisulfite GSM list from SRA. Only used in the cold-start path."""
    conda:
        op.join("..", "envs", "bedtools.yml")
    output:
        gsm = op.join(CRC_DATA, "bisulfites_gsm.txt"),
    log:
        op.join(CRC_DATA, "logs", "download_accessors.log"),
    shell:
        r"""
        mkdir -p $(dirname {output.gsm})
        esearch -db sra -query PRJNA382695 \
          | efetch -format runinfo \
          | cut -f11,13,14,15,30 -d"," \
          | grep Bis \
          | cut -f5 -d"," > {output.gsm} 2> {log}
        """


checkpoint crc_download_bismarks:
    """Cold-start: download per-GSM singleC tarballs from GEO. Skips files
    already present in CRC_RAW so partial downloads can resume cleanly."""
    conda:
        op.join("..", "envs", "bedtools.yml")
    input:
        gsm = op.join(CRC_DATA, "bisulfites_gsm.txt"),
    output:
        flag = op.join(CRC_DATA, "download.flag"),
    params:
        raw = CRC_RAW,
    log:
        op.join(CRC_DATA, "logs", "download_bismarks.log"),
    shell:
        r"""
        mkdir -p {params.raw}
        exec &> {log}
        set -euo pipefail
        while read -r gsm; do
            short="$(echo $gsm | cut -c1-7)"
            url="ftp://ftp.ncbi.nlm.nih.gov/geo/samples/${{short}}nnn/$gsm/suppl/"
            existing=$(find {params.raw} -maxdepth 1 -name "${{gsm}}*" -print -quit)
            if [[ -n "$existing" ]]; then
                echo "skip $gsm: $existing"
            else
                echo "fetch $gsm $url"
                wget --quiet --execute=robots=off --recursive --convert-links \
                     --accept=gz --no-directories \
                     --directory-prefix={params.raw} --timestamping "$url"
            fi
        done < {input.gsm}
        touch {output.flag}
        """


checkpoint crc_make_manifest:
    """Parse CRC_RAW filenames into a cells.tsv keyed by patient/location.
    Checkpoint so the per-(patient, location) DAG fans out from the manifest.

    Filename pattern: GSM<id>_scTrioSeq2Met_<patient>_<location><lane>_<cell>.singleC.txt.gz
    Prototype mode keeps only the patients and locations listed in
    config['crc']. Group column is 'location' (NC/PT/LN); patient is passed
    through as an extra column.
    """
    conda:
        op.join("..", "envs", "r-tools.yml")
    output:
        manifest = op.join(CRC_DATA, "cells.tsv"),
    params:
        raw_dir = CRC_RAW,
        proto_patients = proto_csv("crc", "proto_patients"),
        proto_locations = proto_csv("crc", "proto_locations"),
        cells_per_group = config["prototype"]["cells_per_group"],
        prototype = "true" if config["prototype"]["enabled"] else "false",
    log:
        op.join(CRC_DATA, "logs", "manifest.log"),
    shell:
        """
        Rscript {workflow.basedir}/scripts/make_manifest_crc.R \
            --raw_dir {params.raw_dir} \
            --proto_patients "{params.proto_patients}" \
            --proto_locations "{params.proto_locations}" \
            --cells_per_group {params.cells_per_group} \
            --prototype {params.prototype} \
            --out {output.manifest} &> {log}
        """


rule crc_per_combo_manifest:
    """Sub-manifest for one (patient, location) combo."""
    conda:
        op.join("..", "envs", "python.yml")
    input:
        cells = op.join(CRC_DATA, "cells.tsv"),
    output:
        manifest = op.join(CRC_DATA, "manifests",
                           "{patient}_{location}.tsv"),
    params:
        max_cells = max_cells_per_combo(),
    log:
        op.join(CRC_DATA, "logs", "manifest_{patient}_{location}.log"),
    shell:
        """
        python {workflow.basedir}/scripts/crc_subset_manifest.py \
            --cells {input.cells} \
            --patient {wildcards.patient} \
            --location {wildcards.location} \
            --max-cells {params.max_cells} \
            --out {output.manifest} &> {log}
        """


rule crc_make_windows_bed:
    """Whole-genome fixed-size windows on hg19 (UCSC chrom names match the singleC files)."""
    conda:
        op.join("..", "envs", "bedtools.yml")
    input:
        sizes = op.join(REFS, "hg19_ucsc", "genome.sizes"),
    output:
        bed = op.join(CRC_RUN, "beds", "windows.bed"),
    params:
        win_size = config["crc"]["window_size"],
    log:
        op.join(CRC_RUN, "logs", "make_windows.log"),
    shell:
        r"""
        mkdir -p $(dirname {output.bed})
        bedtools makewindows -g {input.sizes} -w {params.win_size} \
          | awk 'BEGIN{{OFS="\t"}}
                 {{if ($1 ~ /^chr(X|Y|M|MT)$/) next;
                   print $1, $2, $3, "win_"NR}}' \
          | sort -k1,1 -k2,2n > {output.bed} 2> {log}
        """


rule crc_pull_yamet_bed:
    """Materialise CRC_BEDS/<subcat>.<cat>.bed from the right BED source
    (plain BED or .bed.gz). Source mapping in _crc_yamet_bed_path."""
    conda:
        op.join("..", "envs", "bedtools.yml")
    wildcard_constraints:
        subcat = _CRC_SUBCAT_RE,
        cat = _CRC_CAT_RE,
    input:
        bed = lambda w: _crc_yamet_bed_path(w.subcat, w.cat),
    output:
        bed = op.join(CRC_BEDS, "{subcat}.{cat}.bed"),
    log:
        op.join(CRC_DATA, "logs", "pull_yamet_bed_{subcat}_{cat}.log"),
    shell:
        r"""
        mkdir -p $(dirname {output.bed})
        if [[ "{input.bed}" == *.gz ]]; then
            zcat {input.bed} > {output.bed} 2> {log}
        else
            cp {input.bed} {output.bed} 2> {log}
        fi
        """


rule crc_stage_annotation_bed:
    """Stage one annotation BED (already in CRC_BEDS as <subcat>.<cat>.bed)
    into CRC_RUN with feature_id = <subcat>_<index>. Filename convention:
    {subcat}.{cat}.bed. Whole-genome, no chr restriction."""
    conda:
        op.join("..", "envs", "bedtools.yml")
    wildcard_constraints:
        subcat = _CRC_SUBCAT_RE,
        cat = _CRC_CAT_RE,
    input:
        bed = op.join(CRC_BEDS, "{subcat}.{cat}.bed"),
    output:
        bed = op.join(CRC_RUN, "beds", "{subcat}.{cat}.bed"),
    log:
        op.join(CRC_RUN, "logs", "stage_bed_{subcat}_{cat}.log"),
    shell:
        r"""
        mkdir -p $(dirname {output.bed})
        awk -v sc={wildcards.subcat} 'BEGIN{{OFS="\t"; k=0}}
                                       {{if ($1 ~ /^chr(X|Y|M|MT)$/) next;
                                         k++; print $1, $2, $3, sc "_" k}}' \
          {input.bed} > {output.bed} 2> {log}
        """


rule crc_window_annotation_per_pair:
    """Per-window overlap fraction with one annotation BED. For a 4-column BED
    `-a`, `bedtools coverage` appends count, bases_covered, length_A, and
    fraction; the fraction is column 8. Output is a single-column file with
    header `<subcat>_<cat>`."""
    wildcard_constraints:
        subcat = _CRC_SUBCAT_RE,
        cat = _CRC_CAT_RE,
    conda:
        op.join("..", "envs", "bedtools.yml")
    input:
        windows = op.join(CRC_RUN, "beds", "windows.bed"),
        ann = op.join(CRC_RUN, "beds", "{subcat}.{cat}.bed"),
    output:
        frac = temp(op.join(CRC_RUN, "beds",
                            "windows_annotation.{subcat}.{cat}.frac")),
    log:
        op.join(CRC_RUN, "logs",
                "window_annotation_{subcat}_{cat}.log"),
    shell:
        r"""
        mkdir -p $(dirname {output.frac})
        (
          echo "{wildcards.subcat}_{wildcards.cat}"
          bedtools coverage -a {input.windows} -b {input.ann} | cut -f8
        ) > {output.frac} 2> {log}
        """


rule crc_combine_window_annotations:
    """Paste the per-pair window fractions next to the (chrom, start, end,
    feature_id) windows BED. Header line is added so downstream readers can
    use read_tsv directly."""
    conda:
        op.join("..", "envs", "bedtools.yml")
    input:
        windows = op.join(CRC_RUN, "beds", "windows.bed"),
        fracs = expand(
            op.join(CRC_RUN, "beds",
                    "windows_annotation.{subcat}.{cat}.frac"),
            zip,
            subcat = [sc for sc, _ in _CRC_LOCAL_PAIRS],
            cat = [c for _, c in _CRC_LOCAL_PAIRS],
        ),
    output:
        tsv = op.join(CRC_RUN, "beds", "windows_annotation.tsv.gz"),
    log:
        op.join(CRC_RUN, "logs", "combine_window_annotations.log"),
    shell:
        r"""
        mkdir -p $(dirname {output.tsv})
        tmp_header=$(mktemp)
        tmp_body=$(mktemp)
        {{
          echo -e "chrom\tstart\tend\tfeature_id" > $tmp_header
          cat $tmp_header {input.windows} > $tmp_body
          paste $tmp_body {input.fracs} | gzip -c > {output.tsv}
        }} 2> {log}
        rm -f $tmp_header $tmp_body
        """


rule run_amet_on_crc_features:
    """Run amet once per (patient, location) combo across every annotation BED.
    Each BED is passed as a separate --features, so the cell files are parsed
    only once for the whole annotation panel. amet writes a cell_feature,
    feature, and pair_counts file per BED, keyed by the staged BED basename
    <subcat>.<cat>; this rule declares the cell_feature and feature files as
    tracked outputs."""
    wildcard_constraints:
        patient = r"[^_.]+",
        location = r"[^_.]+",
    conda:
        op.join("..", "envs", "bedtools.yml")
    input:
        binary = AMET,
        cells = op.join(CRC_DATA, "manifests",
                        "{patient}_{location}.tsv"),
        genome = op.join(REFS, "hg19_ucsc", "genome.fa"),
        cpg = op.join(REFS, "hg19_ucsc", "genome.fa.cpg"),
        beds = [op.join(CRC_RUN, "beds", f"{sc}.{c}.bed")
                for sc, c in _CRC_LOCAL_PAIRS],
    output:
        cell_feature = [
            op.join(CRC_RUN, "features",
                    "{patient}_{location}." + f"{sc}.{c}.cell_feature.tsv.gz")
            for sc, c in _CRC_LOCAL_PAIRS],
        feature = [
            op.join(CRC_RUN, "features",
                    "{patient}_{location}." + f"{sc}.{c}.feature.tsv.gz")
            for sc, c in _CRC_LOCAL_PAIRS],
    params:
        prefix = op.join(CRC_RUN, "features", "{patient}_{location}"),
        features_flags = lambda w, input: " ".join(
            f"--features {b}" for b in input.beds),
        i_max_lag = config["amet"]["i_max_lag"],
        min_cpgs = config["amet"]["min_cpgs_per_feature"],
        min_cells = min_cells_per_group(),
        thresh = config["amet"]["meth_call_threshold"],
    threads: min(workflow.cores, 4)
    log:
        op.join(CRC_RUN, "logs",
                "amet_features_{patient}_{location}.log"),
    shell:
        """
        mkdir -p $(dirname {params.prefix})
        {input.binary} \
            --genome {input.genome} \
            --cells {input.cells} \
            {params.features_flags} \
            --output-prefix {params.prefix} \
            --i-max-lag {params.i_max_lag} \
            --min-cpgs-per-feature {params.min_cpgs} \
            --min-cells-per-group {params.min_cells} \
            --meth-call-threshold {params.thresh} \
            --threads {threads} &> {log}
        """


rule run_amet_on_crc_windows:
    """Run amet on whole-genome windows for one (patient, location) combo.
    win_size is collapsed to a single hard-coded value."""
    conda:
        op.join("..", "envs", "bedtools.yml")
    input:
        binary = AMET,
        cells = op.join(CRC_DATA, "manifests",
                        "{patient}_{location}.tsv"),
        genome = op.join(REFS, "hg19_ucsc", "genome.fa"),
        cpg = op.join(REFS, "hg19_ucsc", "genome.fa.cpg"),
        bed = op.join(CRC_RUN, "beds", "windows.bed"),
    output:
        cell_feature = op.join(
            CRC_RUN, "windows",
            "{patient}_{location}.cell_feature.tsv.gz"),
        feature = op.join(
            CRC_RUN, "windows",
            "{patient}_{location}.feature.tsv.gz"),
    params:
        prefix = op.join(
            CRC_RUN, "windows", "{patient}_{location}"),
        i_max_lag = config["amet"]["i_max_lag"],
        min_cpgs = config["amet"]["min_cpgs_per_feature"],
        min_cells = min_cells_per_group(),
        thresh = config["amet"]["meth_call_threshold"],
    threads: min(workflow.cores, 4)
    log:
        op.join(CRC_RUN, "logs",
                "amet_windows_{patient}_{location}.log"),
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


def _crc_combos():
    """(patient, location) pairs from cells.tsv after the manifest checkpoint."""
    import csv
    manifest_path = checkpoints.crc_make_manifest.get().output.manifest
    pairs = set()
    with open(manifest_path) as f:
        for row in csv.DictReader(f, delimiter="\t"):
            p = row.get("patient")
            l = row.get("location")
            if p and l:
                pairs.add((p, l))
    return sorted(pairs)


def list_crc_features_outputs(wildcards):
    combos = _crc_combos()
    out = []
    for p, l in combos:
        for sc, c in _CRC_LOCAL_PAIRS:
            out.append(op.join(CRC_RUN, "features",
                               f"{p}_{l}.{sc}.{c}.cell_feature.tsv.gz"))
            out.append(op.join(CRC_RUN, "features",
                               f"{p}_{l}.{sc}.{c}.feature.tsv.gz"))
    return out


def list_crc_windows_outputs(wildcards):
    combos = _crc_combos()
    out = []
    for p, l in combos:
        out.append(op.join(CRC_RUN, "windows",
                           f"{p}_{l}.cell_feature.tsv.gz"))
        out.append(op.join(CRC_RUN, "windows",
                           f"{p}_{l}.feature.tsv.gz"))
    return out


def _crc_render_shell(with_annotation = True):
    helpers = op.join(REPO_ROOT, "workflow", "scripts", "render_logging.R")
    i_max_lag = config["amet"]["i_max_lag"]
    annotation_line = (
        '\n                windows_annotation="{{input.windows_annotation}}",'
        if with_annotation else ""
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
                windows_dir="{{params.windows_dir}}",
                win_bed="{{input.win_bed}}",
                manifest="{{input.manifest}}",
                out_dir="{{params.out_dir}}",{annotation_line}
                log_path="{{log}}",
                threads={{threads}},
                i_max_lag={i_max_lag}),
            quiet=TRUE)' &> {{log}}
        """


## The four analytical Rmds run in this order:
##   crc          (per-feature, independent)
##   crc_windows  -> sce_windows_colon.rds + de_list.rds
##   crc_windows_sce  -> sce_windows_colon_corrected.rds
##   crc_embeddings   -> crc_embeddings_debug.rds, crc_win_varexp.csv, crc_per_cell_summary.csv
## RDS/CSV intermediates are declared as snakemake outputs/inputs so the chain
## is enforced by the file graph, not by html-on-html ordering tricks.


rule render_crc:
    conda:
        op.join("..", "envs", "r-tools.yml")
    input:
        rmd = op.join(REPO_ROOT, "workflow", "Rmd", "crc.Rmd"),
        scripts = RMD_SHARED_SCRIPTS + [DRIVER_UTILS_R],
        features = list_crc_features_outputs,
        win_bed = op.join(CRC_RUN, "beds", "windows.bed"),
        manifest = op.join(CRC_DATA, "cells.tsv"),
    output:
        html = op.join(CRC_RUN, "crc.html"),
        entropy_summaries = op.join(CRC_RUN, "crc_entropy_summaries.rds"),
        driver_sd_range = op.join(CRC_RUN, "crc_driver_sd_range.rds"),
    params:
        rmd_name = "crc",
        out_dir = CRC_RUN,
        features_dir = op.join(CRC_RUN, "features"),
        windows_dir = op.join(CRC_RUN, "windows"),
    log:
        op.join(CRC_RUN, "logs", "render_crc.log"),
    threads: 4
    shell:
        _crc_render_shell(with_annotation = False)


rule render_crc_windows:
    conda:
        op.join("..", "envs", "r-tools.yml")
    input:
        rmd = op.join(REPO_ROOT, "workflow", "Rmd", "crc_windows.Rmd"),
        scripts = RMD_SHARED_SCRIPTS + [DRIVER_UTILS_R, DIFF_TESTING_R],
        features = list_crc_features_outputs,
        windows = list_crc_windows_outputs,
        win_bed = op.join(CRC_RUN, "beds", "windows.bed"),
        windows_annotation = op.join(CRC_RUN, "beds", "windows_annotation.tsv.gz"),
        manifest = op.join(CRC_DATA, "cells.tsv"),
    output:
        html = op.join(CRC_RUN, "crc_windows.html"),
        sce_windows = op.join(CRC_RUN, "sce_windows_colon.rds"),
        de_list = op.join(CRC_RUN, "de_list.rds"),
    params:
        rmd_name = "crc_windows",
        out_dir = CRC_RUN,
        features_dir = op.join(CRC_RUN, "features"),
        windows_dir = op.join(CRC_RUN, "windows"),
    log:
        op.join(CRC_RUN, "logs", "render_crc_windows.log"),
    threads: 4
    shell:
        _crc_render_shell()


rule render_crc_windows_sce:
    conda:
        op.join("..", "envs", "r-tools.yml")
    input:
        rmd = op.join(REPO_ROOT, "workflow", "Rmd", "crc_windows_sce.Rmd"),
        scripts = RMD_SHARED_SCRIPTS + [DIFF_TESTING_R],
        sce_windows = op.join(CRC_RUN, "sce_windows_colon.rds"),
        de_list = op.join(CRC_RUN, "de_list.rds"),
        win_bed = op.join(CRC_RUN, "beds", "windows.bed"),
        windows_annotation = op.join(CRC_RUN, "beds", "windows_annotation.tsv.gz"),
        manifest = op.join(CRC_DATA, "cells.tsv"),
    output:
        html = op.join(CRC_RUN, "crc_windows_sce.html"),
        corrected_sce = op.join(CRC_RUN, "sce_windows_colon_corrected.rds"),
    params:
        rmd_name = "crc_windows_sce",
        out_dir = CRC_RUN,
        features_dir = op.join(CRC_RUN, "features"),
        windows_dir = op.join(CRC_RUN, "windows"),
    log:
        op.join(CRC_RUN, "logs", "render_crc_windows_sce.log"),
    threads: 4
    shell:
        _crc_render_shell()


rule render_crc_embeddings:
    conda:
        op.join("..", "envs", "r-tools.yml")
    input:
        rmd = op.join(REPO_ROOT, "workflow", "Rmd", "crc_embeddings.Rmd"),
        scripts = RMD_SHARED_SCRIPTS + [EMBEDDING_UTILS_R],
        corrected_sce = op.join(CRC_RUN, "sce_windows_colon_corrected.rds"),
        win_bed = op.join(CRC_RUN, "beds", "windows.bed"),
        windows_annotation = op.join(CRC_RUN, "beds", "windows_annotation.tsv.gz"),
        manifest = op.join(CRC_DATA, "cells.tsv"),
    output:
        html = op.join(CRC_RUN, "crc_embeddings.html"),
        embeddings_debug = op.join(CRC_RUN, "crc_embeddings_debug.rds"),
        win_varexp = op.join(CRC_RUN, "crc_win_varexp.csv"),
        per_cell_summary = op.join(CRC_RUN, "crc_per_cell_summary.csv"),
    params:
        rmd_name = "crc_embeddings",
        out_dir = CRC_RUN,
        features_dir = op.join(CRC_RUN, "features"),
        windows_dir = op.join(CRC_RUN, "windows"),
    log:
        op.join(CRC_RUN, "logs", "render_crc_embeddings.log"),
    threads: 4
    shell:
        _crc_render_shell()


rule render_fig_crc:
    """Compact CRC figure (single page panels A-H). Consumes the analytical
    Rmds' entropy/driver/varexp/per-cell summaries plus the embeddings debug
    RDS and the de_list."""
    conda:
        op.join("..", "envs", "r-tools.yml")
    input:
        rmd = op.join(REPO_ROOT, "workflow", "Rmd", "fig_crc.Rmd"),
        scripts = RMD_SHARED_SCRIPTS + [DRIVER_UTILS_R],
        entropy_summaries = op.join(CRC_RUN, "crc_entropy_summaries.rds"),
        driver_sd_range = op.join(CRC_RUN, "crc_driver_sd_range.rds"),
        embeddings_debug = op.join(CRC_RUN, "crc_embeddings_debug.rds"),
        win_varexp = op.join(CRC_RUN, "crc_win_varexp.csv"),
        per_cell_summary = op.join(CRC_RUN, "crc_per_cell_summary.csv"),
        de_list = op.join(CRC_RUN, "de_list.rds"),
        win_bed = op.join(CRC_RUN, "beds", "windows.bed"),
        windows_annotation = op.join(CRC_RUN, "beds", "windows_annotation.tsv.gz"),
        manifest = op.join(CRC_DATA, "cells.tsv"),
    output:
        html = op.join(CRC_RUN, "fig_crc.html"),
    params:
        rmd_name = "fig_crc",
        out_dir = CRC_RUN,
        features_dir = op.join(CRC_RUN, "features"),
        windows_dir = op.join(CRC_RUN, "windows"),
    log:
        op.join(CRC_RUN, "logs", "render_fig_crc.log"),
    threads: 4
    shell:
        _crc_render_shell()


rule render_fig_crc_diffentropy:
    """Differential-entropy CRC figure. Consumes de_list, the corrected SCE
    and the embeddings debug RDS; does not need entropy/driver/varexp/per-cell
    artifacts."""
    conda:
        op.join("..", "envs", "r-tools.yml")
    input:
        rmd = op.join(REPO_ROOT, "workflow", "Rmd", "fig_crc_diffentropy.Rmd"),
        scripts = RMD_SHARED_SCRIPTS,
        de_list = op.join(CRC_RUN, "de_list.rds"),
        embeddings_debug = op.join(CRC_RUN, "crc_embeddings_debug.rds"),
        corrected_sce = op.join(CRC_RUN, "sce_windows_colon_corrected.rds"),
        win_bed = op.join(CRC_RUN, "beds", "windows.bed"),
        windows_annotation = op.join(CRC_RUN, "beds", "windows_annotation.tsv.gz"),
        manifest = op.join(CRC_DATA, "cells.tsv"),
    output:
        html = op.join(CRC_RUN, "fig_crc_diffentropy.html"),
    params:
        rmd_name = "fig_crc_diffentropy",
        out_dir = CRC_RUN,
        features_dir = op.join(CRC_RUN, "features"),
        windows_dir = op.join(CRC_RUN, "windows"),
    log:
        op.join(CRC_RUN, "logs", "render_fig_crc_diffentropy.log"),
    threads: 4
    shell:
        _crc_render_shell()
