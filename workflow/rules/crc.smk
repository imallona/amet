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
  (b) download_from_geo: cold-start path mirroring yamet's GEO rules. Pulls
      the GSM list from SRA (esearch + efetch) and per-sample directories
      from GEO. Used when starting from scratch with no barbara mirror.

Both paths populate the same results/crc/raw/ directory; the manifest rule
reads whatever is there.
"""

CRC_DATA = op.join(RESULTS, "crc")
CRC_RAW = op.join(CRC_DATA, "raw")
CRC_BEDS = op.join(CRC_DATA, "beds")
CRC_RUN = op.join(RESULTS, config["crc"]["run_name"])
CRC_RUN_NAME = config["crc"]["run_name"]

## yamet's CRC ANNOTATIONS dict (rules/crc.smk in yamet). Keys: cat (outer),
## values: list of subcats (inner). Wildcards mirror yamet exactly:
## {subcat}_{cat}_{patient}_{location}.
##
## Locally we only have cpgIslandExt (downloaded from UCSC). The rest live on
## barbara as bed files; they get pulled in when running there. The dict keeps
## yamet's full surface so rules referencing wildcard combos remain stable.
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

## Annotations whose source BEDs are present locally (built from UCSC etc.).
## Extend this set as more annotations are pulled.
_CRC_LOCAL_ANNOTATIONS = {
    "cpgIslandExt": ["cpgIslandExt"],
}

_CRC_LOCAL_PAIRS = [(sc, c) for c, subs in _CRC_LOCAL_ANNOTATIONS.items() for sc in subs]
_CRC_ALL_PAIRS   = [(sc, c) for c, subs in CRC_ANNOTATIONS.items() for sc in subs]
_CRC_SUBCAT_RE = "|".join(sorted({sc for sc, _ in _CRC_ALL_PAIRS}))
_CRC_CAT_RE    = "|".join(sorted({c  for _, c  in _CRC_ALL_PAIRS}))


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
        proto_patients = ",".join(config["crc"]["proto_patients"]),
        proto_locations = ",".join(config["crc"]["proto_locations"]),
        cells_per_group = config["prototype"]["cells_per_group"],
        prototype = "true" if config["prototype"]["enabled"] else "false",
    log:
        op.join(CRC_DATA, "logs", "manifest.log"),
    shell:
        """
        Rscript {workflow.basedir}/scripts/make_manifest_crc.R \
            --raw_dir {params.raw_dir} \
            --proto_patients {params.proto_patients} \
            --proto_locations {params.proto_locations} \
            --cells_per_group {params.cells_per_group} \
            --prototype {params.prototype} \
            --out {output.manifest} &> {log}
        """


rule crc_per_combo_manifest:
    """Sub-manifest for one (patient, location) combo. Mirrors yamet's
    get_harmonized_files."""
    conda:
        op.join("..", "envs", "python.yml")
    input:
        cells = op.join(CRC_DATA, "cells.tsv"),
    output:
        manifest = op.join(CRC_DATA, "manifests",
                           "{patient}_{location}.tsv"),
    params:
        max_cells = config["prototype"]["cells_per_group"],
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
          | awk 'BEGIN{{OFS="\t"}} {{print $1, $2, $3, "win_"NR}}' \
          | sort -k1,1 -k2,2n > {output.bed} 2> {log}
        """


rule crc_make_cgi_bed:
    """UCSC CpG islands on hg19, whole-genome, stamped with feature_id.
    Output filename mirrors yamet's {subcat}.{cat}.bed convention so all
    annotation BEDs share one filename pattern."""
    conda:
        op.join("..", "envs", "bedtools.yml")
    output:
        bed = op.join(CRC_BEDS, "cpgIslandExt.cpgIslandExt.bed"),
    params:
        url = "https://hgdownload.cse.ucsc.edu/goldenpath/hg19/database/cpgIslandExt.txt.gz",
    log:
        op.join(CRC_DATA, "logs", "make_cgi.log"),
    shell:
        r"""
        mkdir -p $(dirname {output.bed})
        curl -sSL {params.url} 2> {log} \
          | gunzip -c \
          | awk 'BEGIN{{OFS="\t"; k=0}}
                 {{k++; print $2, $3, $4, "cpgIslandExt_" k}}' \
          | sort -k1,1 -k2,2n > {output.bed} 2>> {log}
        """


rule crc_stage_annotation_bed:
    """Stage one annotation BED (already in CRC_BEDS as <subcat>.<cat>.bed)
    into CRC_RUN with feature_id = <subcat>_<index>. Mirrors yamet's
    {subcat}.{cat}.bed naming. Whole-genome, no chr restriction."""
    conda:
        op.join("..", "envs", "bedtools.yml")
    wildcard_constraints:
        subcat = _CRC_SUBCAT_RE,
        cat    = _CRC_CAT_RE,
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
                                       {{k++; print $1, $2, $3, sc "_" k}}' \
          {input.bed} > {output.bed} 2> {log}
        """


rule run_amet_on_crc_features:
    """Run amet on one (subcat, cat, patient, location) combo. Mirrors
    yamet's run_yamet_on_separate_features wildcards exactly."""
    wildcard_constraints:
        subcat = _CRC_SUBCAT_RE,
        cat    = _CRC_CAT_RE,
    conda:
        op.join("..", "envs", "bedtools.yml")
    input:
        binary = AMET,
        cells = op.join(CRC_DATA, "manifests",
                        "{patient}_{location}.tsv"),
        genome = op.join(REFS, "hg19_ucsc", "genome.fa"),
        cpg = op.join(REFS, "hg19_ucsc", "genome.fa.cpg"),
        bed = op.join(CRC_RUN, "beds", "{subcat}.{cat}.bed"),
    output:
        cell_feature = op.join(
            CRC_RUN, "features",
            "{subcat}_{cat}_{patient}_{location}.cell_feature.tsv.gz"),
        feature = op.join(
            CRC_RUN, "features",
            "{subcat}_{cat}_{patient}_{location}.feature.tsv.gz"),
    params:
        prefix = op.join(
            CRC_RUN, "features",
            "{subcat}_{cat}_{patient}_{location}"),
        i_max_lag = config["amet"]["i_max_lag"],
        min_cpgs = config["amet"]["min_cpgs_per_feature"],
        min_cells = config["amet"]["min_cells_per_group"],
        thresh = config["amet"]["meth_call_threshold"],
    threads: min(workflow.cores, 4)
    log:
        op.join(CRC_RUN, "logs",
                "amet_{subcat}_{cat}_{patient}_{location}.log"),
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


rule run_amet_on_crc_windows:
    """Run amet on whole-genome windows for one (patient, location) combo.
    Mirrors yamet's run_yamet_on_crc_windows wildcards (collapsed
    win_size to a single hard-coded value)."""
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
        min_cells = config["amet"]["min_cells_per_group"],
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
    for sc, c in _CRC_LOCAL_PAIRS:
        for p, l in combos:
            out.append(op.join(CRC_RUN, "features",
                               f"{sc}_{c}_{p}_{l}.cell_feature.tsv.gz"))
            out.append(op.join(CRC_RUN, "features",
                               f"{sc}_{c}_{p}_{l}.feature.tsv.gz"))
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


def _crc_render_shell():
    return r"""
        mkdir -p {params.out_dir}
        Rscript -e 'rmarkdown::render("{input.rmd}",
            output_file="{wildcards.rmd_name}.html",
            output_dir="{params.out_dir}",
            knit_root_dir="{params.out_dir}",
            params=list(
                features_dir="{params.features_dir}",
                windows_dir="{params.windows_dir}",
                win_bed="{input.win_bed}",
                manifest="{input.manifest}",
                out_dir="{params.out_dir}"),
            quiet=TRUE)' &> {log}
        """


rule render_crc_analytical_rmd:
    """Render one of the four analytical CRC Rmds (crc, _windows,
    _windows_sce, _embeddings). Each writes RDS/CSV intermediates that
    the figure Rmds consume."""
    wildcard_constraints:
        rmd_name = "crc|crc_windows|crc_windows_sce|crc_embeddings",
    conda:
        op.join("..", "envs", "r-tools.yml")
    input:
        rmd = op.join(REPO_ROOT, "workflow", "Rmd", "{rmd_name}.Rmd"),
        features = list_crc_features_outputs,
        windows = list_crc_windows_outputs,
        win_bed = op.join(CRC_RUN, "beds", "windows.bed"),
        manifest = op.join(CRC_DATA, "cells.tsv"),
    output:
        html = op.join(CRC_RUN, "{rmd_name}.html"),
    params:
        out_dir = CRC_RUN,
        features_dir = op.join(CRC_RUN, "features"),
        windows_dir = op.join(CRC_RUN, "windows"),
    log:
        op.join(CRC_RUN, "logs", "render_{rmd_name}.log"),
    shell:
        _crc_render_shell()


rule render_fig_crc_rmd:
    """Render fig_crc.Rmd or fig_crc_diffentropy.Rmd; depends on the four
    analytical Rmds because it loads their RDS intermediates."""
    wildcard_constraints:
        rmd_name = "fig_crc|fig_crc_diffentropy",
    conda:
        op.join("..", "envs", "r-tools.yml")
    input:
        rmd = op.join(REPO_ROOT, "workflow", "Rmd", "{rmd_name}.Rmd"),
        analytical = expand(op.join(CRC_RUN, "{r}.html"),
                            r = ["crc",
                                 "crc_windows",
                                 "crc_windows_sce",
                                 "crc_embeddings"]),
        features = list_crc_features_outputs,
        windows = list_crc_windows_outputs,
        win_bed = op.join(CRC_RUN, "beds", "windows.bed"),
        manifest = op.join(CRC_DATA, "cells.tsv"),
    output:
        html = op.join(CRC_RUN, "{rmd_name}.html"),
    params:
        out_dir = CRC_RUN,
        features_dir = op.join(CRC_RUN, "features"),
        windows_dir = op.join(CRC_RUN, "windows"),
    log:
        op.join(CRC_RUN, "logs", "render_{rmd_name}.log"),
    shell:
        _crc_render_shell()
