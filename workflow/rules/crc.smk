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
CRC_RUN = op.join(RESULTS, config["crc"]["run_name"])
CRC_RUN_NAME = config["crc"]["run_name"]


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


rule crc_make_manifest:
    """Parse CRC_RAW filenames into a cells.tsv keyed by patient/location.

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


rule crc_make_windows_bed:
    """chr19 fixed-size windows on hg19 (UCSC chrom names match the singleC files)."""
    conda:
        op.join("..", "envs", "bedtools.yml")
    input:
        sizes = op.join(REFS, "hg19_ucsc", "chr19.sizes"),
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
    """UCSC CpG islands on hg19, restricted to chr19 and stamped with feature_id."""
    conda:
        op.join("..", "envs", "bedtools.yml")
    output:
        bed = op.join(CRC_RUN, "beds", "cpgIslands.bed"),
    params:
        url = "http://hgdownload.cse.ucsc.edu/goldenpath/hg19/database/cpgIslandExt.txt.gz",
        n_max = config["prototype"]["features_subset"],
    log:
        op.join(CRC_RUN, "logs", "make_cgi.log"),
    shell:
        r"""
        mkdir -p $(dirname {output.bed})
        curl -sSL {params.url} 2> {log} \
          | gunzip -c \
          | awk 'BEGIN{{OFS="\t"}} $2 == "chr19" {{print $2, $3, $4, "cgi"}}' \
          | sort -k1,1 -k2,2n \
          | head -n {params.n_max} \
          | awk 'BEGIN{{OFS="\t"; k=0}} {{k++; print $1, $2, $3, "cgi_" k}}' \
          > {output.bed} 2>> {log}
        """


rule crc_run_amet_on_scope:
    """Run amet on a {scope} BED. scope is 'windows' or 'cpgIslands'."""
    wildcard_constraints:
        scope = "windows|cpgIslands",
    conda:
        op.join("..", "envs", "bedtools.yml")
    input:
        binary = AMET,
        cells = op.join(CRC_DATA, "cells.tsv"),
        genome = op.join(REFS, "hg19_ucsc", "chr19.fa"),
        bed = op.join(CRC_RUN, "beds", "{scope}.bed"),
    output:
        cell_feature = op.join(CRC_RUN, CRC_RUN_NAME + ".{scope}.cell_feature.tsv.gz"),
        feature = op.join(CRC_RUN, CRC_RUN_NAME + ".{scope}.feature.tsv.gz"),
    params:
        prefix = op.join(CRC_RUN, CRC_RUN_NAME + ".{scope}"),
        i_max_lag = config["amet"]["i_max_lag"],
        min_cpgs = config["amet"]["min_cpgs_per_feature"],
        min_cells = config["amet"]["min_cells_per_group"],
        thresh = config["amet"]["meth_call_threshold"],
    threads: min(workflow.cores, 8)
    log:
        op.join(CRC_RUN, "logs", "amet_{scope}.log"),
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


rule crc_render_report:
    """Single paper-style report: feature-level (CGI) and window-level analyses."""
    conda:
        op.join("..", "envs", "r-tools.yml")
    input:
        rmd = op.join(REPO_ROOT, "workflow", "Rmd", "crc.Rmd"),
        feat_cell_feature = op.join(CRC_RUN, CRC_RUN_NAME + ".cpgIslands.cell_feature.tsv.gz"),
        feat_feature = op.join(CRC_RUN, CRC_RUN_NAME + ".cpgIslands.feature.tsv.gz"),
        feat_bed = op.join(CRC_RUN, "beds", "cpgIslands.bed"),
        win_cell_feature = op.join(CRC_RUN, CRC_RUN_NAME + ".windows.cell_feature.tsv.gz"),
        win_feature = op.join(CRC_RUN, CRC_RUN_NAME + ".windows.feature.tsv.gz"),
        win_bed = op.join(CRC_RUN, "beds", "windows.bed"),
        manifest = op.join(CRC_DATA, "cells.tsv"),
    output:
        html = op.join(CRC_RUN, CRC_RUN_NAME + ".html"),
    params:
        out_dir = CRC_RUN,
        run_name = CRC_RUN_NAME,
    log:
        op.join(CRC_RUN, "logs", "render.log"),
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
