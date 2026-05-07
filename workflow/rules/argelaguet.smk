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

ARG_RUN_NAME = config["argelaguet"]["run_name"]
ARG_RUN = op.join(RESULTS, ARG_RUN_NAME)

## Gastrulation-specific features bundled in the scnmt_gastrulation tarball
## (features/genomic_contexts/). H3K27ac distal peaks: GSE125318 companion
## ChIP-seq (see Methods). ESC marks: ENCODE.
_ARGELAGUET_GASTRO_BEDS = [
    "H3K27ac_distal_E7.5_Ect_intersect12",
    "H3K27ac_distal_E7.5_End_intersect12",
    "H3K27ac_distal_E7.5_Mes_intersect12",
    "ESC_p300",
]


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


rule make_argelaguet_manifest:
    """Build amet cells.tsv from the filtered metadata + local cell files."""
    conda:
        op.join("..", "envs", "r-tools.yml")
    input:
        meta = op.join(ARG_DATA, "meta.tsv.gz"),
    output:
        manifest = op.join(ARG_DATA, "cells.tsv"),
    params:
        cells_dir = ARG_CELLS,
        group_col = config["argelaguet"]["group_column"],
    log:
        op.join(ARG_DATA, "logs", "manifest.log"),
    shell:
        """
        Rscript {workflow.basedir}/scripts/make_manifest_argelaguet.R \
            --metadata {input.meta} \
            --cells_dir {params.cells_dir} \
            --group_col {params.group_col} \
            --out {output.manifest} &> {log}
        """


rule argelaguet_filter_mm10_bed:
    """Subset a feature BED to chr19 (prototype scope), cap to N intervals, and
    stamp each peak with a unique feature_id of the form <annotation>_<index>.
    The Rmd recovers the annotation by stripping the trailing _<digits>.
    """
    conda:
        op.join("..", "envs", "bedtools.yml")
    input:
        bed = op.join(ARG_FEATURES_DIR, "{annotation}.bed"),
    output:
        bed = temp(op.join(ARG_RUN, "beds", "{annotation}.bed")),
    params:
        n_max = config["prototype"]["features_subset"],
    log:
        op.join(ARG_RUN, "logs", "filter_bed_{annotation}.log"),
    shell:
        r"""
        mkdir -p $(dirname {output.bed})
        echo "[filter_bed] {wildcards.annotation}: cap=N={params.n_max}, chr19 only" > {log}
        awk -v ann={wildcards.annotation} -v n={params.n_max} '
             BEGIN{{OFS="\t"; k=0}}
             {{ chr=$1; sub(/^chr/, "", chr);
                if (chr == "19" && ++k <= n)
                    print $1, $2, $3, ann "_" k }}' {input.bed} > {output.bed}
        echo "[filter_bed] kept $(wc -l < {output.bed}) intervals" >> {log}
        """


rule argelaguet_prep_gastro_bed:
    """Concatenate per-annotation BEDs into one merged feature BED."""
    conda:
        op.join("..", "envs", "bedtools.yml")
    input:
        beds = expand(op.join(ARG_RUN, "beds", "{annotation}.bed"),
                      annotation = _ARGELAGUET_GASTRO_BEDS),
    output:
        bed = op.join(ARG_RUN, "beds", "gastro_features.bed"),
    log:
        op.join(ARG_RUN, "logs", "prep_gastro_bed.log"),
    shell:
        r"""
        cat {input.beds} | sort -k1,1 -k2,2n > {output.bed}
        echo "[prep_gastro] merged $(wc -l < {output.bed}) intervals from $(echo {input.beds} | wc -w) BEDs" > {log}
        """


rule argelaguet_make_windows:
    """Tiled chr19 windows for the windows/embeddings analysis. Each row is a
    distinct feature_id (chrom:start-end). Whole chr19, prototype-scope size."""
    conda:
        op.join("..", "envs", "bedtools.yml")
    input:
        sizes = op.join(REFS, "mm10_ucsc", "chr19.sizes"),
    output:
        bed = op.join(ARG_RUN, "beds", "windows.bed"),
    params:
        win_size = config["argelaguet"]["window_size"],
    shell:
        r"""
        # Strip 'chr' to match Ensembl-named cells; amet's chrom harmonization
        # would also handle the mismatch, but Stripping here keeps the windows
        # BED self-consistent with the gastro_features BED.
        bedtools makewindows -g {input.sizes} -w {params.win_size} \
          | awk 'BEGIN{{OFS="\t"}} {{sub(/^chr/,"",$1); print $1, $2, $3, "win_"NR}}' \
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


rule run_amet_on_argelaguet_scope:
    """Run amet on a {scope} BED. scope is 'features' or 'windows'."""
    wildcard_constraints:
        scope = "features|windows",
    conda:
        op.join("..", "envs", "bedtools.yml")
    input:
        binary = AMET,
        cells = op.join(ARG_DATA, "cells.tsv"),
        genome = op.join(REFS, "mm10_ucsc", "chr19.fa"),
        bed = lambda w: op.join(
            ARG_RUN, "beds",
            "gastro_features.bed" if w.scope == "features" else "windows.bed",
        ),
    output:
        cell_feature = op.join(ARG_RUN, ARG_RUN_NAME + ".{scope}.cell_feature.tsv.gz"),
        feature = op.join(ARG_RUN, ARG_RUN_NAME + ".{scope}.feature.tsv.gz"),
    params:
        prefix = op.join(ARG_RUN, ARG_RUN_NAME + ".{scope}"),
        i_max_lag = config["amet"]["i_max_lag"],
        min_cpgs = config["amet"]["min_cpgs_per_feature"],
        min_cells = config["amet"]["min_cells_per_group"],
        thresh = config["amet"]["meth_call_threshold"],
    threads: min(workflow.cores, 8)
    log:
        op.join(ARG_RUN, "logs", "amet_{scope}.log"),
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


rule render_argelaguet_report:
    """Single paper report combining feature-level and window-level analyses."""
    conda:
        op.join("..", "envs", "r-tools.yml")
    input:
        rmd = op.join(REPO_ROOT, "workflow", "Rmd", "argelaguet_paper.Rmd"),
        feat_cell_feature = op.join(ARG_RUN, ARG_RUN_NAME + ".features.cell_feature.tsv.gz"),
        feat_feature = op.join(ARG_RUN, ARG_RUN_NAME + ".features.feature.tsv.gz"),
        feat_bed = op.join(ARG_RUN, "beds", "gastro_features.bed"),
        win_cell_feature = op.join(ARG_RUN, ARG_RUN_NAME + ".windows.cell_feature.tsv.gz"),
        win_feature = op.join(ARG_RUN, ARG_RUN_NAME + ".windows.feature.tsv.gz"),
        win_bed = op.join(ARG_RUN, "beds", "windows.bed"),
        manifest = op.join(ARG_DATA, "cells.tsv"),
    output:
        html = op.join(ARG_RUN, ARG_RUN_NAME + ".html"),
    params:
        out_dir = ARG_RUN,
        run_name = ARG_RUN_NAME,
    log:
        op.join(ARG_RUN, "logs", "render.log"),
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
