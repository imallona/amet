"""
Emanuel Sonder's coverage simulations, scored by amet.

simulations_01_sim_data.Rmd (ported from yamet, authored by Emanuel Sonder)
generates per-cell methylation tables across a grid of CpG count, coverage
model and transition matrix. The lowReal coverage regime draws real coverage
patterns from the argelaguet gastrulation cpg_level cells, so the generator
depends on the argelaguet data, mirroring yamet's input/output chaining.

amet then scores each ncpg as a single feature spanning every CpG;
eval_emanuel_coverage.R turns the combined cell_feature table into the
i_total, i_norm and pairwise panels of the simulations report.
"""

EMANUEL_SIM = op.join(SIM, "emanuel")
EMANUEL_SIM_DATA = op.join(EMANUEL_SIM, "sim_data")
EMANUEL_AMET = op.join(SIM, "amet", "emanuel_coverage")
EMANUEL_PARAMS = op.join(REPO_ROOT, "workflow", "resources",
                         "emanuel_parameters.tsv")
## Representative output of the generator; the first parameter-grid row is
## 50 cells x 50 CpGs, rand mode, low coverage, lmr transition matrix.
EMANUEL_SIM_FLAG = op.join(EMANUEL_SIM_DATA,
                           "sim_cell_1_50_50_rand_low_lmr.tsv")


rule generate_emanuel_sim_data:
    """Render the ported simulator: per-cell sim tables and cpgPositions files
    across the coverage grid. Depends on the argelaguet gastrulation cpg_level
    cells for the lowReal coverage regime."""
    conda:
        op.join("..", "envs", "sim_r_emanuel.yml")
    input:
        rmd = op.join(REPO_ROOT, "workflow", "Rmd", "simulations_01_sim_data.Rmd"),
        simpattern = op.join(REPO_ROOT, "workflow", "scripts", "simPattern.R"),
        parameters = EMANUEL_PARAMS,
        argelaguet = op.join(RESULTS, "argelaguet", "cells.tsv"),
    output:
        sim_data = directory(EMANUEL_SIM_DATA),
        flag = EMANUEL_SIM_FLAG,
        html = op.join(EMANUEL_SIM, "simulations_01_sim_data.html"),
    params:
        low_real_dir = op.join(RESULTS, "argelaguet", "cells"),
        out_dir = EMANUEL_SIM_DATA,
        report_dir = EMANUEL_SIM,
    log:
        op.join(SIM, "logs", "generate_emanuel_sim_data.log"),
    shell:
        r"""
        mkdir -p {params.out_dir} $(dirname {log})
        Rscript -e 'rmarkdown::render("{input.rmd}",
            output_file="simulations_01_sim_data.html",
            output_dir="{params.report_dir}",
            knit_root_dir="{params.out_dir}",
            params=list(
                parameters_path="{input.parameters}",
                low_real_dir="{params.low_real_dir}",
                out_dir="{params.out_dir}"),
            quiet=TRUE)' &> {log}
        """


rule run_amet_emanuel_coverage:
    """Convert each simulated cell to allc, run amet once per ncpg with a single
    whole-axis feature, and combine the per-ncpg cell_feature outputs."""
    conda:
        op.join("..", "envs", "bedtools.yml")
    input:
        binary = AMET,
        runner = op.join(REPO_ROOT, "workflow", "scripts", "run_emanuel_coverage.sh"),
        convert = op.join(REPO_ROOT, "workflow", "scripts", "convert_sim.sh"),
        sim_data = EMANUEL_SIM_DATA,
        flag = EMANUEL_SIM_FLAG,
    output:
        cell_feature = op.join(EMANUEL_AMET, "all_cells.tsv.gz"),
    params:
        out_base = EMANUEL_AMET,
        i_max_lag = config["amet"]["i_max_lag"],
        min_cpgs = config["amet"]["min_cpgs_per_feature"],
        min_cells = config["amet"]["min_cells_per_group"],
    threads: 8
    log:
        op.join(SIM, "logs", "run_amet_emanuel_coverage.log"),
    shell:
        """
        mkdir -p {params.out_base} $(dirname {log})
        bash {input.runner} {input.binary} {input.sim_data} {input.convert} \
            {params.out_base} {params.i_max_lag} {params.min_cpgs} \
            {params.min_cells} {threads} &> {log}
        """


rule eval_emanuel_coverage:
    """Render the i_total, i_norm and pairwise panels for the report."""
    conda:
        R_TOOLS_ENV
    input:
        script = op.join(REPO_ROOT, "workflow", "scripts", "eval_emanuel_coverage.R"),
        theme = op.join(REPO_ROOT, "workflow", "scripts", "plot_theme.R"),
        cell_feature = op.join(EMANUEL_AMET, "all_cells.tsv.gz"),
    output:
        i_total = multiext(op.join(SIM, "eval", "emanuel_coverage_i_total"),
                           ".pdf", ".svg", ".csv"),
        i_norm = multiext(op.join(SIM, "eval", "emanuel_coverage_i_norm"),
                          ".pdf", ".svg", ".csv"),
        pairwise = multiext(op.join(SIM, "eval", "emanuel_coverage_pairwise"),
                            ".pdf", ".svg", ".csv"),
    params:
        prefix = op.join(SIM, "eval", "emanuel_coverage"),
    log:
        op.join(SIM, "logs", "eval_emanuel_coverage.log"),
    shell:
        """
        mkdir -p $(dirname {params.prefix}) $(dirname {log})
        Rscript {input.script} \
            --cell_feature {input.cell_feature} \
            --output_prefix {params.prefix} &> {log}
        """
