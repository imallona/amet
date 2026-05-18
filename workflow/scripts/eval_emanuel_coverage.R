## Evaluation: Emanuel Sonder's coverage simulations scored by amet.
##
## Plots the within-cell score i_total and its analytical normalization
## i_norm = i_total / (k_max * H(mean_meth)) against mean methylation, faceted
## by CpG count and coverage regime and coloured by transition matrix. The
## simulation is Emanuel Sonder's; this amet evaluation is the adaptation of
## his yamet plotting Rmd.

suppressPackageStartupMessages({
    library(optparse)
    library(ggplot2)
    library(data.table)
})

.this_dir <- local({
    args <- commandArgs(trailingOnly = FALSE)
    fa <- grep("^--file=", args, value = TRUE)
    if (length(fa) > 0) dirname(sub("^--file=", "", fa[1])) else "."
})
source(file.path(.this_dir, "plot_theme.R"))

opt <- parse_args(OptionParser(option_list = list(
    make_option("--cell_feature", type = "character",
                help = "amet combined cell_feature TSV for the coverage sims"),
    make_option("--output_prefix", type = "character",
                help = "output prefix; writes <prefix>_{i_total,i_norm,pairwise}.{pdf,svg,csv}")
)))

shannon_binary <- function(p) {
    out <- numeric(length(p))
    safe <- !is.na(p) & p > 0 & p < 1
    out[safe] <- -p[safe] * log2(p[safe]) - (1 - p[safe]) * log2(1 - p[safe])
    out[!safe] <- NA_real_
    out
}

cf <- fread(opt$cell_feature)
i_cols <- grep("^i_[0-9]+$", names(cf), value = TRUE)
k_max <- length(i_cols)
cf[, i_norm := i_total / (k_max * shannon_binary(mean_meth))]

## cell_id: sim_cell_<i>_<n_cells>_<n_cpgs>_<mode>_<coverage>_<transMat>
cf[, c("ncpgs", "type", "coverage", "tr") :=
   tstrsplit(cell_id, "_", keep = 5:8, type.convert = TRUE)]
cf <- cf[is.finite(i_total) & is.finite(mean_meth)]

cf[, ncpgs := factor(ncpgs, levels = sort(unique(as.integer(ncpgs))))]
cf[, coverage := factor(coverage,
   levels = intersect(c("low", "lowReal", "real", "medium", "high", "complete"),
                      unique(coverage)))]
cf[, type := factor(type)]
cf[, tr := factor(tr, levels = intersect(c("lmr", "imrCons", "imrRand", "hmr"),
                                         unique(tr)))]

facet <- facet_grid(ncpgs ~ coverage,
    labeller = labeller(ncpgs = function(x) paste("# CpGs:", x),
                        coverage = function(x) paste("cov.:", x)))

scatter <- function(yvar, ylab) {
    ggplot(cf, aes(x = mean_meth, y = .data[[yvar]], colour = tr, shape = type)) +
        geom_point(alpha = 0.6, size = 0.7) +
        facet +
        labs(x = "mean methylation", y = ylab,
             colour = "transition matrix", shape = "coverage model") +
        theme_ng_discrete()
}

p_i_total <- scatter("i_total", expression("i"["total"] * " (lag 1.." * "k MI sum)"))
p_i_norm <- scatter("i_norm", expression("i"["norm"] * " = i"["total"] * " / (k"["max"] * " H(p))"))
p_pairwise <- ggplot(cf, aes(x = i_total, y = i_norm, colour = tr, shape = type)) +
    geom_point(alpha = 0.6, size = 0.7) +
    facet +
    labs(x = expression("i"["total"]), y = expression("i"["norm"]),
         colour = "transition matrix", shape = "coverage model") +
    theme_ng_discrete()

save_eval(p_i_total, cf, paste0(opt$output_prefix, "_i_total"),
          width_mm = 200, height_mm = 170)
save_eval(p_i_norm, cf, paste0(opt$output_prefix, "_i_norm"),
          width_mm = 200, height_mm = 170)
save_eval(p_pairwise, cf, paste0(opt$output_prefix, "_pairwise"),
          width_mm = 200, height_mm = 170)
