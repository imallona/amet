## Feature-variability sweep, designed to give scMET its native data shape:
## many features per cell, with feature-level dispersion as the ground truth.
##
## Each cell carries F features of C CpGs each. Feature f has a true dispersion
## gamma_f drawn from a fixed grid (5 levels, F/5 features per level). For every
## cell, we draw a per-cell-per-feature methylation rate r_{cf} from a Beta
## reparameterised by (mu_f, gamma_f), then emit C iid Bernoulli CpG calls at r_{cf}.
##
## This is the regime scMET was designed for: feature-level overdispersion, no
## within-feature spatial structure. amet's I_norm should sit near zero across
## features. amet's per-feature JSD should rise with gamma_f because cells with
## high gamma_f differ in their per-cell rate. scMET should win on per-feature
## gamma recovery.

suppressPackageStartupMessages({ library(optparse) })
.this_dir <- local({ args <- commandArgs(trailingOnly = FALSE); fa <- grep("^--file=", args, value = TRUE); if (length(fa) > 0) dirname(sub("^--file=", "", fa[1])) else "." })
source(file.path(.this_dir, "write_outputs.R"))

options <- list(
    make_option(c("--out_dir"), type = "character"),
    make_option(c("--n_cells"), type = "integer", default = 50),
    make_option(c("--n_features"), type = "integer", default = 100),
    make_option(c("--n_cpgs_per_feature"), type = "integer", default = 100),
    make_option(c("--mu"), type = "numeric", default = 0.5),
    make_option(c("--gamma_grid"), type = "character",
                default = "0.05,0.10,0.20,0.40,0.80"),
    make_option(c("--seed"), type = "integer", default = 1)
)
opt <- parse_args(OptionParser(option_list = options))

dir.create(opt$out_dir, showWarnings = FALSE, recursive = TRUE)
cells_dir <- file.path(opt$out_dir, "cells")
dir.create(cells_dir, showWarnings = FALSE)

gamma_levels <- as.numeric(strsplit(opt$gamma_grid, ",")[[1]])
if (opt$n_features %% length(gamma_levels) != 0) {
    stop("n_features must be a multiple of |gamma_grid|")
}
features_per_level <- opt$n_features %/% length(gamma_levels)
gamma_per_feature <- rep(gamma_levels, each = features_per_level)

C <- opt$n_cpgs_per_feature
F_n <- opt$n_features

chrom <- "chr1"
total_cpgs <- F_n * C
positions_0based <- 99L + (seq_len(total_cpgs) - 1L) * 100L
write_cpg_reference(file.path(opt$out_dir, "cpgs.tsv"), chrom, positions_0based)

bed_df <- data.frame(
    chrom = chrom,
    start = (seq_len(F_n) - 1L) * C * 100L,
    end = seq_len(F_n) * C * 100L,
    name = sprintf("feat%03d", seq_len(F_n))
)
write.table(bed_df, file.path(opt$out_dir, "feature.bed"),
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

ground_truth <- data.frame(feature_id = bed_df$name, gamma = gamma_per_feature, mu = opt$mu)
write.table(ground_truth, file.path(opt$out_dir, "ground_truth.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)

set.seed(opt$seed)

beta_alpha_beta <- function(mu, gamma) {
    nu <- (1 / gamma) - 1
    list(alpha = mu * nu, beta = (1 - mu) * nu)
}
ab_per_feature <- lapply(gamma_per_feature, function(g) beta_alpha_beta(opt$mu, g))

cell_ids <- character(); groups <- character(); paths <- character()
for (i in seq_len(opt$n_cells)) {
    bits <- integer(total_cpgs)
    for (f in seq_len(F_n)) {
        ab <- ab_per_feature[[f]]
        r <- rbeta(1, ab$alpha, ab$beta)
        offset <- (f - 1L) * C
        bits[offset + seq_len(C)] <- as.integer(runif(C) < r)
    }
    id <- sprintf("cell%03d", i)
    path <- file.path(cells_dir, sprintf("%s.allc.tsv.gz", id))
    write_allc(path, chrom, positions_0based, bits)
    cell_ids <- c(cell_ids, id)
    groups <- c(groups, "all")
    paths <- c(paths, path)
}

write_manifest(file.path(opt$out_dir, "manifest.tsv"), cell_ids, groups, paths)
message(sprintf("[generate_feature_variability_sweep] wrote %d cells, %d features (%d levels)",
                length(cell_ids), F_n, length(gamma_levels)))
