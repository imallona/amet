## Run scMET on amet's per-cell-per-feature output (uses n_ones / n_covered as the
## methylated and total read counts respectively).
##
## scMET fits a Bayesian beta-binomial GLM and reports per-feature:
##   mu      : posterior mean methylation
##   gamma   : overdispersion (the closest existing analog to amet's JSD)
##   epsilon : residual overdispersion after the GLM trend
##
## Input format expected by scmet::scmet():
##   data.frame with columns Feature, Cell, total_reads, met_reads.

suppressPackageStartupMessages({
    library(optparse); library(data.table)
})

options <- list(
    make_option(c("--cell_feature"), type = "character"),
    make_option(c("--output"), type = "character"),
    make_option(c("--iter"), type = "integer", default = 5000),
    make_option(c("--L"), type = "integer", default = 4,
                help = "Number of basis functions for the mu-gamma trend"),
    make_option(c("--threads"), type = "integer", default = 4),
    make_option(c("--seed"), type = "integer", default = 42)
)
opt <- parse_args(OptionParser(option_list = options))

df <- fread(opt$cell_feature, header = TRUE, sep = "\t", na.strings = "NA")
df <- df[!is.na(n_ones) & !is.na(n_zeros) & (n_ones + n_zeros) > 0L]

## In our simulator each cell sits in exactly one (acvi/wcvi, target_p) group, so
## the natural scMET "feature" is the group identifier. With 19 target_p × 10 acvi
## (or 10 wcvi) × ~30 cells, scMET sees ~190 features each with ~30 observations.
scmet_Y <- df[, .(Feature = group,
                  Cell = cell_id,
                  total_reads = as.integer(n_ones + n_zeros),
                  met_reads = as.integer(n_ones))]

message(sprintf("[run_scmet] cells = %d, features = %d, rows = %d",
                length(unique(scmet_Y$Cell)),
                length(unique(scmet_Y$Feature)),
                nrow(scmet_Y)))

fit <- suppressWarnings(scMET::scmet(
    Y = scmet_Y,
    L = opt$L,
    iter = opt$iter,
    n_cores = opt$threads,
    seed = opt$seed
))

## Per-feature posterior point estimates: posterior$mu/gamma/epsilon are
## (n_samples x n_features) matrices; column medians give per-feature scores.
out <- data.frame(
    group = fit$feature_names,
    mu = apply(fit$posterior$mu, 2, median),
    gamma = apply(fit$posterior$gamma, 2, median),
    epsilon = apply(fit$posterior$epsilon, 2, median)
)
fwrite(out, opt$output, sep = "\t")
message(sprintf("[run_scmet] wrote %s (%d groups)", opt$output, nrow(out)))
