## Generate cells with known correlation structure at different lags:
##   - iid:      no comethylation at any lag
##   - markov_p: AR(1)-style first-order persistence (lag-1 dominant)
##   - period_3: deterministic 010010010..., zero lag-1 MI but strong lag-3
##   - period_5: deterministic 0001100011..., zero short-lag MI but strong lag-5

suppressPackageStartupMessages({ library(optparse) })
source(file.path(dirname(sys.frame(1)$ofile), "write_outputs.R"))
source(file.path(dirname(sys.frame(1)$ofile), "simulate_markov_chain.R"))

options <- list(
    make_option(c("--out_dir"), type = "character"),
    make_option(c("--n_cells_per_pattern"), type = "integer", default = 20),
    make_option(c("--n_cpgs"), type = "integer", default = 200),
    make_option(c("--seed"), type = "integer", default = 1)
)
opt <- parse_args(OptionParser(option_list = options))

dir.create(opt$out_dir, showWarnings = FALSE, recursive = TRUE)
cells_dir <- file.path(opt$out_dir, "cells")
dir.create(cells_dir, showWarnings = FALSE)

chrom <- "chr1"
positions_0based <- 99L + (seq_len(opt$n_cpgs) - 1L) * 100L
write_cpg_reference(file.path(opt$out_dir, "cpgs.tsv"), chrom, positions_0based)
write_bed_feature(file.path(opt$out_dir, "feature.bed"), chrom, 0L,
                  positions_0based[length(positions_0based)] + 100L, "region")

iid_cell <- function(n) simulate_markov_cell(n, p00 = 0.5, p11 = 0.5)
markov_cell <- function(n) simulate_markov_cell(n, p00 = 0.85, p11 = 0.85)
period_n_cell <- function(n, period) {
    base <- c(rep(0L, ceiling(period / 2)), rep(1L, floor(period / 2)))
    rep(base, length.out = n)
}

set.seed(opt$seed)
cell_ids <- character(); groups <- character(); paths <- character()

for (pattern in c("iid", "markov", "period3", "period5")) {
    for (i in seq_len(opt$n_cells_per_pattern)) {
        bits <- switch(pattern,
            iid     = iid_cell(opt$n_cpgs),
            markov  = markov_cell(opt$n_cpgs),
            period3 = period_n_cell(opt$n_cpgs, 3L),
            period5 = period_n_cell(opt$n_cpgs, 5L)
        )
        id <- sprintf("%s_seed%d", pattern, i)
        path <- file.path(cells_dir, sprintf("%s.allc.tsv.gz", id))
        write_allc(path, chrom, positions_0based, bits)
        cell_ids <- c(cell_ids, id)
        groups <- c(groups, pattern)
        paths <- c(paths, path)
    }
}

write_manifest(file.path(opt$out_dir, "manifest.tsv"), cell_ids, groups, paths)
message(sprintf("[generate_lag_profile] wrote %d cells across 4 patterns",
                length(cell_ids)))
