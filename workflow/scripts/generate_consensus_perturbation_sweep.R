## Consensus-perturbation sweep, designed to give epiCHAOS its native data shape:
## cells in a group share a consensus methylation pattern, drifted by per-cell
## independent CpG flips at a controlled rate.
##
## Group index pl in 1..10 sets the flip rate:
##   pl = 1  -> flip_rate = 0.00 (cells identical to consensus, low Jaccard distance)
##   pl = 10 -> flip_rate = 0.45 (cells nearly independent of consensus, high distance)
##
## The consensus is drawn iid Bernoulli(p) at marginal p = 0.5, so cells share both
## marginal and methylation positions at low pl. Within-cell spatial structure is
## absent (no Markov persistence), so amet's I_norm should sit near zero across all
## groups while epiCHAOS's eITH grows with pl.

suppressPackageStartupMessages({ library(optparse) })
.this_dir <- local({ args <- commandArgs(trailingOnly = FALSE); fa <- grep("^--file=", args, value = TRUE); if (length(fa) > 0) dirname(sub("^--file=", "", fa[1])) else "." })
source(file.path(.this_dir, "write_outputs.R"))

options <- list(
    make_option(c("--out_dir"), type = "character"),
    make_option(c("--n_cells_per_group"), type = "integer", default = 30),
    make_option(c("--n_cpgs"), type = "integer", default = 2000),
    make_option(c("--marginal"), type = "numeric", default = 0.5),
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

pl_to_flip <- function(pl) (pl - 1L) * 0.05

set.seed(opt$seed)
cell_ids <- character(); groups <- character(); paths <- character()

for (pl in 1:10) {
    consensus <- as.integer(runif(opt$n_cpgs) < opt$marginal)
    flip_rate <- pl_to_flip(pl)
    for (i in seq_len(opt$n_cells_per_group)) {
        flips <- runif(opt$n_cpgs) < flip_rate
        bits <- consensus
        bits[flips] <- 1L - bits[flips]
        id <- sprintf("pl%02d_cell%03d", pl, i)
        path <- file.path(cells_dir, sprintf("%s.allc.tsv.gz", id))
        write_allc(path, chrom, positions_0based, bits)
        cell_ids <- c(cell_ids, id)
        groups <- c(groups, sprintf("pl%02d", pl))
        paths <- c(paths, path)
    }
}

write_manifest(file.path(opt$out_dir, "manifest.tsv"), cell_ids, groups, paths)
message(sprintf("[generate_consensus_perturbation_sweep] wrote %d cells across 10 levels",
                length(cell_ids)))
