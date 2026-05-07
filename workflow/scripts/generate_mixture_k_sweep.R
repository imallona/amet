## Sweep number of mixture components K at fixed component divergence.
## Each "group" is one (K, replicate) combination. Within a group, n_cells cells are
## drawn from a K-component mixture where components have well-separated transition
## matrices.

suppressPackageStartupMessages({ library(optparse) })
.this_dir <- local({ args <- commandArgs(trailingOnly = FALSE); fa <- grep("^--file=", args, value = TRUE); if (length(fa) > 0) dirname(sub("^--file=", "", fa[1])) else "." })

source(file.path(.this_dir, "write_outputs.R"))
source(file.path(.this_dir, "simulate_markov_chain.R"))
source(file.path(.this_dir, "mix_markovs.R"))

build_components <- function(K) {
    ## K components evenly spaced across the persistence axis at fixed pi_1 = 0.5.
    ## c values from 0.30 to 0.90 (anti-persistent -> strongly persistent).
    cs <- seq(0.30, 0.90, length.out = K)
    lapply(cs, function(c) list(p00 = c, p11 = c))
}

options <- list(
    make_option(c("--out_dir"), type = "character"),
    make_option(c("--n_cells_per_k"), type = "integer", default = 60),
    make_option(c("--n_cpgs"), type = "integer", default = 200),
    make_option(c("--k_grid"), type = "character", default = "1,2,4,8"),
    make_option(c("--seed"), type = "integer", default = 1)
)
opt <- parse_args(OptionParser(option_list = options))

dir.create(opt$out_dir, showWarnings = FALSE, recursive = TRUE)
cells_dir <- file.path(opt$out_dir, "cells")
dir.create(cells_dir, showWarnings = FALSE)

ks <- as.integer(strsplit(opt$k_grid, ",")[[1]])
chrom <- "chr1"
positions_0based <- 99L + (seq_len(opt$n_cpgs) - 1L) * 100L
write_cpg_reference(file.path(opt$out_dir, "cpgs.tsv"), chrom, positions_0based)
write_bed_feature(file.path(opt$out_dir, "feature.bed"), chrom, 0L,
                  positions_0based[length(positions_0based)] + 100L, "region")

set.seed(opt$seed)
cell_ids <- character(); groups <- character(); paths <- character()

for (K in ks) {
    components <- build_components(K)
    out <- simulate_mixture_cells(opt$n_cells_per_k, opt$n_cpgs, components,
                                  seed = opt$seed + K)
    for (i in seq_along(out$cells)) {
        id <- sprintf("K%d_cell%d", K, i)
        path <- file.path(cells_dir, sprintf("%s.allc.tsv.gz", id))
        write_allc(path, chrom, positions_0based, out$cells[[i]])
        cell_ids <- c(cell_ids, id)
        groups <- c(groups, sprintf("K%d", K))
        paths <- c(paths, path)
    }
}

write_manifest(file.path(opt$out_dir, "manifest.tsv"), cell_ids, groups, paths)
message(sprintf("[generate_mixture_k_sweep] wrote %d cells across K = %s",
                length(cell_ids), opt$k_grid))
