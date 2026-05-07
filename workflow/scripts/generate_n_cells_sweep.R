## Sweep number of cells per group at fixed mixture (K = 2 with mid divergence).
## Used to characterise JSD finite-sample bias and variance.

suppressPackageStartupMessages({ library(optparse) })
.this_dir <- local({ args <- commandArgs(trailingOnly = FALSE); fa <- grep("^--file=", args, value = TRUE); if (length(fa) > 0) dirname(sub("^--file=", "", fa[1])) else "." })

source(file.path(.this_dir, "write_outputs.R"))
source(file.path(.this_dir, "simulate_markov_chain.R"))
source(file.path(.this_dir, "mix_markovs.R"))

options <- list(
    make_option(c("--out_dir"), type = "character"),
    make_option(c("--n_cells_grid"), type = "character",
                default = "5,10,20,50,100,200"),
    make_option(c("--n_cpgs"), type = "integer", default = 200),
    make_option(c("--seed"), type = "integer", default = 1)
)
opt <- parse_args(OptionParser(option_list = options))

dir.create(opt$out_dir, showWarnings = FALSE, recursive = TRUE)
cells_dir <- file.path(opt$out_dir, "cells")
dir.create(cells_dir, showWarnings = FALSE)

ns <- as.integer(strsplit(opt$n_cells_grid, ",")[[1]])
chrom <- "chr1"
positions_0based <- 99L + (seq_len(opt$n_cpgs) - 1L) * 100L
write_cpg_reference(file.path(opt$out_dir, "cpgs.tsv"), chrom, positions_0based)
write_bed_feature(file.path(opt$out_dir, "feature.bed"), chrom, 0L,
                  positions_0based[length(positions_0based)] + 100L, "region")

components <- list(list(p00 = 0.5, p11 = 0.5),
                   list(p00 = 0.8, p11 = 0.8))

set.seed(opt$seed)
cell_ids <- character(); groups <- character(); paths <- character()

for (n in ns) {
    out <- simulate_mixture_cells(n, opt$n_cpgs, components)
    for (i in seq_along(out$cells)) {
        id <- sprintf("n%d_cell%d", n, i)
        path <- file.path(cells_dir, sprintf("%s.allc.tsv.gz", id))
        write_allc(path, chrom, positions_0based, out$cells[[i]])
        cell_ids <- c(cell_ids, id)
        groups <- c(groups, sprintf("n%d", n))
        paths <- c(paths, path)
    }
}

write_manifest(file.path(opt$out_dir, "manifest.tsv"), cell_ids, groups, paths)
message(sprintf("[generate_n_cells_sweep] wrote %d cells across %d group sizes",
                length(cell_ids), length(ns)))
