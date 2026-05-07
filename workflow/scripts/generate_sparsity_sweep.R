## Sweep sequencing depth at fixed structure. Each cell is generated as a structured
## Markov chain at p = 0.5, then a fragment-based coverage model is applied with a
## variable number of fragments.

suppressPackageStartupMessages({ library(optparse) })
source(file.path(dirname(sys.frame(1)$ofile), "write_outputs.R"))
source(file.path(dirname(sys.frame(1)$ofile), "simulate_markov_chain.R"))
source(file.path(dirname(sys.frame(1)$ofile), "apply_fragment_sparsity.R"))

options <- list(
    make_option(c("--out_dir"), type = "character"),
    make_option(c("--n_cells_per_value"), type = "integer", default = 20),
    make_option(c("--n_cpgs"), type = "integer", default = 500),
    make_option(c("--persistence"), type = "double", default = 0.85),
    make_option(c("--n_fragments_grid"), type = "character",
                default = "10,30,60,100,200,500"),
    make_option(c("--fragment_length_mean"), type = "double", default = 350),
    make_option(c("--fragment_length_sd"), type = "double", default = 150),
    make_option(c("--seed"), type = "integer", default = 1)
)
opt <- parse_args(OptionParser(option_list = options))

dir.create(opt$out_dir, showWarnings = FALSE, recursive = TRUE)
cells_dir <- file.path(opt$out_dir, "cells")
dir.create(cells_dir, showWarnings = FALSE)

n_frag_values <- as.integer(strsplit(opt$n_fragments_grid, ",")[[1]])
chrom <- "chr1"
positions_0based <- 99L + (seq_len(opt$n_cpgs) - 1L) * 100L
write_cpg_reference(file.path(opt$out_dir, "cpgs.tsv"), chrom, positions_0based)
write_bed_feature(file.path(opt$out_dir, "feature.bed"), chrom, 0L,
                  positions_0based[length(positions_0based)] + 100L, "region")

set.seed(opt$seed)
cell_ids <- character(); groups <- character(); paths <- character()

for (nf in n_frag_values) {
    for (i in seq_len(opt$n_cells_per_value)) {
        bits <- simulate_markov_cell(opt$n_cpgs, opt$persistence, opt$persistence)
        kept <- fragment_coverage(positions_0based,
                                  n_fragments = nf,
                                  fragment_length_mean = opt$fragment_length_mean,
                                  fragment_length_sd = opt$fragment_length_sd)
        idx <- kept$indices
        id <- sprintf("nf%d_seed%d", nf, i)
        path <- file.path(cells_dir, sprintf("%s.allc.tsv.gz", id))
        if (length(idx) == 0L) {
            write_allc(path, chrom, integer(0), integer(0))
        } else {
            write_allc(path, chrom, positions_0based[idx], bits[idx])
        }
        cell_ids <- c(cell_ids, id)
        groups <- c(groups, sprintf("nf%d", nf))
        paths <- c(paths, path)
    }
}

write_manifest(file.path(opt$out_dir, "manifest.tsv"), cell_ids, groups, paths)
message(sprintf("[generate_sparsity_sweep] wrote %d cells across %d coverage levels",
                length(cell_ids), length(n_frag_values)))
