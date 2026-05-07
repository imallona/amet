## Sweep feature length (number of CpGs in the simulated region) at fixed structure.

suppressPackageStartupMessages({ library(optparse) })
source(file.path(dirname(sys.frame(1)$ofile), "write_outputs.R"))
source(file.path(dirname(sys.frame(1)$ofile), "simulate_markov_chain.R"))

options <- list(
    make_option(c("--out_dir"), type = "character"),
    make_option(c("--n_cells_per_length"), type = "integer", default = 20),
    make_option(c("--length_grid"), type = "character",
                default = "20,50,100,200,500,1000"),
    make_option(c("--persistence"), type = "double", default = 0.85),
    make_option(c("--seed"), type = "integer", default = 1)
)
opt <- parse_args(OptionParser(option_list = options))

dir.create(opt$out_dir, showWarnings = FALSE, recursive = TRUE)
cells_dir <- file.path(opt$out_dir, "cells")
dir.create(cells_dir, showWarnings = FALSE)

lengths <- as.integer(strsplit(opt$length_grid, ",")[[1]])
chrom <- "chr1"
max_len <- max(lengths)
positions_0based <- 99L + (seq_len(max_len) - 1L) * 100L
write_cpg_reference(file.path(opt$out_dir, "cpgs.tsv"), chrom, positions_0based)

## One feature per length, all sharing the same start.
bed_rows <- character()
for (l in lengths) {
    end <- positions_0based[l] + 100L
    bed_rows <- c(bed_rows,
                  sprintf("%s\t0\t%d\tlen%d", chrom, end, l))
}
writeLines(bed_rows, file.path(opt$out_dir, "feature.bed"))

set.seed(opt$seed)
cell_ids <- character(); groups <- character(); paths <- character()

for (l in lengths) {
    for (i in seq_len(opt$n_cells_per_length)) {
        bits <- simulate_markov_cell(max_len, opt$persistence, opt$persistence)
        id <- sprintf("len%d_seed%d", l, i)
        path <- file.path(cells_dir, sprintf("%s.allc.tsv.gz", id))
        write_allc(path, chrom, positions_0based, bits)
        cell_ids <- c(cell_ids, id)
        groups <- c(groups, sprintf("len%d", l))
        paths <- c(paths, path)
    }
}

write_manifest(file.path(opt$out_dir, "manifest.tsv"), cell_ids, groups, paths)
message(sprintf("[generate_feature_length_sweep] wrote %d cells, %d feature lengths",
                length(cell_ids), length(lengths)))
