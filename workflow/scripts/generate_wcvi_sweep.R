## within-cell Variability Index (wcVI), integer 1..10.
##
## Each cell is an independent realisation of a Markov-with-repeat process at
## marginal target_p. wcVI sets the repeat probability:
##   wcVI = 1  -> p_repeat = 0.95 -> long runs of the same state, very predictable.
##   wcVI = 10 -> p_repeat = 0    -> pure iid Bernoulli at target_p (no within-cell
##                                  structure).
##
## Across-cell variability inside a (wcVI, target_p) group is just finite-sample
## noise from independent draws - fully decoupled from acVI. No bricks are used,
## so increasing wcVI does not enlarge any pool that cells share.

suppressPackageStartupMessages({ library(optparse) })
.this_dir <- local({ args <- commandArgs(trailingOnly = FALSE); fa <- grep("^--file=", args, value = TRUE); if (length(fa) > 0) dirname(sub("^--file=", "", fa[1])) else "." })
source(file.path(.this_dir, "write_outputs.R"))
source(file.path(.this_dir, "bricks.R"))

options <- list(
    make_option(c("--out_dir"), type = "character"),
    make_option(c("--n_cells_per_combo"), type = "integer", default = 30),
    make_option(c("--n_cpgs"), type = "integer", default = 2000),
    make_option(c("--marginals"), type = "character",
                default = "0.05,0.10,0.15,0.20,0.25,0.30,0.35,0.40,0.45,0.50,0.55,0.60,0.65,0.70,0.75,0.80,0.85,0.90,0.95"),
    make_option(c("--seed"), type = "integer", default = 1)
)
opt <- parse_args(OptionParser(option_list = options))

dir.create(opt$out_dir, showWarnings = FALSE, recursive = TRUE)
cells_dir <- file.path(opt$out_dir, "cells")
dir.create(cells_dir, showWarnings = FALSE)

marginals <- as.numeric(strsplit(opt$marginals, ",")[[1]])

chrom <- "chr1"
positions_0based <- 99L + (seq_len(opt$n_cpgs) - 1L) * 100L
write_cpg_reference(file.path(opt$out_dir, "cpgs.tsv"), chrom, positions_0based)
write_bed_feature(file.path(opt$out_dir, "feature.bed"), chrom, 0L,
                  positions_0based[length(positions_0based)] + 100L, "region")

set.seed(opt$seed)
cell_ids <- character(); groups <- character(); paths <- character()

wcvi_to_p_repeat <- function(v) (10L - v) * 0.95 / 9.0

for (wcvi in 1:10) {
    p_repeat <- wcvi_to_p_repeat(wcvi)
    for (p in marginals) {
        for (i in seq_len(opt$n_cells_per_combo)) {
            bits <- simulate_repeat_cell(opt$n_cpgs, p, p_repeat)
            id <- sprintf("wcvi%02d_p%.2f_cell%03d", wcvi, p, i)
            path <- file.path(cells_dir, sprintf("%s.allc.tsv.gz", id))
            write_allc(path, chrom, positions_0based, bits)
            cell_ids <- c(cell_ids, id)
            groups <- c(groups, sprintf("wcvi%02d_p%.2f", wcvi, p))
            paths <- c(paths, path)
        }
    }
}

write_manifest(file.path(opt$out_dir, "manifest.tsv"), cell_ids, groups, paths)
message(sprintf("[generate_wcvi_sweep] wrote %d cells across %d marginals",
                length(cell_ids), length(marginals)))
