## across-cell Variability Index (acVI), integer 1..10.
##
## Uses the same Markov-with-repeat process as wcVI (so the two simulators are
## comparable). The difference is what scales with the index:
##   - wcVI: p_repeat varies WITH the index (per cell, all cells in a group share).
##   - acVI: p_repeat varies WITHIN the group; the per-cell spread scales with
##           acVI. acVI=1 → all cells share p_repeat=0.725; acVI=10 → cells span
##           p_repeat uniformly across [0.5, 0.95].
##
## Each cell is an independent realisation, marginal pinned to target_p exactly
## by construction. target_p is swept 0.05..0.95 in 0.05 steps so the simulation
## covers the full DNA methylation range at every acVI level.

suppressPackageStartupMessages({ library(optparse) })
.this_dir <- local({ args <- commandArgs(trailingOnly = FALSE); fa <- grep("^--file=", args, value = TRUE); if (length(fa) > 0) dirname(sub("^--file=", "", fa[1])) else "." })
source(file.path(.this_dir, "write_outputs.R"))
source(file.path(.this_dir, "bricks.R"))

options <- list(
    make_option(c("--out_dir"), type = "character"),
    make_option(c("--n_cells_per_group"), type = "integer", default = 30),
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

p_centre <- 0.725
p_half_max <- 0.225  # so range is [0.5, 0.95] at acvi=10

for (acvi in 1:10) {
    half_width <- (acvi - 1L) * p_half_max / 9.0
    for (p in marginals) {
        for (i in seq_len(opt$n_cells_per_group)) {
            p_repeat_cell <- p_centre + runif(1, -half_width, half_width)
            p_repeat_cell <- max(0.0, min(0.95, p_repeat_cell))
            bits <- simulate_repeat_cell(opt$n_cpgs, p, p_repeat_cell)
            id <- sprintf("acvi%02d_p%.2f_cell%03d", acvi, p, i)
            path <- file.path(cells_dir, sprintf("%s.allc.tsv.gz", id))
            write_allc(path, chrom, positions_0based, bits)
            cell_ids <- c(cell_ids, id)
            groups <- c(groups, sprintf("acvi%02d_p%.2f", acvi, p))
            paths <- c(paths, path)
        }
    }
}

write_manifest(file.path(opt$out_dir, "manifest.tsv"), cell_ids, groups, paths)
message(sprintf("[generate_acvi_sweep] wrote %d cells across %d marginals",
                length(cell_ids), length(marginals)))
