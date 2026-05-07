## Sweep feature length and target marginal jointly. All cells are generated from
## a Markov-with-repeat process at marginal target_p with moderate structure
## (p_repeat = 0.7), full length (max_len). Each feature in the BED corresponds
## to one length value; the per-cell methylation lookup picks up the appropriate
## subset for each feature.

suppressPackageStartupMessages({ library(optparse) })
.this_dir <- local({ args <- commandArgs(trailingOnly = FALSE); fa <- grep("^--file=", args, value = TRUE); if (length(fa) > 0) dirname(sub("^--file=", "", fa[1])) else "." })

source(file.path(.this_dir, "write_outputs.R"))
source(file.path(.this_dir, "bricks.R"))

options <- list(
    make_option(c("--out_dir"), type = "character"),
    make_option(c("--n_cells_per_combo"), type = "integer", default = 10),
    make_option(c("--length_grid"), type = "character",
                default = "20,50,100,200,500,1000"),
    make_option(c("--p_repeat"), type = "double", default = 0.7),
    make_option(c("--marginals"), type = "character",
                default = "0.10,0.30,0.50,0.70,0.90"),
    make_option(c("--seed"), type = "integer", default = 1)
)
opt <- parse_args(OptionParser(option_list = options))

dir.create(opt$out_dir, showWarnings = FALSE, recursive = TRUE)
cells_dir <- file.path(opt$out_dir, "cells")
dir.create(cells_dir, showWarnings = FALSE)

lengths <- as.integer(strsplit(opt$length_grid, ",")[[1]])
marginals <- as.numeric(strsplit(opt$marginals, ",")[[1]])
chrom <- "chr1"
max_len <- max(lengths)
positions_0based <- 99L + (seq_len(max_len) - 1L) * 100L
write_cpg_reference(file.path(opt$out_dir, "cpgs.tsv"), chrom, positions_0based)

bed_rows <- character()
for (l in lengths) {
    end <- positions_0based[l] + 100L
    bed_rows <- c(bed_rows, sprintf("%s\t0\t%d\tlen%d", chrom, end, l))
}
writeLines(bed_rows, file.path(opt$out_dir, "feature.bed"))

set.seed(opt$seed)
cell_ids <- character(); groups <- character(); paths <- character()

for (p in marginals) {
    for (i in seq_len(opt$n_cells_per_combo)) {
        bits <- simulate_repeat_cell(max_len, p, opt$p_repeat)
        id <- sprintf("p%.2f_seed%d", p, i)
        path <- file.path(cells_dir, sprintf("%s.allc.tsv.gz", id))
        write_allc(path, chrom, positions_0based, bits)
        cell_ids <- c(cell_ids, id)
        groups <- c(groups, sprintf("p%.2f", p))
        paths <- c(paths, path)
    }
}

write_manifest(file.path(opt$out_dir, "manifest.tsv"), cell_ids, groups, paths)
message(sprintf("[generate_feature_length_sweep] wrote %d cells across %d marginals (each scored at %d feature lengths)",
                length(cell_ids), length(marginals), length(lengths)))
