## Sweep sequencing depth at fixed structure across the full methylation range.
## Each cell is a Markov-with-repeat process at marginal target_p with moderate
## structure (p_repeat = 0.7), then a fragment-based coverage model drops random
## CpGs. Sweeps both n_fragments (coverage) and target_p (methylation), so the
## eval can show I_total response to coverage AND methylation in one figure.

suppressPackageStartupMessages({ library(optparse) })
.this_dir <- local({ args <- commandArgs(trailingOnly = FALSE); fa <- grep("^--file=", args, value = TRUE); if (length(fa) > 0) dirname(sub("^--file=", "", fa[1])) else "." })

source(file.path(.this_dir, "write_outputs.R"))
source(file.path(.this_dir, "bricks.R"))
source(file.path(.this_dir, "apply_fragment_sparsity.R"))

options <- list(
    make_option(c("--out_dir"), type = "character"),
    make_option(c("--n_cells_per_combo"), type = "integer", default = 10),
    make_option(c("--n_cpgs"), type = "integer", default = 1000),
    make_option(c("--p_repeat"), type = "double", default = 0.7),
    make_option(c("--n_fragments_grid"), type = "character",
                default = "10,30,60,100,200,500"),
    make_option(c("--marginals"), type = "character",
                default = "0.10,0.30,0.50,0.70,0.90"),
    make_option(c("--fragment_length_mean"), type = "double", default = 350),
    make_option(c("--fragment_length_sd"), type = "double", default = 150),
    make_option(c("--seed"), type = "integer", default = 1)
)
opt <- parse_args(OptionParser(option_list = options))

dir.create(opt$out_dir, showWarnings = FALSE, recursive = TRUE)
cells_dir <- file.path(opt$out_dir, "cells")
dir.create(cells_dir, showWarnings = FALSE)

n_frag_values <- as.integer(strsplit(opt$n_fragments_grid, ",")[[1]])
marginals <- as.numeric(strsplit(opt$marginals, ",")[[1]])
chrom <- "chr1"
positions_0based <- 99L + (seq_len(opt$n_cpgs) - 1L) * 100L
write_cpg_reference(file.path(opt$out_dir, "cpgs.tsv"), chrom, positions_0based)
write_bed_feature(file.path(opt$out_dir, "feature.bed"), chrom, 0L,
                  positions_0based[length(positions_0based)] + 100L, "region")

set.seed(opt$seed)
cell_ids <- character(); groups <- character(); paths <- character()

for (nf in n_frag_values) {
    for (p in marginals) {
        for (i in seq_len(opt$n_cells_per_combo)) {
            bits <- simulate_repeat_cell(opt$n_cpgs, p, opt$p_repeat)
            kept <- fragment_coverage(positions_0based,
                                      n_fragments = nf,
                                      fragment_length_mean = opt$fragment_length_mean,
                                      fragment_length_sd = opt$fragment_length_sd)
            idx <- kept$indices
            id <- sprintf("nf%d_p%.2f_seed%d", nf, p, i)
            path <- file.path(cells_dir, sprintf("%s.allc.tsv.gz", id))
            if (length(idx) == 0L) {
                write_allc(path, chrom, integer(0), integer(0))
            } else {
                write_allc(path, chrom, positions_0based[idx], bits[idx])
            }
            cell_ids <- c(cell_ids, id)
            groups <- c(groups, sprintf("nf%d_p%.2f", nf, p))
            paths <- c(paths, path)
        }
    }
}

write_manifest(file.path(opt$out_dir, "manifest.tsv"), cell_ids, groups, paths)
message(sprintf("[generate_sparsity_sweep] wrote %d cells, %d coverage levels x %d marginals",
                length(cell_ids), length(n_frag_values), length(marginals)))
