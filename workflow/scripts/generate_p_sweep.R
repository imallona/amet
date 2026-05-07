## Generate one simulation experiment: a set of cells with controlled marginal p
## and comethylation, written as allc files plus a manifest, CpG reference, and
## features BED.
##
## Two cell types per condition:
##   - iid: p00 + p11 = 1 (no comethylation)
##   - structured: persistent chain at the same marginal
##
## All files are written under <out_dir>/.

suppressPackageStartupMessages({
    library(optparse)
})

.this_dir <- local({ args <- commandArgs(trailingOnly = FALSE); fa <- grep("^--file=", args, value = TRUE); if (length(fa) > 0) dirname(sub("^--file=", "", fa[1])) else "." })

source(file.path(.this_dir, "write_outputs.R"))
source(file.path(.this_dir, "simulate_markov_chain.R"))

options <- list(
    make_option(c("--out_dir"), type = "character"),
    make_option(c("--n_cells_per_condition"), type = "integer", default = 20),
    make_option(c("--n_cpgs"), type = "integer", default = 200),
    make_option(c("--persistence"), type = "double", default = 0.85,
                help = "Diagonal value for the structured chains (p00=p11=this)"),
    make_option(c("--p_grid"), type = "character",
                default = "0.1,0.3,0.5,0.7,0.9",
                help = "Comma-separated marginal methylation values to sweep"),
    make_option(c("--seed"), type = "integer", default = 1)
)
opt <- parse_args(OptionParser(option_list = options))

dir.create(opt$out_dir, showWarnings = FALSE, recursive = TRUE)
cells_dir <- file.path(opt$out_dir, "cells")
dir.create(cells_dir, showWarnings = FALSE)

p_values <- as.numeric(strsplit(opt$p_grid, ",")[[1]])
chrom <- "chr1"
cpg_step <- 100L
positions_0based <- 99L + (seq_len(opt$n_cpgs) - 1L) * cpg_step

write_cpg_reference(file.path(opt$out_dir, "cpgs.tsv"), chrom, positions_0based)
feat_end <- positions_0based[length(positions_0based)] + 100L
write_bed_feature(file.path(opt$out_dir, "feature.bed"),
                  chrom, 0L, feat_end, "synthetic_region")

set.seed(opt$seed)
cell_ids <- character()
groups <- character()
paths <- character()

for (p in p_values) {
    ## i.i.d.: pi_1 = p, with p00 = 1 - p and p11 = p (so p00 + p11 = 1).
    iid_p00 <- 1 - p
    iid_p11 <- p
    ## structured: same marginal pi_1 = p, with p00 + p11 = 2*persistence,
    ## solving (1 - p00) / (2 - p00 - p11) = p with p00 + p11 = 2*opt$persistence.
    sum_pp <- 2 * opt$persistence
    s_p00 <- 1 - p * (2 - sum_pp)
    s_p11 <- sum_pp - s_p00
    if (s_p00 < 0 || s_p00 > 1 || s_p11 < 0 || s_p11 > 1) {
        message(sprintf("[skip] cannot reach pi_1=%.2f at persistence=%.2f", p, opt$persistence))
        next
    }

    for (i in seq_len(opt$n_cells_per_condition)) {
        for (kind in c("iid", "structured")) {
            if (kind == "iid") {
                bits <- simulate_markov_cell(opt$n_cpgs, iid_p00, iid_p11)
            } else {
                bits <- simulate_markov_cell(opt$n_cpgs, s_p00, s_p11)
            }
            id <- sprintf("%s_p%.2f_seed%d", kind, p, i)
            path <- file.path(cells_dir, sprintf("%s.allc.tsv.gz", id))
            write_allc(path, chrom, positions_0based, bits)
            cell_ids <- c(cell_ids, id)
            groups <- c(groups, sprintf("%s_p%.2f", kind, p))
            paths <- c(paths, path)
        }
    }
}

write_manifest(file.path(opt$out_dir, "manifest.tsv"),
               cell_ids, groups, paths)

message(sprintf("[generate] wrote %d cells across %d marginal values to %s",
                length(cell_ids), length(p_values), opt$out_dir))
