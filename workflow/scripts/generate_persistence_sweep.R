## Sweep comethylation strength at fixed marginal p = 0.5.
## Each persistence c sets p00 = p11 = c. c = 0.5 is i.i.d.; c > 0.5 is comethylated;
## c < 0.5 is anti-persistent.

suppressPackageStartupMessages({ library(optparse) })
source(file.path(dirname(sys.frame(1)$ofile), "write_outputs.R"))
source(file.path(dirname(sys.frame(1)$ofile), "simulate_markov_chain.R"))

options <- list(
    make_option(c("--out_dir"), type = "character"),
    make_option(c("--n_cells_per_value"), type = "integer", default = 20),
    make_option(c("--n_cpgs"), type = "integer", default = 200),
    make_option(c("--persistence_grid"), type = "character",
                default = "0.30,0.40,0.50,0.60,0.70,0.80,0.90"),
    make_option(c("--seed"), type = "integer", default = 1)
)
opt <- parse_args(OptionParser(option_list = options))

dir.create(opt$out_dir, showWarnings = FALSE, recursive = TRUE)
cells_dir <- file.path(opt$out_dir, "cells")
dir.create(cells_dir, showWarnings = FALSE)

cs <- as.numeric(strsplit(opt$persistence_grid, ",")[[1]])
chrom <- "chr1"
positions_0based <- 99L + (seq_len(opt$n_cpgs) - 1L) * 100L
write_cpg_reference(file.path(opt$out_dir, "cpgs.tsv"), chrom, positions_0based)
write_bed_feature(file.path(opt$out_dir, "feature.bed"), chrom, 0L,
                  positions_0based[length(positions_0based)] + 100L, "region")

set.seed(opt$seed)
cell_ids <- character(); groups <- character(); paths <- character()
truth <- data.frame(condition = character(), persistence = numeric(),
                    target_pi1 = numeric(), true_i1 = numeric())

for (c in cs) {
    p00 <- c; p11 <- c
    for (i in seq_len(opt$n_cells_per_value)) {
        bits <- simulate_markov_cell(opt$n_cpgs, p00, p11)
        id <- sprintf("c%.2f_seed%d", c, i)
        path <- file.path(cells_dir, sprintf("%s.allc.tsv.gz", id))
        write_allc(path, chrom, positions_0based, bits)
        cell_ids <- c(cell_ids, id)
        groups <- c(groups, sprintf("c%.2f", c))
        paths <- c(paths, path)
    }
    truth <- rbind(truth, data.frame(
        condition = sprintf("c%.2f", c),
        persistence = c,
        target_pi1 = stationary_marginal(p00, p11),
        true_i1 = markov_i1(p00, p11)
    ))
}

write_manifest(file.path(opt$out_dir, "manifest.tsv"), cell_ids, groups, paths)
write.table(truth, file.path(opt$out_dir, "ground_truth.csv"),
            sep = ",", row.names = FALSE, quote = FALSE)
message(sprintf("[generate_persistence_sweep] wrote %d cells across %d values",
                length(cell_ids), length(cs)))
