## Sweep mixture-component divergence at fixed K = 2.
## Two components: c1 = 0.5 (i.i.d.) and c2 = 0.5 + delta. Larger delta means more
## divergent components, so JSD across cells should increase with delta.

suppressPackageStartupMessages({ library(optparse) })
.this_dir <- local({ args <- commandArgs(trailingOnly = FALSE); fa <- grep("^--file=", args, value = TRUE); if (length(fa) > 0) dirname(sub("^--file=", "", fa[1])) else "." })

source(file.path(.this_dir, "write_outputs.R"))
source(file.path(.this_dir, "simulate_markov_chain.R"))
source(file.path(.this_dir, "mix_markovs.R"))

options <- list(
    make_option(c("--out_dir"), type = "character"),
    make_option(c("--n_cells_per_value"), type = "integer", default = 60),
    make_option(c("--n_cpgs"), type = "integer", default = 200),
    make_option(c("--delta_grid"), type = "character",
                default = "0.0,0.1,0.2,0.3,0.4"),
    make_option(c("--seed"), type = "integer", default = 1)
)
opt <- parse_args(OptionParser(option_list = options))

dir.create(opt$out_dir, showWarnings = FALSE, recursive = TRUE)
cells_dir <- file.path(opt$out_dir, "cells")
dir.create(cells_dir, showWarnings = FALSE)

deltas <- as.numeric(strsplit(opt$delta_grid, ",")[[1]])
chrom <- "chr1"
positions_0based <- 99L + (seq_len(opt$n_cpgs) - 1L) * 100L
write_cpg_reference(file.path(opt$out_dir, "cpgs.tsv"), chrom, positions_0based)
write_bed_feature(file.path(opt$out_dir, "feature.bed"), chrom, 0L,
                  positions_0based[length(positions_0based)] + 100L, "region")

set.seed(opt$seed)
cell_ids <- character(); groups <- character(); paths <- character()

for (delta in deltas) {
    components <- list(
        list(p00 = 0.5, p11 = 0.5),
        list(p00 = 0.5 + delta, p11 = 0.5 + delta)
    )
    out <- simulate_mixture_cells(
        opt$n_cells_per_value,
        opt$n_cpgs,
        components,
        seed = opt$seed + as.integer(round(delta * 100))
    )
    for (i in seq_along(out$cells)) {
        id <- sprintf("d%.2f_cell%d", delta, i)
        path <- file.path(cells_dir, sprintf("%s.allc.tsv.gz", id))
        write_allc(path, chrom, positions_0based, out$cells[[i]])
        cell_ids <- c(cell_ids, id)
        groups <- c(groups, sprintf("d%.2f", delta))
        paths <- c(paths, path)
    }
}

write_manifest(file.path(opt$out_dir, "manifest.tsv"), cell_ids, groups, paths)
message(sprintf("[generate_mixture_divergence_sweep] wrote %d cells across %d divergences",
                length(cell_ids), length(deltas)))
