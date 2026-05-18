#!/usr/bin/env Rscript

## Pivot amet's long windows cell_feature TSV into an HDF5 store of
## windows x cells matrices.
##
## amet writes one row per (cell, window). For a genome-wide windows run over
## thousands of cells that long table is too large to read whole. amet writes
## each cell's rows as one contiguous block with a fixed window order (see
## method/src/main.rs), so this script streams the gzipped input cell block by
## cell block and writes each cell as one column of the output matrices. Peak
## memory stays near one cell plus one read batch.
##
## Output HDF5 datasets:
##   i_total      double matrix, windows x cells
##   meth         double matrix, windows x cells (amet mean_meth)
##   feature_id   character vector, length windows
##   cell_id      character vector, length cells (manifest order)
##   cell_present integer vector, 1 if the cell was found in the input
## Cells absent from the cell_feature input keep an all-NaN column.

suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
  library(rhdf5)
})

opt <- parse_args(OptionParser(option_list = list(
  make_option("--cell-feature", dest = "cell_feature",
              help = "amet windows cell_feature TSV (gzipped)"),
  make_option("--manifest", dest = "manifest",
              help = "cells.tsv manifest with a cell_id column"),
  make_option("--output", dest = "output",
              help = "output HDF5 path"),
  make_option("--batch-lines", dest = "batch_lines", type = "integer",
              default = 2000000L,
              help = "input lines read per streaming batch"),
  make_option("--threads", dest = "threads", type = "integer", default = 4L)
)))

setDTthreads(opt$threads)

manifest <- fread(opt$manifest)
if (!"cell_id" %in% colnames(manifest))
  stop("manifest has no cell_id column")
manifest_cells <- as.character(manifest$cell_id)
n_cells <- length(manifest_cells)
cell_to_col <- setNames(seq_len(n_cells), manifest_cells)
cat("Manifest cells:", n_cells, "\n")

if (file.exists(opt$output)) file.remove(opt$output)
h5createFile(opt$output)

con <- pipe(sprintf("zcat %s", shQuote(opt$cell_feature)), open = "rt")
on.exit(close(con), add = TRUE)

header <- readLines(con, n = 1L)
cols <- strsplit(header, "\t", fixed = TRUE)[[1]]
need <- c("cell_id", "feature_id", "i_total", "mean_meth")
pos <- match(need, cols)
if (anyNA(pos))
  stop("cell_feature missing required columns: ",
       paste(need[is.na(pos)], collapse = ", "))

n_windows <- NA_integer_
feature_ids <- NULL
seen_cols <- integer(0)

flush_cell <- function(dc) {
  cid <- dc$cell_id[1]
  if (is.na(n_windows)) {
    n_windows <<- nrow(dc)
    feature_ids <<- dc$feature_id
    h5createDataset(opt$output, "i_total", dims = c(n_windows, n_cells),
                    storage.mode = "double", chunk = c(n_windows, 1L),
                    level = 4L, fillValue = NaN)
    h5createDataset(opt$output, "meth", dims = c(n_windows, n_cells),
                    storage.mode = "double", chunk = c(n_windows, 1L),
                    level = 4L, fillValue = NaN)
    cat("Windows per cell:", n_windows, "\n")
  } else if (nrow(dc) != n_windows) {
    stop(sprintf("cell %s has %d rows, expected %d", cid, nrow(dc), n_windows))
  } else if (!identical(dc$feature_id, feature_ids)) {
    stop(sprintf("cell %s window order differs from the first cell", cid))
  }
  j <- unname(cell_to_col[cid])
  if (is.na(j)) {
    warning("cell not in manifest, skipped: ", cid)
    return(invisible())
  }
  h5write(matrix(dc$i_total, ncol = 1L), opt$output, "i_total",
          index = list(seq_len(n_windows), j))
  h5write(matrix(dc$meth, ncol = 1L), opt$output, "meth",
          index = list(seq_len(n_windows), j))
  seen_cols <<- c(seen_cols, j)
}

pending <- NULL
repeat {
  lines <- readLines(con, n = opt$batch_lines)
  if (length(lines) == 0L) break
  dt <- fread(text = lines, header = FALSE, sep = "\t")
  dt <- dt[, ..pos]
  setnames(dt, need)
  pending <- if (is.null(pending)) dt else rbindlist(list(pending, dt))
  ## The last cell in pending may continue into the next batch; everything
  ## before it is complete because amet writes each cell as one block.
  last_cid <- pending$cell_id[nrow(pending)]
  complete <- pending[cell_id != last_cid]
  pending <- pending[cell_id == last_cid]
  for (cid in unique(complete$cell_id))
    flush_cell(complete[cell_id == cid])
}
if (!is.null(pending) && nrow(pending) > 0L)
  for (cid in unique(pending$cell_id))
    flush_cell(pending[cell_id == cid])

if (is.na(n_windows)) stop("no cells found in cell_feature input")

cell_present <- integer(n_cells)
cell_present[unique(seen_cols)] <- 1L

h5write(feature_ids, opt$output, "feature_id")
h5write(manifest_cells, opt$output, "cell_id")
h5write(cell_present, opt$output, "cell_present")
H5close()

cat("Cells written:", length(unique(seen_cols)), "/", n_cells, "\n")
cat("Output:", opt$output, "\n")
