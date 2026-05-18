## Loader for the windows HDF5 store built by build_windows_h5.R.
##
## Returns HDF5-backed (DelayedArray) windows x cells matrices for i_total and
## meth, with dimnames set and cells absent from amet's cell_feature output
## dropped. The matrices stay on disk and are read block by block; the store
## exists to avoid reading amet's billion-row long table. Missing values are
## stored as NaN (the HDF5 fill value and genuine missing coverage alike),
## which is.na() and na.rm handling treat as NA.

load_windows_h5 <- function(h5_path) {
  suppressPackageStartupMessages({
    library(rhdf5)
    library(HDF5Array)
  })
  feature_id <- as.character(h5read(h5_path, "feature_id"))
  cell_id <- as.character(h5read(h5_path, "cell_id"))
  keep <- as.logical(h5read(h5_path, "cell_present"))
  i_total <- HDF5Array(h5_path, "i_total")
  meth <- HDF5Array(h5_path, "meth")
  dimnames(i_total) <- list(feature_id, cell_id)
  dimnames(meth) <- list(feature_id, cell_id)
  list(
    i_total = i_total[, keep, drop = FALSE],
    meth = meth[, keep, drop = FALSE],
    feature_id = feature_id,
    cell_id = cell_id[keep]
  )
}
