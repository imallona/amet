## Shared render-time logging and threading for analytical and figure Rmds.
## Called from a `logging_early` chunk at the top of each Rmd.
##
## - `log_path = ""` (or NA) skips file sinking entirely; logs go to the
##   knit console only.
## - Re-entrant: if called twice in the same R process, the second call
##   does not stack another sink on top of the first.
## - `threads` (default 1) caps data.table OpenMP threads and is the worker
##   count for the returned BiocParallel BPPARAM. Set to >1 only with memory
##   headroom -- MulticoreParam forks the parent process for each worker.

.amet_log_state <- new.env(parent = emptyenv())

amet_parse_threads <- function(x, default = 1L) {
  if (is.null(x)) return(default)
  if (is.character(x) && !nzchar(x)) return(default)
  n <- suppressWarnings(as.integer(x))
  if (is.na(n) || n < 1L) default else n
}

amet_make_bpparam <- function(threads) {
  threads <- amet_parse_threads(threads, 1L)
  if (threads <= 1L) BiocParallel::SerialParam()
  else BiocParallel::MulticoreParam(workers = threads)
}

## Shannon entropy of a binary distribution in bits. Vectorized; returns NA
## at p in {0, 1, NA}.
shannon_binary <- function(p) {
  out <- numeric(length(p))
  safe <- !is.na(p) & p > 0 & p < 1
  out[safe] <- -p[safe] * log2(p[safe]) - (1 - p[safe]) * log2(1 - p[safe])
  out[!safe] <- NA_real_
  out
}

## Canonical i_norm: i_total normalised by its marginal-entropy ceiling.
## i_total / (k_max * H(p_hat)). Matches the formula used in the eval scripts
## and simulations_report.Rmd; k_max is amet's --i-max-lag.
compute_i_norm <- function(i_total, mean_meth, k_max) {
  k_max <- amet_parse_threads(k_max, 1L)
  denom <- k_max * shannon_binary(mean_meth)
  out <- i_total / denom
  out[!is.finite(out)] <- NA_real_
  out
}

## Load a windows_annotation.tsv.gz produced by the snakemake helper rules and
## attach annotation columns to a data.table keyed by feature_id. The TSV has
## chrom/start/end/feature_id followed by one numeric coverage-fraction column
## per annotation. Returns a data.table with feature_id and the per-annotation
## columns, or NULL if path is empty/missing. Caller is responsible for
## merging into a SCE / wide data.frame; the helper just handles I/O.
amet_load_annotation_matrix <- function(path) {
  if (is.null(path) || is.na(path) || !nzchar(path)) return(NULL)
  if (!file.exists(path)) {
    message("[annotation] file not found: ", path)
    return(NULL)
  }
  if (!requireNamespace("data.table", quietly = TRUE))
    stop("amet_load_annotation_matrix requires data.table")
  m <- data.table::fread(path)
  if (nrow(m) == 0L) return(NULL)
  expected <- c("chrom", "start", "end", "feature_id")
  if (!all(expected %in% names(m))) {
    message("[annotation] missing chrom/start/end/feature_id columns; got: ",
            paste(names(m), collapse = ", "))
    return(NULL)
  }
  ann_cols <- setdiff(names(m), expected)
  if (length(ann_cols) == 0L) {
    message("[annotation] no annotation columns in ", path)
    return(NULL)
  }
  m
}

amet_setup_render_logging <- function(log_path, threads = 1L) {
  threads <- amet_parse_threads(threads, 1L)
  if (requireNamespace("data.table", quietly = TRUE)) {
    data.table::setDTthreads(threads)
  }
  if (!is.null(log_path) && !is.na(log_path) && nzchar(log_path)) {
    if (isTRUE(.amet_log_state$sink_active)) {
      message("[render_logging] sinks already active; skipping re-sink")
    } else {
      log_con <- file(log_path, open = "at")
      ## split = TRUE so knitr's stdout capture still sees chunk output
      ## (results = "asis" and cat() emit Markdown into the rendered HTML).
      sink(log_con, split = TRUE)
      sink(log_con, type = "message")
      .amet_log_state$sink_active <- TRUE
      .amet_log_state$log_con <- log_con
    }
  }
  starts <- new.env(parent = emptyenv())
  knitr::knit_hooks$set(progress = function(before, options, envir) {
    label <- if (is.null(options$label) || !nzchar(options$label)) "<unnamed>" else options$label
    if (before) {
      starts[[label]] <- Sys.time()
      message("[chunk start] ", label)
    } else {
      t0 <- starts[[label]]
      if (is.null(t0)) {
        message("[chunk end]   ", label, " (unknown duration)")
      } else {
        elapsed <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
        message("[chunk end]   ", label, " (", round(elapsed, 2), "s)")
      }
    }
  })
  knitr::knit_hooks$set(error = function(x, options) {
    label <- if (is.null(options$label) || !nzchar(options$label)) "<unnamed>" else options$label
    formatted <- tryCatch(paste(format(x), collapse = "\n"),
                          error = function(e) as.character(x))
    message("knitr error in chunk '", label, "':\n", formatted)
    x
  })
  knitr::opts_chunk$set(progress = TRUE)
  invisible(amet_make_bpparam(threads))
}
