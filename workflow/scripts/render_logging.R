## Shared render-time logging for analytical and figure Rmds.
## Called from a `logging_early` chunk at the top of each Rmd.

amet_setup_render_logging <- function(log_path) {
  if (nzchar(log_path)) {
    log_con <- file(log_path, open = "at")
    sink(log_con)
    sink(log_con, type = "message")
  }
  starts <- new.env(parent = emptyenv())
  knitr::knit_hooks$set(progress = function(before, options, envir) {
    label <- if (is.null(options$label)) "<unnamed>" else options$label
    if (before) {
      starts[[label]] <- Sys.time()
      message("[chunk start] ", label)
    } else {
      elapsed <- as.numeric(Sys.time() - starts[[label]])
      message("[chunk end]   ", label, " (", round(elapsed, 2), "s)")
    }
  })
  knitr::knit_hooks$set(error = function(x, options) {
    message("knitr error in chunk '", options$label, "':\n", x)
    knitr::knit_exit()
  })
  knitr::opts_chunk$set(progress = TRUE)
}
