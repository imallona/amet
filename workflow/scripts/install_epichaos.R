## Install epiCHAOS from GitHub into the active conda env's R library.
## Idempotent: skips if already present. Writes a sentinel file so snakemake
## can use it as the rule output.

suppressPackageStartupMessages({ library(optparse) })

options <- list(
    make_option(c("--sentinel"), type = "character",
                help = "Path to write a marker file when install succeeds")
)
opt <- parse_args(OptionParser(option_list = options))

target_lib <- .libPaths()[1]
message(sprintf("[install_epichaos] target library: %s", target_lib))
options(repos = c(CRAN = "https://cloud.r-project.org"))

if (!requireNamespace("remotes", quietly = TRUE)) {
    install.packages("remotes", lib = target_lib)
}

if (!requireNamespace("jaccard", quietly = TRUE)) {
    message("[install_epichaos] installing jaccard (CRAN)")
    install.packages("jaccard", lib = target_lib)
}

if (!requireNamespace("epiCHAOS", quietly = TRUE)) {
    message("[install_epichaos] installing CompEpigen/epiCHAOS")
    remotes::install_github("CompEpigen/epiCHAOS",
                            lib = target_lib,
                            upgrade = "never",
                            quiet = FALSE)
} else {
    message("[install_epichaos] already present")
}

stopifnot(requireNamespace("epiCHAOS", quietly = TRUE))

dir.create(dirname(opt$sentinel), showWarnings = FALSE, recursive = TRUE)
writeLines(format(Sys.time()), opt$sentinel)
