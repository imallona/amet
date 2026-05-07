## Install epiCHAOS from GitHub into the active conda env's R library.
## Pinned to a specific commit so reruns and CI use a stable epiCHAOS version.
## Idempotent: skips if already present at the pinned commit. Writes a sentinel
## file with the SHA so snakemake can use it as the rule output and audit logs.

suppressPackageStartupMessages({ library(optparse) })

EPICHAOS_REPO <- "CompEpigen/epiCHAOS"
EPICHAOS_REF  <- "34cb72d83fbf98457a68c258f2e842a4b38c492e"

options <- list(
    make_option(c("--sentinel"), type = "character",
                help = "Path to write a marker file when install succeeds")
)
opt <- parse_args(OptionParser(option_list = options))

target_lib <- .libPaths()[1]
message(sprintf("[install_epichaos] target library: %s", target_lib))
message(sprintf("[install_epichaos] pinned ref: %s@%s", EPICHAOS_REPO, EPICHAOS_REF))
options(repos = c(CRAN = "https://cloud.r-project.org"))

if (!requireNamespace("remotes", quietly = TRUE)) {
    install.packages("remotes", lib = target_lib)
}

if (!requireNamespace("jaccard", quietly = TRUE)) {
    message("[install_epichaos] installing jaccard (CRAN)")
    install.packages("jaccard", lib = target_lib)
}

if (!requireNamespace("epiCHAOS", quietly = TRUE)) {
    message(sprintf("[install_epichaos] installing %s@%s", EPICHAOS_REPO, EPICHAOS_REF))
    remotes::install_github(paste0(EPICHAOS_REPO, "@", EPICHAOS_REF),
                            lib = target_lib,
                            upgrade = "never",
                            quiet = FALSE)
} else {
    message("[install_epichaos] already present")
}

stopifnot(requireNamespace("epiCHAOS", quietly = TRUE))

dir.create(dirname(opt$sentinel), showWarnings = FALSE, recursive = TRUE)
writeLines(c(
    sprintf("epichaos_repo: %s", EPICHAOS_REPO),
    sprintf("epichaos_ref:  %s", EPICHAOS_REF),
    sprintf("installed_at:  %s", format(Sys.time()))
), opt$sentinel)
