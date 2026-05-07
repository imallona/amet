## Generate a single-cell binary methylation sequence from a 2-state Markov chain.
##
## Transition matrix M[s, s'] = P(next = s' | current = s).
## Parameters: p11 (P(1|1)) and p00 (P(0|0)); off-diagonals follow.
##
## Stationary marginal:
##   pi_1 = (1 - p00) / (2 - p00 - p11)
##
## Persistence (comethylation): p11 + p00 > 1 is persistent (ordered);
## p11 + p00 < 1 is anti-persistent; p11 + p00 = 1 is i.i.d.

suppressPackageStartupMessages({
    library(optparse)
})

.this_dir <- local({ args <- commandArgs(trailingOnly = FALSE); fa <- grep("^--file=", args, value = TRUE); if (length(fa) > 0) dirname(sub("^--file=", "", fa[1])) else "." })

source(file.path(.this_dir, "write_outputs.R"))

simulate_markov_cell <- function(n, p00, p11, seed = NA) {
    if (!is.na(seed)) set.seed(seed)
    if (n < 1) return(integer(0))
    pi1 <- (1 - p00) / (2 - p00 - p11)
    state <- if (runif(1) < pi1) 1L else 0L
    out <- integer(n)
    out[1] <- state
    if (n == 1) return(out)
    for (i in seq.int(2, n)) {
        if (state == 0L) {
            state <- if (runif(1) < (1 - p00)) 1L else 0L
        } else {
            state <- if (runif(1) < p11) 1L else 0L
        }
        out[i] <- state
    }
    out
}

stationary_marginal <- function(p00, p11) {
    (1 - p00) / (2 - p00 - p11)
}

## True per-step entropy rate of the chain (bits per CpG).
markov_entropy_rate <- function(p00, p11) {
    pi1 <- stationary_marginal(p00, p11)
    pi0 <- 1 - pi1
    safe_log2 <- function(x) ifelse(x > 0, log2(x), 0)
    h_row0 <- -p00 * safe_log2(p00) - (1 - p00) * safe_log2(1 - p00)
    h_row1 <- -p11 * safe_log2(p11) - (1 - p11) * safe_log2(1 - p11)
    pi0 * h_row0 + pi1 * h_row1
}

## True I_1 of the stationary chain: I_1 = H(pi) - h_rate.
markov_i1 <- function(p00, p11) {
    pi1 <- stationary_marginal(p00, p11)
    pi0 <- 1 - pi1
    h_marg <- if (pi1 == 0 || pi1 == 1) 0 else -pi1 * log2(pi1) - pi0 * log2(pi0)
    h_marg - markov_entropy_rate(p00, p11)
}

if (sys.nframe() == 0) {
    options <- list(
        make_option(c("--n_cpgs"), type = "integer", default = 100),
        make_option(c("--p00"), type = "double", default = 0.8),
        make_option(c("--p11"), type = "double", default = 0.8),
        make_option(c("--chrom"), type = "character", default = "chr1"),
        make_option(c("--start_pos"), type = "integer", default = 100),
        make_option(c("--cpg_step"), type = "integer", default = 100),
        make_option(c("--seed"), type = "integer", default = 1),
        make_option(c("--output"), type = "character", default = "cell.allc.tsv.gz")
    )
    opt <- parse_args(OptionParser(option_list = options))
    bits <- simulate_markov_cell(opt$n_cpgs, opt$p00, opt$p11, opt$seed)
    positions_0based <- opt$start_pos + (seq_along(bits) - 1) * opt$cpg_step
    write_allc(opt$output, opt$chrom, positions_0based, bits)
    message(sprintf(
        "[markov] wrote %d CpGs to %s (true pi1=%.3f, true I_1=%.4f bits)",
        length(bits), opt$output,
        stationary_marginal(opt$p00, opt$p11),
        markov_i1(opt$p00, opt$p11)
    ))
}
