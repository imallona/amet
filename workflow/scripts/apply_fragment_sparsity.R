## Two sparsity models for synthetic single-cell methylation data.
##
## bernoulli_dropout: each CpG is independently kept with probability p_keep.
##   The trivial baseline; ignores any spatial structure of sequencing coverage.
##
## fragment_coverage: simulate sequencing as random fragments laid down along the
##   genome and call a CpG observed if any fragment covers its position.
##   Controlled by:
##     n_fragments: how many fragments are sequenced from this cell.
##     fragment_length_mean / fragment_length_sd: log-normal distribution for fragment
##         lengths (typical short-read methylation libraries: mean ~300-500 bp).
##     genome_size_bp: total span of the genomic region the cell is read from
##         (used to place fragment starts uniformly).
##
##   This models coverage as a physical sequencing process rather than as a statistical
##   run-length pattern. CpGs near the centre of a long fragment are more likely covered
##   than CpGs at fragment boundaries; the resulting sparsity has spatial correlation
##   that emerges from fragment overlap, not from a hand-tuned run-length distribution.

bernoulli_dropout <- function(bits, p_keep, seed = NA) {
    if (!is.na(seed)) set.seed(seed)
    keep <- runif(length(bits)) < p_keep
    list(bits = bits[keep], indices = which(keep))
}

fragment_coverage <- function(positions_0based,
                              n_fragments,
                              fragment_length_mean = 350,
                              fragment_length_sd = 150,
                              genome_size_bp = NULL,
                              seed = NA) {
    if (!is.na(seed)) set.seed(seed)
    n_cpg <- length(positions_0based)
    if (is.null(genome_size_bp)) {
        genome_size_bp <- max(positions_0based) - min(positions_0based) + 1L
    }
    genome_start <- min(positions_0based)

    log_mean <- log(fragment_length_mean^2 /
                    sqrt(fragment_length_sd^2 + fragment_length_mean^2))
    log_sd <- sqrt(log(1 + (fragment_length_sd / fragment_length_mean)^2))
    fragment_lengths <- pmax(1L, round(rlnorm(n_fragments, log_mean, log_sd)))
    fragment_starts <- floor(runif(n_fragments,
                                   min = genome_start,
                                   max = genome_start + genome_size_bp))
    fragment_ends <- fragment_starts + fragment_lengths - 1L

    covered <- logical(n_cpg)
    if (n_fragments > 0) {
        ord <- order(fragment_starts)
        fs <- fragment_starts[ord]
        fe <- fragment_ends[ord]
        ## Walk fragments in start order; for each, mark CpGs in [fs, fe].
        for (j in seq_along(fs)) {
            lo <- which(positions_0based >= fs[j])
            if (length(lo) == 0L) next
            lo <- lo[1]
            hi <- which(positions_0based[lo:n_cpg] <= fe[j])
            if (length(hi) == 0L) next
            hi <- lo + max(hi) - 1L
            covered[lo:hi] <- TRUE
        }
    }
    list(indices = which(covered),
         n_fragments = n_fragments,
         fragment_length_mean_realised = mean(fragment_lengths))
}
