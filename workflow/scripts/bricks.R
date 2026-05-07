## Templated brick library for the acVI sweep, plus a Markov-with-repeat process
## for the wcVI sweep. The two share no machinery so wcVI and acVI are decoupled
## by construction:
##   - wcVI cells come from `simulate_repeat_cell` (independent realisations of a
##     Markov-with-repeat process at marginal target_p, regularity = p_repeat).
##     Each cell is independent -> across-cell variability is finite-sample noise,
##     not a function of wcVI.
##   - acVI cells use templated bricks (one brick per cell, repeated). Templates
##     have very different L-mer signatures so JSD has dynamic range across the
##     acVI ladder.

## ---------- Markov-with-repeat process for wcVI ----------

## Each CpG either copies the previous one (prob p_repeat) or is sampled iid
## Bernoulli at target_p. Stationary marginal = target_p exactly.
##   p_repeat = 0   -> pure iid Bernoulli at target_p (no within-cell structure).
##   p_repeat = 0.95 -> long runs of the same state.
simulate_repeat_cell <- function(n_cpgs, target_p, p_repeat) {
    bits <- integer(n_cpgs)
    bits[1] <- as.integer(runif(1) < target_p)
    if (n_cpgs == 1) return(bits)
    iid_draws <- as.integer(runif(n_cpgs - 1) < target_p)
    repeats <- runif(n_cpgs - 1) < p_repeat
    for (i in 2:n_cpgs) {
        bits[i] <- if (repeats[i - 1]) bits[i - 1] else iid_draws[i - 1]
    }
    bits
}

## ---------- Templated bricks for acVI ----------
## Each template takes (L, k) where k = number of methylated CpGs in the brick
## (k = round(target_p * L)). All templates produce a length-L binary vector with
## exactly k ones, but with very different intrinsic patterns. Templates are
## indexed 1..10 so the acVI pool grows from 1 distinct pattern (acVI=1) to 10
## distinct patterns (acVI=10).

template_block_right <- function(L, k) {
    c(rep(0L, L - k), rep(1L, k))
}
template_block_left <- function(L, k) {
    c(rep(1L, k), rep(0L, L - k))
}
## Evenly spaced positions guaranteed unique: (0..k-1) * L %/% k + 1.
template_alternating_a <- function(L, k) {
    if (k == 0) return(rep(0L, L))
    if (k == L) return(rep(1L, L))
    pos <- (seq.int(0L, k - 1L) * L %/% k) + 1L
    b <- integer(L); b[pos] <- 1L; b
}
template_alternating_b <- function(L, k) {
    if (k == 0) return(rep(0L, L))
    if (k == L) return(rep(1L, L))
    pos <- (seq.int(0L, k - 1L) * L %/% k) + 1L
    pos <- ((pos %% L) + 1L)  # shift by 1
    b <- integer(L); b[pos] <- 1L; b
}
## Two clumps centred around L/4 and 3L/4. Falls back to a single clump if the
## two would overlap; always returns a brick of length L with exactly k ones.
template_two_clumps_a <- function(L, k) {
    if (k == 0) return(rep(0L, L))
    if (k == L) return(rep(1L, L))
    b <- integer(L)
    h1 <- k %/% 2L; h2 <- k - h1
    q <- L %/% 4L
    s1 <- max(1L, q - h1 %/% 2L + 1L)
    e1 <- min(L, s1 + h1 - 1L)
    s2 <- max(e1 + 1L, 3L * q - h2 %/% 2L + 1L)
    e2 <- min(L, s2 + h2 - 1L)
    if (e1 - s1 + 1L < h1) {
        e1 <- min(L, s1 + h1 - 1L)
    }
    if (e2 - s2 + 1L < h2) {
        s2 <- max(1L, L - h2 + 1L)
        e2 <- L
    }
    b[seq.int(s1, e1)] <- 1L
    b[seq.int(s2, e2)] <- 1L
    actual_k <- sum(b)
    if (actual_k < k) {
        zeros <- which(b == 0L)
        b[zeros[seq_len(k - actual_k)]] <- 1L
    } else if (actual_k > k) {
        ones <- which(b == 1L)
        b[ones[seq_len(actual_k - k)]] <- 0L
    }
    b
}
## Two clumps shifted toward the edges of the brick.
template_two_clumps_b <- function(L, k) {
    if (k == 0) return(rep(0L, L))
    if (k == L) return(rep(1L, L))
    b <- integer(L)
    h1 <- k %/% 2L; h2 <- k - h1
    o <- max(1L, L %/% 8L)
    s1 <- o
    e1 <- min(L, s1 + h1 - 1L)
    s2 <- max(e1 + 1L, L - o - h2 + 1L)
    e2 <- min(L, s2 + h2 - 1L)
    b[seq.int(s1, e1)] <- 1L
    if (s2 <= e2) b[seq.int(s2, e2)] <- 1L
    actual_k <- sum(b)
    if (actual_k < k) {
        zeros <- which(b == 0L)
        b[zeros[seq_len(k - actual_k)]] <- 1L
    } else if (actual_k > k) {
        ones <- which(b == 1L)
        b[ones[seq_len(actual_k - k)]] <- 0L
    }
    b
}
template_scattered <- function(L, k, seed) {
    if (k == 0) return(rep(0L, L))
    if (k == L) return(rep(1L, L))
    set.seed(seed * 17L + L * 1000L + k * 7L)
    pos <- sample(L, k)
    b <- integer(L); b[pos] <- 1L; b
}

## Order templates so each new entry adds a genuinely new L-mer signature when
## repeated to fill a cell. Entries 1-4 are the four primary signatures (block,
## alternating, clumped, scattered); entries 5-10 are minor variants that fill in
## the curve between primary jumps. block_left and the alternating shift are NOT
## near the start because they produce the same L-mer histograms as their
## counterparts (mirror / shift symmetries when repeated).
template_funcs <- list(
    function(L, k) template_block_right(L, k),
    function(L, k) template_alternating_a(L, k),
    function(L, k) template_two_clumps_a(L, k),
    function(L, k) template_scattered(L, k, seed = 1L),
    function(L, k) template_scattered(L, k, seed = 2L),
    function(L, k) template_two_clumps_b(L, k),
    function(L, k) template_scattered(L, k, seed = 3L),
    function(L, k) template_alternating_b(L, k),
    function(L, k) template_scattered(L, k, seed = 4L),
    function(L, k) template_block_left(L, k)
)

## Build the acVI brick pool of size `pool_size` at marginal target_p, using the
## first `pool_size` templates from the list.
build_acvi_pool <- function(L, target_p, pool_size) {
    k <- min(L, max(0L, as.integer(round(target_p * L))))
    n <- min(pool_size, length(template_funcs))
    lapply(seq_len(n), function(i) template_funcs[[i]](L, k))
}

## Cell built by repeating one brick from the pool across all slots.
build_acvi_cell <- function(pool, n_slots) {
    chosen <- pool[[sample(length(pool), 1L)]]
    rep(chosen, length.out = n_slots * length(chosen))
}
