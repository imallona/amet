## Generate K cells from a mixture of K Markov chains, each chain with its own
## (p00, p11). Cell i is sampled from component i mod K, so the mixture has K
## equally-weighted components.

source(file.path(dirname(sys.frame(1)$ofile), "write_outputs.R"))
source(file.path(dirname(sys.frame(1)$ofile), "simulate_markov_chain.R"))

simulate_mixture_cells <- function(n_cells, n_cpgs, components, seed = 1) {
    set.seed(seed)
    k <- length(components)
    cells <- vector("list", n_cells)
    component_assignments <- integer(n_cells)
    for (i in seq_len(n_cells)) {
        comp_idx <- ((i - 1L) %% k) + 1L
        comp <- components[[comp_idx]]
        cells[[i]] <- simulate_markov_cell(n_cpgs, comp$p00, comp$p11)
        component_assignments[i] <- comp_idx
    }
    list(cells = cells, components = component_assignments)
}
