## Multi-panel benchmark figure for amet's scores (I_norm, JSD) and the
## wcVI/acVI simulators. Panels:
##   A : example cells at increasing wcVI (within-cell variability)
##   B : example cells at increasing acVI (across-cell variability)
##   C : I_norm per cell vs mean methylation, coloured by wcVI
##   D : JSD per group vs mean methylation, coloured by acVI
##   E : orthogonality: JSD response to wcVI (should be flat for true decoupling)
##   F : orthogonality: I_norm response to acVI (should be flat for true decoupling)
##   G : NMI of I_norm recovery of wcVI vs cells per group (subsampled)
##   H : NMI of JSD recovery of acVI vs cells per group (subsampled)

suppressPackageStartupMessages({
    library(optparse); library(ggplot2); library(dplyr); library(patchwork); library(tidyr)
})
.this_dir <- local({ args <- commandArgs(trailingOnly = FALSE); fa <- grep("^--file=", args, value = TRUE); if (length(fa) > 0) dirname(sub("^--file=", "", fa[1])) else "." })
source(file.path(.this_dir, "plot_theme.R"))
source(file.path(.this_dir, "bricks.R"))

shannon <- function(probs) { probs <- probs[probs > 0]; -sum(probs * log2(probs)) }
shannon_binary <- function(p) {
    out <- numeric(length(p)); safe <- !is.na(p) & p > 0 & p < 1
    out[safe] <- -p[safe] * log2(p[safe]) - (1 - p[safe]) * log2(1 - p[safe])
    out[!safe] <- NA_real_; out
}
nmi <- function(x, y) {
    tab <- table(x, y); n <- sum(tab); p_xy <- tab / n
    p_x <- rowSums(p_xy); p_y <- colSums(p_xy)
    h_x <- shannon(p_x); h_y <- shannon(p_y); h_xy <- shannon(as.vector(p_xy))
    mi <- h_x + h_y - h_xy
    if (h_x + h_y == 0) 0 else 2 * mi / (h_x + h_y)
}
score_to_bins <- function(s, n_bins = 10) {
    keep <- is.finite(s)
    breaks <- unique(quantile(s[keep], probs = seq(0, 1, length.out = n_bins + 1), na.rm = TRUE))
    if (length(breaks) < 2) return(rep(1L, length(s)))
    as.integer(cut(s, breaks = breaks, include.lowest = TRUE))
}
multi_jsd_pmf <- function(pmf_list) {
    if (length(pmf_list) < 2) return(0)
    h_avg <- mean(sapply(pmf_list, function(p) { p <- p[p > 0]; -sum(p * log2(p)) }))
    mix <- Reduce("+", pmf_list) / length(pmf_list)
    h_mix <- { m <- mix[mix > 0]; -sum(m * log2(m)) }
    max(0, h_mix - h_avg)
}
cell_conditional <- function(c00, c01, c10, c11, min_row_count = 20) {
    if ((c00 + c01) < min_row_count || (c10 + c11) < min_row_count) return(NULL)
    sm <- c(c00, c01, c10, c11) + 1
    r0 <- sm[1] + sm[2]; r1 <- sm[3] + sm[4]
    c(sm[1] / r0, sm[2] / r0, sm[3] / r1, sm[4] / r1) / 2
}

options <- list(
    make_option(c("--wcvi_cell_feature"), type = "character"),
    make_option(c("--wcvi_feature"), type = "character"),
    make_option(c("--wcvi_pair_counts"), type = "character"),
    make_option(c("--acvi_cell_feature"), type = "character"),
    make_option(c("--acvi_feature"), type = "character"),
    make_option(c("--acvi_pair_counts"), type = "character"),
    make_option(c("--output_prefix"), type = "character")
)
opt <- parse_args(OptionParser(option_list = options))

read_amet <- function(cf_path, ft_path, pc_path) {
    cf <- read.table(gzfile(cf_path), header = TRUE, sep = "\t", na.strings = "NA",
                     stringsAsFactors = FALSE)
    ft <- read.table(gzfile(ft_path), header = TRUE, sep = "\t", na.strings = "NA",
                     stringsAsFactors = FALSE)
    pc <- read.table(gzfile(pc_path), header = TRUE, sep = "\t", stringsAsFactors = FALSE)
    pc <- pc[pc$lag == 1, ]
    pc_split <- split(pc, list(pc$feature_id, pc$group), drop = TRUE)
    ## pair_counts available but JSD_cond is not used in this figure; JSD is the
    ## within-marginal (no-conditional) version reported by amet.
    list(cf = cf, ft = ft)
}

wc <- read_amet(opt$wcvi_cell_feature, opt$wcvi_feature, opt$wcvi_pair_counts)
ac <- read_amet(opt$acvi_cell_feature, opt$acvi_feature, opt$acvi_pair_counts)

wc$cf$wcvi <- as.integer(sub("wcvi([0-9]+)_.*", "\\1", wc$cf$cell_id))
wc$cf$target_p <- as.numeric(sub(".*_p([0-9.]+)_.*", "\\1", wc$cf$cell_id))
wc$ft$wcvi <- as.integer(sub("wcvi([0-9]+)_.*", "\\1", wc$ft$group))
wc$ft$target_p <- as.numeric(sub(".*_p([0-9.]+)$", "\\1", wc$ft$group))
i_cols <- grep("^i_[0-9]+$", names(wc$cf), value = TRUE)
k_max <- length(i_cols)
wc$cf$h_marg <- shannon_binary(wc$cf$mean_meth)
wc$cf$i_norm <- wc$cf$i_total / (k_max * wc$cf$h_marg)

ac$cf$acvi <- as.integer(sub("acvi([0-9]+)_.*", "\\1", ac$cf$cell_id))
ac$cf$target_p <- as.numeric(sub(".*_p([0-9.]+)_.*", "\\1", ac$cf$cell_id))
ac$ft$acvi <- as.integer(sub("acvi([0-9]+)_.*", "\\1", ac$ft$group))
ac$ft$target_p <- as.numeric(sub(".*_p([0-9.]+)$", "\\1", ac$ft$group))
ac$cf$h_marg <- shannon_binary(ac$cf$mean_meth)
ac$cf$i_norm <- ac$cf$i_total / (k_max * ac$cf$h_marg)

## ===== Panels A and B: simulator visualisations =====
make_wcvi_examples <- function(p = 0.5, n_cpgs = 200, levels = c(1, 3, 5, 7, 10)) {
    set.seed(2026)
    wcvi_to_p_repeat <- function(v) (10L - v) * 0.95 / 9.0
    do.call(rbind, lapply(levels, function(v) {
        bits <- simulate_repeat_cell(n_cpgs, p, wcvi_to_p_repeat(v))
        data.frame(wcvi = v, pos = seq_along(bits), bit = bits)
    }))
}
## Five cells per acVI level so across-cell variability is visible.
make_acvi_examples <- function(p = 0.5, n_cpgs = 200, levels = c(1, 3, 5, 7, 10),
                                 n_per_level = 5) {
    set.seed(2026)
    p_centre <- 0.725; p_half_max <- 0.225
    do.call(rbind, lapply(levels, function(v) {
        half_width <- (v - 1L) * p_half_max / 9.0
        do.call(rbind, lapply(seq_len(n_per_level), function(i) {
            p_repeat_cell <- p_centre + runif(1, -half_width, half_width)
            p_repeat_cell <- max(0, min(0.95, p_repeat_cell))
            bits <- simulate_repeat_cell(n_cpgs, p, p_repeat_cell)
            data.frame(acvi = v, cell = i, pos = seq_along(bits), bit = bits)
        }))
    }))
}
wcvi_ex <- make_wcvi_examples()
acvi_ex <- make_acvi_examples()

## Lollipop-style: filled circle = methylated, open circle = unmethylated.
lollipop_fill <- scale_fill_manual(values = c("0" = "white", "1" = "black"),
                                   guide = "none")
pA <- ggplot(wcvi_ex, aes(x = pos, y = factor(wcvi), fill = factor(bit))) +
    geom_point(shape = 21, size = 1.2, stroke = 0.25, colour = "black") +
    lollipop_fill +
    labs(x = "CpG position", y = "wcVI", title = "A. wcVI cell examples (p = 0.5)") +
    theme_ng() + theme(aspect.ratio = NULL)
pB <- ggplot(acvi_ex, aes(x = pos, y = factor(cell), fill = factor(bit))) +
    geom_point(shape = 21, size = 1.2, stroke = 0.25, colour = "black") +
    lollipop_fill +
    facet_wrap(~ acvi, ncol = 1, scales = "free_y", labeller = label_both) +
    labs(x = "CpG position", y = "cell",
         title = "B. acVI cell examples (p = 0.5, 5 cells per level)") +
    theme_ng() + theme(aspect.ratio = NULL,
                       axis.text.y = element_blank(),
                       axis.ticks.y = element_blank())

## 4-corner sanity panel: cells at the extremes of (wcVI, acVI).
make_corner_examples <- function(p = 0.5, n_cpgs = 200, n_per_corner = 6) {
    set.seed(2027)
    wcvi_to_p_repeat <- function(v) (10L - v) * 0.95 / 9.0
    p_centre <- 0.725; p_half_max <- 0.225
    corners <- list(
        list(label = "wcVI=1, acVI=1",  wcvi = 1,  acvi = 1),
        list(label = "wcVI=1, acVI=10", wcvi = 1,  acvi = 10),
        list(label = "wcVI=10, acVI=1",  wcvi = 10, acvi = 1),
        list(label = "wcVI=10, acVI=10", wcvi = 10, acvi = 10)
    )
    do.call(rbind, lapply(seq_along(corners), function(ci) {
        co <- corners[[ci]]
        do.call(rbind, lapply(seq_len(n_per_corner), function(i) {
            base_p_repeat <- wcvi_to_p_repeat(co$wcvi)
            half_width <- (co$acvi - 1L) * p_half_max / 9.0
            p_repeat_cell <- max(0, min(0.95, base_p_repeat + runif(1, -half_width, half_width)))
            bits <- simulate_repeat_cell(n_cpgs, p, p_repeat_cell)
            data.frame(corner = co$label, corner_idx = ci, cell = i,
                       pos = seq_along(bits), bit = bits)
        }))
    }))
}
corner_ex <- make_corner_examples()
corner_ex$corner <- factor(corner_ex$corner,
                           levels = c("wcVI=1, acVI=1", "wcVI=1, acVI=10",
                                      "wcVI=10, acVI=1", "wcVI=10, acVI=10"))

pCorner <- ggplot(corner_ex, aes(x = pos, y = factor(cell), fill = factor(bit))) +
    geom_point(shape = 21, size = 1.2, stroke = 0.25, colour = "black") +
    lollipop_fill +
    facet_wrap(~ corner, ncol = 2) +
    labs(x = "CpG position", y = "cell",
         title = "Sanity check: cells at the four (wcVI, acVI) corners (p = 0.5)") +
    theme_ng() + theme(aspect.ratio = NULL,
                       axis.text.y = element_blank(),
                       axis.ticks.y = element_blank())

## ===== Panels C and D: scores vs methylation =====
hot <- scale_colour_viridis_c(option = "inferno", limits = c(1, 10),
                              breaks = c(1, 5, 10))

pC <- ggplot(wc$cf, aes(x = mean_meth, y = i_norm, colour = wcvi)) +
    geom_point(alpha = 0.4, size = 0.4) +
    scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25)) +
    hot +
    labs(x = "mean methylation", y = expression(I[norm]), colour = "wcVI",
         title = "C. I_norm per cell vs methylation") +
    theme_ng() + theme(aspect.ratio = NULL)

pD <- ggplot(ac$ft, aes(x = mean_meth_mean, y = jsd, colour = acvi)) +
    geom_point(size = 1) +
    scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25)) +
    hot +
    labs(x = "group mean methylation", y = "JSD", colour = "acVI",
         title = "D. JSD per group vs methylation") +
    theme_ng() + theme(aspect.ratio = NULL)

## ===== Panels E and F: orthogonality =====
pE <- ggplot(wc$ft, aes(x = factor(wcvi), y = jsd, colour = mean_meth_mean)) +
    geom_jitter(width = 0.15, alpha = 0.7, size = 1) +
    scale_colour_viridis_c(option = "inferno", limits = c(0, 1)) +
    stat_summary(fun = mean, geom = "line", aes(group = 1), colour = "black", linewidth = 0.4) +
    stat_summary(fun = mean, geom = "point", colour = "black", size = 1.5) +
    labs(x = "wcVI", y = "JSD", colour = "group mean methylation",
         title = "E. JSD per group vs wcVI (orthogonality)") + theme_ng()
pF <- ggplot(ac$cf[is.finite(ac$cf$i_norm), ],
             aes(x = factor(acvi), y = i_norm, colour = target_p)) +
    geom_jitter(width = 0.25, alpha = 0.4, size = 0.4) +
    scale_colour_viridis_c(option = "inferno", limits = c(0, 1)) +
    stat_summary(fun = mean, geom = "line", aes(group = 1), colour = "black", linewidth = 0.4) +
    stat_summary(fun = mean, geom = "point", colour = "black", size = 1.5) +
    labs(x = "acVI", y = expression(I[norm]), colour = "target p",
         title = "F. I_norm per cell vs acVI (orthogonality)") + theme_ng()

## ===== Panels G and H: NMI vs cells per group (subsampled) =====
nmi_at_subsample <- function(df, true_col, score_col, n_target, n_reps = 5, seed = 7) {
    set.seed(seed)
    dat <- df[is.finite(df[[score_col]]), ]
    truth <- dat[[true_col]]
    score <- dat[[score_col]]
    res <- numeric(n_reps)
    for (r in seq_len(n_reps)) {
        idx <- unlist(lapply(split(seq_len(nrow(dat)), truth), function(g) {
            sample(g, min(length(g), n_target))
        }))
        sub_truth <- truth[idx]
        sub_score <- score[idx]
        res[r] <- nmi(sub_truth, score_to_bins(sub_score))
    }
    mean(res)
}

n_grid <- c(3, 5, 10, 20, 30)

g_data <- expand.grid(score = "I_norm", n_cells = n_grid, stringsAsFactors = FALSE)
g_data$NMI <- mapply(function(s, n) nmi_at_subsample(wc$cf, "wcvi", "i_norm", n_target = n),
                     g_data$score, g_data$n_cells)

h_data <- expand.grid(score = "JSD", n_cells = n_grid, stringsAsFactors = FALSE)
h_data$NMI <- mapply(function(s, n) nmi_at_subsample(ac$ft, "acvi", "jsd", n_target = n),
                     h_data$score, h_data$n_cells)

pG <- ggplot(g_data, aes(x = n_cells, y = NMI, colour = score, group = score)) +
    geom_line() + geom_point(size = 1.5) +
    scale_x_log10(breaks = n_grid) +
    scale_colour_manual(values = c(I_norm = "#d95f02")) +
    coord_cartesian(ylim = c(0, 1)) +
    labs(x = "cells per group (log)", y = "NMI", colour = NULL,
         title = "G. wcVI recovery NMI vs n_cells") + theme_ng()
pH <- ggplot(h_data, aes(x = n_cells, y = NMI, colour = score, group = score)) +
    geom_line() + geom_point(size = 1.5) +
    scale_x_log10(breaks = n_grid) +
    scale_colour_manual(values = c(JSD = "#7570b3")) +
    coord_cartesian(ylim = c(0, 1)) +
    labs(x = "cells per group (log)", y = "NMI", colour = NULL,
         title = "H. acVI recovery NMI vs n_cells") + theme_ng()

## ===== Combine =====
combined <- (pA | pB) / pCorner / pC / pD / (pE | pF) / (pG | pH)
metrics <- rbind(
    cbind(panel = "G", g_data),
    cbind(panel = "H", h_data)
)
save_eval(combined, metrics, opt$output_prefix, width_mm = 200, height_mm = 380)
message("[eval_benchmark_summary] done")
