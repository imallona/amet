## acVI recovery: per-group JSD against the index, per-cell I_norm
## as the orthogonality check.
##
## Layout (4 panels):
##   A   JSD vs acVI                        (per group, colour=mean_meth_mean)
##   B   JSD_norm vs mean_meth              (per group, faceted by acVI)
##   C   I_norm vs acVI                     (per cell, colour=target_p)
##   D   I_norm vs mean_meth                (per cell, faceted by acVI)
##
## All methylation axes use the full (0, 1) range.

suppressPackageStartupMessages({
    library(optparse); library(ggplot2); library(dplyr); library(patchwork)
})
.this_dir <- local({ args <- commandArgs(trailingOnly = FALSE); fa <- grep("^--file=", args, value = TRUE); if (length(fa) > 0) dirname(sub("^--file=", "", fa[1])) else "." })
source(file.path(.this_dir, "plot_theme.R"))

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
recovery <- function(true_int, score) {
    keep <- is.finite(score) & is.finite(true_int)
    if (sum(keep) < 5) return(data.frame(NMI = NA, Spearman = NA, Kendall = NA))
    s <- score[keep]; t <- true_int[keep]
    bin <- as.integer(cut(s,
        breaks = unique(quantile(s, probs = seq(0, 1, length.out = 11), na.rm = TRUE)),
        include.lowest = TRUE))
    data.frame(NMI = nmi(t, bin),
               Spearman = cor(t, s, method = "spearman"),
               Kendall = cor(t, s, method = "kendall"))
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
    make_option(c("--cell_feature"), type = "character"),
    make_option(c("--feature"), type = "character"),
    make_option(c("--pair_counts"), type = "character"),
    make_option(c("--output_prefix"), type = "character")
)
opt <- parse_args(OptionParser(option_list = options))

cf <- read.table(gzfile(opt$cell_feature), header = TRUE, sep = "\t",
                 na.strings = "NA", stringsAsFactors = FALSE)
cf$acvi <- as.integer(sub("acvi([0-9]+)_.*", "\\1", cf$cell_id))
cf$target_p <- as.numeric(sub(".*_p([0-9.]+)_.*", "\\1", cf$cell_id))
cf <- cf[!is.na(cf$i_total) & !is.na(cf$acvi), ]
i_cols <- grep("^i_[0-9]+$", names(cf), value = TRUE)
k_max <- length(i_cols)
cf$h_marg <- shannon_binary(cf$mean_meth)
cf$i_norm <- cf$i_total / (k_max * cf$h_marg)

ft <- read.table(gzfile(opt$feature), header = TRUE, sep = "\t",
                 na.strings = "NA", stringsAsFactors = FALSE)
ft$acvi <- as.integer(sub("acvi([0-9]+)_.*", "\\1", ft$group))
ft$target_p <- as.numeric(sub(".*_p([0-9.]+)$", "\\1", ft$group))
ft <- ft[!is.na(ft$acvi), ]

pc <- read.table(gzfile(opt$pair_counts), header = TRUE, sep = "\t",
                 stringsAsFactors = FALSE)
pc <- pc[pc$lag == 1, ]
pc_split <- split(pc, list(pc$feature_id, pc$group), drop = TRUE)
jsd_cond_per_group <- do.call(rbind, lapply(pc_split, function(sub) {
    pmfs <- lapply(seq_len(nrow(sub)),
                   function(i) cell_conditional(sub$n00[i], sub$n01[i], sub$n10[i], sub$n11[i]))
    pmfs <- pmfs[!sapply(pmfs, is.null)]
    data.frame(feature_id = sub$feature_id[1], group = sub$group[1],
               jsd_cond = multi_jsd_pmf(pmfs),
               n_cells_used = length(pmfs))
}))
ft <- merge(ft, jsd_cond_per_group, by = c("feature_id", "group"), all.x = TRUE)
ft$jsd_norm <- ft$jsd / (2 * shannon_binary(ft$mean_meth_mean))
ft$jsd_norm[!is.finite(ft$jsd_norm)] <- NA

eval_score <- function(score_name, true_int, score, scope_label) {
    cbind(scope = scope_label, score = score_name, recovery(true_int, score))
}
metrics <- rbind(
    eval_score("JSD_raw",  ft$acvi, ft$jsd,      "overall (per group)"),
    eval_score("JSD_norm", ft$acvi, ft$jsd_norm, "overall (per group)"),
    eval_score("I_norm",   cf$acvi, cf$i_norm,   "overall (per cell)")
)
rownames(metrics) <- NULL

annot <- paste(sprintf("%-10s %s: NMI=%.2f  Spearman=%.2f",
                       metrics$score, metrics$scope, metrics$NMI, metrics$Spearman),
               collapse = "\n")

x_meth <- scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25))
meth_grad <- scale_colour_viridis_c(option = "inferno", limits = c(0, 1))

pA <- ggplot(ft, aes(x = factor(acvi), y = jsd)) +
    geom_jitter(aes(colour = mean_meth_mean), width = 0.15, alpha = 0.7, size = 1) +
    meth_grad +
    stat_summary(fun = mean, geom = "line", aes(group = 1), colour = "black", linewidth = 0.4) +
    stat_summary(fun = mean, geom = "point", size = 1.5, colour = "black") +
    scale_y_sqrt() +
    labs(x = "acVI", y = "JSD (sqrt)", colour = "group mean methylation",
         title = "A. JSD per group vs acVI") + theme_ng() +
    theme(aspect.ratio = NULL)
pB <- ggplot(ft[is.finite(ft$jsd_norm), ], aes(x = mean_meth_mean, y = jsd_norm)) +
    geom_point(size = 1, alpha = 0.6, colour = "grey40") +
    geom_smooth(method = "loess", se = FALSE, colour = "darkred", linewidth = 0.5, span = 0.8) +
    facet_wrap(~ acvi, ncol = 5, labeller = label_both) +
    x_meth +
    labs(x = "group mean methylation",
         y = expression(JSD[norm] ~ "= JSD / (2 H(p))"),
         title = "B. JSD_norm per group vs methylation, faceted by acVI (inverted-U flattened)") +
    theme_ng() + theme(aspect.ratio = NULL)
pC <- ggplot(cf[is.finite(cf$i_norm), ],
             aes(x = factor(acvi), y = i_norm)) +
    geom_jitter(aes(colour = target_p), width = 0.25, alpha = 0.4, size = 0.4) +
    meth_grad +
    stat_summary(fun = mean, geom = "line", aes(group = 1), colour = "black", linewidth = 0.4) +
    stat_summary(fun = mean, geom = "point", size = 1.5, colour = "black") +
    labs(x = "acVI", y = expression(I[norm]), colour = "target p",
         title = "C. I_norm per cell vs acVI (orthogonality)") + theme_ng() +
    theme(aspect.ratio = NULL)
pD <- ggplot(cf[is.finite(cf$i_norm), ], aes(x = mean_meth, y = i_norm)) +
    geom_point(alpha = 0.35, size = 0.4, colour = "grey40") +
    geom_smooth(method = "loess", se = FALSE, colour = "darkred", linewidth = 0.5, span = 0.6) +
    facet_wrap(~ acvi, ncol = 5, labeller = label_both) +
    x_meth +
    labs(x = "mean methylation", y = expression(I[norm]),
         title = "D. I_norm per cell vs methylation, faceted by acVI") + theme_ng() +
    theme(aspect.ratio = NULL)

combined <- pA / pB / pC / pD +
    plot_annotation(subtitle = annot, theme = theme(plot.subtitle = element_text(size = 6)))
save_eval(combined, metrics, opt$output_prefix, width_mm = 220, height_mm = 280)
message("[eval_acvi_recovery] done")
