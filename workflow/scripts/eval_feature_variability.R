## Recovery of per-feature dispersion gamma_f on the feature_variability_sweep.
##
## Each feature has a true gamma_f; cells in the (single) group draw their
## per-cell-per-feature methylation rate from Beta(mu, gamma_f), then emit iid
## CpG calls. The strongest recovery should come from scMET, which is designed
## around per-feature beta-binomial overdispersion. amet's per-feature JSD
## should also rise with gamma_f. amet's I_norm should be flat (no within-feature
## spatial structure). epiCHAOS's per-feature eITH should pick up the dispersion
## via Jaccard distance among cells at each feature.

suppressPackageStartupMessages({
    library(optparse); library(ggplot2); library(dplyr); library(patchwork); library(tidyr)
})
.this_dir <- local({ args <- commandArgs(trailingOnly = FALSE); fa <- grep("^--file=", args, value = TRUE); if (length(fa) > 0) dirname(sub("^--file=", "", fa[1])) else "." })
source(file.path(.this_dir, "plot_theme.R"))

shannon_binary <- function(p) {
    out <- numeric(length(p)); safe <- !is.na(p) & p > 0 & p < 1
    out[safe] <- -p[safe] * log2(p[safe]) - (1 - p[safe]) * log2(1 - p[safe])
    out[!safe] <- NA_real_; out
}
spearman <- function(true_v, score) {
    keep <- is.finite(score) & is.finite(true_v)
    if (sum(keep) < 5) return(NA_real_)
    cor(true_v[keep], score[keep], method = "spearman")
}

options <- list(
    make_option(c("--cell_feature"), type = "character"),
    make_option(c("--feature"), type = "character"),
    make_option(c("--ground_truth"), type = "character"),
    make_option(c("--scmet"), type = "character"),
    make_option(c("--epichaos"), type = "character"),
    make_option(c("--amet_bench"), type = "character"),
    make_option(c("--scmet_bench"), type = "character"),
    make_option(c("--epichaos_bench"), type = "character"),
    make_option(c("--output_prefix"), type = "character")
)
opt <- parse_args(OptionParser(option_list = options))

cf <- read.table(gzfile(opt$cell_feature), header = TRUE, sep = "\t",
                 na.strings = "NA", stringsAsFactors = FALSE)
ft <- read.table(gzfile(opt$feature), header = TRUE, sep = "\t",
                 na.strings = "NA", stringsAsFactors = FALSE)
gt <- read.table(opt$ground_truth, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
names(gt)[names(gt) == "gamma"] <- "gamma_truth"
names(gt)[names(gt) == "mu"] <- "mu_truth"

i_cols <- grep("^i_[0-9]+$", names(cf), value = TRUE); k_max <- length(i_cols)
cf$i_norm <- cf$i_total / (k_max * shannon_binary(cf$mean_meth))
per_feature_amet <- cf %>%
    group_by(feature_id) %>%
    summarise(i_norm = mean(i_norm, na.rm = TRUE), .groups = "drop")
ft$jsd_norm <- ft$jsd / (2 * shannon_binary(ft$mean_meth_mean))

scmet <- read.table(opt$scmet, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
names(scmet)[names(scmet) == "group"] <- "feature_id"
epi <- read.table(opt$epichaos, header = TRUE, sep = "\t", stringsAsFactors = FALSE)

joined <- gt %>%
    inner_join(per_feature_amet, by = "feature_id") %>%
    inner_join(ft, by = "feature_id", suffix = c("", ".ft")) %>%
    inner_join(scmet, by = "feature_id") %>%
    inner_join(epi, by = "feature_id")

scores <- c("i_norm", "jsd", "jsd_norm", "mu", "gamma", "epsilon", "eITH")
metrics <- do.call(rbind, lapply(scores, function(s)
    data.frame(score = s, Spearman = spearman(joined$gamma_truth, joined[[s]]))))
metrics$score <- factor(metrics$score, levels = scores)

score_pal <- c(i_norm = "#1b9e77", jsd = "#1b9e77", jsd_norm = "#1b9e77",
               mu = "#d95f02", gamma = "#d95f02", epsilon = "#d95f02",
               eITH = "#7570b3")
tool_pal <- c(amet = "#1b9e77", scMET = "#d95f02", epiCHAOS = "#7570b3")

bar_recovery <- ggplot(metrics, aes(x = score, y = abs(Spearman), fill = score)) +
    geom_col() +
    geom_text(aes(label = sprintf("%.2f", Spearman)), vjust = -0.3, size = 2) +
    scale_fill_manual(values = score_pal, guide = "none") +
    coord_cartesian(ylim = c(0, 1.05)) +
    labs(x = NULL, y = "|Spearman rho|",
         title = "A. Recovery of per-feature dispersion gamma") +
    theme_ng() + theme(aspect.ratio = NULL,
                       axis.text.x = element_text(angle = 30, hjust = 1))

long <- joined %>%
    pivot_longer(cols = all_of(scores), names_to = "score", values_to = "value")
long$score <- factor(long$score, levels = scores)
long$family <- ifelse(long$score %in% c("i_norm", "jsd", "jsd_norm"), "amet",
               ifelse(long$score %in% c("mu", "gamma", "epsilon"), "scMET", "epiCHAOS"))

trace <- ggplot(long, aes(x = gamma_truth, y = value, colour = score)) +
    geom_jitter(width = 0.0, alpha = 0.7, size = 0.8) +
    geom_smooth(aes(group = score), method = "loess", se = FALSE,
                linewidth = 0.4, span = 1.0) +
    scale_colour_manual(values = score_pal) +
    facet_wrap(~ family, scales = "free_y", ncol = 3) +
    labs(x = "ground-truth gamma (per feature)",
         y = "score (per feature)",
         title = "B. Score response to gamma") +
    theme_ng() + theme(aspect.ratio = NULL)

read_bench <- function(p, tool) {
    if (!file.exists(p)) return(data.frame(tool = tool, cpu_seconds = NA_real_, mem_mb = NA_real_))
    d <- read.table(p, header = TRUE, sep = "\t")
    data.frame(tool = tool, cpu_seconds = d$cpu_time[1], mem_mb = d$max_rss[1])
}
bench <- rbind(
    read_bench(opt$amet_bench,     "amet"),
    read_bench(opt$scmet_bench,    "scMET"),
    read_bench(opt$epichaos_bench, "epiCHAOS")
)
bar_cpu <- ggplot(bench, aes(x = tool, y = cpu_seconds, fill = tool)) +
    geom_col() + scale_fill_manual(values = tool_pal, guide = "none") +
    labs(x = NULL, y = "CPU seconds", title = "C. CPU time") +
    theme_ng() + theme(aspect.ratio = NULL)
bar_mem <- ggplot(bench, aes(x = tool, y = mem_mb, fill = tool)) +
    geom_col() + scale_fill_manual(values = tool_pal, guide = "none") +
    labs(x = NULL, y = "max RSS (MB)", title = "D. Peak memory") +
    theme_ng() + theme(aspect.ratio = NULL)

combined <- bar_recovery / trace / (bar_cpu | bar_mem)

metrics_csv <- rbind(
    cbind(table = "recovery", score = as.character(metrics$score),
          metric = "Spearman", value = metrics$Spearman),
    cbind(table = "benchmark", score = bench$tool,
          metric = "cpu_seconds", value = bench$cpu_seconds),
    cbind(table = "benchmark", score = bench$tool,
          metric = "mem_mb", value = bench$mem_mb)
)
save_eval(combined, metrics_csv, opt$output_prefix, width_mm = 220, height_mm = 240)
message("[eval_feature_variability] done")
