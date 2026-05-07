## Head-to-head comparison of amet against scMET and epiCHAOS on the wcVI/acVI
## sweeps. Reports recovery (NMI / Spearman) per score and runtime / memory
## footprint per tool from snakemake benchmark TSVs.

suppressPackageStartupMessages({
    library(optparse); library(ggplot2); library(dplyr); library(patchwork); library(tidyr)
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

options <- list(
    make_option(c("--wcvi_cell_feature"), type = "character"),
    make_option(c("--wcvi_feature"), type = "character"),
    make_option(c("--wcvi_scmet"), type = "character"),
    make_option(c("--wcvi_epichaos"), type = "character"),
    make_option(c("--acvi_cell_feature"), type = "character"),
    make_option(c("--acvi_feature"), type = "character"),
    make_option(c("--acvi_scmet"), type = "character"),
    make_option(c("--acvi_epichaos"), type = "character"),
    make_option(c("--wcvi_amet_bench"), type = "character"),
    make_option(c("--acvi_amet_bench"), type = "character"),
    make_option(c("--wcvi_scmet_bench"), type = "character"),
    make_option(c("--acvi_scmet_bench"), type = "character"),
    make_option(c("--wcvi_epichaos_bench"), type = "character"),
    make_option(c("--acvi_epichaos_bench"), type = "character"),
    make_option(c("--output_prefix"), type = "character")
)
opt <- parse_args(OptionParser(option_list = options))

read_amet_per_group <- function(cf_path, ft_path) {
    cf <- read.table(gzfile(cf_path), header = TRUE, sep = "\t",
                     na.strings = "NA", stringsAsFactors = FALSE)
    ft <- read.table(gzfile(ft_path), header = TRUE, sep = "\t",
                     na.strings = "NA", stringsAsFactors = FALSE)
    i_cols <- grep("^i_[0-9]+$", names(cf), value = TRUE); k_max <- length(i_cols)
    cf$i_norm <- cf$i_total / (k_max * shannon_binary(cf$mean_meth))
    per_group <- cf %>%
        group_by(group) %>%
        summarise(i_total = mean(i_total, na.rm = TRUE),
                  i_norm = mean(i_norm, na.rm = TRUE),
                  mean_meth = mean(mean_meth, na.rm = TRUE),
                  .groups = "drop")
    ft$jsd_norm <- ft$jsd / (2 * shannon_binary(ft$mean_meth_mean))
    list(per_group = per_group, ft = ft)
}

wc <- read_amet_per_group(opt$wcvi_cell_feature, opt$wcvi_feature)
ac <- read_amet_per_group(opt$acvi_cell_feature, opt$acvi_feature)
wc_scmet <- read.table(opt$wcvi_scmet, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
ac_scmet <- read.table(opt$acvi_scmet, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
wc_epi <- read.table(opt$wcvi_epichaos, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
ac_epi <- read.table(opt$acvi_epichaos, header = TRUE, sep = "\t", stringsAsFactors = FALSE)

wc_joined <- wc$per_group %>%
    inner_join(wc$ft, by = "group", suffix = c("", ".ft")) %>%
    inner_join(wc_scmet, by = "group") %>%
    inner_join(wc_epi, by = "group") %>%
    mutate(wcvi = as.integer(sub("wcvi([0-9]+)_.*", "\\1", group)))
ac_joined <- ac$per_group %>%
    inner_join(ac$ft, by = "group", suffix = c("", ".ft")) %>%
    inner_join(ac_scmet, by = "group") %>%
    inner_join(ac_epi, by = "group") %>%
    mutate(acvi = as.integer(sub("acvi([0-9]+)_.*", "\\1", group)))

eval_set <- function(j, true_col, score_cols) {
    do.call(rbind, lapply(score_cols, function(s) {
        r <- recovery(j[[true_col]], j[[s]])
        cbind(score = s, r)
    }))
}

scores_all <- c("i_norm", "jsd_norm",
                "mu", "gamma", "epsilon", "eITH")
wc_metrics <- cbind(axis = "wcVI", eval_set(wc_joined, "wcvi", scores_all))
ac_metrics <- cbind(axis = "acVI", eval_set(ac_joined, "acvi", scores_all))
metrics <- rbind(wc_metrics, ac_metrics)
metrics$score <- factor(metrics$score, levels = scores_all)

read_bench <- function(p, tool, axis) {
    if (!file.exists(p)) return(data.frame(tool = tool, axis = axis,
                                            cpu_seconds = NA_real_, mem_mb = NA_real_))
    d <- read.table(p, header = TRUE, sep = "\t")
    data.frame(tool = tool, axis = axis,
               cpu_seconds = d$cpu_time[1],
               mem_mb = d$max_rss[1])
}
bench <- rbind(
    read_bench(opt$wcvi_amet_bench, "amet", "wcVI"),
    read_bench(opt$acvi_amet_bench, "amet", "acVI"),
    read_bench(opt$wcvi_scmet_bench, "scMET", "wcVI"),
    read_bench(opt$acvi_scmet_bench, "scMET", "acVI"),
    read_bench(opt$wcvi_epichaos_bench, "epiCHAOS", "wcVI"),
    read_bench(opt$acvi_epichaos_bench, "epiCHAOS", "acVI")
)

tool_pal <- c(amet = "#1b9e77", scMET = "#d95f02", epiCHAOS = "#7570b3")
score_pal <- c(i_norm = "#1b9e77", jsd_norm = "#1b9e77",
               mu = "#d95f02", gamma = "#d95f02", epsilon = "#d95f02",
               eITH = "#7570b3")

bar_recovery <- ggplot(metrics, aes(x = score, y = abs(Spearman), fill = score)) +
    geom_col() +
    geom_text(aes(label = sprintf("%.2f", Spearman)), vjust = -0.3, size = 2) +
    facet_wrap(~ axis, ncol = 2) +
    scale_fill_manual(values = score_pal, guide = "none") +
    coord_cartesian(ylim = c(0, 1.05)) +
    labs(x = NULL, y = "|Spearman ρ|",
         title = "A. Recovery: amet (green), scMET (orange), epiCHAOS (purple)") +
    theme_ng() + theme(aspect.ratio = NULL,
                       axis.text.x = element_text(angle = 30, hjust = 1))

bar_cpu <- ggplot(bench, aes(x = tool, y = cpu_seconds, fill = tool)) +
    geom_col() +
    facet_wrap(~ axis, ncol = 2) +
    scale_fill_manual(values = tool_pal, guide = "none") +
    labs(x = NULL, y = "CPU seconds",
         title = "B. CPU time per tool per sweep") +
    theme_ng() + theme(aspect.ratio = NULL,
                       axis.text.x = element_text(angle = 0))

bar_mem <- ggplot(bench, aes(x = tool, y = mem_mb, fill = tool)) +
    geom_col() +
    facet_wrap(~ axis, ncol = 2) +
    scale_fill_manual(values = tool_pal, guide = "none") +
    labs(x = NULL, y = "max RSS (MB)",
         title = "C. Peak memory per tool per sweep") +
    theme_ng() + theme(aspect.ratio = NULL,
                       axis.text.x = element_text(angle = 0))

combined <- bar_recovery / (bar_cpu | bar_mem)
out_csv <- list(metrics = metrics, bench = bench)
metrics_for_csv <- rbind(
    cbind(table = "recovery", axis = metrics$axis, score = as.character(metrics$score),
          metric = "NMI",      value = metrics$NMI),
    cbind(table = "recovery", axis = metrics$axis, score = as.character(metrics$score),
          metric = "Spearman", value = metrics$Spearman),
    cbind(table = "recovery", axis = metrics$axis, score = as.character(metrics$score),
          metric = "Kendall",  value = metrics$Kendall),
    cbind(table = "benchmark", axis = bench$axis, score = bench$tool,
          metric = "cpu_seconds", value = bench$cpu_seconds),
    cbind(table = "benchmark", axis = bench$axis, score = bench$tool,
          metric = "mem_mb",      value = bench$mem_mb)
)
save_eval(combined, metrics_for_csv, opt$output_prefix, width_mm = 200, height_mm = 200)
message("[eval_tool_comparison] done")
