## Recovery of the perturbation level (pl) on the consensus_perturbation_sweep.
##
## Each group has a consensus methylation pattern; cells deviate from it by
## independent CpG flips at rate (pl - 1) * 0.05. The Jaccard distance between
## cells should grow monotonically with pl, so this is the regime where epiCHAOS
## is at its strongest. amet's JSD is expected to be non-monotone (rises with
## drift, then drops as cells converge to iid noise at high flip rates).

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
spearman <- function(true_int, score) {
    keep <- is.finite(score) & is.finite(true_int)
    if (sum(keep) < 5) return(NA_real_)
    cor(true_int[keep], score[keep], method = "spearman")
}

options <- list(
    make_option(c("--cell_feature"), type = "character"),
    make_option(c("--feature"), type = "character"),
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
i_cols <- grep("^i_[0-9]+$", names(cf), value = TRUE); k_max <- length(i_cols)
cf$i_norm <- cf$i_total / (k_max * shannon_binary(cf$mean_meth))
per_group <- cf %>%
    group_by(group) %>%
    summarise(i_norm = mean(i_norm, na.rm = TRUE),
              mean_meth = mean(mean_meth, na.rm = TRUE),
              .groups = "drop")
ft$jsd_norm <- ft$jsd / (2 * shannon_binary(ft$mean_meth_mean))

scmet <- read.table(opt$scmet, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
epi   <- read.table(opt$epichaos, header = TRUE, sep = "\t", stringsAsFactors = FALSE)

joined <- per_group %>%
    inner_join(ft, by = "group", suffix = c("", ".ft")) %>%
    inner_join(scmet, by = "group") %>%
    inner_join(epi,   by = "group") %>%
    mutate(pl = as.integer(sub("pl([0-9]+)", "\\1", group)),
           flip_rate = (pl - 1L) * 0.05)

scores <- c("i_norm", "jsd", "jsd_norm", "mu", "gamma", "epsilon", "eITH")
metrics <- do.call(rbind, lapply(scores, function(s)
    data.frame(score = s, Spearman = spearman(joined$pl, joined[[s]]))))
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
         title = "A. Recovery of perturbation level pl") +
    theme_ng() + theme(aspect.ratio = NULL,
                       axis.text.x = element_text(angle = 30, hjust = 1))

long <- joined %>%
    pivot_longer(cols = all_of(scores), names_to = "score", values_to = "value")
long$score <- factor(long$score, levels = scores)
long$family <- ifelse(long$score %in% c("i_norm", "jsd", "jsd_norm"), "amet",
               ifelse(long$score %in% c("mu", "gamma", "epsilon"), "scMET", "epiCHAOS"))

trace <- ggplot(long, aes(x = pl, y = value, colour = score, group = score)) +
    geom_line(linewidth = 0.4) + geom_point(size = 1) +
    scale_colour_manual(values = score_pal) +
    facet_wrap(~ family, scales = "free_y", ncol = 3) +
    scale_x_continuous(breaks = 1:10) +
    labs(x = "perturbation level pl (flip rate = (pl-1)*0.05)",
         y = "score (per group)",
         title = "B. Score trajectories vs pl") +
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
message("[eval_consensus_perturbation] done")
