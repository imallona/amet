## Headline contrast: I_norm vs the marginal-only Shannon baseline H(p_obs).
## A single panel overlays both curves for both iid and structured cells, so the
## contrast between an axis-specific score (I_norm) and a marginal-only score
## (H(p)) is visible at a glance.

suppressPackageStartupMessages({
    library(optparse); library(ggplot2); library(dplyr); library(tidyr)
})
.this_dir <- local({ args <- commandArgs(trailingOnly = FALSE); fa <- grep("^--file=", args, value = TRUE); if (length(fa) > 0) dirname(sub("^--file=", "", fa[1])) else "." })

source(file.path(.this_dir, "plot_theme.R"))

shannon_h <- function(p) {
    out <- numeric(length(p))
    safe <- !is.na(p) & p > 0 & p < 1
    out[safe] <- -p[safe] * log2(p[safe]) - (1 - p[safe]) * log2(1 - p[safe])
    out
}

options <- list(
    make_option(c("--cell_feature"), type = "character"),
    make_option(c("--output_prefix"), type = "character")
)
opt <- parse_args(OptionParser(option_list = options))

df <- read.table(gzfile(opt$cell_feature), header = TRUE, sep = "\t",
                 na.strings = "NA", stringsAsFactors = FALSE)
df$shannon_marginal <- shannon_h(df$mean_meth)
i_cols <- grep("^i_[0-9]+$", names(df), value = TRUE); k_max <- length(i_cols)
df$i_norm <- df$i_total / (k_max * shannon_h(df$mean_meth))
df <- df %>% filter(is.finite(i_norm))
df$structure <- ifelse(grepl("^iid_", df$cell_id), "iid", "structured")

agg <- df %>%
    mutate(p_bin = round(mean_meth, 1)) %>%
    group_by(p_bin, structure) %>%
    summarise(mean_meth = mean(mean_meth),
              I_norm = mean(i_norm),
              `H(p)` = mean(shannon_marginal),
              .groups = "drop")

long <- agg %>%
    pivot_longer(cols = c(I_norm, `H(p)`), names_to = "score", values_to = "value")
long$score <- factor(long$score, levels = c("I_norm", "H(p)"))

p <- ggplot(long, aes(x = mean_meth, y = value,
                      colour = structure, linetype = score, shape = score)) +
    geom_line(linewidth = 0.4) +
    geom_point(size = 1.5) +
    scale_colour_manual(values = c(iid = "grey40", structured = "firebrick")) +
    scale_linetype_manual(values = c(I_norm = "solid", `H(p)` = "dashed")) +
    scale_shape_manual(values = c(I_norm = 16, `H(p)` = 1)) +
    scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25)) +
    labs(x = "mean methylation", y = "score",
         colour = "structure", linetype = "score", shape = "score",
         title = "amet I_norm (solid) vs marginal-only Shannon H(p) (dashed)") +
    theme_ng() + theme(aspect.ratio = NULL)

save_eval(p, agg, opt$output_prefix, width_mm = 110, height_mm = 80)
message(sprintf("[eval_vs_marginal_baseline] wrote %s.{pdf,svg,csv}", opt$output_prefix))
