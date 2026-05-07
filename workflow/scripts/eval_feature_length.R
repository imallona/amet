## I_norm mean and spread vs feature length, single panel coloured by methylation.

suppressPackageStartupMessages({
    library(optparse); library(ggplot2); library(dplyr)
})
.this_dir <- local({ args <- commandArgs(trailingOnly = FALSE); fa <- grep("^--file=", args, value = TRUE); if (length(fa) > 0) dirname(sub("^--file=", "", fa[1])) else "." })
source(file.path(.this_dir, "plot_theme.R"))

options <- list(
    make_option(c("--cell_feature"), type = "character"),
    make_option(c("--output_prefix"), type = "character")
)
opt <- parse_args(OptionParser(option_list = options))

shannon_binary <- function(p) {
    out <- numeric(length(p)); safe <- !is.na(p) & p > 0 & p < 1
    out[safe] <- -p[safe] * log2(p[safe]) - (1 - p[safe]) * log2(1 - p[safe])
    out[!safe] <- NA_real_; out
}
df <- read.table(gzfile(opt$cell_feature), header = TRUE, sep = "\t",
                 na.strings = "NA", stringsAsFactors = FALSE)
df$length <- as.integer(sub("len", "", df$feature_id))
df$target_p <- as.numeric(sub("p([0-9.]+)_seed.*", "\\1", df$cell_id))
i_cols <- grep("^i_[0-9]+$", names(df), value = TRUE); k_max <- length(i_cols)
df$i_norm <- df$i_total / (k_max * shannon_binary(df$mean_meth))

agg <- df %>%
    filter(is.finite(i_norm)) %>%
    group_by(length, target_p) %>%
    summarise(i_norm_mean = mean(i_norm),
              i_norm_sd = sd(i_norm),
              n = dplyr::n(), .groups = "drop")

p <- ggplot(agg, aes(x = length, y = i_norm_mean,
                     colour = target_p, group = factor(target_p))) +
    geom_errorbar(aes(ymin = i_norm_mean - i_norm_sd,
                      ymax = i_norm_mean + i_norm_sd),
                  width = 0, alpha = 0.6) +
    geom_line(linewidth = 0.4) +
    geom_point(size = 1.4) +
    scale_x_log10() +
    scale_colour_viridis_c(option = "inferno", limits = c(0, 1)) +
    labs(x = "feature length (CpGs, log)", y = expression(I[norm]),
         colour = "target p") +
    theme_ng() + theme(aspect.ratio = NULL)

save_eval(p, agg, opt$output_prefix, width_mm = 100, height_mm = 75)
message(sprintf("[eval_feature_length] wrote %s.{pdf,svg,csv}", opt$output_prefix))
