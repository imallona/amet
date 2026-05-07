## I_k profile across patterns: iid -> flat near zero, AR(1) -> I_1 high then decay,
## period_k -> spike at lag k.

suppressPackageStartupMessages({
    library(optparse); library(ggplot2); library(dplyr); library(tidyr)
})
.this_dir <- local({ args <- commandArgs(trailingOnly = FALSE); fa <- grep("^--file=", args, value = TRUE); if (length(fa) > 0) dirname(sub("^--file=", "", fa[1])) else "." })

source(file.path(.this_dir, "plot_theme.R"))

options <- list(
    make_option(c("--cell_feature"), type = "character"),
    make_option(c("--output_prefix"), type = "character")
)
opt <- parse_args(OptionParser(option_list = options))

df <- read.table(gzfile(opt$cell_feature), header = TRUE, sep = "\t",
                 na.strings = "NA", stringsAsFactors = FALSE)
df$pattern <- df$group

i_cols <- grep("^i_[0-9]+$", names(df), value = TRUE)
long <- df %>%
    filter(!is.na(i_total)) %>%
    pivot_longer(all_of(i_cols), names_to = "lag", values_to = "i_k") %>%
    mutate(lag = as.integer(sub("i_", "", lag)))

agg <- long %>%
    group_by(pattern, lag) %>%
    summarise(i_k_mean = mean(i_k, na.rm = TRUE),
              i_k_sd = sd(i_k, na.rm = TRUE),
              .groups = "drop")

p <- ggplot(agg, aes(x = lag, y = i_k_mean, colour = pattern)) +
    geom_errorbar(aes(ymin = i_k_mean - i_k_sd, ymax = i_k_mean + i_k_sd),
                  width = 0) +
    geom_line() +
    geom_point(size = 1) +
    scale_x_continuous(breaks = unique(agg$lag)) +
    labs(x = "lag k", y = expression(I[k] ~ "(bits)"),
         colour = "pattern") +
    theme_ng()

save_eval(p, agg, opt$output_prefix, width_mm = 80, height_mm = 70)
message(sprintf("[eval_lag_profile] wrote %s.{pdf,svg,csv}", opt$output_prefix))
