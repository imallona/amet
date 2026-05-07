## JSD finite-sample bias and variance vs cells per group.

suppressPackageStartupMessages({
    library(optparse); library(ggplot2); library(dplyr)
})
.this_dir <- local({ args <- commandArgs(trailingOnly = FALSE); fa <- grep("^--file=", args, value = TRUE); if (length(fa) > 0) dirname(sub("^--file=", "", fa[1])) else "." })

source(file.path(.this_dir, "plot_theme.R"))

options <- list(
    make_option(c("--feature"), type = "character"),
    make_option(c("--output_prefix"), type = "character")
)
opt <- parse_args(OptionParser(option_list = options))

df <- read.table(gzfile(opt$feature), header = TRUE, sep = "\t",
                 na.strings = "NA", stringsAsFactors = FALSE)
df$n <- as.integer(sub("n", "", df$group))

p <- ggplot(df, aes(x = n, y = jsd)) +
    geom_point(size = 1) +
    geom_line() +
    scale_x_log10(breaks = unique(df$n)) +
    labs(x = "n cells per group (log scale)", y = "JSD (bits)") +
    theme_ng()

save_eval(p, df, opt$output_prefix, width_mm = 70, height_mm = 70)
message(sprintf("[eval_n_cells] wrote %s.{pdf,svg,csv}", opt$output_prefix))
