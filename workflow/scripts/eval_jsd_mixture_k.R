## JSD should be ~0 at K=1 (homogeneous) and increase with K.

suppressPackageStartupMessages({
    library(optparse); library(ggplot2); library(dplyr); library(patchwork)
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
df$K <- as.integer(sub("K", "", df$group))

p1 <- ggplot(df, aes(x = K, y = jsd)) +
    geom_point(size = 1) + geom_line() +
    scale_x_log10(breaks = unique(df$K)) +
    labs(x = "K (number of mixture components)", y = "JSD (bits)") +
    theme_ng()
p2 <- ggplot(df, aes(x = mean_meth_mean, y = jsd, colour = K)) +
    geom_point(size = 1.4) +
    scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
    scale_colour_gradient(low = "navy", high = "darkred", trans = "log10",
                          breaks = unique(df$K)) +
    labs(x = "group mean methylation", y = "JSD (bits)", colour = "K") +
    theme_ng()

save_eval(p1 + p2, df, opt$output_prefix, width_mm = 120, height_mm = 60)
message(sprintf("[eval_jsd_mixture_k] wrote %s.{pdf,svg,csv}", opt$output_prefix))
