## Evaluation: I_norm should be flat against marginal p when there is no comethylation,
## and stay positive when there is. Reads amet's per-cell-per-feature output and plots
## I_norm vs mean_meth, faceted by ground-truth structure label encoded in cell_id.

suppressPackageStartupMessages({
    library(optparse)
    library(ggplot2)
    library(dplyr)
})

.this_dir <- local({ args <- commandArgs(trailingOnly = FALSE); fa <- grep("^--file=", args, value = TRUE); if (length(fa) > 0) dirname(sub("^--file=", "", fa[1])) else "." })

source(file.path(.this_dir, "plot_theme.R"))

options <- list(
    make_option(c("--cell_feature"), type = "character",
                help = "Path to amet's cell_feature.tsv.gz"),
    make_option(c("--output_prefix"), type = "character",
                help = "Output prefix (writes .pdf, .svg, .csv)")
)
opt <- parse_args(OptionParser(option_list = options))

df <- read.table(gzfile(opt$cell_feature), header = TRUE, sep = "\t",
                 na.strings = "NA", stringsAsFactors = FALSE)

## We expect cell_id to encode the ground-truth label: e.g., "iid_p0.3_seed5" or
## "markov_p11-0.8_p00-0.8_seed5". Strip the seed suffix to get the condition label.
df$condition <- sub("_seed[0-9]+$", "", df$cell_id)
df$structure <- ifelse(grepl("^iid_", df$condition), "iid", "structured")

shannon_binary <- function(p) {
    out <- numeric(length(p)); safe <- !is.na(p) & p > 0 & p < 1
    out[safe] <- -p[safe] * log2(p[safe]) - (1 - p[safe]) * log2(1 - p[safe])
    out[!safe] <- NA_real_; out
}
i_cols <- grep("^i_[0-9]+$", names(df), value = TRUE); k_max <- length(i_cols)
df$i_norm <- df$i_total / (k_max * shannon_binary(df$mean_meth))

agg <- df %>%
    filter(is.finite(i_norm)) %>%
    group_by(condition, structure) %>%
    summarise(
        mean_meth = mean(mean_meth, na.rm = TRUE),
        i_norm_mean = mean(i_norm),
        i_norm_sd = sd(i_norm),
        n_cells = dplyr::n(),
        .groups = "drop"
    )

p <- ggplot(agg, aes(x = mean_meth, y = i_norm_mean, colour = structure)) +
    geom_point(size = 1) +
    geom_errorbar(aes(ymin = i_norm_mean - i_norm_sd,
                      ymax = i_norm_mean + i_norm_sd),
                  width = 0) +
    scale_colour_manual(values = c(iid = "grey40", structured = "firebrick")) +
    labs(x = "mean methylation", y = expression(I[norm]),
         colour = NULL) +
    theme_ng()

save_eval(p, agg, opt$output_prefix, width_mm = 70, height_mm = 70)
message(sprintf("[eval_p_decoupling] wrote %s.{pdf,svg,csv}", opt$output_prefix))
