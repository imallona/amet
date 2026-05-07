## Evaluation: I_total should be flat against marginal p when there is no comethylation,
## and stay positive when there is. Reads amet's per-cell-per-feature output and plots
## I_total vs mean_meth, faceted by ground-truth structure label encoded in cell_id.

suppressPackageStartupMessages({
    library(optparse)
    library(ggplot2)
    library(dplyr)
})

source(file.path(dirname(sys.frame(1)$ofile), "plot_theme.R"))

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

agg <- df %>%
    filter(!is.na(i_total)) %>%
    group_by(condition, structure) %>%
    summarise(
        mean_meth = mean(mean_meth, na.rm = TRUE),
        i_total_mean = mean(i_total),
        i_total_sd = sd(i_total),
        n_cells = dplyr::n(),
        .groups = "drop"
    )

p <- ggplot(agg, aes(x = mean_meth, y = i_total_mean, colour = structure)) +
    geom_point(size = 1) +
    geom_errorbar(aes(ymin = i_total_mean - i_total_sd,
                      ymax = i_total_mean + i_total_sd),
                  width = 0) +
    scale_colour_manual(values = c(iid = "grey40", structured = "firebrick")) +
    labs(x = "mean methylation", y = expression(I[total] ~ "(bits)"),
         colour = NULL) +
    theme_ng()

save_eval(p, agg, opt$output_prefix, width_mm = 70, height_mm = 70)
message(sprintf("[eval_p_decoupling] wrote %s.{pdf,svg,csv}", opt$output_prefix))
