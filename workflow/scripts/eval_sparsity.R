## I_total robustness as coverage drops. Plots I_total against the realised number of
## covered CpGs and against mean methylation.

suppressPackageStartupMessages({
    library(optparse); library(ggplot2); library(dplyr); library(patchwork)
})
source(file.path(dirname(sys.frame(1)$ofile), "plot_theme.R"))

options <- list(
    make_option(c("--cell_feature"), type = "character"),
    make_option(c("--output_prefix"), type = "character")
)
opt <- parse_args(OptionParser(option_list = options))

df <- read.table(gzfile(opt$cell_feature), header = TRUE, sep = "\t",
                 na.strings = "NA", stringsAsFactors = FALSE)
df$n_fragments <- as.integer(sub("nf([0-9]+)_seed.*", "\\1", df$cell_id))

agg <- df %>%
    filter(!is.na(i_total)) %>%
    group_by(n_fragments) %>%
    summarise(n_covered_mean = mean(n_covered),
              mean_meth = mean(mean_meth, na.rm = TRUE),
              i_total_mean = mean(i_total),
              i_total_sd = sd(i_total),
              n = dplyr::n(), .groups = "drop")

p1 <- ggplot(agg, aes(x = n_covered_mean, y = i_total_mean)) +
    geom_errorbar(aes(ymin = i_total_mean - i_total_sd,
                      ymax = i_total_mean + i_total_sd), width = 0) +
    geom_point(size = 1) +
    labs(x = "mean covered CpGs", y = expression(I[total] ~ "(bits)")) +
    theme_ng()
p2 <- ggplot(agg, aes(x = mean_meth, y = i_total_mean)) +
    geom_errorbar(aes(ymin = i_total_mean - i_total_sd,
                      ymax = i_total_mean + i_total_sd), width = 0) +
    geom_point(size = 1) +
    labs(x = "mean methylation", y = expression(I[total] ~ "(bits)")) +
    theme_ng()

save_eval(p1 + p2, agg, opt$output_prefix, width_mm = 120, height_mm = 60)
message(sprintf("[eval_sparsity] wrote %s.{pdf,svg,csv}", opt$output_prefix))
