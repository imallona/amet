## I_total should rise monotonically with comethylation persistence c at fixed marginal.
## Two panels: I_total vs persistence; I_total vs realised mean methylation (must be flat).

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
df$persistence <- as.numeric(sub("c([0-9.]+)_seed.*", "\\1", df$cell_id))

agg <- df %>%
    filter(!is.na(i_total)) %>%
    group_by(persistence) %>%
    summarise(mean_meth = mean(mean_meth, na.rm = TRUE),
              i_total_mean = mean(i_total),
              i_total_sd = sd(i_total),
              n = dplyr::n(), .groups = "drop")

p1 <- ggplot(agg, aes(x = persistence, y = i_total_mean)) +
    geom_errorbar(aes(ymin = i_total_mean - i_total_sd,
                      ymax = i_total_mean + i_total_sd), width = 0) +
    geom_point(size = 1) +
    labs(x = "persistence c (= p00 = p11)",
         y = expression(I[total] ~ "(bits)")) +
    theme_ng()
p2 <- ggplot(agg, aes(x = mean_meth, y = i_total_mean)) +
    geom_errorbar(aes(ymin = i_total_mean - i_total_sd,
                      ymax = i_total_mean + i_total_sd), width = 0) +
    geom_point(size = 1) +
    labs(x = "mean methylation", y = expression(I[total] ~ "(bits)")) +
    theme_ng()

save_eval(p1 + p2, agg, opt$output_prefix, width_mm = 120, height_mm = 60)
message(sprintf("[eval_persistence_response] wrote %s.{pdf,svg,csv}", opt$output_prefix))
