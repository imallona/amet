## Sanity-check the data generators by reading amet's per-cell output and comparing
## the realised marginal methylation against the target encoded in cell_id, plus the
## empirical I_1 against the analytic value for Markov chains.

suppressPackageStartupMessages({
    library(optparse); library(ggplot2); library(dplyr); library(patchwork)
})
source(file.path(dirname(sys.frame(1)$ofile), "plot_theme.R"))
source(file.path(dirname(sys.frame(1)$ofile), "simulate_markov_chain.R"))

options <- list(
    make_option(c("--cell_feature"), type = "character"),
    make_option(c("--ground_truth"), type = "character",
                help = "Optional ground_truth.csv from a generator that emits one"),
    make_option(c("--output_prefix"), type = "character")
)
opt <- parse_args(OptionParser(option_list = options))

df <- read.table(gzfile(opt$cell_feature), header = TRUE, sep = "\t",
                 na.strings = "NA", stringsAsFactors = FALSE)

## Try to extract a target marginal from the cell_id, e.g. "iid_p0.30_seed5".
df$target_p <- as.numeric(sub(".*p([0-9.]+).*", "\\1",
                              ifelse(grepl("p[0-9.]+", df$cell_id), df$cell_id, NA)))

agg <- df %>%
    filter(!is.na(mean_meth)) %>%
    group_by(cell_id) %>%
    summarise(mean_meth = mean(mean_meth, na.rm = TRUE),
              i_1 = mean(i_1, na.rm = TRUE),
              target_p = first(target_p),
              .groups = "drop")

p_marg <- if (any(!is.na(agg$target_p))) {
    ggplot(agg %>% filter(!is.na(target_p)),
           aes(x = target_p, y = mean_meth)) +
        geom_abline(slope = 1, intercept = 0, linetype = 2, colour = "grey50") +
        geom_point(size = 1) +
        labs(x = "target marginal p", y = "realised mean methylation",
             title = "Marginal calibration") +
        theme_ng()
} else {
    ggplot() + labs(title = "Marginal calibration: no target_p in cell_id") + theme_ng()
}

p_i1 <- ggplot(agg, aes(x = mean_meth, y = i_1)) +
    geom_point(size = 1) +
    labs(x = "mean methylation", y = expression(I[1] ~ "(bits)"),
         title = "I_1 across cells") +
    theme_ng()

combined <- p_marg + p_i1 + plot_layout(guides = "collect")
save_eval(combined, agg, opt$output_prefix, width_mm = 120, height_mm = 60)
message(sprintf("[eval_simulator_diagnostics] wrote %s.{pdf,svg,csv}", opt$output_prefix))
