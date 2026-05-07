## Headline contrast: I_total vs a marginal-only Shannon baseline.
## The Shannon baseline H(p_obs) is computed from the per-cell mean methylation; this
## is what a metric coupled to the marginal would look like. I_total should be flat
## against mean methylation (no inflation in the middle, no suppression at edges)
## whereas the Shannon baseline is curved.

suppressPackageStartupMessages({
    library(optparse); library(ggplot2); library(dplyr); library(patchwork)
})
source(file.path(dirname(sys.frame(1)$ofile), "plot_theme.R"))

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
df <- df %>% filter(!is.na(i_total))
df$structure <- ifelse(grepl("^iid_", df$cell_id), "iid", "structured")

agg <- df %>%
    mutate(p_bin = round(mean_meth, 1)) %>%
    group_by(p_bin, structure) %>%
    summarise(mean_meth = mean(mean_meth),
              i_total_mean = mean(i_total),
              shannon_mean = mean(shannon_marginal),
              .groups = "drop")

p1 <- ggplot(agg, aes(x = mean_meth, y = i_total_mean, colour = structure)) +
    geom_point(size = 1) +
    scale_colour_manual(values = c(iid = "grey40", structured = "firebrick")) +
    labs(x = "mean methylation",
         y = expression(I[total] ~ "(bits)"),
         title = "amet I_total") +
    theme_ng()
p2 <- ggplot(agg, aes(x = mean_meth, y = shannon_mean, colour = structure)) +
    geom_point(size = 1) +
    scale_colour_manual(values = c(iid = "grey40", structured = "firebrick")) +
    labs(x = "mean methylation", y = "H(p) (bits)",
         title = "marginal-only Shannon baseline") +
    theme_ng()

save_eval(p1 + p2 + plot_layout(guides = "collect"), agg, opt$output_prefix,
          width_mm = 130, height_mm = 65)
message(sprintf("[eval_vs_marginal_baseline] wrote %s.{pdf,svg,csv}", opt$output_prefix))
