## Plotting theme: square panels, no grid lines, black borders.
## Used across all evaluation scripts and Rmd reports.

suppressPackageStartupMessages({
    library(ggplot2)
    library(patchwork)
})

theme_ng <- function(base_size = 8, base_family = "Helvetica") {
    theme_classic(base_size = base_size, base_family = base_family) %+replace%
        theme(
            panel.grid = element_blank(),
            panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.4),
            axis.line = element_blank(),
            axis.ticks = element_line(linewidth = 0.3, colour = "black"),
            axis.ticks.length = unit(1.5, "mm"),
            axis.text = element_text(size = base_size, colour = "black"),
            axis.title = element_text(size = base_size, colour = "black"),
            strip.background = element_blank(),
            strip.text = element_text(size = base_size, colour = "black"),
            legend.background = element_blank(),
            legend.key = element_blank(),
            legend.key.size = unit(3, "mm"),
            legend.text = element_text(size = base_size - 1),
            legend.title = element_text(size = base_size),
            plot.title = element_text(size = base_size + 1, face = "plain", hjust = 0),
            plot.subtitle = element_text(size = base_size, colour = "grey30"),
            plot.background = element_blank(),
            aspect.ratio = 1,
            complete = TRUE
        )
}

save_eval <- function(plot, data, prefix, width_mm = 60, height_mm = 60) {
    dir.create(dirname(prefix), showWarnings = FALSE, recursive = TRUE)
    ggsave(paste0(prefix, ".pdf"), plot,
           width = width_mm, height = height_mm, units = "mm",
           device = grDevices::cairo_pdf)
    ggsave(paste0(prefix, ".svg"), plot,
           width = width_mm, height = height_mm, units = "mm")
    write.table(data, paste0(prefix, ".csv"),
                sep = ",", row.names = FALSE, quote = FALSE)
}

## Two-panel score plot: score vs swept parameter, score vs mean methylation.
## The second panel exposes any unintended marginal coupling.
two_panel_score <- function(data, x_var, y_var, colour_var = NULL,
                             x_label = NULL, y_label = NULL) {
    if (is.null(x_label)) x_label <- x_var
    if (is.null(y_label)) y_label <- y_var
    p1 <- ggplot(data, aes(x = .data[[x_var]], y = .data[[y_var]])) +
        geom_point(size = 1) +
        labs(x = x_label, y = y_label) +
        theme_ng()
    p2 <- ggplot(data, aes(x = mean_meth, y = .data[[y_var]])) +
        geom_point(size = 1) +
        labs(x = "mean methylation", y = y_label) +
        theme_ng()
    if (!is.null(colour_var)) {
        p1 <- p1 + aes(colour = .data[[colour_var]])
        p2 <- p2 + aes(colour = .data[[colour_var]])
    }
    p1 + p2 + plot_layout(guides = "collect")
}
