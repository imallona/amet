## Plotting theme: square panels, no grid lines, black borders.
## Used across all evaluation scripts and Rmd reports.
## theme_ng provides a 45-degree x-axis text variant; ng_fig_size and save_ng
## are shared helpers used by the Rmd reports.

suppressPackageStartupMessages({
    library(ggplot2)
    library(viridis)
    library(patchwork)
})

cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73",
                "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

theme_ng <- function(base_size = 8, base_family = "Helvetica") {
    theme_classic(base_size = base_size, base_family = base_family) %+replace%
        theme(
            panel.grid = element_blank(),
            panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.4),
            axis.line = element_blank(),
            axis.ticks = element_line(linewidth = 0.3, colour = "black"),
            axis.ticks.length = unit(1.5, "mm"),
            axis.text = element_text(size = base_size, colour = "black"),
            axis.text.x = element_text(size = base_size, colour = "black",
                                       angle = 45, hjust = 1, vjust = 1),
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

## Compute figure dimensions so each panel is ~panel_mm wide (and tall, given
## aspect.ratio = 1).
ng_fig_size <- function(ncol = 1, nrow = 1,
                        panel_mm  = 40,
                        legend_mm = 25,
                        strip_mm  = 6) {
    w <- (ncol * panel_mm + legend_mm + strip_mm) / 25.4
    h <- (nrow * panel_mm + strip_mm  * nrow)     / 25.4
    list(w = round(w, 1), h = round(h, 1))
}

guide_x_nolap <- function() guide_axis(check.overlap = TRUE)

theme_ng_discrete <- function(base_size = 8, base_family = "Helvetica") {
    list(
        theme_ng(base_size = base_size, base_family = base_family),
        scale_color_manual(values = cbbPalette)
    )
}

theme_ng_continuous <- function(base_size = 8, base_family = "Helvetica") {
    list(
        theme_ng(base_size = base_size, base_family = base_family),
        scale_fill_viridis_c()
    )
}

save_ng <- function(plot, file, width_mm = 85, height_mm = 60) {
    ggsave(paste0(file, ".png"), plot,
           width = width_mm, height = height_mm, units = "mm", dpi = 600)
    ggsave(paste0(file, ".svg"), plot,
           width = width_mm, height = height_mm, units = "mm")
    invisible(plot)
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
