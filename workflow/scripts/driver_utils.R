## Shared driver categorization for amet reports.
##
## Classifies genomic annotations as "across-cell driven", "within-cell driven",
## "both", or "neither" based on how much median_jsd (across-cell) and
## median_i_total (within-cell) vary across biological groups.
##
## For each annotation, compute the SD of the group-level medians. If one SD
## is at least 1.5x the other, that component dominates. If both SDs sit
## below the 30th percentile of all annotations, the annotation is "neither".
## Otherwise it is "both".
##
## Requires dplyr.

categorize_drivers <- function(grp_df, group_col) {
  stopifnot(
    all(c("annotation", "median_i_total", "median_jsd", group_col) %in% names(grp_df))
  )

  grp_df <- grp_df[!is.na(grp_df[[group_col]]), , drop = FALSE]

  var_df <- grp_df %>%
    dplyr::group_by(annotation) %>%
    dplyr::summarise(
      jsd_sd = sd(median_jsd, na.rm = TRUE),
      i_total_sd = sd(median_i_total, na.rm = TRUE),
      .groups = "drop"
    )

  var_df$jsd_sd[!is.finite(var_df$jsd_sd)] <- 0
  var_df$i_total_sd[!is.finite(var_df$i_total_sd)] <- 0

  jsd_thr <- quantile(var_df$jsd_sd, 0.3)
  i_total_thr <- quantile(var_df$i_total_sd, 0.3)

  var_df$driver <- dplyr::case_when(
    var_df$jsd_sd < jsd_thr & var_df$i_total_sd < i_total_thr ~ "neither",
    var_df$jsd_sd >= var_df$i_total_sd * 1.5 ~ "across-cell driven",
    var_df$i_total_sd >= var_df$jsd_sd * 1.5 ~ "within-cell driven",
    TRUE ~ "both"
  )

  var_df
}

plot_driver_scatter <- function(driver_df,
                                x_label = "SD of median jsd across groups",
                                y_label = "SD of median i_total across groups") {
  ggplot2::ggplot(driver_df,
                  ggplot2::aes(x = jsd_sd, y = i_total_sd,
                               color = driver, shape = driver,
                               label = annotation)) +
    ggplot2::geom_point(size = 3) +
    ggrepel::geom_text_repel(size = 3, max.overlaps = Inf,
                             box.padding = 0.5, point.padding = 0.2,
                             min.segment.length = 0, force = 2,
                             segment.size = 0.2, segment.alpha = 0.5) +
    ggplot2::scale_color_manual(values = driver_pal) +
    ggplot2::scale_shape_manual(values = driver_shapes) +
    ggplot2::labs(x = x_label, y = y_label, color = "driver", shape = "driver") +
    theme_ng() +
    ggplot2::theme(plot.margin = ggplot2::margin(3, 6, 3, 6, "mm"))
}
