# figures.R

clrs <- c(
  "21%"    = "#E41A1C",
  "0.5%"   = "#377EB8",
  "0.2%"   = "#08306b",
  "DMSO"   = "#4DAF4A",
  "BAY"    = "#984EA3",
  "siCTL"  = "#999999",
  "siMYC"  = "#F781BF",
  "N.DMSO" = "#b2df8a",
  "H.DMSO" = "#33a02c",
  "N.BAY" =  "#cab2d6",
  "H.BAY" =  "#6a3d9a"
)

theme_plots <- function() {
  list(
    wmo::theme_wmo(
      base_family = "Calibri",
      base_size = 8
    ),
    ggplot2::theme(
      # panel.border = ggplot2::element_rect(size = 0.1),
      axis.line = ggplot2::element_line(
        colour = "black",
        linewidth = 0.25,
        lineend = "square"
      ),
      panel.border = ggplot2::element_blank(),
      plot.margin = ggplot2::margin(5, 5, 5, 5),
      plot.tag = ggplot2::element_text(face = "bold"),
      axis.title.y.left = ggplot2::element_text(margin = ggplot2::margin(r = 3))
    ),
    ggplot2::coord_cartesian(clip = "off")
  )
}

theme_patchwork <- function(
    design = NULL,
    widths = NULL,
    heights = NULL,
    tags = "A",
    ...
) {
  list(
    patchwork::plot_layout(
      design = design,
      widths = widths,
      heights = heights,
      ...
    ),
    patchwork::plot_annotation(
      tag_levels = tags,
      theme = ggplot2::theme(plot.margin = ggplot2::margin(0, 0, 0, 0))
    )
  )
}

write_figures <- function(plot, filename, path = "manuscript/figs") {
  gtab <- patchwork::patchworkGrob(plot)

  overall_width <-
    grid::convertWidth(
      sum(gtab$widths),
      unitTo = "in",
      valueOnly = TRUE
    )

  overall_height <-
    grid::convertHeight(
      sum(gtab$heights),
      unitTo = "in",
      valueOnly = TRUE
    )

  ggplot2::ggsave(
    filename = filename,
    plot = plot,
    device = ragg::agg_png,
    path = path,
    width = overall_width,
    height = overall_height,
    units = "in",
    res = 300
  )

  if (file.exists("Rplots.pdf")) unlink("Rplots.pdf")

  stringr::str_c(path, "/", filename)
}


plot_time_lines <- function(
    df,
    y,
    ylab,
    clr = c("oxygen", "treatment", "group")
) {
  ggplot2::ggplot(df) +
    ggplot2::aes(
      x = .data$time,
      y = .data[[y]],
      color = .data[[clr]],
      fill = .data[[clr]]
    ) +
    # ggplot2::geom_line(
    #   ggplot2::aes(group = interaction(date, group)),
    #   size = 0.25,
    #   alpha = 0.25,
    #   show.legend = FALSE
    # ) +
    # ggplot2::stat_summary(
    #   geom = "linerange",
    #   fun.data = ggplot2::mean_se,
    #   size = 0.5,
    #   show.legend = FALSE
  # ) +
  ggplot2::stat_summary(
    geom = "errorbar",
    fun.data = ggplot2::mean_se,
    color = "black",
    width = 2,
    linewidth = 0.25,
    show.legend = FALSE
  ) +
    ggplot2::stat_summary(
      geom = "line",
      fun.data = ggplot2::mean_se,
      size = 0.5,
      show.legend = FALSE
    ) +
    ggplot2::stat_summary(
      geom = "point",
      fun = "mean",
      pch = 21,
      color = "white",
      size = 1.5,
      stroke = 0.4,
      show.legend = TRUE
    ) +
    ggplot2::labs(
      x = "Time (h)",
      y = ylab,
      color = NULL,
      fill = NULL
    ) +
    ggplot2::scale_x_continuous(breaks = seq(0, 72, 24)) +
    ggplot2::scale_y_continuous(limits = c(0, NA)) +
    ggplot2::scale_color_manual(values = clrs) +
    ggplot2::scale_fill_manual(values = clrs) +
    ggplot2::coord_cartesian(ylim = c(0, NA)) +
    theme_plots()
}

plot_growth_curve <- function(
    df,
    cell = c("lf", "df", "pasmc"),
    exper = c("02", "05", "bay")
) {
  df |>
    dplyr::filter(
      .data$cell_type %in% cell &
        .data$experiment %in% exper &
        .data$metabolite == "cells" &
        .data$time < 96
    ) |>
    dplyr::group_by(
      .data$date,
      .data$group,
      .data$time
    ) |>
    dplyr::summarize(count = mean(.data$conc, na.rm = TRUE)) |>
    plot_time_lines(
      y = .data$count,
      ylab = "Cell count",
      clr = "group"
    ) +
    ggplot2::scale_y_continuous(
      labels = scales::label_number(scale_cut = scales::cut_short_scale())
    ) +
    ggplot2::coord_cartesian(
      ylim = c(0, NA),
      clip = "off"
    )
}

plot_growth_rates <- function(
    df,
    cell = c("lf", "pasmc"),
    exper = c("02", "05", "bay")
) {
  x <-
    df |>
    dplyr::filter(.data$cell_type %in% cell & .data$experiment %in% exper)

  annot <-
    lmerTest::lmer(mu ~ group + (1 | date), data = x) |>
    emmeans::emmeans(~ group) |>
    graphics::pairs(adjust = "mvt") |>
    broom::tidy() |>
    dplyr::mutate(
      group = stringr::str_extract(.data$contrast, "(?<= - ).*"),
      x = 1.5,
      y = Inf,
      vjust = 1,
      label = annot_p(.data$p.value)
    )

  ggplot2::ggplot(x) +
    ggplot2::aes(
      x = .data$group,
      y = .data$mu
    ) +
    ggplot2::stat_summary(
      ggplot2::aes(fill = .data$group),
      geom = "col",
      fun = "mean",
      # width = 0.6,
      show.legend = FALSE,
      alpha = 0.5
    ) +
    ggbeeswarm::geom_beeswarm(
      ggplot2::aes(fill = .data$group),
      method = "center",
      pch = 21,
      size = 1,
      stroke = 0.2,
      cex = 4,
      color = "white",
      show.legend = FALSE
    ) +
    ggplot2::stat_summary(
      geom = "errorbar",
      fun.data = ggplot2::mean_se,
      width = 0.2,
      linewidth = 0.25,
      show.legend = FALSE
    ) +
    ggplot2::geom_text(
      data = annot,
      ggplot2::aes(
        x = .data$x,
        label = .data$label,
        y = .data$y,
        vjust = .data$vjust
      ),
      family = "Calibri",
      size = 8/ggplot2::.pt
    ) +
    ggplot2::labs(
      x = "Treatment",
      y = "Growth rate (/h)"
    ) +
    ggplot2::scale_fill_manual(
      values = clrs,
      limits = force
    ) +
    theme_plots() +
    ggplot2::coord_cartesian(
      # ylim = c(0, NA),
      clip = "off"
    ) +
    NULL
}

plot_expression <- function(
    df,
    exp,
    prot = c("ldha", "hif1a"),
    ylab
) {
  df |>
    dplyr::filter(.data$experiment %in% exp) |>
    dplyr::filter(.data$protein == prot) |>
    plot_time_lines(y = "fold_change", ylab = ylab, clr = "group")
}


arrange_growth <- function(p1, p2, p3) {
  p1 + p2 + p3 +
    theme_patchwork(
      widths = ggplot2::unit(c(1.5, 1, 1), "in"),
      heights = ggplot2::unit(1, "in")
    ) &
    ggplot2::theme(
      legend.position = "bottom",
      legend.key.width = ggplot2::unit(0.5, "lines"),
      legend.key.height = ggplot2::unit(0.5, "lines"),
      legend.box.margin = ggplot2::margin(t = -10)
    )
}

arrange_viability <- function(p1) {
  p1 +
    theme_patchwork(
      widths = ggplot2::unit(1, "in"),
      heights = ggplot2::unit(1, "in"),
      tags = NULL
    ) &
    ggplot2::theme(
      legend.position = "bottom",
      legend.key.width = ggplot2::unit(0.5, "lines"),
      legend.key.height = ggplot2::unit(0.5, "lines"),
      legend.box.margin = ggplot2::margin(t = -10)
    )
}

arrange_hif_targets <- function(p1, p2, p3, p4, p5) {
  layout <- "
  ab#
  cde
  "

  p1 + p2 + p3 + p4 + p5 +
    theme_patchwork(
      design = layout,
      widths = ggplot2::unit(1, "in"),
      heights = ggplot2::unit(1, "in"),
      guides = "collect"
    ) &
    ggplot2::theme(
      legend.position = "bottom",
      legend.key.width = ggplot2::unit(0.5, "lines"),
      legend.key.height = ggplot2::unit(0.5, "lines"),
      legend.box.margin = ggplot2::margin(t = -10)
    )
}
