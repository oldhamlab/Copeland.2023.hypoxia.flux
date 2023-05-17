# figures.R

clrs <- c(
  "21%"    = "#E41A1C",
  "0.5%"   = "#377EB8",
  "0.2%"   = "#08306b",
  "DMSO"   = "#4DAF4A",
  "BAY"    = "#984EA3",
  "siCTL"  = "#999999",
  "siMYC"  = "#F781BF",
  "siPHD2" = "#984EA3",
  "N.DMSO" = "#b2df8a",
  "H.DMSO" = "#33a02c",
  "N.BAY"  = "#cab2d6",
  "H.BAY"  = "#6a3d9a"
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
    filename = stringr::str_c(filename, ".png"),
    plot = plot,
    device = ragg::agg_png,
    # device = ragg::agg_tiff,
    # device = cairo_pdf,
    path = path,
    width = overall_width,
    height = overall_height,
    units = "in",
    # res = 300
  )

  if (file.exists("Rplots.pdf")) unlink("Rplots.pdf")

  stringr::str_c(path, "/", filename, ".png")
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
      show.legend = FALSE
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
      y = "count",
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

arrange_f1_s1 <- function(p1, p2, p3) {
  layout <- "
  ab#
  ccc
  "

  p1 + p2 + p3 +
    theme_patchwork(
      design = layout,
      widths = ggplot2::unit(c(1, 1, 1.25), "in"),
      heights = ggplot2::unit(1, "in")
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

plot_high_fluxes <- function(
    df,
    cell = c("lf", "pasmc"),
    exper = c("02", "05", "bay")
) {

  x <-
    df |>
    dplyr::filter(
      cell_type %in% cell &
        experiment %in% exper &
        metabolite %in% c("lactate", "glucose")
    )

  annot <-
    x |>
    dplyr::group_by(abbreviation) |>
    tidyr::nest() |>
    dplyr::mutate(
      model = purrr::map(
        data,
        ~lmerTest::lmer(flux ~ group + (1 | date), data = .x) |>
          emmeans::emmeans(~ group) |>
          pairs() |>
          broom::tidy()
      )
    ) |>
    tidyr::unnest(c(model)) |>
    dplyr::mutate(
      group = stringr::str_extract(contrast, "(?<= - ).*"),
      y = Inf,
      vjust = 1,
      label = annot_p(p.value)
    )

  ggplot2::ggplot(x) +
    ggplot2::aes(
      x = reorder(toupper(abbreviation), flux),
      y = flux
    ) +
    ggplot2::stat_summary(
      ggplot2::aes(fill = group),
      geom = "col",
      fun = "mean",
      # width = 0.6,
      position = ggplot2::position_dodge2(),
      show.legend = TRUE,
      alpha = 0.5
    ) +
    ggplot2::geom_hline(
      yintercept = 0,
      color = "black",
      lwd = 0.25
    ) +
    ggbeeswarm::geom_beeswarm(
      ggplot2::aes(fill = group),
      method = "center",
      dodge.width = 0.9,
      pch = 21,
      size = 1,
      stroke = 0.2,
      cex = 4,
      color = "white",
      show.legend = FALSE
    ) +
    ggplot2::stat_summary(
      ggplot2::aes(group = group),
      geom = "errorbar",
      fun.data = ggplot2::mean_se,
      position = ggplot2::position_dodge(width = 0.9),
      width = 0.3,
      size = 0.25,
      show.legend = FALSE
    ) +
    ggplot2::geom_text(
      data = annot,
      ggplot2::aes(
        x = abbreviation,
        label = label,
        y = y,
        vjust = vjust
      ),
      family = "Calibri",
      size = 8/ggplot2::.pt,
      inherit.aes = FALSE
    ) +
    ggplot2::labs(
      x = NULL,
      y = "Flux (fmol/cell/h)",
      fill = NULL
    ) +
    ggplot2::scale_fill_manual(values = clrs, limits = force) +
    ggplot2::scale_y_continuous(
      expand = ggplot2::expansion(mult = 0.2),
      breaks = scales::extended_breaks(n = 7)
    ) +
    ggplot2::guides(
      fill = ggplot2::guide_legend(override.aes = list(alpha = 1))
    ) +
    theme_plots() +
    ggplot2::theme(
      legend.key.width = ggplot2::unit(0.5, "lines"),
      legend.key.height = ggplot2::unit(0.5, "lines"),
      legend.position = "bottom",
      legend.box.margin = ggplot2::margin(t = -10)
    )
}

plot_low_fluxes <- function(
    df,
    cell = c("lf", "pasmc"),
    exper = c("02", "05", "bay")
) {

  x <-
    df |>
    dplyr::filter(
      cell_type %in% cell &
        experiment %in% exper &
        metabolite %nin% c("lactate", "glucose")
    ) |>
    dplyr::filter(
      !(metabolite == "glutamine" & flux > 0)
    ) |>
    dplyr::mutate(width = 0.3)

  annot <-
    x |>
    dplyr::group_by(abbreviation) |>
    tidyr::nest() |>
    dplyr::mutate(
      model = purrr::map(
        data,
        ~lmerTest::lmer(flux ~ group + (1 | date), data = .x) |>
          emmeans::emmeans(~ group) |>
          pairs() |>
          broom::tidy()
      )
    ) |>
    tidyr::unnest(c(model)) |>
    dplyr::mutate(
      group = stringr::str_extract(contrast, "(?<= - ).*"),
      y = Inf,
      vjust = 1,
      label = annot_p(p.value)
    )

  ggplot2::ggplot(x) +
    ggplot2::aes(
      x = reorder(toupper(abbreviation), flux),
      y = flux
    ) +
    ggplot2::stat_summary(
      ggplot2::aes(fill = group),
      geom = "col",
      fun = "mean",
      position = ggplot2::position_dodge2(width = 1),
      alpha = 0.5,
      show.legend = FALSE
    ) +
    ggplot2::geom_hline(
      yintercept = 0,
      color = "black",
      lwd = 0.25
    ) +
    ggplot2::stat_summary(
      ggplot2::aes(
        group = group,
        width = width
      ),
      geom = "errorbar",
      fun.data = ggplot2::mean_se,
      position = ggplot2::position_dodge(width = 0.9),
      width = 0.35,
      linewidth = 0.25,
      show.legend = FALSE
    ) +
    ggplot2::geom_text(
      data = annot,
      ggplot2::aes(
        x = abbreviation,
        y = y,
        vjust = vjust,
        label = label,
        group = group
      ),
      family = "Calibri",
      size = 8/ggplot2::.pt
    ) +
    ggplot2::labs(
      x = NULL,
      y = "Flux (fmol/cell/h)",
      fill = NULL
    ) +
    ggplot2::scale_fill_manual(values = clrs, limits = force) +
    ggplot2::scale_y_continuous(
      trans = ggallin::pseudolog10_trans,
      breaks = c(-100, -10, 0, 10),
      limits = c(-250, 50)
    ) +
    ggplot2::guides(
      fill = ggplot2::guide_legend(override.aes = list(alpha = 1))
    ) +
    theme_plots() +
    ggplot2::theme(
      panel.grid.major.x = ggplot2::element_line(color = "gray80", size = 0.1),
      legend.key.width = ggplot2::unit(0.5, "lines"),
      legend.key.height = ggplot2::unit(0.5, "lines"),
      legend.position = "bottom",
      legend.box.margin = ggplot2::margin(t = -10)
    )
}

arrange_fluxes <- function(p1, p2, p3, p4, p5, p6, p7, p8, p9, p10) {
  layout <- "
  abc
  def
  ghi
  jjj
  "

  p1 + p2 + p3 + p4 + p5 + p6 + p7 + p8 + p9 + p10 +
    theme_patchwork(
      design = layout,
      widths = unit(1, "in"),
      heights = unit(1, "in"),
      guides = "collect"
    ) &
    theme(
      legend.position = "bottom",
      legend.key.width = ggplot2::unit(0.5, "lines"),
      legend.key.height = ggplot2::unit(0.5, "lines"),
      legend.box.margin = ggplot2::margin(t = -10)
    )
}

plot_glc6_curve <- function(conc_std) {
  conc_std |>
    dplyr::ungroup() |>
    dplyr::filter(experiment == "substrate" & metabolite == "lactate") |>
    dplyr::select(date, data) |>
    tidyr::unnest(c("data")) |>
    ggplot2::ggplot() +
    ggplot2::aes(
      x = conc / 1000,
      y = value
    ) +
    ggplot2::geom_smooth(
      method = "lm",
      formula = y ~ x,
      color = clrs[[2]],
      linewidth = 0.25,
      se = FALSE
    ) +
    ggplot2::geom_point(
      pch = 21,
      color = "white",
      fill = "black",
      size = 1,
      stroke = 0.25,
      show.legend = FALSE
    ) +
    ggplot2::labs(
      x = "Lactate (mM)",
      y = "Peak area ratio"
    ) +
    theme_plots() +
    ggplot2::coord_cartesian(xlim = c(0, NA), clip = "off") +
    NULL
}

plot_glc6_mass <- function(flux_measurements) {
  flux_measurements |>
    dplyr::filter(experiment == "substrate" & metabolite == "lactate") |>
    dplyr::group_by(dplyr::across(cell_type:time)) |>
    wmo::remove_nested_outliers("nmol", remove = TRUE) |>
    dplyr::summarise(umol = mean(nmol) / 1000) |>
    plot_time_lines(y = "umol", ylab = bquote("[U-"^13 * "C"[3] * "]-Lactate (μmol)"), clr = "oxygen")
}

plot_glc6_fluxes <- function(df) {
  x <-
    df |>
    dplyr::filter(experiment == "substrate" & metabolite == "lactate")

  annot <-
    lmerTest::lmer(flux ~ group + (1 | date), data = x) |>
    emmeans::emmeans(~ group) |>
    pairs() |>
    broom::tidy()|>
    dplyr::mutate(
      group = stringr::str_extract(contrast, "(?<= - ).*"),
      y = Inf,
      vjust = 1,
      label = annot_p(p.value),
      abbreviation = "LAC"
    )

  ggplot2::ggplot(x) +
    ggplot2::aes(
      x = group,
      y = flux
    ) +
    ggplot2::stat_summary(
      ggplot2::aes(fill = group),
      geom = "col",
      fun = "mean",
      # width = 0.6,
      position = ggplot2::position_dodge2(),
      show.legend = FALSE,
      alpha = 0.5
    ) +
    ggbeeswarm::geom_beeswarm(
      ggplot2::aes(fill = group),
      method = "center",
      dodge.width = 0.9,
      pch = 21,
      size = 1,
      stroke = 0.2,
      cex = 4,
      color = "white",
      show.legend = FALSE
    ) +
    ggplot2::stat_summary(
      ggplot2::aes(group = group),
      geom = "errorbar",
      fun.data = ggplot2::mean_se,
      position = ggplot2::position_dodge(width = 0.9),
      width = 0.3,
      size = 0.25,
      show.legend = FALSE
    ) +
    ggplot2::geom_text(
      data = annot,
      ggplot2::aes(
        x = 1.5,
        label = label,
        y = y,
        vjust = vjust
      ),
      family = "Calibri",
      size = 8/ggplot2::.pt,
      inherit.aes = FALSE
    ) +
    ggplot2::labs(
      x = "Treatment",
      title = bquote("[U-"^13 * "C"[6] * "]-GLC → [U-"^13 * "C"[3] * "]-LAC"),
      y = "Flux (fmol/cell/h)",
      fill = NULL
    ) +
    ggplot2::scale_fill_manual(values = clrs, limits = force) +
    ggplot2::scale_y_continuous(
      breaks = scales::extended_breaks(n = 7)
    ) +
    ggplot2::guides(
      fill = ggplot2::guide_legend(override.aes = list(alpha = 1))
    ) +
    theme_plots() +
    ggplot2::theme(
      legend.key.width = ggplot2::unit(0.5, "lines"),
      legend.key.height = ggplot2::unit(0.5, "lines"),
      legend.position = "bottom",
      legend.box.margin = ggplot2::margin(t = -10),
      plot.title = ggplot2::element_text(size = 8, hjust = 1),
      plot.title.position = "plot"
    )
}

arrange_glc6_flux <- function(p1, p2, p3) {
  layout <- "abc"

  p1 + p2 + p3 +
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

plot_evap_data <- function(evap_clean) {
  evap_clean |>
    dplyr::filter(experiment == "05" & cell_type == "lf") |>
    dplyr::mutate(
      oxygen = factor(oxygen, levels = c("21%", "0.5%"))
    ) |>
    ggplot2::ggplot() +
    ggplot2::aes(
      x = time,
      y = volume,
      color = oxygen,
      fill = oxygen
    ) +
    ggplot2::geom_smooth(
      method = "lm",
      formula = y ~ x,
      se = FALSE,
      size = 0.5,
      show.legend = FALSE
    ) +
    # ggplot2::stat_summary(
    #   geom = "linerange",
    #   fun.data = "mean_se",
    #   size = 0.5,
    #   show.legend = FALSE
    # ) +
    ggplot2::stat_summary(
      geom = "errorbar",
      fun.data = ggplot2::mean_se,
      color = "black",
      width = 2,
      size = 0.25,
      show.legend = FALSE
    ) +
    ggplot2::stat_summary(
      geom = "point",
      fun = "mean",
      pch = 21,
      color = "white",
      size = 1.5,
      stroke = 0.2,
      show.legend = FALSE
    ) +
    ggplot2::labs(
      x = "Time (h)",
      y = "Volume (mL)"
    ) +
    ggplot2::scale_x_continuous(breaks = seq(0, 96, 24)) +
    ggplot2::scale_color_manual(values = clrs) +
    ggplot2::scale_fill_manual(values = clrs) +
    theme_plots()
}

plot_k <- function(degradation_rates, k) {
  annot <-
    k |>
    dplyr::select(-k) |>
    dplyr::mutate(
      label = "*",
      ypos = Inf,
      vjust = 1
    )

  degradation_rates |>
    dplyr::mutate(
      group = dplyr::case_when(
        oxygen == "21%" & treatment == "None" ~ "21%",
        oxygen == "0.5%" ~ "0.5%",
        treatment == "DMSO" ~ "DMSO"
      ),
      group = factor(group, levels = c("21%", "0.5%", "DMSO"))
    ) |>
    dplyr::left_join(annot, by = c("metabolite", "oxygen", "treatment")) |>
    ggplot2::ggplot() +
    ggplot2::aes(
      x = reorder(toupper(abbreviation), k),
      y = k
    ) +
    ggplot2::stat_summary(
      ggplot2::aes(fill = group),
      geom = "col",
      fun = "mean",
      position = ggplot2::position_dodge(width = 0.6),
      width = 0.6,
      show.legend = TRUE
    ) +
    ggplot2::geom_hline(
      yintercept = 0,
      color = "black",
      lwd = 0.25
    ) +
    ggplot2::stat_summary(
      ggplot2::aes(group = group),
      geom = "errorbar",
      fun.data = ggplot2::mean_se,
      position = ggplot2::position_dodge(width = 0.6),
      width = 0.2,
      size = 0.25,
      show.legend = FALSE
    ) +
    ggplot2::geom_text(
      ggplot2::aes(
        color = group,
        label = label,
        y = ypos,
        vjust = vjust
      ),
      family = "Calibri",
      size = 6/ggplot2::.pt,
      position = ggplot2::position_dodge(width = 0.6),
      show.legend = FALSE
    ) +
    ggplot2::labs(
      x = "Metabolite",
      y = "Rate constant (/h)",
      fill = NULL,
      color = NULL
    ) +
    ggplot2::scale_color_manual(
      values = clrs,
      limits = force
    ) +
    ggplot2::scale_fill_manual(
      values = clrs,
      limits = force
    ) +
    ggplot2::guides(
      fill = ggplot2::guide_legend(override.aes = list(alpha = 1))
    ) +
    # ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = 0.2)) +
    theme_plots() +
    ggplot2::theme(
      panel.grid.major.x = ggplot2::element_line(color = "gray90"),
      legend.key.width = ggplot2::unit(0.5, "lines"),
      legend.key.height = ggplot2::unit(0.5, "lines"),
      legend.position = "bottom",
      legend.box.margin = ggplot2::margin(t = -10)
    )
}

plot_substrate <- function(df, cell, title) {
  x <-
    df |>
    dplyr::filter(experiment == "substrate" & treatment != "GLC6") |>
    dplyr::filter(cell_type == cell) |>
    dplyr::ungroup() |>
    dplyr::mutate(grand_mean = mean(.data$mu)) |>
    dplyr::group_by(date) |>
    dplyr::mutate(
      exp_mean = mean(.data$mu),
      cf = grand_mean - exp_mean,
      mu_corr = mu + cf
    ) |>
    dplyr::group_by(oxygen, treatment) |>
    wmo::remove_nested_outliers("mu_corr", remove = TRUE)

  annot <-
    lmerTest::lmer(mu ~ oxygen * treatment + (1 | date), data = x) |>
    emmeans::emmeans(~ oxygen * treatment) |>
    graphics::pairs(simple = "each", combine = TRUE, adjust = "mvt") |>
    broom::tidy() |>
    dplyr::select("contrast", "treatment", "oxygen", pval = "adj.p.value") |>
    dplyr::mutate(
      lab = annot_p(.data$pval),
      y = Inf,
      vjust = 1
    )

  annot1 <-
    dplyr::filter(annot, .data$oxygen != ".") |>
    dplyr::mutate(
      treatment = rep(c("-GLC", "-GLN", "-GLN"), 2),
      treatment = factor(.data$treatment, levels = c("None", "-GLC", "-GLN")),
      oxygen = factor(.data$oxygen, levels = c("21%", "0.5%")),

    )

  annot2 <-
    dplyr::filter(annot, .data$treatment != ".") |>
    dplyr::mutate(
      oxygen = "21%",
      treatment = factor(.data$treatment, levels = c("None", "-GLC", "-GLN")),
      oxygen = factor(.data$oxygen, levels = c("21%", "0.5%"))
    )

  ggplot2::ggplot(x) +
    ggplot2::aes(
      x = .data$treatment,
      y = .data$mu_corr,
      fill = .data$oxygen
    ) +
    ggplot2::stat_summary(
      geom = "col",
      fun = "mean",
      # width = 0.6,
      position = ggplot2::position_dodge2(),
      show.legend = TRUE,
      alpha = 0.5
    ) +
    ggbeeswarm::geom_beeswarm(
      method = "center",
      dodge.width = 0.9,
      pch = 21,
      size = 1,
      stroke = 0.2,
      cex = 4,
      color = "white",
      show.legend = FALSE
    ) +
    ggplot2::stat_summary(
      # ggplot2::aes(group = treatment),
      geom = "errorbar",
      fun.data = ggplot2::mean_se,
      position = ggplot2::position_dodge(width = 0.9),
      width = 0.2,
      linewidth = 0.25,
      show.legend = FALSE
    ) +
    ggplot2::geom_text(
      data = annot1,
      ggplot2::aes(
        color = .data$oxygen,
        y = .data$y,
        label = .data$lab,
        vjust = .data$vjust
      ),
      family = "Calibri",
      size = 8/ggplot2::.pt,
      position = ggplot2::position_dodge(width = 0.9),
      show.legend = FALSE
    ) +
    ggplot2::geom_text(
      data = annot2,
      ggplot2::aes(
        y = .data$y,
        label = .data$lab,
        vjust = .data$vjust
      ),
      family = "Calibri",
      size = 8/ggplot2::.pt,
      color = "black",
      show.legend = FALSE
    ) +
    ggplot2::scale_y_continuous(
      limits = c(0, 0.03)
    ) +
    ggplot2::scale_color_manual(
      values = clrs,
      limits = force,
      aesthetics = c("fill", "color")
    ) +
    ggplot2::labs(
      x = "Treatment",
      y = "Growth rate (/h)",
      fill = NULL,
      color = NULL,
      title = title
    ) +
    ggplot2::guides(
      fill = ggplot2::guide_legend(override.aes = list(alpha = 1))
    ) +
    theme_plots() +
    ggplot2::theme(
      legend.key.width = ggplot2::unit(0.5, "lines"),
      legend.key.height = ggplot2::unit(0.5, "lines"),
      legend.position = "bottom",
      legend.box.margin = ggplot2::margin(t = -10)
    )
}

arrange_substrate <- function(p1, p2) {
  layout <- "ab"

  p1 + p2 +
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

plot_all_mids <- function(df, cell, t = 72, title) {
  tracer_labels <- c(
    expression(paste("[1,2-"^13, "C"[2], "]-glucose")),
    expression(paste("[U-"^13, "C"[6], "]-glucose")),
    expression(paste("[U-"^13, "C"[5], "]-glutamine")),
    expression(paste("[U-"^13, "C"[3], "]-lactate"))
  )
  tracer_levels <- c("glc2", "glc6", "q5", "lac3")
  metab <- c(
    "FBP",
    "3PG",
    "PYR",
    # "ALA",
    # "SER",
    "LAC",
    "CIT",
    "AKG",
    # "GLN",
    "GLU",
    "MAL"
    # "ASP"
  )

  x <-
    df |>
    dplyr::filter(
      cell_type == cell &
        time == t &
        metabolite %in% metab
    ) |>
    dplyr::filter(
      !(metabolite == "CIT" & isotope == "M6")
    ) |>
    dplyr::mutate(
      metabolite = factor(metabolite, levels = metab),
      metabolite = forcats::fct_recode(metabolite, "`3PG`" = "3PG"),
      tracer = factor(tracer, levels = tracer_levels, labels = tracer_labels),
      isotope = stringr::str_replace(isotope, "M", ""),
      group = dplyr::case_when(
        oxygen == "21%" & treatment == "None" ~ "21%",
        oxygen == "0.5%" & treatment == "None" ~ "0.5%",
        oxygen == "21%" & treatment == "DMSO" ~ "DMSO",
        oxygen == "21%" & treatment == "BAY" ~ "BAY"
      ),
      group = factor(group, levels = c("21%", "0.5%", "DMSO", "BAY"))
    )

  annot <-
    x |>
    dplyr::group_by(tracer, metabolite) |>
    tidyr::nest() |>
    dplyr::mutate(annot = purrr::map(data, annot_mids_main)) |>
    tidyr::unnest(c(annot))

  ggplot2::ggplot(x) +
    ggplot2::aes(
      x = isotope,
      y = mid
    ) +
    ggplot2::facet_grid(
      tracer ~ metabolite,
      labeller = ggplot2::label_parsed,
      # switch = "y",
      scales = "free_x",
      space = "free_x"
    ) +
    ggplot2::stat_summary(
      ggplot2::aes(fill = group),
      geom = "col",
      fun = "mean",
      position = ggplot2::position_dodge(width = 0.6),
      width = 0.6,
      show.legend = TRUE
    ) +
    ggplot2::geom_hline(
      yintercept = 0,
      color = "black",
      lwd = 0.1
    ) +
    ggplot2::stat_summary(
      ggplot2::aes(group = group),
      geom = "errorbar",
      fun.data = "mean_se",
      position = ggplot2::position_dodge(width = 0.6),
      width = 0.2,
      size = 0.25,
      show.legend = FALSE
    ) +
    ggplot2::geom_text(
      data = annot,
      ggplot2::aes(
        x = -Inf,
        y = Inf,
        vjust = 1.5,
        hjust = -0.3,
        label = annot
      ),
      family = "Calibri",
      size = 8/ggplot2::.pt
    ) +
    ggplot2::labs(
      x = "Isotope",
      y = "Mole fraction",
      fill = NULL,
      title = title
    ) +
    ggplot2::scale_y_continuous(
      breaks = seq(0, 1, 0.25),
      limits = c(0, 1.1)
    ) +
    ggplot2::theme(
      strip.placement = "outside",
      legend.title = ggplot2::element_blank()
    ) +
    ggplot2::scale_fill_manual(values = clrs, limits = force) +
    theme_plots() +
    ggplot2::theme(
      panel.grid.major.y = ggplot2::element_line(color = "grey80", size = 0.1),
      legend.key.width = ggplot2::unit(0.5, "lines"),
      legend.key.height = ggplot2::unit(0.5, "lines"),
      legend.position = "bottom",
      legend.box.margin = ggplot2::margin(t = -10),
      axis.line = ggplot2::element_blank(),
      panel.border = ggplot2::element_rect(color = "black", fill = NA, size = 0.25)
    ) +
    theme_patchwork(
      widths = unit(9.5, "in"),
      heights = unit(7.5, "in"),
      tags = NULL
    )
}

annot_mids_main <- function(a, formula) {
  lmerTest::lmer(
    mid ~ group * isotope + (1 | date),
    data = a
  ) |>
    emmeans::emmeans(~ group * isotope) |>
    emmeans::mvcontrast("pairwise", mult.name = "isotope") |>
    tibble::as_tibble() |>
    dplyr::select(contrast, tidyselect::contains("p.value")) |>
    dplyr::rename_with(~ "pval", .cols = tidyselect::contains("p.value")) |>
    dplyr::mutate(
      label = dplyr::case_when(
        pval < 0.05 & contrast == "21% - 0.5%" ~ "*",
        pval < 0.05 & contrast == "21% - 0.2%" ~ "*",
        pval < 0.05 & contrast == "DMSO - BAY" ~ "†",
        pval < 0.05 & contrast == "0.5% - BAY" ~ "‡"
      ),
      order = dplyr::case_when(
        contrast == "21% - 0.5%" ~ 1,
        contrast == "21% - 0.2%" ~ 1,
        contrast == "DMSO - BAY" ~ 2,
        contrast == "0.5% - BAY" ~ 3
      )
    ) |>
    dplyr::filter(!is.na(label)) |>
    dplyr::arrange(order) |>
    dplyr::pull(label) |>
    stringr::str_c(collapse = " ")
}

plot_m5_citrate <- function(df) {
  x <-
    df |>
    dplyr::filter(
      metabolite == "CIT" &
        tracer == "q5" &
        treatment == "None" &
        isotope == "M5" &
        method == "sim"
    ) |>
    dplyr::filter(
      (cell_type == "lf" & time == 48) |
        (cell_type == "pasmc" & time == 36)
    )

  annot <-
    x |>
    dplyr::group_by(cell_type) |>
    tidyr::nest() |>
    dplyr::mutate(
      model = purrr::map(
        data,
        ~lmerTest::lmer(mid ~ oxygen + (1 | date), data = .x) |>
          emmeans::emmeans(~ oxygen) |>
          pairs() |>
          broom::tidy()
      )
    ) |>
    tidyr::unnest(c(model)) |>
    dplyr::mutate(
      x = 1.5,
      y = Inf,
      vjust = 1.5,
      label = annot_p(p.value)
    )

  ggplot2::ggplot(x) +
    ggplot2::aes(
      x = oxygen,
      y = mid
    ) +
    ggplot2::facet_wrap(
      ~cell_type,
      labeller = ggplot2::as_labeller(toupper)
    ) +
    ggplot2::stat_summary(
      ggplot2::aes(fill = oxygen),
      geom = "col",
      fun = "mean",
      # width = 0.6,
      show.legend = FALSE,
      alpha = 0.5
    ) +
    ggbeeswarm::geom_beeswarm(
      ggplot2::aes(fill = oxygen),
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
      size = 0.25,
      show.legend = FALSE
    ) +
    ggplot2::geom_text(
      data = annot,
      ggplot2::aes(
        x = x,
        label = label,
        y = y,
        vjust = vjust
      ),
      family = "Calibri",
      size = 8/ggplot2::.pt
    ) +
    ggplot2::labs(
      x = "Oxygen",
      y = "M5 Citrate fraction"
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
    ggplot2::theme(
      axis.line = ggplot2::element_blank(),
      panel.border = ggplot2::element_rect(color = "black", fill = NA, size = 0.25)
    ) +
    NULL
}

plot_mids <- function(df, cell, metab, t = 72, track) {
  tracer_labels <- c(
    glc2 = expression(paste("[1,2-"^13, "C"[2], "]-GLC")),
    glc6 = bquote("[U-"^13 * "C"[6] * "]-GLC →" ~ .(metab)),
    q5 = bquote("[U-"^13 * "C"[5] * "]-GLN →" ~ .(metab)),
    lac3 = expression(paste("[U-"^13, "C"[3], "]-LAC"))
  )
  tracer_levels <- c("glc2", "glc6", "q5", "lac3")

  x <-
    df |>
    dplyr::filter(
      cell_type == cell &
        metabolite == metab &
        time == t &
        tracer == track
    ) |>
    dplyr::filter(
      !(metabolite == "CIT" & isotope == "M6")
    ) |>
    dplyr::mutate(
      metabolite = factor(metabolite, levels = metab),
      tracer = factor(tracer, levels = tracer_levels, labels = tracer_labels),
      group = dplyr::case_when(
        oxygen == "21%" & treatment == "None" ~ "21%",
        oxygen == "0.5%" & treatment == "None" ~ "0.5%",
        oxygen == "21%" & treatment == "DMSO" ~ "DMSO",
        oxygen == "21%" & treatment == "BAY" ~ "BAY"
      ),
      group = factor(group, levels = c("21%", "0.5%", "DMSO", "BAY"))
    )

  annot <-
    x |>
    dplyr::group_by(tracer, metabolite) |>
    tidyr::nest() |>
    dplyr::mutate(annot = purrr::map(data, annot_mids_main)) |>
    tidyr::unnest(c(annot))

  ggplot2::ggplot(x) +
    ggplot2::aes(
      x = isotope,
      y = mid
    ) +
    ggplot2::stat_summary(
      ggplot2::aes(fill = group),
      geom = "col",
      fun = "mean",
      position = ggplot2::position_dodge(width = 0.6),
      width = 0.6,
      show.legend = TRUE
    ) +
    ggplot2::geom_hline(
      yintercept = 0,
      color = "black",
      lwd = 0.1
    ) +
    ggplot2::stat_summary(
      ggplot2::aes(group = group),
      geom = "errorbar",
      fun.data = "mean_se",
      position = ggplot2::position_dodge(width = 0.6),
      width = 0.2,
      size = 0.25,
      show.legend = FALSE
    ) +
    ggplot2::geom_text(
      data = annot,
      ggplot2::aes(
        x = -Inf,
        y = Inf,
        vjust = 1.5,
        hjust = -0.3,
        label = annot
      ),
      family = "Calibri",
      size = 8/ggplot2::.pt
    ) +
    ggplot2::labs(
      x = "Isotope",
      y = "Mole fraction",
      fill = NULL,
      title = tracer_labels[track]
    ) +
    ggplot2::scale_y_continuous(
      breaks = seq(0, 1, 0.25),
      limits = c(0, 1.1)
    ) +
    ggplot2::theme(
      strip.placement = "outside",
      legend.title = ggplot2::element_blank()
    ) +
    ggplot2::scale_fill_manual(values = clrs, limits = force) +
    theme_plots() +
    ggplot2::theme(
      panel.grid.major.y = ggplot2::element_line(color = "grey80", size = 0.1),
      legend.key.width = ggplot2::unit(0.5, "lines"),
      legend.key.height = ggplot2::unit(0.5, "lines"),
      legend.position = "bottom",
      legend.box.margin = ggplot2::margin(t = -10)
    )
}

arrange_f4 <- function(p1, p2, p3, p4) {
  layout <- "
  ab
  cd
  "

  p1 + p2 + p3 + p4 +
    theme_patchwork(
      design = layout,
      widths = unit(1.75, "in"),
      heights = unit(1, "in"),
      guides = "collect"
    ) &
    theme(
      legend.position = "bottom",
      legend.key.width = ggplot2::unit(0.5, "lines"),
      legend.key.height = ggplot2::unit(0.5, "lines"),
      legend.box.margin = ggplot2::margin(t = -10)
    )
}

format_time_course_mids <- function(model_mids) {
  tracer_labels <- c(
    expression(paste("[1,2-"^13, "C"[2], "]-glucose")),
    expression(paste("[U-"^13, "C"[6], "]-glucose")),
    expression(paste("[U-"^13, "C"[5], "]-glutamine")),
    expression(paste("[U-"^13, "C"[3], "]-lactate"))
  )
  tracer_levels <- c("glc2", "glc6", "q5", "lac3")

  model_mids |>
    tidyr::unnest(c(data)) |>
    dplyr::filter(metabolite %in% c("FBP", "PYR", "CIT", "MAL")) |>
    dplyr::filter(tracer != "lac3") |>
    # dplyr::filter(isotope != "M6") |>
    dplyr::mutate(
      tracer = factor(
        tracer,
        levels = tracer_levels,
        labels = tracer_labels
      ),
      metabolite = factor(metabolite, levels = c("FBP", "PYR", "CIT", "MAL"))
    )
}

plot_mid_time_course <- function(df, cells, o2, treat, color) {
  if (cells == "lf") {
    brks <- seq(0, 72, 24)
  } else if (cells == "pasmc") {
    brks <- seq(0, 48, 12)
  }

  df |>
    dplyr::filter(
      cell_type == cells &
        oxygen %in% o2 &
        treatment == treat &
        isotope == "M0"
    ) |>
    ggplot2::ggplot() +
    ggplot2::aes(
      x = time,
      y = 1 - mean,
      color = oxygen,
      fill = oxygen
    ) +
    ggplot2::facet_grid(
      metabolite ~ tracer,
      labeller = ggplot2::label_parsed
    ) +
    ggplot2::geom_linerange(
      ggplot2::aes(
        ymin = 1 - mean - se,
        ymax = 1 - mean + se
      ),
      color = "black",
      size = 0.25,
      show.legend = FALSE
    ) +
    ggplot2::geom_line(
      size = 0.5,
      show.legend = FALSE
    ) +
    ggplot2::geom_point(
      pch = 21,
      color = "white",
      size = 1.5,
      stroke = 0.2,
      show.legend = TRUE
    ) +
    ggplot2::labs(
      x = "Time (h)",
      y = "Labeled fraction",
      color = NULL,
      fill = NULL
    ) +
    ggplot2::scale_x_continuous(breaks = brks) +
    ggplot2::scale_color_manual(
      values = clrs,
      limits = force,
      aesthetics = c("fill", "color")
    ) +
    ggplot2::coord_cartesian(ylim = c(0, NA)) +
    theme_plots() +
    ggplot2::theme(
      panel.grid.major.y = ggplot2::element_line(color = "grey80", size = 0.1),
      legend.key.size = ggplot2::unit(0.5, units = "lines"),
      axis.line = ggplot2::element_blank(),
      panel.border = ggplot2::element_rect(color = "black", fill = NA, size = 0.25),
      legend.position = "bottom",
      legend.box.margin = ggplot2::margin(t = -10)
    ) +
    theme_patchwork(
      widths = unit(3, "in"),
      heights = unit(3, "in"),
      tags = NULL
    )
}

arrange_f5_s2 <- function(p1, p2) {
  layout <- "ab"

  p1 + p2 +
    theme_patchwork(
      design = layout,
      widths = unit(3.5, "in"),
      heights = unit(3.5, "in"),
      guides = "collect"
    )
}

format_flux_table <- function(
    flux_differences,
    cell = c("lf", "pasmc"),
    experiment = c("0.5%", "BAY"),
    ssr_ctl = NULL,
    ssr_exp = NULL
) {

  big_border <- flextable::fp_border_default(color = "black", width = 1)
  small_border <- flextable::fp_border_default(color = "black", width = 0.25)

  conditions <- c(unique(flux_differences$ctl), unique(flux_differences$exp))

  df <-
    flux_differences |>
    dplyr::ungroup() |>
    dplyr::filter(normalization == "none" & cell_type == cell & exp == experiment) |>
    dplyr::mutate(
      equation = stringr::str_replace_all(.data$equation, "<->", "↔︎"),
      equation = stringr::str_replace_all(.data$equation, "->", "→"),
      equation = stringr::str_replace_all(.data$equation, "<-", "←")
    )

  conditions <- c(unique(df$ctl), unique(df$exp))

  df |>
    dplyr::select(-c(normalization, cell_type, index, ctl, exp)) |>
    dplyr::arrange(desc(type)) |>
    dplyr::mutate(
      pathway = stringr::str_to_sentence(pathway),
      type = toupper(type),
      dplyr::across(tidyselect::matches("ctl|exp"), ~scales::scientific(.x))
    ) |>
    dplyr::select("type", "pathway", tidyselect::everything()) |>
    flextable::flextable() |>
    flextable::merge_v(j = c("type", "pathway")) |>
    flextable::border_remove() |>
    # flextable::bold(j = 1, i = ~ !is.na(pathway), bold = TRUE, part = "body") |>
    # flextable::bold(j = 2, i = ~ !is.na(type), bold = TRUE, part = "body") |>
    flextable::valign(j = 1:2, valign = "top", part = "all") |>
    flextable::add_header_row(
      values = c("", conditions, ""),
      colwidths = c(4, 3, 3, 1)
    ) |>
    flextable::compose(
      i = 1,
      j = 5:7,
      part = "header",
      value = flextable::as_paragraph(conditions[[1]], flextable::as_sup("a"))
    ) |>
    flextable::compose(
      i = 1,
      j = 8:10,
      part = "header",
      value = flextable::as_paragraph(conditions[[2]], flextable::as_sup("b"))
    ) |>
    flextable::set_header_labels(
      type = "Type",
      pathway = "Pathway",
      id = "ID",
      equation = "Reaction",
      flux_ctl = "Flux",
      lb_ctl = "LB",
      ub_ctl = "UB",
      flux_exp = "Flux",
      lb_exp = "LB",
      ub_exp = "UB",
      ratio = "Ratio"
    ) |>
    flextable::align(i = 1, part = "header", align = "center") |>
    flextable::align(i = 2, j = 5:10, part = "header", align = "center") |>
    flextable::merge_h(part = "header") |>
    flextable::bold(part = "header", bold = TRUE) |>
    flextable::colformat_double(
      digits = 2,
      big.mark = ""
    ) |>
    flextable::hline_top(part = "header", border = big_border) |>
    flextable::hline_bottom(part = "all", border = big_border) |>
    flextable::hline(i = ~ !is.na(type), border = small_border) |>
    flextable::hline(i = 1, j = c(5:7, 8:10), border = small_border, part = "header") |>
    flextable::add_footer_lines(c("a", "b")) |>
    flextable::compose(
      i = 1,
      part = "footer",
      value = flextable::as_paragraph(flextable::as_sup("a"), ssr_ctl)
    ) |>
    flextable::compose(
      i = 2,
      part = "footer",
      value = flextable::as_paragraph(flextable::as_sup("b"), ssr_exp)
    ) |>
    flextable::font(fontname = "Calibri", part = "all") |>
    flextable::fontsize(size = 8, part = "all") |>
    flextable::set_table_properties(layout = "autofit")
}

plot_exch_flux <- function(df, enzyme) {
  x <-
    df |>
    dplyr::filter(id == enzyme & type == "exch") |>
    dplyr::filter(treatment %in% c("21%", "0.5%"))

  eq <- unique(x$equation)

  lhs <- stringr::str_extract(eq, "\\w+")
  rhs <- stringr::str_extract(eq, "(?<=-> ).*")

  title <- stringr::str_c(rhs, " → ", lhs)

  ggplot2::ggplot(x) +
    ggplot2::aes(
      x = toupper(cell_type),
      y = flux,
      fill = treatment
    ) +
    ggplot2::geom_errorbar(
      ggplot2::aes(
        ymin = lb,
        ymax = ub
      ),
      position = ggplot2::position_dodge(width = 0.5),
      width = 0.2,
      size = 0.25
    ) +
    ggplot2::geom_point(
      pch = 21,
      position = ggplot2::position_dodge(width = 0.5),
      color = "white",
      size = 1.5,
      stroke = 0.2
    ) +
    ggplot2::scale_fill_manual(values = clrs, limits = force) +
    ggplot2::labs(
      x = NULL,
      y = "Exchange flux\n(fmol/cell/h)",
      fill = NULL,
      title = title
    ) +
    ggplot2::coord_cartesian(ylim = c(0, NA)) +
    theme_plots() +
    ggplot2::theme(
      panel.grid.major.y = ggplot2::element_line(color = "grey80", size = 0.1),
      legend.key.width = ggplot2::unit(0.5, "lines"),
      legend.key.height = ggplot2::unit(0.5, "lines"),
      legend.position = "bottom",
      legend.box.margin = ggplot2::margin(t = -10)
    )
}

arrange_rc <- function(p1, p2) {
  layout <- "ab"

  p1 + p2 +
    theme_patchwork(
      design = layout,
      widths = unit(c(1.25, 1), "in"),
      heights = unit(1, "in")
    ) &
    ggplot2::theme(
      legend.position = "bottom",
      legend.key.width = ggplot2::unit(0.5, "lines"),
      legend.key.height = ggplot2::unit(0.5, "lines"),
      legend.box.margin = ggplot2::margin(t = -10)
    )
}

arrange_graphs <- function(p1) {
  p1 +
    theme_patchwork(
      widths = ggplot2::unit(3.5, "in"),
      heights = ggplot2::unit(3.5, "in"),
      tags = NULL,
      guides = "collect"
    )
}

plot_lactate_mids <- function(df, cell, t = 72) {
  tracer_labels <- c(
    expression(paste("[1,2-"^13, "C"[2], "]-glucose")),
    expression(paste("[U-"^13, "C"[6], "]-glucose")),
    expression(paste("[U-"^13, "C"[5], "]-glutamine")),
    expression(paste("[U-"^13, "C"[3], "]-lactate"))
  )
  tracer_levels <- c("glc2", "glc6", "q5", "lac3")
  metab <- c(
    "LAC",
    "FBP",
    "PYR",
    "CIT",
    "AKG",
    "MAL"
  )

  x <-
    df |>
    dplyr::filter(
      cell_type == cell &
        time == t &
        metabolite %in% metab &
        tracer == "lac3"
    ) |>
    dplyr::filter(
      !(metabolite == "CIT" & isotope == "M6")
    ) |>
    dplyr::mutate(
      metabolite = factor(metabolite, levels = metab),
      metabolite = forcats::fct_recode(metabolite, "`3PG`" = "3PG"),
      tracer = factor(tracer, levels = tracer_levels, labels = tracer_labels),
      isotope = stringr::str_replace(isotope, "M", ""),
      group = dplyr::case_when(
        oxygen == "21%" & treatment == "None" ~ "21%",
        oxygen == "0.5%" & treatment == "None" ~ "0.5%",
        oxygen == "21%" & treatment == "DMSO" ~ "DMSO",
        oxygen == "21%" & treatment == "BAY" ~ "BAY"
      ),
      group = factor(group, levels = c("21%", "0.5%", "DMSO", "BAY"))
    )

  annot <-
    x |>
    dplyr::group_by(tracer, metabolite) |>
    tidyr::nest() |>
    dplyr::mutate(annot = purrr::map(data, annot_mids_main)) |>
    tidyr::unnest(c(annot))

  ggplot2::ggplot(x) +
    ggplot2::aes(
      x = isotope,
      y = mid
    ) +
    ggplot2::facet_grid(
      tracer ~ metabolite,
      labeller = ggplot2::label_parsed,
      # switch = "y",
      scales = "free_x",
      space = "free_x"
    ) +
    ggplot2::stat_summary(
      ggplot2::aes(fill = group),
      geom = "col",
      fun = "mean",
      position = ggplot2::position_dodge(width = 0.6),
      width = 0.6,
      show.legend = TRUE
    ) +
    ggplot2::geom_hline(
      yintercept = 0,
      color = "black",
      lwd = 0.1
    ) +
    ggplot2::stat_summary(
      ggplot2::aes(group = group),
      geom = "errorbar",
      fun.data = ggplot2::mean_se,
      position = ggplot2::position_dodge(width = 0.6),
      width = 0.2,
      size = 0.25,
      show.legend = FALSE
    ) +
    ggplot2::geom_text(
      data = annot,
      ggplot2::aes(
        x = -Inf,
        y = Inf,
        vjust = 1.5,
        hjust = -0.3,
        label = annot
      ),
      family = "Calibri",
      size = 8/ggplot2::.pt
    ) +
    ggplot2::labs(
      x = "Isotope",
      y = "Mole fraction",
      fill = NULL
    ) +
    ggplot2::scale_y_continuous(
      breaks = seq(0, 1, 0.25),
      limits = c(0, 1.1)
    ) +
    ggplot2::theme(
      strip.placement = "outside",
      legend.title = ggplot2::element_blank()
    ) +
    ggplot2::scale_fill_manual(values = clrs, limits = force) +
    theme_plots() +
    ggplot2::theme(
      panel.grid.major.y = ggplot2::element_line(color = "grey80", size = 0.1),
      legend.key.width = ggplot2::unit(0.5, "lines"),
      legend.key.height = ggplot2::unit(0.5, "lines"),
      legend.position = "bottom",
      legend.box.margin = ggplot2::margin(t = -10),
      axis.line = ggplot2::element_blank(),
      panel.border = ggplot2::element_rect(color = "black", fill = NA, size = 0.25)
    )
}

arrange_f7 <- function(p1, p2) {
  layout <- "ab"

  p1 + p2 +
    theme_patchwork(
      design = layout,
      widths = unit(c(1, 3), "in"),
      heights = unit(1, "in")
    ) &
    ggplot2::theme(
      legend.position = "bottom",
      legend.key.width = ggplot2::unit(0.5, "lines"),
      legend.key.height = ggplot2::unit(0.5, "lines"),
      legend.box.margin = ggplot2::margin(t = -10)
    )
}

plot_hyp_bay_fluxes <- function(df, annot, metab, ylab) {
  z <-
    df |>
    dplyr::ungroup() |>
    dplyr::filter(metabolite == metab) |>
    dplyr::mutate(grand_mean = mean(flux)) |>
    dplyr::group_by(date) |>
    dplyr::mutate(
      exp_mean = mean(flux),
      adj = grand_mean - exp_mean,
      flux_corr = flux + adj
    )

  ggplot2::ggplot(z) +
    ggplot2::aes(
      x = oxygen,
      y = flux_corr
    ) +
    ggplot2::stat_summary(
      ggplot2::aes(fill = treatment),
      geom = "col",
      fun = "mean",
      # width = 0.6,
      position = ggplot2::position_dodge2(),
      show.legend = TRUE,
      alpha = 0.5
    ) +
    ggbeeswarm::geom_beeswarm(
      ggplot2::aes(fill = treatment),
      method = "center",
      dodge.width = 0.9,
      pch = 21,
      size = 1,
      stroke = 0.2,
      cex = 4,
      color = "white",
      show.legend = FALSE
    ) +
    ggplot2::stat_summary(
      ggplot2::aes(group = treatment),
      geom = "errorbar",
      fun.data = ggplot2::mean_se,
      position = ggplot2::position_dodge(width = 0.9),
      width = 0.2,
      size = 0.25,
      show.legend = FALSE
    ) +
    ggplot2::geom_text(
      data = dplyr::filter(annot, metabolite == metab),
      ggplot2::aes(
        x = oxygen,
        y = y_pos,
        vjust = vjust,
        label = lab,
      ),
      family = "Calibri",
      size = 8/ggplot2::.pt
    ) +
    ggplot2::labs(
      x = "Oxygen",
      y = ylab,
      fill = NULL
    ) +
    ggplot2::scale_fill_manual(
      values = clrs,
      limits = force
    ) +
    ggplot2::scale_y_continuous(
      expand = ggplot2::expansion(mult = c(0.05, 0.1)),
      breaks = scales::pretty_breaks(n = 6)
    ) +
    ggplot2::guides(
      fill = ggplot2::guide_legend(override.aes = list(alpha = 1))
    ) +
    theme_plots() +
    ggplot2::theme(
      legend.key.width = ggplot2::unit(0.5, "lines"),
      legend.key.height = ggplot2::unit(0.5, "lines"),
      legend.position = "bottom",
      legend.box.margin = ggplot2::margin(t = -10)
    ) +
    NULL
}

arrange_f8 <- function(p1, p2, p3) {
  layout <- "abc"

  p1 + p2 + p3 +
    theme_patchwork(
      design = layout,
      widths = unit(c(1, 1), "in"),
      heights = unit(1, "in"),
      guides = "collect"
    ) &
    ggplot2::theme(
      legend.position = "bottom",
      legend.key.width = ggplot2::unit(0.5, "lines"),
      legend.key.height = ggplot2::unit(0.5, "lines"),
      legend.box.margin = ggplot2::margin(t = -10)
    )
}

plot_siphd <- function(df, annot, prot, ylab) {
  ggplot2::ggplot(df) +
    ggplot2::aes(
      x = oxygen,
      y = fold_change
    ) +
    ggplot2::stat_summary(
      ggplot2::aes(fill = treatment),
      geom = "col",
      fun = "mean",
      # width = 0.6,
      position = ggplot2::position_dodge2(),
      show.legend = TRUE,
      alpha = 0.5
    ) +
    ggbeeswarm::geom_beeswarm(
      ggplot2::aes(fill = treatment),
      method = "center",
      dodge.width = 0.9,
      pch = 21,
      size = 1,
      stroke = 0.2,
      cex = 4,
      color = "white",
      show.legend = FALSE
    ) +
    ggplot2::stat_summary(
      ggplot2::aes(group = treatment),
      geom = "errorbar",
      fun.data = ggplot2::mean_se,
      position = ggplot2::position_dodge(width = 0.9),
      width = 0.2,
      size = 0.25,
      show.legend = FALSE
    ) +
    ggplot2::geom_text(
      data = dplyr::filter(annot, !is.na(treatment)),
      ggplot2::aes(
        x = oxygen,
        y = y_pos,
        color = treatment,
        vjust = vjust,
        label = lab,
      ),
      position = ggplot2::position_dodge(width = 0.9),
      family = "Calibri",
      size = 8/ggplot2::.pt,
      show.legend = FALSE
    ) +
    ggplot2::geom_text(
      data = dplyr::filter(annot, is.na(treatment)),
      ggplot2::aes(
        x = oxygen,
        y = y_pos,
        vjust = vjust,
        label = lab,
      ),
      family = "Calibri",
      size = 8/ggplot2::.pt,
      show.legend = FALSE
    ) +
    ggplot2::labs(
      x = "Oxygen",
      y = ylab,
      fill = NULL
    ) +
    ggplot2::scale_fill_manual(
      values = clrs,
      limits = force,
      aesthetics = c("fill", "color")
    ) +
    ggplot2::scale_y_continuous(
      expand = ggplot2::expansion(mult = c(0.05, 0.1)),
      breaks = scales::pretty_breaks(n = 6)
    ) +
    ggplot2::guides(
      fill = ggplot2::guide_legend(override.aes = list(alpha = 1)),
      color = "none"
    ) +
    theme_plots() +
    ggplot2::theme(
      legend.key.width = ggplot2::unit(0.5, "lines"),
      legend.key.height = ggplot2::unit(0.5, "lines"),
      legend.position = "bottom",
      legend.box.margin = ggplot2::margin(t = -10)
    ) +
    NULL
}

plot_simyc <- function(df, annot, prot, ylab) {
  ggplot2::ggplot(df) +
    ggplot2::aes(
      x = treatment,
      y = fold_change
    ) +
    ggplot2::stat_summary(
      ggplot2::aes(fill = oxygen),
      geom = "col",
      fun = "mean",
      # width = 0.6,
      position = ggplot2::position_dodge2(),
      show.legend = TRUE,
      alpha = 0.5
    ) +
    ggbeeswarm::geom_beeswarm(
      ggplot2::aes(fill = oxygen),
      method = "center",
      dodge.width = 0.9,
      pch = 21,
      size = 1,
      stroke = 0.2,
      cex = 4,
      color = "white",
      show.legend = FALSE
    ) +
    ggplot2::stat_summary(
      ggplot2::aes(group = oxygen),
      geom = "errorbar",
      fun.data = ggplot2::mean_se,
      position = ggplot2::position_dodge(width = 0.9),
      width = 0.2,
      size = 0.25,
      show.legend = FALSE
    ) +
    ggplot2::geom_text(
      data = dplyr::filter(annot, !is.na(treatment)),
      ggplot2::aes(
        x = treatment,
        y = y_pos,
        vjust = vjust,
        label = lab,
      ),
      position = ggplot2::position_dodge(width = 0.9),
      family = "Calibri",
      size = 8/ggplot2::.pt,
      show.legend = FALSE
    ) +
    ggplot2::geom_text(
      data = dplyr::filter(annot, is.na(treatment)),
      ggplot2::aes(
        x = treatment,
        y = y_pos,
        color = oxygen,
        vjust = vjust,
        label = lab,
      ),
      family = "Calibri",
      size = 8/ggplot2::.pt,
      show.legend = FALSE
    ) +
    ggplot2::labs(
      x = "Oxygen",
      y = ylab,
      fill = NULL
    ) +
    ggplot2::scale_fill_manual(
      values = clrs,
      limits = force,
      aesthetics = c("fill", "color")
    ) +
    ggplot2::scale_y_continuous(
      expand = ggplot2::expansion(mult = c(0.05, 0.1)),
      breaks = scales::pretty_breaks(n = 6)
    ) +
    ggplot2::guides(
      fill = ggplot2::guide_legend(override.aes = list(alpha = 1)),
      color = "none"
    ) +
    theme_plots() +
    ggplot2::theme(
      legend.key.width = ggplot2::unit(0.5, "lines"),
      legend.key.height = ggplot2::unit(0.5, "lines"),
      legend.position = "bottom",
      legend.box.margin = ggplot2::margin(t = -10)
    ) +
    NULL
}

arrange_f8_s1 <- function(p1, p2, p3, p4, p5, p6) {
  layout <- "
  abc
  def
  "

  p1 + p2 + p3 + p4 + p5 + p6 +
    theme_patchwork(
      design = layout,
      widths = unit(c(1, 1), "in"),
      heights = unit(1, "in"),
      guides = "collect"
    ) &
    ggplot2::theme(
      legend.position = "bottom",
      legend.key.width = ggplot2::unit(0.5, "lines"),
      legend.key.height = ggplot2::unit(0.5, "lines"),
      legend.box.margin = ggplot2::margin(t = -10)
    )
}

arrange_f9 <- function(p1, p2, p3, p4) {
  layout <- "
  ab
  cd
  "

  p1 + p2 + p3 + p4 +
    theme_patchwork(
      design = layout,
      widths = unit(2.5, "in"),
      heights = unit(1.5, "in"),
      guides = "collect"
    ) &
    ggplot2::theme(
      legend.position = "bottom",
      legend.key.width = ggplot2::unit(0.5, "lines"),
      legend.key.height = ggplot2::unit(0.5, "lines"),
      legend.box.margin = ggplot2::margin(t = -10)
    )
}

arrange_f9_s1 <- function(p1, p2, p3, p4, p5) {
  layout <- "
  ad
  be
  ce
  "

  p1 + p2 + p3 + p4 + p5 +
    theme_patchwork(
      design = layout,
      widths = unit(2.5, "in"),
      heights = unit(c(1.5), "in")
    )
}

arrange_f9_s2 <- function(p1, p2, p3, p4, p5, p6) {
  layout <- "
  abc
  def
  "

  p1 + p2 + p3 + p4 + p5 + p6 +
    theme_patchwork(
      design = layout,
      widths = unit(1, "in"),
      heights = unit(1, "in"),
      guides = "collect"
    ) &
    ggplot2::theme(
      legend.position = "bottom",
      legend.key.width = ggplot2::unit(0.5, "lines"),
      legend.key.height = ggplot2::unit(0.5, "lines"),
      legend.box.margin = ggplot2::margin(t = -10)
    )
}

arrange_f10 <- function(p1, p2, p3, p4, p5, p6, p7) {
  layout <- "
  acdf
  bce#
  "

  p1 + p2 + p3 + p4 + p5 + p6 +
    theme_patchwork(
      design = layout,
      widths = unit(1.5, "in"),
      heights = unit(c(1.5), "in")
    )
}

arrange_f10_s1 <- function(p1, p2, p3, p4, p5, p6, p7, p8, p9) {
  layout <- "
  aeg
  beh
  #e#
  cfi
  df#
  "

  p1 + p2 + p3 + p4 + p5 + p6 + p7 + p8 + p9 +
    theme_patchwork(
      design = layout,
      widths = unit(2.5, "in"),
      heights = unit(c(1.5, 1.5, 1.1, 1.5, 1.5), "in")
    )
}

plot_cdkn1a <- function(dds) {
  idx <- which(SummarizedExperiment::rowData(dds)$hgnc_symbol == "CDKN1A")
  df <-
    dds[idx, ] |>
    SummarizedExperiment::assay() |>
    tibble::as_tibble(rownames = "entrez") |>
    tidyr::pivot_longer(
      -"entrez",
      names_to = "id",
      values_to = "count"
    ) |>
    dplyr::left_join(
      dplyr::select(
        tibble::as_tibble(SummarizedExperiment::colData(dds[idx, ])),
        "id":"group"
      ),
      by = "id",
      copy = TRUE
    ) |>
    dplyr::left_join(
      dplyr::select(
        tibble::as_tibble(
          SummarizedExperiment::rowData(dds[idx, ]),
          rownames = "entrez"
        ),
        "entrez",
        symbol = "hgnc_symbol"
      ),
      by = "entrez",
      copy = TRUE
    ) |>
    dplyr::mutate(
      oxygen = factor(.data$oxygen, levels = c("N", "H"), labels = c("21%", "0.5%"))
    )

  annot <-
    df |>
    dplyr::group_by(symbol) |>
    tidyr::nest() |>
    dplyr::mutate(
      m = purrr::map(
        data, \(x) lmerTest::lmer(
          count ~ oxygen * treatment + (1 | experiment),
          data = x
        )
      ),
      res = purrr::map(m, \(x) emmeans::emmeans(
        x,
        "pairwise" ~ oxygen * treatment,
        simple = "each",
        adjust = "mvt",
        combine = TRUE
      )[["contrasts"]]
      ),
      out = purrr::map(res, broom::tidy)
    ) |>
    tidyr::unnest(c(out)) |>
    dplyr::select(symbol, oxygen, treatment, adj.p.value) |>
    dplyr::mutate(
      oxygen = replace(oxygen, oxygen == ".", "0.5%"),
      oxygen = factor(oxygen, levels = c("21%", "0.5%")),
      treatment = factor(treatment, levels = c("DMSO", "BAY")),
      y_pos = Inf,
      vjust = 1,
      lab = dplyr::case_when(
        adj.p.value < 0.05 ~ "*",
        TRUE ~ NA_character_
      )
    )

  ggplot2::ggplot(df) +
    ggplot2::aes(
      x = .data$oxygen,
      y = .data$count/1000
    ) +
    ggplot2::stat_summary(
      ggplot2::aes(fill = .data$treatment),
      geom = "col",
      fun = "mean",
      position = ggplot2::position_dodge2(),
      show.legend = TRUE,
      alpha = 0.5
    ) +
    ggbeeswarm::geom_beeswarm(
      ggplot2::aes(fill = .data$treatment),
      method = "center",
      dodge.width = 0.9,
      pch = 21,
      size = 1,
      stroke = 0.2,
      cex = 4,
      color = "white",
      show.legend = FALSE
    ) +
    ggplot2::stat_summary(
      ggplot2::aes(group = treatment),
      geom = "errorbar",
      fun.data = ggplot2::mean_se,
      position = ggplot2::position_dodge(width = 0.9),
      width = 0.2,
      size = 0.25,
      show.legend = FALSE
    ) +
    ggplot2::geom_text(
      data = dplyr::filter(annot, !is.na(treatment)),
      ggplot2::aes(
        x = oxygen,
        y = y_pos,
        color = treatment,
        vjust = vjust,
        label = lab,
      ),
      position = ggplot2::position_dodge(width = 0.9),
      family = "Calibri",
      size = 8/ggplot2::.pt,
      show.legend = FALSE
    ) +
    ggplot2::geom_text(
      data = dplyr::filter(annot, is.na(treatment)),
      ggplot2::aes(
        x = oxygen,
        y = y_pos,
        vjust = vjust,
        label = lab,
      ),
      family = "Calibri",
      size = 8/ggplot2::.pt,
      show.legend = FALSE
    ) +
    ggplot2::scale_fill_manual(
      values = clrs,
      limits = force,
      aesthetics = c("fill", "color")
    ) +
    ggplot2::labs(
      x = "Oxygen",
      y = "CDKN1A<br>(count ×10<sup>3</sup>)",
      fill = NULL
    ) +
    ggplot2::guides(
      fill = ggplot2::guide_legend(override.aes = list(alpha = 1)),
      color = "none"
    ) +
    ggplot2::scale_y_continuous(
      expand = ggplot2::expansion(mult = c(0.05, 0.1)),
      breaks = scales::pretty_breaks(n = 6)
    ) +
    theme_plots() +
    ggplot2::theme(
      legend.key.width = ggplot2::unit(0.5, "lines"),
      legend.key.height = ggplot2::unit(0.5, "lines"),
      legend.position = "bottom",
      legend.box.margin = ggplot2::margin(t = -10),
      axis.title.y.left = ggtext::element_markdown(margin = ggplot2::margin(r = 3))
    ) +
    NULL
}

arrange_f10_s2 <- function(p1, p2) {
  layout <- "ab"

  p1 + p2 +
    theme_patchwork(
      design = layout,
      widths = unit(1, "in"),
      heights = unit(c(1), "in")
    )
}

plot_hyp_bay_densities <- function(df, annot, prot, ylab) {
  annot1 <-
    dplyr::filter(annot, !is.na(treatment)) |>
    dplyr::mutate(
      treatment = factor(treatment, levels = c("DMSO", "BAY")),
      oxygen = factor(oxygen, levels = c("21%", "0.5%"))
    )
  annot2 <-
    dplyr::filter(annot, is.na(treatment)) |>
    dplyr::mutate(
      treatment = factor(treatment, levels = c("DMSO", "BAY")),
      oxygen = factor(oxygen, levels = c("21%", "0.5%"))
    )

  df |>
    dplyr::mutate(
      treatment = factor(treatment, levels = c("DMSO", "BAY")),
      oxygen = factor(oxygen, levels = c("21%", "0.5%"))
    ) |>
    ggplot2::ggplot() +
    ggplot2::aes(
      x = oxygen,
      y = fold_change
    ) +
    ggplot2::stat_summary(
      ggplot2::aes(fill = treatment),
      geom = "col",
      fun = "mean",
      # width = 0.6,
      position = ggplot2::position_dodge2(),
      show.legend = TRUE,
      alpha = 0.5
    ) +
    ggbeeswarm::geom_beeswarm(
      ggplot2::aes(fill = treatment),
      method = "center",
      dodge.width = 0.9,
      pch = 21,
      size = 1,
      stroke = 0.2,
      cex = 4,
      color = "white",
      show.legend = FALSE
    ) +
    ggplot2::stat_summary(
      ggplot2::aes(group = treatment),
      geom = "errorbar",
      fun.data = ggplot2::mean_se,
      position = ggplot2::position_dodge(width = 0.9),
      width = 0.2,
      size = 0.25,
      show.legend = FALSE
    ) +
    ggplot2::geom_text(
      data = annot1,
      ggplot2::aes(
        color = treatment,
        y = y_pos,
        label = lab,
        vjust = vjust
      ),
      family = "Calibri",
      size = 8/ggplot2::.pt,
      position = ggplot2::position_dodge(width = 0.9),
      show.legend = FALSE
    ) +
    ggplot2::geom_text(
      data = annot2,
      ggplot2::aes(
        y = y_pos,
        label = lab,
        vjust = vjust
      ),
      family = "Calibri",
      size = 8/ggplot2::.pt,
      color = "black",
      show.legend = FALSE
    ) +
    ggplot2::labs(
      x = "Oxygen",
      y = ylab,
      fill = NULL,
      color = NULL
    ) +
    ggplot2::scale_fill_manual(
      values = clrs,
      limits = force,
      aesthetics = c("fill", "color")
    ) +
    ggplot2::scale_y_continuous(
      expand = ggplot2::expansion(mult = c(0.05, 0.1)),
      breaks = scales::pretty_breaks(n = 6)
    ) +
    ggplot2::guides(
      fill = ggplot2::guide_legend(override.aes = list(alpha = 1))
    ) +
    theme_plots() +
    ggplot2::theme(
      legend.key.width = ggplot2::unit(0.5, "lines"),
      legend.key.height = ggplot2::unit(0.5, "lines"),
      legend.position = "bottom",
      legend.box.margin = ggplot2::margin(t = -10)
    ) +
    NULL
}

plot_oemyc <- function(df, annot, prot, ylab) {
  ggplot2::ggplot(df) +
    ggplot2::aes(
      x = virus,
      y = fold_change
    ) +
    ggplot2::stat_summary(
      ggplot2::aes(fill = treatment),
      geom = "col",
      fun = "mean",
      # width = 0.6,
      position = ggplot2::position_dodge2(),
      show.legend = TRUE,
      alpha = 0.5
    ) +
    ggbeeswarm::geom_beeswarm(
      ggplot2::aes(fill = treatment),
      method = "center",
      dodge.width = 0.9,
      pch = 21,
      size = 1,
      stroke = 0.2,
      cex = 4,
      color = "white",
      show.legend = FALSE
    ) +
    ggplot2::stat_summary(
      ggplot2::aes(group = treatment),
      geom = "errorbar",
      fun.data = ggplot2::mean_se,
      position = ggplot2::position_dodge(width = 0.9),
      width = 0.2,
      size = 0.25,
      show.legend = FALSE
    ) +
    ggplot2::geom_text(
      data = dplyr::filter(annot, !is.na(treatment)),
      ggplot2::aes(
        x = virus,
        y = y_pos,
        color = treatment,
        vjust = vjust,
        label = lab,
      ),
      position = ggplot2::position_dodge(width = 0.9),
      family = "Calibri",
      size = 8/ggplot2::.pt,
      show.legend = FALSE
    ) +
    ggplot2::geom_text(
      data = dplyr::filter(annot, is.na(treatment)),
      ggplot2::aes(
        x = virus,
        y = y_pos,
        vjust = vjust,
        label = lab,
      ),
      family = "Calibri",
      size = 8/ggplot2::.pt,
      show.legend = FALSE
    ) +
    ggplot2::labs(
      x = "Treatment",
      y = ylab,
      fill = NULL
    ) +
    ggplot2::scale_fill_manual(
      values = clrs,
      limits = force,
      aesthetics = c("fill", "color")
    ) +
    ggplot2::scale_y_continuous(
      expand = ggplot2::expansion(mult = c(0.05, 0.1)),
      breaks = scales::pretty_breaks(n = 6)
    ) +
    ggplot2::guides(
      fill = ggplot2::guide_legend(override.aes = list(alpha = 1)),
      color = "none"
    ) +
    theme_plots() +
    ggplot2::theme(
      legend.key.width = ggplot2::unit(0.5, "lines"),
      legend.key.height = ggplot2::unit(0.5, "lines"),
      legend.position = "bottom",
      legend.box.margin = ggplot2::margin(t = -10)
    ) +
    NULL
}

arrange_f11 <- function(p1, p2, p3, p4, p5, p6, p7, p8, p9) {
  layout <- "
  abcc
  def#
  ghi#
  "

  p1 + p2 + p3 + p4 + p5 + p6 + p7 + p8 + p9 +
    theme_patchwork(
      design = layout,
      widths = unit(1, "in"),
      heights = unit(c(1), "in")
    )
}

arrange_f11_s1 <- function(p1, p2) {
  layout <- "ab"

  p1 + p2 +
    theme_patchwork(
      design = layout,
      widths = unit(c(1, 1), "in"),
      heights = unit(1, "in")
    ) &
    ggplot2::theme(
      legend.position = "bottom",
      legend.key.width = ggplot2::unit(0.5, "lines"),
      legend.key.height = ggplot2::unit(0.5, "lines"),
      legend.box.margin = ggplot2::margin(t = -10)
    )
}

arrange_f11_s2 <- function(p1, p2, p3, p4, p5, p6, p7, p8, p9, p10) {
  layout <- "
  abbbb
  cdddd
  "

  row1 <- p2 + p3 + p4 + p5 + patchwork::plot_layout(nrow = 1, guides = "collect")
  row2 <- p7 + p8 + p9 + p10 + patchwork::plot_layout(nrow = 1, guides = "collect")

  p1 + row1 + p6 + row2 +
    theme_patchwork(
      widths = unit(c(1, 5), "in"),
      heights = unit(c(1), "in")
    ) &
    ggplot2::theme(
      legend.position = "bottom",
      legend.key.width = ggplot2::unit(0.5, "lines"),
      legend.key.height = ggplot2::unit(0.5, "lines"),
      legend.box.margin = ggplot2::margin(t = -10)
    )
}

arrange_f12 <- function(p1, p2) {
  layout <- "ab"

  p1 + p2 +
    theme_patchwork(
      design = layout,
      widths = unit(1, "in"),
      heights = unit(c(1), "in"),
      guides = "collect"
    ) &
    ggplot2::theme(
      legend.position = "bottom",
      legend.key.width = ggplot2::unit(0.5, "lines"),
      legend.key.height = ggplot2::unit(0.5, "lines"),
      legend.box.margin = ggplot2::margin(t = -10)
    )
}

create_resources <- function(path) {
  readr::read_csv(path, show_col_types = FALSE) |>
    flextable::as_grouped_data(groups = c("category")) |>
    flextable::as_flextable(hide_grouplabel = TRUE) |>
    flextable::bold(j = 1, i = ~ !is.na(category), bold = TRUE, part = "body") |>
    flextable::bold(part = "header", bold = TRUE) |>
    flextable::colformat_double(
      i = ~ is.na(category),
      j = "Reagent/Resource",
      digits = 0,
      big.mark = ""
    ) |>
    flextable::compose(
      i = 11,
      j = 1,
      part = "body",
      value = flextable::as_paragraph("[1,2-", flextable::as_sup("13"), "C", flextable::as_sub("2"), "]-glucose")
    ) |>
    flextable::compose(
      i = 12,
      j = 1,
      part = "body",
      value = flextable::as_paragraph("[U-", flextable::as_sup("13"), "C", flextable::as_sub("6"), "]-glucose")
    ) |>
    flextable::compose(
      i = 13,
      j = 1,
      part = "body",
      value = flextable::as_paragraph("[U-", flextable::as_sup("13"), "C", flextable::as_sub("5"), "]-glutamine")
    ) |>
    flextable::compose(
      i = 14,
      j = 1,
      part = "body",
      value = flextable::as_paragraph("[U-", flextable::as_sup("13"), "C", flextable::as_sub("3"), "]-lactate")
    ) |>
    flextable::font(fontname = "Calibri", part = "all") |>
    flextable::fontsize(size = 9, part = "all") |>
    flextable::set_table_properties(layout = "autofit")
}
