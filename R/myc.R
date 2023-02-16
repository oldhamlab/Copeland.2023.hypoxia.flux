# myc.R

combine_fluxes <- function(growth_rates, fluxes, exp) {
  growth_rates |>
    dplyr::filter(.data$experiment == exp) |>
    dplyr::rename(flux = "mu") |>
    dplyr::select(-"X0") |>
    dplyr::mutate(metabolite = "growth") |>
    dplyr::bind_rows(dplyr::filter(fluxes, .data$experiment == exp)) |>
    dplyr::filter(.data$treatment %nin% c("siHIF1A", "siHIF2A")) |>
    dplyr::mutate(
      flux = ifelse(
        .data$experiment == "bay-myc" &
          .data$metabolite == "lactate" &
          .data$date %in% c("2021-09-21", "2021-11-01"),
        .data$flux / 3,
        .data$flux
      )
    ) |>
    dplyr::group_by(.data$oxygen, .data$treatment, .data$virus) |>
    wmo::remove_nested_outliers("flux", TRUE)
}

annot_myc_fluxes <- function(df, intervention) {
  if (intervention == "treatment") {
    lvls <- c("siCTL", "siMYC")
    x <- "oxygen"
  } else if (intervention == "virus") {
    lvls <- c("YFP", "MYC")
    x <- "virus"
  }
  fo1 <- stats::as.formula(paste("flux ~ ", x, " * treatment + (1 | date)"))
  fo2 <- stats::as.formula(paste("pairwise ~ ", x, "* treatment"))

  int <- rlang::sym(intervention)

  df |>
    dplyr::group_by(.data$metabolite) |>
    tidyr::nest() |>
    dplyr::mutate(
      m = purrr::map(.data$data, \(x) lmerTest::lmer(fo1, data = x)),
      res = purrr::map(
        .data$m,
        \(x) emmeans::emmeans(
          x,
          fo2,
          simple = "each",
          adjust = "mvt",
          combine = TRUE
        )[["contrasts"]]
      ),
      out = purrr::map(.data$res, broom::tidy)
    ) |>
    tidyr::unnest(c("out")) |>
    dplyr::filter(!!int != ".") |>
    dplyr::select("metabolite", {{ intervention }}, "adj.p.value") |>
    dplyr::mutate(
      group = factor(!!int, levels = lvls),
      y_pos = Inf,
      vjust = 1.5,
      lab = annot_p(.data$adj.p.value)
    ) |>
    dplyr::select(-!!int) |>
    dplyr::rename(!!int := "group")
}

plot_myc <- function(df, annot, metab, ylab, x, fill) {
  z <-
    df |>
    dplyr::ungroup() |>
    dplyr::filter(.data$metabolite == metab) |>
    dplyr::mutate(grand_mean = mean(.data$flux)) |>
    dplyr::group_by(.data$date) |>
    dplyr::mutate(
      exp_mean = mean(.data$flux),
      adj = .data$grand_mean - .data$exp_mean,
      flux_corr = .data$flux + .data$adj
    )

  ggplot2::ggplot(z) +
    ggplot2::aes(
      x = {{x}},
      y = .data$flux_corr
    ) +
    ggplot2::stat_summary(
      ggplot2::aes(
        fill = {{fill}}
      ),
      geom = "col",
      fun = "mean",
      position = ggplot2::position_dodge2(),
      alpha = 0.5,
      show.legend = TRUE
    ) +
    ggbeeswarm::geom_beeswarm(
      ggplot2::aes(fill = {{fill}}),
      dodge.width = 0.9,
      pch = 21,
      size = 1,
      stroke = 0.2,
      cex = 4,
      color = "white",
      show.legend = FALSE
    ) +
    ggplot2::stat_summary(
      ggplot2::aes(
        group = {{fill}}
      ),
      geom = "errorbar",
      fun.data = ggplot2::mean_se,
      position = ggplot2::position_dodge(width = 0.9),
      width = 0.2,
      linewidth = 0.25,
      show.legend = FALSE,
      color = "black"
    ) +
    ggplot2::geom_text(
      data = dplyr::filter(annot, .data$metabolite == metab),
      ggplot2::aes(
        x = {{x}},
        y = .data$y_pos,
        vjust = .data$vjust,
        label = .data$lab,
      ),
      position = ggplot2::position_dodge(width = 0.9),
      family = "Calibri",
      color = "black",
      size = 8/ggplot2::.pt,
      show.legend = FALSE
    ) +
    ggplot2::labs(
      x = NULL,
      y = ylab,
      color = NULL,
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
