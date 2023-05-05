# cell-dna.R

clean_dna_per_cell <- function(path) {
  path |>
    read_multi_excel() |>
    purrr::map(clean_technical_replicates) |>
    dplyr::bind_rows(.id = "id") |>
    tidyr::separate(.data$id, c("cell_type", "volume"), sep = "_", convert = TRUE) |>
    dplyr::mutate(cells = 1000 * .data$cells)
}

calculate_cells_per_dna <- function(df) {
  df |>
    dplyr::filter(!(.data$cell_type == "lf" & .data$volume == "100" & .data$cells == 300000)) |>
    dplyr::filter(!(.data$cell_type == "pasmc" & .data$volume == "200" & .data$cells == 400000)) |>
    dplyr::group_by(.data$cell_type, .data$volume) |>
    tidyr::nest() |>
    dplyr::mutate(
      model = purrr::map(.data$data, \(x) stats::lm(conc ~ 0 + cells, data = x, na.action = modelr::na.warn)),
      glance = purrr::map(.data$model, broom::tidy)
    ) |>
    tidyr::unnest(c("glance")) |>
    dplyr::select("cell_type", "volume", slope = "estimate") |>
    dplyr::mutate(slope = 1 / .data$slope) |>
    dplyr::ungroup()
}

plot_cells_per_dna <- function(dna_per_cell_clean, cell = c("lf", "pasmc")) {
  dna_per_cell_clean |>
    dplyr::filter(cell_type == cell) |>
    dplyr::filter(.data$volume == 200 & cells < 400000) |>
    ggplot2::ggplot() +
    ggplot2::facet_wrap(~cell_type, labeller = ggplot2::as_labeller(toupper)) +
    ggplot2::aes(
      x = cells,
      y = conc
    ) +
    ggplot2::geom_smooth(
      formula = y ~ 0 + x,
      method = "lm",
      color = clrs[[2]],
      size = 0.5,
      se = FALSE
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
      width = 7500,
      size = 0.25,
      show.legend = FALSE
    ) +
    ggplot2::stat_summary(
      geom = "point",
      fun = "mean",
      pch = 21,
      color = "white",
      fill = "black",
      size = 1.5,
      show.legend = FALSE,
      stroke = 0.2
    ) +
    ggplot2::labs(
      x = "Cell count",
      y = "DNA (ng)"
    ) +
    ggplot2::scale_x_continuous(
      labels = scales::label_number(scale_cut = scales::cut_short_scale())
    ) +
    theme_plots() +
    ggplot2::coord_cartesian(xlim = c(0, NA), clip = "off") +
    NULL
}
