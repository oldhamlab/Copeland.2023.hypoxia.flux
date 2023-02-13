# nadp.R

clean_nad <- function(nad_data) {
  nad_data |>
    dplyr::bind_rows(.id = "metabolite") |>
    dplyr::rename(rep = "replicate") |>
    clean_technical_replicates() |>
    tidyr::separate(.data$experiment, c(NA, "date"), "_") |>
    dplyr::filter(!(.data$metabolite == "nadp" & .data$conc > 0.2) | is.na(.data$conc)) |>
    dplyr::filter(!(.data$metabolite == "nad" & .data$conc > 1) | is.na(.data$conc))
}

finalize_nad <- function(nad_interp, cells_per_dna) {
  x <-
    cells_per_dna |>
    dplyr::filter(.data$cell_type == "lf" & .data$volume == 200) |>
    dplyr::pull(.data$slope)

  df <-
    nad_interp |>
    dplyr::filter(!(.data$date == "2023-01-11" & .data$oxygen == "21%" & .data$treatment == "BAY")) |>
    dplyr::group_by(
      .data$metabolite,
      .data$date,
      .data$oxygen,
      .data$treatment,
      .data$nucleotide
    ) |>
    dplyr::summarise(conc = mean(.data$conc)) |>
    dplyr::mutate(
      conc = dplyr::case_when(
        .data$metabolite == "dna" ~ .data$conc * x,
        .data$metabolite == "nad" ~ .data$conc * 720,
        .data$metabolite == "nadp" ~ .data$conc * 480,
      )
    ) |>
    dplyr::ungroup()

  counts <-
    df |>
    dplyr::filter(is.na(.data$nucleotide)) |>
    dplyr::select("date":"treatment", count = "conc")

  df |>
    dplyr::filter(!is.na(.data$nucleotide)) |>
    dplyr::mutate(
      nucleotide = dplyr::case_when(
        stringr::str_detect(.data$nucleotide, "H") ~ "red",
        TRUE ~ "ox"
      )
    ) |>
    tidyr::pivot_wider(names_from = "nucleotide", values_from = "conc") |>
    dplyr::left_join(counts, by = c("date", "oxygen", "treatment")) |>
    dplyr::mutate(
      oxygen = factor(.data$oxygen, levels = c("21%", "0.5%")),
      treatment = factor(
        .data$treatment,
        levels = c("none", "DMSO", "BAY"),
        labels = c("None", "DMSO", "BAY")
      ),
      ratio = .data$red / .data$ox,
      dplyr::across(c("red", "ox"), \(x) x / .data$count * 1000)
    ) |>
    dplyr::select(-"count") |>
    tidyr::pivot_longer(
      c("ox", "red", "ratio"),
      names_to = "measurement",
      values_to = "value"
    ) |>
    dplyr::mutate(
      measurement = dplyr::case_when(
        .data$metabolite == "nad" & .data$measurement == "ox" ~ "NAD",
        .data$metabolite == "nad" & .data$measurement == "red" ~ "NADH",
        .data$metabolite == "nad" & .data$measurement == "ratio" ~ "NADH/NAD",
        .data$metabolite == "nadp" & .data$measurement == "ox" ~ "NADP",
        .data$metabolite == "nadp" & .data$measurement == "red" ~ "NADPH",
        .data$metabolite == "nadp" & .data$measurement == "ratio" ~ "NADPH/NADP"
      )
    ) |>
    dplyr::select(-"metabolite") |>
    dplyr::arrange(.data$measurement, .data$oxygen, .data$treatment)
}

annot_nad <- function(nad_final) {
  nad_final |>
    dplyr::filter(.data$treatment != "None") |>
    dplyr::group_by(.data$measurement) |>
    tidyr::nest() |>
    dplyr::mutate(
      m = purrr::map(
        .data$data,
        \(x) lmerTest::lmer(value ~ treatment * oxygen + (1 | date), data = x)
      ),
      s = purrr::map(
        .data$m,
        \(x) emmeans::emmeans(x, ~ treatment * oxygen) |>
          graphics::pairs(simple = "each", combine = TRUE) |>
          broom::tidy()
      )
    ) |>
    tidyr::unnest(c("s")) |>
    dplyr::select("measurement", "oxygen", "treatment", pval = "adj.p.value") |>
    dplyr::mutate(
      lab = annot_p(.data$pval),
      y = Inf,
      vjust = 1
    )
}

plot_nad <- function(df, annot, metab, ylab) {
  annot <- dplyr::filter(annot, .data$measurement == metab)
  annot1 <-
    dplyr::filter(annot, .data$oxygen != ".") |>
    dplyr::mutate(
      treatment = "BAY",
      treatment = factor(.data$treatment, levels = c("DMSO", "BAY")),
      oxygen = factor(.data$oxygen, levels = c("21%", "0.5%"))
    )
  annot2 <-
    dplyr::filter(annot, .data$treatment != ".") |>
    dplyr::mutate(
      oxygen = "21%",
      treatment = factor(.data$treatment, levels = c("DMSO", "BAY")),
      oxygen = factor(.data$oxygen, levels = c("21%", "0.5%"))
    )

  z <-
    df |>
    dplyr::ungroup() |>
    dplyr::filter(.data$treatment != "None" & .data$measurement == metab) |>
    dplyr::mutate(grand_mean = mean(.data$value)) |>
    dplyr::group_by(date) |>
    dplyr::mutate(
      exp_mean = mean(.data$value),
      adj = .data$grand_mean - .data$exp_mean,
      value_corr = .data$value + .data$adj
    )

  ggplot2::ggplot(z) +
    ggplot2::aes(
      x = .data$treatment,
      y = .data$value_corr,
      fill = .data$oxygen
    ) +
    ggplot2::stat_summary(
      # ggplot2::aes(fill = treatment),
      geom = "col",
      fun = "mean",
      # width = 0.6,
      position = ggplot2::position_dodge2(),
      show.legend = TRUE,
      alpha = 0.5
    ) +
    ggbeeswarm::geom_beeswarm(
      # ggplot2::aes(fill = treatment),
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
      size = 0.25,
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
    ggplot2::scale_color_manual(
      values = clrs,
      limits = force,
      aesthetics = c("fill", "color")
    ) +
    ggplot2::labs(
      x = "Treatment",
      y = ylab,
      fill = NULL,
      color = NULL
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
