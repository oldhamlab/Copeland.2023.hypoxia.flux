# model.R

clean_model_fluxes <- function(map_flux_files, model_reactions) {
  pathways <-
    model_reactions |>
    dplyr::select(-"equation") |>
    dplyr::add_row(
      name = "BIOMASS",
      pathway = factor("biomass"),
      index = 0
    )

  map_flux_files[stringr::str_detect(map_flux_files, "(lf|pasmc)_.*_model\\.csv")] |>
    {\(x) rlang::set_names(x, stringr::str_extract(basename(x), pattern = ".*(?=_model\\.csv)"))}() |>
    purrr::map_dfr(
      readr::read_csv,
      show_col_types = FALSE,
      .id = "experiment"
    ) |>
    tidyr::separate(
      "experiment",
      c("cell_type", "treatment"),
      sep = "_"
    ) |>
    dplyr::mutate(
      treatment = dplyr::case_when(
        .data$treatment == "21" ~ "21%",
        .data$treatment == "05" ~ "0.5%",
        .data$treatment == "dmso" ~ "DMSO",
        .data$treatment == "bay" ~ "BAY"
      ),
      treatment = factor(.data$treatment, levels = c("21%", "0.5%", "DMSO", "BAY"))
    ) |>

    # identify net and exchange fluxes
    tidyr::separate("id", c("id", "type"), sep = " ", fill = "right") |>
    dplyr::mutate(type = replace(.data$type, is.na(.data$type), "net")) |>
    dplyr::select(-"se") |>

    # add pathway info
    dplyr::left_join(pathways, by = c("id" = "name")) |>
    dplyr::select(
      "cell_type",
      "treatment",
      "pathway",
      "index",
      "id",
      "type",
      "equation",
      flux = "value",
      "lb",
      "ub"
    ) |>
    dplyr::group_by(.data$cell_type, .data$treatment) |>
    dplyr::arrange(.data$treatment, .data$pathway)
}

calculate_flux_differences <- function(map_fluxes, cell, control, experiment) {
  map_fluxes |>
    dplyr::filter(.data$cell_type == cell) |>
    dplyr::filter(.data$treatment %in% c(control, experiment)) |>
    dplyr::mutate(
      treatment = dplyr::case_when(
        .data$treatment == control ~ "ctl",
        .data$treatment == experiment ~ "exp"
      )
    ) |>
    tidyr::pivot_wider(
      names_from = "treatment",
      values_from = c("flux", "lb", "ub")
    ) |>
    dplyr::mutate(
      ratio = dplyr::if_else(
        pmax(.data$lb_ctl, .data$lb_exp) - pmin(.data$ub_ctl, .data$ub_exp) > 0,
        .data$flux_exp / .data$flux_ctl,
        NA_real_
      ),
      ctl = control,
      exp = experiment
    ) |>
    dplyr::select(
      "cell_type":"equation",
      "ctl",
      "exp",
      tidyselect::contains("ctl"),
      tidyselect::contains("exp"),
      "ratio"
    )
}

assemble_flux_differences <- function(map_fluxes) {
  glc_norm <- normalize_fluxes(map_fluxes, "GLUT")
  growth_norm <- normalize_fluxes(map_fluxes, "BIOMASS")

  tibble::tibble(
    map_fluxes = list(
      map_fluxes,
      map_fluxes,
      map_fluxes,
      glc_norm,
      glc_norm,
      glc_norm,
      growth_norm,
      growth_norm,
      growth_norm
    ),
    cell = rep(c("lf", "lf", "pasmc"), 3),
    control = rep(c("21%", "DMSO", "21%"), 3),
    experiment = rep(c("0.5%", "BAY", "0.5%"), 3)
  ) |>
    purrr::pmap_dfr(calculate_flux_differences, .id = "normalization") |>
    dplyr::mutate(
      normalization = dplyr::case_when(
        .data$normalization %in% 1:3 ~ "none",
        .data$normalization %in% 4:6 ~ "glucose",
        .data$normalization %in% 7:9 ~ "growth"
      )
    )
}

normalize_fluxes <- function(map_fluxes, reference_flux) {
  map_fluxes |>
    dplyr::group_by(.data$cell_type, .data$treatment) |>
    dplyr::mutate(
      lb = .data$lb / .data$flux[.data$id == reference_flux],
      ub = .data$ub / .data$flux[.data$id == reference_flux],
      flux = .data$flux / .data$flux[.data$id == reference_flux]
    )
}
