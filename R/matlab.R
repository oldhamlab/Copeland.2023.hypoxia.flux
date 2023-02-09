# matlab.R

format_reactions <- function(reactions_file) {
  readr::read_csv(
    reactions_file,
    show_col_types = FALSE
  ) |>
    dplyr::filter(!is.na(.data$name)) |>
    dplyr::mutate(
      pathway = factor(
        .data$pathway,
        levels = c(
          "transport",
          "glycolysis",
          "pentose phosphate pathway",
          "anaplerosis",
          "tricarboxylic acid cycle",
          "amino acid metabolism",
          "biomass",
          "mixing",
          "dilution"
        )
      )
    )
}

format_fluxes <- function(growth_rates, fluxes) {
  growth <-
    growth_rates |>
    dplyr::filter(.data$experiment %in% c("05", "02", "bay")) |>
    dplyr::mutate(metabolite = "biomass") |>
    dplyr::select(-"X0") |>
    dplyr::rename(flux = "mu")

  fluxes_final <-
    fluxes |>
    dplyr::filter(.data$experiment %in% c("05", "02", "bay")) |>
    dplyr::filter(!(.data$experiment == "02" & .data$metabolite == "pyruvate")) |>
    dplyr::filter(!(.data$experiment == "02" & .data$oxygen == "0.2%")) |>
    dplyr::select(-"abbreviation") |>
    dplyr::bind_rows(growth) |>
    dplyr::group_by(.data$metabolite, .data$cell_type, .data$oxygen, .data$treatment) |>
    wmo::remove_nested_outliers("flux", remove = TRUE)

  flux_names <-
    tibble::tribble(
      ~ metabolite, ~ name,
      "biomass", "BIOMASS",
      "alanine", "ALAR",
      "glucose", "GLUT",
      "glutamine", "GLNR",
      "glutamate", "GLUR",
      "lactate", "MCT",
      "pyruvate", "PYRR",
      "serine", "SERR",
      "cystine", "CYSR",
      "glycine", "GLYR",
      "aspartate", "ASPR"
    )

  model_fluxes <-
    c(
      "biomass",
      "alanine",
      "glucose",
      "glutamine",
      "glutamate",
      "lactate",
      "pyruvate",
      "serine",
      "cystine",
      "glycine",
      "aspartate"
    )

  fluxes_final |>
    dplyr::filter(.data$metabolite %in% model_fluxes) |>
    dplyr::group_by(.data$metabolite, .data$cell_type, .data$oxygen, .data$treatment) |>
    dplyr::summarise(
      mean = mean(.data$flux, na.rm = TRUE),
      se = stats::sd(.data$flux, na.rm = TRUE) / sqrt(dplyr::n())
    ) |>
    dplyr::filter(!is.nan(.data$mean)) |>
    dplyr::left_join(flux_names, by = "metabolite") |>
    dplyr::group_by(.data$cell_type) |>
    tidyr::nest()
}

format_mids <- function(mids) {
  model_metabolites <-
    c(
      "pyruvate",
      "lactate",
      "alanine",
      "2OG",
      "malate",
      "aspartate",
      "glutamate",
      "glutamine",
      "citrate",
      "serine",
      "FBP",
      "3PG"
    )

  mids_filtered <-
    mids |>
    dplyr::filter(.data$metabolite %in% model_metabolites) |>
    dplyr::filter(.data$time != 96)

  mids_lf_new <-
    mids_filtered |>
    dplyr::filter(
      .data$cell_type == "lf" &
        .data$date %in% c(
          "2018-11-15", "2018-11-20", "2018-11-25", "2018-12-16",
          "2018-11-11", "2018-11-16",
          "2019-05-06", "2019-05-10", "2019-05-14"
        )
    ) |>
    dplyr::filter(
      .data$method == "sim" |
        (.data$method == "fs" & .data$metabolite %in% c("FBP", "3PG"))
    )

  mids_pasmc_05 <-
    mids_filtered |>
    dplyr::filter(
      .data$cell_type == "pasmc" &
        .data$oxygen %in% c("21%", "0.5%")
    ) |>
    dplyr::filter(
      .data$method == "sim" |
        (.data$method == "fs" & .data$metabolite %in% c("FBP", "3PG"))
    )

  dplyr::bind_rows(mids_lf_new, mids_pasmc_05) |>
    dplyr::group_by(
      .data$method,
      .data$cell_type,
      .data$tracer,
      .data$oxygen,
      .data$treatment,
      .data$metabolite,
      .data$time,
      .data$isotope
    ) |>
    wmo::remove_nested_outliers("mid", remove = TRUE) |>
    dplyr::mutate(
      metabolite = replace(.data$metabolite, .data$metabolite == "glutamine", "GLN"),
      metabolite = replace(.data$metabolite, .data$metabolite == "2OG", "AKG"),
      metabolite = stringr::str_sub(.data$metabolite, 1, 3),
      metabolite = toupper(.data$metabolite)
    ) |>
    dplyr::select(-"method")
}

summarize_mids <- function(df) {
  df |>
    dplyr::summarise(
      mean = mean(.data$mid, na.rm = TRUE),
      se = stats::sd(.data$mid, na.rm = TRUE) / sqrt(dplyr::n())
    ) |>
    dplyr::ungroup() |>
    dplyr::arrange(
      .data$metabolite,
      .data$tracer,
      .data$oxygen,
      .data$treatment,
      .data$time,
      .data$isotope
    ) |>
    dplyr::group_by(.data$cell_type) |>
    tidyr::nest()
}
