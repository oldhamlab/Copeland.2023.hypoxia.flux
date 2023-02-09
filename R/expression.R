# expression.R

read_data <- function(data_files) {
  data_files[stringr::str_detect(data_files, "\\.csv$")] |>
    {\(x) rlang::set_names(
      x,
      stringr::str_extract(x, "(lf|pasmc)_(02|05-bay|05|bay|)")
    )}() |>
    purrr::map_dfr(readr::read_csv, .id = "experiment") |>
    dplyr::mutate(
      oxygen = factor(
        .data$oxygen,
        levels = c("21%", "0.5%", "0.2%"),
        ordered = TRUE
      ),
      treatment = factor(
        .data$treatment,
        levels = c("None", "DMSO", "BAY"),
        ordered = TRUE
      )
    )
}

normalize_densities <- function(blot_raw) {
  blot_raw |>
    dplyr::filter(.data$time < 96) |>
    dplyr::filter(!(.data$experiment == "lf_05-bay" & .data$gel %in% c("b", "e", "f"))) |>
    tidyr::pivot_longer(
      "blot":tidyselect::last_col(),
      names_to = "protein",
      values_to = "value",
      values_drop_na = TRUE
    ) |>
    dplyr::group_by(.data$experiment, .data$gel, .data$protein) |>
    dplyr::mutate(norm = .data$value / mean(.data$value, na.rm = TRUE)) |>
    dplyr::group_by(dplyr::across("experiment":"time")) |>
    dplyr::mutate(density = .data$norm / .data$norm[.data$protein == "blot"]) |>
    dplyr::filter(.data$protein != "blot") |>
    dplyr::group_by(.data$experiment, .data$protein) |>
    dplyr::mutate(
      fold_change = .data$density /
        mean(.data$density[
          .data$oxygen == min(.data$oxygen) &
            .data$treatment %in% c("None", "DMSO") &
            .data$time == min(.data$time)
        ]),
      fold_change = replace(
        .data$fold_change,
        .data$experiment == "lf_bay" & .data$protein == "hif1a",
        sqrt(.data$fold_change)
      ),
      group = dplyr::case_when(
        .data$experiment %in% c("lf_02", "lf_05", "lf_bay", "pasmc_05") & .data$treatment == "None" & .data$oxygen == "21%" ~ "21%",
        .data$experiment %in% c("lf_02", "lf_05", "lf_bay", "pasmc_05") & .data$treatment == "DMSO" ~ "DMSO",
        .data$experiment %in% c("lf_02", "lf_05", "lf_bay", "pasmc_05") & .data$treatment == "BAY" ~ "BAY",
        .data$experiment %in% c("lf_02", "lf_05", "lf_bay", "pasmc_05") & .data$oxygen == "0.5%" ~ "0.5%",
        .data$experiment %in% c("lf_02", "lf_05", "lf_bay", "pasmc_05") & .data$oxygen == "0.2%" ~ "0.2%",
      ),
      group = factor(.data$group, levels = c("21%", "0.5%", "0.2%", "DMSO", "BAY"))
    ) |>
    dplyr::group_by(
      .data$experiment,
      .data$oxygen,
      .data$treatment,
      .data$group,
      .data$time,
      .data$protein
    ) |>
    # wmo::remove_nested_outliers(fold_change, remove = TRUE) |>
    dplyr::relocate("group", .after = "treatment")
}

normalize_qpcr <- function(raw_mrna) {
  raw_mrna |>
    dplyr::mutate(gene = tolower(.data$gene)) |>
    dplyr::group_by(dplyr::across(c("experiment":"gene"))) |>
    dplyr::summarize(ct = mean(.data$ct, na.rm = TRUE)) |>
    tidyr::pivot_wider(
      names_from = "gene",
      values_from = "ct"
    ) |>
    dplyr::ungroup() |>
    dplyr::mutate(
      dplyr::across(-c("experiment":"time", "actin"), \(x) x - .data$actin)
    ) |>
    dplyr::select(-"actin") |>
    tidyr::pivot_longer(
      -c("experiment":"time"),
      names_to = "protein",
      values_to = "dct"
    ) |>
    dplyr::filter(!is.na(.data$dct)) |>
    dplyr::group_by(.data$experiment, .data$protein) |>
    dplyr::mutate(
      ddct = .data$dct -
        mean(.data$dct[
          .data$oxygen == min(.data$oxygen) &
            .data$treatment %in% c("None", "DMSO") &
            .data$time == 0
        ]),
      fold_change = 2 ^ -.data$ddct,
      group = dplyr::case_when(
        .data$experiment %in% c("lf_02", "lf_05", "lf_bay", "pasmc_05") & .data$treatment == "None" & .data$oxygen == "21%" ~ "21%",
        .data$experiment %in% c("lf_02", "lf_05", "lf_bay", "pasmc_05") & .data$treatment == "DMSO" ~ "DMSO",
        .data$experiment %in% c("lf_02", "lf_05", "lf_bay", "pasmc_05") & .data$treatment == "BAY" ~ "BAY",
        .data$experiment %in% c("lf_02", "lf_05", "lf_bay", "pasmc_05") & .data$oxygen == "0.5%" ~ "0.5%",
        .data$experiment %in% c("lf_02", "lf_05", "lf_bay", "pasmc_05") & .data$oxygen == "0.2%" ~ "0.2%",
      ),
      group = factor(.data$group, levels = c("21%", "0.5%", "0.2%", "DMSO", "BAY"))
    ) |>
    dplyr::group_by(
      .data$experiment,
      .data$oxygen,
      .data$treatment,
      .data$group,
      .data$time,
      .data$protein
    ) |>
    wmo::remove_nested_outliers("fold_change", remove = TRUE) |>
    dplyr::relocate("group", .after = "treatment")
}
