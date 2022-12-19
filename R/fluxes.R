# fluxes.R

clean_flux_meta <- function(file_list) {
  file_list |>
    {\(x) rlang::set_names(x, stringr::str_extract(basename(x), "^.+(?=_)"))}() |>
    purrr::map_dfr(
      readr::read_csv,
      col_types = "dccccdc",
      .id = "experiment"
    ) |>
    tidyr::separate("experiment", c("cell_type", "experiment"), "_")
}

assemble_flux_data <- function(file_list) {
  nms <- sub("\\..*$", "", basename(file_list))
  sheets <- unique(unlist(purrr::map(file_list, readxl::excel_sheets)))
  data_list <- purrr::map(file_list, read_multi_excel)

  purrr::map(sheets, ~purrr::map(data_list, .x)) |>
    rlang::set_names(sheets) |>
    purrr::map(rlang::set_names, nms) |>
    purrr::map(dplyr::bind_rows, .id = "experiment")
}

clean_fluxes <- function(data_list, cells_per_dna){
  df <-
    data_list[c("dna", "glc", "lac", "pyr")] |>
    dplyr::bind_rows(.id = "metabolite") |>
    dplyr::filter(!is.na(.data$a)) |>
    dplyr::select(tidyselect::where(\(x) any(!is.na(x)))) |>
    clean_technical_replicates() |>
    tidyr::separate(.data$experiment, c("cell_type", "experiment", "batch", "date"), "_")

  dna <-
    df |>
    dplyr::filter(.data$metabolite == "dna") |>
    dplyr::mutate(volume = dplyr::if_else(.data$date <= "2018-05-25", 100, 200)) |>
    dplyr::left_join(cells_per_dna, by = c("cell_type", "volume")) |>
    dplyr::mutate(
      conc = .data$conc * .data$slope,
      metabolite = "dna",
      detector = "picogreen"
    ) |>
    dplyr::select(-c("volume", "slope"))

  others <-
    df |>
    dplyr::filter(.data$metabolite != "dna") |>
    dplyr::mutate(
      conc = dplyr::case_when(
        metabolite == "lac" & batch == "a" ~ .data$conc * 10.5,
        metabolite == "lac" & batch != "a" ~ .data$conc * 10,
        metabolite == "glc" ~ .data$conc * 555.074,
        metabolite == "pyr" ~ .data$conc * 20
      ),
      metabolite = dplyr::case_when(
        metabolite == "lac" ~ "lactate",
        metabolite == "glc" ~ "glucose",
        metabolite == "pyr" ~ "pyruvate"
      ),
      detector = "enzyme"
    )

  pyr <-
    data_list[["pyr"]] |>
    dplyr::filter(!is.na(.data$pyruvate)) |>
    tidyr::separate("experiment", c("cell_type", "experiment", "batch", "date"), "_") |>
    dplyr::mutate(istd = dplyr::coalesce(.data$KV, .data$`d8-valine`)) |>
    dplyr::mutate(
      detector = dplyr::case_when(
        !is.na(.data$KV) ~ "hplc",
        !is.na(.data$`d8-valine`) ~ "lcms",
        TRUE ~ "enzyme"
      ),
      istd = dplyr::case_when(
        experiment == "05" & batch == "a" & run == "a" & !is.na(.data$conc) ~ istd * 25,
        TRUE ~ .data$istd
      ),
      value = .data$pyruvate / .data$istd,
      metabolite = "pyruvate"
    ) |>
    dplyr::select("metabolite", "cell_type":"conc", "value", "detector")

  istds <- c("Norvaline", "Sarcosine")
  secondary_aa <- c("Hydroxyproline", "Proline")

  aa <-
    data_list[["aa"]] |>
    dplyr::mutate(
      dplyr::across(
        c(tidyselect::contains("1") & !tidyselect::contains(c(istds, secondary_aa))),
        ~ . /.data$`1 Norvaline`
      ),
      dplyr::across(
        c(tidyselect::contains("2") & !tidyselect::contains(c(istds, secondary_aa))),
        ~ . /.data$`2 Norvaline`
      ),
      dplyr::across(
        c(tidyselect::contains("1") & tidyselect::contains(secondary_aa)),
        ~ . /.data$`1 Sarcosine`
      ),
      dplyr::across(
        c(tidyselect::contains("2") & tidyselect::contains(secondary_aa)),
        ~ . /.data$`2 Sarcosine`
      ),
    ) |>
    tidyr::pivot_longer(
      tidyselect::matches("\\d .*"),
      names_to = "metabolite",
      values_to = "value"
    ) |>
    tidyr::separate(.data$metabolite, c("detector", "metabolite"), " ") |>
    tidyr::separate(.data$experiment, c("cell_type", "experiment", "batch", "date"), "_") |>
    dplyr::mutate(
      metabolite = tolower(.data$metabolite),
      detector = dplyr::if_else(.data$detector == 1, "mwd", "fld"),
      conc = replace(
        .data$conc,
        .data$experiment == "bay" &
          .data$date %in% c("2018-11-06", "2018-11-11") &
          .data$metabolite %in% c("asparagine", "glutamine", "tryptophan") &
          .data$conc == 225,
        22.5),
      conc = dplyr::case_when(
        .data$batch == "a" ~ 200 / 180 * .data$conc,
        TRUE ~ 200 / 190 * .data$conc
      )
    ) |>
    dplyr::filter(!(.data$cell_type == "pasmc" & .data$metabolite == "glutamine"))

  gln <-
    data_list[["gln"]] |>
    dplyr::mutate(
      value = .data$`1 Glutamine` / .data$`1 Norvaline`,
      detector = "mwd",
      conc = 20 * .data$conc,
      metabolite = "glutamine"
    ) |>
    dplyr::select("experiment":"conc", "detector", "value", "metabolite") |>
    tidyr::separate("experiment", c("cell_type", "experiment", "batch", "date"), "_")

  dplyr::bind_rows(dna, others, pyr, aa, gln) |>
    dplyr::relocate("detector", .after = "metabolite") |>
    dplyr::arrange(
      .data$metabolite,
      .data$detector,
      .data$cell_type,
      .data$experiment,
      .data$batch,
      .data$date,
      .data$run,
      .data$conc,
      .data$id
    ) |>
    dplyr::mutate(detector = tidyr::replace_na(.data$detector, "na")) |>
    dplyr::filter(!is.na(.data$value)) |>
    dplyr::filter(.data$metabolite %nin% c(tolower(istds), "hydroxyproline"))
}

format_glc6_raw <- function(path) {
  read_ms_excel(path) |>
    dplyr::filter(!stringr::str_detect(.data$id, "B")) |>
    tidyr::separate(
      .data$sample,
      c("type", "number"),
      sep = "(?<=[[:alpha:]])(?=[[:digit:]])",
      fill = "left",
      convert = TRUE
    ) |>
    dplyr::mutate(
      conc = ifelse(is.na(.data$type), .data$number, NA),
      conc = as.numeric(.data$conc),
      number = ifelse(!is.na(.data$conc), NA, .data$number),
      type = dplyr::case_when(
        .data$type == "QC" ~ "qc",
        is.na(.data$type) ~ "std",
        TRUE ~ tolower(.data$type)
      )
    ) |>
    new_tbl_se(
      a_data = "area",
      f_names = "metabolite",
      f_data = c("metabolite", "mz", "rt"),
      s_names = "id",
      s_data = c(
        "type",
        "number",
        "conc"
      )
    ) |>
    tbl_to_se()
}

clean_glc6_fluxes <- function(df) {
  df |>
    dplyr::select("metabolite", "id", "type", "number", "conc", "area") |>
    dplyr::group_by(.data$id) |>
    dplyr::mutate(area = .data$area / .data$area[.data$metabolite == "VAL D8"]) |>
    dplyr::filter(
      (.data$metabolite == "LAC M0" & .data$type == "std") |
        (.data$metabolite == "LAC M3" & .data$type %in% letters[1:3])
    ) |>
    tidyr::separate(.data$metabolite, c("metabolite", "isotope"), sep = " ") |>
    dplyr::mutate(metabolite = "lactate") |>
    dplyr::left_join(qbias_correction_factors, by = c("metabolite", "isotope" = "M")) |>
    dplyr::filter(.data$batch == "b" & .data$isotope %in% c("M0", "M3")) |>
    dplyr::mutate(value = .data$area * .data$cf) |>
    dplyr::select("metabolite", "type", id = "number", "conc", "value") |>
    dplyr::mutate(
      detector = "lcms",
      cell_type = "lf",
      experiment = "substrate",
      batch = "b",
      date = dplyr::case_when(
        .data$type == "a" ~ "2022-11-12",
        .data$type == "b" ~ "2022-11-17",
        .data$type == "c" ~ "2022-11-22"
      ),
      run = "a",
      conc = .data$conc * 1000
    ) |>
    dplyr::select(-"type")
}
