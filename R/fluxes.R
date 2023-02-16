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

preclean_fluxes <- function(data_list, metabolite, name = metabolite) {
  data_list[[metabolite]] |>
    dplyr::mutate(metabolite = name) |>
    dplyr::filter(!is.na(.data$a)) |>
    dplyr::select(tidyselect::where(\(x) any(!is.na(x))) | letters[1:3]) |>
    clean_technical_replicates() |>
    tidyr::separate(.data$experiment, c("cell_type", "experiment", "batch", "date"), "_")
}

clean_dna_fluxes <- function(data_list, cells_per_dna) {
  preclean_fluxes(data_list, "dna") |>
    dplyr::mutate(volume = dplyr::if_else(.data$date <= "2018-05-25", 100, 200)) |>
    dplyr::left_join(cells_per_dna, by = c("cell_type", "volume")) |>
    dplyr::mutate(
      conc = .data$conc * .data$slope,
      metabolite = "dna",
      detector = "picogreen"
    ) |>
    dplyr::select(-c("volume", "slope"))
}

clean_glc_fluxes <- function(data_list) {
  preclean_fluxes(data_list, "glc", "glucose") |>
    dplyr::mutate(
      conc = .data$conc * 555.074,
      detector = "enzyme"
    )
}

clean_lac_fluxes <- function(data_list) {
  preclean_fluxes(data_list, "lac", "lactate") |>
    dplyr::mutate(
      conc = dplyr::case_when(
        batch == "a" ~ .data$conc * 10.5,
        batch != "a" ~ .data$conc * 10,
      ),
      detector = "enzyme"
    )
}

clean_pyr_fluxes <- function(data_list) {
  pyr_1 <-
    preclean_fluxes(data_list, "pyr", "pyruvate") |>
    dplyr::mutate(
      conc = .data$conc * 20,
      detector = "enzyme"
    )

  pyr_2 <-
    data_list[["pyr"]] |>
    dplyr::filter(!is.na(.data$pyruvate)) |>
    tidyr::separate(
      "experiment",
      c("cell_type", "experiment", "batch", "date"),
      sep = "_"
    ) |>
    dplyr::mutate(istd = dplyr::coalesce(.data$KV, .data$`d8-valine`)) |>
    dplyr::mutate(
      detector = dplyr::case_when(
        !is.na(.data$KV) ~ "hplc",
        !is.na(.data$`d8-valine`) ~ "lcms"
      ),
      istd = dplyr::case_when(
        experiment == "05" & batch == "a" & run == "a" & !is.na(.data$conc) ~ istd * 25,
        TRUE ~ .data$istd
      ),
      value = .data$pyruvate / .data$istd,
      metabolite = "pyruvate"
    ) |>
    dplyr::select("metabolite", "cell_type":"conc", "value", "detector")

  dplyr::bind_rows(pyr_1, pyr_2)
}

clean_gln_fluxes <- function(data_list) {
  data_list[["gln"]] |>
    dplyr::mutate(
      value = .data$`1 Glutamine` / .data$`1 Norvaline`,
      detector = "mwd",
      conc = 20 * .data$conc,
      metabolite = "glutamine"
    ) |>
    dplyr::select("experiment":"conc", "detector", "value", "metabolite") |>
    tidyr::separate("experiment", c("cell_type", "experiment", "batch", "date"), "_")
}

clean_aa_fluxes <- function(data_list) {
  istds <- c("Norvaline", "Sarcosine")
  secondary_aa <- c("Hydroxyproline", "Proline")

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

clean_glc6_fluxes <- function(df, cf) {
  new_id <-
    outer(10:12, 12 * 0:7, `+`) |>
    as.vector()

  x <-
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
    dplyr::left_join(cf, by = c("metabolite", "isotope" = "M"), multiple = "all") |>
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
      conc = .data$conc * 1000,
      id = replace(.data$id, !is.na(.data$id), new_id)
    ) |>
    dplyr::select(-"type")

  samples <-
    x |>
    dplyr::filter(is.na(.data$conc))

  std <-
    x |>
    dplyr::filter(is.na(.data$id))

  dplyr::bind_rows(std, std, std) |>
    dplyr::mutate(date = rep(c("2022-11-12", "2022-11-17", "2022-11-22"), each = 12)) |>
    dplyr::bind_rows(samples) |>
    dplyr::arrange(.data$date)
}

clean_fluxes <- function(...) {
  dplyr::bind_rows(...) |>
    dplyr::relocate(c("metabolite", "detector"), .before = 1) |>
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
    dplyr::filter(.data$metabolite %nin% c("norvaline", "sarcosine", "hydroxyproline"))
}

clean_flux_std <- function(df) {
  outliers <-
    tibble::tribble(
      ~metabolite, ~experiment, ~date, ~batch, ~run, ~detector, ~conc,
      "lactate", "bay", "2018-06-02", "a", "a", "enzyme", 10500,
      "pyruvate", "bay", "2018-11-02", "b", "a", "enzyme", 1500,
      "glycine", "05", "2017-11-06", "a", "b", "fld", 10
    )

  df |>
    dplyr::filter(!(.data$detector == "fld" & .data$conc > 900)) |>
    dplyr::filter(!(.data$experiment == "substrate" & .data$conc == 10000)) |>
    dplyr::anti_join(
      outliers,
      by = c(
        "metabolite",
        "experiment",
        "date",
        "batch",
        "run",
        "detector",
        "conc"
      )
    ) |>
    make_std_curves()
}

fill_missing_fluxes <- function(df, meta) {

  df_meta <-
    dplyr::left_join(df, meta, by = c("cell_type", "experiment", "id"))

  # https://github.com/tidyverse/tidyr/issues/971
  metabolite <- type <- detector <- time <- well <- NULL

  missing_data <-
    df_meta |>
    dplyr::group_by(.data$cell_type, .data$experiment, .data$batch) |>
    tidyr::complete(
      .data$date,
      .data$treatment,
      .data$oxygen,
      .data$virus,
      tidyr::nesting(
        metabolite,
        type,
        detector,
        time,
        well
      )
    ) |>
    dplyr::filter(is.na(.data$conc))

  empty_05_t0 <-
    missing_data |>
    dplyr::filter(
      .data$experiment == "05" &
        .data$oxygen == "0.5%" &
        .data$time == 0 &
        .data$type == "empty"
    ) |>
    dplyr::select(-c("run", "id", "conc")) |>
    dplyr::mutate(oxygen = forcats::fct_recode(.data$oxygen, "21%" = "0.5%")) |>
    dplyr::left_join(
      df_meta,
      by = c(
        "cell_type",
        "experiment",
        "batch",
        "date",
        "virus",
        "treatment",
        "oxygen",
        "metabolite",
        "type",
        "detector",
        "time",
        "well"
      )
    ) |>
    dplyr::mutate(oxygen = forcats::fct_recode(.data$oxygen, "0.5%" = "21%"))

  empty_simyc_t0 <-
    missing_data |>
    dplyr::filter(.data$experiment == "05-simyc" & .data$time == 0) |>
    dplyr::select(-c("run", "id", "conc")) |>
    dplyr::left_join(
      df_meta,
      by = c(
        "cell_type",
        "experiment",
        "batch",
        "date",
        "oxygen",
        "virus",
        "metabolite",
        "type",
        "detector",
        "time",
        "well"
      ),
      multiple = "all"
    ) |>
    dplyr::select(-"treatment.y") |>
    dplyr::rename(treatment = "treatment.x")

  dplyr::bind_rows(df_meta, empty_05_t0, empty_simyc_t0)
}

filter_assays <- function(df) {
  mwd <-
    c(
      "cystine",
      "glutamine",
      "isoleucine",
      "leucine",
      "lysine",
      "valine"
    )

  df |>
    dplyr::mutate(
      keep = dplyr::case_when(
        .data$metabolite %in% mwd & .data$detector == "mwd" ~ TRUE,
        .data$metabolite %nin% mwd & .data$detector == "fld" ~ TRUE,
        .data$detector %in% c("enzyme", "picogreen") ~ TRUE,
        .data$metabolite == "pyruvate" & .data$experiment == "02" ~ TRUE,
        .data$metabolite == "lactate" & .data$detector == "lcms" ~ TRUE,
        TRUE ~ FALSE
      )
    ) |>
    dplyr::filter(.data$keep) |>
    dplyr::select(-"keep") |>
    dplyr::ungroup()
}

assemble_evap_data <- function(data_list) {
  data_list[["evap"]] |>
    dplyr::group_by(.data$experiment, .data$oxygen) |>
    dplyr::mutate(plate_mass = .data$mass - .data$mass[[1]]) |>
    dplyr::filter(.data$time != -24) |>
    dplyr::mutate(volume = 2 * .data$plate_mass / .data$plate_mass[[1]]) |>
    dplyr::filter(!is.na(.data$volume)) |>
    tidyr::nest() |>
    dplyr::mutate(
      model = purrr::map(.data$data, ~stats::lm(volume ~ time, data = .x))
    ) |>
    dplyr::mutate(pred_vol = purrr::map2(.data$model, .data$data, stats::predict)) |>
    dplyr::select(-"model") |>
    tidyr::unnest(c("data", "pred_vol")) |>
    dplyr::select(-c("volume", "plate_mass", "mass")) |>
    dplyr::rename(volume = "pred_vol") |>
    tidyr::separate("experiment", c("cell_type", "experiment", "batch", "date"), "_")
}

fill_missing_evap <- function(evap, samples) {
  evap_dup_bay <-
    evap |>
    dplyr::filter(.data$experiment == "bay") |>
    dplyr::group_by(.data$oxygen, .data$time) |>
    dplyr::summarize(volume = mean(.data$volume, na.rm = TRUE))

  evap_bay_a <-
    samples |>
    dplyr::filter(.data$experiment == "bay" & .data$batch == "a") |>
    dplyr::select(
      "cell_type",
      "experiment",
      "batch",
      "date",
      "oxygen"
    ) |>
    dplyr::distinct() |>
    dplyr::left_join(evap_dup_bay, by = "oxygen", multiple = "all")

  evap_dup_hyp <-
    evap |>
    dplyr::filter(.data$experiment == "05") |>
    dplyr::group_by(.data$oxygen, .data$time) |>
    dplyr::summarize(volume = mean(.data$volume, na.rm = TRUE)) |>
    dplyr::mutate(oxygen = replace(.data$oxygen, .data$oxygen == "0.5%", "0.2%"))

  evap_02 <-
    samples |>
    dplyr::filter(.data$experiment == "02") |>
    dplyr::select(
      "cell_type",
      "experiment",
      "batch",
      "date",
      "oxygen"
    ) |>
    dplyr::distinct() |>
    dplyr::left_join(evap_dup_hyp, by = "oxygen", multiple = "all")

  dplyr::bind_rows(evap, evap_02, evap_bay_a)
}

assemble_flux_measurements <- function(conc_clean, evap_clean) {
  abbreviations <-
    tibble::tibble(metabolite = unique(conc_clean$metabolite)) |>
    dplyr::mutate(
      abbreviation = dplyr::case_when(
        .data$metabolite == "dna" ~ "cells",
        .data$metabolite == "glucose" ~ "glc",
        .data$metabolite == "asparagine" ~ "asn",
        .data$metabolite == "cystine" ~ "cyx",
        .data$metabolite == "glutamine" ~ "gln",
        .data$metabolite == "isoleucine" ~ "ile",
        .data$metabolite == "tryptophan" ~ "trp",
        TRUE ~ stringr::str_extract(.data$metabolite, "^[a-z]{3}")
      )
    )

  conc_clean |>
    dplyr::left_join(
      evap_clean,
      by = c(
        "cell_type",
        "experiment",
        "batch",
        "date",
        "oxygen",
        "time"
      )
    ) |>
    dplyr::left_join(abbreviations, by = "metabolite") |>
    dplyr::filter(!is.na(.data$conc)) |>
    dplyr::select(
      "cell_type",
      "experiment",
      "batch",
      "date",
      "metabolite",
      "abbreviation",
      "detector",
      "type",
      "oxygen",
      "virus",
      "treatment",
      "time",
      "well",
      "conc",
      "volume"
    ) |>
    dplyr::mutate(
      treatment = factor(
        .data$treatment,
        levels = c(
          "none",
          "DMSO",
          "BAY",
          "siCTL",
          "siMYC",
          "siHIF1A",
          "siHIF2A",
          "siPHD2",
          "-GLC",
          "-GLN",
          "GLC6"
        ),
        labels = c(
          "None",
          "DMSO",
          "BAY",
          "siCTL",
          "siMYC",
          "siHIF1A",
          "siHIF2A",
          "siPHD2",
          "-GLC",
          "-GLN",
          "GLC6"
        )
      ),
      oxygen = factor(.data$oxygen, levels = c("21%", "0.5%", "0.2%")),
      virus = factor(.data$virus, levels = c("YFP", "MYC", "none"), labels = c("YFP", "MYC", "None")),
      group = dplyr::case_when(
        .data$experiment %in% c("02", "05", "bay") & .data$treatment == "None" & .data$oxygen == "21%" ~ "21%",
        .data$experiment %in% c("02", "05", "bay") & .data$treatment == "DMSO" ~ "DMSO",
        .data$experiment %in% c("02", "05", "bay") & .data$treatment == "BAY" ~ "BAY",
        .data$experiment %in% c("02", "05", "bay") & .data$oxygen == "0.5%" ~ "0.5%",
        .data$experiment %in% c("02", "05", "bay") & .data$oxygen == "0.2%" ~ "0.2%",
      ),
      group = factor(.data$group, levels = c("21%", "0.5%", "0.2%", "DMSO", "BAY")),
      metabolite = replace(.data$metabolite, .data$metabolite == "dna", "cells"),
      nmol = .data$conc * .data$volume,
      abbreviation = toupper(.data$abbreviation)
    ) |>
    dplyr::relocate("group", .before = "time") |>
    dplyr::filter(!(.data$experiment == "05-simyc" & .data$time > 48)) |>
    dplyr::filter(!(.data$experiment == "bay-myc" & .data$time > 48)) |>
    dplyr::filter(!(.data$experiment == "05-siphd" & .data$time > 48))
}

plot_raw_curves <- function(data, title, y, xlab = "Time (h)", ylab, ...) {
  ggplot2::ggplot(data) +
    ggplot2::aes(
      x = .data$time,
      y = {{y}},
      color = interaction(.data$oxygen, .data$treatment, .data$virus, sep = " | ")
    ) +
    ggplot2::geom_point(
      ggplot2::aes(shape = .data$well),
      size = 3,
      alpha = 0.3
    ) +
    geom_fit(...) +
    ggplot2::stat_summary(
      fun = "mean",
      size = 4,
      geom = "point",
      alpha = 0.8
    ) +
    ggplot2::labs(
      title = title,
      x = xlab,
      y = ylab,
      color = "condition"
    ) +
    wmo::theme_wmo()
}

geom_fit <- function(fit = "line", fo = NULL, method = NULL, ...) {
  if (fit == "reg"){
    ggplot2::geom_smooth(
      method = method,
      formula = fo,
      se = FALSE,
      ...
    )
  } else if (fit == "line"){
    ggplot2::stat_summary(
      fun = "mean",
      geom = "line",
      ...
    )
  }
}

plot_masses <- function(df, plot_function, ...) {
  df |>
    dplyr::group_by(
      .data$cell_type,
      .data$experiment,
      .data$batch,
      .data$date,
      .data$metabolite,
      .data$abbreviation
    ) |>
    tidyr::nest() |>
    dplyr::mutate(
      title = stringr::str_c(
        .data$metabolite,
        .data$cell_type,
        .data$experiment,
        .data$batch,
        .data$date,
        sep = "_"),
      plots = purrr::map2(
        .data$data,
        .data$title,
        plot_function,
        ...
      )
    )
}

plot_growth_curves <- function(flux_measurements) {
  flux_measurements |>
    dplyr::filter(.data$metabolite == "cells") |>
    plot_masses(
      plot_raw_curves,
      y = .data$conc,
      fit = "line",
      ylab = "Cell count"
    )
}

calculate_growth_rates <- function(growth_curves) {
  growth_m <- function(df){
    fit <- MASS::rlm(log(conc) ~ time, data = df, maxit = 1000)
    names(fit$coefficients) <- c("X0", "mu")
    fit
  }

  growth_curves |>
    dplyr::select(-c("title", "plots")) |>
    tidyr::unnest(c("data")) |>
    dplyr::filter(.data$time < 96) |>
    dplyr::group_by(
      .data$cell_type,
      .data$experiment,
      .data$batch,
      .data$date,
      .data$oxygen,
      .data$treatment,
      .data$virus
    ) |>
    tidyr::nest() |>
    dplyr::mutate(
      model = purrr::map(.data$data, growth_m),
      summary = purrr::map(.data$model, broom::tidy)
    ) |>
    tidyr::unnest(c("summary")) |>
    dplyr::select(-c("std.error", "statistic")) |>
    tidyr::pivot_wider(names_from = "term", values_from = "estimate") |>
    dplyr::mutate(
      X0 = exp(.data$X0),
      group = dplyr::case_when(
        .data$experiment %in% c("02", "05", "bay") & .data$treatment == "None" & .data$oxygen == "21%" ~ "21%",
        .data$experiment %in% c("02", "05", "bay") & .data$treatment == "DMSO" ~ "DMSO",
        .data$experiment %in% c("02", "05", "bay") & .data$treatment == "BAY" ~ "BAY",
        .data$experiment %in% c("02", "05", "bay") & .data$oxygen == "0.5%" ~ "0.5%",
        .data$experiment %in% c("02", "05", "bay") & .data$oxygen == "0.2%" ~ "0.2%",
      ),
      group = factor(.data$group, levels = c("21%", "0.5%", "0.2%", "DMSO", "BAY"))
    ) |>
    dplyr::arrange(dplyr::desc(.data$experiment), .data$oxygen, .data$treatment) |>
    dplyr::select(-c("data", "model")) |>
    dplyr::relocate("group", .before = "X0")
}

plot_degradation_curves <- function(flux_measurements) {
  flux_measurements |>
    dplyr::filter(.data$type == "empty") |>
    plot_masses(
      plot_raw_curves,
      y = log(.data$nmol),
      ylab = "ln(Mass (nmol))",
      fit = "reg",
      method = MASS::rlm,
      method.args = list(maxit = 100),
      fo = y ~ x
    )
}

calculate_degradation_rates <- function(df){
  degradation_m <- function(df){
    fit <- MASS::rlm(log(nmol) ~ time, data = df, maxit = 1000)
    names(fit$coefficients) <- c("intercept", "k")
    fit
  }

  df |>
    dplyr::select(-c("title", "plots")) |>
    tidyr::unnest(c("data")) |>
    dplyr::group_by(
      .data$cell_type,
      .data$experiment,
      .data$batch,
      .data$date,
      .data$metabolite,
      .data$abbreviation,
      .data$oxygen,
      .data$virus,
      .data$treatment
    ) |>
    tidyr::nest() |>
    dplyr::mutate(
      model = purrr::map(.data$data, degradation_m),
      summary = purrr::map(.data$model, broom::tidy)
    ) |>
    tidyr::unnest(c("summary")) |>
    dplyr::filter(.data$term == "k") |>
    dplyr::select(
      -c(
        "term",
        "model",
        "data",
        "std.error",
        "statistic"
      )
    ) |>
    dplyr::rename(k = "estimate")
}

clean_degradation_rates <- function(degradation_rates) {
  k <-
    degradation_rates |>
    dplyr::group_by(
      .data$metabolite,
      .data$oxygen,
      .data$treatment,
      .data$virus
    ) |>
    wmo::remove_nested_outliers("k", remove = TRUE) |>
    tidyr::nest() |>
    dplyr::mutate(
      ttest = purrr::map(.data$data, ~ t.test(.x$k, mu = 0)),
      summary = purrr::map(.data$ttest, broom::tidy)
    ) |>
    tidyr::unnest(c("summary")) |>
    dplyr::filter(.data$p.value < 0.01) |>
    dplyr::select(
      "metabolite",
      "oxygen",
      "virus",
      "treatment",
      k = "estimate"
    )

  hyp_02 <-
    k |>
    dplyr::filter(.data$oxygen == "0.5%") |>
    dplyr::mutate(oxygen = forcats::fct_recode(.data$oxygen, "0.2%" = "0.5%"))

  bay <-
    k |>
    dplyr::filter(.data$treatment == "DMSO") |>
    dplyr::mutate(treatment = forcats::fct_recode(.data$treatment, "BAY" = "DMSO"))

  dplyr::bind_rows(k, hyp_02, bay) |>
    dplyr::mutate(k = -.data$k) |>
    dplyr::arrange(.data$metabolite, .data$oxygen, .data$virus, .data$treatment) |>
    dplyr::ungroup()
}

plot_mass_curves <- function(flux_measurements) {
  flux_measurements |>
    dplyr::filter(.data$type == "cells" & .data$metabolite != "cells") |>
    plot_masses(
      plot_raw_curves,
      y = log(.data$nmol),
      ylab = "ln(Mass (nmol))",
      fit = "line"
    )
}

plot_flux_curves <- function(mass_curves, k, growth_rates) {
  mass_curves |>
    dplyr::select(-c("title", "plots")) |>
    tidyr::unnest(c("data")) |>
    dplyr::left_join(k, by = c("metabolite", "oxygen", "treatment", "virus")) |>
    dplyr::left_join(
      growth_rates,
      by = c(
        "cell_type",
        "experiment",
        "batch",
        "date",
        "oxygen",
        "treatment",
        "virus"
      )
    ) |>
    dplyr::mutate(
      k = tidyr::replace_na(.data$k, 0),
      x = exp((.data$mu + .data$k) * .data$time) - 1,
      y = .data$nmol * exp(.data$k * .data$time)
    ) |>
    plot_masses(
      plot_raw_curves,
      y = .data$y,
      xlab = expression(e^{(mu+k)*t} - 1),
      ylab = expression(M*e^{k*t}),
      fit = "reg",
      method = MASS::rlm,
      method.args = list(maxit = 100),
      fo = y ~ x
    )
}

calculate_fluxes <- function(flux_curves) {
  flux_m <- function(df){
    fit <- MASS::rlm(y ~ x, data = df, maxit = 1000)
    names(fit$coefficients) <- c("M0", "m")
    fit
  }

  flux_curves |>
    dplyr::select(-c("title", "plots")) |>
    tidyr::unnest(c("data")) |>
    dplyr::filter(.data$time < 96) |>
    dplyr::group_by(
      .data$cell_type,
      .data$experiment,
      .data$batch,
      .data$date,
      .data$metabolite,
      .data$abbreviation,
      .data$oxygen,
      .data$treatment,
      .data$virus,
      .data$k,
      .data$X0,
      .data$mu
    ) |>
    tidyr::nest() |>
    dplyr::mutate(
      model = purrr::map(.data$data, flux_m),
      summary = purrr::map(.data$model, broom::tidy)
    ) |>
    tidyr::unnest(c("summary")) |>
    dplyr::ungroup() |>
    dplyr::filter(.data$term == "m") |>
    dplyr::select(
      -c(
        "term",
        "model",
        "data",
        "std.error",
        "statistic"
      )
    ) |>
    dplyr::rename(m = "estimate") |>
    dplyr::mutate(
      flux = .data$m * (.data$mu + .data$k) / .data$X0 * 1E6,
      # group = dplyr::case_when(
      #   .data$experiment %in% c("02", "05", "bay") & .data$treatment == "None" & .data$oxygen == "21%" ~ "21%",
      #   .data$experiment %in% c("02", "05", "bay") & .data$treatment == "DMSO" ~ "DMSO",
      #   .data$experiment %in% c("02", "05", "bay") & .data$treatment == "BAY" ~ "BAY",
      #   .data$experiment %in% c("02", "05", "bay") & .data$oxygen == "0.5%" ~ "0.5%",
      #   .data$experiment %in% c("02", "05", "bay") & .data$oxygen == "0.2%" ~ "0.2%",
      # ),
      # group = factor(.data$group, levels = c("21%", "0.5%", "0.2%", "DMSO", "BAY"))
    ) |>
    dplyr::select(-c("k", "X0", "mu", "m")) |>
    dplyr::relocate("metabolite", "abbreviation") |>
    # dplyr::relocate("group", .after = "treatment") |>
    dplyr::arrange(
      .data$metabolite,
      .data$oxygen,
      .data$treatment,
      .data$date
    )
}
