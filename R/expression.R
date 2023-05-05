# expression.R

read_data <- function(data_files) {
  data_files[stringr::str_detect(data_files, "\\.csv$")] |>
    {\(x) rlang::set_names(
      x,
      stringr::str_extract(x, "(lf|pasmc)_(02|bay-myc|bay|05-bay|05-siphd|05-simyc|05)")
    )}() |>
    purrr::map_dfr(readr::read_csv, .id = "experiment", show_col_types = FALSE) |>
    dplyr::mutate(
      oxygen = factor(
        .data$oxygen,
        levels = c("21%", "0.5%", "0.2%"),
        ordered = TRUE
      ),
      virus = replace(.data$virus, is.na(.data$virus), "None"),
      virus = factor(
        .data$virus,
        levels = c("None", "YFP", "MYC"),
        ordered = TRUE
      ),
      treatment = replace(.data$treatment, is.na(.data$treatment), "None"),
      treatment = factor(
        .data$treatment,
        levels = c("None", "DMSO", "BAY", "siCTL", "siMYC", "siPHD2"),
        ordered = TRUE
      )
    ) |>
    dplyr::relocate("virus", .after = "oxygen")
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
            .data$treatment == min(.data$treatment) &
            .data$virus == min(.data$virus) &
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
      .data$virus,
      .data$treatment,
      .data$group,
      .data$time,
      .data$protein
    ) |>
    # wmo::remove_nested_outliers(fold_change, remove = TRUE) |>
    dplyr::relocate("group", .after = "treatment")
}

normalize_qpcr <- function(raw_mrna) {
  dct <-
    raw_mrna |>
    dplyr::mutate(gene = tolower(.data$gene)) |>
    dplyr::group_by(dplyr::across(c("experiment":"gene"))) |>
    dplyr::summarize(ct = mean(.data$ct, na.rm = TRUE)) |>
    dplyr::mutate(dct = .data$ct - .data$ct[.data$gene == "actin"]) |>
    dplyr::filter(.data$gene != "actin") |>
    dplyr::group_by(experiment) |>
    dplyr::arrange(
      .data$gene,
      .data$oxygen,
      .data$virus,
      .data$treatment,
      .data$time,
      .data$date,
      .by_group = TRUE
    )

  ctl <-
    dct |>
    dplyr::filter(
      .data$oxygen == min(.data$oxygen) &
        .data$treatment == min(.data$treatment) &
        .data$virus == min(.data$virus) &
        .data$time == min(.data$time)
    ) |>
    dplyr::group_by(dplyr::across(c("experiment", "plate", "oxygen":"gene"))) |>
    dplyr::summarise(ctl = mean(.data$dct)) |>
    dplyr::ungroup() |>
    dplyr::select("experiment", "plate", "gene", "ctl")

  dplyr::left_join(dct, ctl, by = c("experiment", "plate", "gene")) |>
    dplyr::mutate(
      ddct = .data$dct - .data$ctl,
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
      .data$virus,
      .data$treatment,
      .data$group,
      .data$time,
      .data$gene
    ) |>
    wmo::remove_nested_outliers("fold_change", remove = TRUE) |>
    dplyr::relocate("group", .after = "treatment") |>
    dplyr::rename(protein = "gene")
}

analyze_siphd_expression <- function(x, prot, exp, take_out = FALSE) {
  df <-
    x |>
    dplyr::filter(experiment == exp) |>
    dplyr::filter(protein == prot) |>
    dplyr::group_by(protein, oxygen, treatment)

  if ("date" %in% names(df)) {
    fo <- as.formula(fold_change ~ oxygen * treatment + (1 | date))
  } else if ("gel" %in% names(df)) {
    fo <- as.formula(fold_change ~ oxygen * treatment + (1 | gel))
  }

  if (take_out) {
    df <- wmo::remove_nested_outliers(df, fold_change, remove = TRUE)
  }

  annot <-
    df |>
    dplyr::group_by(protein) |>
    tidyr::nest() |>
    dplyr::mutate(
      m = purrr::map(data, ~lmerTest::lmer(fo, data = .x)),
      res = purrr::map(m, ~emmeans::emmeans(
        .x,
        "pairwise" ~ oxygen * treatment,
        simple = "each",
        adjust = "mvt",
        combine = TRUE
      )[["contrasts"]]
      ),
      out = purrr::map(res, broom::tidy)
    ) |>
    tidyr::unnest(c(out)) |>
    dplyr::select(protein, oxygen, treatment, adj.p.value) |>
    dplyr::mutate(
      oxygen = replace(oxygen, oxygen == ".", "0.5%"),
      oxygen = factor(oxygen, levels = c("21%", "0.5%")),
      treatment = factor(treatment, levels = c("DMSO", "BAY", "siCTL", "siPHD2", "siMYC")),
      y_pos = Inf,
      vjust = 1,
      lab = dplyr::case_when(
        adj.p.value < 0.05 ~ "*",
        TRUE ~ NA_character_
      )
    )

  list(data = df, annot = annot)
}

analyze_hyp_bay_densities <- function(x, prot) {
  df <-
    x |>
    dplyr::filter(experiment == "lf_05-bay") |>
    dplyr::filter(protein == prot) |>
    dplyr::group_by(protein, oxygen, treatment) |>
    # wmo::remove_nested_outliers(fold_change, remove = TRUE) |>
    identity()

  annot <-
    df |>
    dplyr::group_by(protein) |>
    tidyr::nest() |>
    dplyr::mutate(
      m = purrr::map(data, ~lmerTest::lmer(fold_change ~ oxygen * treatment + (1 | gel), data = .x)),
      res = purrr::map(m, ~emmeans::emmeans(
        .x,
        "pairwise" ~ oxygen * treatment,
        simple = "each",
        adjust = "mvt",
        combine = TRUE
      )[["contrasts"]]
      ),
      out = purrr::map(res, broom::tidy)
    ) |>
    tidyr::unnest(c(out)) |>
    dplyr::select(protein, oxygen, treatment, adj.p.value) |>
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

  list(data = df, annot = annot)
}

analyze_oemyc_expression <- function(x, prot, take_out = FALSE) {
  df <-
    x |>
    dplyr::filter(experiment == "lf_bay-myc") |>
    dplyr::filter(protein == prot) |>
    dplyr::group_by(protein, treatment, virus)

  if ("date" %in% names(df)) {
    fo <- as.formula(fold_change ~ virus * treatment + (1 | date))
  } else if ("gel" %in% names(df)) {
    fo <- as.formula(fold_change ~ virus * treatment + (1 | gel))
  }

  if (take_out) {
    df <- wmo::remove_nested_outliers(df, fold_change, remove = TRUE)
  }

  annot <-
    df |>
    dplyr::group_by(protein) |>
    tidyr::nest() |>
    dplyr::mutate(
      m = purrr::map(data, ~lmerTest::lmer(fo, data = .x)),
      res = purrr::map(m, ~emmeans::emmeans(
        .x,
        "pairwise" ~ virus * treatment,
        simple = "each",
        adjust = "mvt",
        combine = TRUE
      )[["contrasts"]]
      ),
      out = purrr::map(res, broom::tidy)
    ) |>
    tidyr::unnest(c(out)) |>
    dplyr::select(protein, virus, treatment, adj.p.value) |>
    dplyr::mutate(
      virus = replace(virus, virus == ".", "MYC"),
      virus = factor(virus, levels = c("YFP", "MYC")),
      treatment = factor(treatment, levels = c("DMSO", "BAY")),
      y_pos = Inf,
      vjust = 1,
      lab = dplyr::case_when(
        adj.p.value < 0.05 ~ "*",
        TRUE ~ NA_character_
      )
    )

  list(data = df, annot = annot)
}
