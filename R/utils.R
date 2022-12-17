# utils.R

"%nin%" <- function(x, table) {
  match(x, table, nomatch = 0L) == 0L
}

raw_data_path <- function(nm) {
  dir(
    path = "data-raw",
    pattern = nm,
    all.files = TRUE,
    full.names = TRUE,
    recursive = TRUE,
    include.dirs = TRUE
  )
}

write_data <- function(...) {
  targets::tar_load(...)
  usethis::use_data(..., overwrite = TRUE)
}

report_path <- function(nm) {
  stringr::str_c("analysis/", nm)
}

read_multi_excel <- function(excel_file) {
  sheets <- readxl::excel_sheets(excel_file)
  purrr::map(sheets, ~readxl::read_excel(excel_file, sheet = .x)) |>
    rlang::set_names(sheets)
}

clean_technical_replicates <- function(tbl) {
  tidyr::pivot_longer(
    data = tbl,
    cols = "a":"c",
    names_to = "replicate",
    values_to = "value"
  ) |>
    dplyr::group_by(
      dplyr::across(-c("replicate", "value"))
    ) |>
    dplyr::mutate(value = replace_outliers(.data$value)) |>
    dplyr::summarise(value = mean(.data$value, na.rm = TRUE)) |>
    dplyr::ungroup()
}

replace_outliers <- function(vec) {
  if (stats::mad(vec, na.rm = TRUE) == 0) return (vec)
  replace(
    vec,
    abs(vec - stats::median(vec, na.rm = TRUE)) / stats::mad(vec, na.rm = TRUE) > 2,
    NA
  )
}

make_std_curves <- function(df, fo = NULL) {
  if (is.null(fo)){
    fo <- ~lm(value ~ conc, data = .x, na.action = modelr::na.warn)
  }

  df |>
    dplyr::filter(!is.na(.data$conc)) |>
    dplyr::select(tidyselect::where(\(x) all(!is.na(x)))) |>
    dplyr::group_by(dplyr::across(-c("conc", "value"))) |>
    tidyr::nest() |>
    {\(x) dplyr::mutate(
      x,
      title = stringr::str_c(!!!rlang::syms(dplyr::groups(x)), sep = "_")
    )}() |>
    dplyr::ungroup() |>
    dplyr::mutate(
      model = furrr::future_map(.data$data, fo),
      summary = furrr::future_map(.data$model, ~broom::glance(.x)),
      plots = furrr::future_map2(.data$data, .data$title, make_std_plots)
    ) |>
    dplyr::group_by(
      dplyr::across(
        -c("data", "title", "model", "summary", "plots")
      )
    )
}

make_std_plots <- function(df, title = NULL) {
  ggplot2::ggplot(df) +
    ggplot2::aes(
      x = "conc",
      y = "value"
    ) +
    ggplot2::geom_smooth(
      method = stats::lm,
      formula = y ~ x,
      color = "gray20",
      se = FALSE
    ) +
    ggplot2::geom_point(
      size = 3,
      alpha = 0.3,
      color = "blue"
    ) +
    ggplot2::stat_summary(
      fun = "mean",
      size = 4,
      geom = "point",
      alpha = 0.8,
      color = "blue"
    ) +
    ggplot2::labs(
      x = "Concentration",
      y = "Value",
      title = title
    )
}

interp_data <- function(tbl, std) {
  tbl |>
    dplyr::filter(is.na(.data$conc)) |>
    dplyr::select(-"conc") |>
    dplyr::group_by(dplyr::across(dplyr::group_vars(std))) |>
    tidyr::nest() |>
    dplyr::left_join(dplyr::select(std, "model")) |>
    dplyr::mutate(conc = purrr::map2(.data$data, .data$model, wmo::interpolate)) |>
    tidyr::unnest(c("data", "conc")) |>
    dplyr::select(-c("model", "value"))
}

my_kable <- function(data, ...) {
  kableExtra::kable(data, booktabs = TRUE, linesep = "", ...) |>
    kableExtra::kable_styling(
      latex_options = c("hold_position"),
      font_size = 9
    )
}
