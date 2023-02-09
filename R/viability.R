# viability.R

clean_viability <- function(viability_file) {
  readr::read_csv(
    viability_file,
    show_col_types = FALSE
  ) |>
    dplyr::mutate(
      viability = 100 * .data$live / (.data$dead + .data$live),
      oxygen = factor(.data$oxygen, levels = c("21%", "0.5%"))
    ) |>
    dplyr::group_by(.data$time, .data$oxygen) |>
    wmo::remove_nested_outliers("viability", remove = TRUE)
}
