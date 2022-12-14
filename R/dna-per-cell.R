# dna-per-cell.R

clean_dna_per_cell <- function(path) {
  path |>
    read_multi_excel() |>
    purrr::map(clean_technical_replicates) |>
    dplyr::bind_rows(.id = "id") |>
    tidyr::separate("id", c("cell_type", "volume"), sep = "_", convert = TRUE) |>
    dplyr::mutate(cells = 1000 * .data$cells)
}

calculate_cells_per_dna <- function(df) {
  df |>
    dplyr::filter(!(.data$cell_type == "lf" & .data$volume == "100" & .data$cells == 300000)) |>
    dplyr::filter(!(.data$cell_type == "pasmc" & .data$volume == "200" & .data$cells == 400000)) |>
    dplyr::group_by(.data$cell_type, .data$volume) |>
    tidyr::nest() |>
    dplyr::mutate(
      model = purrr::map(
        .data$data,
        ~stats::lm(
          conc ~ 0 + cells,
          data = .x,
          na.action = modelr::na.warn
        )
      ),
      glance = purrr::map(.data$model, broom::tidy)
    ) |>
    tidyr::unnest(c("glance")) |>
    dplyr::select("cell_type", "volume", slope = "estimate") |>
    dplyr::mutate(slope = 1 / .data$slope) |>
    dplyr::ungroup()
}

