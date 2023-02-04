# nadp.R

clean_nad <- function(nad_data) {
  nad_data |>
    dplyr::bind_rows(.id = "metabolite") |>
    dplyr::rename(rep = replicate) |>
    clean_technical_replicates() |>
    tidyr::separate(.data$experiment, c(NA, "date"), "_")
}
