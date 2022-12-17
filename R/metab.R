# metab.R

read_ms_excel <- function(path, orientation = "long") {
  if (orientation == "long") {
    readxl::read_excel(
      path = path,
      sheet = 2
    ) |>
      dplyr::select(
        id = "Raw File Name",
        sample = "Sample ID",
        metabolite = "Compound Name",
        mz = "Detected Mass",
        rt = "RT",
        area = "Peak Area"
      )
  }
}

new_tbl_se <- function(
    tbl,
    a_data,
    f_names,
    f_data = NULL,
    s_names,
    s_data = NULL
) {
  structure(
    tbl,
    class = c("tbl_se", class(tbl)),
    a_data = a_data,
    f_names = f_names,
    s_names = s_names,
    f_data = f_data,
    s_data = s_data
  )
}

tbl_to_se <- function(tbl_se, assay_name){
  assay_data <-
    tbl_se |>
    dplyr::select(
      attr(tbl_se, "f_names"),
      attr(tbl_se, "s_names"),
      attr(tbl_se, "a_data")
    ) |>
    tidyr::pivot_wider(
      names_from = attr(tbl_se, "s_names"),
      values_from = attr(tbl_se, "a_data")
    ) |>
    tibble::column_to_rownames(attr(tbl_se, "f_names"))

  feature_data <-
    tbl_se |>
    dplyr::select(
      attr(tbl_se, "f_names"),
      attr(tbl_se, "f_data")
    ) |>
    dplyr::group_by(!!rlang::sym(attr(tbl_se, "f_names"))) |>
    dplyr::summarise(
      metabolite = unique(.data$metabolite),
      mz = mean(.data$mz, na.rm = TRUE),
      mz_min = min(.data$mz, na.rm = TRUE),
      mz_max = max(.data$mz, na.rm = TRUE),
      rt = mean(.data$rt, na.rm = TRUE),
      rt_min = min(.data$rt, na.rm = TRUE),
      rt_max = max(.data$rt, na.rm = TRUE)
    ) |>
    tibble::column_to_rownames(attr(tbl_se, "f_names")) |>
    {\(x) x[match(rownames(assay_data), rownames(x)), ]}()

  sample_data <-
    tbl_se |>
    dplyr::select(
      attr(tbl_se, "s_names"),
      attr(tbl_se, "s_data")
    ) |>
    dplyr::distinct() |>
    tibble::column_to_rownames(attr(tbl_se, "s_names")) |>
    {\(x) x[match(colnames(assay_data), rownames(x)), ]}()

  SummarizedExperiment::SummarizedExperiment(
    assays = assay_data,
    rowData = feature_data,
    colData = sample_data
  )
}

remove_missing_metab <- function(raw){
  qc <- SummarizedExperiment::assay(raw[, raw$type == "qc"])
  missing <- names(which(apply(qc, 1, function(x) sum(is.na(x))) > 0))
  raw[rownames(raw) %nin% missing, raw$type %nin% c("water", "blank")]
}

prepare_assay_data <- function(se){
  SummarizedExperiment::assay(se) |>
    tibble::rownames_to_column("hmdb") |>
    tidyr::pivot_longer(-"hmdb", names_to = "sample", values_to = "value") |>
    dplyr::mutate(
      value = log(.data$value),
      run_order = as.numeric(stringr::str_extract(.data$sample, "\\d{2}"))
    ) |>
    dplyr::group_by(.data$hmdb) |>
    tidyr::nest()
}

correct_drift <- function(missing){
  models <-
    missing[, missing$type == "qc"] |>
    prepare_assay_data() |>
    dplyr::mutate(
      model = purrr::map(
        .data$data,
        \(x) stats::smooth.spline(
          x = x$run_order,
          y = x$value,
          spar = 0.2)
      ),
      mean = purrr::map_dbl(
        .data$data,
        \(x) mean(x$value))
    ) |>
    dplyr::select(-"data")

  corrected <-
    missing |>
    prepare_assay_data() |>
    dplyr::left_join(models, by = "hmdb") |>
    dplyr::mutate(pred = purrr::map2(
      .data$model,
      .data$data, \(x, y) stats::predict(x, y$run_order)$y)
    ) |>
    tidyr::unnest(c(.data$data, .data$pred)) |>
    dplyr::mutate(corr = .data$value + .data$mean - .data$pred) |>
    dplyr::select("hmdb", "sample", "corr") |>
    tidyr::pivot_wider(names_from = "sample", values_from = "corr") |>
    tibble::column_to_rownames("hmdb") |>
    exp()

  SummarizedExperiment::assay(missing) <- corrected
  missing
}

se_to_tbl <- function(se, rownames) {
  samples <-
    SummarizedExperiment::colData(se) |>
    tibble::as_tibble(rownames = "id")

  features <-
    SummarizedExperiment::rowData(se) |>
    tibble::as_tibble(rownames = rownames)

  SummarizedExperiment::assay(se) |>
    tibble::as_tibble(rownames = rownames) |>
    tidyr::pivot_longer(
      -tidyselect::all_of(rownames),
      names_to = "id",
      values_to = "area"
    ) |>
    dplyr::left_join(samples, by = "id") |>
    dplyr::left_join(features, by = rownames)
}
