# model.R

clean_model_fluxes <- function(map_flux_files, model_reactions) {
  pathways <-
    model_reactions |>
    dplyr::select(-"equation") |>
    dplyr::add_row(
      name = "BIOMASS",
      pathway = factor("biomass"),
      index = 0
    )

  map_flux_files[stringr::str_detect(map_flux_files, "(lf|pasmc)_.*_model\\.csv")] |>
    {\(x) rlang::set_names(x, stringr::str_extract(basename(x), pattern = ".*(?=_model\\.csv)"))}() |>
    purrr::map_dfr(
      readr::read_csv,
      show_col_types = FALSE,
      .id = "experiment"
    ) |>
    tidyr::separate(
      "experiment",
      c("cell_type", "treatment"),
      sep = "_"
    ) |>
    dplyr::mutate(
      treatment = dplyr::case_when(
        .data$treatment == "21" ~ "21%",
        .data$treatment == "05" ~ "0.5%",
        .data$treatment == "dmso" ~ "DMSO",
        .data$treatment == "bay" ~ "BAY"
      ),
      treatment = factor(.data$treatment, levels = c("21%", "0.5%", "DMSO", "BAY"))
    ) |>

    # identify net and exchange fluxes
    tidyr::separate("id", c("id", "type"), sep = " ", fill = "right") |>
    dplyr::mutate(type = replace(.data$type, is.na(.data$type), "net")) |>
    dplyr::select(-"se") |>

    # add pathway info
    dplyr::left_join(pathways, by = c("id" = "name")) |>
    dplyr::select(
      "cell_type",
      "treatment",
      "pathway",
      "index",
      "id",
      "type",
      "equation",
      flux = "value",
      "lb",
      "ub"
    ) |>
    dplyr::group_by(.data$cell_type, .data$treatment) |>
    dplyr::arrange(.data$treatment, .data$pathway)
}

assemble_flux_differences <- function(map_fluxes) {
  glc_norm <- normalize_fluxes(map_fluxes, "GLUT")
  growth_norm <- normalize_fluxes(map_fluxes, "BIOMASS")

  tibble::tibble(
    map_fluxes = list(
      map_fluxes,
      map_fluxes,
      map_fluxes,
      glc_norm,
      glc_norm,
      glc_norm,
      growth_norm,
      growth_norm,
      growth_norm
    ),
    cell = rep(c("lf", "lf", "pasmc"), 3),
    control = rep(c("21%", "DMSO", "21%"), 3),
    experiment = rep(c("0.5%", "BAY", "0.5%"), 3)
  ) |>
    purrr::pmap_dfr(calculate_flux_differences, .id = "normalization") |>
    dplyr::mutate(
      normalization = dplyr::case_when(
        .data$normalization %in% 1:3 ~ "none",
        .data$normalization %in% 4:6 ~ "glucose",
        .data$normalization %in% 7:9 ~ "growth"
      )
    )
}

calculate_flux_differences <- function(map_fluxes, cell, control, experiment) {
  map_fluxes |>
    dplyr::filter(.data$cell_type == cell) |>
    dplyr::filter(.data$treatment %in% c(control, experiment)) |>
    dplyr::mutate(
      treatment = dplyr::case_when(
        .data$treatment == control ~ "ctl",
        .data$treatment == experiment ~ "exp"
      )
    ) |>
    tidyr::pivot_wider(
      names_from = "treatment",
      values_from = c("flux", "lb", "ub")
    ) |>
    dplyr::mutate(
      ratio = dplyr::if_else(
        pmax(.data$lb_ctl, .data$lb_exp) - pmin(.data$ub_ctl, .data$ub_exp) > 0,
        .data$flux_exp / .data$flux_ctl,
        NA_real_
      ),
      ctl = control,
      exp = experiment
    ) |>
    dplyr::select(
      "cell_type":"equation",
      "ctl",
      "exp",
      tidyselect::contains("ctl"),
      tidyselect::contains("exp"),
      "ratio"
    )
}

make_cell_ratio_graph <- function(map_fluxes, nodes) {
  edges <-
    map_fluxes |>
    dplyr::filter(.data$treatment == "21%") |>
    dplyr::mutate(
      cell_type = dplyr::case_when(
        .data$cell_type == "lf" ~ "ctl",
        .data$cell_type == "pasmc" ~ "exp"
      ),
      equation = stringr::str_replace_all(.data$equation, "\\d+\\.?\\d+\\*", "")
    ) |>
    tidyr::pivot_wider(
      names_from = "cell_type",
      values_from = c("flux", "lb", "ub")
      ) |>
    dplyr::mutate(
      ratio = dplyr::if_else(
        pmax(.data$lb_ctl, .data$lb_exp) - pmin(.data$ub_ctl, .data$ub_exp) > 0,
        .data$flux_exp / .data$flux_ctl,
        NA_real_
      ),
      treatment_ctl = "lf",
      treatment_exp = "pasmc"
    ) |>
    dplyr::ungroup() |>
    dplyr::select(
      "pathway":"equation",
      "treatment_ctl",
      "treatment_exp",
      tidyselect::contains("ctl"),
      tidyselect::contains("exp"),
      "ratio"
    ) |>
    tidyr::pivot_longer(
      tidyselect::matches("_(ctl|exp)"),
      names_to = c(".value", NA),
      names_sep = "_"
    ) |>
    dplyr::filter(.data$treatment == "lf") |>
    parse_eq() |>
    dplyr::select(
      from = "substrate",
      to = "product",
      "flux",
      "ratio"
    ) |>
    dplyr::mutate(ratio = suppressWarnings(log(.data$ratio, base = 2))) |>
    dplyr::distinct()

  tidygraph::tbl_graph(
    nodes = nodes,
    edges = edges,
    directed = TRUE
  )
}

normalize_fluxes <- function(map_fluxes, reference_flux) {
  map_fluxes |>
    dplyr::group_by(.data$cell_type, .data$treatment) |>
    dplyr::mutate(
      lb = .data$lb / .data$flux[.data$id == reference_flux],
      ub = .data$ub / .data$flux[.data$id == reference_flux],
      flux = .data$flux / .data$flux[.data$id == reference_flux]
    )
}

make_graph <- function(flux_differences, nodes, cell, treat, normalizer) {
  edges <-
    flux_differences |>
    dplyr::rename(
      treatment_ctl = "ctl",
      treatment_exp = "exp",
    ) |>
    tidyr::pivot_longer(
      tidyselect::matches("_(ctl|exp)"),
      names_to = c(".value", NA),
      names_sep = "_"
    ) |>
    dplyr::filter(
      .data$cell_type == cell &
        .data$normalization == normalizer &
        .data$treatment == treat
    ) |>
    parse_eq() |>
    dplyr::select(
      "cell_type",
      "treatment",
      from = "substrate",
      to = "product",
      "flux",
      "ratio"
    ) |>
    dplyr::mutate(ratio = suppressWarnings(log(.data$ratio, base = 2))) |>
    dplyr::distinct()

  tidygraph::tbl_graph(
    nodes = nodes,
    edges = edges,
    directed = TRUE
  )
}

make_grid <- function(arrow, lhs, rhs) {
  forward <- tidyr::expand_grid(substrate = lhs, product = rhs)
  reverse <- tidyr::expand_grid(substrate = rhs, product = lhs)

  if (arrow == "->") {
    forward
  } else {
    reverse
  }
}

parse_eq <- function(df) {
  pattern <- "\\w?[^0-9*\\s\\.]\\w+(\\.[a-z]*)?"

  df |>

    # select relevant fluxes
    dplyr::filter(.data$type == "net" & !stringr::str_detect(.data$equation, ".ms")) |>

    # parse equation
    dplyr::mutate(
      arrow = stringr::str_extract(.data$equation, pattern = "(<->)|(<-)|(->)"),
      halved = stringr::str_split(.data$equation, pattern = .data$arrow),
      lhs = purrr::map(
        .data$halved,
        \(x) unlist(stringr::str_extract_all(x[1], pattern = pattern))
      ),
      rhs = purrr::map(
        .data$halved,
        \(x) unlist(stringr::str_extract_all(x[2], pattern = pattern))
      )
    ) |>

    # deal with directionality
    dplyr::mutate(
      arrow = dplyr::case_when(
        .data$flux >= 0 ~ "->",
        .data$flux < 0 ~ "<-"
      ),
      lb = dplyr::case_when(
        .data$flux >= 0 ~ .data$lb,
        .data$flux < 0 ~ -.data$ub
      ),
      ub = dplyr::case_when(
        .data$flux >= 0 ~ .data$ub,
        .data$flux < 0 ~ -.data$lb
      ),
      flux = abs(.data$flux),
      grid = purrr::pmap(list(.data$arrow, .data$lhs, .data$rhs), make_grid)
    ) |>
    tidyr::unnest(c("grid")) |>

    # eliminate duplicated fluxes
    dplyr::distinct()
}

plot_ratio_network <- function(graph, caption, edges = TRUE) {
  fold <- c(Inf, 5, 3, 2, 1.5, 1.1, 0.91, 0.67, 0.5, 0.33, 0.2, 0)

  fold_cuts <- log(fold, base = 2)

  lvls <-
    tidygraph::activate(graph, edges) |>
    dplyr::pull(.data$ratio) |>
    cut(
      fold_cuts,
      include.lowest = TRUE
    )

  clrs <-
    RColorBrewer::brewer.pal(11, "PRGn") |>
    rlang::set_names(levels(lvls))

  # clrs[ceiling(length(clrs)/2)] <- "grey90"

  labs <-
    c(
      "> +5.0",
      "+3.0",
      "+2.0",
      "+1.5",
      "+1.1",
      "0.0",
      "-1.1",
      "-1.5",
      "-2.0",
      "-3.0",
      "< -5.0"
    ) |>
    rev() |>
    rlang::set_names(names(clrs))

  graph <-
    graph |>
    tidygraph::activate(edges) |>
    dplyr::mutate(Flux = dplyr::case_when(
      flux > 1000 ~ "> 1000",
      flux > 100 ~ "> 100",
      flux > 10 ~ "> 10",
      TRUE ~ "> 0"
    ))

  set_gr_style()

  ggraph::ggraph(
    graph,
    circular = FALSE,
    x = tidygraph::.N()$x,
    y = tidygraph::.N()$y
  ) +
    ggraph::geom_edge_parallel(
      ggplot2::aes(
        color = cut(.data$ratio, fold_cuts, include.lowest = TRUE),
        start_cap = ggraph::label_rect(
          .data$node1.name,
          fontfamily = "Calibri",
          fontsize = 6,
          padding = ggplot2::margin(1, 1, 1, 1, "mm")
        ),
        end_cap = ggraph::label_rect(
          .data$node2.name,
          fontfamily = "Calibri",
          fontsize = 6,
          padding = ggplot2::margin(1, 1, 1, 1, "mm")
        ),
        edge_width = if (edges) .data$Flux
      ),
      arrow = ggplot2::arrow(
        length = ggplot2::unit(1, "mm"),
        angle = 30,
        type = "closed"
      ),
      linejoin = "mitre",
      lineend = "butt",
      linemitre = 10,
      key_glyph = "custom",
      show.legend = TRUE
    ) +
    ggraph::geom_node_text(
      ggplot2::aes(label = .data$name),
      size = 6/ggplot2::.pt
    ) +
    ggraph::scale_edge_color_manual(
      name = "Fold change",
      values = clrs,
      labels = labs,
      na.translate = TRUE,
      na.value = "grey90",
      drop = FALSE
    ) +
    ggraph::scale_edge_width_manual(
      values = c(0.25, 0.5, 1, 1.5),
      name = "Flux",
      guide = "legend",
      limits = c("> 0", "> 10", "> 100", "> 1000")
    ) +
    ggplot2::labs(
      x = NULL,
      y = NULL,
      title = caption
    ) +
    theme_plots() +
    ggplot2::theme(
      legend.position = c(0.07, 0.2),
      plot.title = ggplot2::element_text(
        size = ggplot2::rel(1),
        hjust = 0.5,
        face = "bold",
        margin = ggplot2::margin(b = 0)
      ),
      legend.key.height = ggplot2::unit(0.6, "lines"),
      plot.tag = ggplot2::element_text(face = "bold"),
      axis.ticks = ggplot2::element_blank(),
      panel.border = ggplot2::element_blank(),
      axis.title = ggplot2::element_blank(),
      axis.text = ggplot2::element_blank(),
      axis.text.y = ggplot2::element_blank(),
      axis.line = ggplot2::element_blank(),
      legend.text.align = 1
    ) +
    NULL
}

plot_network <- function(gr, caption) {
  set_gr_style()

  gr |>
    tidygraph::activate("edges") |>
    dplyr::filter(.data$flux > 0.1) |>
    ggraph::ggraph(
      circular = FALSE,
      x = tidygraph::.N()$x,
      y = tidygraph::.N()$y
    ) +
    ggraph::geom_edge_parallel(
      ggplot2::aes(
        # color = log(flux, base = 10),
        color = .data$flux,
        start_cap = ggraph::label_rect(
          .data$node1.name,
          fontfamily = "Calibri",
          fontsize = 6,
          padding = ggplot2::margin(1, 1, 1, 1, "mm")
        ),
        end_cap = ggraph::label_rect(
          .data$node2.name,
          fontfamily = "Calibri",
          fontsize = 6,
          padding = ggplot2::margin(1, 1, 1, 1, "mm")
        )
      ),
      edge_width = 0.5,
      arrow = ggplot2::arrow(
        length = ggplot2::unit(1, "mm"),
        type = "closed",
        angle = 30
      ),
      linejoin = "mitre",
      lineend = "butt",
      linemitre = 10,
      key_glyph = "custom"
    ) +
    ggraph::geom_node_text(
      ggplot2::aes(label = .data$name),
      size = 6/ggplot2::.pt
    ) +
    ggraph::scale_edge_color_viridis(
      name = "Flux",
      trans = "log10",
      end = 0.95,
      limits = c(1, 1000),
      na.value = "grey90",
      oob = scales::squish,
      # breaks = c(0.01, 0.1, 1, 10, 100, 1000),
      # labels = c(1, 10, 100, 1000),
      option = "plasma"
    ) +
    ggplot2::guides(
      edge_color = ggraph::guide_edge_colorbar(barwidth = 0.5, barheight = 3)
    ) +
    ggplot2::labs(
      title = caption,
      y = NULL
    ) +
    theme_plots() +
    ggplot2::theme(
      legend.position = c(0.07, 0.14),
      plot.title = ggplot2::element_text(
        size = ggplot2::rel(1),
        hjust = 0.5,
        face = "bold",
        margin = ggplot2::margin(b = 0)
      ),
      legend.key.height = ggplot2::unit(0.6, "lines"),
      plot.tag = ggplot2::element_text(face = "bold"),
      axis.ticks = ggplot2::element_blank(),
      panel.border = ggplot2::element_blank(),
      axis.title = ggplot2::element_blank(),
      axis.text = ggplot2::element_blank(),
      axis.text.y = ggplot2::element_blank(),
      axis.line = ggplot2::element_blank(),
      legend.text.align = 1
    ) +
    NULL
}

set_gr_style <- function() {
  ggraph::set_graph_style(
    family = "Calibri",
    face = "plain",
    size = 8,
    text_size = 6,
    text_colour = "black",
    title_face = "plain",
    plot_margin = ggplot2::margin(5, 5, 5, 5)
  )
}

draw_key_custom = function(data, params, size) {
  grid::segmentsGrob(
    0.1, 0.5, 0.9, 0.5,
    gp = grid::gpar(
      col = scales::alpha(data$edge_colour, data$edge_alpha),
      fill = scales::alpha(data$edge_colour, data$edge_alpha),
      lwd = data$edge_width * ggplot2::.pt,
      lty = data$edge_linetype,
      lineend = "butt",
      linejoin = "mitre"
    ),
    arrow = params$arrow
  )
}
