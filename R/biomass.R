# biomass.R

clean_biomass <- function(biomass_file) {
  readr::read_csv(
    biomass_file,
    show_col_types = FALSE
  ) |>
    dplyr::mutate(
      diff = .data$post - .data$pre,
      cell_mass = .data$diff / .data$cell_number * 1E12
    )
}

cell_composition <-
  tibble::tribble(
    ~metabolite, ~umol_per_gDW,
    "ALA", 0.6,
    "ARG", 0.377,
    "ASP", 0.359,
    "ASN", 0.288,
    "CYS", 0.145,
    "GLN", 0.322,
    "GLU", 0.386,
    "GLY", 0.538,
    "HIS", 0.143,
    "ILE", 0.324,
    "LEU", 0.564,
    "LYS", 0.57,
    "MET", 0.138,
    "PHE", 0.219,
    "PRO", 0.313,
    "SER", 0.43,
    "THR", 0.386,
    "TRP", 0.044,
    "TYR", 0.182,
    "VAL", 0.416,
    "glycogen", 0.279,
    "dAMP", 0.0148,
    "dGMP", 0.0099,
    "dCMP", 0.0099,
    "dTMP", 0.0148,
    "AMP", 0.033,
    "GMP", 0.0624,
    "CMP", 0.0551,
    "UMP", 0.033,
    "cholesterol", 0.018,
    "phosphatidylcholine", 0.069,
    "phosphatidylethanolamine", 0.026,
    "phosphatidylinositol", 0.01,
    "phosphatidylserine", 0.003,
    "phosphatidylglycerol", 0.001,
    "sphingomyelin", 0.008,
    "cardiolipin", 0.003
  ) |>
  dplyr::mutate(
    class = dplyr::case_when(
      .data$metabolite %in% c("dAMP", "dGMP", "AMP", "GMP") ~ "purine",
      .data$metabolite %in% c("dCMP", "CMP", "UMP") ~ "pyrimidine",
      .data$metabolite == "dTMP" ~ "thymine",
      TRUE ~ .data$metabolite
    )
  )

metabolite_composition <-
  tibble::tribble(
    ~class, ~component, ~stoichiometry,
    "purine", c("P5P", "GLY", "CO2", "MEETHF"), c(1, 1, 1, 2),
    "pyrimidine", c("P5P", "CO2", "ASP"), c(1, 1, 1),
    "thymine", c("P5P", "CO2", "ASP", "MEETHF"), c(1, 1, 1, 1),
    "glycogen", "G6P", 1,
    "cholesterol", "AcCoA", 18,
    "phosphatidylcholine", c("AcCoA", "DHAP"), c(17.43, 1),
    "phosphatidylethanolamine", c("AcCoA", "DHAP"), c(17.43, 1),
    "phosphatidylserine", c("AcCoA", "DHAP", "SER"), c(17.43, 1, 1),
    "phosphatidylinositol", c("AcCoA", "DHAP", "G6P"), c(17.43, 1, 1),
    "phosphatidylglycerol", c("AcCoA", "DHAP"), c(17.43, 2),
    "sphingomyelin", c("AcCoA", "SER"), c(17.43, 1),
    "cardiolipin", c("AcCoA", "DHAP"), c(34.86, 2)
  ) |>
  tidyr::unnest(c("component", "stoichiometry"))

umol_per_mass <-
  dplyr::left_join(cell_composition, metabolite_composition, by = "class") |>
  dplyr::mutate(
    component = dplyr::if_else(is.na(.data$component), .data$metabolite, .data$component),
    stoichiometry = tidyr::replace_na(stoichiometry, 1),
    umol_per_g = .data$umol_per_gDW * .data$stoichiometry
  ) |>
  dplyr::group_by(.data$component) |>
  dplyr::summarise(umol_per_gDW = sum(.data$umol_per_g)) |>
  dplyr::rename(metabolite = "component")

calculate_biomass <- function(biomass_clean) {
  biomass_clean |>
    dplyr::group_by(.data$cell_type, .data$date) |>
    dplyr::mutate(cell_mass = replace_outliers(.data$cell_mass)) |>
    dplyr::summarise(cell_mass = mean(.data$cell_mass, na.rm = TRUE)) |>
    dplyr::mutate(cell_mass = replace_outliers(.data$cell_mass)) |>
    dplyr::summarise(biomass = mean(.data$cell_mass, na.rm = TRUE))
}

calculate_biomass_equations <- function(biomass) {
  biomass_coefs <- function(biomass, metabolites) {
    umol_per_mass |>
      dplyr::mutate(coefficient = .data$umol_per_gDW * biomass) |>
      dplyr::filter(.data$metabolite %in% metabolites) |>
      dplyr::select("metabolite", "coefficient") |>
      dplyr::mutate(
        metabolite =
          dplyr::case_when(
            .data$metabolite %in% c("AcCoA", "ASP") ~ stringr::str_c(.data$metabolite, ".c"),
            TRUE ~ .data$metabolite
          )
      )
  }

  biomass_eq <- function(biomass_coefs) {
    biomass_coefs |>
      dplyr::mutate(
        combined = stringr::str_c(
          as.character(format(.data$coefficient, digits = 1)),
          .data$metabolite,
          sep = " ")
      ) |>
      dplyr::pull(.data$combined) |>
      stringr::str_trim() |>
      stringr::str_c(collapse = " + ") |>
      stringr::str_c(" -> Biomass")
  }

  metabs <- c(
    "ALA",
    "ASP",
    "GLU",
    "GLN",
    "P5P",
    "AcCoA",
    # "LEU",
    # "ILE",
    # "VAL",
    "SER",
    "CYS",
    "GLY",
    "MEETHF",
    "CO2",
    "DHAP",
    # "ASN",
    # "ARG",
    "G6P"
  )

  biomass |>
    dplyr::mutate(
      coefs = purrr::map(.data$biomass, biomass_coefs, metabs),
      eq = purrr::map(.data$coefs, biomass_eq)
    )
}

write_matlab_input <- function(x, column, suffix){
  path <- report_path("modeling/matlab-input")
  column <- rlang::ensym(column)

  purrr::walk2(
    x[[column]],
    x[["cell_type"]],
    \(x, y) readr::write_csv(x, file = file.path(path, stringr::str_c(y, suffix)))
  )

  purrr::map_chr(
    x[["cell_type"]],
    \(x) file.path(path, stringr::str_c(x, suffix))
  )
}
