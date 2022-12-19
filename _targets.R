# _targets.R

# setup -------------------------------------------------------------------

devtools::load_all()
library(targets)
library(tarchetypes)

# set target options
tar_option_set(
  packages = c("tidyverse"),
  format = "qs"
)

# source R scripts
tar_source()

# parallel processing
options(clustermq.scheduler = "multicore")
future::plan(future.callr::callr(workers = future::availableCores() - 1))

# quiet packages
options(
  dplyr.summarise.inform = FALSE
)


# targets -----------------------------------------------------------------

list(

  # picogreen ---------------------------------------------------------------

  tar_target(
    dna_per_cell_file,
    raw_data_path("dna-per-cell.xlsx"),
    format = "qs"
  ),
  tar_target(
    dna_per_cell_raw,
    clean_dna_per_cell(dna_per_cell_file)
  ),
  tar_target(
    dna_per_cell_std,
    make_std_curves(dna_per_cell_raw)
  ),
  tar_target(
    dna_per_cell_clean,
    interp_data(dna_per_cell_raw, dna_per_cell_std)
  ),
  tar_target(
    cells_per_dna,
    calculate_cells_per_dna(dna_per_cell_clean)
  ),
  tar_target(
    cells_per_dna_data,
    write_data(cells_per_dna),
    format = "qs"
  ),
  tar_quarto(
    dna_per_cell_report,
    path = report_path("dna-per-cell.qmd"),
    extra_files = c("_quarto.yml")
  ),

  # qbias -------------------------------------------------------------------

  tar_target(
    qbias_files,
    raw_data_path("q-bias-correction_[ab]\\.csv"),
    format = "file"
  ),
  tar_target(
    qbias_ratios,
    read_qbias(qbias_files)
  ),
  tar_target(
    qbias_correction_factors,
    calculate_correction_factors(qbias_ratios, predicted_ratios)
  ),
  tar_target(
    qbias_correction_factors_data,
    write_data(qbias_correction_factors)
  ),
  tar_quarto(
    qbias_report,
    path = report_path("qbias.qmd"),
    extra_files = c("_quarto.yml")
  ),

  # fluxes ------------------------------------------------------------------

  tar_target(
    fluxes_meta_files,
    raw_data_path("(lf|pasmc)_.*_meta\\.csv"),
    format = "file",
    cue = tar_cue("always")
  ),
  tar_target(
    fluxes_meta,
    clean_flux_meta(fluxes_meta_files)
  ),
  tar_target(
    fluxes_data_files,
    raw_data_path("(lf|pasmc)_.*_[a-z]_\\d{4}-\\d{2}-\\d{2}\\.xlsx"),
    format = "file",
    cue = tar_cue("always")
  ),
  tar_target(
    fluxes_data,
    assemble_flux_data(fluxes_data_files)
  ),
  tar_target(
    fluxes_dna_raw,
    clean_dna_fluxes(fluxes_data, cells_per_dna)
  ),
  tar_target(
    conc_raw,
    clean_fluxes(fluxes_data, cells_per_dna)
  ),
  tar_target(
    fluxes_glc6_files,
    raw_data_path("lf_substrate_b_lactate\\.xlsx")
  ),
  tar_target(
    fluxes_glc6_se,
    format_glc6_raw(fluxes_glc6_files)
  ),
  tar_target(
    fluxes_glc6_drift,
    correct_drift(fluxes_glc6_se) |>
      se_to_tbl(rownames = "metabolite")
  ),
  tar_target(
    fluxes_glc6_raw,
    clean_glc6_fluxes(fluxes_glc6_drift)
  ),
  tar_target(
    conc_std,
    make_std_curves(dplyr::bind_rows(conc_raw, fluxes_glc6_raw))
  ),
  tar_target(
    conc_std_plots,
    print_plots(conc_std$plots, conc_std$title, "fluxes/01_standard_curves"),
    format = "file"
  ),
  NULL
)
