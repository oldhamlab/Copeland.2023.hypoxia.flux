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
    raw_data_path("qbias_[ab]\\.csv"),
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
    fluxes_raw_data,
    assemble_flux_data(fluxes_data_files)
  ),
  tar_target(
    fluxes_raw_dna,
    clean_dna_fluxes(fluxes_raw_data, cells_per_dna)
  ),
  tar_target(
    fluxes_raw_glc,
    clean_glc_fluxes(fluxes_raw_data)
  ),
  tar_target(
    fluxes_raw_lac,
    clean_lac_fluxes(fluxes_raw_data)
  ),
  tar_target(
    fluxes_raw_pyr,
    clean_pyr_fluxes(fluxes_raw_data)
  ),
  tar_target(
    fluxes_raw_gln,
    clean_gln_fluxes(fluxes_raw_data)
  ),
  tar_target(
    fluxes_raw_aa,
    clean_aa_fluxes(fluxes_raw_data)
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
    fluxes_raw_glc6,
    clean_glc6_fluxes(fluxes_glc6_drift, qbias_correction_factors)
  ),
  tar_target(
    conc_raw,
    clean_fluxes(
      fluxes_raw_dna,
      fluxes_raw_glc,
      fluxes_raw_lac,
      fluxes_raw_pyr,
      fluxes_raw_gln,
      fluxes_raw_aa,
      fluxes_raw_glc6
    )
  ),
  tar_target(
    conc_std,
    make_std_curves(conc_raw)
  ),
  tar_target(
    conc_std_plots,
    print_plots(conc_std$plots, conc_std$title, "fluxes/01_standard_curves"),
    format = "file"
  ),
  tar_target(
    conc_std_clean_fld,
    make_std_curves(dplyr::filter(conc_raw, !(detector == "fld" & conc > 900)))
  ),
  tar_target(
    conc_std_clean,
    clean_flux_std(conc_raw)
  ),
  tar_target(
    conc_interp,
    interp_data(conc_raw, conc_std_clean)
  ),
  tar_target(
    conc_with_missing,
    fill_missing_fluxes(conc_interp, fluxes_meta)
  ),
  tar_target(
    conc_clean,
    filter_assays(conc_with_missing)
  ),
  tar_target(
    evap_raw,
    assemble_evap_data(fluxes_raw_data)
  ),
  tar_target(
    evap_clean,
    fill_missing_evap(evap_raw, conc_clean)
  ),
  tar_target(
    flux_measurements,
    assemble_flux_measurements(conc_clean, evap_clean)
  ),
  tar_target(
    flux_measurements_data,
    write_data(flux_measurements)
  ),
  tar_target(
    growth_curves,
    plot_growth_curves(flux_measurements)
  ),
  tar_target(
    growth_curve_plots,
    print_plots(growth_curves$plots, growth_curves$title, "fluxes/02_growth_curves"),
    format = "file"
  ),
  tar_target(
    growth_rates,
    calculate_growth_rates(growth_curves)
  ),
  tar_target(
    growth_rates_data,
    write_data(growth_rates)
  ),
  tar_target(
    degradation_curves,
    plot_degradation_curves(flux_measurements)
  ),
  tar_target(
    degradation_curve_plots,
    print_plots(degradation_curves$plots, degradation_curves$title, "fluxes/03_degradation_curves"),
    format = "file"
  ),
  tar_target(
    degradation_rates,
    calculate_degradation_rates(degradation_curves)
  ),
  tar_target(
    k,
    clean_degradation_rates(degradation_rates)
  ),
  tar_target(
    k_data,
    write_data(k)
  ),
  tar_target(
    mass_curves,
    plot_mass_curves(flux_measurements)
  ),
  tar_target(
    mass_curve_plots,
    print_plots(mass_curves$plots, mass_curves$title, "fluxes/04_mass_curves"),
    format = "file"
  ),
  tar_target(
    flux_curves,
    plot_flux_curves(mass_curves, k, growth_rates)
  ),
  tar_target(
    flux_curve_plots,
    print_plots(flux_curves$plots, flux_curves$title, "fluxes/05_flux_curves"),
    format = "file"
  ),
  tar_target(
    fluxes,
    calculate_fluxes(flux_curves)
  ),
  tar_target(
    fluxes_data,
    write_data(fluxes)
  ),
  tar_quarto(
    glc6_report,
    path = report_path("glc6-lactate.qmd"),
    extra_files = c("_quarto.yml")
  ),
  tar_quarto(
    fluxes_report,
    path = report_path("fluxes.qmd"),
    extra_files = c("_quarto.yml")
  ),

# mids --------------------------------------------------------------------

  NULL
)
