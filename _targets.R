# _targets.R

# setup -------------------------------------------------------------------

devtools::load_all()
library(targets)
library(tarchetypes)

# set target options
tar_option_set(
  packages = c(
    "tidyverse",
    "scales"
  ),
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
    print_plots(
      growth_curves$plots,
      growth_curves$title,
      "fluxes/02_growth_curves"
    ),
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
    print_plots(
      degradation_curves$plots,
      degradation_curves$title,
      "fluxes/03_degradation_curves"
    ),
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

  tar_target(
    mid_files,
    raw_data_path("(a|b)_(fs|sim)_(lf|pasmc)_.*\\.csv"),
    format = "file",
    cue = tar_cue("always")
  ),
  tar_target(
    mid_clean,
    clean_mids(mid_files)
  ),
  tar_target(
    mid_correct,
    correct_mid(mid_clean)
  ),
  tar_target(
    mids,
    remove_mid_outliers(mid_correct)
  ),
  tar_target(
    mids_data,
    write_data(mids)
  ),
  tar_target(
    mid_curves,
    plot_mid_curves(mids)
  ),
  tar_target(
    mid_curve_plots,
    print_plots(mid_curves$plots, mid_curves$title, "mids"),
    format = "file"
  ),
  tar_quarto(
    mids_report,
    path = report_path("mids.qmd"),
    extra_files = c("_quarto.yml")
  ),

  # biomass -----------------------------------------------------------------

  tar_target(
    biomass_file,
    raw_data_path("cell-mass.csv"),
    format = "file"
  ),
  tar_target(
    biomass_clean,
    clean_biomass(biomass_file)
  ),
  tar_target(
    biomass,
    calculate_biomass(biomass_clean)
  ),
  tar_target(
    biomass_equations,
    calculate_biomass_equations(biomass)
  ),
  tar_target(
    biomass_equations_out,
    write_matlab_input(biomass_equations, coefs, "_biomass.csv"),
    format = "file"
  ),
  tar_quarto(
    biomass_report,
    path = report_path("biomass.qmd"),
    extra_files = c("_quarto.yml")
  ),

  # matlab ------------------------------------------------------------------

  tar_target(
    reactions_file,
    report_path("modeling/matlab-input/reactions.csv"),
    format = "file"
  ),
  tar_target(
    model_reactions,
    format_reactions(reactions_file)
  ),
  tar_target(
    model_reactions_data,
    write_data(model_reactions)
  ),
  tar_target(
    model_fluxes,
    format_fluxes(growth_rates, fluxes)
  ),
  tar_target(
    model_fluxes_data,
    write_data(model_fluxes)
  ),
  tar_target(
    model_fluxes_out,
    write_matlab_input(model_fluxes, data, "_fluxes.csv"),
    format = "file"
  ),
  tar_target(
    pruned_mids,
    format_mids(mids)
  ),
  tar_target(
    model_mids,
    summarize_mids(pruned_mids)
  ),
  tar_target(
    model_mids_out,
    write_matlab_input(model_mids, data, "_mids.csv"),
    format = "file"
  ),

  # viability ---------------------------------------------------------------

  tar_target(
    viability_file,
    raw_data_path("cell-viability.csv"),
    format = "file"
  ),
  tar_target(
    viability,
    clean_viability(viability_file)
  ),

  # blots -------------------------------------------------------------------

  tar_target(
    blot_files,
    raw_data_path("blots"),
    format = "file",
    cue = tar_cue("always")
  ),
  tar_target(
    blot_raw,
    read_data(blot_files)
  ),
  tar_target(
    blot_norm,
    normalize_densities(blot_raw)
  ),

  # mrna --------------------------------------------------------------------

  tar_target(
    mrna_files,
    raw_data_path("mrna"),
    format = "file",
    cue = tar_cue("always")
  ),
  tar_target(
    mrna_raw,
    read_data(mrna_files)
  ),
  tar_target(
    mrna_norm,
    normalize_qpcr(mrna_raw)
  ),

  # model -------------------------------------------------------------------

  tar_target(
    graph_flux_files,
    raw_data_path("model\\.csv"),
    format = "file"
  ),
  tar_target(
    graph_fluxes,
    clean_model_fluxes(graph_flux_files, model_reactions)
  ),
  tar_target(
    flux_differences,
    assemble_flux_differences(graph_fluxes)
  ),
  tar_target(
    flux_differences_data,
    write_data(flux_differences)
  ),
  tar_quarto(
    flux_differences_report,
    path = report_path("flux-differences.qmd"),
    extra_files = c("_quarto.yml")
  ),
  tar_target(
    node_file,
    raw_data_path("nodes\\.csv"),
    format = "file"
  ),
  tar_target(
    nodes,
    readr::read_csv(node_file, show_col_types = FALSE)
  ),
  tar_map(
    values = list(
      names = list(
        "lf_norm_none",
        "lf_hyp_none",
        "lf_norm_growth",
        "lf_hyp_growth",
        "lf_norm_glc",
        "lf_hyp_glc",
        "pasmc_norm_none",
        "pasmc_hyp_none",
        "pasmc_norm_growth",
        "pasmc_hyp_growth",
        "pasmc_norm_glc",
        "pasmc_hyp_glc",
        "lf_dmso_none",
        "lf_bay_none",
        "lf_dmso_growth",
        "lf_bay_growth",
        "lf_dmso_glc",
        "lf_bay_glc"
      ),
      cell = rep(c("lf", "pasmc", "lf"), each = 6),
      treat = c(rep(c("21%", "0.5%"), 6), rep(c("DMSO", "BAY"), 3)),
      normalizer = rep(c("none", "growth", "glucose"), each = 2, 3)
    ),
    names = names,
    tar_target(
      graph,
      make_graph(flux_differences, nodes, cell = cell, treat = treat, normalizer = normalizer)
    )
  ),
  tar_target(
    graph_cells_norm_none,
    make_cell_ratio_graph(graph_fluxes, nodes)
  ),
  tar_map(
    values = list(
      names = list(
        "lf_hyp_ratio",
        "lf_hyp_growth_ratio",
        "pasmc_hyp_ratio",
        "lf_bay_ratio",
        "cells_norm_none_ratio"
      ),
      graphs = rlang::syms(
        c(
          "graph_lf_hyp_none",
          "graph_lf_hyp_growth",
          "graph_pasmc_hyp_none",
          "graph_lf_bay_none",
          "graph_cells_norm_none"
        )
      ),
      captions = list(
        "Hypoxia/Normoxia",
        "Hypoxia/Normoxia\nGrowth Rate Normalized",
        "PASMC\nHypoxia/Normoxia",
        "BAY/DMSO",
        "PASMC/LF"
      ),
      edges = list(
        TRUE,
        FALSE,
        TRUE,
        TRUE,
        TRUE
      )
    ),
    names = names,
    tar_target(
      graph_ratio,
      plot_ratio_network(graphs, captions, edges)
    )
  ),
  tar_map(
    values = list(
      names = list(
        "lf_norm",
        "pasmc_norm"
      ),
      gr = rlang::syms(
        c(
          "graph_lf_norm_none",
          "graph_pasmc_norm_none"
        )
      ),
      caption = list(
        "LF\nNormoxia",
        "PASMC\nNormoxia"
      )
    ),
    names = names,
    tar_target(
      graph_raw,
      plot_network(gr, caption)
    )
  ),

  # nadp --------------------------------------------------------------------

  tar_target(
    nad_files,
    raw_data_path("nadp?_.*\\.xlsx"),
    format = "file",
    cue = tar_cue("always")
  ),
  tar_target(
    nad_data,
    assemble_flux_data(nad_files)
  ),
  tar_target(
    nad_raw,
    clean_nad(nad_data)
  ),
  tar_target(
    nad_conc_std,
    make_std_curves(
      nad_raw,
      fo = \(x) MASS::rlm(value ~ poly(conc, 2, raw = TRUE), data = x, maxit = 1000)
    )
  ),
  tar_target(
    nad_conc_std_plots,
    print_plots(nad_conc_std$plots, nad_conc_std$title, "nad/01_standard_curves"),
    format = "file"
  ),
  tar_target(
    nad_interp,
    interp_data(nad_raw, nad_conc_std)
  ),
  tar_target(
    nad_final,
    finalize_nad(nad_interp, cells_per_dna)
  ),
  tar_target(
    nad_annot,
    annot_nad(nad_final)
  ),
  tar_map(
    values = list(
      names = list(
        "nad",
        "nadh",
        "nadh_ratio",
        "nadp",
        "nadph",
        "nadph_ratio"
      ),
      metab = list(
        "NAD",
        "NADH",
        "NADH/NAD",
        "NADP",
        "NADPH",
        "NADPH/NADP"
      ),
      ylab = list(
        "NAD\n(nmol/cell)",
        "NADH\n(nmol/cell)",
        "NADH/NAD ratio",
        "NADP\n(nmol/cell)",
        "NADPH\n(nmol/cell)",
        "NADPH/NADP ratio"
      )
    ),
    names = names,
    tar_target(
      plot,
      plot_nad(nad_final, nad_annot, metab, ylab)
    )
  ),

  # rnaseq ------------------------------------------------------------------

  tar_target(
    dds,
    count_rnaseq()
  ),
  tar_target(
    pca_data,
    vst_rnaseq(dds)
  ),
  tar_target(
    rnaseq_pca,
    plot_rnaseq_pca(pca_data)
  ),
  tar_target(
    hallmark_pathways,
    get_msigdb_pathways(category = "H")
  ),
  tar_target(
    tfea,
    run_tfea(dds)
  ),
  tar_target(
    tfea_fit,
    fit_tfea(dds, tfea)
  ),
  tar_map(
    values = list(
      names = c("hyp", "bay", "hyp_bay", "int"),
      comp = list(
        "h.dmso - n.dmso",
        "n.bay - n.dmso",
        "h.bay - n.bay",
        "(h.dmso - n.dmso) - (n.bay - n.dmso)"
      ),
      xlab = c(
        "Hypoxia/Normoxia",
        "BAY/DMSO",
        "Hypoxia/Normoxia in BAY",
        "ΔHypoxia/ΔBAY"
      ),
      colors = list(
        clrs[c("0.5%", "21%")],
        clrs[c("BAY", "DMSO")],
        clrs[c("0.5%", "21%")],
        clrs[c(2, 4)]
      ),
      filename = stringr::str_c("gsea_", c("hyp", "bay", "hyp-bay", "int"))
    ),
    names = names,
    tar_target(
      deg,
      identify_deg(dds, comp)
    ),
    tar_target(
      rnaseq_vol,
      plot_rnaseq_volcano(deg, xlab = xlab)
    ),
    tar_target(
      gsea,
      run_gsea(deg, hallmark_pathways)
    ),
    tar_target(
      gsea_table,
      plot_gsea_table(gsea, title = xlab, clr = colors)
    ),
    tar_target(
      gsea_file,
      write_table(gsea_table, path = "analysis/figures/gsea/", filename),
      format = "file"
    ),
    tar_target(
      gsea_plot,
      plot_table(gsea_file)
    ),
    tar_target(
      tfea_res,
      index_tfea(tfea_fit, names)
    ),
    tar_target(
      tfea_plot,
      plot_tfea_volcanoes(tfea_res, xlab = xlab, colors = colors, nudge = 10)
    ),
    NULL
  ),
  tar_target(
    rnaseq_venn,
    plot_rnaseq_venn(deg_hyp, deg_bay)
  ),
  tar_target(
    gsea_venn,
    plot_gsea_venn(gsea_hyp, gsea_bay)
  ),
  tar_target(
    tfea_venn,
    plot_tfea_venn(tfea_res_hyp, tfea_res_bay)
  ),
  tar_quarto(
    rnaseq_report,
    path = report_path("rnaseq.qmd"),
    extra_files = c("_quarto.yml")
  ),


  NULL
)
