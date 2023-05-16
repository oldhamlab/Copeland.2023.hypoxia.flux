# _targets.R

# setup -------------------------------------------------------------------

devtools::load_all()
library(targets)
library(tarchetypes)

# set target options
tar_option_set(
  packages = c(
    "tidyverse",
    "patchwork",
    "scales",
    "grid"
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

extrafont::loadfonts(quiet = TRUE)


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
  tar_map(
    values = list(
      cell = c("lf", "lf", "pasmc"),
      experiment = c("0.5%", "BAY", "0.5%"),
      ssr1 = c(
        " SSR 391.7 [311.2-416.6] (95% CI, 362 DOF)",
        " SSR 393.5 [311.2-416.6] (95% CI, 362 DOF)",
        " SSR 575.6 [499.1-630.6] (95% CI, 563 DOF)"
      ),
      ssr2 = c(
        " SSR 334.3 [311.2-416.6] (95% CI, 362 DOF)",
        " SSR 392.4 [308.4-413.4] (95% CI, 359 DOF)",
        " SSR 521.3 [482.2-611.6] (95% CI, 545 DOF)"
      ),
      names = c(
        "lf_hyp",
        "lf_bay",
        "pasmc_hyp"
      )
    ),
    names = names,
    tar_target(
      table,
      format_flux_table(flux_differences, cell, experiment, ssr1, ssr2)
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
        # "NAD\n(nmol/cell)",
        "NAD<sup>+</sup><br>(nmol/cell)",
        "NADH\n(nmol/cell)",
        # "NADH/NAD ratio",
        "NADH/NAD<sup>+</sup> ratio",
        # "NADP\n(nmol/cell)",
        "NADP<sup>+</sup><br>(nmol/cell)",
        "NADPH\n(nmol/cell)",
        # "NADPH/NADP ratio",
        "NADPH/NADP<sup>+</sup> ratio"
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
        clrs[c("0.5%", "BAY")]
      ),
      binx = c(50, 50, 50, 50),
      biny = c(2, 2, 6, 6),
      filename = stringr::str_c("gsea_", c("hyp", "bay", "hyp-bay", "int"))
    ),
    names = names,
    tar_target(
      deg,
      identify_deg(dds, comp)
    ),
    tar_target(
      rnaseq_vol,
      plot_rnaseq_volcano(deg, xlab = xlab, binx = binx, biny = biny)
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

  # metab -------------------------------------------------------------------

  tar_target(
    metab_tar_file,
    raw_data_path("lf_05-bay_metabolomics-targeted.xlsx"),
    format = "file"
  ),
  tar_target(
    metab_tar_raw,
    format_metab_tar(metab_tar_file)
  ),
  tar_target(
    metab_tar_clean,
    remove_missing_metab(metab_tar_raw) |>
      correct_drift() |>
      quality_control() |>
      impute_missing() |>
      pqn() |>
      log_transform() |>
      annot_metabs()
  ),
  tar_target(
    metab_tar_pca,
    calc_metab_pca(metab_tar_clean)
  ),
  tar_target(
    metab_tar_pca_plot,
    plot_metab_pca(metab_tar_clean, metab_tar_pca)
  ),
  tar_target(
    metab_tar_limma,
    fit_metab_limma(metab_tar_clean)
  ),
  tar_target(
    metab_pathways_file,
    raw_data_path("metabolites.tab"),
    format = "file"
  ),
  tar_target(
    metab_pathways,
    read_metab_pathways(metab_pathways_file)
  ),
  tar_map(
    values = list(
      names = list(
        "hyp",
        "bay",
        "hyp_bay",
        "int"
      ),
      colors = list(
        clrs[c("0.5%", "21%")],
        clrs[c("BAY", "DMSO")],
        clrs[c("0.5%", "21%")],
        clrs[c("0.5%", "BAY")]
      ),
      xlab = list(
        "Hypoxia/Normoxia",
        "BAY/DMSO",
        "Hypoxia/Normoxia in BAY",
        "ΔHypoxia/ΔBAY"
      ),
      title = list(
        "Hypoxia/Normoxia",
        "BAY/DMSO",
        "Hypoxia/Normoxia in BAY",
        "ΔHypoxia/ΔBAY"
      ),
      filename = stringr::str_c("msea_", list("hyp", "bay", "hyp-bay", "int"))
    ),
    names = names,
    tar_target(
      metab_tar_res,
      index_metab_limma(metab_tar_clean, metab_tar_limma, names)
    ),
    tar_target(
      metab_tar_vol,
      plot_metab_volcano(metab_tar_res, colors = colors, xlab = xlab)
    ),
    tar_target(
      metab_tar_msea,
      run_msea(metab_tar_res, metab_pathways)
    ),
    tar_target(
      metab_tar_msea_table,
      plot_msea_table(metab_tar_msea, title, colors, names)
    ),
    tar_target(
      msea_tar_file,
      write_table(metab_tar_msea_table, path = "analysis/figures/msea/", filename),
      format = "file"
    ),
    tar_target(
      msea_tar_plot,
      plot_table(msea_tar_file)
    ),
    NULL
  ),
  tar_target(
    metab_venn,
    plot_metab_venn(metab_tar_res_hyp, metab_tar_res_bay)
  ),
  tar_target(
    tca_leading_edge,
    plot_leading_edge(
      metab_tar_res_hyp_bay,
      metab_pathways,
      "KEGG | Citrate cycle (TCA cycle) - Homo sapiens (human)"
    )
  ),
  tar_quarto(
    metab_report,
    path = report_path("metab.qmd"),
    extra_files = c("_quarto.yml")
  ),

  # myc ---------------------------------------------------------------------

  tar_map(
    values = list(
      names = c("simyc", "oemyc"),
      exp = list("05-simyc", "bay-myc"),
      intervention = list("treatment", "virus"),
      x = rlang::syms(c("treatment", "virus")),
      y = rlang::syms(c("oxygen", "treatment"))
    ),
    names = names,
    tar_target(
      myc_fluxes,
      combine_fluxes(growth_rates, fluxes, exp = exp)
    ),
    tar_target(
      myc_fluxes_annot,
      annot_myc_fluxes(myc_fluxes, intervention)
    ),
    tar_target(
      myc_growth_plot,
      plot_myc(myc_fluxes, myc_fluxes_annot, "growth", "Growth rate (/h)", x = x, fill = y)
    ),
    tar_target(
      myc_lactate_plot,
      plot_myc(myc_fluxes, myc_fluxes_annot, "lactate", "Lactate\n(fmol/cell/h)", x = x, fill = y)
    ),
    NULL
  ),

  # figure 1 ----------------------------------------------------------------

  tar_map(
    values = list(
      names = list(
        "lf_hyp_05",
        "lf_hyp_02",
        "pasmc_hyp_05",
        "lf_bay"
      ),
      timeline_path = list(
        "manuscript/figs-raw/lf_hyp_05_timeline.png",
        "manuscript/figs-raw/lf_hyp_02_timeline.png",
        "manuscript/figs-raw/pasmc_hyp_05_timeline.png",
        "manuscript/figs-raw/lf_bay_timeline.png"
      ),
      blot_path = list(
        "manuscript/figs-raw/lf_05_hif1a-ldha-blots.png",
        "manuscript/figs-raw/lf_02_hif1a-ldha-blots.png",
        "manuscript/figs-raw/pasmc_05_hif1a-ldha-blots.png",
        "manuscript/figs-raw/lf_bay_hif1a-ldha-blots.png"
      ),
      cell = c("lf", "lf", "pasmc", "lf"),
      exp1 = c("05", "02", "05", "bay"),
      exp2 = c("lf_05", "lf_02", "pasmc_05", "lf_bay"),
      filename = list(
        "Figure 01.pdf",
        "Figure 01 - figure supplement 3.pdf",
        "Figure 01 - figure supplement 4.pdf",
        "Figure 02.pdf"
      )
    ),
    names = names,
    tar_target(
      timeline_png,
      system.file(
        timeline_path,
        package = "Copeland.2023.hypoxia.flux"
      ),
      format = "file"
    ),
    tar_target(
      timeline,
      plot_image(timeline_png, scale = 1.6, hjust = 0.2, vjust = 0.1)
    ),
    tar_target(
      growth_curve,
      plot_growth_curve(flux_measurements, cell = cell, exp = exp1)
    ),
    tar_target(
      growth_rate,
      plot_growth_rates(growth_rates, cell = cell, exp = exp1)
    ),
    tar_target(
      blot_png,
      system.file(
        blot_path,
        package = "Copeland.2023.hypoxia.flux"
      ),
      format = "file"
    ),
    tar_target(
      blot,
      plot_image(blot_png, scale = 1.3, hjust = 0.2, vjust = 0)
    ),
    tar_target(
      hif1a_prot,
      plot_expression(blot_norm, exp2, "hif1a", "HIF-1α protein\n(normalized)")
    ),
    tar_target(
      ldha_prot,
      plot_expression(blot_norm, exp2, "ldha", "LDHA protein\n(normalized)")
    ),
    tar_target(
      glut1_rna,
      plot_expression(mrna_norm, exp2, "glut1", "GLUT1 mRNA\n(normalized)")
    ),
    tar_target(
      ldha_rna,
      plot_expression(mrna_norm, exp2, "ldha", "LDHA mRNA\n(normalized)")
    ),
    tar_target(
      high,
      plot_high_fluxes(fluxes, cell, exp1)
    ),
    tar_target(
      low,
      plot_low_fluxes(fluxes, cell, exp1)
    ),
    tar_target(
      fluxes_panel,
      arrange_fluxes(
        timeline,
        growth_curve,
        growth_rate,
        blot,
        hif1a_prot,
        glut1_rna,
        ldha_rna,
        ldha_prot,
        high,
        low
      )
    ),
    tar_target(
      flux_figures,
      write_figures(fluxes_panel, filename = filename),
      format = "file"
    ),
    NULL
  ),
  tar_target(
    viability_plot,
    plot_time_lines(viability, y = "viability", ylab = "Cell viability (%)", clr = "oxygen")
  ),
  tar_target(
    evap_plot,
    plot_evap_data(evap_clean)
  ),
  tar_target(
    k_plot,
    plot_k(degradation_rates, k)
  ),
  tar_target(
    f1_s1_plot,
    arrange_f1_s1(
      viability_plot,
      evap_plot,
      k_plot
    )
  ),
  tar_target(
    f1_s1_figure,
    write_figures(f1_s1_plot, "Figure 01 - figure supplement 1.pdf"),
    format = "file"
  ),
  tar_target(
    glc6_std_curve,
    plot_glc6_curve(conc_std)
  ),
  tar_target(
    glc6_mass,
    plot_glc6_mass(flux_measurements)
  ),
  tar_target(
    glc6_flux,
    plot_glc6_fluxes(fluxes)
  ),
  tar_target(
    f1_lcms_lactate,
    arrange_glc6_flux(
      glc6_std_curve,
      glc6_mass,
      glc6_flux
    )
  ),
  tar_target(
    f1_s2_figure,
    write_figures(f1_lcms_lactate, "Figure 01 - figure supplement 2.pdf"),
    format = "file"
  ),

  # figure 3 ----------------------------------------------------------------

  tar_map(
    values = list(
      cell = c("lf", "pasmc"),
      title = c("LFs", "PASMCs")
    ),
    names = cell,
    tar_target(
      substrate_plot,
      plot_substrate(growth_rates, cell, title)
    )
  ),
  tar_target(
    f3_substrate,
    arrange_substrate(substrate_plot_lf, substrate_plot_pasmc)
  ),
  tar_target(
    f3_figure,
    write_figures(f3_substrate, "Figure 03.pdf"),
    format = "file"
  ),

  # figure 4 ----------------------------------------------------------------

  tar_map(
    values = list(
      metab = c("PYR", "CIT", "CIT"),
      track = c("glc6", "glc6", "q5"),
      names = c("pyr_glc6", "cit_glc6", "cit_q5")
    ),
    names = names,
    tar_target(
      mid,
      plot_mids(pruned_mids, "lf", metab, track = track)
    )
  ),
  tar_target(
    mid_cit_q5_m5,
    plot_m5_citrate(pruned_mids)
  ),
  tar_target(
    f4_mids,
    arrange_f4(
      mid_pyr_glc6,
      mid_cit_glc6,
      mid_cit_q5,
      mid_cit_q5_m5
    )
  ),
  tar_target(
    f4_mids_figure,
    write_figures(f4_mids, "Figure 04.pdf"),
    format = "file"
  ),
  tar_map(
    values = list(
      cell = c("lf", "pasmc"),
      time = c(72, 36),
      title = c("LFs", "PASMCs"),
      filename = c(
        "Figure 04 - figure supplement 1.pdf",
        "Figure 04 - figure supplement 2.pdf"
      )
    ),
    names = cell,
    tar_target(
      mids,
      plot_all_mids(pruned_mids, cell, time, title)
    ),
    tar_target(
      mids_figure,
      write_figures(mids, filename),
      format = "file"
    )
  ),

  # figure 5 ----------------------------------------------------------------

  tar_target(
    time_course_mids,
    format_time_course_mids(model_mids)
  ),
  tar_target(
    mid_time_course,
    plot_mid_time_course(time_course_mids, "lf", c("21%", "0.5%"), "None")
  ),
  tar_target(
    f5_s1_figure,
    write_figures(mid_time_course, "Figure 05 - figure supplement 1.pdf"),
    format = "file"
  ),
  tar_target(
    f5_s2,
    arrange_f5_s2(graph_raw_lf_norm, graph_raw_pasmc_norm)
  ),
  tar_target(
    f5_s2_figure,
    write_figures(f5_s2, "Figure 05 - figure supplement 2.pdf")
  ),
  tar_target(
    f5,
    arrange_f5_s2(graph_ratio_lf_hyp_ratio, graph_ratio_lf_bay_ratio)
  ),
  tar_target(
    f5_figure,
    write_figures(f5, "Figure 05.pdf")
  ),
  tar_target(
    f5_s3_figure,
    arrange_graphs(graph_ratio_cells_norm_none_ratio) |>
      write_figures("Figure 05 - figure supplement 3.pdf")
  ),
  tar_target(
    f5_s4_figure,
    arrange_graphs(graph_ratio_pasmc_hyp_ratio) |>
      write_figures("Figure 05 - figure supplement 4.pdf")
  ),
  tar_target(
    f5_s5_figure,
    arrange_graphs(graph_ratio_lf_hyp_growth_ratio) |>
      write_figures("Figure 05 - figure supplement 5.pdf")
  ),

  # figure 6 ----------------------------------------------------------------

  tar_target(
    rc_scheme_png,
    system.file(
      "manuscript/figs-raw/13c.png",
      package = "Copeland.2023.hypoxia.flux"
    ),
    format = "file"
  ),
  tar_target(
    rc_scheme,
    plot_image(rc_scheme_png, scale = 1.3, hjust = 0, vjust = 0)
  ),
  tar_target(
    rc_fluxes,
    plot_exch_flux(graph_fluxes, "IDH")
  ),
  tar_target(
    f6_rc,
    arrange_rc(rc_scheme, rc_fluxes)
  ),
  tar_target(
    f6_rc_figure,
    write_figures(f6_rc, "Figure 06.pdf")
  ),

  # figure 7 ----------------------------------------------------------------

  tar_target(
    mct_fluxes,
    plot_exch_flux(graph_fluxes, "MCT")
  ),
  tar_target(
    lactate_mids,
    plot_lactate_mids(pruned_mids, "lf")
  ),
  tar_target(
    f7_lactate_ox,
    arrange_f7(
      mct_fluxes,
      lactate_mids
    )
  ),
  tar_target(
    f7_lactate_ox_figure,
    write_figures(f7_lactate_ox, "Figure 07.pdf")
  ),

  # figure 8 ----------------------------------------------------------------

  tar_map(
    values = list(
      names = c("05_bay", "05_siphd"),
      exp = c("05-bay", "05-siphd")
    ),
    names = names,
    tar_target(
      fluxes_stats,
      analyze_hyp_bay_fluxes(growth_rates, fluxes, exp)
    )
  ),
  tar_map(
    values = list(
      df = rlang::syms(
        c(
          rep("fluxes_stats_05_bay", 3),
          rep("fluxes_stats_05_siphd", 2)
        )
      ),
      metab = c("growth", "glucose", "lactate", "growth", "lactate"),
      ylab = c(
        "Growth Rate (/h)",
        "Glucose\n(fmol/cell/h)",
        "Lactate\n(fmol/cell/h)",
        "Growth Rate (/h)",
        "Lactate\n(fmol/cell/h)"
      ),
      names = c(
        "bay_growth",
        "bay_glucose",
        "bay_lactate",
        "phd_growth",
        "phd_lactate"
      )
    ),
    names = names,
    tar_target(
      hyp_bay_fluxes,
      plot_hyp_bay_fluxes(df$data, df$annot, metab, ylab)
    )
  ),
  tar_target(
    f8_hyp_bay,
    arrange_f8(
      hyp_bay_fluxes_bay_growth,
      hyp_bay_fluxes_bay_glucose + ggplot2::scale_y_reverse(),
      hyp_bay_fluxes_bay_lactate
    )
  ),
  tar_target(
    f8_hyp_bay_figure,
    write_figures(f8_hyp_bay, "Figure 08.pdf")
  ),
  tar_map(
    values = list(
      df = rlang::syms(
        c(
          "blot_norm",
          "blot_norm",
          "mrna_norm",
          "mrna_norm",
          "mrna_norm"
        )
      ),
      protein = c(
        "ldha",
        "myc",
        "glut1",
        "ldha",
        "cdkn1a"
      ),
      label = c(
        "LDHA protein\n(normalized)",
        "MYC protein\n(normalized)",
        "GLUT1 mRNA\n(normalized)",
        "LDHA mRNA\n(normalized)",
        "CDKN1A mRNA\n(normalized)"
      ),
      names = c(
        "prot_ldha",
        "prot_myc",
        "mrna_glut1",
        "mrna_ldha",
        "mrna_cdkn1a"
      )
    ),
    names = names,
    tar_target(
      siphd_annot,
      analyze_siphd_expression(df, exp = "lf_05-siphd", prot = protein)
    ),
    tar_target(
      siphd_plot,
      plot_siphd(siphd_annot$data, siphd_annot$annot, prot = protein, ylab = label)
    )
  ),
  tar_target(
    siphd_blot_1_png,
    system.file(
      "manuscript/figs-raw/lf_05-siphd_phd-hif-ldha-blots.png",
      package = "Copeland.2023.hypoxia.flux"
    ),
    format = "file"
  ),
  tar_target(
    siphd_blot_1,
    plot_image(siphd_blot_1_png, scale = 1.25, hjust = 0.1, vjust = 0)
  ),
  tar_target(
    f8_s1_figure,
    arrange_f8_s1(
      hyp_bay_fluxes_phd_growth,
      siphd_blot_1,
      siphd_plot_prot_ldha,
      siphd_plot_mrna_glut1,
      siphd_plot_mrna_ldha,
      hyp_bay_fluxes_phd_lactate
    ) |>
      write_figures("Figure 08 - figure supplement 1.pdf")
  ),

  # figure 9 ----------------------------------------------------------------

  tar_target(
    f9_figure,
    arrange_f9(
      metab_tar_pca_plot,
      metab_tar_vol_hyp_bay,
      msea_tar_plot_hyp_bay,
      tca_leading_edge
    ) |>
      write_figures("Figure 09.pdf")
  ),
  tar_target(
    f9_s1_figure,
    arrange_f9_s1(
      metab_tar_vol_hyp,
      metab_tar_vol_bay,
      metab_venn,
      msea_tar_plot_hyp,
      msea_tar_plot_bay
    ) |>
      write_figures("Figure 09 - figure supplement 1.pdf")
  ),
  tar_target(
    f9_s2_figure,
    arrange_f9_s2(
      plot_nad,
      plot_nadh,
      plot_nadh_ratio,
      plot_nadp,
      plot_nadph,
      plot_nadph_ratio
    ) |>
      write_figures("Figure 09 - figure supplement 2.pdf")
  ),

  # figure 10 ---------------------------------------------------------------

  tar_map(
    values = list(
      sets = list(
        c("HALLMARK_E2F_TARGETS", "HALLMARK_G2M_CHECKPOINT"),
        c("HALLMARK_MYC_TARGETS_V1", "HALLMARK_MYC_TARGETS_V2"),
        "HALLMARK_OXIDATIVE_PHOSPHORYLATION"
      ),
      titles = list(
        "E2F Targets and\nG2/M Checkpoint",
        "MYC Targets",
        "Oxidative Phosphorylation"
      ),
      names = list(
        "e2f",
        "myc",
        "oxphos")
    ),
    names = names,
    tar_target(
      gsea_deg,
      plot_pathway_volcanoes(
        deg_hyp_bay,
        hallmark_pathways,
        sets = sets,
        title = titles
      )
    )
  ),
  tar_target(
    f10_figure,
    arrange_f10(
      rnaseq_pca,
      rnaseq_vol_hyp_bay,
      gsea_plot_hyp_bay,
      gsea_deg_e2f,
      gsea_deg_myc,
      tfea_plot_hyp_bay
    ) |>
      write_figures("Figure 10.pdf")
  ),
  tar_target(
    f10_s1_figure,
    arrange_f10_s1(
      rnaseq_vol_hyp,
      rnaseq_vol_bay,
      rnaseq_venn,
      gsea_venn,
      gsea_plot_hyp,
      gsea_plot_bay,
      tfea_plot_hyp,
      tfea_plot_bay,
      tfea_venn
    ) |>
      write_figures("Figure 10 - figure supplement 1.pdf")
  ),
  tar_target(
    cdkn1a,
    plot_cdkn1a(dds)
  ),
  tar_target(
    f10_s2_figure,
    arrange_f10_s2(
      cdkn1a,
      siphd_plot_mrna_cdkn1a
    ) |>
      write_figures("Figure 10 - figure supplement 2.pdf")
  ),

  # figure 11 ---------------------------------------------------------------

  tar_target(
    myc_blot_file,
    manuscript_path("figs-raw/lf_05-bay_myc-blots.png"),
    format = "file"
  ),
  tar_target(
    myc_blot,
    plot_image(myc_blot_file)
  ),
  tar_target(
    myc_blot_quant,
    analyze_hyp_bay_densities(blot_norm, "myc")
  ),
  tar_target(
    myc_blot_plot,
    plot_hyp_bay_densities(
      myc_blot_quant$data,
      myc_blot_quant$annot,
      "myc",
      "MYC protein\n(normalized)"
    )
  ),
  tar_target(
    model_image_file,
    manuscript_path("figs-raw/working-model_2.png"),
    format = "file"
  ),
  tar_target(
    model_image,
    plot_image(model_image_file, scale = 1.7, hjust = 0.2, vjust = 0.2)
  ),
  tar_target(
    simyc_image_file,
    manuscript_path("figs-raw/lf_05-simyc_myc-blots.png"),
    format = "file"
  ),
  tar_target(
    simyc_image,
    plot_image(simyc_image_file)
  ),
  tar_target(
    oemyc_image_file,
    manuscript_path("figs-raw/lf_bay-myc_myc-blots.png"),
    format = "file"
  ),
  tar_target(
    oemyc_image,
    plot_image(oemyc_image_file)
  ),
  tar_target(
    f11_figure,
    arrange_f11(
      myc_blot,
      myc_blot_plot,
      model_image,
      simyc_image,
      myc_growth_plot_simyc + ggplot2::geom_hline(yintercept = 0, size = 0.25),
      myc_lactate_plot_simyc,
      oemyc_image,
      myc_growth_plot_oemyc,
      myc_lactate_plot_oemyc
    ) |>
      write_figures("Figure 11.pdf")
  ),
  tar_target(
    simyc_blot_2_png,
    system.file(
      "manuscript/figs-raw/lf_05-simyc_hif-ldha-blots.png",
      package = "Copeland.2023.hypoxia.flux"
    ),
    format = "file"
  ),
  tar_target(
    simyc_blot_2,
    plot_image(simyc_blot_2_png, scale = 1.25, hjust = -0.15, vjust = 0)
  ),
  tar_map(
    values = list(
      df = rlang::syms(
        c(
          "blot_norm",
          "blot_norm",
          "mrna_norm",
          "mrna_norm"
        )
      ),
      protein = c(
        "hif1a",
        "ldha",
        "glut1",
        "ldha"
      ),
      label = c(
        "HIF-1α protein\n(normalized)",
        "LDHA protein\n(normalized)",
        "GLUT1 mRNA\n(normalized)",
        "LDHA mRNA\n(normalized)"
      ),
      names = c(
        "simyc_prot_hif1a",
        "simyc_prot_ldha",
        "simyc_mrna_glut1",
        "simyc_mrna_ldha"
      )
    ),
    names = names,
    tar_target(
      stats,
      analyze_siphd_expression(df, exp = "lf_05-simyc", prot = protein)
    ),
    tar_target(
      plot,
      plot_simyc(stats$data, stats$annot, prot = protein, ylab = label)
    )
  ),
  tar_target(
    siphd_blot_2_png,
    system.file(
      "manuscript/figs-raw/lf_05-siphd_myc-blots.png",
      package = "Copeland.2023.hypoxia.flux"
    ),
    format = "file"
  ),
  tar_target(
    siphd_blot_2,
    plot_image(siphd_blot_2_png, scale = 1.25, hjust = -0.15, vjust = 0)
  ),
  tar_target(
    f11_s1_figure,
    arrange_f11_s1(
      siphd_blot_2,
      siphd_plot_prot_myc
    ) |>
      write_figures("Figure 11 - figure supplement 1.pdf")
  ),
  tar_target(
    oemyc_blot_2_png,
    system.file(
      "manuscript/figs-raw/lf_bay-myc_hif-ldha-blots.png",
      package = "Copeland.2023.hypoxia.flux"
    ),
    format = "file"
  ),
  tar_target(
    oemyc_blot_2,
    plot_image(oemyc_blot_2_png, scale = 1.2, hjust = -0.07, vjust = 0)
  ),
  tar_map(
    values = list(
      df = rlang::syms(
        c(
          "blot_norm",
          "blot_norm",
          "mrna_norm",
          "mrna_norm"
        )
      ),
      protein = c(
        "hif1a",
        "ldha",
        "glut1",
        "ldha"
      ),
      label = c(
        "HIF-1α protein\n(normalized)",
        "LDHA protein\n(normalized)",
        "GLUT1 mRNA\n(normalized)",
        "LDHA mRNA\n(normalized)"
      ),
      names = c(
        "oemyc_prot_hif1a",
        "oemyc_prot_ldha",
        "oemyc_mrna_glut1",
        "oemyc_mrna_ldha"
      )
    ),
    names = names,
    tar_target(
      stats,
      analyze_oemyc_expression(df, prot = protein)
    ),
    tar_target(
      plot,
      plot_oemyc(stats$data, stats$annot, prot = protein, ylab = label)
    )
  ),
  tar_target(
    f11_s2_figure,
    arrange_f11_s2(
      simyc_blot_2,
      plot_simyc_prot_hif1a,
      plot_simyc_prot_ldha,
      plot_simyc_mrna_glut1,
      plot_simyc_mrna_ldha,
      oemyc_blot_2,
      plot_oemyc_prot_hif1a,
      plot_oemyc_prot_ldha,
      plot_oemyc_mrna_glut1,
      plot_oemyc_mrna_ldha
    ) |>
      write_figures("Figure 11 - figure supplement 2.pdf")
  ),

  # figure 12 ---------------------------------------------------------------

  tar_map(
    values = list(cell = c("lf", "pasmc")),
    tar_target(
      dna_curve,
      plot_cells_per_dna(dna_per_cell_clean, cell)
    )
  ),
  tar_target(
    f12_figure,
    arrange_f12(
      dna_curve_lf,
      dna_curve_pasmc
    ) |>
      write_figures("Figure 12.pdf")
  ),

  # resources table ---------------------------------------------------------

  tar_target(
    resources_file,
    system.file(
      "manuscript/figs-raw/resources.csv",
      package = "Copeland.2023.hypoxia.flux"
    )
  ),
  tar_target(
    resources_table,
    create_resources(resources_file)
  ),

  # manuscript --------------------------------------------------------------

  tar_target(
    template,
    system.file(
      "manuscript/template.docx",
      package = "Copeland.2023.hypoxia.flux"
    ),
    format = "file"
  ),
  tar_target(
    pkg_citations,
    write_pkg_cites(),
    cue = tar_cue(mode = "always")
  ),
  tar_target(
    csl,
    system.file(
      "manuscript/elife.csl",
      package = "Copeland.2023.hypoxia.flux"
    ),
    format = "file"
  ),
  tar_target(
    refs,
    rbbt::bbt_update_bib(
      path = "manuscript/manuscript.qmd",
      ignore = c("R-base"),
      path_bib = "manuscript/library.bib"
    ),
    cue = tar_cue("always")
  ),
  tar_quarto(
    manuscript,
    path = manuscript_path("manuscript.qmd"),
  ),
  NULL
)
