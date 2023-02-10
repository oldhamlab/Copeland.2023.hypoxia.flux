# data.R

#' Cell count per DNA
#'
#' A data set containing the slopes to interpolate cell number from DNA
#' measurements.
#'
#' \describe{
#'   \item{cell_type}{
#'   `lf`    = lung fibroblasts \cr
#'   `pasmc` = pulmonary artery smooth muscle cells}
#'   \item{volume}{volume of lysis buffer used to extract DNA}
#'   \item{slope}{cells / ng DNA}
#'   }
"cells_per_dna"

#' Quadrupole bias correction factors
#'
#' A data set containing experimentally determined correction factors to adjust
#' peak areas for isotopes detected in selected ion monitoring mode. Isotope peak
#' areas are multiplied by the correction factor. The corrected peak area is
#' subsequently used to determine the mass isotope distribution.
#'
#'  \describe{
#'    \item{batch}{correction factors associated with specific experiments}
#'    \item{metabolite}{name of measurement}
#'    \item{M}{isotope}
#'    \item{cf}{correction factor}
#'   }
"qbias_correction_factors"

#' Measurements for extracellular flux calculations
#'
#' A data set containing the interpolated metabolite concentrations, cell counts,
#' and evaporation volumes for extracellular flux determinations.
#'
#' \describe{
#'   \item{cell_type}{
#'   `lf`    = lung fibroblasts \cr
#'   `pasmc` = pulmonary artery smooth muscle cells}
#'   \item{experiment}{
#'   `02`        = 0.2% oxygen for hypoxia \cr
#'   `05`        = 0.5% oxygen for hypoxia \cr
#'   `bay`       = molidustat treatment \cr
#'   `05-bay`    = 0.5% oxygen plus molidustat \cr
#'   `05-simyc`  = 0.5% oxygen plus siMYC \cr
#'   `05-siphd`  = 0.5% oxygen plus siPHD2 \cr
#'   `bay-myc`   = molidustat plus MYC overexpression \cr
#'   `substrate` = glucose or glutamine dropout and <sup>13</sup>C-glucose labeling}
#'   \item{batch}{a group of biological replicates analyzed similarly}
#'   \item{date}{start date of an experiment}
#'   \item{metabolite}{name of measurement}
#'   \item{abbreviation}{abbreviated metabolite name}
#'   \item{detector}{
#'   `fld`       = HPLC fluorescence detector \cr
#'   `mwd`       = HPLC multi-wavelength detector \cr
#'   `enzyme`    = enzymatic assay \cr
#'   `hplc`      = OPD-derivatized HPLC detection \cr
#'   `lcms`      = liquid chromatography-mass spectrometry assay \cr
#'   `picogreen` = fluorescence dye labeling of DNA}
#'   \item{type}{
#'   `cells` = conditioned medium \cr
#'   `empty` = unconditioned medium}
#'   \item{oxygen}{ambient oxygen level}
#'   \item{virus}{
#'   `None` = no virus \cr
#'   `YFP`  = YFP control virus \cr
#'   `MYC`  = MYC overexpression virus}
#'   \item{treatment}{
#'   `None`    = no treatment \cr
#'   `DMSO`    = 0.1% DMSO \cr
#'   `BAY`     = 10 μM molidustat \cr
#'   `siCTL`   = non-targeting control siRNA \cr
#'   `siHIF1A` = siRNA targeting HIF1A \cr
#'   `siHIF2A` = siRNA targeting HIF2A \cr
#'   `siMYC`   = siRNA targeting MYC \cr
#'   `siPHD2`  = siRNA targeting PHD2 \cr
#'   `-GLC`    = MCDB131 medium lacking glucose \cr
#'   `-GLN`    = MCDB131 medium lacking glutamine \cr
#'   `GLC6`    = MDCB131 medium supplemented with \[U-<sup>13</sup>C<sub>6</sub>\]-glucose}
#'   \item{time}{hours}
#'   \item{well}{denotes technical replicates in each experiment}
#'   \item{conc}{measured in number for cells and μM for metabolites}
#'   \item{volume}{measured in mL, extrapolated from evaporation controls}
#'   \item{nmol}{metabolite mass}
#'   }
"flux_measurements"

#' Growth rate parameters
#'
#' A data set containing the growth rate and X0 values from linear fitting of
#' cell count data.
#'
#' \describe{
#'   \item{cell_type}{
#'   `lf`    = lung fibroblasts \cr
#'   `pasmc` = pulmonary artery smooth muscle cells}
#'   \item{experiment}{
#'   `02`        = 0.2% oxygen for hypoxia \cr
#'   `05`        = 0.5% oxygen for hypoxia \cr
#'   `bay`       = molidustat treatment \cr
#'   `05-bay`.   = 0.5% oxygen plus molidustat \cr
#'   `05-simyc`  = 0.5% oxygen plus siMYC \cr
#'   `05-siphd`  = 0.5% oxygen plus siPHD2 \cr
#'   `bay-myc`   = molidustat plus MYC overexpression \cr
#'   `substrate` = glucose or glutamine dropout and <sup>13</sup>C-glucose labeling}
#'   \item{batch}{a group of biological replicates analyzed similarly}
#'   \item{date}{start date of an experiment}
#'   \item{oxygen}{ambient oxygen level}
#'   \item{virus}{
#'   `None`    = no virus \cr
#'   `YFP`     = YFP control virus \cr
#'   `MYC`     = MYC overexpression virus}
#'   \item{treatment}{
#'   `None`    = no treatment \cr
#'   `DMSO`    = 0.1% DMSO \cr
#'   `BAY`     = 10 μM molidustat \cr
#'   `siCTL`   = non-targeting control siRNA \cr
#'   `siHIF1A` = siRNA targeting HIF1A \cr
#'   `siHIF2A` = siRNA targeting HIF2A \cr
#'   `siMYC`   = siRNA targeting MYC \cr
#'   `siPHD2`  = siRNA targeting PHD2 \cr
#'   `-GLC`    = MCDB131 medium lacking glucose \cr
#'   `-GLN`.   = MCDB131 medium lacking glutamine \cr
#'   `GLC6`    = MDCB131 medium supplemented with \[U-<sup>13</sup>C<sub>6</sub>\]-glucose}
#'   \item{X0}{cell count at time 0}
#'   \item{mu}{growth rate per hour}
#'   }
"growth_rates"

#' Spontaneous degradation and accumulation rates
#'
#' A data set containing the rates of significant metabolite accumulation or
#' degradation in cell-free medium.
#'
#' \describe{
#'   \item{metabolite}{name}
#'   \item{oxygen}{ambient oxygen level}
#'   \item{treatment}{
#'   `None` = no treatment \cr
#'   `DMSO` = 0.1% DMSO \cr
#'   `BAY`  = 10 μM molidustat}
#'   \item{k}{accumulation or decay rate per hour, positive values indicate decay}
#'   }
"k"

#' Calculated fluxes
#'
#' A data set containing metabolite fluxes.
#'
#'  \describe{
#'   \item{metabolite}{name of measurement}
#'   \item{abbreviation}{abbreviated metabolite name}
#'   \item{cell_type}{
#'   `lf`    = lung fibroblasts \cr
#'   `pasmc` = pulmonary artery smooth muscle cells}
#'   \item{experiment}{
#'   `02`        = 0.2% oxygen for hypoxia \cr
#'   `05`        = 0.5% oxygen for hypoxia \cr
#'   `bay`       = molidustat treatment \cr
#'   `05-bay`    = 0.5% oxygen plus molidustat \cr
#'   `05-simyc`  = 0.5% oxygen plus siMYC \cr
#'   `05-siphd`  = 0.5% oxygen plus siPHD2 \cr
#'   `bay-myc`   = molidustat plus MYC overexpression \cr
#'   `substrate` = glucose or glutamine dropout and <sup>13</sup>C-glucose labeling}
#'   \item{batch}{a group of biological replicates analyzed similarly}
#'   \item{date}{start date of an experiment}
#'   \item{oxygen}{ambient oxygen level}
#'   \item{virus}{
#'   `None` = no virus \cr
#'   `YFP`  = YFP control virus \cr
#'   `MYC`  = MYC overexpression virus}
#'   \item{treatment}{
#'   `None`    = no treatment \cr
#'   `DMSO`    = 0.1% DMSO \cr
#'   `BAY`     = 10 μM molidustat \cr
#'   `siCTL`   = non-targeting control siRNA \cr
#'   `siHIF1A` = siRNA targeting HIF1A \cr
#'   `siHIF2A` = siRNA targeting HIF2A \cr
#'   `siMYC`   = siRNA targeting MYC \cr
#'   `siPHD2`  = siRNA targeting PHD2 \cr
#'   `-GLC`    = MCDB131 medium lacking glucose \cr
#'   `-GLN`    = MCDB131 medium lacking glutamine \cr
#'   `GLC6`    = MDCB131 medium supplemented with \[U-<sup>13</sup>C<sub>6</sub>\]-glucose}
#'   \item{group}{
#'   }
#'   \item{flux}{fmol / cell / h, positive fluxes indicate secretion, negative fluxes
#'   indicate uptake}
#'   }
"fluxes"

#' Mass isotope distributions
#'
#' A data set containing the mass isotope distributions corrected for natural
#' abundance for all experiments.
#'
#'  \describe{
#'   \item{method}{
#'   `fs`  = full scan \cr
#'   `sim` = selected ion monitoring}
#'   \item{cell_type}{
#'   `lf`    = lung fibroblasts \cr
#'   `pasmc` = pulmonary artery smooth muscle cells}
#'   \item{date}{start date of an experiment}
#'   \item{tracer}{
#'   `glc2` = \[1,2-<sup>13</sup>C<sub>2</sub>\]-glucose \cr
#'   `glc6` = \[U-<sup>13</sup>C<sub>6</sub>\]-glucose \cr
#'   `q5`   = \[U-<sup>13</sup>C<sub>5</sub>\]-glutamine \cr
#'   `lac3` = \[U-<sup>13</sup>C<sub>3</sub>\]-lactate}
#'   \item{oxygen}{ambient oxygen level}
#'   \item{treatment}{
#'   `None` = no treatment \cr
#'   `DMSO` = 0.1% DMSO \cr
#'   `BAY`  = 10 μM molidustat \cr}
#'   \item{time}{hours}
#'   \item{metabolite}{name of measurement}
#'   \item{isotope}{}
#'   \item{mid}{mole fraction of the isotope}
#'   }
"mids"

#' Model reactions
#'
#' A data table of model reactions formatted for INCA analysis.
#'
#'  \describe{
#'   \item{index}{reaction identifier}
#'   \item{name}{reaction name}
#'   \item{equation}{equation for INCA input}
#'   \item{pathway}{metabolic pathway}
#'   }
"model_reactions"

#' Model fluxes
#'
#' A data table of experimentally determined extracellular fluxes used to
#' constrain the metabolic model.
#'
#'  \describe{
#'   \item{cell_type}{
#'   `lf`    = lung fibroblasts \cr
#'   `pasmc` = pulmonary artery smooth muscle cells}
#'   \item{data}{
#'   `metabolite` = metabolite name \cr
#'   `oxygen` = ambient oxygen level \cr
#'   `treatment` = `None`, `DMSO`, `BAY` \cr
#'   `mean` = average flux from multiple biological replicates \cr
#'   `se` = standard error of flux from multiple biological replicates \cr
#'   `name` = corresponds to name of model reaction}
#'   }
"model_fluxes"

#' Flux differences
#'
#' A data table of summarized fluxes and significant differences. This also
#' includes data for flux differences normalized to growth rate or glucose
#' consumptions.
#'
#'  \describe{
#'   \item{normalization}{
#'   `none` = no normalization \cr
#'   `glucose` = fluxes normalized to GLUT flux \cr
#'   `growth` = fluxed normalized to biomass flux}
#'   \item{cell_type}{
#'   `lf`    = lung fibroblasts \cr
#'   `pasmc` = pulmonary artery smooth muscle cells}
#'   \item{pathway}{metabolic pathway}
#'   \item{index}{numeric ID for a metabolic reaction}
#'   \item{id}{name for a metabolic reaction}
#'   \item{type}{
#'   `net` = net flux \cr
#'   `exch` = exchange flux}
#'   \item{equation}{flux reaction equation}
#'   \item{ctl}{control condition}
#'   \item{exp}{experimental condition}
#'   \item{flux_ctl}{flux in the control condition}
#'   \item{lb_ctl}{lower bound of the control flux}
#'   \item{ub_ctl}{upper bound of the control flux}
#'   \item{flux_exp}{flux in the experimental condition}
#'   \item{lb_exp}{lower bound of the experimental flux}
#'   \item{ub_flux}{upper bound of the experimental flux}
#'   \item{ratio}{ratio of the experimental flux to the control flux}
#'   }
"flux_differences"


