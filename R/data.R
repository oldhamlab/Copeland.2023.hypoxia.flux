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
