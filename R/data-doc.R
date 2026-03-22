#' Reference gene expression matrix
#'
#' Packaged reference/core expression matrix with genes in rows and samples
#' in columns.
#'
#' @format A numeric matrix.
"exp_core_g"

#' Reference sample annotations
#'
#' Packaged sample annotation data frame for the reference/core cohort.
#' Row names match `colnames(exp_core_g)` and the data frame contains a `CTS`
#' column.
#'
#' @format A data frame.
"core_samples"
