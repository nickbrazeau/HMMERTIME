#------------------------------------------------
#' @title polyIBD
#'
#' @description Infer blocks of identity by descent between polyclonal samples using unphased haplotype data
#'
#' @docType package
#' @name polyIBD
NULL

#------------------------------------------------
# link to Rcpp
#' @useDynLib polyIBD, .registration = TRUE
#' @importFrom Rcpp sourceCpp
NULL

#------------------------------------------------
# unload dll when package is unloaded
#' @noRd
.onUnload <- function(libpath) {
  library.dynam.unload("polyIBD", libpath)  # nocov
}
