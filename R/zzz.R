#------------------------------------------------
#' @title HMMERTIME
#'
#' @description Infer blocks of identity by descent between polyclonal samples using unphased haplotype data
#'
#' @docType package
#' @name HMMERTIME
NULL

#------------------------------------------------
# link to Rcpp
#' @useDynLib HMMERTIME, .registration = TRUE
#' @importFrom Rcpp sourceCpp
NULL

#------------------------------------------------
# unload dll when package is unloaded
#' @noRd
.onUnload <- function(libpath) {
  library.dynam.unload("HMMERTIME", libpath)  # nocov
}
