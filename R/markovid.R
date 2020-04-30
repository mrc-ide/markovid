#------------------------------------------------
#' @title MCMC inference of durations and flows from UK hospital data
#'
#' @description Re-purpose of drjacoby package, tailored to COVID application.
#'
#' @docType package
#' @name markovid
NULL

#------------------------------------------------
#' @useDynLib markovid, .registration = TRUE
#' @importFrom Rcpp sourceCpp
#' @importFrom magrittr %>%
NULL

#------------------------------------------------
# unload DLL when package is unloaded
#' @noRd
.onUnload <- function(libpath) {
  library.dynam.unload("markovid", libpath)  # nocov
}
