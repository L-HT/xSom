#' xSom: A Parallel Version of the Som Algorithm for Extendable Datasets
#'
#' A parallelized version of the self-organizing map.
#' It also supports the extension of data sets where new columns
#' are appended without changing the resulting codebook matrix in the
#' old columns.
#'
#' @docType package
#' @name xSom
NULL

#' @useDynLib xSom, .registration = TRUE
#' @importFrom Rcpp sourceCpp
#' @importFrom RcppParallel RcppParallelLibs
NULL

.onUnload <- function (libpath) {
  library.dynam.unload("xSom", libpath)
}
