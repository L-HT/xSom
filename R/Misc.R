
#' Set Number of Threads
#'
#' This function uses the function \code{setThreadOptions()} from the RcppParallel package.
#'
#' See also the function \code{RcppParallel::setThreadOptions()} if you
#' want to change the stack size für the worker threads.
#'
#' @param n Number of threads.
#' @export
#'
setNumberOfThreads <- function(n){
  RcppParallel::setThreadOptions(numThreads=n)
}

#' Get Default Number of Threads
#'
#' The default value depends on the system.
#'
#' @return The default number of threads.
#' @export
getDefaultNumberOfThreads <- function(){
  return(RcppParallel::defaultNumThreads())
}


