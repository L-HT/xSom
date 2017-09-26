#ifndef NEURONSUCHE_CPP
#define NEURONSUCHE_CPP

#include <Rcpp.h>
// [[Rcpp::depends(RcppParallel)]]
#include <RcppParallel.h>
#include <iostream>


Rcpp::List findWinningNeuron(Rcpp::NumericMatrix weightMatrix, Rcpp::NumericVector x, Rcpp::LogicalVector oldColumns);

int findMinimumIndex(Rcpp::NumericVector input);

int findMinimumIndexWithR(Rcpp::NumericVector input);

#endif
