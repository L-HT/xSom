#ifndef SIMPLEMATRIXOPERATIONS_HPP
#define SIMPLEMATRIXOPERATIONS_HPP

#include <Rcpp.h>
// [[Rcpp::depends(RcppParallel)]]
#include <RcppParallel.h>
#include <iostream>

Rcpp::NumericMatrix matrixTimesScalar(Rcpp::NumericMatrix m, double a);

Rcpp::NumericMatrix elementwiseMultiplication(Rcpp::NumericMatrix m1, Rcpp::NumericMatrix m2);
Rcpp::NumericMatrix elementwiseAddition(Rcpp::NumericMatrix m1, Rcpp::NumericMatrix m2);

#endif

