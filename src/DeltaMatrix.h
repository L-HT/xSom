#ifndef DELTAMATRIX_HPP
#define DELTAMATRIX_HPP

#include <Rcpp.h>

Rcpp::NumericMatrix calculateDelta(Rcpp::NumericMatrix inputMatrix, Rcpp::NumericVector inputVector, bool naExist);

bool isInvalidNumber(double x);

#endif 