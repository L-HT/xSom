#ifndef DELTAMATRIX_HPP
#define DELTAMATRIX_HPP

#include <Rcpp.h>

void calculateDelta(const Rcpp::NumericMatrix& inputMatrix, const Rcpp::NumericVector& inputVector, bool naExist, Rcpp::NumericMatrix& resultDelta);

bool isInvalidNumber(double x);

#endif
