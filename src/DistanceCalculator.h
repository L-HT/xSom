#ifndef DISTANCECALCULATOR_HPP
#define DISTANCECALCULATOR_HPP

#include <Rcpp.h>

void calculateEuclidianDistances(const Rcpp::NumericMatrix& deltaMatrix, const Rcpp::LogicalVector& oldColumns, Rcpp::NumericVector& resultEuclidianDistances2);

#endif
