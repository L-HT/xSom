#ifndef DISTANCECALCULATOR_HPP
#define DISTANCECALCULATOR_HPP

#include <Rcpp.h>

Rcpp::NumericVector calculateEuclidianDistances(Rcpp::NumericMatrix deltaMatrix, Rcpp::LogicalVector oldColumns);

#endif 