#ifndef SOMFUNCTIONS_HPP
#define SOMFUNCTIONS_HPP

#include <Rcpp.h>
// [[Rcpp::depends(RcppParallel)]]
#include <RcppParallel.h>


Rcpp::NumericMatrix somTwoStep(Rcpp::NumericMatrix dataSet, Rcpp::NumericMatrix weightMatrix, 
                               Rcpp::LogicalVector oldColumns,
                               Rcpp::IntegerVector rlen, 
                               Rcpp::NumericVector initLearnRate, 
                               Rcpp::NumericVector initRadius,
                               double radiusReduction, double learnRateReduction,
                               int normType, int sampling, bool naExist,
                               bool updateParametersPerEpoch);
                               
Rcpp::List somWithMapping(Rcpp::NumericMatrix dataSet, Rcpp::NumericMatrix weightMatrix, 
                         Rcpp::LogicalVector oldColumns,
                         Rcpp::IntegerVector rlen, 
                         Rcpp::NumericVector initLearnRate, 
                         Rcpp::NumericVector initRadius,
                         double radiusReduction, double learnRateReduction,
                         int normType, int sampling, bool naExist,
                         bool updateParametersPerEpoch);
  
#endif