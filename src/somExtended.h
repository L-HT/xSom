#ifndef SOMEXTENDED_HPP
#define SOMEXTENDED_HPP

#include <Rcpp.h>
// [[Rcpp::depends(RcppParallel)]]
#include <RcppParallel.h>


Rcpp::NumericMatrix learnCyclesExtended(Rcpp::NumericMatrix dataSet, Rcpp::NumericMatrix weightMatrix, 
                                        Rcpp::LogicalVector oldColumns,
                                        unsigned int cycles, double initLearnRate, 
                                        double learnRateReduction, double initRadius, double radiusReduction,
                                        int normType, int sampling, bool naExist, bool updateParametersPerEpoch);
#endif