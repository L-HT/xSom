#ifndef NEIGHBORHOODMATRIX_HPP
#define NEIGHBORHOODMATRIX_HPP

#include <Rcpp.h>
// [[Rcpp::depends(RcppParallel)]]
#include <RcppParallel.h>
#include <iostream>
#include <functional>
#include <cmath>


Rcpp::NumericVector calculateNeighborhoodTable(int somSize, double radius);
Rcpp::NumericMatrix tableToCodebookMatrix(int somSize, int winnerNeuronR, int xDim, Rcpp::NumericVector lookupTable);
Rcpp::NumericVector calculateNeighborhoodMatrix(int winnerNeuronR, int somSize, double radius);
Rcpp::NumericMatrix matrixToCodebookMatrix(Rcpp::NumericVector matrix, int xDim, Rcpp::NumericMatrix& result);

#endif
