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

void calculateNeighborhoodMatrix(const int& winnerNeuronR, const int& somSize, const double& radius, Rcpp::NumericVector& resultVector);
//void calculateNeighborhoodMatrix(int winnerNeuronR, int somSize, double radius, Rcpp::NumericVector& result);

void matrixToCodebookMatrix(const Rcpp::NumericVector& matrix, Rcpp::NumericMatrix& result);

#endif
