/*
 * Berechne daraus einen Vektor euclidianDistances2 mit
 *
 * euclidianDistances2_i = Delta_i^T * Delta_i
 *
 */

#include <Rcpp.h>
// [[Rcpp::depends(RcppParallel)]]
#include <RcppParallel.h>
#include <iostream>
#include <functional>

#include "DistanceCalculator.h"

// l?uft ?ber Zeilen der Matrix und bildet da das Skalarprodukt mit sich selbst
struct EuclidianDistanceCalculator: RcppParallel::Worker{
  const RcppParallel::RMatrix<double> inputMatrix_;
  RcppParallel::RVector<double> resultED_;
  RcppParallel::RVector<int> oldColumns_;


  EuclidianDistanceCalculator(const Rcpp::NumericMatrix inputMatrix,
                              const Rcpp::NumericVector resultED, const Rcpp::LogicalVector oldColumns)
    : inputMatrix_(inputMatrix), resultED_(resultED), oldColumns_(oldColumns){

  }

  void operator()(std::size_t begin, std::size_t end){
    for (std::size_t i = begin; i < end; i++){
      RcppParallel::RMatrix<double>::Row row = inputMatrix_.row(i);

      double partialEuclidianDistance = 0;
      for (unsigned int j = 0; j < row.length(); j++){
        if (oldColumns_[j]){
          partialEuclidianDistance += row[j] * row[j];
        }
      }
      resultED_[i] = partialEuclidianDistance;
    }
  }
};

// [[Rcpp::export]]
Rcpp::NumericVector calculateEuclidianDistances(Rcpp::NumericMatrix deltaMatrix, Rcpp::LogicalVector oldColumns){

  int length = deltaMatrix.nrow();
  Rcpp::NumericVector resultEuclidianDistances2(length);

  EuclidianDistanceCalculator edc(deltaMatrix, resultEuclidianDistances2, oldColumns);
  RcppParallel::parallelFor(0, length, edc);

  // return resultDelta;
  return resultEuclidianDistances2;
}

