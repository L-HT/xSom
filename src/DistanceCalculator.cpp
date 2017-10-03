#include <Rcpp.h>
// [[Rcpp::depends(RcppParallel)]]
#include <RcppParallel.h>
//#include "DistanceCalculator.h"

struct EuclidianDistanceCalculator: RcppParallel::Worker{
  const RcppParallel::RMatrix<double> inputMatrix_;
  RcppParallel::RVector<double> resultED_;
  const RcppParallel::RVector<int> oldColumns_;

  EuclidianDistanceCalculator(const Rcpp::NumericMatrix& inputMatrix,
                              Rcpp::NumericVector& resultED, const Rcpp::LogicalVector& oldColumns)
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
void calculateEuclidianDistances(const Rcpp::NumericMatrix& deltaMatrix, const Rcpp::LogicalVector& oldColumns, Rcpp::NumericVector& resultEuclidianDistances2){
  int length = deltaMatrix.nrow();
  EuclidianDistanceCalculator edc(deltaMatrix, resultEuclidianDistances2, oldColumns);
  RcppParallel::parallelFor(0, length, edc);
}

