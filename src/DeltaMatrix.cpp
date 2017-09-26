/*
 * Berechne die Delta-Matrix parallelisiert mit den Zeilen
 *
 * Delta_i = Goal - W_i
 */

#include <Rcpp.h>
// [[Rcpp::depends(RcppParallel)]]
#include <RcppParallel.h>
#include <iostream>
#include <functional>

#include "DeltaMatrix.h"

double addSquares(double x, double y){
  return x + y*y;
}


bool isInvalidNumber(double x){
  bool a = Rcpp::traits::is_na<REALSXP>(x);
  bool b = Rcpp::traits::is_nan<REALSXP>(x);
  bool c = Rcpp::traits::is_infinite<REALSXP>(x);
  return a || b || c;
}

//das urspr?ngliche Soft-Minus (x-Vektor - Gewichtsvektor)
double softMinus(double x, double y){
  double result = 0;
  if (!isInvalidNumber(x) && !isInvalidNumber(y)){
    result = x - y;
  }
  return result;
}

//neues Soft-Minus (wegen der Toleranz gegen?ber neuen Spalten
//muss man als Ausgangspunkt nicht x, sondern w nehmen und davon
//die L?nge von w), berechnet also -(w-x)
double softMinus2(double x, double y){
  double result = 0;
  if (!isInvalidNumber(x) && !isInvalidNumber(y)){
    result = y - x;
  }
  return result;
}

//normales Minus
double nonSoftMinus(double x, double y){
  return y - x;

}

struct DeltaMatrixCalculator: RcppParallel::Worker{
  const RcppParallel::RMatrix<double> inputMatrix_;
  const RcppParallel::RVector<double> inputVector_;
  RcppParallel::RMatrix<double> resultDelta_;

  DeltaMatrixCalculator(const Rcpp::NumericMatrix inputMatrix,
                        const Rcpp::NumericVector inputVector,
                        const Rcpp::NumericMatrix resultDelta)
    : inputMatrix_(inputMatrix), inputVector_(inputVector), resultDelta_(resultDelta){
  }

  void operator()(std::size_t begin, std::size_t end){
    for (std::size_t i = begin; i < end; i++){
      RcppParallel::RMatrix<double>::Row row = inputMatrix_.row(i);

      //neuer Code: Rechnet nun nur bis zur Breite der Gewichtsmatrix, Rest wird ignoriert
      std::transform(
        row.begin(), row.end(),
        inputVector_.begin(),
        resultDelta_.row(i).begin(),
        softMinus2
      );

    }
  }
};

struct DeltaMatrixCalculatorNoNA: RcppParallel::Worker{
  const RcppParallel::RMatrix<double> inputMatrix_;
  const RcppParallel::RVector<double> inputVector_;
  RcppParallel::RMatrix<double> resultDelta_;

  DeltaMatrixCalculatorNoNA(const Rcpp::NumericMatrix inputMatrix,
                        const Rcpp::NumericVector inputVector,
                        const Rcpp::NumericMatrix resultDelta)
    : inputMatrix_(inputMatrix), inputVector_(inputVector), resultDelta_(resultDelta){
  }

  void operator()(std::size_t begin, std::size_t end){
    for (std::size_t i = begin; i < end; i++){
      RcppParallel::RMatrix<double>::Row row = inputMatrix_.row(i);

      //neuer Code: Rechnet nun nur bis zur Breite der Gewichtsmatrix, Rest wird ignoriert
      // std::transform(
      //   row.begin(), row.end(),
      //   inputVector_.begin(),
      //   resultDelta_.row(i).begin(),
      //   nonSoftMinus
      // );
      std::transform(
        inputVector_.begin(), inputVector_.end(),
        row.begin(),
        resultDelta_.row(i).begin(),
        std::minus<double>()
      );

    }
  }
};

/*** R
replaceNAByMeanInVector <- function(x){
  result <- x
  result[is.na(result)] = mean(x, na.rm=T)
  return(result)
}

replaceNAByMean <- function(x, byRow=F){
  if (byRow == F){
    result <- apply(a, MARGIN=2, replaceNAByMeanInVector)
  }else{
    result <- t(apply(a, MARGIN=2, replaceNAByMeanInVector))
  }
  return(result)
}
*/

// [[Rcpp::export]]
Rcpp::NumericMatrix calculateDelta(Rcpp::NumericMatrix inputMatrix, Rcpp::NumericVector inputVector, bool naExist){
  Rcpp::NumericMatrix resultDelta(inputMatrix.nrow(), inputMatrix.ncol());
  Rcpp::NumericVector resultEuclidianDistances2(inputMatrix.nrow());
  if (naExist){
    DeltaMatrixCalculator dmc(inputMatrix, inputVector, resultDelta);
    RcppParallel::parallelFor(0, inputMatrix.nrow(), dmc);
  }else{
    DeltaMatrixCalculatorNoNA dmc(inputMatrix, inputVector, resultDelta);
    RcppParallel::parallelFor(0, inputMatrix.nrow(), dmc);
  }
  return resultDelta;
}

