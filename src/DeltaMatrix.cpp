#include <Rcpp.h>
// [[Rcpp::depends(RcppParallel)]]
#include <RcppParallel.h>

//#include "DeltaMatrix.h"

bool isInvalidNumber(const double& x){
  bool a = Rcpp::traits::is_na<REALSXP>(x);
  bool b = Rcpp::traits::is_nan<REALSXP>(x);
  bool c = Rcpp::traits::is_infinite<REALSXP>(x);
  return a || b || c;
}

double softMinus(const double& x, const double& y){
  double result = 0;
  if (!isInvalidNumber(x) && !isInvalidNumber(y)){
    result = x - y;
  }
  return result;
}

//neues Soft-Minus (wegen der Toleranz gegen?ber neuen Spalten
//muss man als Ausgangspunkt nicht x, sondern w nehmen und davon
//die L?nge von w), berechnet also -(w-x)
double softMinus2(const double& x, const double& y){
  double result = 0;
  if (!isInvalidNumber(x) && !isInvalidNumber(y)){
    result = y - x;
  }
  return result;
}


struct DeltaMatrixCalculator: RcppParallel::Worker{
  const RcppParallel::RMatrix<double> inputMatrix_;
  const RcppParallel::RVector<double> inputVector_;
  RcppParallel::RMatrix<double> resultDelta_;

  DeltaMatrixCalculator(const Rcpp::NumericMatrix& inputMatrix,
                        const Rcpp::NumericVector& inputVector,
                        Rcpp::NumericMatrix& resultDelta)
    : inputMatrix_(inputMatrix), inputVector_(inputVector), resultDelta_(resultDelta){

  }

  void operator()(std::size_t begin, std::size_t end){
    for (std::size_t i = begin; i < end; i++){
      RcppParallel::RMatrix<double>::Row row = inputMatrix_.row(i);
      // RcppParallel::RMatrix<double>::Row resultRow = resultDelta_.row(i);
      // for (unsigned int k = 0; k < row.length(); k++){
      //   resultRow[k] = softMinus2(row[k], inputVector_[k]);
      // }
      //neuer Code: Rechnet nun nur bis zur Breite der Gewichtsmatrix, Rest wird ignoriert
      std::transform(
        row.begin(), row.end(),
        //inputMatrix_.row(i)).begin(), inputMatrix_.row(i).end(),
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

  DeltaMatrixCalculatorNoNA(const Rcpp::NumericMatrix& inputMatrix,
                        const Rcpp::NumericVector& inputVector,
                        Rcpp::NumericMatrix& resultDelta)
    : inputMatrix_(inputMatrix), inputVector_(inputVector), resultDelta_(resultDelta){

  }

  void operator()(std::size_t begin, std::size_t end){
    for (std::size_t i = begin; i < end; i++){
      RcppParallel::RMatrix<double>::Row row = inputMatrix_.row(i);
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
void calculateDelta(const Rcpp::NumericMatrix& inputMatrix, const Rcpp::NumericVector& inputVector, const bool naExist, Rcpp::NumericMatrix& resultDelta){
  if (naExist){
    DeltaMatrixCalculator dmc(inputMatrix, inputVector, resultDelta);
    RcppParallel::parallelFor(0, inputMatrix.nrow(), dmc);
  }else{
    DeltaMatrixCalculatorNoNA dmc(inputMatrix, inputVector, resultDelta);
    RcppParallel::parallelFor(0, inputMatrix.nrow(), dmc);
  }
}

