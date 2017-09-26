#include <Rcpp.h>
// [[Rcpp::depends(RcppParallel)]]
#include <RcppParallel.h>
#include <iostream>

#include "Neuronsuche.h"
#include "DeltaMatrix.h"
#include "DistanceCalculator.h"

double quadraticDifference(double x, double y){
  return (x-y)*(x-y);
}

struct WinningNeuronSearch : public RcppParallel::Worker{
  int resultIndex_;
  double resultDistance_;
  double tempDistance_;
  const RcppParallel::RMatrix<double> inputMatrix_;
  const RcppParallel::RVector<double> x_;
  RcppParallel::RVector<double> distanceVector_;
  WinningNeuronSearch(const Rcpp::NumericMatrix inputMatrix, const Rcpp::NumericVector x,
                      const Rcpp::NumericVector distanceVector)
      : resultIndex_(-1), resultDistance_(-1), tempDistance_(0), inputMatrix_(inputMatrix), x_(x),
        distanceVector_(distanceVector){


  }
  void operator()(std::size_t begin, std::size_t end){
    //f?r alle diesem Arbeiter zugewiesenen Matrixreihen
    for (std::size_t i = begin; i < end; i++) {
      RcppParallel::RMatrix<double>::Row row = inputMatrix_.row(i);
      std::vector<double> tempDifferenceVector(row.length());
      std::transform(row.begin(), row.end(), x_.begin(), tempDifferenceVector.begin(), quadraticDifference);
      tempDistance_ = std::accumulate(tempDifferenceVector.begin(), tempDifferenceVector.end(), 0.0);

      if (tempDistance_ < resultDistance_ || resultIndex_ == -1){
         resultIndex_ = i+1;
         resultDistance_ = tempDistance_;
      }
    }
  }
};

//'
//' Determine the Best Matching Unit (BMU) for a data vector
//'
//' This function takes a vector \code{x} and determines the row of the weight matrix which
//' minimizes the euclidian distance.
//' @param weightMatrix The weight matrix (codebook matrix).
//' @param x The vector.
//' @param oldColumns A boolean vector specifying the columns of the data vector which were already
//' available before the data set has been extended. In other words, the only
//' these columns are used for the calculation of the euclidian distance (for the BMU)
//' while the other other columns are ignored during this step.
//'
//' @return A list containing the following elements: The number of the row which is closest
//' to \code{x}, the x-coordinate and the y-coordinate of the found neuron in the quadratic
//' SOM grid.
//'
//' @export
// [[Rcpp::export]]
Rcpp::List findWinningNeuron(Rcpp::NumericMatrix weightMatrix, Rcpp::NumericVector x, Rcpp::LogicalVector oldColumns){
  Rcpp::NumericMatrix resultDelta = calculateDelta(weightMatrix, x, true);

  Rcpp::NumericVector euclidianDistances = calculateEuclidianDistances(resultDelta, oldColumns);

  int somSize = std::sqrt(weightMatrix.nrow());
  int minIndex = findMinimumIndex(euclidianDistances);
  int neuronX = 1 + (minIndex-1) % somSize;
  int neuronY = 1 + (int) (minIndex-1) / somSize;

  // @param useRearrangedRows A boolean value specifying whether a different numbering is
  // used for the neurons in the SOM grid. See the 'Details' of the documentation for
  // the function \code{somWithMappingR} for more information.
  // if (useRearrangedRows){
  //   neuronY = somSize - neuronY + 1;
  // }

  return Rcpp::List::create(
    Rcpp::Named("neuronNumber") = minIndex,
    Rcpp::Named("neuronX") = neuronX,
    Rcpp::Named("neuronY") = neuronY
  );
}


/////////////////////////////////////////////////
/////////////////////////////////////////////////

struct MinimumVectorSearch : RcppParallel::Worker{
  const RcppParallel::RVector<double> inputVector_;
  int result_;

  MinimumVectorSearch(Rcpp::NumericVector inputVector)
    : inputVector_(inputVector), result_(0){

  }
  MinimumVectorSearch(const MinimumVectorSearch& mvs, RcppParallel::Split)
    : inputVector_(mvs.inputVector_), result_(0){

  }

  void operator()(std::size_t begin, std::size_t end){
    for (std::size_t i = begin; i < end; i++){
      if (inputVector_[i] < inputVector_[result_]){
        result_ = i;
      }
    }
  }
  void join(const MinimumVectorSearch& mvs){
    if (inputVector_[mvs.result_] < inputVector_[result_]){
      result_ = mvs.result_;
    }
  }
};

int findMinimumIndex(Rcpp::NumericVector input){
  MinimumVectorSearch mvs(input);
  RcppParallel::parallelReduce(0, input.length(), mvs);
  return mvs.result_+1; //+1 wegen R
}

int findMinimumIndexWithR(Rcpp::NumericVector input){
  return Rcpp::which_min(input)+1;
}
