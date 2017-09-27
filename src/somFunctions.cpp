#include <Rcpp.h>
// [[Rcpp::depends(RcppParallel)]]
#include <RcppParallel.h>
#include <iostream>
#include <functional>

#include "somExtended.h"
#include "DeltaMatrix.h"
#include "somFunctions.h"
#include "Neuronsuche.h"

Rcpp::NumericMatrix somTwoStep(Rcpp::NumericMatrix dataSet, Rcpp::NumericMatrix weightMatrix,
                                Rcpp::LogicalVector oldColumns,
                                Rcpp::IntegerVector rlen = Rcpp::IntegerVector::create(0),
                                Rcpp::NumericVector initLearnRate = Rcpp::NumericVector::create(0),
                                Rcpp::NumericVector initRadius = Rcpp::NumericVector::create(0.0),
                                double radiusReduction = -1.0, double learnRateReduction = 0.0,
                                int normType = 2, int sampling = 0, bool naExist = true,
                                bool updateParametersPerEpoch = true){

  Rcpp::NumericMatrix result(weightMatrix);
  int somSize = std::sqrt(weightMatrix.nrow());
  if (rlen[0] == 0){
    rlen = Rcpp::IntegerVector::create(2, 10);
  }else{
    if (rlen.length() == 1){
      rlen = Rcpp::IntegerVector::create(rlen[0], rlen[0]*10);
    }
  }
  if (initLearnRate[0] == 0.0){
    initLearnRate = Rcpp::NumericVector::create(0.05, 0.02);
  }else{
    if (initLearnRate.length() == 1){
      initLearnRate = Rcpp::NumericVector::create(initLearnRate[0], initLearnRate[0]/2.0);
    }
  }
  if (initRadius[0] == 0.0){
    initRadius = Rcpp::NumericVector::create(somSize, (somSize < 3 ? somSize : 3) );
  }else{
    if (initRadius.length() == 1){
      initRadius = Rcpp::NumericVector::create(initRadius[0],  (initRadius[0] > 2 ? initRadius[0] : 2.0 ) );
    }
  }
  Rcpp::Rcout << "Erste Trainingsstufe (rlen = " << rlen[0] << " Zyklen, initLearnRate = " << initLearnRate[0] << ", initRadius = " <<
    initRadius[0] << ")" << std::endl;
  result = learnCyclesExtended(dataSet, weightMatrix, oldColumns, rlen[0],
                              initLearnRate[0], learnRateReduction, initRadius[0],
                              radiusReduction, normType, sampling, naExist,
                              updateParametersPerEpoch
                              );

  if (rlen[1] != 0){
    Rcpp::Rcout << "Zweite Trainingsstufe (rlen = " << rlen[1] << " Zyklen, initLearnRate = " << initLearnRate[1] << ", initRadius = " <<
      initRadius[1] << ")" << std::endl;
    result = learnCyclesExtended(dataSet, result, oldColumns, rlen[1],
                                 initLearnRate[1], learnRateReduction, initRadius[1],
                                 radiusReduction, normType, sampling, naExist,
                                 updateParametersPerEpoch
                                 );
  }
  return result;
}

//' A SOM algorithm that handles NA values.
//'
//' If the argument \code{naExist} is set to \code{FALSE}, the first 100 lines of the data set will be
//' scanned for NA values. Depending on whether NA values are found slightly different
//' SOM algorithms are executed:
//' If NA values exist, the algorithm will ignore these values during the calculation of the
//' Best Matching Unit (BMU). If no NA values exist, this step is omitted which improves
//' the runtime of the algorithm.
//'
//' @param dataSet The data set.
//' @param weightMatrix The initial weight matrix (also known as codebook matrix). A possible
//' initialization is provided with the function \code{som.init.extended}. Note that the
//' weight matrix must have N rows where N is a squared number because this SOM algorithm
//' arranges the N codebook neurons in a two-dimensional quadratic grid.
//' @param oldColumns A boolean vector specifying the columns of the data set which were already
//' available before the data set has been extended. In other words, the SOM algorithm only
//' uses these columns for the calculation of Best Matching Units while the other other columns
//' are ignored during this step, but still modified during the update of the weight matrix
//' (codebook matrix).
//' @param rlen A vector with two integers specifying how many cycles are used for the two training
//' steps. Default value is \code{(2,10)}. If only one value x is provided, it will be
//' interpreted as \code{(x, 10*x)}. This is based on the \code{som} package written by Jun Yan.
//' @param initLearnRate A vector with two numbers specifying the initial learn rate for the
//' two training steps. The default is \code{(0.05,0.02)}. If only one value \code{x} given, it will be
//' interpreted as \code{(x, x/2.0)}.
//' This is based on the \code{som} package written by Jun Yan.
//' @param initRadius A vector with two numbers specifying the initial radii for the
//' gaussian neighborhood function during the two training steps.
//' The default is \code{(somSize,min(somSize,3))}, where somSize refers to the size of the
//' quadratic grid which is used for the codebook neurons (see the argument \code{weightMatrix}.
//' If only one value \code{x} given, it will be interpreted as \code{(x, max(x,2))}.
//' This is based on the \code{som} package written by Jun Yan.
//' @param radiusReduction A number specifying how the radius for the gaussian neighborhood function
//' is reduced during each parameter update (see also \code{updateParametersPerEpoch}).
//' @param learnRateReduction A number specifying how the radius for the gaussian neighborhood function
//' is reduced during each parameter update (see also \code{updateParametersPerEpoch}).
//' @param normType This must be either 1 or 2. It specifies how the distance between neurons
//' in the quadratic grid is calculated. The value 1 uses the 1-norm (absolute value) and leads
//' to rectangular structures in the resulting SOM while the value 2 uses the 2-norm (root of the
//' sum of squares, also known as euclidian norm) which leads to circular structures.
//' @param sampling A non-negative integer specifying how the data set is run through during
//' a training cycle. The value 0 makes the algorithm run through the data set without any
//' changes. The value 1 makes it use a random order while any value K greater than 1
//' makes it randomly select K rows from the data sets (there it is possible that the
//' same row is selected multiple times)
//' @param naExist A boolean value indicating whether NAs are present. If the argument
//' \code{naExist} is set to \code{FALSE}, the first 100 lines of the data set will be
//' scanned for NA values. Depending on whether NA values are found slightly different
//' SOM algorithms are executed:
//' If NA values exist, the algorithm will ignore these values during the calculation of the
//' Best Matching Unit (BMU). If no NA values exist, this step is omitted which improves
//' the runtime of the algorithm.
//' @param updateParametersPerEpoch A boolean value specifying when the learning rate and the
//' radius are updated during the SOM training. If it is set to \code{TRUE}, the update is done
//' after the algorithm went through the data set (or the sampled rows if \code{sampling} > 1).
//' The \code{som} algorithm written bei Jun Yan updates the learn rate and radius after every row.
//' This behavior is replicated by setting \code{updateParametersPerEpoch} to \code{FALSE}.
//' @return A weight matrix (codebook matrix) after the training has been completed.
//' @export
//' @examples
//' # generate data and initialize weight matrix
//' dataSet <- matrix(as.numeric(1:400),ncol=2)
//' weightMatrix <- som.init.extended(dataSet, somSize=2, oldColumns=c(TRUE,TRUE))
//'
//' # apply the algorithm
//' result <- somWithMapping(dataSet, weightMatrix, oldColumns=c(TRUE,TRUE))
// [[Rcpp::export]]
Rcpp::NumericMatrix somCheckNa(Rcpp::NumericMatrix dataSet, Rcpp::NumericMatrix weightMatrix,
                               Rcpp::LogicalVector oldColumns,
                               Rcpp::IntegerVector rlen = Rcpp::IntegerVector::create(0),
                               Rcpp::NumericVector initLearnRate = Rcpp::NumericVector::create(0),
                               Rcpp::NumericVector initRadius = Rcpp::NumericVector::create(0.0),
                               double radiusReduction = -1.0, double learnRateReduction = 0.0,
                               int normType = 2, int sampling = 0, bool naExist = false,
                               bool updateParametersPerEpoch = true){

  int somSize = std::sqrt(weightMatrix.nrow());
  Rcpp::NumericMatrix result(weightMatrix);
  if (rlen[0] == 0){
    rlen = Rcpp::IntegerVector::create(2, 10);
  }else{
    if (rlen.length() == 1){
      rlen = Rcpp::IntegerVector::create(rlen[0], rlen[0]*10);
    }
  }
  if (initLearnRate[0] == 0.0){
    initLearnRate = Rcpp::NumericVector::create(0.05, 0.02);
  }else{
    if (initLearnRate.length() == 1){
      initLearnRate = Rcpp::NumericVector::create(initLearnRate[0], initLearnRate[0]/2.0);
    }
  }
  if (initRadius[0] == 0.0){
    initRadius = Rcpp::NumericVector::create(somSize, (somSize < 3 ? somSize : 3) );
  }else{
    if (initRadius.length() == 1){
      initRadius = Rcpp::NumericVector::create(initRadius[0],  (initRadius[0] > 2 ? initRadius[0] : 2.0 ) );
    }
  }
  if (naExist){
    result = somTwoStep(dataSet, weightMatrix,  oldColumns,rlen, initLearnRate, initRadius,radiusReduction,
                         learnRateReduction, normType, sampling, naExist, updateParametersPerEpoch);
  }else{
    int n = dataSet.nrow();
    if (n > 100){
      n = 100;
    }
    Rcpp::NumericMatrix subset = dataSet(Rcpp::Range(0,n-1), Rcpp::_);
    int subRow = subset.nrow();
    int subCol = subset.ncol();
    bool naFound = false;
    for (int i = 0; i < subRow; i++){
      Rcpp::NumericMatrix::Row tempRow = subset.row(i);
      for (int k = 0; k < subCol; k++){
        if (isInvalidNumber(tempRow[k])){
          naFound = true;
        }
      }
    }
    if (naFound){
      Rcpp::Rcout << " -- NAs in den ersten 100 Zeilen erkannt. Es wird Training mit NAs genutzt. -- " << std::endl;
    }
    result = somTwoStep(dataSet, weightMatrix,  oldColumns,rlen, initLearnRate, initRadius,radiusReduction,
                        learnRateReduction, normType, sampling, naFound, updateParametersPerEpoch);
  }
  return result;
}

//' A SOM algorithm that handles NA values.
//'
//' This function is identical to \code{somCheckNa}, but returns a list containing the following
//' elements: (insert here)
//'
//' ...
//'
//' @inheritParams somCheckNa
//' @return A list described above.
//' //@export
//' @examples
//' # generate data and initialize weight matrix
//' dataSet <- matrix(as.numeric(1:400),ncol=2)
//' weightMatrix <- som.init.extended(dataSet, somSize=2, oldColumns=c(TRUE,TRUE))
//'
//' # apply the algorithm
//' result <- somWithMapping(dataSet, weightMatrix, oldColumns=c(TRUE,TRUE))
// [[Rcpp::export]]
Rcpp::List somWithMapping(Rcpp::NumericMatrix dataSet, Rcpp::NumericMatrix weightMatrix,
                               Rcpp::LogicalVector oldColumns,
                               Rcpp::IntegerVector rlen = Rcpp::IntegerVector::create(0),
                               Rcpp::NumericVector initLearnRate = Rcpp::NumericVector::create(0),
                               Rcpp::NumericVector initRadius = Rcpp::NumericVector::create(0.0),
                               double radiusReduction = -1.0, double learnRateReduction = 0.0,
                               int normType = 2, int sampling = 0, bool naExist = true,
                               bool updateParametersPerEpoch = true){

  Rcpp::NumericMatrix result = somCheckNa(dataSet, weightMatrix, oldColumns, rlen, initLearnRate, initRadius,
                                          radiusReduction, learnRateReduction, normType, sampling,
                                          naExist, updateParametersPerEpoch);
  Rcpp::NumericVector dataCounter(dataSet.nrow());
  Rcpp::NumericVector neuronNumberVector(dataSet.nrow());
  Rcpp::NumericVector neuronXVector(dataSet.nrow());
  Rcpp::NumericVector neuronYVector(dataSet.nrow());
  Rcpp::List neuronCoordinates;

  for (int i=0; i<dataSet.nrow(); i++){
    dataCounter[i] = i+1;
    neuronCoordinates = findWinningNeuron(result, dataSet.row(i), oldColumns);
    neuronNumberVector[i] = neuronCoordinates["neuronNumber"];
    neuronXVector[i] = neuronCoordinates["neuronX"];
    neuronYVector[i] = neuronCoordinates["neuronY"];
  }

  return Rcpp::List::create(
    Rcpp::Named("codebook") = result,
    Rcpp::Named("mapping") = Rcpp::DataFrame::create(
      Rcpp::Named("n") = dataCounter,
      Rcpp::Named("neuronNumber") = neuronNumberVector,
      Rcpp::Named("neuronX") = neuronXVector,
      Rcpp::Named("neuronYVector") = neuronYVector
    )
  );
}

// Umsortieren der Zeilen, damit es zu anderen Nummerierungen der Neuronen im Gitter passt
// aktuell ist
//   1 2 3
//   4 5 6
//   7 8 9
// gewollt ist:
//   7 8 9
//   4 5 6
//   1 2 3
// [[Rcpp::export]]
Rcpp::List rearrangeNeuronCoordinates(Rcpp::List neuronCoordinates, int somSize){
  int newY = somSize - (int) neuronCoordinates["neuronY"] + 1;
  return Rcpp::List::create(
    Rcpp::Named("neuronNumber") = (newY - 1) * somSize + (int) neuronCoordinates["neuronX"],
    Rcpp::Named("neuronX") = neuronCoordinates["neuronX"],
    Rcpp::Named("neuronY") = newY
  );
}

