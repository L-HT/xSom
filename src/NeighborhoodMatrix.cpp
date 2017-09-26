#include <Rcpp.h>
// [[Rcpp::depends(RcppParallel)]]
#include <RcppParallel.h>
#include <iostream>
#include <functional>
#include <cmath>

#include "NeighborhoodMatrix.h"

typedef struct xyCoordinates {
  int x_;
  int y_;
} xyCoordinates;


struct NeighborhoodTableCalculator : RcppParallel::Worker{
  const double radius_;
  RcppParallel::RVector<double> resultTable_;

  NeighborhoodTableCalculator(double radius, Rcpp::NumericVector resultTable)
      : radius_(radius), resultTable_(resultTable){

  }
  void operator()(std::size_t begin, std::size_t end){
    for (std::size_t i = begin; i < end; i++){
      resultTable_[i] = std::exp( ((double)i)*i / (-2 * radius_ * radius_));
    }
  }
};

// [[Rcpp::export]]
Rcpp::NumericVector calculateNeighborhoodTable(int somSize, double radius){
  Rcpp::NumericVector resultTable(2*somSize-1);//2*N-1, weil das erste Element die 1 selbst ist?
  NeighborhoodTableCalculator ntc(radius, resultTable);
  parallelFor(0, 2*somSize-1, ntc);
  return resultTable;
}

struct TableToCodebookMatrixConverter : RcppParallel::Worker{
  int xWinner_;
  int yWinner_;
  const int somSize_;
  const unsigned int xDim_;
  RcppParallel::RMatrix<double> resultMatrix_;
  const RcppParallel::RVector<double> lookupTable_;

  TableToCodebookMatrixConverter(int xWinner, int yWinner, int somSize, int xDim,
                                 Rcpp::NumericMatrix resultMatrix, Rcpp::NumericVector lookupTable)
    : xWinner_(xWinner), yWinner_(yWinner), somSize_(somSize), xDim_(xDim), resultMatrix_(resultMatrix),
      lookupTable_(lookupTable){

  }
  void operator()(std::size_t begin, std::size_t end){
    int distance = 0;
    double lookupResult = 0;
    for (std::size_t i = begin; i < end; i++){
      int yTemp = (int) i / somSize_;
      int xTemp = i % somSize_;
      distance = std::abs(xTemp-xWinner_)+std::abs(yTemp-yWinner_);
      lookupResult = lookupTable_[distance];
      for (std::size_t j=0; j < xDim_; j++){
        resultMatrix_.row(i)[j] = lookupResult;
      }
    }
  }
};

// [[Rcpp::export]]
Rcpp::NumericMatrix tableToCodebookMatrix(int somSize, int winnerNeuronR, int xDim, Rcpp::NumericVector lookupTable){
  Rcpp::NumericMatrix result(somSize*somSize, xDim);

  int winnerNeuron = winnerNeuronR - 1;
  int yWinner = (int) winnerNeuron / somSize;
  int xWinner = winnerNeuron % somSize;

  TableToCodebookMatrixConverter ttcmc(xWinner, yWinner, somSize, xDim, result, lookupTable);
  parallelFor(0, somSize*somSize, ttcmc);
  return result;
}

///////////////////////////////////////////
///////////////////////////////////////////


inline xyCoordinates relativeCoordinates(xyCoordinates ref, xyCoordinates given){
  xyCoordinates result;
  result.x_ = ref.x_ - given.x_;
  result.y_ = ref.y_ - given.y_;
  return result;
}

inline int xyToRowNumber(int x, int y, int xdim){
  return x + y*xdim;
}
struct NeighborhoodMatrixCalculator : RcppParallel::Worker{
  int xWinner_;
  int yWinner_;
  int somSize_;

  unsigned int largerLimit_;
  unsigned int smallerLimit_;

  double radius_;
  RcppParallel::RVector<double> result_;

  NeighborhoodMatrixCalculator(int xWinner, int yWinner, int somSize,
                               int largerLimit, int smallerLimit, double radius, Rcpp::NumericVector result)
      : xWinner_(xWinner), yWinner_(yWinner), somSize_(somSize), largerLimit_(largerLimit),
        smallerLimit_(smallerLimit), radius_(radius), result_(result){

  }

  void operator()(std::size_t begin, std::size_t end){
    int tempX = 0;
    int tempY = 0;
    double neighborhoodDistance = 0;
    for (std::size_t i = begin; i < end; i++){ //ich muss nicht durch somSize laufen, nur entlang des Quadranten
      for (std::size_t j = i; j <= largerLimit_; j++){ //das gr??ere Limit innen
        //if(i>4){std::cout<< i << " " << j << std::endl;}
        neighborhoodDistance = std::exp( (i*i + j*j) / (-2 * radius_*radius_));

        // Symmetrien ausnutzen (insgesamt 8...)
        tempX = xWinner_ + i;
        tempY = yWinner_ + j;

        if (tempX < somSize_ && tempY < somSize_){
          result_[xyToRowNumber(tempX, tempY, somSize_)] = neighborhoodDistance;
        }
        {
          tempX = xWinner_ + j;
          tempY = yWinner_ + i;
          if (tempX < somSize_ && tempY < somSize_){
            result_[xyToRowNumber(tempX, tempY, somSize_)] = neighborhoodDistance;
          }

          tempX = xWinner_ - i;
          tempY = yWinner_ - j;
          if (tempX >= 0 && tempY >= 0){
            result_[xyToRowNumber(tempX, tempY, somSize_)] = neighborhoodDistance;
          }

          tempX = xWinner_ - j;
          tempY = yWinner_ - i;
          if (tempX >= 0 && tempY >= 0){
            result_[xyToRowNumber(tempX, tempY, somSize_)] = neighborhoodDistance;
          }

          tempX = xWinner_ - i;
          tempY = yWinner_ + j;
          if (tempX >= 0 && tempY < somSize_){
            result_[xyToRowNumber(tempX, tempY, somSize_)] = neighborhoodDistance;
          }

          tempX = xWinner_ - j;
          tempY = yWinner_ + i;
          if (tempX >= 0 && tempY < somSize_){
            result_[xyToRowNumber(tempX, tempY, somSize_)] = neighborhoodDistance;
          }

          tempX = xWinner_ + i;
          tempY = yWinner_ - j;
          if (tempX < somSize_ && tempY >= 0){
            result_[xyToRowNumber(tempX, tempY, somSize_)] = neighborhoodDistance;
          }

          tempX = xWinner_ + j;
          tempY = yWinner_ - i;
          if (tempX < somSize_ && tempY >= 0){
            result_[xyToRowNumber(tempX, tempY, somSize_)] = neighborhoodDistance;
          }
        }
      }
    }
  }
};

// [[Rcpp::export]]
Rcpp::NumericVector calculateNeighborhoodMatrix(int winnerNeuronR, int somSize, double radius){
  //Neighborhood-Matrix als Tabelle (gibt (somSize x somSize)-Matrix zur?ck)
  Rcpp::NumericVector result(somSize*somSize);
  int winnerNeuron = winnerNeuronR - 1;
  int yWinner = (int) winnerNeuron / somSize;
  int xWinner = winnerNeuron % somSize;
  int xLimit = 0;
  int yLimit = 0;
  int largestQuadrant = 0; //zweistellige Zahl xy; 11: oben rechts, 01: oben links, 10: unten rechts; 00: unten links

  //std::cout << winnerNeuron << " " << yWinner << " " << xWinner << std::endl;
  //gr??ter Quadrant
  if (yWinner > somSize/2){
    largestQuadrant = 1;
  }
  if (xWinner < somSize/2){
    largestQuadrant += 10;
  }


  switch(largestQuadrant){ //am meisten Platz ist...
  case(11): //oben rechts
    xLimit = somSize - xWinner;
    yLimit = yWinner;
    break;
  case(10): //unten rechts
    xLimit = somSize - xWinner;
    yLimit = somSize - yWinner;
    break;
  case(01): //oben links
    xLimit = xWinner;
    yLimit = yWinner;
    break;
  case(00): //unten links
    xLimit = xWinner;
    yLimit = somSize - yWinner;
    break;
  }

  //bestimmt, wie der gr??te Quadrant zu durchlaufen ist (zeilenweise oder spaltenweise)
  int* largerLimit = &xLimit;
  int* smallerLimit = &yLimit;
  if (xLimit < yLimit){
    largerLimit = &yLimit;
    smallerLimit = &xLimit;
  }

 // std::cout << "nmc " << *largerLimit << " " << *smallerLimit << std::endl;
  NeighborhoodMatrixCalculator nmc(xWinner, yWinner, somSize, *largerLimit, *smallerLimit, radius, result);
 // std::cout << "nmc ";
  parallelFor(0, *smallerLimit+1, nmc);
 // std::cout << "nmc ";
  return result;
}

/*** R
# eine richtige Matrix (nicht nur die Tabelle) berechnen
calculateNeighborhoodMatrixT <- function(winner, somSize, radius){
  nbMatrix <- calculateNeighborhoodMatrix(winner, somSize, radius=1.0)
  return(matrix(nbMatrix,ncol=somSize, byrow=T))

}
*/

struct MatrixToCodebookMatrixConverter : RcppParallel::Worker{

  RcppParallel::RVector<double> inputMatrix_;
  unsigned int xDim_;
  RcppParallel::RMatrix<double> result_;

  MatrixToCodebookMatrixConverter(Rcpp::NumericVector inputMatrix, int xDim, Rcpp::NumericMatrix result)
    : inputMatrix_(inputMatrix), xDim_(xDim), result_(result){

  }

  void operator()(std::size_t begin, std::size_t end){
    for (std::size_t i = begin; i < end; i++){
      //std::transform(result_.row(i).begin(), result_.row(i).end(),
      //               result_.row(i).begin(), [](double d){ return inputMatrix_[i]; });
      for (std::size_t j = 0; j < xDim_; j++){
        //result_[i,j] = inputMatrix_[i];
        result_.row(i)[j] = inputMatrix_[i];
      }
    }
  }
};

// [[Rcpp::export]]
Rcpp::NumericMatrix matrixToCodebookMatrix(Rcpp::NumericVector matrix, int xDim){

  Rcpp::NumericMatrix result(matrix.length(), xDim);
  MatrixToCodebookMatrixConverter mtcmc(matrix, xDim, result);
  parallelFor(0, matrix.length(), mtcmc);

  return result;
}
