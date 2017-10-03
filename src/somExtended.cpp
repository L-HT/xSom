#include <Rcpp.h>
// [[Rcpp::depends(RcppParallel)]]
#include <RcppParallel.h>
#include <iostream>

#include "somExtended.h"
#include "DeltaMatrix.h"
#include "DistanceCalculator.h"
#include "Neuronsuche.h"
#include "NeighborhoodMatrix.h"
#include "SimpleMatrixOperations.h"

/*
 * allokiere eine Matrix so groß wie W
 * In die wird die Delta-Matrix geschrieben.
 * genause eine Matrix so groß wie W: da rein kommt die matrixToCodebookMatrix
 *
 * Dann: calculateEuclidianDistances könnte const ref nehmen
 */

// [[Rcpp::export]]
Rcpp::NumericMatrix learnCyclesExtended(Rcpp::NumericMatrix dataSet, Rcpp::NumericMatrix weightMatrix,
                                Rcpp::LogicalVector oldColumns,
                                unsigned int cycles = 1, double initLearnRate = 0.01,
                                double learnRateReduction = 0.0, double initRadius = 1.0, double radiusReduction = -1.0,
                                int normType = 2, int sampling = 1, bool naExist = true, bool updateParametersPerEpoch = true){

  int newColumns = weightMatrix.ncol();
  int somSize = std::sqrt(weightMatrix.nrow());
  int cycleIntern = 0;
  int maxCycleIntern = cycles * dataSet.nrow();
  double currentRadius = initRadius;
  double learnRate = initLearnRate;
  double inverseLearnRateC = -learnRateReduction;

  Rcpp::NumericMatrix result(weightMatrix);
  Rcpp::IntegerVector data_order;
  Rcpp::NumericVector x;
  Rcpp::NumericMatrix resultDelta(weightMatrix.nrow(), weightMatrix.ncol());
  Rcpp::NumericMatrix neighborhoodMatrix(weightMatrix.nrow(), weightMatrix.ncol());;
  Rcpp::NumericVector euclidianDistances2(weightMatrix.nrow());
  Rcpp::NumericVector nbMatrix(weightMatrix.nrow());

  int progressStep = 1;

  // if (cycles >= 100){
  //   progressStep = (int) (cycles / 100);
  //   Rcpp::Rcout << "1%                                                                                              100%" << std::endl;
  // }else{
  //   if (!updateParametersPerEpoch && maxCycleIntern > 100){
  //     progressStep = (int) (maxCycleIntern / 100);
  //     Rcpp::Rcout << "1%                                                                                              100%" << std::endl;
  //   }
  // }
  if (maxCycleIntern > 100){
    progressStep = (int) (maxCycleIntern / 100);
    Rcpp::Rcout << "1%                                                                                              100%" << std::endl;
  }
  for (unsigned int cycle = 0; cycle < cycles; cycle++){
    // if (cycle % progressStep == 0 && cycles >= 100 && updateParametersPerEpoch){
    //   Rcpp::Rcout << ".";
    // }
    Rcpp::IntegerVector data_order;
    if (sampling == 0){
      data_order = Rcpp::seq_len(dataSet.nrow());
    }else{
      if (sampling == 1){
        data_order = Rcpp::sample(dataSet.nrow(), dataSet.nrow());
      }else{
        data_order = Rcpp::sample(dataSet.nrow(), sampling, true);
      }
    }
    for (Rcpp::IntegerVector::iterator it = data_order.begin(); it < data_order.end(); it++){
      cycleIntern++;
      if (!updateParametersPerEpoch){
        if (radiusReduction < 0){
          currentRadius = 1+(initRadius-1)*(maxCycleIntern-cycleIntern)/(maxCycleIntern);
        }else{
          currentRadius -= radiusReduction;
        }
        if (learnRateReduction == 0){
          // Rcpp::Rcout << currentRadius << "--" << maxCycleIntern << "--" << cycleIntern << "--" << learnRate << std::endl;

          learnRate = initLearnRate*(1.0-((double)cycleIntern)/maxCycleIntern);
        }else{
          if (learnRateReduction > 0){
            learnRate -= learnRateReduction;
          }else{
            learnRate = initLearnRate * inverseLearnRateC / (inverseLearnRateC + cycleIntern);
          }
        }
      }
      if (cycleIntern % progressStep == 0 && maxCycleIntern >= 100){
        Rcpp::Rcout << ".";
      }
      int chosenIndex = *it - 1;
      x = dataSet.row(chosenIndex);
      try{
        calculateDelta(result, x, naExist, resultDelta);
        calculateEuclidianDistances(resultDelta, oldColumns, euclidianDistances2);
        int winner = findMinimumIndex(euclidianDistances2);

        if (normType == 1){
          Rcpp::NumericVector nbTable = calculateNeighborhoodTable(somSize, currentRadius);
          neighborhoodMatrix = tableToCodebookMatrix(somSize, winner, newColumns, nbTable);
        }else{
          //Rcpp::NumericVector nbMatrix = calculateNeighborhoodMatrix(winner, somSize, currentRadius);
          calculateNeighborhoodMatrix(winner, somSize, currentRadius, nbMatrix);
          //nbMatrix.attr("dim") = Rcpp::Dimension(somSize, somSize);
          matrixToCodebookMatrix(nbMatrix, neighborhoodMatrix); //ehemalig noch newColumns als 2. Argument

        }
        result = elementwiseAddition(
          result,
          matrixTimesScalar(
            elementwiseMultiplication(
              resultDelta,
              neighborhoodMatrix
            ),
            learnRate
          )
        );
      }catch(std::exception &ex){
        forward_exception_to_r(ex);
      }
    }
    if (updateParametersPerEpoch){
      currentRadius = 1+(initRadius-1)*(cycles-(cycle+1.0))/(cycles);
      if (learnRateReduction >= 0){
        learnRate = initLearnRate*(1.0-(cycle+1.0)/cycles);
      }else{
        learnRate = initLearnRate * inverseLearnRateC / (inverseLearnRateC + (cycle+1));
      }
    }
    Rcpp::checkUserInterrupt();
  }
  Rcpp::Rcout << std::endl;
  return result;
}
