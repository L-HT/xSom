learnCyclesExtendedR <- function(dataSet, weightMatrix, oldColumns,cycles = 1, initLearnRate = 0.01,
                                 learnRateReduction = 0.0,  initRadius = 1.0,  radiusReduction = -1.0,
                                 normType = 2, sampling = 1, naExist = T, updateParametersPerEpoch = T){

  nrowDataSet <- nrow(dataSet)
  newColumns <- ncol(weightMatrix)
  somSize <- sqrt(nrow(weightMatrix))
  cycleIntern <- 0
  maxCycleIntern <- cycles * nrowDataSet
  currentRadius <- initRadius;
  learnRate <- initLearnRate;
  inverseLearnRateC <- -learnRateReduction;
  progressStep = 1;

  # if (cycles >= 100){
  #   progressStep <-  (cycles / 100.0);
  #     cat("1%                                                                                              100%\n")
  # }else{
  #   if (!updateParametersPerEpoch && maxCycleIntern > 100){
  #     progressStep <- (maxCycleIntern %/% 100.0);
  #     cat("1%                                                                                              100%\n")
  #   }
  # }
  if (maxCycleIntern > 100){
    progressStep <- (maxCycleIntern %/% 100.0);
    cat("1%                                                                                              100%\n")
  }
  for (cycle in 1:cycles){
    # if (cycle %% progressStep == 0 && cycles >= 100 && updateParametersPerEpoch){
    #   cat(".")
    # }
    data_order = calculateDataOrder(nrowDataSet, sampling)
    for (it in data_order){
      cycleIntern <- cycleIntern + 1
      if (!updateParametersPerEpoch){
        if (radiusReduction < 0){
          currentRadius <- 1+(initRadius-1)*(maxCycleIntern-cycleIntern)/(maxCycleIntern)
        }else{
          currentRadius <- currentRadius - radiusReduction;
        }
        if (learnRateReduction == 0){
          # print(paste(currentRadius,maxCycleIntern,cycleIntern,learnRate,sep="--"))

          learnRate <- initLearnRate*(1.0-(cycleIntern)/maxCycleIntern)
        }else{
          if (learnRateReduction > 0){
            learnRate <- learnRate - learnRateReduction;
          }else{
            learnRate <- initLearnRate * inverseLearnRateC / (inverseLearnRateC + cycleIntern)
          }
        }
      }
      # print(paste(currentRadius,maxCycleIntern,cycleIntern,learnRate,sep="--"))
      if (cycleIntern %% progressStep == 0 && maxCycleIntern >= 100){
        cat(".")
      }
      x = dataSet[it,]
      resultDelta = calculateDelta(weightMatrix, x, naExist);
      euclidianDistances2 = calculateEuclidianDistances(resultDelta, oldColumns)
      winner = which.min(euclidianDistances2);
      if (normType == 1){
        weightMatrix <- weightMatrix + learnRate * (resultDelta *
                        tableToCodebookMatrix(somSize, winner, newColumns, calculateNeighborhoodTable(somSize, currentRadius))
                        )
      }else{
        nbMatrix = calculateNeighborhoodMatrix(winner, somSize, currentRadius);
        dim(nbMatrix) = c(somSize, somSize);
        weightMatrix <- weightMatrix + learnRate * (resultDelta * matrixToCodebookMatrix(nbMatrix, newColumns))
      }
    }
    if (updateParametersPerEpoch){
      currentRadius <- 1+(initRadius-1)*(cycles-cycle)/(cycles);
      if (learnRateReduction >= 0){
        learnRate <- initLearnRate*(1.0-(cycle)/cycles);
      }else{
        learnRate <- initLearnRate * inverseLearnRateC / (inverseLearnRateC + cycle);
      }
    }
  }
  cat("\n")
  return(weightMatrix)
}
calculateDataOrder <- function(nrowDataSet,sampling){
  if (sampling == 0){
    return(seq_len(nrowDataSet))
  }else{
    if (sampling == 1){
      return(sample(seq_len(nrowDataSet), nrowDataSet))
    }else{
      return(sample(seq_len(nrowDataSet), sampling, T))
    }
  }
}

#' Codebook Matrix Initialization using PCA
#'
#' This function provides an initial weight matrix for the SOM algorithm by applying
#' the principal component analysis (PCA) to a given data set, based on the
#' \code{som.init} function from the \code{som} package written by Jun Yan.
#' It also supports NA values and extended data sets.
#'
#' @param data The data set.
#' @param somSize The size of the resulting (quadratic) SOM grid of size \code{(somSize, somSize)}.
#' @param oldColumns A boolean vector specifying the columns of the data set which were already
#' available before the data set has been extended. In other words, the SOM initialization only
#' uses these columns for the calculation of Best Matching Units while the other other columns
#' are ignored during this step, but still modified during the update of the weight matrix
#  (codebook matrix).
#'
#' @return The initialized codebook matrix.
#' @export
som.init.extended <- function(data, somSize, oldColumns){
  oldData <- data[,oldColumns]
  rowsWithNoNA <- !apply(oldData,MARGIN = 1, function(x){any(is.na(x))})
  oldDataNoNA <- oldData[rowsWithNoNA,]
  oldW = som::som.init(oldDataNoNA, somSize, somSize, init="linear")
  resultW <- oldW

  if (any(!oldColumns)){
    newData <- data[,!oldColumns]
    rowsWithNoNA <- !apply(data,MARGIN = 1, function(x){any(is.na(x))})#newData,MARGIN = 1, function(x){any(is.na(x))})
    newDataNoNA <- newData[rowsWithNoNA,]
    newW = som::som.init(newDataNoNA, somSize, somSize, init="linear")
    resultW <- cbind(oldW, newW)
  }
  return(resultW)
}

