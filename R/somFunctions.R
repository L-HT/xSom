somTwoStepR <- function(dataSet, weightMatrix,
                       oldColumns,
                       rlen = c(2,10),
                       initLearnRate = c(0.05,0.02),
                       initRadius = c(),
                       radiusReduction = -1.0,
                       learnRateReduction = 0.0,
                       normType = 2, sampling = 0, naExist = T,
                       updateParametersPerEpoch = T){

  result <- weightMatrix
  somSize = sqrt(nrow(weightMatrix));
  str <- paste("Erste Trainingsstufe (rlen = ", rlen[1], " Zyklen, initLearnRate = ", initLearnRate[1], ", initRadius = ",
    initRadius[1], ")", sep="")
  print(str)
  result <- learnCyclesExtendedR(dataSet, weightMatrix, oldColumns, rlen[1],
                              initLearnRate[1], learnRateReduction, initRadius[1],
                              radiusReduction, normType, sampling, naExist,
                              updateParametersPerEpoch
                              )

  if (rlen[2] != 0){
    str <- paste("Zweite Trainingsstufe (rlen = ", rlen[2], " Zyklen, initLearnRate = ", initLearnRate[2], ", initRadius = ",
                 initRadius[2], ")", sep="")
    print(str)
    result <- learnCyclesExtendedR(dataSet, result, oldColumns, rlen[2],
                                 initLearnRate[2], learnRateReduction, initRadius[2],
                                 radiusReduction, normType, sampling, naExist,
                                 updateParametersPerEpoch
                                 )
  }
  return(result)
}

#' A SOM algorithm that handles NA values.
#'
#' In terms of functionality this function is identical to \code{somCheckNaR},
#' but uses more code based on R which might change the runtime of the function.
#' It is recommended to try both \code{somCheckNa} and
#' \code{somCheckNaR} and see which function is faster.
#'
#' If the argument \code{naExist} is set to \code{FALSE}, the first 100 lines of the data set will be
#' scanned for NA values. Depending on whether NA values are found slightly different
#' SOM algorithms are executed:
#' If NA values exist, the algorithm will ignore these values during the calculation of the
#' Best Matching Unit (BMU). If no NA values exist, this step is omitted which improves
#' the runtime of the algorithm.
#' @inheritParams somCheckNa
#' @export
#' @examples
#' # generate data and initialize weight matrix
#' dataSet <- matrix(as.numeric(1:20),ncol=2)
#' weightMatrix <- som.init.extended(dataSet, somSize=2, oldColumns=c(TRUE,TRUE))
#'
#' # apply the algorithm
#' result <- somCheckNaR(dataSet, weightMatrix, oldColumns=c(TRUE,TRUE))
somCheckNaR <- function(dataSet,  weightMatrix,
                        oldColumns,
                        rlen = c(),
                       initLearnRate = c(),
                       initRadius = c(),
                       radiusReduction = -1.0, learnRateReduction = 0.0,
                       normType = 2, sampling = 0, naExist = T,
                       updateParametersPerEpoch = T){

  somSize <- sqrt(nrow(weightMatrix));
  if (length(rlen) == 0){
    rlen <- c(2,10)
  }
  if (length(rlen) == 1){
    rlen <- c(rlen[1], rlen[1]*10)
  }
  if (length(initLearnRate) == 0){
    initLearnRate = c(0.05, 0.02)
  }else{
    if (length(initLearnRate) == 1){
      initLearnRate = c(initLearnRate[1], initLearnRate[1]/2.0);
    }
  }
  if (length(initRadius) == 0.0){
    initRadius = c(somSize, min(somSize,3) );
  }else{
    if (length(initRadius) == 1){
      initRadius = c(initRadius[1],  max(initRadius[1],2) );
    }
  }
  if (naExist){
    result <- somTwoStepR(dataSet, weightMatrix,  oldColumns,rlen, initLearnRate, initRadius,radiusReduction,
               learnRateReduction, normType, sampling, naExist = T, updateParametersPerEpoch)
  }else{
    n <- min(100, nrow(dataSet))
    subset <- dataSet[1:n,]
    if (any(is.na(subset))){
      print(" -- NAs erkannt. Es wird Training mit NAs genutzt.")
      result <- somTwoStepR(dataSet, weightMatrix,  oldColumns,rlen, initLearnRate, initRadius,radiusReduction,
                 learnRateReduction, normType, sampling, naExist = T, updateParametersPerEpoch)
    }else{
      result <- somTwoStepR(dataSet, weightMatrix,  oldColumns,rlen, initLearnRate, initRadius,radiusReduction,
                 learnRateReduction, normType, sampling, naExist = F, updateParametersPerEpoch)
    }
  }
  return(result)

}


#' A SOM algorithm that handles NA values.
#'
#' This function is identical to \code{somCheckNaR}, but returns a list containing the following
#' elements:
#' \describe{
#'   \item{\code{codebook}}{The resulting weight matrix (codebook matrix).}
#'   \item{\code{feature.mapping}}{A data frame specifying the SOM node mapped to each row of the data set
#'   as well as the position of the nodes on the quadratic SOM grid.}
#'   \item{\code{node.summary}}{A data frame which shows how many rows of the data set
#'   are mapped to each node on the quadratic SOM grid.}
#' }
#'
#' In terms of functionality this function is similar to \code{somWithMapping},
#' but \code{somWithMapping} lacks some smaller features (namely the node.summary), so
#' as of now, it is recommended to use this function instead of \code{somWithMapping}.
#'
#' @inheritParams somCheckNa
#' @return A list described above.
#' @export
#' @examples
#' # generate data and initialize weight matrix
#' dataSet <- matrix(as.numeric(1:20),ncol=2)
#' weightMatrix <- som.init.extended(dataSet, somSize=2, oldColumns=c(TRUE,TRUE))
#'
#' # apply the algorithm
#' resultList <- somWithMappingR(dataSet, weightMatrix, oldColumns=c(TRUE,TRUE))
somWithMappingR <- function(dataSet, weightMatrix, oldColumns,
                               rlen = c(),
                               initLearnRate = c(),
                               initRadius = c(),
                               radiusReduction = -1.0,
                               learnRateReduction = 0.0,
                               normType = 2, sampling = 0, naExist = T,
                               updateParametersPerEpoch = T
                               ){

  result = somCheckNaR(dataSet, weightMatrix,oldColumns, rlen, initLearnRate, initRadius,
                                          radiusReduction, learnRateReduction, normType, sampling,
                                          naExist, updateParametersPerEpoch);
  N <- nrow(dataSet)
  somSize <- sqrt(nrow(weightMatrix))
  dataCounter <- seq_len(N);
  neuronNumberVector <- rep_len(0, N);
  neuronXVector <- rep_len(0, N);
  neuronYVector <- rep_len(0, N);

  node.summary <- cbind.data.frame(x = rep(1:somSize, somSize), y = rep(1:somSize, each=somSize),
                                   n.features = rep(0, somSize*somSize))

  for (i in 1:nrow(dataSet)){
    neuronCoordinates = findWinningNeuron(result, dataSet[i,], oldColumns);
    neuronNumberVector[i] = neuronCoordinates[["neuronNumber"]];
    neuronXVector[i] = neuronCoordinates[["neuronX"]];
    neuronYVector[i] = neuronCoordinates[["neuronY"]];
    # if (rearrangeRows){
    #   neuronCoordinates <- rearrangeNeuronCoordinates(neuronCoordinates, somSize)
    # }

    # Zähltabelle
    oldValue <- node.summary[node.summary$x == neuronCoordinates["neuronX"] & node.summary$y ==
                   neuronCoordinates["neuronY"], "n.features"]
    node.summary[node.summary$x == neuronCoordinates["neuronX"] & node.summary$y ==
                   neuronCoordinates["neuronY"], "n.features"] <- oldValue + 1
  }

  # if (rearrangeRows){
  #   result <- result[calculateNewRowOrder(somSize),]
  # }
  return(
    list(
      "codebook" = result,
      "feature.mapping" = data.frame(
        "dataRowNumber" = dataCounter,
        "neuronNumber" = neuronNumberVector,
        "neuronX" = neuronXVector,
        "neuronY" = neuronYVector
      ),
      "node.summary" = node.summary
    )
  );
}

# About \code{rearrangeRows}:
#
# In the default setting, the neurons in the SOM grid are numbered as follows:
# \preformatted{1 2 3
#
# 4 5 6
#
# 7 8 9}
#
# By enabling \code{rearrangeRows}, the weight matrix's rows and the resulting
# mapping will be regarded as follows:
# \preformatted{7 8 9
#
# 4 5 6
#
# 1 2 3}
#
# @param useRearrangeRows A boolean vector specifying whether a different numbering of the neurons
# in the SOM grid is used for the argument \code{weightMatrix}. More details can be found below. The default setting is \code{FALSE}.

# calculateNewRowOrder <- function(somSize){
#   k <- seq_len(somSize * somSize)
#   result <- somSize^2 - ((k-1) %/% somSize) * somSize + (k-1) %% somSize - (somSize - 1)
#   return(result)
# }

