% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/somFunctions.R
\name{calculateMapping}
\alias{calculateMapping}
\title{Generate a SOM mapping from a given codebook matrix}
\usage{
calculateMapping(dataSet, weightMatrix, oldColumns)
}
\arguments{
\item{dataSet}{The data set.}

\item{weightMatrix}{A given weight matrix (codebook matrix).}

\item{oldColumns}{A boolean vector specifying the columns of the data set which were already
available before the data set has been extended. In other words, the SOM algorithm only
uses these columns for the calculation of Best Matching Units while the other other columns
are ignored during this step, but still modified during the update of the weight matrix
(codebook matrix).}
}
\value{
A list described above.
}
\description{
This function is meant to be used in conjunction with \code{learnCyclesExtended} or
\code{learnCyclesExtendedR}. It takes a given codebook matrix and calculates a list with a structure
similar to \code{somWithMappingR}:
\describe{
  \item{\code{codebook}}{The resulting weight matrix (codebook matrix).}
  \item{\code{feature.mapping}}{A data frame specifying the SOM node mapped to each row of the data set
  as well as the position of the nodes on the quadratic SOM grid.}
  \item{\code{node.summary}}{A data frame which shows how many rows of the data set
  are mapped to each node on the quadratic SOM grid.}
}
}
\examples{
#generate Data
dataSet <- matrix(as.numeric(1:400),ncol=2)
weightMatrix <- som.init.extended(dataSet, somSize=3, oldColumns=c(TRUE,TRUE))

#apply the algorithm
result <- learnCyclesExtendedR(dataSet, weightMatrix,
     oldColumns=c(TRUE,TRUE), currentTrainingStep=1)
result <- learnCyclesExtendedR(dataSet, result,
     oldColumns=c(TRUE,TRUE), currentTrainingStep=2)
#result now contains the final result

#calculate the mapping
resultList <- calculateMapping(dataSet, result, oldColumns=c(TRUE,TRUE))
}
