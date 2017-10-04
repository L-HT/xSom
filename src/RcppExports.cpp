// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// calculateDelta
void calculateDelta(const Rcpp::NumericMatrix& inputMatrix, const Rcpp::NumericVector& inputVector, const bool naExist, Rcpp::NumericMatrix& resultDelta);
RcppExport SEXP xSom_calculateDelta(SEXP inputMatrixSEXP, SEXP inputVectorSEXP, SEXP naExistSEXP, SEXP resultDeltaSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::NumericMatrix& >::type inputMatrix(inputMatrixSEXP);
    Rcpp::traits::input_parameter< const Rcpp::NumericVector& >::type inputVector(inputVectorSEXP);
    Rcpp::traits::input_parameter< const bool >::type naExist(naExistSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix& >::type resultDelta(resultDeltaSEXP);
    calculateDelta(inputMatrix, inputVector, naExist, resultDelta);
    return R_NilValue;
END_RCPP
}
// calculateEuclidianDistances
void calculateEuclidianDistances(const Rcpp::NumericMatrix& deltaMatrix, const Rcpp::LogicalVector& oldColumns, Rcpp::NumericVector& resultEuclidianDistances2);
RcppExport SEXP xSom_calculateEuclidianDistances(SEXP deltaMatrixSEXP, SEXP oldColumnsSEXP, SEXP resultEuclidianDistances2SEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::NumericMatrix& >::type deltaMatrix(deltaMatrixSEXP);
    Rcpp::traits::input_parameter< const Rcpp::LogicalVector& >::type oldColumns(oldColumnsSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector& >::type resultEuclidianDistances2(resultEuclidianDistances2SEXP);
    calculateEuclidianDistances(deltaMatrix, oldColumns, resultEuclidianDistances2);
    return R_NilValue;
END_RCPP
}
// calculateNeighborhoodTable
Rcpp::NumericVector calculateNeighborhoodTable(int somSize, double radius);
RcppExport SEXP xSom_calculateNeighborhoodTable(SEXP somSizeSEXP, SEXP radiusSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type somSize(somSizeSEXP);
    Rcpp::traits::input_parameter< double >::type radius(radiusSEXP);
    rcpp_result_gen = Rcpp::wrap(calculateNeighborhoodTable(somSize, radius));
    return rcpp_result_gen;
END_RCPP
}
// tableToCodebookMatrix
Rcpp::NumericMatrix tableToCodebookMatrix(int somSize, int winnerNeuronR, int xDim, Rcpp::NumericVector lookupTable);
RcppExport SEXP xSom_tableToCodebookMatrix(SEXP somSizeSEXP, SEXP winnerNeuronRSEXP, SEXP xDimSEXP, SEXP lookupTableSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type somSize(somSizeSEXP);
    Rcpp::traits::input_parameter< int >::type winnerNeuronR(winnerNeuronRSEXP);
    Rcpp::traits::input_parameter< int >::type xDim(xDimSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type lookupTable(lookupTableSEXP);
    rcpp_result_gen = Rcpp::wrap(tableToCodebookMatrix(somSize, winnerNeuronR, xDim, lookupTable));
    return rcpp_result_gen;
END_RCPP
}
// calculateNeighborhoodMatrix
void calculateNeighborhoodMatrix(const int& winnerNeuronR, const int& somSize, const double& radius, Rcpp::NumericVector& resultVector);
RcppExport SEXP xSom_calculateNeighborhoodMatrix(SEXP winnerNeuronRSEXP, SEXP somSizeSEXP, SEXP radiusSEXP, SEXP resultVectorSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const int& >::type winnerNeuronR(winnerNeuronRSEXP);
    Rcpp::traits::input_parameter< const int& >::type somSize(somSizeSEXP);
    Rcpp::traits::input_parameter< const double& >::type radius(radiusSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector& >::type resultVector(resultVectorSEXP);
    calculateNeighborhoodMatrix(winnerNeuronR, somSize, radius, resultVector);
    return R_NilValue;
END_RCPP
}
// matrixToCodebookMatrix
void matrixToCodebookMatrix(const Rcpp::NumericVector& matrixAsVector, Rcpp::NumericMatrix& result);
RcppExport SEXP xSom_matrixToCodebookMatrix(SEXP matrixAsVectorSEXP, SEXP resultSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::NumericVector& >::type matrixAsVector(matrixAsVectorSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix& >::type result(resultSEXP);
    matrixToCodebookMatrix(matrixAsVector, result);
    return R_NilValue;
END_RCPP
}
// findWinningNeuron
Rcpp::List findWinningNeuron(Rcpp::NumericMatrix weightMatrix, Rcpp::NumericVector x, Rcpp::LogicalVector oldColumns);
RcppExport SEXP xSom_findWinningNeuron(SEXP weightMatrixSEXP, SEXP xSEXP, SEXP oldColumnsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type weightMatrix(weightMatrixSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< Rcpp::LogicalVector >::type oldColumns(oldColumnsSEXP);
    rcpp_result_gen = Rcpp::wrap(findWinningNeuron(weightMatrix, x, oldColumns));
    return rcpp_result_gen;
END_RCPP
}
// learnCyclesExtended
Rcpp::NumericMatrix learnCyclesExtended(Rcpp::NumericMatrix dataSet, Rcpp::NumericMatrix weightMatrix, Rcpp::LogicalVector oldColumns, unsigned int cycles, double initLearnRate, double learnRateReduction, double initRadius, double radiusReduction, int normType, int sampling, bool naExist, bool updateParametersPerEpoch, unsigned int currentTrainingStep);
RcppExport SEXP xSom_learnCyclesExtended(SEXP dataSetSEXP, SEXP weightMatrixSEXP, SEXP oldColumnsSEXP, SEXP cyclesSEXP, SEXP initLearnRateSEXP, SEXP learnRateReductionSEXP, SEXP initRadiusSEXP, SEXP radiusReductionSEXP, SEXP normTypeSEXP, SEXP samplingSEXP, SEXP naExistSEXP, SEXP updateParametersPerEpochSEXP, SEXP currentTrainingStepSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type dataSet(dataSetSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type weightMatrix(weightMatrixSEXP);
    Rcpp::traits::input_parameter< Rcpp::LogicalVector >::type oldColumns(oldColumnsSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type cycles(cyclesSEXP);
    Rcpp::traits::input_parameter< double >::type initLearnRate(initLearnRateSEXP);
    Rcpp::traits::input_parameter< double >::type learnRateReduction(learnRateReductionSEXP);
    Rcpp::traits::input_parameter< double >::type initRadius(initRadiusSEXP);
    Rcpp::traits::input_parameter< double >::type radiusReduction(radiusReductionSEXP);
    Rcpp::traits::input_parameter< int >::type normType(normTypeSEXP);
    Rcpp::traits::input_parameter< int >::type sampling(samplingSEXP);
    Rcpp::traits::input_parameter< bool >::type naExist(naExistSEXP);
    Rcpp::traits::input_parameter< bool >::type updateParametersPerEpoch(updateParametersPerEpochSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type currentTrainingStep(currentTrainingStepSEXP);
    rcpp_result_gen = Rcpp::wrap(learnCyclesExtended(dataSet, weightMatrix, oldColumns, cycles, initLearnRate, learnRateReduction, initRadius, radiusReduction, normType, sampling, naExist, updateParametersPerEpoch, currentTrainingStep));
    return rcpp_result_gen;
END_RCPP
}
// somCheckNa
Rcpp::NumericMatrix somCheckNa(Rcpp::NumericMatrix dataSet, Rcpp::NumericMatrix weightMatrix, Rcpp::LogicalVector oldColumns, Rcpp::IntegerVector rlen, Rcpp::NumericVector initLearnRate, Rcpp::NumericVector initRadius, double radiusReduction, double learnRateReduction, int normType, int sampling, bool naExist, bool updateParametersPerEpoch);
RcppExport SEXP xSom_somCheckNa(SEXP dataSetSEXP, SEXP weightMatrixSEXP, SEXP oldColumnsSEXP, SEXP rlenSEXP, SEXP initLearnRateSEXP, SEXP initRadiusSEXP, SEXP radiusReductionSEXP, SEXP learnRateReductionSEXP, SEXP normTypeSEXP, SEXP samplingSEXP, SEXP naExistSEXP, SEXP updateParametersPerEpochSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type dataSet(dataSetSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type weightMatrix(weightMatrixSEXP);
    Rcpp::traits::input_parameter< Rcpp::LogicalVector >::type oldColumns(oldColumnsSEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type rlen(rlenSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type initLearnRate(initLearnRateSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type initRadius(initRadiusSEXP);
    Rcpp::traits::input_parameter< double >::type radiusReduction(radiusReductionSEXP);
    Rcpp::traits::input_parameter< double >::type learnRateReduction(learnRateReductionSEXP);
    Rcpp::traits::input_parameter< int >::type normType(normTypeSEXP);
    Rcpp::traits::input_parameter< int >::type sampling(samplingSEXP);
    Rcpp::traits::input_parameter< bool >::type naExist(naExistSEXP);
    Rcpp::traits::input_parameter< bool >::type updateParametersPerEpoch(updateParametersPerEpochSEXP);
    rcpp_result_gen = Rcpp::wrap(somCheckNa(dataSet, weightMatrix, oldColumns, rlen, initLearnRate, initRadius, radiusReduction, learnRateReduction, normType, sampling, naExist, updateParametersPerEpoch));
    return rcpp_result_gen;
END_RCPP
}
// somWithMapping
Rcpp::List somWithMapping(Rcpp::NumericMatrix dataSet, Rcpp::NumericMatrix weightMatrix, Rcpp::LogicalVector oldColumns, Rcpp::IntegerVector rlen, Rcpp::NumericVector initLearnRate, Rcpp::NumericVector initRadius, double radiusReduction, double learnRateReduction, int normType, int sampling, bool naExist, bool updateParametersPerEpoch);
RcppExport SEXP xSom_somWithMapping(SEXP dataSetSEXP, SEXP weightMatrixSEXP, SEXP oldColumnsSEXP, SEXP rlenSEXP, SEXP initLearnRateSEXP, SEXP initRadiusSEXP, SEXP radiusReductionSEXP, SEXP learnRateReductionSEXP, SEXP normTypeSEXP, SEXP samplingSEXP, SEXP naExistSEXP, SEXP updateParametersPerEpochSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type dataSet(dataSetSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type weightMatrix(weightMatrixSEXP);
    Rcpp::traits::input_parameter< Rcpp::LogicalVector >::type oldColumns(oldColumnsSEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type rlen(rlenSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type initLearnRate(initLearnRateSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type initRadius(initRadiusSEXP);
    Rcpp::traits::input_parameter< double >::type radiusReduction(radiusReductionSEXP);
    Rcpp::traits::input_parameter< double >::type learnRateReduction(learnRateReductionSEXP);
    Rcpp::traits::input_parameter< int >::type normType(normTypeSEXP);
    Rcpp::traits::input_parameter< int >::type sampling(samplingSEXP);
    Rcpp::traits::input_parameter< bool >::type naExist(naExistSEXP);
    Rcpp::traits::input_parameter< bool >::type updateParametersPerEpoch(updateParametersPerEpochSEXP);
    rcpp_result_gen = Rcpp::wrap(somWithMapping(dataSet, weightMatrix, oldColumns, rlen, initLearnRate, initRadius, radiusReduction, learnRateReduction, normType, sampling, naExist, updateParametersPerEpoch));
    return rcpp_result_gen;
END_RCPP
}
// rearrangeNeuronCoordinates
Rcpp::List rearrangeNeuronCoordinates(Rcpp::List neuronCoordinates, int somSize);
RcppExport SEXP xSom_rearrangeNeuronCoordinates(SEXP neuronCoordinatesSEXP, SEXP somSizeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::List >::type neuronCoordinates(neuronCoordinatesSEXP);
    Rcpp::traits::input_parameter< int >::type somSize(somSizeSEXP);
    rcpp_result_gen = Rcpp::wrap(rearrangeNeuronCoordinates(neuronCoordinates, somSize));
    return rcpp_result_gen;
END_RCPP
}
