#include <Rcpp.h>
// [[Rcpp::depends(RcppParallel)]]
#include <RcppParallel.h>
#include <iostream>

#include "SimpleMatrixOperations.h"

struct MatrixTimesScalarMultiplier : RcppParallel::Worker{
  double scalar_;
  RcppParallel::RMatrix<double> input_;
  RcppParallel::RMatrix<double> output_;

  MatrixTimesScalarMultiplier(double scalar, Rcpp::NumericMatrix input, Rcpp::NumericMatrix output)
    : scalar_(scalar), input_(input), output_(output){

  }

  void operator()(std::size_t begin, std::size_t end){
    std::transform(
        input_.begin() + begin,
        input_.begin() + end,
        output_.begin() + begin,
        std::bind1st(std::multiplies<double>(),scalar_)
    );
  }
};

Rcpp::NumericMatrix matrixTimesScalar(Rcpp::NumericMatrix m, double a){
  Rcpp::NumericMatrix result(m.nrow(), m.ncol());
  MatrixTimesScalarMultiplier mtsm(a, m, result);
  RcppParallel::parallelFor(0, m.length(), mtsm);
  return result;
}

//////////////////////

struct ElementwiseMultiplier : RcppParallel::Worker{
  RcppParallel::RMatrix<double> m1_;
  RcppParallel::RMatrix<double> m2_;
  RcppParallel::RMatrix<double> output_;

  ElementwiseMultiplier(Rcpp::NumericMatrix m1, Rcpp::NumericMatrix m2, Rcpp::NumericMatrix output)
    :m1_(m1), m2_(m2), output_(output){

  }

  void operator()(std::size_t begin, std::size_t end){
    std::transform(
      m1_.begin() + begin,
      m1_.begin() + end,
      m2_.begin() + begin,
      output_.begin() + begin,
      std::multiplies<double>()
    );
  }
};

Rcpp::NumericMatrix elementwiseMultiplication(Rcpp::NumericMatrix m1, Rcpp::NumericMatrix m2){
  Rcpp::NumericMatrix result(m1.nrow(), m2.ncol());
  if (m1.nrow() != m2.nrow() || m1.ncol() != m2.ncol()){
    std::cerr << "Matrixdimensionen passen nicht!" << std::endl;
  }
  ElementwiseMultiplier em(m1, m2, result);
  RcppParallel::parallelFor(0, m1.length(), em);
  return result;
}

///////////////////////////////////////////////////////////
struct ElementwiseAdder : RcppParallel::Worker{
  RcppParallel::RMatrix<double> m1_;
  RcppParallel::RMatrix<double> m2_;
  RcppParallel::RMatrix<double> output_;

  ElementwiseAdder(Rcpp::NumericMatrix m1, Rcpp::NumericMatrix m2, Rcpp::NumericMatrix output)
    :m1_(m1), m2_(m2), output_(output){

  }

  void operator()(std::size_t begin, std::size_t end){
    //std::cout << begin << " , " << end;
    std::transform(
      m1_.begin() + begin,
      m1_.begin() + end,
      m2_.begin() + begin,
      output_.begin() + begin,
      std::plus<double>()
    );
  }
};

Rcpp::NumericMatrix elementwiseAddition(Rcpp::NumericMatrix m1, Rcpp::NumericMatrix m2){
  Rcpp::NumericMatrix result(m1.nrow(), m2.ncol());
  if (m1.nrow() != m2.nrow() || m1.ncol() != m2.ncol()){
    std::cerr << "Matrixdimensionen passen nicht!" << std::endl;
  }
  ElementwiseAdder ea(m1, m2, result);
  RcppParallel::parallelFor(0, m1.length(), ea);
  return result;
}

