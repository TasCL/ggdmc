#include <RcppArmadillo.h>

bool checkDDM(std::vector<double> pVec) ;

Rcpp::NumericMatrix getAccumulatorMatrix(Rcpp::NumericVector pVec,
  std::string cell, Rcpp::NumericVector model, bool n1order) ;

std::vector<double> ddmc(Rcpp::DataFrame x, Rcpp::NumericVector pVec,
  double precision, double minLike) ;

std::vector<double> ddmc_parallel(Rcpp::DataFrame x, Rcpp::NumericVector pVec,
  double precision, double minLike) ;
