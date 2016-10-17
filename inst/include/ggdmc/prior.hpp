#include <RcppArmadillo.h>

Rcpp::NumericVector dprior(Rcpp::NumericVector& parameterVector,
  Rcpp::List& priorList) ;

Rcpp::NumericVector rprior_scalar(Rcpp::List priorList) ;

Rcpp::NumericMatrix rprior(int n, Rcpp::List priorList) ;

