// These functions are based on Jonathan Olmsted's <jpolmsted@gmail.com>;
// RcppTN 0.1-8 source codes (https://github.com/olmjo/RcppTN); Correct JO's
// codes based on Christopher Jackson's <chris.jackson@mrc-bsu.cam.ac.uk>
// R codes

#include <RcppArmadillo.h>

 inline bool CheckSimple(const double lower,  const double upper)
{
  // Check if simpler subalgorithm is appropriate.
  // Init Values Used in Inequality of Interest
  double val1 = (2 * sqrt(exp(1))) / (lower + sqrt(pow(lower, 2) + 4));
  double val2 = exp((pow(lower, 2) - lower * sqrt(pow(lower, 2) + 4)) / (4)) ;
  //

  // Test if Simple is Preferred
  if (upper > lower + val1 * val2) {
    return true ;
  } else {
    return false ;
  }
}

// tn_check and rtn_check are the same
bool tn_check(const double mean, const double sd,
  const double lower, const double upper)
{
  bool flag = true ;
  if (sd <= 0) {
    flag = false ;
  }
  if ((sd == R_NegInf) | (sd == R_PosInf)) {
    flag = false ;
  }
  if ((mean == R_NegInf) | (mean == R_PosInf)) {
    flag = false ;
  }
  if (lower >= upper) {
    flag = false ;
  }
  return(flag) ;
}

// Naive Accept-Reject algorithm; CJ's standard method (0)
inline double alg1(const double lower, const double upper)
{
  int valid = 0 ;   // Init Valid Flag
  double z = 0.0 ;  // Init Draw Storage

  // Loop Until Valid Draw
  while (valid == 0)
  {
    z = Rf_rnorm(0.0, 1.0) ;
    if (z <= upper && z >= lower) { valid = 1 ; }
  }
  return z ;
}

// Accept-Reject Algorithm 2; CJ's alg==1 ie expl
inline double alg2(const double lower, const double upper)
{
  // Init flag and values
  int valid  = 0 ;  // p 122, right column
  const double alphaStar = (lower + sqrt(pow(lower, 2.0) + 4.0)) / 2.0 ;
  double z   = 0 ;
  double u   = 0 ;
  double rho = 0 ;

  // Loop until valid draw
  while (valid == 0)
  {
    z = Rf_rexp(alphaStar) + lower; // control lower boundary
    u = Rf_runif(0, 1) ;
    rho = exp(-pow(z - alphaStar, 2) / 2) ;
    if (u <= rho && z <= upper) { valid = 1 ; } // Keep Successes
  }
  return z ;
}

// Accept-Reject Algorithm 3; CJ's alg==2 ie expu
inline double alg3(const double lower, const double upper)
{
  // Init flag and values
  int valid  = 0 ;  // p 122, right column
  const double alphaStar = (-upper + sqrt(pow(upper, 2.0) + 4.0)) / 2.0 ;
  double z   = 0 ;
  double u   = 0 ;
  double rho = 0 ;

  // Loop until valid draw
  while (valid == 0)
  {
    z = Rf_rexp(alphaStar) - upper; // control upper boundary
    u = Rf_runif(0, 1) ;
    rho = exp(-pow(z - alphaStar, 2) / 2) ;
    if (u <= rho && z <= -lower) { valid = 1 ; } // Keep Successes
  }
  return z ;
}

// Accept-Reject Algorithm 4; page 123. 2.2. Two-sided truncated normal dist.
inline double alg4(const double lower, const double upper)
{
  int valid = 0 ;
  double z = 0 ;
  double u = 0 ;
  double rho = 0 ;

  while (valid == 0)  // Loop Until Valid Draw
  {
    z = Rf_runif(lower, upper) ;
    if (lower > 0)
    {
      rho = exp( (pow(lower, 2) - pow(z, 2)) / 2 ) ;
    }
    else if (upper < 0)
    {
      rho = exp( (pow(upper, 2) - pow(z, 2)) / 2 ) ;
    } else  // 0 belongs to [lower, upper]
    {
      rho = exp( -pow(z, 2) / 2) ;
    }
    u = Rf_runif(0, 1) ;
    if (u <= rho) { valid = 1 ; }
  }
  return z ;
}

//' @export
double rtn_scalar(const double mean,  const double sd, const double lower,
  const double upper)
{
  // Standardised upper and lower
  double stdLower = (lower - mean) / sd ; // Algorithm works on mean 0, sd 1 scale
  double stdUpper = (upper - mean) / sd ; // bounds are standardised
  double draw ;

/* Different algorithms depending on where upper/lower limits lie.
 * CJ's -1. return NaN if lower > upper
 * CJ's  0. standard "simulate from normal and reject if outside limits" method.
            Use if bounds are wide.
 * CJ's  1. rejection sampling with exponential proposal. Use if lower >> mean
 * CJ's  2. rejection sampling with exponential proposal. Use if upper << mean.
 * CJ's  3. rejection sampling with uniform proposal. Use if bounds are narrow and
 */

  // CJ's -1
  bool a0 = stdLower > stdUpper ;
  // CJ's 0
  bool a1 = (stdLower < 0 && stdUpper==INFINITY)  ||
    (stdLower == -INFINITY  && stdUpper > 0) ||
    (std::isfinite(stdLower) && std::isfinite(stdUpper) &&
    (stdLower < 0) && (stdUpper > 0) &&
    (stdUpper - stdLower > sqrt(2*PI))) ;

  // CJ's 1
  double term1_a2 = stdLower ;
  double term2_a2 = 2*sqrt(exp(1)) / (stdLower+sqrt(pow(stdLower, 2.0) + 4.0)) ;
  double term3_a2 = exp( (stdLower*2 - stdLower*sqrt(pow(stdLower, 2.0) + 4.0)) / 4) ;
  double eq_a2 = term1_a2 + term2_a2 * term3_a2;
  bool a2 = (stdLower >= 0) && (stdUpper > eq_a2 ) ;

  // CJ's 2
  double term1_a3 = -stdUpper ;
  double term2_a3 = 2*sqrt(exp(1)) / (-stdUpper+sqrt(pow(stdUpper, 2.0) + 4.0)) ;
  double term3_a3 = exp( (stdUpper*2 - stdUpper*sqrt(pow(stdUpper, 2.0) + 4.0)) / 4) ;
  double eq_a3 = term1_a3 + term2_a3 * term3_a3;
  bool a3 = (stdUpper <= 0) && ( -stdLower > eq_a3) ;


  // CJ's 3 otherwise
    if (a0) { draw = NAN ;}  // nan
    else if (a1)  // no
    {
      draw = alg1(stdLower, stdUpper) ;
    }
    else if (a2)  // expl
    {
      draw = alg2(stdLower, stdUpper) ;
    }
    else if (a3)  // expu
    {
      draw = alg3(stdLower, stdUpper) ;
    }
    else   // u
    {
      draw = alg4(stdLower, stdUpper) ;
    }

    double out = draw * sd + mean ;
    return out ;
}

//' @export
double dtn_scalar(const double x, const double mean, const double sd,
  const double lower, const double upper, int islog)
{
  double dens, numer, denom;
  // Default value for dens either -infinity or 0
  // similar to C. Jackson's
  // ret[x < lower | x > upper] <- if (log) -Inf else 0

  if ((x >= lower) && (x <= upper))
  {
    // Use Chris Jackson's method. Jonathan Olmsted's division method will
    // return infinity when denominator is 0.
    // lower.tail (lt) = 1; log.p (lg) = 0
    denom = R::pnorm(upper, mean, sd, 1, 0) - R::pnorm(lower, mean, sd, 1, 0);
    numer = R::dnorm(x, mean, sd, islog);
    dens = islog ? (numer - log(denom)) : (numer/denom);
  }
  else
  {
    dens = islog ? -INFINITY : 0;
  }
  return(dens) ;
}

//' @export
void rtn(Rcpp::NumericVector &mean, Rcpp::NumericVector &sd,
  Rcpp::NumericVector &lower, Rcpp::NumericVector &upper,
  Rcpp::NumericVector &draws) {

  Rcpp::NumericVector::iterator itM = mean.begin() ;
  Rcpp::NumericVector::iterator itS = sd.begin() ;
  Rcpp::NumericVector::iterator itL = lower.begin() ;
  Rcpp::NumericVector::iterator itH = upper.begin() ;
  Rcpp::NumericVector::iterator itD = draws.begin() ;

  while (itD != draws.end())    // Draw from TN
  {
    if (tn_check(*itM, *itS, *itL, *itH)) {
      *itD = rtn_scalar(*itM, *itS, *itL, *itH) ;
    } else {
      *itD = NA_REAL ;
    }
    itD++ ;
    itM++ ;
    itS++ ;
    itL++ ;
    itH++ ;
  }
}

//' @export
void dtn(Rcpp::NumericVector &x,
  Rcpp::NumericVector &mean,  Rcpp::NumericVector &sd,
  Rcpp::NumericVector &lower, Rcpp::NumericVector &upper,
  Rcpp::LogicalVector &islog, Rcpp::NumericVector &dens) {

  // Initialise containers
  Rcpp::NumericVector::iterator itX     = x.begin() ;
  Rcpp::NumericVector::iterator itM     = mean.begin() ;
  Rcpp::NumericVector::iterator itS     = sd.begin() ;
  Rcpp::NumericVector::iterator itL     = lower.begin() ;
  Rcpp::NumericVector::iterator itU     = upper.begin() ;
  Rcpp::LogicalVector::iterator itIsLog = islog.begin() ;
  Rcpp::NumericVector::iterator itD     = dens.begin() ;

  // Based on Jonathan Olmsted's code. Change values internally at the memory
  // locations without returning
  while (itD != dens.end())
  {
    if (tn_check(*itM, *itS, *itL, *itU)) {
      *itD = dtn_scalar(*itX, *itM, *itS, *itL, *itU, *itIsLog) ;
    } else {
      *itD = NA_REAL ;
    }
    itX++ ;
    itM++ ;
    itS++ ;
    itL++ ;
    itU++ ;
    itD++ ;
    itIsLog++ ;
  }
}

//' @export
RcppExport SEXP rtn_wrapper(const SEXP mean_, const SEXP sd_,
  const SEXP lower_, const SEXP upper_)
{
  Rcpp::NumericVector Mean(mean_) ;
  Rcpp::NumericVector Sd(sd_) ;
  Rcpp::NumericVector Lower(lower_) ;
  Rcpp::NumericVector Upper(upper_) ;

  Rcpp::NumericVector Draws(Mean.size(), NA_REAL) ;
  Rcpp::RNGScope scope ;  // To replicate exactly the same result control this
  rtn(Mean, Sd, Lower, Upper, Draws) ;
  return Draws ;
}

//' @export
RcppExport SEXP dtn_wrapper(const SEXP x_, const SEXP mean_, const SEXP sd_,
  const SEXP lower_, const SEXP upper_, const SEXP islog_)
{
  // Convert input to Rcpp's NumericVector and LogicalVector
  // pass to Rcpp's dtn to get probability density
  Rcpp::NumericVector X(x_) ;
  Rcpp::NumericVector Mean(mean_) ;
  Rcpp::NumericVector Sd(sd_) ;
  Rcpp::NumericVector Lower(lower_) ;
  Rcpp::NumericVector Upper(upper_) ;
  Rcpp::LogicalVector Islog(islog_) ;
  Rcpp::NumericVector Dens(X.size(), 0.0) ;
  dtn(X, Mean, Sd, Lower, Upper, Islog, Dens) ;
  return Dens;
}


