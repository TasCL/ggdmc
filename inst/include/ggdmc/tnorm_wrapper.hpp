inline bool CheckSimple(const double lower,  const double upper) ;

bool tn_check(const double mean, const double sd,
  const double lower, const double upper) ;

// Naive Accept-Reject algorithm; CJ's standard method (0)
inline double alg1(const double lower, const double upper) ;

// Accept-Reject Algorithm 2; CJ's alg==1 ie expl
inline double alg2(const double lower, const double upper) ;

// Accept-Reject Algorithm 3; CJ's alg==2 ie expu
inline double alg3(const double lower, const double upper) ;

// Accept-Reject Algorithm 4; page 123. 2.2. Two-sided truncated normal dist.
inline double alg4(const double lower, const double upper) ;

double rtn_scalar(const double mean,  const double sd, const double lower,
  const double upper) ;

double dtn_scalar(const double x, const double mean, const double sd,
  const double lower, const double upper, int islog) ;

void rtn(Rcpp::NumericVector &mean, Rcpp::NumericVector &sd,
  Rcpp::NumericVector &lower, Rcpp::NumericVector &upper,
  Rcpp::NumericVector &draws) ;

void dtn(Rcpp::NumericVector &x,
  Rcpp::NumericVector &mean,  Rcpp::NumericVector &sd,
  Rcpp::NumericVector &lower, Rcpp::NumericVector &upper,
  Rcpp::LogicalVector &islog, Rcpp::NumericVector &dens) ;

RcppExport SEXP rtn_wrapper(const SEXP mean_, const SEXP sd_,
  const SEXP lower_, const SEXP upper_) ;
RcppExport SEXP dtn_wrapper(const SEXP x_, const SEXP mean_, const SEXP sd_,
  const SEXP lower_, const SEXP upper_, const SEXP islog_) ;


