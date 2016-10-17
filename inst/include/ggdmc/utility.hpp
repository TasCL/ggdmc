Rcpp::NumericMatrix na_matrix(int nr, int nc) ;

arma::uvec getArmaIdx(Rcpp::NumericVector x, int trueOrFalse) ;

Rcpp::List assign_pp_pLists(Rcpp::List& samples, Rcpp::List& usePhi) ;

void assign_pp_pList(Rcpp::List& phi, Rcpp::List& pList) ;

Rcpp::List swap(Rcpp::NumericVector& pVec, Rcpp::List& pList) ;

Rcpp::List isConstant(Rcpp::List& ppList) ;

Rcpp::List resize( const Rcpp::List& x, int n ) ;

arma::cube get_ps(Rcpp::List& samples) ;

arma::cube get_ps_use(Rcpp::List samples) ;

void add_nmc0 (Rcpp::List& samples, Rcpp::List& pLists) ;

double summed_log_likelihood(arma::vec& pVec, Rcpp::List& data)  ;

double summed_log_likelihood_parallel(arma::vec& pVec, Rcpp::List& data)  ;

double summed_log_likelihood_hyper(arma::mat& data_hyper,
  Rcpp::List& pList_hyper) ;

double summed_log_prior(arma::vec& pVec, Rcpp::List& pList) ;

double summed_log_prior_hyper(Rcpp::List& pVec, Rcpp::List& ppList) ;

arma::vec update_useLogPrior(arma::mat& theta, Rcpp::List& pLists) ;
