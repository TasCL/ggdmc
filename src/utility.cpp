//   Copyright (C) <2016>  <Yi-Shin Lin>
//
//   This program is free software; you can redistribute it and/or modify
//   it under the terms of the GNU General Public License as published by
//   the Free Software Foundation; version 2
//
//   This program is distributed in the hope that it will be useful,
//   but WITHOUT ANY WARRANTY; without even the implied warranty of
//   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//   GNU General Public License for more details.
//
//   You should have received a copy of the GNU General Public License along
//   with this program; if not, write to the Free Software Foundation, Inc.,
//   51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
#include <RcppArmadillo.h>
#include <ggdmc.hpp>

Rcpp::NumericMatrix na_matrix(int nr, int nc) {
  Rcpp::NumericMatrix m(nr, nc) ;
  std::fill( m.begin(), m.end(), Rcpp::NumericVector::get_na() ) ;
  return m ;
}

arma::uvec getArmaIdx(Rcpp::NumericVector x, int trueOrFalse) {
  arma::vec xArma(x.begin(), x.size(), false) ;
  arma::uvec index = find(xArma == trueOrFalse) ;
  return index ;
}

/* This function uses 1st subject's p.prior as a template and replaces
 the location and scale for each parameter with those generated from
hyper distributions. The return is a nChains list with each has its own
p.prior */
// [[Rcpp::export]]
Rcpp::List assign_pp_pLists(Rcpp::List& samples_s1, Rcpp::List& usePhi) {
  Rcpp::List pList         = samples_s1["p.prior"] ;
  int nChains              = samples_s1["n.chains"] ;
  arma::mat usePhiLocation = usePhi[0] ; // nChains x npars
  arma::mat usePhiScale    = usePhi[1] ;
  Rcpp::List pLists(nChains) ; //p.priors (plural)
  Rcpp::List dist_setting ;

  // let usePhi follow normal order eg 0, 1, ..., 23
  // But shuffle the chain sequence
  Rcpp::IntegerVector chainSeq = Rcpp::seq_len(nChains) ; /* eg 1:24 */

  // obtain a time-based seed:
  unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
  std::shuffle(chainSeq.begin(), chainSeq.end(), std::default_random_engine(seed));

  for(int i =0; i<nChains; i++)   // eg 0 to 23
  {
    int ii = chainSeq[i] - 1 ;
    Rcpp::List pList_cp(clone(pList)) ; // Make hard copies several times
    for(Rcpp::List::iterator it2=pList_cp.begin(); it2!=pList_cp.end(); ++it2)
    {
      int j = std::distance(pList_cp.begin(), it2) ; // j is ordered
      dist_setting    = *it2 ; // pList_cp's memory location is altered
      dist_setting[0] = as_scalar(usePhiLocation.row(i).col(j)) ;
      dist_setting[1] = as_scalar(usePhiScale.row(i).col(j)) ;
    }
    pLists[ii] = pList_cp ; // ii is random pList_cp follows the order
  }
return pLists ;
}

/* been used in summed_log_likelihood_hyper. This function changes
 * memory location directly; input pList will be modified */
void assign_pp_pList(Rcpp::List& phi, Rcpp::List& pList) {
  // To replace 1st assign.pp; locat and scale has size of npar
  arma::vec location = phi[0] ;
  arma::vec scale    = phi[1] ;
  for(Rcpp::List::iterator it = pList.begin(); it != pList.end(); ++it)
  { // Loop through a, v, z, sv, ...
    int i = std::distance(pList.begin(), it) ;
    Rcpp::List pList_same_mem_loc = *it ;
    pList_same_mem_loc[0] = location(i) ;
    pList_same_mem_loc[1] = scale(i) ;
  }
}

Rcpp::List swap(Rcpp::NumericVector& pVec, Rcpp::List& pList) {
  std::vector<std::string> pVecNames  = pVec.names() ;
  std::vector<std::string> pListNames = pList.names() ;
  int pVecSize  = pVec.size() ;
  int pListSize = pList.size() ;
  Rcpp::List out(pList.size()) ;

  if (pVecSize != pListSize)
  {
    Rcpp::stop("p.vector does not match parameter prior list.") ;
  }

  for (std::vector<std::string>::iterator it1 = pVecNames.begin();
    it1 != pVecNames.end(); ++it1)
  {
    int i = std::distance(pVecNames.begin(), it1) ;
    if (pVecNames[i] != pListNames[i])
    {
      Rcpp::Rcout << "Found an element in p.prior and p.vector is ordered "
                  << "differently. p.prior will be rearranged to follow "
                  << "p.vector.\n" ;
      for (Rcpp::List::iterator it2 = pList.begin(); it2 != pList.end(); ++it2)
      {
        int j = std::distance(pList.begin(), it2) ;
        if (pListNames[j] == pVecNames[i])
        {
          out[i] = pList[j];
        }
      }
    } else
    {
      out[i] = pList[i];
    }
  }

  out.names() = pVecNames ;
  return out ;
}

Rcpp::List isConstant(Rcpp::List& ppList) {
  Rcpp::List location = ppList[0] ;
  Rcpp::List scale    = ppList[1] ;
  Rcpp::LogicalVector location_out(location.size()) ;
  Rcpp::LogicalVector scale_out(scale.size()) ;
  for(int i=0; i < location.size(); i++)
  {
    Rcpp::List a = location[i] ;
    std::string b = a.attr("dist") ;
    location_out[i] = b == "constant" ;
  }

  for(int j=0; j < scale.size(); j++)
  {
    Rcpp::List a = scale[j] ;
    std::string b = a.attr("dist") ;
    scale_out[j] = b == "constant" ;
  }

  Rcpp::List out = Rcpp::List::create(
    Rcpp::Named("location") = location_out,
    Rcpp::Named("scale")    = scale_out) ;

  return out;
}

Rcpp::List resize( const Rcpp::List& x, int n ) {
  int oldsize = x.size() ;
  Rcpp::List y(n) ;
  for( int i=0; i<oldsize; i++) y[i] = x[i] ;
  return y ;
}

arma::cube get_ps(Rcpp::List& samples) {
  /* Collect 1st nmc and reshape ps cube as nSubjects x npars x nChains
   * Andrew's ps array is nChains x nSubjects x npars; but Armadillo's cube
   * uses slice order first, so chain dimension has to be in slice */
  Rcpp::List hyper = samples.attr("hyper") ;
  int nChains      = hyper["n.chains"] ;
  int nsubs        = samples.size() ;
  int npars        = hyper["n.pars"] ;
  int start_R      = hyper["start"] ; // Note start is R index

  arma::cube ps(nsubs, npars, nChains) ;
  for(Rcpp::List::iterator it = samples.begin(); it != samples.end(); ++it)
  {
    int it_i = std::distance(samples.begin(), it) ;
    Rcpp::List subject_i = *it ;
    arma::cube theta_s   = subject_i["theta"] ;        // nChain x npars x nmc
    arma::mat theta_nmc1 = theta_s.slice(start_R - 1); // take 1st nmc; nChain x npars
    for(int i=0; i<nChains; i++)
    {
      ps.slice(i).row(it_i) = theta_nmc1.row(i) ; // nsubs x npars x nchains
    }
  }
  return ps ;
}

arma::cube get_ps_use(Rcpp::List samples) {
  Rcpp::List samples1 = samples[0] ;
  Rcpp::List use1     = samples1.attr("use") ;
  arma::mat theta1    = use1["theta"] ;
  int nChains         = theta1.n_rows ;
  int nsubs           = samples.size() ;
  int npars           = theta1.n_cols ;

  arma::cube ps(nsubs, npars, nChains) ;
  for(Rcpp::List::iterator it = samples.begin(); it != samples.end(); ++it)
  {
    int it_i = std::distance(samples.begin(), it) ;
    Rcpp::List subject_i = *it ;
    Rcpp::List use       = subject_i.attr("use") ;
    arma::mat theta_i    = use["theta"] ;        // nChain x npars
    for(int i=0; i<nChains; i++)
    { // nsubs x npars x nchains
      ps.slice(i).row(it_i) = theta_i.row(i) ;
    }
  }
  return ps ;
}

void add_nmc0 (Rcpp::List& samples, Rcpp::List& pLists) {
  for(Rcpp::List::iterator it=samples.begin(); it!=samples.end(); ++it)
  {
    int it_i = std::distance(samples.begin(), it) ;
    Rcpp::List samples_i  = *it ;
    arma::cube theta_i    = samples_i["theta"] ; // nChains x npars x nmc
    arma::mat logprior_i  = samples_i["summed_log_prior"] ; // nmc x nChains
    arma::mat loglike_i   = samples_i["log_likelihoods"] ;  // nmc x nChains
    int start_R  = samples_i["start"] ; // R index
    int start_C  = start_R - 1 ;        // C index
    arma::mat theta    = theta_i.slice(start_C) ; // nChains x npars
    arma::vec logprior = vectorise(logprior_i.row(start_C)) ;   // nChains x 1
    arma::vec loglike  = vectorise(loglike_i.row(start_C)) ;   // nChains x 1
    Rcpp::List use = Rcpp::List::create(
      Rcpp::Named("theta")    = theta,    // nChains x npars
      Rcpp::Named("logprior") = logprior, // nChains x 1
      Rcpp::Named("loglike")  = loglike,  // nChains x 1
      Rcpp::Named("store_i")  = start_C,
      Rcpp::Named("p.priors") = pLists) ;
    samples_i.attr("use") = use ;
    samples[it_i] = samples_i ;
  }
}

//' Sum and Log Probability Density of a EAM model
//'
//' Get log likelihood summed over a model data instance.  The input data has
//' to be a data frame carrying a model specification, which is usually
//' created by \code{model.data.dmc} function.
//' \code{summed_log_likelihood_parallel} does calculation in parallel. Use it
//' when only the data set is big.
//'
//' @param pVec a parameter vector
//' @param data a model data instance
//' argument.
//' @return a double scalar
//' @export
//' @examples
//' m1 <- model.dmc(
//'   p.map     = list(a="1", v="1", z="1", d="1", sz="1", sv="1", t0="1", st0="1"),
//'   constants = c(st0=0, d=0),
//'   match.map = list(M=list(s1="r1", s2="r2")),
//'   factors   = list(S=c("s1", "s2")),
//'   responses = c("r1", "r2"),
//'   type      = "rd")
//'
//' pVec <- c(a=1, v=1, z=0.5, sz=0.25, sv=0.2,t0=.15)
//'
//' ## Set up a model-data instance
//' dat <- simulate(m1, nsim=1e2, p.vector=pVec)
//' mdi <- data.model.dmc(dat, m1)
//' summed_log_likelihood(pVec, mdi)
//' ## [1] 0.3796048
// [[Rcpp::export]]
double summed_log_likelihood(arma::vec& pVec, Rcpp::List& data) {
    // TODO: check pVec_NV and pVecNA names match
    // Should find Rcpp sugar match
  Rcpp::NumericVector model    = data.attr("model") ;
  Rcpp::NumericVector pVecNA   = model.attr("p.vector") ;
  Rcpp::CharacterVector pNames = pVecNA.names() ;
  Rcpp::NumericVector pVec_NV  = Rcpp::as<Rcpp::NumericVector>(Rcpp::wrap(pVec)) ;
  pVec_NV.names()              = pNames ;
  std::vector<double> density  = ddmc(data, pVec_NV, 2.5, 1e-10) ;
  double summedLogLike         = 0 ;

  for(std::vector<double>::iterator it=density.begin(); it!=density.end();
  ++it)
  { summedLogLike += log(*it) ; }
  return summedLogLike ;
}


//' @rdname summed_log_likelihood
//' @export
// [[Rcpp::export]]
double summed_log_likelihood_parallel(arma::vec& pVec, Rcpp::List& data) {
  Rcpp::NumericVector model    = data.attr("model") ;
  Rcpp::NumericVector pVecNA   = model.attr("p.vector") ;
  Rcpp::CharacterVector pNames = pVecNA.names() ;
  Rcpp::NumericVector pVec_NV  = Rcpp::as<Rcpp::NumericVector>(Rcpp::wrap(pVec)) ;
  pVec_NV.names()              = pNames ;
  double summedLogLike         = 0 ;

  std::vector<double> density = ddmc_parallel(data, pVec_NV, 2.5, 1e-10) ;

  for(std::vector<double>::iterator it=density.begin(); it!=density.end();
  ++it)
  {
    summedLogLike += log(*it) ;
  }
  return summedLogLike ;
}

double summed_log_likelihood_hyper(arma::mat& data_hyper,
  Rcpp::List& pList_hyper) {
  /* hyper log likelihood sums over all DDM parameters and over all
  participants replacing the location and scale parameters for each
  DDM parameter with the updated phi vectors
  * data_hyper is from data_hyper_ii, one chain from migrate subchains */
  Rcpp::CharacterVector pNames = pList_hyper.names() ;
  int nSubjects = data_hyper.n_rows ; // (nSubjects x npars)

  double sumOverSubjects = 0;
  for(int i=0; i<nSubjects; i++)
  {
    arma::mat pVec_i = data_hyper.row(i) ; // 1 x npars mat
    Rcpp::NumericVector pVecNV_i =
      Rcpp::as<Rcpp::NumericVector>(Rcpp::wrap(pVec_i)) ;
    pVecNV_i.names() = pNames ;
    double sumOverPars = sum(dprior(pVecNV_i, pList_hyper)) ;
    sumOverSubjects += sumOverPars ;
  }

  return sumOverSubjects ;
}


//' Sum and Log Prior Density of a EAM model
//'
//' Get log likelihood summed over all prior parameters
//'
//' @param pVec a parameter vector
//' @param pPrior p.prior
//' @export
//' @examples
//' pVec <- c(a=1, v=1, z=0.5, sz=0.25, sv=0.2,t0=.15)
//'
//' p.prior  <- prior.p.dmc(
//'   dists = rep("tnorm", 6),
//'   p1    = c(a=2,   v=2.5, z=0.5, sz=0.3, sv=1,  t0=0.3),
//'   p2    = c(a=0.5, v=.5,  z=0.1, sz=0.1, sv=.3, t0=0.05) * 5,
//'   lower = c(0,-5, 0, 0, 0, 0),
//'   upper = c(5, 7, 2, 2, 2, 2))
//'
//' summed_log_prior(pVec, p.prior)
//' ## -3.224406
// [[Rcpp::export]]
double summed_log_prior(arma::vec& pVec, Rcpp::List& pPrior) {
  Rcpp::NumericVector pVec_NV  = Rcpp::as<Rcpp::NumericVector>(Rcpp::wrap(pVec)) ;
  std::vector<std::string> pNames = pPrior.names() ;
  pVec_NV.names() = pNames ;
  return sum(dprior(pVec_NV, pPrior)) ;
}

double summed_log_prior_hyper(Rcpp::List& pVec, Rcpp::List& ppList) {
  arma::mat phiVecLocation = pVec[0] ;
  arma::mat phiVecScale    = pVec[1] ;
  Rcpp::NumericVector location =
    Rcpp::as<Rcpp::NumericVector>(Rcpp::wrap(phiVecLocation)) ;
  Rcpp::NumericVector scale    =
    Rcpp::as<Rcpp::NumericVector>(Rcpp::wrap(phiVecScale)) ;
  Rcpp::List ppList_location            = ppList[0] ;
  Rcpp::List ppList_scale               = ppList[1] ;
  Rcpp::CharacterVector pNames_location = ppList_location.names() ;
  Rcpp::CharacterVector pNames_scale    = ppList_scale.names() ;
  location.names() = pNames_location ;
  scale.names()    = pNames_scale ;
  double out = sum(dprior(location, ppList_location)) +
    sum(dprior(scale, ppList_scale)) ;
  return out ;
}

arma::vec update_useLogPrior(arma::mat& theta, Rcpp::List& pLists)
{
  int nChains = theta.n_rows ;   // nChains x npars
  arma::vec out(nChains) ;
  for(int i=0; i<nChains; i++)
  {
    arma::vec pVec   = vectorise(theta.row(i)) ;
    Rcpp::List pList = pLists[i] ; // new priors
    out(i)           = summed_log_prior(pVec, pList) ;
  }
  return out ;
}

