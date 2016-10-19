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
#include <ggdmc.hpp>

class sampler {

private:

  /* use "double gammaMult" to differentiate migrate and crossover
   * migrate=migrate_primitive, migrate.h=migrate_data, h.migrate=migrate_hyper */
  inline std::vector<int> pick_2chains (int k, Rcpp::IntegerVector chains) {
    chains.erase(chains.begin() + k) ;  // Remove current chain; ie [-k] R-C
    std::random_shuffle(chains.begin(), chains.end());
    Rcpp::IntegerVector tmp = chains[Rcpp::Range(0, 1)] ; // ie index
    std::vector<int>pickedChains(2) ;
    pickedChains[0] = tmp[0] - 1 ; // R index corrected to C index
    pickedChains[1] = tmp[1] - 1 ;
    return pickedChains ;
  } ;

  inline std::vector<int> get_subchains (int nChains) {
    Rcpp::IntegerVector chainSeq = Rcpp::seq_len(nChains) ; /* eg 1:24 */
  std::random_shuffle(chainSeq.begin(), chainSeq.end()) ;
  int nSubchains = chainSeq[0] ; // how many groups to work with; eg 3

  // shuffle again; lnum2; 3, 4, 9
  std::random_shuffle(chainSeq.begin(), chainSeq.end());
  Rcpp::IntegerVector subchains = chainSeq[Rcpp::Range(0, nSubchains-1)] ;
  std::vector<int> out(nSubchains) ;
  std::sort(subchains.begin(), subchains.end()) ;

  // Convert to C index
  for(Rcpp::IntegerVector::iterator it=subchains.begin(); it!=subchains.end();
  it++)
  {
    int idx  = std::distance(subchains.begin(), it) ;
    out[idx] = *it - 1;
  }
  return out ;
  } ;

  inline std::vector<int> shuffle_chains (int nChains) {
        Rcpp::IntegerVector chainSeq = Rcpp::seq_len(nChains) ; /* eg 1:24 */
        std::random_shuffle(chainSeq.begin(), chainSeq.end()) ;
        std::vector<int> out(nChains) ;

        // Convert to C index
        for(Rcpp::IntegerVector::iterator it=chainSeq.begin(); it!=chainSeq.end();
            it++)
        {
            int idx  = std::distance(chainSeq.begin(), it) ;
            out[idx] = *it - 1;
        }
        return out ;
    } ;

public:
  arma::mat migrate_primitive(arma::mat& useTheta,
    arma::vec& useLogPrior, arma::vec& useLogLike,
    Rcpp::List& pPrior, Rcpp::List& data, double rp) {
    /* Oizuki */
    int nChains = useTheta.n_rows ;
    int npars   = useTheta.n_cols ;
    std::vector<int> subchains    = get_subchains(nChains) ; // C index
    int nSubchains                = subchains.size() ;

    /* Gyakuzuki the useTheta ------------------------------------------------- */
    arma::mat thetaSet(nSubchains, npars) ;
    arma::vec currentLogPrior(nSubchains) ;
    arma::vec proposeLogPrior(nSubchains) ;
    arma::vec currentLogLike(nSubchains) ;
    arma::vec proposeLogLike(nSubchains) ;

    arma::vec currwLogLike(nSubchains) ;
    arma::vec propwLogLike(nSubchains) ;
    arma::vec currwLogPrior(nSubchains) ;
    arma::vec propwLogPrior(nSubchains) ;

    Rcpp::NumericVector model       = data.attr("model") ;
    Rcpp::NumericVector pVecNA      = model.attr("p.vector") ;
    thetaSet.fill(NA_REAL) ;

    for(int i=0; i<nSubchains; i++)
    {
      arma::rowvec perturbation  = Rcpp::runif(npars, -rp, rp) ;
      // create a set of particles to swap
      int ii              = subchains[i] ; // ii is non-continuous chain index
      thetaSet.row(i)     = useTheta.row(ii) + perturbation ;
      currentLogPrior(i)  = useLogPrior(ii) ; // nChains x 1
      currentLogLike(i)   = useLogLike(ii) ;
      arma::vec pVec      = vectorise(thetaSet.row(i)) ; // Proposed new prior
      proposeLogPrior(i)  = summed_log_prior(pVec, pPrior) ;
      proposeLogLike(i)   = summed_log_likelihood(pVec, data) ;

      propwLogPrior(i) = proposeLogPrior(i) ;
      propwLogLike(i)  = proposeLogLike(i) ;

      currwLogPrior(i) = currentLogPrior(i) ;
      currwLogLike(i)  = currentLogLike(i) ;
    }

    /* ppLogLike stands for "proposed posterior log likelihood".
    * cpLogLike stands for "current  posterior log likelihood". */
    double ppLogLike = proposeLogLike(nSubchains-1) + proposeLogPrior(nSubchains-1) ;
    double cpLogLike = currentLogLike(0) + currentLogPrior(0) ;
    if (std::isnan(ppLogLike)) { ppLogLike = -INFINITY ; }
    double rho = exp(ppLogLike - cpLogLike) ;

    if (rho > Rf_runif(0, 1)) {
      useTheta.row(subchains[0]) = thetaSet.row(nSubchains-1) ;
      useLogPrior(subchains[0])  = proposeLogPrior(nSubchains-1) ;
      useLogLike(subchains[0])   = proposeLogLike(nSubchains-1) ;
    }

    if ( nSubchains != 1 ) {
      for(int k=0; k<(nSubchains-2); k++)  // i & j were used before
      {
        // If the current selected chain is more probable than its follower,
        // replace its follower with the current chain.
        double ppLogLike = proposeLogLike(k) + proposeLogPrior(k) ;
        double cpLogLike = currentLogLike(k+1) + currentLogPrior(k+1) ;
        if (std::isnan(ppLogLike)) { ppLogLike = -INFINITY ; }
        double rho = exp (ppLogLike - cpLogLike) ;

        if ( rho > Rf_runif(0, 1) )
        {
          useTheta.row(subchains[k+1]) = thetaSet.row(k) ;
          useLogPrior(subchains[k+1])  = proposeLogPrior(k) ;
          useLogLike(subchains[k+1])   = proposeLogLike(k) ;
        }
      }
    }

    arma::mat out (nChains, 2+npars);
    out.col(0) = useLogPrior ;
    out.col(1) = useLogLike ;
    out.cols(2, 1+npars) = useTheta ;

    return out ;
  } ;

  arma::mat migrate_primitive_parallel(arma::mat& useTheta,
    arma::vec& useLogPrior, arma::vec& useLogLike,
    Rcpp::List& pPrior, Rcpp::List& data, double rp) {
    /* Oizuki */
    int nChains = useTheta.n_rows ;
    int npars   = useTheta.n_cols ;
    std::vector<int> subchains    = get_subchains(nChains) ; // C index
    int nSubchains                = subchains.size() ;

    /* Gyakuzuki the useTheta ------------------------------------------------- */
    arma::mat thetaSet(nSubchains, npars) ;
    arma::vec currentLogPrior(nSubchains) ;
    arma::vec proposeLogPrior(nSubchains) ;
    arma::vec currentLogLike(nSubchains) ;
    arma::vec proposeLogLike(nSubchains) ;

    arma::vec currwLogLike(nSubchains) ;
    arma::vec propwLogLike(nSubchains) ;
    arma::vec currwLogPrior(nSubchains) ;
    arma::vec propwLogPrior(nSubchains) ;

    Rcpp::NumericVector model       = data.attr("model") ;
    Rcpp::NumericVector pVecNA      = model.attr("p.vector") ;
    thetaSet.fill(NA_REAL) ;

    for(int i=0; i<nSubchains; i++)
    {
      arma::rowvec perturbation  = Rcpp::runif(npars, -rp, rp) ;
      // create a set of particles to swap
      int ii              = subchains[i] ; // ii is non-continuous chain index
      thetaSet.row(i)     = useTheta.row(ii) + perturbation ;
      currentLogPrior(i)  = useLogPrior(ii) ; // nChains x 1
      currentLogLike(i)   = useLogLike(ii) ;
      arma::vec pVec      = vectorise(thetaSet.row(i)) ; // Proposed new prior
      proposeLogPrior(i)  = summed_log_prior(pVec, pPrior) ;
      proposeLogLike(i)   = summed_log_likelihood_parallel(pVec, data) ;

      propwLogPrior(i) = proposeLogPrior(i) ;
      propwLogLike(i)  = proposeLogLike(i) ;

      currwLogPrior(i) = currentLogPrior(i) ;
      currwLogLike(i)  = currentLogLike(i) ;
    }

    /* ppLogLike stands for "proposed posterior log likelihood".
     * cpLogLike stands for "current  posterior log likelihood". */
    double ppLogLike = proposeLogLike(nSubchains-1) + proposeLogPrior(nSubchains-1) ;
    double cpLogLike = currentLogLike(0) + currentLogPrior(0) ;
    if (std::isnan(ppLogLike)) { ppLogLike = -INFINITY ; }
    double rho = exp(ppLogLike - cpLogLike) ;

    if (rho > Rf_runif(0, 1)) {
      useTheta.row(subchains[0]) = thetaSet.row(nSubchains-1) ;
      useLogPrior(subchains[0])  = proposeLogPrior(nSubchains-1) ;
      useLogLike(subchains[0])   = proposeLogLike(nSubchains-1) ;
    }

    if ( nSubchains != 1 ) {
      for(int k=0; k<(nSubchains-2); k++)  // i & j were used before
      {
        // If the current selected chain is more probable than its follower,
        // replace its follower with the current chain.
        double ppLogLike = proposeLogLike(k) + proposeLogPrior(k) ;
        double cpLogLike = currentLogLike(k+1) + currentLogPrior(k+1) ;
        if (std::isnan(ppLogLike)) { ppLogLike = -INFINITY ; }
        double rho = exp (ppLogLike - cpLogLike) ;

        if ( rho > Rf_runif(0, 1) )
        {
          useTheta.row(subchains[k+1]) = thetaSet.row(k) ;
          useLogPrior(subchains[k+1])  = proposeLogPrior(k) ;
          useLogLike(subchains[k+1])   = proposeLogLike(k) ;
        }
      }
    }

    arma::mat out (nChains, 2+npars);
    out.col(0) = useLogPrior ;
    out.col(1) = useLogLike ;
    out.cols(2, 1+npars) = useTheta ;

    return out ;
  } ;

  /* data level migrate with different priors for each chain, pLists=p.priors */
  arma::mat migrate_data(arma::mat& useTheta, // nChains x npars
    arma::vec& useLogPrior,                   // nChains x 1
    arma::vec& useLogLike,                    // nChains x 1
    Rcpp::List& pLists,
    Rcpp::List& data, double rp) {
    int nChains = useTheta.n_rows ; // the chain numbers
    int npars   = useTheta.n_cols ; // npars
    std::vector<int> subchains    = get_subchains(nChains) ;
    int nSubchains                = subchains.size() ;

    /* Gyakuzuki the useTheta ------------------------------------------------- */
    arma::mat thetaSet(nSubchains, npars) ;
    arma::vec proposeLogPrior(nSubchains) ;
    arma::vec currentLogPrior(nSubchains) ;
    arma::vec proposeLogLike(nSubchains) ;
    arma::vec currentLogLike(nSubchains) ;
    Rcpp::NumericVector model    = data.attr("model") ;
    Rcpp::NumericVector pVecNA   = model.attr("p.vector") ;
    Rcpp::CharacterVector pNames = pVecNA.names() ;
    thetaSet.fill(NA_REAL);

    for(int i=0; i<nSubchains; i++)
    {
      arma::rowvec perturbation  = Rcpp::runif(npars, -rp, rp) ;
      int ii              = subchains[i] ; // ii is non-continuous chain index
      thetaSet.row(i)     = useTheta.row(ii) + perturbation ;
      currentLogPrior(i)  = useLogPrior(ii) ; // nChains x 1
      currentLogLike(i)   = useLogLike(ii) ;
      arma::vec pVec      = vectorise(thetaSet.row(i)) ; // Proposed new prior
      Rcpp::List pList_ii = pLists[ii] ;
      proposeLogPrior(i)  = summed_log_prior(pVec, pList_ii) ;
      proposeLogLike(i)   = summed_log_likelihood(pVec, data) ;
    }

    /* M-H algorithm.
    * - ppLogLike == "proposed posterior log likelihood".
    * - cpLogLike == "current  posterior log likelihood".
    * - nSubchains - 1 again, 'cause R-C index complex */
    double ppLogLike = proposeLogLike(nSubchains-1) + proposeLogPrior(nSubchains-1) ;
    double cpLogLike = currentLogLike(0) + currentLogPrior(0) ;
    //double logprior = arma::as_scalar(proposeLogPrior(nSubchains-1)) ;
    if (std::isnan(ppLogLike)) { ppLogLike = -INFINITY ; }
    //if (std::isnan(logprior)) { proposeLogPrior(nSubchains-1) = -INFINITY ; }

    double rho = exp(ppLogLike - cpLogLike) ;

    // useTheta nChains x npars
    //if ( std::isnormal(rho) && rho > Rf_runif(0, 1) ) {
    if (rho > Rf_runif(0, 1)) {
      useTheta.row(subchains[0]) = thetaSet.row(nSubchains-1) ;
      useLogPrior(subchains[0])  = proposeLogPrior(nSubchains-1) ;
      useLogLike(subchains[0])   = proposeLogLike(nSubchains-1) ;
    }
    if ( nSubchains != 1 ) {
      for(int k=0; k<(nSubchains-2); k++)  // i & j were used before
      {
        // If the current selected chain is more probable than its follower,
        // replace its follower with the current chain.
        double ppLogLike = proposeLogLike(k) + proposeLogPrior(k) ;
        double cpLogLike = currentLogLike(k+1) + currentLogPrior(k+1) ;
        if (std::isnan(ppLogLike)) { ppLogLike = -INFINITY ; }
        if (std::isnan(proposeLogPrior(k))) { proposeLogPrior(k) = -INFINITY ; }
        double rho = exp (ppLogLike - cpLogLike) ;

        if ( rho > Rf_runif(0, 1) )
        {
          useTheta.row(subchains[k+1]) = thetaSet.row(k) ;
          useLogPrior(subchains[k+1])  = proposeLogPrior(k) ;
          useLogLike(subchains[k+1])   = proposeLogLike(k) ;
        }
      }
    }

    arma::mat out (nChains, 2+npars);
    out.col(0) = useLogPrior ;
    out.col(1) = useLogLike ;
    out.cols(2, 1+npars) = useTheta ;
    return out ;
  } ;

  arma::mat migrate_data_parallel(arma::mat& useTheta, // nChains x npars
    arma::vec& useLogPrior,                   // nChains x 1
    arma::vec& useLogLike,                    // nChains x 1
    Rcpp::List& pLists,
    Rcpp::List& data, double rp) {
    int nChains = useTheta.n_rows ; // the chain numbers
    int npars   = useTheta.n_cols ; // npars
    std::vector<int> subchains    = get_subchains(nChains) ;
    int nSubchains                = subchains.size() ;

    /* Gyakuzuki the useTheta ------------------------------------------------- */
    arma::mat thetaSet(nSubchains, npars) ;
    arma::vec proposeLogPrior(nSubchains) ;
    arma::vec currentLogPrior(nSubchains) ;
    arma::vec proposeLogLike(nSubchains) ;
    arma::vec currentLogLike(nSubchains) ;
    Rcpp::NumericVector model    = data.attr("model") ;
    Rcpp::NumericVector pVecNA   = model.attr("p.vector") ;
    Rcpp::CharacterVector pNames = pVecNA.names() ;
    thetaSet.fill(NA_REAL);

    for(int i=0; i<nSubchains; i++)
    {
      arma::rowvec perturbation  = Rcpp::runif(npars, -rp, rp) ;
      int ii              = subchains[i] ; // ii is non-continuous chain index
      thetaSet.row(i)     = useTheta.row(ii) + perturbation ;
      currentLogPrior(i)  = useLogPrior(ii) ; // nChains x 1
      currentLogLike(i)   = useLogLike(ii) ;
      arma::vec pVec      = vectorise(thetaSet.row(i)) ; // Proposed new prior
      Rcpp::List pList_ii = pLists[ii] ;
      proposeLogPrior(i)  = summed_log_prior(pVec, pList_ii) ;
      proposeLogLike(i)   = summed_log_likelihood_parallel(pVec, data) ;
    }

    /* M-H algorithm.
    * - ppLogLike == "proposed posterior log likelihood".
    * - cpLogLike == "current  posterior log likelihood".
    * - nSubchains - 1 again, 'cause R-C index complex */
    double ppLogLike = proposeLogLike(nSubchains-1) + proposeLogPrior(nSubchains-1) ;
    double cpLogLike = currentLogLike(0) + currentLogPrior(0) ;
    //double logprior = arma::as_scalar(proposeLogPrior(nSubchains-1)) ;
    if (std::isnan(ppLogLike)) { ppLogLike = -INFINITY ; }
    //if (std::isnan(logprior)) { proposeLogPrior(nSubchains-1) = -INFINITY ; }

    double rho = exp(ppLogLike - cpLogLike) ;

    // useTheta nChains x npars
    //if ( std::isnormal(rho) && rho > Rf_runif(0, 1) ) {
    if (rho > Rf_runif(0, 1)) {
      useTheta.row(subchains[0]) = thetaSet.row(nSubchains-1) ;
      useLogPrior(subchains[0])  = proposeLogPrior(nSubchains-1) ;
      useLogLike(subchains[0])   = proposeLogLike(nSubchains-1) ;
    }
    if ( nSubchains != 1 ) {
      for(int k=0; k<(nSubchains-2); k++)  // i & j were used before
      {
        // If the current selected chain is more probable than its follower,
        // replace its follower with the current chain.
        double ppLogLike = proposeLogLike(k) + proposeLogPrior(k) ;
        double cpLogLike = currentLogLike(k+1) + currentLogPrior(k+1) ;
        if (std::isnan(ppLogLike)) { ppLogLike = -INFINITY ; }
        if (std::isnan(proposeLogPrior(k))) { proposeLogPrior(k) = -INFINITY ; }
        double rho = exp (ppLogLike - cpLogLike) ;

        if ( rho > Rf_runif(0, 1) )
        {
          useTheta.row(subchains[k+1]) = thetaSet.row(k) ;
          useLogPrior(subchains[k+1])  = proposeLogPrior(k) ;
          useLogLike(subchains[k+1])   = proposeLogLike(k) ;
        }
      }
    }

    arma::mat out (nChains, 2+npars);
    out.col(0) = useLogPrior ;
    out.col(1) = useLogLike ;
    out.cols(2, 1+npars) = useTheta ;
    return out ;
  } ;

  /* h.migrate; doing DEMCMC migrate set, all chains, hyper level */
  arma::mat migrate_hyper(Rcpp::List& usePhi,
    arma::vec& useLogPrior,
    arma::vec& useLogLike,
    Rcpp::List pList,  // definitively make a copy of pList
    arma::cube& data_hyper,
    Rcpp::List& ppList,
    double rp) {
    // de is a list with two vectors negative below
    // a  v.f1  v.f2     z    sz sv.f1 sv.f2    t0
    // TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE
    // a  v.f1  v.f2     z    sz sv.f1 sv.f2    t0
    // TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE
    arma::mat usePhiLocation   = usePhi[0] ; // nChains x npars
    arma::mat usePhiScale      = usePhi[1] ;
    int nChains                = usePhiLocation.n_rows ;
    int npars                  = usePhiLocation.n_cols ;
    std::vector<int> subchains = get_subchains(nChains) ;
    int nSubchains             = subchains.size() ;

    /* Gyakuzuki various containers--------------------- */
    arma::mat phiSetLocation(nSubchains, npars) ;
    arma::mat phiSetScale(nSubchains, npars)  ;
    arma::vec proposeSetLogPrior(nSubchains) ;
    arma::vec proposeSetLogLike(nSubchains) ;
    arma::vec currentSet(nSubchains) ;
    arma::vec proposeSet(nSubchains) ;
    phiSetLocation.fill(NA_REAL) ;
    phiSetScale.fill(NA_REAL) ;

    for(int i=0; i<nSubchains; i++) // 0, 1, 2
    {
      arma::rowvec perturbation  = Rcpp::runif(npars, -rp, rp) ;
      int ii = subchains[i] ; // eg 3, 4, 9
      phiSetLocation.row(i)    = usePhiLocation.row(ii) + perturbation ;
      phiSetScale.row(i)       = usePhiScale.row(ii) + perturbation ;
      arma::vec phiVecLocation = vectorise(phiSetLocation.row(i)) ;
      arma::vec phiVecScale    = vectorise(phiSetScale.row(i)) ;
      Rcpp::List pVec = Rcpp::List::create(
        Rcpp::Named("location") = phiVecLocation, // nSubchains x npars
        Rcpp::Named("scale")    = phiVecScale) ;
      arma::mat data_hyper_ii = data_hyper.slice(ii) ; // nSubjects x npars
      assign_pp_pList(pVec, pList) ; // pList_cp
      proposeSetLogPrior(i)   = summed_log_prior_hyper(pVec, ppList) ;
      proposeSetLogLike(i)    = summed_log_likelihood_hyper(data_hyper_ii, pList) ;
      proposeSet(i) = proposeSetLogPrior(i) + proposeSetLogLike(i) ;
      currentSet(i) = useLogPrior(ii) + useLogLike(ii) ;
    }

    double ppLogLike = proposeSet(nSubchains-1)  ;
    double cpLogLike = currentSet(0) ;
    if (std::isnan(ppLogLike)) { ppLogLike = -INFINITY ; }
    double rho = exp(ppLogLike - cpLogLike) ;
    if ( rho > Rf_runif(0, 1) ) {
      usePhiLocation.row(subchains[0]) = phiSetLocation.row(nSubchains-1) ;
      usePhiScale.row(subchains[0])    = phiSetScale.row(nSubchains-1) ;
      useLogPrior(subchains[0])        = proposeSetLogPrior(nSubchains-1) ;
      useLogLike(subchains[0])         = proposeSetLogLike(nSubchains-1) ;
    }

    /* Compare all selected chains, except the last ------------ */
    if( nSubchains != 1) // if the selected group contains more than 1 chain
    {
      for(int k=0; k<(nSubchains-2); k++)
      {
        double ppLogLike = proposeSet(k) ;
        double cpLogLike = currentSet(k+1) ;
        if (std::isnan(ppLogLike)) { ppLogLike = -INFINITY ; }
        double rho = exp(ppLogLike - cpLogLike) ;
        if (rho > Rf_runif(0, 1))
        {
          usePhiLocation.row(subchains[k+1]) = phiSetLocation.row(k) ;
          usePhiScale.row(subchains[k+1])    = phiSetScale.row(k) ;
          useLogPrior(subchains[k+1])        = proposeSetLogPrior(k);
          useLogLike(subchains[k+1])         = proposeSetLogLike(k);
        }
      }
    }

    /* Yasume: Construct a matrix to return-------------------- */
    arma::mat out (nChains, 2 + 2*npars);
    out.col(0)                   = useLogPrior ;
    out.col(1)                   = useLogLike ;
    out.cols(2, 1+npars)         = usePhiLocation ;
    out.cols(2+npars, 1+2*npars) = usePhiScale ; // Maybe this usePhiLocation
    return out ;                               //error
  } ;

  arma::mat crossover_primitive(arma::mat& useTheta,
    arma::vec& useLogPrior, // nChains x 1
    arma::vec& useLogLike,
    Rcpp::List& pList,
    Rcpp::List& data,
    double rp,
    double gammaMult) {
    int nChains = useTheta.n_rows ; // useTheta
    int npars   = useTheta.n_cols;  // npars
    double gamma = std::isnan(gammaMult) ? Rf_runif(0.5, 1.0) :
      gammaMult/sqrt(2*npars) ;
    arma::vec gamma_vec = Rcpp::rep(gamma, npars) ;  // make it 1 x 8
    arma::mat out(nChains, 2+npars) ;
    std::vector<int> shuffledChains = shuffle_chains(nChains) ;

    for(std::vector<int>::iterator it=shuffledChains.begin(); it!=shuffledChains.end();
        it++) {
        // for(int i=0; i<nChains; i++) {
      Rcpp::IntegerVector chainSeq  = Rcpp::seq_len(nChains) ; // from 1 (not 0)
      std::vector<int> pickedChains = pick_2chains(*it, chainSeq) ;

      arma::vec theta0 = vectorise(useTheta.row(*it)) ;
      arma::vec theta1 = vectorise(useTheta.row(pickedChains[0])) ;
      arma::vec theta2 = vectorise(useTheta.row(pickedChains[1])) ;
      arma::vec perturbation  = Rcpp::runif(npars, -rp, rp) ;
      // This iterator implements Storn and Price's (1995) DE1
      // % is armadillo's element-wise *
      arma::vec proposal_theta = theta0 + perturbation + gamma_vec % (theta1 - theta2) ;
      // Note No Rcpp::List pList = pLists[i] ;
      // If force == TRUE
      // double clogPrior = summed_log_prior(theta0, pList) ;
      // double clogLike  = summed_log_likelihood(theta0, data) ;
      // double cpLogLike = clogPrior + clogLike ;
      double cpLogLike =  useLogLike(*it) + useLogPrior(*it) ;

      double logPrior  = summed_log_prior(proposal_theta, pList) ;
      double logLike   = summed_log_likelihood(proposal_theta, data) ;
      double ppLogLike = logPrior + logLike ;

      /* Metropolis ratio=rho */
      if (std::isnan(ppLogLike)) { ppLogLike = -INFINITY; }
      double rho = exp(ppLogLike - cpLogLike) ;
      if ( Rf_runif(0, 1) < rho) {
        useTheta.row(*it) = proposal_theta.t() ;
        useLogPrior(*it)  = logPrior ;
        useLogLike(*it)   = logLike ;
      }

      /* Yasume -------------------- */
      arma::vec update_a_chain(2+npars) ;
      update_a_chain(0) = useLogPrior(*it) ;
      update_a_chain(1) = useLogLike(*it) ;
      update_a_chain.subvec(2, npars+1) = useTheta.row(*it).t() ; // nChains x npars
      out.row(*it) = update_a_chain.t() ;
    }
    return out ;
  } ;

  arma::mat crossover_primitive_parallel(arma::mat& useTheta,
    arma::vec& useLogPrior, // nChains x 1
    arma::vec& useLogLike,
    Rcpp::List& pList,
    Rcpp::List& data,
    double rp,
    double gammaMult) {
    int nChains = useTheta.n_rows ; // useTheta
    int npars   = useTheta.n_cols;  // npars
    double gamma = std::isnan(gammaMult) ? Rf_runif(0.5, 1.0) :
      gammaMult/sqrt(2*npars) ;
    arma::vec gamma_vec = Rcpp::rep(gamma, npars) ;  // make it 1 x 8
    arma::mat out(nChains, 2+npars) ;
    std::vector<int> shuffledChains = shuffle_chains(nChains) ;

    for(std::vector<int>::iterator it=shuffledChains.begin(); it!=shuffledChains.end();
    it++) {
      // for(int i=0; i<nChains; i++) {
      Rcpp::IntegerVector chainSeq  = Rcpp::seq_len(nChains) ; // from 1 (not 0)
      std::vector<int> pickedChains = pick_2chains(*it, chainSeq) ;

      arma::vec theta0 = vectorise(useTheta.row(*it)) ;
      arma::vec theta1 = vectorise(useTheta.row(pickedChains[0])) ;
      arma::vec theta2 = vectorise(useTheta.row(pickedChains[1])) ;
      arma::vec perturbation  = Rcpp::runif(npars, -rp, rp) ;
      // This iterator implements Storn and Price's (1995) DE1
      // % is armadillo's element-wise *
      arma::vec proposal_theta = theta0 + perturbation + gamma_vec % (theta1 - theta2) ;
      // Note No Rcpp::List pList = pLists[i] ;
      // If force == TRUE
      // double clogPrior = summed_log_prior(theta0, pList) ;
      // double clogLike  = summed_log_likelihood(theta0, data) ;
      // double cpLogLike = clogPrior + clogLike ;
      double cpLogLike =  useLogLike(*it) + useLogPrior(*it) ;

      double logPrior  = summed_log_prior(proposal_theta, pList) ;
      double logLike   = summed_log_likelihood_parallel(proposal_theta, data) ;
      double ppLogLike = logPrior + logLike ;

      /* Metropolis ratio=rho */
      if (std::isnan(ppLogLike)) { ppLogLike = -INFINITY; }
      double rho = exp(ppLogLike - cpLogLike) ;
      if ( Rf_runif(0, 1) < rho) {
        useTheta.row(*it) = proposal_theta.t() ;
        useLogPrior(*it)  = logPrior ;
        useLogLike(*it)   = logLike ;
      }

      /* Yasume -------------------- */
      arma::vec update_a_chain(2+npars) ;
      update_a_chain(0) = useLogPrior(*it) ;
      update_a_chain(1) = useLogLike(*it) ;
      update_a_chain.subvec(2, npars+1) = useTheta.row(*it).t() ; // nChains x npars
      out.row(*it) = update_a_chain.t() ;
    }
    return out ;
  } ;


  arma::mat crossover_primitive_original(arma::mat& useTheta,
    arma::vec& useLogPrior, // nChains x 1
    arma::vec& useLogLike,
    Rcpp::List& pList,
    Rcpp::List& data,
    double rp,
    double gammaMult) {
    int nChains = useTheta.n_rows ; // useTheta
    int npars   = useTheta.n_cols;  // npars
    double gamma = std::isnan(gammaMult) ? Rf_runif(0.5, 1.0) :
      gammaMult/sqrt(2*npars) ;
    arma::vec gamma_vec = Rcpp::rep(gamma, npars) ;  // make it 1 x 8
    arma::mat out(nChains, 2+npars) ;

    for(int i=0; i<nChains; i++) {
      Rcpp::IntegerVector chainSeq  = Rcpp::seq_len(nChains) ; // from 1 (not 0)
      std::vector<int> pickedChains = pick_2chains(i, chainSeq) ;

      arma::vec theta0 = vectorise(useTheta.row(i)) ;
      arma::vec theta1 = vectorise(useTheta.row(pickedChains[0])) ;
      arma::vec theta2 = vectorise(useTheta.row(pickedChains[1])) ;
      arma::vec perturbation  = Rcpp::runif(npars, -rp, rp) ;
      // This iterator implements Storn and Price's (1995) DE1
      // % is armadillo's element-wise *
      arma::vec proposal_theta = theta0 + perturbation + gamma_vec % (theta1 - theta2) ;
      // Note No Rcpp::List pList = pLists[i] ;
      // If force == TRUE
      // double clogPrior = summed_log_prior(theta0, pList) ;
      // double clogLike  = summed_log_likelihood(theta0, data) ;
      // double cpLogLike = clogPrior + clogLike ;
      double cpLogLike =  useLogLike(i) + useLogPrior(i) ;

      double logPrior  = summed_log_prior(proposal_theta, pList) ;
      double logLike   = summed_log_likelihood(proposal_theta, data) ;
      double ppLogLike = logPrior + logLike ;

      /* Metropolis ratio=rho */
      if (std::isnan(ppLogLike)) { ppLogLike = -INFINITY; }
      double rho = exp(ppLogLike - cpLogLike) ;
      if ( Rf_runif(0, 1) < rho) {
        useTheta.row(i) = proposal_theta.t() ;
        useLogPrior(i)  = logPrior ;
        useLogLike(i)   = logLike ;
      }

      /* Yasume -------------------- */
      arma::vec update_a_chain(2+npars) ;
      update_a_chain(0) = useLogPrior(i) ;
      update_a_chain(1) = useLogLike(i) ;
      update_a_chain.subvec(2, npars+1) = useTheta.row(i).t() ; // nChains x npars
      out.row(i) = update_a_chain.t() ;
    }
    return out ;
  } ;


  /* This is crossover.h */
  arma::mat crossover_data(arma::mat& useTheta, // nChains x npars
    arma::vec& useLogPrior, // nChains x 1
    arma::vec& useLogLike,  // nChains x 1
    Rcpp::List& pLists,
    Rcpp::List& data,
    double rp,
    double gammaMult) {
    int nChains  = useTheta.n_rows ;
    int npars    = useTheta.n_cols ;
    double gamma = std::isnan(gammaMult) ? Rf_runif(0.5, 1.0) :
      gammaMult/std::sqrt(2*npars) ;
    arma::vec gamma_vec = Rcpp::rep(gamma, npars) ;
    arma::mat out (nChains, 2 + npars) ;
    std::vector<int> shuffledChains = shuffle_chains(nChains) ;

    for(std::vector<int>::iterator it=shuffledChains.begin(); it!=shuffledChains.end();
    it++) {
    //for(int i=0; i<nChains; i++) {
      Rcpp::IntegerVector chainSeq = Rcpp::seq_len(nChains) ; // from 1 (not 0)
      arma::vec theta0 = vectorise(useTheta.row(*it)) ;
      std::vector<int> pickedChains = pick_2chains(*it, chainSeq) ;
      arma::vec theta1 = vectorise(useTheta.row(pickedChains[0])) ;
      arma::vec theta2 = vectorise(useTheta.row(pickedChains[1])) ;
      arma::vec perturbation  = Rcpp::runif(npars, -rp, rp) ;
      // Update mu. ie chain0
      arma::vec proposal_theta = theta0 + perturbation +
        gamma_vec % (theta1 - theta2) ;

      Rcpp::List pList = pLists[*it] ; // !!!The only different from primitive!!!
      double logPrior  = summed_log_prior(proposal_theta, pList) ;
      double logLike   = summed_log_likelihood(proposal_theta, data) ;
      double ppLogLike = logPrior + logLike ;
      double cpLogLike = useLogPrior(*it) + useLogLike(*it) ;

      // Metropolis step
      if (std::isnan(ppLogLike)) { ppLogLike = -INFINITY ; }
      double rho = exp(ppLogLike - cpLogLike) ;
      //if (std::isnormal(rho) && rho > Rf_runif(0, 1)) {
      if (rho > Rf_runif(0, 1)) {
        useTheta.row(*it) = proposal_theta.t() ; // 8 x 1 to 1x8
        useLogPrior(*it)  = logPrior ;
        useLogLike(*it)   = logLike ;
      }

      /* Yasume -------------------- */
      arma::vec update_a_chain (2+npars) ; // column vector
      update_a_chain(0) = useLogPrior(*it) ;
      update_a_chain(1) = useLogLike(*it) ;
      update_a_chain.rows(2, npars+1) = useTheta.row(*it).t() ;
      out.row(*it) = update_a_chain.t() ; // out: nChains, 2 + npars
    }
    return out ;
  } ;

  arma::mat crossover_data_parallel(arma::mat& useTheta, // nChains x npars
    arma::vec& useLogPrior, // nChains x 1
    arma::vec& useLogLike,  // nChains x 1
    Rcpp::List& pLists,
    Rcpp::List& data,
    double rp,
    double gammaMult) {
    int nChains  = useTheta.n_rows ;
    int npars    = useTheta.n_cols ;
    double gamma = std::isnan(gammaMult) ? Rf_runif(0.5, 1.0) :
      gammaMult/std::sqrt(2*npars) ;
    arma::vec gamma_vec = Rcpp::rep(gamma, npars) ;
    arma::mat out (nChains, 2 + npars) ;
    std::vector<int> shuffledChains = shuffle_chains(nChains) ;

    for(std::vector<int>::iterator it=shuffledChains.begin(); it!=shuffledChains.end();
    it++) {
    //for(int i=0; i<nChains; i++) {
      Rcpp::IntegerVector chainSeq = Rcpp::seq_len(nChains) ; // from 1 (not 0)
      arma::vec theta0 = vectorise(useTheta.row(*it)) ;
      std::vector<int> pickedChains = pick_2chains(*it, chainSeq) ;
      arma::vec theta1 = vectorise(useTheta.row(pickedChains[0])) ;
      arma::vec theta2 = vectorise(useTheta.row(pickedChains[1])) ;
      arma::vec perturbation  = Rcpp::runif(npars, -rp, rp) ;
      // Update mu. ie chain0
      arma::vec proposal_theta = theta0 + perturbation +
        gamma_vec % (theta1 - theta2) ;

      Rcpp::List pList = pLists[*it] ; // !!!The only different from primitive!!!
      double logPrior  = summed_log_prior(proposal_theta, pList) ;
      double logLike   = summed_log_likelihood_parallel(proposal_theta, data) ;
      double ppLogLike = logPrior + logLike ;
      double cpLogLike = useLogPrior(*it) + useLogLike(*it) ;

      // Metropolis step
      if (std::isnan(ppLogLike)) { ppLogLike = -INFINITY ; }
      double rho = exp(ppLogLike - cpLogLike) ;
      //if (std::isnormal(rho) && rho > Rf_runif(0, 1)) {
      if (rho > Rf_runif(0, 1)) {
        useTheta.row(*it) = proposal_theta.t() ; // 8 x 1 to 1x8
        useLogPrior(*it)  = logPrior ;
        useLogLike(*it)   = logLike ;
      }

      /* Yasume -------------------- */
      arma::vec update_a_chain (2+npars) ; // column vector
      update_a_chain(0) = useLogPrior(*it) ;
      update_a_chain(1) = useLogLike(*it) ;
      update_a_chain.rows(2, npars+1) = useTheta.row(*it).t() ;
      out.row(*it) = update_a_chain.t() ; // out: nChains, 2 + npars
    }
    return out ;
  } ;

  /* This is h.crossover inside blocked.h.crossover */
  arma::mat crossover_hyper(Rcpp::List& usePhi,
    arma::vec& useLogPrior, // nChains x 1
    arma::vec& useLogLike,  // nChains x 1
    Rcpp::List pList,
    arma::cube& data_hyper, // nSubjects x npars x nChains
    Rcpp::List& ppList,
    double rp, double gammaMult)
  {
    arma::mat usePhiLocation = usePhi[0] ; // nChains x npars
    arma::mat usePhiScale    = usePhi[1] ; // nChains x npars
    int nChains              = usePhiLocation.n_rows ;
    int npars                = usePhiLocation.n_cols ;

    // Brandon's tunning parameter (gamma) page 372; step size
    double hgamma = std::isnan(gammaMult) ? Rf_runif(0.5, 1.0) :
      gammaMult/std::sqrt(2*npars*2) ; // extra *2 as p1 and p2
    arma::vec hgamma_vec = Rcpp::rep(hgamma, npars) ;
    arma::mat out (nChains, 2 + 2*npars) ;
    std::vector<int> shuffledChains = shuffle_chains(nChains) ;

    for(std::vector<int>::iterator it=shuffledChains.begin(); it!=shuffledChains.end();
    it++) {
    //for(int i=0; i<nChains; i++) {
      Rcpp::IntegerVector chainSeq = Rcpp::seq_len(nChains) ; // from 1 (not 0)
      // Update use.loglike for new ps, nChains x npars
      arma::vec usePhiLocation0 = vectorise(usePhiLocation.row(*it)) ;
      arma::vec usePhiScale0    = vectorise(usePhiScale.row(*it)) ;
      Rcpp::List pVec = Rcpp::List::create(
        Rcpp::Named("location") = usePhiLocation0,
        Rcpp::Named("scale")    = usePhiScale0) ;
      // data_hyper_k is nSubjects x npars
      arma::mat data_hyper_k = data_hyper.slice(*it) ;
      assign_pp_pList(pVec, pList) ; // use pVec to modifiy pList

      // Use chain i's hyper data to update useLogLike
      // Get updated cpLogLike based on updated useLogLike and current logPrior
      useLogLike(*it) = summed_log_likelihood_hyper(data_hyper_k, pList) ;
      double cpLogLike =  useLogPrior(*it) + useLogLike(*it) ;

      // DE step: 1. Pick two other chains; usePhi_location: nChains x npars
      std::vector<int> pickedChains = pick_2chains(*it, chainSeq) ;

      arma::vec phi_location1 = vectorise(usePhiLocation.row(pickedChains[0])) ;
      arma::vec phi_location2 = vectorise(usePhiLocation.row(pickedChains[1])) ;
      arma::vec phi_scale1    = vectorise(usePhiScale.row(pickedChains[0])) ;
      arma::vec phi_scale2    = vectorise(usePhiScale.row(pickedChains[1])) ;
      arma::vec perturbation  = Rcpp::runif(npars, -rp, rp) ;

      // Update mu and sigma
      arma::vec proposal_location = usePhiLocation0 + perturbation +
        hgamma_vec % (phi_location1 - phi_location2) ;
      arma::vec proposal_scale = usePhiScale0 + perturbation +
        hgamma_vec % (phi_scale1 - phi_scale2) ;
      Rcpp::List proposal_pVec  = Rcpp::List::create(
        Rcpp::Named("location") = proposal_location,
        Rcpp::Named("scale")    = proposal_scale) ;
      assign_pp_pList(proposal_pVec, pList) ;       // modify pList again

      // Get new post
      double logPrior  = summed_log_prior_hyper(proposal_pVec, ppList) ;
      double logLike   = summed_log_likelihood_hyper(data_hyper_k, pList) ;
      double ppLogLike = logPrior + logLike ;

      // Metropolis step; nChains x npars
      if (std::isnan(ppLogLike)) { ppLogLike = -INFINITY ; }
      double rho = exp(ppLogLike - cpLogLike) ;
      //if (std::isnormal(rho) && rho > Rf_runif(0, 1)) {
      if (rho > Rf_runif(0, 1)) {
        usePhiLocation.row(*it) = proposal_location.t() ; // 8x1 to 1x8
        usePhiScale.row(*it)    = proposal_scale.t() ;
        useLogPrior(*it)        = logPrior ;
        useLogLike(*it)         = logLike ;
      }

      /* usePhi_location: nChains x npars -------------------- */
      arma::vec update_a_chain (2+2*npars) ;
      update_a_chain(0) = useLogPrior(*it) ;
      update_a_chain(1) = useLogLike(*it) ;
      update_a_chain.rows(2, npars+1)         = usePhiLocation.row(*it).t() ;
      update_a_chain.rows(2+npars, 2*npars+1) = usePhiScale.row(*it).t() ;
      out.row(*it) = update_a_chain.t() ;
    }
    return out ;
  } ;

} ;

//' Run a Bayesian EAM Model for Fixed-effect or Random-effect
//'
//' These are the C++ functions hierarchical Bayesian models. The user usually
//' should not call this function directly. Use \code{h.run.dmc} or
//' \code{run.dmc} wrapper instead.  Two parallel versions of \code{run} uses
//' Open MPI to calculate diffusion probability density. Use them by setting
//' \code{cores} greater than 2 in \code{h.run.dmc}
//'
//' @param samples a DMC sample/object
//' @param setting a container to store all settings for running model fitting
//' @param debug a debugging switch to test chain randomisation. Default false
//' @return a DMC sample
//' @export
// [[Rcpp::export]]
Rcpp::List run_data(Rcpp::List& samples, Rcpp::List& setting, bool debug=false) {
  /* Extract samples arguments ---------------------------------------*/
  Rcpp::List samples_in(clone(samples)) ; // so R' original samples stays
  int nThins                   = samples_in["thin"];
  int nmc                      = samples_in["nmc"] ;
  int start_R                  = samples_in["start"] ;
  int start_C                  = start_R - 1 ;
  int nsamp                    = 1 + (nmc - start_R) * nThins ;
  double rp                    = samples_in["rp"] ;
  Rcpp::List pList             = samples_in["p.prior"];
  Rcpp::List data              = samples_in["data"] ;
  Rcpp::NumericVector RT       = data["RT"] ;
  Rcpp::NumericVector model    = data.attr("model") ;
  Rcpp::NumericVector pVecNA   = model.attr("p.vector") ;
  std::vector<std::string> pNames = pVecNA.names() ;

  /* Zenkutsu-dachi: Extract Bayes data------------------------------------*/
  arma::cube theta      = samples_in["theta"] ; // nChains x npars x nmc
  arma::mat logPrior    = samples_in["summed_log_prior"] ; // nmc x nChains
  arma::mat logLike     = samples_in["log_likelihoods"] ;  // nmc x nChains
  arma::mat useTheta    = theta.slice(start_C) ; // nChains x npars
  arma::vec useLogPrior = vectorise(logPrior.row(start_C)) ; // nChains x 1
  arma::vec useLogLike  = vectorise(logLike.row(start_C)) ;  // nChains x 1
  int npars             = pNames.size() ;
  int store_i           = start_C ;
  double pMigrate       = setting["p.migrate"] ;
  double gammaMult      = setting["gamma.mult"] ;
  int report            = setting["report"] ;
  arma::mat evolution ;   // nChains x 2+npars;
  sampler DE ;

  // Start from 2nd nsamp, so i=1, instead of 0 (C index)
  for(int i=1; i<nsamp; i++) {
    if (Rf_runif(0, 1) < pMigrate) {
      evolution = DE.migrate_primitive(useTheta, useLogPrior, useLogLike,
        pList, data, rp) ;
    } else if (debug==true) {
        evolution = DE.crossover_primitive_original(useTheta, useLogPrior, useLogLike,
                                           pList, data, rp, gammaMult) ;
    } else {
      evolution = DE.crossover_primitive(useTheta, useLogPrior, useLogLike,
        pList, data, rp, gammaMult) ;
    }

    useLogPrior = evolution.col(0) ; // nChains x 1
    useLogLike  = evolution.col(1) ;
    useTheta    = evolution.cols(2, 1+npars) ;

    if (i % nThins == 0) {
      store_i++ ;
      if ((store_i+1) % report == 0) { Rcpp::Rcout << store_i+1 << " "; }

      // arma::rowvec abc = logPrior.row(store_i) ;
      // for (int j=1; j<abc.n_elem; j++)
      // {
      //   double aaa = abc(j);
      //   Rcpp::Rcout <<  aaa << "  ";
      // }

      logPrior.row(store_i) = useLogPrior.t() ;  // nmc x nChains

      // std::cout <<"After" << std::endl;
      // arma::rowvec abc2 = logPrior.row(store_i) ;
      // for (int j=1; j<abc2.n_elem; j++)
      // {
      //   double aaa2 = abc2(j);
      //   Rcpp::Rcout <<  aaa2 << "  ";
      // }


      logLike.row(store_i)  = useLogLike.t() ;   // nmc x nChains
      theta.slice(store_i)  = useTheta ;
    }
  }

  /* ------------------Output ------------------------------------- */
  samples_in["summed_log_prior"] = logPrior ; // nmc x nChains
  samples_in["log_likelihoods"]  = logLike ;  // nmc x nChains
  samples_in["theta"]            = theta ;    // nChains x npars x nmc
  return samples_in ;
}


//' @rdname run_data
//' @export
// [[Rcpp::export]]
Rcpp::List run_data_parallel(Rcpp::List& samples, Rcpp::List& setting,
  bool debug=false) {
  /* Extract samples arguments ---------------------------------------*/
  Rcpp::List samples_in(clone(samples)) ; // so R' original samples stays
  int nThins                   = samples_in["thin"];
  int nmc                      = samples_in["nmc"] ;
  int start_R                  = samples_in["start"] ;
  int start_C                  = start_R - 1 ;
  int nsamp                    = 1 + (nmc - start_R) * nThins ;
  double rp                    = samples_in["rp"] ;
  Rcpp::List pList             = samples_in["p.prior"];
  Rcpp::List data              = samples_in["data"] ;
  Rcpp::NumericVector RT       = data["RT"] ;
  Rcpp::NumericVector model    = data.attr("model") ;
  Rcpp::NumericVector pVecNA   = model.attr("p.vector") ;
  std::vector<std::string> pNames = pVecNA.names() ;

  /* Zenkutsu-dachi: Extract Bayes data------------------------------------*/
  arma::cube theta      = samples_in["theta"] ; // nChains x npars x nmc
  arma::mat logPrior    = samples_in["summed_log_prior"] ; // nmc x nChains
  arma::mat logLike     = samples_in["log_likelihoods"] ;  // nmc x nChains
  arma::mat useTheta    = theta.slice(start_C) ; // nChains x npars
  arma::vec useLogPrior = vectorise(logPrior.row(start_C)) ; // nChains x 1
  arma::vec useLogLike  = vectorise(logLike.row(start_C)) ;  // nChains x 1
  int npars             = pNames.size() ;
  int store_i           = start_C ;
  double pMigrate       = setting["p.migrate"] ;
  double gammaMult      = setting["gamma.mult"] ;
  int report            = setting["report"] ;
  // int cores             = setting["cores"] ;  // leave for openmp
  arma::mat evolution ;   // nChains x 2+npars;
  sampler DE ;

  // Start from 2nd nsamp, so i=1, instead of 0 (C index)
  for(int i=1; i<nsamp; i++) {
    if (Rf_runif(0, 1) < pMigrate) {
      evolution = DE.migrate_primitive_parallel(useTheta, useLogPrior, useLogLike,
        pList, data, rp) ;
    } else {
      //Rcpp::Rcout << "double shuffle" << std::endl;
      evolution = DE.crossover_primitive_parallel(useTheta, useLogPrior, useLogLike,
        pList, data, rp, gammaMult) ;
    }

    useLogPrior = evolution.col(0) ; // nChains x 1
    useLogLike  = evolution.col(1) ;
    useTheta    = evolution.cols(2, 1+npars) ;

    if (i % nThins == 0) {
      store_i++ ;
      if ((store_i+1) % report == 0) { Rcpp::Rcout << store_i+1 << " "; }
      logPrior.row(store_i) = useLogPrior.t() ;  // nmc x nChains
      logLike.row(store_i)  = useLogLike.t() ;   // nmc x nChains
      theta.slice(store_i)  = useTheta ;

    }
  }


  /* ------------------Output ------------------------------------- */

  samples_in["summed_log_prior"] = logPrior ; // nmc x nChains
  samples_in["log_likelihoods"]  = logLike ;  // nmc x nChains
  samples_in["theta"]            = theta ;    // nChains x npars x nmc
  return samples_in ;
}

//' @rdname run_data
//' @export
// [[Rcpp::export]]
Rcpp::List run_hyper(Rcpp::List samples, Rcpp::List& setting) {
  /* Yoi: clone makes sure the samples is not altered----------------- */
  Rcpp::List samples_in(clone(samples)) ;
  Rcpp::List hyper, blocks ;
  double pMigrate, pMigrate_hyper, gammaMult ;
  hyper          = samples_in.attr("hyper") ;
  int npars      = hyper["n.pars"] ;
  int nmc        = hyper["nmc"] ;
  int nThins     = hyper["thin"] ;
  double rp      = hyper["rp"] ;    // rp is defined in initialise
  int start_R    = hyper["start"] ; // start_R == 1;
  //Rcpp::Rcout << "start_R " << start_R << std::endl;
  int start_C    = start_R - 1 ;    // start_C == 0;
  //Rcpp::Rcout << "start_C " << start_C << std::endl;
  int nsamp      = 1 + (nmc - start_R) * nThins ;
  //Rcpp::Rcout << "nsamp " << nsamp << std::endl;
  pMigrate       = setting["p.migrate"] ;
  pMigrate_hyper = setting["h.p.migrate"] ;
  gammaMult      = setting["gamma.mult"] ;
  int report     = setting["report"] ;

  /* Kihon: Get hyper-level data.
  * First iteration data are from initialised 1st nmc, collected from
  the theta per chain per participant. The C++'s data_hyper mirrors
  (nSubjects x npars x nChains) R's "ps" array (nChains x nSubjects x npars),
  'cause I use Armadillo chain/slice major. */
  Rcpp::List phi             = hyper["phi"] ;
  arma::mat hlogPrior        = hyper["h_summed_log_prior"] ; // nmc x nChains
  arma::mat hlogLike         = hyper["h_log_likelihoods"] ;  // nmc x nChains
  Rcpp::List ppList          = hyper["pp.prior"] ;
  arma::cube data_hyper      = get_ps(samples_in) ; // 1st nmc theta
  arma::cube phiLocation     = phi[0] ; // nChains x npars x nmc
  arma::cube phiScale        = phi[1] ;
  arma::vec usehLogPrior     = vectorise(hlogPrior.row(start_C)) ;
  arma::vec usehLogLike      = vectorise(hlogLike.row(start_C)) ;
  arma::mat usePhiLocation   = phiLocation.slice(start_C) ;
  arma::mat usePhiScale      = phiScale.slice(start_C) ; // nChains x npars
  Rcpp::List usePhi          = Rcpp::List::create(
    Rcpp::Named("location")  = usePhiLocation,
    Rcpp::Named("scale")     = usePhiScale) ;
  Rcpp::List ppList_location = ppList[0] ;
  Rcpp::List ppList_scale    = ppList[1] ;
  int store_i                = start_C ; // store_i == 0;
  Rcpp::List ppList_constant = isConstant(ppList) ; // two logical vectors

  /* Setup for data level, get "p.priors" and start taking samples */
  Rcpp::List samples_s1        = samples_in[0] ;
  Rcpp::List pList             = samples_s1["p.prior"] ;
  Rcpp::List pLists            = assign_pp_pLists(samples_s1, usePhi) ;
  Rcpp::CharacterVector pNames = pList.names() ;
  Rcpp::CharacterVector pars   = pNames ;
  add_nmc0(samples_in, pLists) ;   // Add "use" for data level
  arma::mat update_theta, evolution, evolution_data ;
  arma::vec update_logprior, update_loglike ;
  Rcpp::List updated_use ;
  sampler DE, DE_data ;

  for(int i=1; i<nsamp; i++) {
    /* Update usePhi, usehLogPrior, and usehLogLike */
    if ( Rf_runif(0, 1) < pMigrate_hyper ) {
      evolution = DE.migrate_hyper(usePhi, usehLogPrior, usehLogLike,
        pList, data_hyper, ppList, rp) ;
    } else {
      evolution = DE.crossover_hyper(usePhi, usehLogPrior, usehLogLike,
        pList, data_hyper, ppList, rp, gammaMult) ;
    }

    usehLogPrior              = evolution.col(0) ; // nChains x 1
    usehLogLike               = evolution.col(1) ; //
    arma::mat usePhiLocation  = evolution.cols(2, 1+npars) ;
    arma::mat usePhiScale     = evolution.cols(2+npars, 1+2*npars) ;

    usePhi                    = Rcpp::List::create(  // nChains x npars
      Rcpp::Named("location") = usePhiLocation,
      Rcpp::Named("scale")    = usePhiScale) ;
    pLists = assign_pp_pLists(samples_s1, usePhi) ;

    /*  Update data level and hyper data */
    for(Rcpp::List::iterator it=samples_in.begin(); it!=samples_in.end(); ++it)
    {
      Rcpp::List subject_j     = *it ;
      Rcpp::List data_j        = subject_j["data"] ;
      Rcpp::List use           = subject_j.attr("use") ;
      int store_j              = use["store_i"] ;
      arma::mat useTheta_j     = use["theta"] ; // nChains x npars
      arma::vec useLogPrior_j  = update_useLogPrior(useTheta_j, pLists) ;
      arma::vec useLogLike_j   = use["loglike"] ; // nChains
      if ( Rf_runif(0, 1) < pMigrate )
      {
        evolution_data = DE_data.migrate_data(useTheta_j, useLogPrior_j,
          useLogLike_j, pLists, data_j, rp) ;
      } else {
        evolution_data = DE_data.crossover_data(useTheta_j, useLogPrior_j,
          useLogLike_j, pLists, data_j, rp, gammaMult) ;
      }
      //arma::mat temp2 = shuffle(evolution_data, 0);
      update_logprior = evolution_data.col(0) ;
      update_loglike  = evolution_data.col(1) ;
      update_theta    = evolution_data.cols(2, 1+npars) ;

      // Harvest results; updated_use points at the same location as
      // attr(samples_i, "use")
      updated_use    = Rcpp::List::create(
        Rcpp::Named("theta")    = update_theta,    // nChains x npars
        Rcpp::Named("logprior") = update_logprior, // nChains x 1
        Rcpp::Named("loglike")  = update_loglike,  // nChains x 1
        Rcpp::Named("store_i")  = store_j,
        Rcpp::Named("p.priors") = pLists) ;
      subject_j.attr("use")  = updated_use ;

      if ( i % nThins == 0 ) {
        store_j++ ;  // store_ii starts from 0, nsamp starts from 1
        //Rcpp::Rcout << "\n" << store_ii << " . ";
        //if ((store_ii+1) % report == 0) { Rcpp::Rcout << store_ii+1 << " . "; }
        arma::mat tmp1       = subject_j["summed_log_prior"] ; // nmc x nChains
        arma::mat tmp2       = subject_j["log_likelihoods"] ;  // nmc x nChains
        arma::cube tmp3      = subject_j["theta"] ; // nChains x npars x nmc
        tmp1.row(store_j)    = update_logprior.t() ;
        tmp2.row(store_j)    = update_loglike.t() ;
        tmp3.slice(store_j)  = update_theta ;
        subject_j["summed_log_prior"] = tmp1 ;
        subject_j["log_likelihoods"]  = tmp2 ;
        subject_j["theta"]            = tmp3 ;
        // No need to feed back to subject_i !!!For DMC sake, still uses store_i
        updated_use["store_i"] = store_j ;
      }
    }

    data_hyper = get_ps_use(samples_in) ; // get new hyper data

    if ( i % nThins == 0 )
    {
      store_i++ ; // store_i starts from 0, nsamp starts from 1
      if ((store_i+1) % report == 0) { Rcpp::Rcout << store_i+1 << " "; }
      hlogPrior.row(store_i)      = usehLogPrior.t() ;  // nmc x nChains
      hlogLike.row(store_i)       = usehLogLike.t() ;   // nmc x nChains
      phiLocation.slice(store_i)  = evolution.cols(2, 1+npars) ; // nChains x npars x nmc
      phiScale.slice(store_i)     = evolution.cols(2+npars, 1+2*npars) ;
      //phiLocation.slice(store_i)  = usePhiLocation ;
      //phiScale.slice(store_i)     = usePhiScale ;
      Rcpp::List updated_phi      = Rcpp::List::create(
        Rcpp::Named("location")   = phiLocation,
        Rcpp::Named("scale")      = phiScale) ;
      hyper["h_log_likelihoods"]  = hlogLike ;
      hyper["h_summed_log_prior"] = hlogPrior ;
      hyper["phi"]                = updated_phi ;
      hyper["nmc"]                = nmc;
    }
  }

  /* Yame */
  Rcpp::Rcout << std::endl ;
  samples_in.attr("hyper") = hyper ;
  return samples_in ;
}

//' @rdname run_data
//' @export
// [[Rcpp::export]]
Rcpp::List run_hyper_parallel(Rcpp::List samples, Rcpp::List& setting) {
  /* Yoi: clone makes sure the samples is not altered----------------- */
  Rcpp::List samples_in(clone(samples)) ;
  Rcpp::List hyper, blocks ;
  hyper       = samples_in.attr("hyper") ;
  int npars   = hyper["n.pars"] ;
  int nmc     = hyper["nmc"] ;
  int nThins  = hyper["thin"] ;
  double rp   = hyper["rp"] ;    // rp is defined in initialise
  int start_R = hyper["start"] ; // start_R == 1;
  int start_C = start_R - 1 ;    // start_C == 0;
  int nsamp   = 1 + (nmc - start_R) * nThins ;
  double pMigrate, pMigrate_hyper, gammaMult ;
  //double gammaMult_hyper ;
  int report ;
  pMigrate       = setting["p.migrate"] ;
  pMigrate_hyper = setting["h.p.migrate"] ;
  gammaMult      = setting["gamma.mult"] ;
  //gammaMult_hyper= setting["h.gamma.mult"] ;
  report         = setting["report"] ;

  /* Kihon: Get hyper-level data.
  * First iteration data are from initialised 1st nmc, collected from
  the theta per chain per participant. The C++'s data_hyper mirrors
  (nSubjects x npars x nChains) R's "ps" array (nChains x nSubjects x npars),
  but I use Armadillo chain/slice major. */
  Rcpp::List phi             = hyper["phi"] ;
  arma::mat hlogPrior        = hyper["h_summed_log_prior"] ; // nmc x nChains
  arma::mat hlogLike         = hyper["h_log_likelihoods"] ;  // nmc x nChains
  Rcpp::List ppList          = hyper["pp.prior"] ;
  arma::cube data_hyper      = get_ps(samples_in) ; // 1st nmc theta
  arma::cube phiLocation     = phi[0] ; // nChains x npars x nmc
  arma::cube phiScale        = phi[1] ;
  arma::vec usehLogPrior     = vectorise(hlogPrior.row(start_C)) ;
  arma::vec usehLogLike      = vectorise(hlogLike.row(start_C)) ;
  arma::mat usePhiLocation   = phiLocation.slice(start_C) ;
  arma::mat usePhiScale      = phiScale.slice(start_C) ; // nChains x npars
  Rcpp::List usePhi          = Rcpp::List::create(
    Rcpp::Named("location")  = usePhiLocation,
    Rcpp::Named("scale")     = usePhiScale) ;
  Rcpp::List ppList_location = ppList[0] ;
  Rcpp::List ppList_scale    = ppList[1] ;
  int store_i                = start_C ; // store_i == 0;
  Rcpp::List ppList_constant = isConstant(ppList) ; // two logical vectors

  /* Setup for data level, get "p.priors" and start taking samples */
  Rcpp::List samples_s1        = samples_in[0] ;
  Rcpp::List pList             = samples_s1["p.prior"] ;
  Rcpp::List pLists            = assign_pp_pLists(samples_s1, usePhi) ;
  Rcpp::CharacterVector pNames = pList.names() ;
  Rcpp::CharacterVector pars   = pNames ;
  add_nmc0(samples_in, pLists) ;   // Add "use" for data level
  arma::mat update_theta, evolution, evolution_data ;
  arma::vec update_logprior, update_loglike ;
  Rcpp::List updated_use ;
  sampler DE, DE_data ;

  for(int i=1; i<nsamp; i++) {
    //Rcpp::Rcout << i << "  " ;
    /* Update usePhi, usehLogPrior, and usehLogLike */
    if ( Rf_runif(0, 1) < pMigrate_hyper ) {
      evolution = DE.migrate_hyper(usePhi, usehLogPrior, usehLogLike,
        pList, data_hyper, ppList, rp) ;
    } else {
      evolution = DE.crossover_hyper(usePhi, usehLogPrior, usehLogLike,
        pList, data_hyper, ppList, rp, gammaMult) ;
    }
    // nChains x npars
    // arma::mat temp1           = shuffle(evolution, 0);
    arma::mat temp1           = evolution ;
    usehLogPrior              = temp1.col(0) ; // nChains x 1
    usehLogLike               = temp1.col(1) ; //
    arma::mat usePhiLocation  = temp1.cols(2, 1+npars) ;
    arma::mat usePhiScale     = temp1.cols(2+npars, 1+2*npars) ;
    usePhi                    = Rcpp::List::create(  // nChains x npars
      Rcpp::Named("location") = usePhiLocation,
      Rcpp::Named("scale")    = usePhiScale) ;
    pLists = assign_pp_pLists(samples_s1, usePhi) ;

    /*  Update data level and hyper data */
    for(Rcpp::List::iterator it=samples_in.begin(); it!=samples_in.end(); ++it)
    {
      Rcpp::List subject_j     = *it ;
      Rcpp::List data_j        = subject_j["data"] ;
      Rcpp::List use           = subject_j.attr("use") ;
      int store_j              = use["store_i"] ;
      arma::mat useTheta_j     = use["theta"] ; // nChains x npars
      arma::vec useLogPrior_j  = update_useLogPrior(useTheta_j, pLists) ;
      arma::vec useLogLike_j   = use["loglike"] ; // nChains
      if ( Rf_runif(0, 1) < pMigrate )
      {
        evolution_data = DE_data.migrate_data_parallel(useTheta_j, useLogPrior_j,
          useLogLike_j, pLists, data_j, rp) ;
      } else {
        evolution_data = DE_data.crossover_data_parallel(useTheta_j, useLogPrior_j,
          useLogLike_j, pLists, data_j, rp, gammaMult) ;
      }
      // arma::mat temp2 = shuffle(evolution_data, 0);
      arma::mat temp2 = evolution_data ;
      update_logprior = temp2.col(0) ;
      update_loglike  = temp2.col(1) ;
      update_theta    = temp2.cols(2, 1+npars) ;

      // Harvest results; updated_use points at the same location as
      // attr(samples_i, "use")
      updated_use    = Rcpp::List::create(
        Rcpp::Named("theta")    = update_theta,    // nChains x npars
        Rcpp::Named("logprior") = update_logprior, // nChains x 1
        Rcpp::Named("loglike")  = update_loglike,  // nChains x 1
        Rcpp::Named("store_i")  = store_j,
        Rcpp::Named("p.priors") = pLists) ;
      subject_j.attr("use")  = updated_use ;

      if ( i % nThins == 0 ) {
        store_j++ ;  // store_ii starts from 0, nsamp starts from 1
        //Rcpp::Rcout << "\n" << store_ii << " . ";
        //if ((store_ii+1) % report == 0) { Rcpp::Rcout << store_ii+1 << " . "; }
        arma::mat tmp1       = subject_j["summed_log_prior"] ; // nmc x nChains
        arma::mat tmp2       = subject_j["log_likelihoods"] ;  // nmc x nChains
        arma::cube tmp3      = subject_j["theta"] ; // nChains x npars x nmc
        tmp1.row(store_j)    = update_logprior.t() ;
        tmp2.row(store_j)    = update_loglike.t() ;
        tmp3.slice(store_j)  = update_theta ;
        subject_j["summed_log_prior"] = tmp1 ;
        subject_j["log_likelihoods"]  = tmp2 ;
        subject_j["theta"]            = tmp3 ;
        // No need to feed back to subject_i !!!For DMC sake, still uses store_i
        updated_use["store_i"] = store_j ;
      }
    }

    data_hyper = get_ps_use(samples_in) ; // get new hyper data

    if ( i % nThins == 0 )
    {
      store_i++ ; // store_i starts from 0, nsamp starts from 1
      if ((store_i+1) % report == 0) { Rcpp::Rcout << store_i+1 << " "; }
      hlogPrior.row(store_i)      = usehLogPrior.t() ;  // nmc x nChains
      hlogLike.row(store_i)       = usehLogLike.t() ;   // nmc x nChains
      phiLocation.slice(store_i)  = temp1.cols(2, 1+npars) ; // nChains x npars x nmc
      phiScale.slice(store_i)     = temp1.cols(2+npars, 1+2*npars) ;
      Rcpp::List updated_phi      = Rcpp::List::create(
        Rcpp::Named("location")   = phiLocation,
        Rcpp::Named("scale")      = phiScale) ;
      hyper["h_summed_log_prior"] = hlogPrior ;
      hyper["h_log_likelihoods"]  = hlogLike ;
      hyper["phi"]                = updated_phi ;
    }
  }

  /* Yame */
  Rcpp::Rcout << std::endl ;
  samples_in.attr("hyper") = hyper ;
  return samples_in ;
}
