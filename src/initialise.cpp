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

arma::cube get_cps(Rcpp::List& samples) {
  int nSubjects = samples.size() ;
  Rcpp::List samples_s1 = samples[0] ;
  Rcpp::List data_s1          = samples_s1["data"] ;
  arma::cube modelCube        = data_s1.attr("model") ;
  Rcpp::NumericVector model   = data_s1.attr("model") ;
  Rcpp::NumericVector pVecNA  = model.attr("p.vector") ;
  Rcpp::CharacterVector pNames= pVecNA.names() ;
  int npars       = pNames.size() ;
  int chainFamily = 3 ;
  int nChains     = chainFamily*npars ;
  arma::cube cps(nChains, npars, nSubjects) ;
  for(Rcpp::List::iterator it=samples.begin(); it!=samples.end(); ++it)
  {
    int it_i = std::distance(samples.begin(), it) ;
    Rcpp::List subject_i = *it ;
    arma::cube theta_i = subject_i["theta"] ; // nChains x npars x nmc
    cps.slice(it_i) = theta_i.slice(0) ; // Take 1st nmc so index is 0
  }
  return cps ;
}

Rcpp::LogicalVector testHyper(Rcpp::List x, Rcpp::List y) {
  std::vector<std::string> a = x.names() ;
  std::vector<std::string> b = y.names() ;
  Rcpp::LogicalVector c(x.size()) ;
  for(std::vector<std::string>::iterator ita = a.begin();
    ita != a.end(); ++ ita)
  {
    int i = std::distance(a.begin(), ita) ;
    for(std::vector<std::string>::iterator itb = b.begin();
      itb != b.end(); ++ itb)
    {
      if(*ita == *itb) {c[i] = true;}
    }
  }
  return c ;
}

//' Set up a DMC Sample for a Participant
//'
//' This is C++ function to initialse a DMC sample. The user usually should
//' not directly use it.
//'
//' @param nmc number of Markov Chain Monte Carlo iteration
//' @param pList prior distribution setting. This is R's p.prior
//' @param data a model data instance created by \code{data.model.dmc}
//' @param samples a DMC posterior sample/object
//' @param theta1 A user supplied initial theta cube
//' @param startPList A user supplied (different) prior distribution setting.
//' This is R's start.prior.
//' @param setting a list container to store all DMC setting
//' @export
// [[Rcpp::export]]
Rcpp::List initialise_data(int nmc,
  Rcpp::Nullable<Rcpp::List> pList           = R_NilValue,
  Rcpp::Nullable<Rcpp::List> data            = R_NilValue,
  Rcpp::Nullable<Rcpp::List> samples         = R_NilValue,
  Rcpp::Nullable<Rcpp::NumericMatrix> theta1 = R_NilValue,
  Rcpp::Nullable<Rcpp::List> startPList      = R_NilValue,
  Rcpp::Nullable<Rcpp::List> setting         = R_NilValue)
{
  // Set up variables
  int nThins, restart, startFrom, nChains, npars ;
  int chainFamily = 3;
  bool add, verbose ;
  double rp ;
  Rcpp::List data_in, pList_in, startPList_in, samples_out, oPList, oStartPList,
  samples_in;
  Rcpp::NumericMatrix theta1_in ;
  arma::cube theta ;
  arma::mat summed_log_prior_mat, log_likelihoods ;
  Rcpp::NumericVector pVecNA, model ;
  Rcpp::CharacterVector pNames ;
  Rcpp::Nullable<Rcpp::NumericVector> remove ;

  // Test all Nullable arguments
  if (samples.isNull()) {
    data_in = Rcpp::as< Rcpp::List >(data) ;
    arma::cube modelCube            = data_in.attr("model") ;
    model       = data_in.attr("model") ;
    pVecNA      = model.attr("p.vector") ;
    pNames = pVecNA.names() ;
    npars = pNames.size() ;
    nChains = chainFamily * npars ;
  } else if (data.isNull()) {
    samples_in = Rcpp::as< Rcpp::List >(samples) ;
    data_in = samples_in["data"] ;   // Get pNames and npars
    arma::cube modelCube = data_in.attr("model") ;
    model       = data_in.attr("model") ;
    pVecNA      = model.attr("p.vector") ;
    pNames = pVecNA.names() ;
    npars = pNames.size() ;
    nChains = chainFamily * npars ;
  } else {
    Rcpp::stop("Neither samples nor data is found\n") ;
  }

  if (pList.isNull()) { Rcpp::stop("Must specify a p.prior argument\n") ; }
  else {
    pList_in = Rcpp::as< Rcpp::List >(pList) ;
    oPList = swap(pVecNA, pList_in) ; // ordered pList
  }

  if (setting.isNull()) {
    rp = 0.001 ;
    nThins = 1 ;
    startFrom = 1 ;
    add = false ;
    restart = true ;
    verbose = false;
    if (verbose)
    {
      Rcpp::Rcout << "uses the setting for: thin=1, rp=0.001, \n" ;
      Rcpp::Rcout << "start.from=1, " << " n.chains=" << nChains ;
      Rcpp::Rcout << " restart and verbose are true" << std::endl ;
    }
  } else {
    Rcpp::List setting_in = Rcpp::as< Rcpp::List >(setting) ;
    Rcpp::CharacterVector setting_names = setting_in.names();
    Rcpp::CharacterVector expected_names(7) ;
    expected_names[0] = "add" ;
    expected_names[1] = "rp"  ;
    expected_names[2] = "thin" ;
    expected_names[3] = "start.from" ;
    expected_names[4] = "restart" ;
    expected_names[5] = "verbose" ;
    expected_names[6] = "remove" ;
    // SEXP Rf_match(SEXP itable, SEXP ix, int nmatch)
    Rcpp::LogicalVector isSet = Rf_match(setting_names, expected_names, 0) ;
    add       = isSet[0] ? setting_in["add"] : false ;
    rp        = isSet[1] ? setting_in["rp"] : 0.001 ;
    nThins    = isSet[2] ? setting_in["thin"] : 1 ;
    startFrom = isSet[3] ? setting_in["start.from"] : 1 ;
    restart   = isSet[4] ? setting_in["restart"] : true ;
    verbose   = isSet[5] ? setting_in["verbose"] : false ;
    remove    = isSet[6] ? setting_in["remove"] : R_NilValue ;
    if(startFrom < 0) {Rcpp::stop("start.from must be an positive integer.\n") ; }
    if (verbose)
    {
      Rcpp::Rcout << "uses the setting: thin=" << nThins ;
      Rcpp::Rcout << ", rp=" << rp <<  ", start.from=" << startFrom ;
      Rcpp::Rcout << ", n.chains=" << nChains << ", restart=" << restart ;
      Rcpp::Rcout << " verbose=" << verbose << std::endl ;
    }
  }

  if (theta1.isNotNull())
  {
    Rcpp::Rcout << "User entered theta1 will be used. " << "\n" ;
    theta1_in = Rcpp::as< Rcpp::NumericMatrix >(theta1) ;
  }

  Rcpp::NumericVector initPVec ;

  /* Set up samples ---------------------------------------------------------*/
  if (samples.isNotNull())
  {
    if(verbose) {Rcpp::Rcout << "restarts at the samples location " ;}
    samples_in = Rcpp::as< Rcpp::List >(samples) ;
    Rcpp::List samples_in_cp(clone(samples_in)) ;

    arma::cube theta               = samples_in_cp["theta"] ;
    arma::mat summed_log_prior_mat = samples_in_cp["summed_log_prior"] ;
    arma::mat log_likelihoods      = samples_in_cp["log_likelihoods"] ;
    int nmc_in                     = samples_in_cp["nmc"] ;

    // 1. If the user requests to remove some samples
    if(remove.isNotNull())
    {
      Rcpp::NumericVector remove_in = Rcpp::as< Rcpp::NumericVector >(remove) ;
      int nRemove = remove_in.size() ;
      int remove_switch = 1;
      if(nRemove == 1) { remove_switch = Rcpp::as<int>(remove) ; }
      if (nRemove > 0 && remove_switch != 0)  // Remove switch
      {
        if (restart)
        {
          Rcpp::stop("restart must be FALSE, if remove is TRUE.\n") ;
        }
        if (std::isnan(remove_in[0]))
        {
          Rcpp::stop("Use FALSE to turn off remove, instead of NA.\n") ;
        }
        Rcpp::Rcout << "remove elements from " << remove_in[0] << " to "
                    << remove_in[nRemove-1] << std::endl ;
        for (int i=0; i<nRemove; i++)
        {
          int removeIdx = remove_in[i] - 1 ;
          theta.slice(removeIdx).fill(NA_REAL) ;
          summed_log_prior_mat.row(removeIdx).fill(NA_REAL) ;
          log_likelihoods.row(removeIdx).fill(NA_REAL) ;
        }
        samples_in_cp["theta"]            = theta ;
        samples_in_cp["summed_log_prior"] = summed_log_prior_mat ;
        samples_in_cp["log_likelihoods"]  = log_likelihoods ;
        samples_in_cp["thin"]             = nThins ;
      }
    }

    // 2. If the user wants to add new samples, we allocate spaces for it
    if (add)           // add on switch
    {
      if (verbose) {
        Rcpp::Rcout << "Add " << nmc << " new samples onto an existing one.\n" ;
      }
      arma::cube theta_add = resize(theta, nChains, npars, nmc+nmc_in) ;
      arma::mat summed_log_prior_add = resize(summed_log_prior_mat, nmc+nmc_in, nChains) ;
      arma::mat log_likelihoods_add  = resize(log_likelihoods,  nmc+nmc_in, nChains) ;
      theta_add.slices(nmc_in, nmc+nmc_in-1).fill(NA_REAL) ;
      summed_log_prior_add.rows(nmc_in, nmc+nmc_in-1).fill(-INFINITY) ;
      log_likelihoods_add.rows(nmc_in, nmc+nmc_in-1).fill(-INFINITY) ;

      samples_in_cp["theta"]            = theta_add ;
      samples_in_cp["summed_log_prior"] = summed_log_prior_add ;
      samples_in_cp["log_likelihoods"]  = log_likelihoods_add ;
      samples_in_cp["start"]            = nmc_in ; // should not add 1
      samples_in_cp["nmc"]              = nmc+nmc_in ;
      samples_in_cp["thin"]             = nThins ;
    } else if (restart && startFrom > 1) {
      if (startFrom > nmc_in)
      {
        Rcpp::Rcout << "start.from mustn't be greater than " << nmc_in << "\n";
        Rcpp::stop("Check your start.from.\n") ;
      }

      if(verbose) { Rcpp::Rcout << startFrom << std::endl; }

      arma::mat theta_new_start = theta.slice(startFrom-1) ;
      arma::mat summed_log_prior_new_start = summed_log_prior_mat.row(startFrom-1) ;
      arma::mat log_likelihoods_new_start  = log_likelihoods.row(startFrom-1) ;

      arma::cube theta_restart = resize(theta, nChains, npars, nmc) ;
      arma::mat summed_log_prior_restart = resize(summed_log_prior_mat, nmc, nChains) ;
      arma::mat log_likelihoods_restart  = resize(log_likelihoods,  nmc, nChains) ;

      theta_restart.fill(NA_REAL) ;
      summed_log_prior_restart.fill(-INFINITY) ;
      log_likelihoods_restart.fill(-INFINITY) ;


      theta_restart.slice(0) = theta_new_start ;
      summed_log_prior_restart.row(0) = summed_log_prior_new_start ;
      log_likelihoods_restart.row(0) = log_likelihoods_new_start ;

      samples_in_cp["theta"] = theta_restart ;
      samples_in_cp["summed_log_prior"] = summed_log_prior_restart ;
      samples_in_cp["log_likelihoods"]  = log_likelihoods_restart ;

      samples_in_cp["start"] = 1 ;
      samples_in_cp["thin"] = nThins ;
      samples_in_cp["nmc"]  = nmc ;
    }
    else {
      Rcpp::Rcout << "Neither add nor restart is true, but found samples!\n" ;
      Rcpp::Rcout << "Assuming you want to restart from the previous " <<
        "iteration\n";
      startFrom = nmc_in - 1 ;
      Rcpp::Rcout << "Change start.from to "<< startFrom+1 << " ";
      arma::mat theta_new_start = theta.slice(startFrom-1) ;
      arma::mat summed_log_prior_new_start = summed_log_prior_mat.row(startFrom-1) ;
      arma::mat log_likelihoods_new_start  = log_likelihoods.row(startFrom-1) ;

      theta.fill(NA_REAL) ;
      summed_log_prior_mat.fill(-INFINITY) ;
      log_likelihoods.fill(-INFINITY) ;

      arma::cube theta_restart = resize(theta, nChains, npars, nmc) ;
      arma::mat summed_log_prior_restart = resize(summed_log_prior_mat, nmc, nChains) ;
      arma::mat log_likelihoods_restart  = resize(log_likelihoods,  nmc, nChains) ;

      theta_restart.slice(0) = theta_new_start ;
      summed_log_prior_restart.row(0) = summed_log_prior_new_start ;
      log_likelihoods_restart.row(0) = log_likelihoods_new_start ;

      samples_in_cp["theta"] = theta_restart ;
      samples_in_cp["summed_log_prior"] = summed_log_prior_restart ;
      samples_in_cp["log_likelihoods"]  = log_likelihoods_restart ;

      samples_in_cp["start"] = 1 ;
      samples_in_cp["thin"] = nThins ;
      samples_in_cp["nmc"]  = nmc ;

    }

    samples_in_cp["p.prior"]  = oPList ;
    samples_out = samples_in_cp ;
  } else
  {
    if(verbose){ Rcpp::Rcout << "sets up a new samples\n" ; }
    summed_log_prior_mat = arma::mat(nmc, nChains).fill(-INFINITY) ;
    log_likelihoods      = arma::mat(nmc, nChains).fill(-INFINITY) ;
    theta                = arma::cube(nChains, npars, nmc).fill(NA_REAL) ;

    if (startPList.isNotNull())
    {
      startPList_in = Rcpp::as< Rcpp::List >(startPList) ;
      oStartPList   = swap(pVecNA, startPList_in) ;
    }

    if (theta1.isNotNull())
    {
      if( theta1_in.ncol() != npars || theta1_in.nrow() != nChains )
      {
        Rcpp::stop("The dimension of theta1 does not match theta dimension") ;
      }
    }


    // Start fresh from a new (NULL) samples
    if (verbose) Rcpp::Rcout << "Generating start points for each chain: " ;
    Rcpp::List oo_pList ;

    for (int i=0; i<nChains; i++)
    {
      if (verbose) Rcpp::Rcout << "." ;
      int j = 1;

      while( !summed_log_prior_mat.row(0).col(i).is_finite() ||
             !log_likelihoods.row(0).col(i).is_finite() )
      {
        if (!theta1.isNotNull())       // Generate parameters from prior
        {
          if (!startPList.isNotNull()) // sample from prior
          {
            theta.slice(0).row(i) = Rcpp::as<arma::rowvec>(rprior_scalar(oPList)) ;
            oo_pList = oPList ;
          } else                       // sample start.prior
          {
            theta.slice(0).row(i) = Rcpp::as<arma::rowvec>(rprior_scalar(oStartPList)) ;
            oo_pList = oStartPList ;
          }
        } else                         // Copy parameters from theta 1
        { // theta1 is a nChains x npar
          arma::mat theta1_arma(theta1_in.begin(), nChains, npars, false) ;
          theta.slice(0).row(i) = theta1_arma.row(i) ;
          oo_pList = oPList ;
        }

        arma::vec tmp = vectorise(theta.slice(0).row(i)) ;

//         initPVec = Rcpp::as<Rcpp::NumericVector>(Rcpp::wrap(tmp)) ;
//         initPVec.names() = pNames ;
//         Rcpp::CharacterVector tmp0 = initPVec.names();
//
//         double sumLogLike = 0 ;
//         std::vector<double> density = ddmc(data_in, initPVec) ;
//         for(std::vector<double>::iterator it1=density.begin();
//             it1!=density.end(); ++it1)
//         { sumLogLike += log(*it1) ; }
//
//         log_likelihoods.row(0).col(i) = sumLogLike ;
//         summed_log_prior_mat.row(0).col(i) = sum(dprior(initPVec, oPList)) ;
        log_likelihoods.row(0).col(i) = summed_log_likelihood(tmp, data_in) ;
        summed_log_prior_mat.row(0).col(i)= summed_log_prior(tmp, oo_pList) ;

        j++ ;
        if (j > 1e4)
        {
          Rcpp::Rcout << "Oops, fail to sample valid data-level parameter "
                      << "after 10,000 attempts. Please check the p.prior "
                      << "you entered" ;
          Rcpp::stop("Initiation fails.") ;
        }
      }

      samples_out = Rcpp::List::create(
        Rcpp::Named("theta")            = theta,
        Rcpp::Named("summed_log_prior") = summed_log_prior_mat,
        Rcpp::Named("log_likelihoods")  = log_likelihoods,
        Rcpp::Named("data")             = data,
        Rcpp::Named("p.prior")          = oPList,
        Rcpp::Named("start")            = 1,
        Rcpp::Named("n.pars")           = npars,
        Rcpp::Named("p.names")          = pNames,
        Rcpp::Named("rp")               = rp,
        Rcpp::Named("nmc")              = nmc,
        Rcpp::Named("thin")             = nThins,
        Rcpp::Named("n.chains")         = nChains) ;
    }
  }
  return samples_out ;
}

//' Set up a DMC Sample for Multiple Participants
//'
//' This is C++ function to initialse a DMC sample with multiple participants.
//' The user usually should not directly use it.
//'
//' @param nmc number of Markov Chain Monte Carlo iteration
//' @param pList prior distribution setting. This is R's p.prior
//' @param data a model data instance created by \code{data.model.dmc}
//' @param samples a DMC posterior sample
//' @param theta1 A user supplied initial theta cube
//' @param startPList A user supplied (different) prior distribution setting.
//' This is R's start.prior.
//' @param phi1 A user supplied initial phi cube
//' @param ppList prior distribution setting, hyper level
//' @param hStartPList A user supplied (different) hyper-prior distribution
//' setting
//' @param setting a list container to store all DMC setting
//' @export
// [[Rcpp::export]]
Rcpp::List initialise_hyper(int nmc,
  Rcpp::Nullable<Rcpp::List> pList           = R_NilValue,
  Rcpp::Nullable<Rcpp::List> data            = R_NilValue,
  Rcpp::Nullable<Rcpp::List> samples         = R_NilValue,
  Rcpp::Nullable<Rcpp::NumericMatrix> theta1 = R_NilValue,
  Rcpp::Nullable<Rcpp::List> startPList      = R_NilValue,
  Rcpp::Nullable<Rcpp::List> phi1            = R_NilValue,
  Rcpp::Nullable<Rcpp::List> ppList          = R_NilValue,
  Rcpp::Nullable<Rcpp::List> hStartPList     = R_NilValue,
  Rcpp::Nullable<Rcpp::List> setting         = R_NilValue)
{
  /* Yoi: Set up variables  -------------------------------------------------*/
  double rp ;
  bool add, verbose, restart ;
  int nThins, startFrom, nChains, npars, nSubjects ;
  int chainFamily = 3;
  Rcpp::List data_in, data_s1, pList_in, startPList_in, samples_in, newhyper,
  samples_cp, oPList, oStartPList, phi1_in, ppList_in,
  setting_in, phi, hyper;

  Rcpp::NumericMatrix theta1_in ;
  Rcpp::NumericVector pVecNA, remove, model ;
  Rcpp::LogicalVector hasHyper, hasSigma ;
  Rcpp::CharacterVector sNames, pNames ; // subject and parameter names
  arma::cube theta, phiLoc, phiSca ;
  arma::mat summed_log_prior_mat, log_likelihoods, h_summed_log_prior,
  h_log_likelihoods, consts ;

  /* Kokutsu-dachi: Test 4 R_NilValues -----------------------------------*/
  if (setting.isNull()) {
    add        = false ;
    rp         = 0.001 ;
    nThins     = 1 ;
    startFrom  = 1 ;
    restart    = false ;
    verbose    = false ;
    setting_in = Rcpp::List::create(
      Rcpp::Named("add")        = add,
      Rcpp::Named("rp")         = rp,
      Rcpp::Named("thin")       = nThins,
      Rcpp::Named("start.from") = startFrom,
      Rcpp::Named("verbose")    = verbose,
      Rcpp::Named("restart")    = restart,
      Rcpp::Named("remove")     = R_NilValue ) ;
  } else {
    setting_in = Rcpp::as< Rcpp::List >(setting) ;
    Rcpp::CharacterVector setting_names = setting_in.names();
    Rcpp::CharacterVector expected_names(7) ;
    expected_names[0] = "add" ;
    expected_names[1] = "rp"  ;
    expected_names[2] = "thin" ;
    expected_names[3] = "start.from" ;
    expected_names[4] = "restart" ;
    expected_names[5] = "verbose" ;
    expected_names[6] = "remove" ;
    // SEXP Rf_match(SEXP itable, SEXP ix, int nmatch)
    Rcpp::LogicalVector isSet = Rf_match(setting_names, expected_names, 0) ;
    add       = isSet[0] ? setting_in["add"] : false ;
    rp        = isSet[1] ? setting_in["rp"] : 0.001 ;
    nThins    = isSet[2] ? setting_in["thin"] : 1 ;
    startFrom = isSet[3] ? setting_in["start.from"] : 1 ;
    restart   = isSet[4] ? setting_in["restart"] : false ;
    verbose   = isSet[5] ? setting_in["verbose"] : true ;
    remove    = isSet[6] ? setting_in["remove"] : R_NilValue ;
    if(startFrom < 0) {Rcpp::stop("start.from must be an positive integer.\n") ; }
    setting_in = Rcpp::List::create(
      Rcpp::Named("add")     = add,
      Rcpp::Named("rp")      = rp,
      Rcpp::Named("thin")    = nThins,
      Rcpp::Named("start.from") = startFrom,
      Rcpp::Named("restart") = restart,
      Rcpp::Named("verbose") = verbose,
      Rcpp::Named("remove")  = R_NilValue ) ;
  }

  if(data.isNull()) {
    if (samples.isNull()) { Rcpp::stop("Must specify 'samples'."); }
    if (verbose) Rcpp::Rcout << "Found samples. Using samples to initialise.\n" ;
    samples_in = Rcpp::as< Rcpp::List >(samples) ;
    samples_cp = clone(samples_in) ;
    nSubjects  = samples_cp.size() ;
    sNames     = samples_cp.names() ; // subject names
    Rcpp::List samples_s1 = samples_cp[0] ;
    hyper                 = samples_cp.attr("hyper") ;
    pList_in              = samples_s1["p.prior"] ;
    nThins                = samples_s1["thin"] ;
    data_s1               = samples_s1["data"] ;
    arma::cube modelCube  = data_s1.attr("model") ;
    model   = data_s1.attr("model") ;
    pVecNA  = model.attr("p.vector") ;
    pNames  = pVecNA.names() ;
    npars   = pNames.size() ;
    nChains = chainFamily * npars ;
  } else if (samples.isNull()) {
    if (data.isNull()) { Rcpp::stop("Must specify 'data'."); }
    if (verbose) { Rcpp::Rcout << "Found data. Using data to initialise \n"; }

    data_in    = Rcpp::as< Rcpp::List >(data) ;
    nSubjects  = data_in.size() ;
    sNames     = data_in.names() ; // subject names
    data_s1    = data_in[0] ;
    arma::cube modelCube      = data_s1.attr("model") ;
    model   = data_s1.attr("model") ;
    pVecNA  = model.attr("p.vector") ;
    pNames  = pVecNA.names() ;
    npars   = pNames.size() ;
    nChains = chainFamily * npars ;

    if (pList.isNull()) {
      Rcpp::stop("Must specify a parameter prior list (p.prior)\n") ;
    } else {
      pList_in = Rcpp::as< Rcpp::List >(pList) ;
    }
  } else {
    Rcpp::stop("Neither data nor samples are found!") ;
  }


  if (theta1.isNull()) {
    if(verbose) { Rcpp::Rcout << "theta1 not found;  " ; }
  } else {
    theta1_in = Rcpp::as< Rcpp::NumericMatrix >(theta1) ;
  }

  if (startPList.isNull()) {
    if(verbose) { Rcpp::Rcout << "start.prior not found;  " ; }
  } else {
    startPList_in = Rcpp::as< Rcpp::List >(startPList) ;
    //oStartPList   = swap(pVecNA, startPList_in) ;
  }


  Rcpp::List samples_out(nSubjects) ;
  Rcpp::List constantPrior_List ; // two npars LogicalVector

  /* Initialise non-hierarchical model -------------------------------*/
  if (ppList.isNull()) {
    if(verbose) {
      Rcpp::Rcout << "pp.prior not found. Initialise non-hierarchical model\n";
    }
    if (samples.isNull()) {// No samples is found
      Rcpp::Rcout << "Generating start points for each participant ('.'= 1): " ;
      for(Rcpp::List::iterator it = data_in.begin(); it != data_in.end(); ++it)
      {
        int subjIdx = std::distance(data_in.begin(), it) ;
        Rcpp::Rcout << "." ;
        // "initialise" arguments are
        // (1) nmc, (2) p.prior, (3) data, (4) samples, (5) theta1,
        // (6) start.prior, (7) setting
        // -- samples.isNull, so send a R_NilValue for (4).
        // -- add must be FALSE, so send R_NilValue for theta1 and start.prior;
        // -- setting_in inherits from initialise.R
        samples_out[subjIdx] = initialise_data(nmc, pList_in, *it, R_NilValue,
          R_NilValue, R_NilValue, setting_in) ;
        samples_out.names() = sNames ;
      }
    }
  } else {
    /* Initialise hierarchical model -------------------------------
    * Found hyper prior parameter list (ie pp.prior)
    * checks for pp.prior's are in initialise.R */
    Rcpp::List loc, sca, oLoc, oSca ; // raw and ordered location and scale pList
    // Check if h.start.prior is set; Assign loc, sca etc outside chain loop

    if (hStartPList.isNotNull()) {
      Rcpp::Rcout << "Use hstart.prior (ie use's pp.prior) to initialise \n" ;
      Rcpp::Rcout << "hierarchical model.\n";
      ppList_in = Rcpp::as< Rcpp::List >(hStartPList) ;
      loc  = ppList_in[0] ;
      sca  = ppList_in[1] ;
      oLoc = swap(pVecNA, loc) ;
      oSca = swap(pVecNA, sca) ; // check dist attribute to see if it is constant
      constantPrior_List = isConstant(ppList_in) ;
      hasSigma = testHyper(loc, sca) ;
      hasHyper = testHyper(pList_in, sca) ;
    } else {
      if (verbose) Rcpp::Rcout << "Use pp.prior to initialise hierarchical model\n";
      ppList_in = Rcpp::as< Rcpp::List >(ppList) ; // a necessary redundancy
      loc  = ppList_in[0] ;
      sca  = ppList_in[1] ;
      oLoc = swap(pVecNA, loc) ;
      oSca = swap(pVecNA, sca) ;
      constantPrior_List = isConstant(ppList_in) ;
      hasSigma = testHyper(loc, sca) ;
      hasHyper = testHyper(pList_in, sca) ;
    }

    /* If samples is not provided, add=true is impossible-------------------*/
    if (samples.isNull()) {
      if (add) {Rcpp::stop("No samples, so add cannot be true!") ; }
      h_summed_log_prior = arma::mat(nmc, nChains).fill(-INFINITY) ;
      h_log_likelihoods  = arma::mat(nmc, nChains).fill(-INFINITY) ;
      phiLoc             = arma::cube(nChains, npars, nmc).fill(NA_REAL) ;
      phiSca             = arma::cube(nChains, npars, nmc).fill(NA_REAL) ;

      /*1. --If phi1 is NULL, use rprior to generate phi-Loc/Sca.slice(0)
      --Then, use dprior to get log prior prob. density
      --Next, sum over data-level parameters (eg, a, v, z, etc ... ) and
      over hyper-level parameters (ie location and scale)
      -- This step produces h_summed_log_prior ---------------------*/
      if (phi1.isNull()) {
        if (verbose) Rcpp::Rcout << "Use rprior to generating hyper-start points for each chain: " ;
        for (int i=0; i<nChains; i++) {
          if (verbose) Rcpp::Rcout << i+1 << " " ; // one chain per dot
          // Use rprior to generate random initial hyper-parameters
          // initialise.R has checked location and scale prior
          // If scale prior is not set, its dist attribute will be "constant",
          // hence rprior_scalar will return 1
          phiLoc.slice(0).row(i) = Rcpp::as<arma::rowvec>(rprior_scalar(oLoc)) ;
          phiSca.slice(0).row(i) = Rcpp::as<arma::rowvec>(rprior_scalar(oSca)) ;

          // Use dprior to get probability densities
          Rcpp::NumericVector initPVec_loc ;
          Rcpp::NumericVector initPVec_sca ;
          arma::mat tmp_loc = phiLoc.slice(0).row(i) ;
          arma::mat tmp_sca = phiSca.slice(0).row(i) ;
          initPVec_loc = Rcpp::as<Rcpp::NumericVector>(Rcpp::wrap(tmp_loc)) ;
          initPVec_sca = Rcpp::as<Rcpp::NumericVector>(Rcpp::wrap(tmp_sca)) ;
          initPVec_loc.names() = pNames ;
          initPVec_sca.names() = pNames ;
          // constant's den value will be 0, because log == TRUE (default)
          h_summed_log_prior.row(0).col(i) =
            sum(dprior(initPVec_loc, oLoc)) + sum(dprior(initPVec_sca, oSca)) ;
        } // Finish nChains loop
      } else {// else branch for the phi1 test; if phi1 is found
        Rcpp::Rcout << "Found phi1, using it to generating hyper-start points\n";
        Rcpp::Rcout << "for each chain: " ;
        phi1_in = Rcpp::as< Rcpp::List >(phi1) ;
        // Assume phi1_in[0] and phi1_in[1] are nChains x npars matrices
        arma::mat phi1_loc = phi1_in[0] ;
        arma::mat phi1_sca = phi1_in[1] ;
        for(int i=0; i<nChains; i++) {// overuse i same thing as in phi1.isNull
          Rcpp::Rcout << "." ; // one chain per dot
          Rcpp::NumericVector initPVec_loc ;
          Rcpp::NumericVector initPVec_sca ;
          arma::mat tmp_loc = phi1_loc.row(i) ;
          arma::mat tmp_sca = phi1_sca.row(i) ;
          initPVec_loc = Rcpp::as<Rcpp::NumericVector>(Rcpp::wrap(tmp_loc)) ;
          initPVec_sca = Rcpp::as<Rcpp::NumericVector>(Rcpp::wrap(tmp_sca)) ;
          initPVec_loc.names() = pNames ;
          initPVec_sca.names() = pNames ;
          h_summed_log_prior.row(0).col(i) =
            sum(dprior(initPVec_loc, oLoc)) + sum(dprior(initPVec_sca, oSca)) ;
        }
      }

      /* 2. If theta1 is NULL, ...
      * -- This step produces samples_out -------------------------------*/
      if (theta1.isNull()) {
        if (verbose) {
          Rcpp::Rcout << "\nGenerating start points for each participant ('.' = "
        << "1 participant): " ;
        }
        for(Rcpp::List::iterator it = data_in.begin(); it != data_in.end(); ++it)
        {
          Rcpp::Rcout << "." ; // one subj per dot
          int it_i = std::distance(data_in.begin(), it) ;
          // "initialise" arguments are
          // (1) nmc, (2) p.prior, (3) data, (4) samples, (5) theta1,
          // (6) start.prior, (7) setting
          // -- still under samples.isNull, so send a R_NilValue for (4).
          // -- add only work when samples is supplied. theta1.isNUll, so send
          //    R_NilValue for (5) and assume start.prior == NULL;
          // -- setting_in inherits from initialise.R
          samples_out[it_i] = initialise_data(nmc, pList_in, *it, R_NilValue,
            R_NilValue, R_NilValue, setting_in) ;
          samples_out.names() = sNames ;
        }
      } else {  // If the user does provide theta1
        theta1_in = Rcpp::as< Rcpp::NumericMatrix >(theta1) ;
        if (verbose) Rcpp::Rcout << "\nGenerate from theta1 ('.' = 1 participant): " ;
        for(Rcpp::List::iterator it = data_in.begin(); it != data_in.end(); ++it)
        {
          int it_i = std::distance(data_in.begin(), it) ;
          Rcpp::Rcout << "." ; // one subj per dot
          // "initialise" arguments are
          // (1) nmc, (2) p.prior, (3) data, (4) samples, (5) theta1,
          // (6) start.prior, (7) setting
          // -- still under samples.isNull, so send a R_NilValue for (4).
          // -- add only work when samples is supplied. theta1.isNotNUll, so
          //    (5)==theta1_in, and assume start.prior == NULL;
          // -- setting_in inherits from initialise.R
          samples_out[it_i] = initialise_data(nmc, pList_in, *it, R_NilValue,
            theta1_in, R_NilValue, setting_in) ;
          samples_out.names() = sNames ;
        }
      }

      /* 3. This step produces h_log_likelihoods */
      double sumOverSubject = 0; // cps refers to nChains x npars x nSubjects cube
      arma::cube cps = get_cps(samples_out) ;
      // phiLoc is a cube with nChains x npars x nmc
      for(int k=0; k<nChains; k++) {
        arma::vec phiLoc1nmc = phiLoc.slice(0).row(k).t() ;
        arma::vec phiSca1nmc = phiSca.slice(0).row(k).t() ;
        Rcpp::List phi_i = Rcpp::List::create(
          Rcpp::Named("location") = phiLoc1nmc ,
          Rcpp::Named("scale")    = phiSca1nmc ) ;
        assign_pp_pList(phi_i, pList_in) ; // pList_in is modified here

        for(int j=0; j<nSubjects; j++) {
          // cps is a nChains x npars x nSubjects cube from the 1st nmc of theta
          arma::mat ps_si          = cps.slice(j).row(k) ;
          Rcpp::NumericVector pVec = Rcpp::as<Rcpp::NumericVector>(Rcpp::wrap(ps_si)) ;
          pVec.names() = pNames ;
          // First sum over parameters and then sum over subjects
          sumOverSubject += sum(dprior(pVec, pList_in)) ;
        }
        h_log_likelihoods.row(0).col(k) = sumOverSubject ;
      }

      // Check if all elements/chains are finite; // nmc x nChains
      if (h_log_likelihoods.row(0).has_inf()) {
        Rcpp::Rcout << "Sampling start points was valid for data but not hyper"
                    << "level.\n" ;
        Rcpp::stop("Try a tighter pp.prior or hstart.prior\n") ;
      }

      // Correct data level prior given hypers
      // phiLoc is a cube with nChains x npars x nmc
      Rcpp::LogicalVector constantPrior = constantPrior_List[0] ; //location dist attr
      int constantSize = Rcpp::sum(constantPrior);
      Rcpp::List ppList_constant(constantSize) ; // pp.priors
      arma::mat consts(nChains, constantSize) ;
      Rcpp::CharacterVector consts_names(constantSize) ;
      int consts_i = 0;
      if (constantSize > 0) {
        for(Rcpp::LogicalVector::iterator it2 = constantPrior.begin();
          it2 != constantPrior.end() ; ++it2)
        {
          int it2_i = std::distance(constantPrior.begin(), it2) ;
          if(*it2)
          { // nCahins x constantSize
            consts.col(consts_i) = phiLoc.slice(0).col(it2_i) ;
            consts_names(consts_i) = pNames[it2_i] ;
            for(int k=0; k<ppList_constant.size(); k++)
            {
              ppList_constant[k] = oLoc[it2_i] ; // pp.priors
            }
            consts_i++;
          }
        }
      }

      Rcpp::List subject1 = samples_out[0] ;
      Rcpp::List pList_l  = subject1["p.prior"] ;
      for (int l =0; l<nChains; l++) {
        arma::vec phiLoc1nmc = phiLoc.slice(0).row(l).t() ;
        arma::vec phiSca1nmc = phiSca.slice(0).row(l).t() ;
        Rcpp::List phi_l = Rcpp::List::create(
          Rcpp::Named("location") = phiLoc1nmc ,
          Rcpp::Named("scale")    = phiSca1nmc );
        assign_pp_pList(phi_l, pList_l) ;

        // Loop through all subjects; this takes quite a while
        for(Rcpp::List::iterator it4 = samples_out.begin();
          it4 != samples_out.end() ; ++it4) {
          Rcpp::List subject_i = *it4 ;
          arma::cube subject_i_theta = subject_i["theta"] ;
          arma::mat pVec_arma1 = subject_i_theta.slice(0).row(l) ;
          arma::mat pVec_arma2 = consts.row(l) ;
          Rcpp::NumericVector pVec1 = Rcpp::as<Rcpp::NumericVector>(Rcpp::wrap(pVec_arma1)) ;
          Rcpp::NumericVector pVec2 = Rcpp::as<Rcpp::NumericVector>(Rcpp::wrap(pVec_arma2)) ;
          pVec1.names() = pNames ;
          pVec2.names() = consts_names ;
          double tmp ;
          if (constantSize==0) {
            tmp = sum(dprior(pVec1, pList_l));
          } else {
            Rcpp::Rcout << "Found constant hyperprior\n" ;
            tmp = sum(dprior(pVec1, pList_l)) + sum(dprior(pVec2, ppList_constant)) ;
          }
          arma::mat summed_log_prior_tmp = subject_i["summed_log_prior"] ;
          summed_log_prior_tmp.row(0).col(l) = tmp ;
          subject_i["summed_log_prior"] = summed_log_prior_tmp ;
        }
      }

      phi = Rcpp::List::create(
        Rcpp::Named("location") = phiLoc,
        Rcpp::Named("scale")    = phiSca) ;

      Rcpp::List ppList_out = Rcpp::List::create(
        Rcpp::Named("location") = oLoc,
        Rcpp::Named("scale")    = oSca) ;

      hyper = Rcpp::List::create(
        Rcpp::Named("phi")                = phi,
        Rcpp::Named("h_log_likelihoods")  = h_log_likelihoods,
        Rcpp::Named("h_summed_log_prior") = h_summed_log_prior,
        Rcpp::Named("pp.prior")           = ppList_out,
        Rcpp::Named("start")              = startFrom,
        Rcpp::Named("n.pars")             = npars,
        Rcpp::Named("p.names")            = pNames,
        Rcpp::Named("rp")                 = setting_in["rp"],
        Rcpp::Named("nmc")                = nmc,
        Rcpp::Named("n.chains")           = nChains,
        Rcpp::Named("thin")               = setting_in["thin"],
        Rcpp::Named("has.sigma")          = hasSigma,
        Rcpp::Named("has.hyper")          = hasHyper) ;

      samples_out.attr("hyper") = hyper ;
      /* If samples exists, restart from previous sampling------------------*/
    } else {
      if (verbose) Rcpp::Rcout << "Participant " ;
      for(Rcpp::List::iterator it = samples_cp.begin(); it != samples_cp.end(); ++it)
      {
        int it_i = std::distance(samples_cp.begin(), it) ;
        if (verbose) Rcpp::Rcout << it_i+1 << " " ;
        Rcpp::List samples_i = *it ;
        // "initialise" arguments are (1) nmc, (2) p.prior, (3) data,
        // (4) samples, (5) theta1, (6) start.prior, (7) setting
        // -- !!!samples exit!!!, so (3) == NULL and (4) == samples_i
        // -- add == TRUE is handled later and theta1.isNull, so send
        //    R_NilValue for (5) and assume start.prior == NULL;
        // -- setting_in inherits from initialise.R
        samples_out[it_i] = initialise_data(nmc, pList_in, R_NilValue,
          samples_i, R_NilValue, R_NilValue, setting_in) ;
      }
      samples_out.names() = sNames ;


      /* Handle add only work when samples is supplied ------------------*/
      if (add) {
        int previous_nmc = hyper["nmc"] ;
        int newnmc = previous_nmc + nmc ;
        phiLoc = arma::cube(nChains, npars, newnmc).fill(NA_REAL) ;
        phiSca = arma::cube(nChains, npars, newnmc).fill(NA_REAL) ;
        h_summed_log_prior = arma::mat(newnmc, nChains).fill(-INFINITY) ;
        h_log_likelihoods  = arma::mat(newnmc, nChains).fill(-INFINITY) ;

        Rcpp::List p_phi = hyper["phi"] ;
        arma::mat p_h_summed_log_prior = hyper["h_summed_log_prior"] ;
        arma::mat p_h_log_likelihoods  = hyper["h_log_likelihoods"] ;
        arma::cube p_phiLoc = p_phi[0] ;
        arma::cube p_phiSca = p_phi[1] ;
        for(int i=0; i<previous_nmc; i++) {
          phiLoc.slice(i) = p_phiLoc.slice(i) ;
          phiSca.slice(i) = p_phiSca.slice(i) ; // Correct p_phiLoc to p_phiSca
          h_summed_log_prior.row(i) = p_h_summed_log_prior.row(i) ;
          h_log_likelihoods.row(i)  = p_h_log_likelihoods.row(i) ;
        }

        Rcpp::List newphi = Rcpp::List::create(
          Rcpp::Named("location") = phiLoc,
          Rcpp::Named("scale")    = phiSca) ;

        newhyper = Rcpp::List::create(
          Rcpp::Named("phi")                = newphi,
          Rcpp::Named("h_log_likelihoods")  = h_log_likelihoods,
          Rcpp::Named("h_summed_log_prior") = h_summed_log_prior,
          Rcpp::Named("pp.prior")           = ppList_in,
          Rcpp::Named("start")              = previous_nmc, // no +1
          Rcpp::Named("n.pars")             = npars,
          Rcpp::Named("p.names")            = pNames,
          Rcpp::Named("rp")                 = setting_in["rp"],
          Rcpp::Named("nmc")                = newnmc,
          Rcpp::Named("n.chains")           = nChains,
          Rcpp::Named("thin")               = setting_in["thin"],
          Rcpp::Named("has.sigma")          = hasSigma,
          Rcpp::Named("has.hyper")          = hasHyper) ;
      } else { // add else branch
        phiLoc = arma::cube(nChains, npars, nmc).fill(NA_REAL) ;
        phiSca = arma::cube(nChains, npars, nmc).fill(NA_REAL) ;
        h_summed_log_prior = arma::mat(nmc, nChains).fill(-INFINITY) ;
        h_log_likelihoods  = arma::mat(nmc, nChains).fill(-INFINITY) ;

        Rcpp::List p_phi = hyper["phi"] ;
        arma::mat p_h_summed_log_prior = hyper["h_summed_log_prior"] ;
        arma::mat p_h_log_likelihoods  = hyper["h_log_likelihoods"] ;
        arma::cube previous_phiLoc = p_phi[0] ;
        arma::cube previous_phiSca = p_phi[1] ;

        int old_nmc = p_h_log_likelihoods.n_rows ;  // correct for R index
        phiLoc.slice(0) = previous_phiLoc.slice(old_nmc-1) ;
        phiSca.slice(0) = previous_phiSca.slice(old_nmc-1) ;
        h_summed_log_prior.row(0) = p_h_summed_log_prior.row(old_nmc-1) ;
        h_log_likelihoods.row(0)  = p_h_log_likelihoods.row(old_nmc-1) ;

        Rcpp::List newphi = Rcpp::List::create(
          Rcpp::Named("location") = phiLoc,
          Rcpp::Named("scale")    = phiSca) ;
        newhyper = Rcpp::List::create(
          Rcpp::Named("phi")                = newphi,
          Rcpp::Named("h_log_likelihoods")  = h_log_likelihoods,
          Rcpp::Named("h_summed_log_prior") = h_summed_log_prior,
          Rcpp::Named("pp.prior")           = ppList_in,
          Rcpp::Named("start")              = 1,
          Rcpp::Named("n.pars")             = npars,
          Rcpp::Named("p.names")            = pNames,
          Rcpp::Named("rp")                 = setting_in["rp"],
          Rcpp::Named("nmc")                = nmc,
          Rcpp::Named("n.chains")           = nChains,
          Rcpp::Named("thin")               = setting_in["thin"],
          Rcpp::Named("has.sigma")          = hasSigma,
          Rcpp::Named("has.hyper")          = hasHyper) ;
      }

      samples_out.attr("hyper") = newhyper ;
    }
  }
  return samples_out ;
}
