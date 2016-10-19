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

//' Calculate Prior Probability Density for an EAM
//'
//' \code{dprior} is C++ function. It matches 5 string types: \code{tnorm},
//' \code{beta_lu}, \code{gamma_l}, \code{lnorm_l}, and \code{constant} to
//' determine which density functions to call (via R API). For truncated normal
//' density, \code{dprior} calls \code{dtn_scalar}, an internal Rcpp function
//' built specific for \pkg{ggdmc}. Whetehr log the probability density is
//' determined by the boolean \code{log} sent in via \code{p.prior}. This
//' function is akin to DMC's \code{log.prior.dmc}.
//'
//' @param pVec the user's supplied parameter vector or a sampler supplied
//' theta/phi vector.
//' @param pPrior a list of list usually created by prior.p.dmc to store the
//' prior parameter setting.
//' @return a double vector with probability density for each model parameter
//' @export
//' @examples
//' ## Use Drift-diffusion model as an example
//' ddm.prior <- prior.p.dmc(
//' dists = c("tnorm", "tnorm", "beta", "tnorm", "beta", "beta"),
//' p1    = c(a = 1, v = 0, z = 1, sz = 1, sv = 1, t0 = 1),
//' p2    = c(a = 1, v = 2, z = 1, sz = 1, sv = 1, t0 = 1),
//' lower = c(0,-5, NA, NA, 0, NA),
//' upper = c(2, 5, NA, NA, 2, NA))
//'
//' view(ddm.prior)
//' ##      mean sd lower upper log    dist  untrans
//' ##   a     1  1     0     2   1   tnorm identity
//' ##   v     0  2    -5     5   1   tnorm identity
//' ##   z     1  1     0     1   1 beta_lu identity
//' ##   sz    1  1  -Inf   Inf   1   tnorm identity
//' ##   sv    1  1     0     2   1 beta_lu identity
//' ##   t0    1  1     0     1   1 beta_lu identity
//'
//' ddm.pVec <- c(a=1.15, v=-0.10, z=0.74, sz=1.23, sv=0.11, t0=0.87)
//' dprior(ddm.pVec, ddm.prior)
//' ##         a          v          z         sz         sv         t0
//' ##-0.5484734 -1.6008386  0.0000000 -0.9453885  0.0000000  0.0000000
//'
//' ## Use LBA model as an example
//' lba.prior <- prior.p.dmc(
//' dists = c("tnorm", "tnorm", "tnorm", "tnorm", "tnorm", "tnorm"),
//' p1    = c(A=.4, B=.6, mean_v.true=1,  mean_v.false=0,  sd_v.true=.5, t0=.3),
//' p2    = c(A=.1, B=.1, mean_v.true=.2, mean_v.false=.2, sd_v.true=.1, t0=.05),
//' lower = c(0,   0, NA, NA, 0, .1),
//' upper = c(NA, NA, NA, NA, NA, 1))
//'
//' view(lba.prior)
//' ##              mean   sd lower upper log  dist  untrans
//' ## A             0.4  0.1     0   Inf   1 tnorm identity
//' ## B             0.6  0.1     0   Inf   1 tnorm identity
//' ## mean_v.true     1  0.2  -Inf   Inf   1 tnorm identity
//' ## mean_v.false    0  0.2  -Inf   Inf   1 tnorm identity
//' ## sd_v.true     0.5  0.1     0   Inf   1 tnorm identity
//' ## t0            0.3 0.05   0.1     1   1 tnorm identity
//'
//' lba.pVec <- c(A=0.398, B=0.614, mean_v.true=1.040,
//'               mean_v.false=-0.032, sd_v.true=0.485, t0=0.271)
//' dprior(lba.pVec, lba.prior)
//' ##          A            B  mean_v.true mean_v.false    sd_v.true
//' ##  1.3834782    1.3738466    0.6704994    0.6776994    1.3723968
// [[Rcpp::export]]
Rcpp::NumericVector dprior(Rcpp::NumericVector& pVec, Rcpp::List& pPrior)
{
  Rcpp::List l1;   // a container to loop through inner list
  std::vector<std::string> pNames = pVec.names() ;
  Rcpp::NumericVector out = Rcpp::NumericVector(pVec.size());

  std::string distType1 ("tnorm");  // Available pdf' in DMC
  std::string distType2 ("beta_lu");
  std::string distType3 ("gamma_l");
  std::string distType4 ("lnorm_l");
  std::string distType5 ("constant");

  int islog;
  double x, p1, p2, lower, upper, den;

  for (Rcpp::NumericVector::iterator it = pVec.begin(); it != pVec.end(); ++it)
  {
    int i = std::distance(pVec.begin(), it) ;
    l1 = pPrior[pNames[i]];
    std::string distName = l1.attr("dist");

    x  = pVec[i];
    p1 = Rcpp::as<double>(l1[0]);  // mean; shape1; shape; meanlog
    p2 = Rcpp::as<double>(l1[1]);  // sd;   shape2; scale; sdlog

    // Do do.call
    if (distName.compare(distType1) == 0) {         // tnorm
      lower = Rcpp::as<double>(l1[2]);
      upper = Rcpp::as<double>(l1[3]);
      islog = Rcpp::as<int>(l1[4]);
      den = dtn_scalar(x, p1, p2, lower, upper, islog);
    } else if (distName.compare(distType2) == 0) {  // beta_ul
      lower = Rcpp::as<double>(l1[2]);
      upper = Rcpp::as<double>(l1[3]);
      islog = Rcpp::as<int>(l1[4]);
      den = R::dbeta((x-lower) / (upper-lower), p1, p2, islog);
    } else if (distName.compare(distType3) == 0) {  // gamma_l
      lower = Rcpp::as<double>(l1[2]);
      islog = Rcpp::as<int>(l1[3]);
      den = R::dgamma(x-lower, p1, p2, islog);
    } else if (distName.compare(distType4) == 0) {  // lnorm_l
      lower = Rcpp::as<double>(l1[2]);
      islog = Rcpp::as<int>(l1[3]);
      den = R::dlnorm(x-lower, p1, p2, islog);
    } else if (distName.compare(distType5) == 0) {  // constant
      islog = Rcpp::as<int>(l1[1]);
      den = islog ? 0 : 1;
    } else {
      Rcpp::Rcout << "Distribution type not yet supported" << "\n";
      den = 0;
    }
    out[i] = den;
  }

  out.attr("names") = pNames;
  return out;
}

Rcpp::NumericVector rprior_scalar(Rcpp::List priorList)
{
  Rcpp::List l1;   // a list to hold each parameter's setting inside priorList
  Rcpp::NumericVector out = Rcpp::NumericVector(priorList.size());
  std::vector<std::string> parameterNames = priorList.attr("names");

  std::string distType1 ("tnorm");  // Available pdf's Andrew has implemented
  std::string distType2 ("beta_lu");
  std::string distType3 ("gamma_l");
  std::string distType4 ("lnorm_l");
  std::string distType5 ("constant");

  double p1, p2, lower, upper, x;

  // Diff. EAMs has diff. numbers of parameter and each item in the priorList
  // contains the dist. setting for a parameter, so I use priorList.size() to
  // get the size of the priorList; thereby each of them is calculated one after
  // another.
  for(int i=0; i < priorList.size(); i++)
  {
    l1 = priorList[parameterNames[i]];
    std::string distName = l1.attr("dist");

    p1 = Rcpp::as<double>(l1[0]);  // parameter1: mean; shape1; shape; meanlog
    p2 = Rcpp::as<double>(l1[1]);  // parameter2: sd;   shape2; scale; sdlog

    if (distName.compare(distType1) == 0) {         // tnorm
      lower = Rcpp::as<double>(l1[2]);
      upper = Rcpp::as<double>(l1[3]);
      x = rtn_scalar(p1, p2, lower, upper);
    } else if (distName.compare(distType2) == 0) {  // beta_ul
      x = R::rbeta(p1, p2);
    } else if (distName.compare(distType3) == 0) {  // gamma_l
      x = R::rgamma(p1, p2);
    } else if (distName.compare(distType4) == 0) {  // lnorm_l
      x = R::rlnorm(p1, p2);
    } else if (distName.compare(distType5) == 0) {  // constant
      x = 1;
    } else {
      Rcpp::Rcout << "Distribution type not yet supported" << "\n";
      x = 0;
    }
    out[i] = x;
  }

  out.attr("names") = parameterNames;
  return out;
}

//' Generate Random Numbers from Prior Probability Distribution
//'
//' \code{dprior} matches 5 types of string: \code{tnorm}, \code{beta_lu},
//' \code{gamma_l}, \code{lnorm_l}, and \code{constant} to determine which
//' density functions to call. \code{dprior} calls beta, gamma, log-normal
//' density functions via R API.  For truncated normal density, \code{dprior}
//' calls \code{dtn_scalar}, an internal Rcpp function built specific for
//' ggdmc. Whetehr log the probability density is determined by the boolean
//' \code{log} sent in via p.prior.  This function is akin to DMC's
//' \code{log.prior.dmc}.
//'
//' @param pPrior a p.prior list
//' @param n how many random number to generate
//' @return a double matrix with nrow equal to the number of random numbers and
//' ncol equal to the number of EAM parameters (npar)
//' @export
//' @examples
//' ddm.prior <- prior.p.dmc(
//' dists = c("tnorm", "tnorm", "beta", "tnorm", "beta", "beta"),
//'   p1    = c(a = 1, v = 0, z = 1, sz = 1, sv = 1, t0 = 1),
//'   p2    = c(a = 1, v = 2, z = 1, sz = 1, sv = 1, t0 = 1),
//'   lower = c(0,-5, NA, NA, 0, NA),
//'   upper = c(2, 5, NA, NA, 2, NA))
//'
//' view(ddm.prior)
//' ##      mean sd lower upper log    dist  untrans
//' ##   a     1  1     0     2   1   tnorm identity
//' ##   v     0  2    -5     5   1   tnorm identity
//' ##   z     1  1     0     1   1 beta_lu identity
//' ##   sz    1  1  -Inf   Inf   1   tnorm identity
//' ##   sv    1  1     0     2   1 beta_lu identity
//' ##   t0    1  1     0     1   1 beta_lu identity
//'
//' rprior(ddm.prior, 9)
//' ##               a           v         z         sz        sv         t0
//' ## [1,] 0.97413686  0.78446178 0.9975199 -0.5264946 0.5364492 0.55415052
//' ## [2,] 0.72870190  0.97151662 0.8516604  1.6008591 0.3399731 0.96520848
//' ## [3,] 1.63153685  1.96586939 0.9260939  0.7041254 0.4138329 0.78367440
//' ## [4,] 1.55866180  1.43657110 0.6152371  0.1290078 0.2957604 0.23027759
//' ## [5,] 1.32520281 -0.07328408 0.2051155  2.4040387 0.9663111 0.06127237
//' ## [6,] 0.49628528 -0.19374770 0.5142829  2.1452972 0.4335482 0.38410626
//' ## [7,] 0.03655549  0.77223432 0.1739831  1.4431507 0.6257398 0.63228368
//' ## [8,] 0.71197612 -1.15798082 0.8265523  0.3813370 0.4465184 0.23955415
//' ## [9,] 0.38049166  3.32132034 0.9888108  0.9684292 0.8437480 0.13502154
// [[Rcpp::export]]
Rcpp::NumericMatrix rprior(Rcpp::List pPrior, int n)
{
   int nrow = n ;
   int ncol = pPrior.size() ;

   Rcpp::NumericMatrix m = na_matrix(nrow, ncol) ;
   Rcpp::NumericVector a ;

  for (int i = 0; i < nrow; i++)
  {
    a = rprior_scalar(pPrior);
    for (int j = 0; j < ncol; j++)
    {
      m(i, j) = a[j] ;
    }
  }
  Rcpp::CharacterVector pNames = pPrior.names() ;
  colnames(m) = pNames;
  return m ;
}

