#' Supersonic DMC
#'
#' \pkg{ggdmc} implements graphic grammar, Rcpp, and Open MPI for
#' Andrew Heathcote's Dynamic Models of Choice (DMC), an algorithm for fitting
#' various evidence accumulation models via differential evolution MCMC
#' sampler. Currently, ggdmc implements only drift-diffusion model.
#'
#' @keywords package
#' @name ggdmc
#' @docType package
#' @author  Yi-Shin Lin <yishin.lin@utas.edu.au> \cr
#' Andrew Heathcote <andrew.heathcote@utas.edu.au>
#' @references Turner, B. M., & Sederberg P. B. (2012). Approximate
#' Bayesian computation with differential evolution,
#' \emph{Journal of Mathematical Psychology}, \bold{56}, 375--385. \cr\cr
#' Ter Braak (2006). A Markov Chain Monte Carlo version of the genetic
#' algorithm Differential Evolution: easy Bayesian computing for real
#' parameter spaces. \emph{Statistics and Computing}, \bold{16}, 239-249.
#' @importFrom Rcpp evalCpp
#' @useDynLib ggdmc
NULL
