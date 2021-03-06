% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dmc_sampling.R
\name{prior.p.dmc}
\alias{prior.p.dmc}
\title{Makes a list of prior distribution parameters.}
\usage{
prior.p.dmc(p1, p2, lower = rep(NA, length(p1)), upper = rep(NA,
  length(p1)), dists = rep("tnorm", length(p1)), untrans = rep("identity",
  length(p1)), dist.types = c("tnorm", "beta", "gamma", "lnorm", "constant"))
}
\arguments{
\item{p1}{the values of location parameters for each prior distribution, set
as a double vector}

\item{p2}{ditto for scale parameter vector}

\item{lower}{lower support boundary}

\item{upper}{upper support boundary}

\item{dists}{indicate which prior distribution, e.g., uniform, beta etc.}

\item{untrans}{whether do log transformation or not. Default is identity,
namely not to transform}

\item{dist.types}{allowed prior distributions in current version of DMC}
}
\description{
\code{prior.p.dmc} creates a list of prior distribution an array
object ("model") with a set of attributes
specifying a particular model and parameterization. Call \pkg{coda} to
summarise the model parameters in a DMC samples with multiple participants
at the hyper level.
}
\keyword{prior.p.dmc}

