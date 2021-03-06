% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dmc_sampling.R
\name{dtnorm}
\alias{dtnorm}
\alias{rtnorm}
\title{Truncated Normal Distribution}
\usage{
dtnorm(x = 0, mean = 0, sd = 1, lower = -Inf, upper = Inf, log = 0,
  checks = FALSE)

rtnorm(n = 1, mean = 0, sd = 1, lower = -Inf, upper = Inf,
  checks = TRUE)
}
\arguments{
\item{x, n}{x in \code{dtnorm} is a vector of quantiles; n in
\code{rtnorm} indicates how many random number to generate}

\item{mean}{vector of means}

\item{sd}{vector of standard deviations}

\item{lower}{lower truncation point. Default as -Inf}

\item{upper}{upper truncation point. Default as Inf.}

\item{log}{default 0 (FALSE). Enter 1 to make it TRUE. Whether to
calculate density. Use only in dtnorm.}

\item{checks}{a switch to turn on check functions to check inputs and ouput.
Default to FALSE for the sake of speed.}
}
\description{
The two wrapper functions call Rcpp functions to speed up computation. The
codes are based on Christopher Jackson's \pkg{msm} (1.5) package and Jonathan
Olmsted's \pkg{RcppTN} (0.1-8) package.
}
\details{
\code{dtnorm} calculates probability density for a truncated normal
distribution with mean equal to mean and standard deviation equal to
sd before truncation, and truncated on the interval [lower, upper].

\code{rtnorm} generates a random number from a truncated Normal
distribution with mean equal to mean and standard deviation equal
to sd before truncation, and truncated on the interval [lower, upper].
}
\examples{
## Use dtnorm and rtnorm with their default values
dtnorm()
rtnorm()

## A similar curve plotting example extracted from dnorm functions
plot(function(x) dnorm(x, log = FALSE), -2.5, 2.5,
     main = "Normal Distribution", ylim=c(0,0.45), ylab="Density")
curve(dtnorm(x, lower=-2, upper=2), add=TRUE, col="tomato", lwd=2)
mtext("dnorm(x)", adj = 0)
mtext("dtnorm(x)", col = "tomato", adj = 1)
}

