% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dmc_analysis.R
\name{summary.dmc}
\alias{summary.dmc}
\title{Summarise a DMC Sample with One Participant}
\usage{
\method{summary}{dmc}(object, start = 1, end = NA, ...)
}
\arguments{
\item{object}{a model samples}

\item{start}{summarise from which MCMC iteration. Default uses the first
iteration.}

\item{end}{summarise to the end of MCMC iteration. For example, set
\code{start=101} and \code{end=1000}, instructs the function to calculate
from 101 to 1000 iteration. Default uses the last iteration.}

\item{...}{other arguments}
}
\description{
Call coda package to summarise the model parameters in a DMC samples
}
\examples{
m1 <- model.dmc(
p.map=list(a="1",v="1",z="1",d="1",sz="1",sv="1", t0="1",st0="1"),
constants=c(st0=0,d=0),
match.map=list(M=list(s1="r1",s2="r2")),
factors=list(S=c("s1","s2")),
responses=c("r1","r2"),
type="rd")

p.prior <- prior.p.dmc(
  dists = rep("tnorm", 6),
  p1=c(a=2,   v=2.5, z=0.5, sz=0.3, sv=1,  t0=0.3),
  p2=c(a=0.5, v=.5,  z=0.1, sz=0.1, sv=.3, t0=0.05),
  lower=c(0,-5, 0, 0, 0, 0),
  upper=c(5, 7, 2, 2, 2, 2))

pVec <- c(a=1,v=1, z=0.5, sz=0.25, sv=0.2,t0=.15)

dat1 <- simulate(m1, nsim=30, p.vector=pVec)
mdi1 <- data.model.dmc(dat1, m1)

samples0 <- samples.dmc(nmc=100, p.prior=p.prior, data=mdi1)
samples0 <- run.dmc(samples0, p.migrate=.05)
gelman.diag.dmc(samples0)
class(samples0)
## [1] "dmc"
summary(samples0)
}
\seealso{
\code{\link{summary.dmc.list}, \link{summary.hyper}}
}
\keyword{summary}

