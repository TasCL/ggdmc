% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dmc_hierarchical.R
\name{h.run.dmc}
\alias{h.run.dmc}
\title{Fit an EAM with Multiple Participants}
\usage{
h.run.dmc(samples, report = 100, cores = 1, blocks = NA, p.migrate = 0,
  h.p.migrate = 0, force.hyper = NA, force.data = NA, gamma.mult = 2.38,
  h.gamma.mult = NA, setting = NULL, verbose = FALSE)
}
\arguments{
\item{samples}{a DMC samples, usually generated by calling
\code{h.samples.dmc}.}

\item{report}{the progress report. Default is every 100 iterations report
once}

\item{cores}{a switch for parallel computation.  In the case of fixed-effect
model, \code{h.run.dmc} calls \pkg{snow} (via \pkg{snowfall}) or
\pkg{parallel}. When setting for example \code{cores=4}, \code{h.run.dmc}
requests 4 parallel R processes. In the case of random-effect model,
\code{h.run.dmc} calls OpenMP library, which will automatically decide the
number of cores to use, so setting core number higher than 1 has the same
effect.  Use OpenMP only when the data set is large (e.g., more than 1000
trial per condition).}

\item{blocks}{a DEMCMC parameter for blocking. Not implemented yet. Default
is NA}

\item{p.migrate}{the probability of applying migration sampler. Set it
greater than will trigger the migration sampler. For example
\code{p.migrate=0.05} will use migration sampler about 5\% chance.}

\item{h.p.migrate}{similar to p.migrate for hyper level}

\item{force.hyper}{This argument has no funciton. It is here only for
compatibility (with DMC) reason.}

\item{force.data}{This argument has no funciton. It is here only for
compatibility (with DMC) reason.}

\item{gamma.mult}{a DE-MCMC tuning parameter, affecting the size of jump.
Default is 2.38.}

\item{h.gamma.mult}{Similar with gamm.mult.  This argument has no funciton.
It is here only for compatibility reason.}

\item{setting}{a list to store setting arguments, such as p.migrate,
gamma.mult, etc, sending to C++ function. The user should not use this
argument}

\item{verbose}{a swtich to see debugging information}
}
\description{
\code{h.run.dmc} calls \code{run_data} in the case of fixed-effect model, or
\code{run_hyper} in the case of random-effect model. At the moment, it fits
either fixed-effect or random-effect model by looking for the
\code{hyper} attribute in a DMC sample/object.
}
\examples{
m1 <- model.dmc(
p.map     = list(a="1",v="1",z="1",d="1",sz="1",sv="1",t0="1",st0="1"),
match.map = list(M=list(s1="r1", s2="r2")),
factors   = list(S=c("s1", "s2")),
constants = c(st0=0, d=0),
responses = c("r1","r2"),
type      = "rd")

## Population distribution
pop.prior <- prior.p.dmc(
  dists = rep("tnorm", 6),
  p1    = c(a=2,  v=2.5, z=.5, sz=.3, sv=1,  t0=.3),
  p2    = c(a=.5, v=.5,  z=.1, sz=.1, sv=.3, t0=.05),
  lower = c(0,-5, 0, 0, 0, 0),
  upper = c(5, 7, 2, 2, 2, 2))

dat <- h.simulate.dmc(m1, nsim=30, p.prior=pop.prior, ns=8)
mdi <- data.model.dmc(dat, m1)
p.prior  <- prior.p.dmc(
  dists = rep("tnorm", 6),
  p1    = c(a=2,  v=2.5, z=.5, sz=.3, sv=1,  t0=.3),
  p2    = c(a=.5, v=.5,  z=.1, sz=.1, sv=.3, t0=.05) * 5,
  lower = c(0,-5, 0, 0, 0, 0),
  upper = c(5, 7, 2, 2, 2, 2))

## Fixed-effect model
samplesInit <- h.samples.dmc(nmc=20, p.prior=p.prior, data=mdi, thin=1)
samples0    <- h.run.dmc(samples=samplesInit, report=10)

## Make a hyper-prior list
mu.prior <- prior.p.dmc(
  dists = rep("tnorm", 6),
  p1    = c(a=2,  v=2.5, z=.5, sz=.3, sv=1,  t0=.3),
  p2    = c(a=.5, v=.5,  z=.1, sz=.1, sv=.3, t0=.05) * 5,
  lower = c(0,-5, 0, 0, 0, 0),
  upper = c(5, 7, 2, 2, 2, 2))

sigma.prior <- prior.p.dmc(
  dists = rep("beta", 6),
  p1    = c(a=1, v=1, z=1, sz=1, sv=1, t0=1),
  p2    = c(1,1,1,1,1,1),
  upper = c(2,2,2,2,2,2))

pp.prior <- list(mu.prior, sigma.prior)

## Random-effect model
hsamplesInit <- h.samples.dmc(nmc=20, p.prior=p.prior, pp.prior=pp.prior,
  data=mdi, thin=1)
hsamples0 <- h.run.dmc(samples=hsamplesInit, report=10, p.migrate=.05,
  h.p.migrate=.05)
}
\keyword{h.run.dmc}

