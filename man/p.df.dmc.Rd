% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dmc_model.R
\name{p.df.dmc}
\alias{p.df.dmc}
\title{Gets Parameter Data Frame}
\usage{
p.df.dmc(model, p.vector, cell, n1order = TRUE)
}
\arguments{
\item{model}{a model, created by model.dmc}

\item{p.vector}{a parameter vector}

\item{cell}{a string indicating a design cell, e.g., \code{s1.f1.r1}}

\item{n1order}{a boolean switch to use specific LBA ordering for its
parameters}
}
\value{
rows in natural (r1, r2 etc., used by simulate.dmc) or "n1" order
(used by likelihood.dmc)
}
\description{
Gets parameter data frame (one row for each accumulator) for a design cell
(specified by name or index) with model picking out the appropriate
elements of par, and function transform.par transforms them appropriately
for model specified in model.
}
\examples{
m1 <- model.dmc(
  p.map     = list(a="1",v="F",z="1",d="1",sz="1",sv="1", t0="1",st0="1"),
  constants = c(st0=0,d=0),
  match.map = list(M=list(s1="r1",s2="r2")),
  factors   = list(S=c("s1","s2"), F=c("f1", "f2")),
  responses = c("r1","r2"),
  type      = "rd")
pVec <- c(a=1, v.f1=1, v.f2=1.5, z=0.5, sz=0.25, sv=0.2,t0=.15)
p.df.dmc(m1, pVec, "s1.f1.r1", n1order=TRUE)
}

