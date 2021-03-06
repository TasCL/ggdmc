% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dmc_plotting.R
\name{plot_priors}
\alias{plot_priors}
\title{Plot Prior Probability Density}
\usage{
plot_priors(p.prior, save.dat = FALSE, ...)
}
\arguments{
\item{p.prior}{a prior parameter list, created usually by prior.p.dmc}

\item{save.dat}{a boolean switch for save data frame, if the user wants to
modifiy the plotting parameter}

\item{...}{other arguments}
}
\description{
Plot each the prior density distribution for each EAM parameter and then
plot them all together in a canvas.
}
\examples{
pop.prior <- prior.p.dmc(
             dists = rep("tnorm", 7),
             p1=c(a=2,  v.f1=4, v.f2=3, z=0.5,sv=1, sz=0.3, t0=0.3),
             p2=c(a=0.5,v.f1=.5,v.f2=.5,z=0.1,sv=.3,sz=0.1,t0=0.05),
             lower=c(0,-5, -5, 0, 0, 0, 0),
             upper=c(5, 7,  7, 1, 2, 1, 1))

view(pop.prior)
plot_priors(pop.prior)
d <- plot_priors(pop.prior, save.dat=TRUE)

require(ggplot2)
p2 <- ggplot(d, aes(x = xpos, y = ypos)) + geom_line() +
      facet_wrap(~gpvar, scales="free") + theme_bw(base_size =14)
}
\keyword{plot_priors}

