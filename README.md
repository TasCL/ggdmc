# DMC with Graphic Grammar, Parallel Computing and C++ Capabilities

ggdmc implements hierarchical Bayesian model of evidence accumulation 
model (EAM). This release includes drift-diffusion model (Ratcliff, 1978).  
This DDM implements Voss, Rothermund, and Voss's (2004) equations based on 
their fast-dm 30.2 C functions as C++ functions. 

ggdmc is derived from Andrew Heathcote's Dynamiuc Model of Choice (DMC), 
which has also implemented numerous other EAMs.  Similar to DMC, ggdmc uses 
differential evolution Markov Chain Monte Carlo sampler to search optimal 
theta and phi that maximise posterior likelihood. See NEWS for new features

## Quick Getting Started
Below shows a simple example extract from Andrew Heathcote's DMC workshop 
materials. For further details, 

```
require(ggdmc) 

## Use a 6-parameter drift-diffusion model  
## This only experimental factor is Stimulus type
m1 <- model.dmc(
  p.map     = list(a="1", v="1", z="1", d="1", sz="1", sv="1", t0="1",st0="1"),
  constants = c(st0=0, d=0),  ## Fixed st0 and d at 0
  match.map = list(M = list(s1="r1", s2="r2")),
  factors   = list(S=c("s1", "s2")),  ## This only experimental factor is 
  responses = c("r1","r2"),           ## Two response types
  type      = "rd")                   ## rd stands for Ratcliff's diffusion

## Set up a prior probability distribution setting
p.prior <- prior.p.dmc(
  dists = rep("tnorm", 6),
  p1    = c(a=2,  v=2.5, z=.5, sz=.3, sv=1,  t0=.3),
  p2    = c(a=.5, v=.5,  z=.1, sz=.1, sv=.3, t0=.05),
  lower = c(0,-5, 0, 0, 0, 0),
  upper = c(5, 7, 2, 2, 2, 2))

## Assume a true DDM parameter vector
p.vector <- c(a=1, v=1, z=.5, sz=.25, sv=0.2,t0=.15)

## Use simulate function to simulate choice-RT data based on a
dat1 <- simulate(m1, nsim=1e2, p.vector=p.vector)

## Set up a data model instace. This bind the data with the model set-up
mdi1 <- data.model.dmc(dat1, m1)

## Initialise a small sample 
## (1) iteration # == 250 
## (2) thin length == 1
## (3) prior distributions are listed in p.prior 
## (4) data == mdi1, an assumed model and a simulated/empirical data frame
samples0 <- samples.dmc(nmc=250, p.prior=p.prior, data=mdi1, thin=1)

## Fit the Bayesian model 
samples0 <- run.dmc(samples0)


## Check if model converged, etc.
gelman.diag.dmc(samples0)
## Potential scale reduction factors:
## 
##    Point est. Upper C.I.
## a        1.16       1.29
## v        1.20       1.36
## z        1.18       1.33
## sz       1.12       1.21
## sv       1.24       1.48
## t0       1.16       1.34

plot(samples0)  ## Check traceplot to see the fit

## Further details and more arguments, please see
?gelman.diag.dmc
?samples.dmc
?run.dmc
?plot.dmc

```


### Prerequisities
 -  R (>= 3.0.2)
 -  Rcpp (>= 0.12.3), RcppArmadillo (>= 0.6.700.6.0), ggplot2 (>= 2.1.0),
 rtdists (>= 0.6-6), gridExtra (>= 2.2-1), ggmcmc (>= 0.7.3), 
 ggthemes (>= 3.0.1),  snowfall (>= 1.84.6-1), stats (>= 3.2.2), 
 hypergeo (>= 1.2-13), pracma (>= 1.8.8), statmod (>= 1.4.24), loo (>= 0.1.6) 
 - coda (>= 0.16-1)


### Installing

```
install.packages("ggdmc_0.1.3.5.tar.gz", repos = NULL, type="source")

```

## Contributors

The ggdmc C++ codes are developed by [Yi-Shin Lin](http://www.tascl.org/yi-shin-lin.html). 


The DMC R codes are developed by [Andrew Heathcote](http://www.tascl.org/andrew-heathcote.html) 

## License

GPL-2 

## Acknowledgments

* density.cpp is based on Voss & Voss's (2012) density.c in fast-dm 30.2. Any 
identical parts in density.cpp belong to their copyright.
* The tnorm_wrapper.cpp and tnorm.R are based on Jonathan Olmsted's
<jpolmsted@gmail.com> RcppTN 0.1-8 (https://github.com/olmjo/RcppTN) 
and Christopher Jackson's <chris.jackson@mrc-bsu.cam.ac.uk> R codes in msm. 
Any identical parts in tnorm_wrapper.cpp and tnrom.R belong to their copyright.
* C++ codes in this package depend largely on Dirk Eddelbuettel, Romain 
Francois and Doug Bates's Rcpp, and RcppArmadillo.  
* Armadillo is a collection of C++ library for performing linear
algebra <http://arma.sourceforge.net/>, authored by Conrad Sanderson. 
* DEMCMC sampler is based on the following papers, Turner & Sederberg (2012),
Turner et al (2013) and Ter Braak (2006).
* Thanks to Matthew Gretton's consulation about rtdists's intenal mechanism
