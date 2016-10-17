# DMC with Graphic Grammar, Parallel Computing and C++ Capabilities

ggdmc implements the computation of hierarchical Bayesian model fitting 
evidence accumulation model (EAM). This release includes only drift-diffusion
model (Ratcliff, 1978) and LBA model (Brown & Heathcote, 2008).  This former
replements Voss, Rothermund, and Voss's (2004) equations and their 
fast-dm 30.2 C functions as C++ functions. 

ggdmc is derived from Andrew Heathcote's Dynamiuc Model of Choice (DMC), 
which has implemented other types of EAMs.  Similar to DMC, ggdmc uses 
differential evolution Markov Chain Monte Carlo sampler to search optimal 
theta and phi sets that maximise posterior likelihood. 


Version 0.1.3.4 New features:
  See NEWs
  
Version 0.1.2.0 New features:
  * random-effect fit with DE sampler
  * Dynamic choose of samplers, using C++ class
  * R convenient function, view, to see p.prior and pp.prior easier  
  * example scripts in examples folder
  
Version 0.1.1.3 New features:
  * traceplot.dmc to get trace and density plots at one go
  * run function replaces run.dmc to run a bit faster 
  * truncated normal functions, _dtn_ and _rtn_ run in C++ mode. _dtn_ and _rtn_
  run faster than msm's dtnorm and rtnorm. 
  * constant drift diffusion functions, "g_plus" and "g_minus". You can call
  ddiffusion's equations directly without loading rtdists.  Note g_plus and
  g_minus expect a parameter vector (p.vector+precision).   
  * ggdmc now has its own initiation function (_initialise_). It acts as 
  DMC's sample.dmc. Further, it helps you to order parameters in p.prior 
  (ggdmc calls it pList, as prior parameter list).  

## Getting Started
Below is an example of fitting fixed-effect model. ggdmc has not incorporated  
DMC model building and simulation functions, so this example still loads msm
and rtdists packages. That is, DMC requires them.

Firstly, follow the usual procedure to set up a model and simulate data.  
Note the naming sequence in defining pop.mean, pop.scale.  ggdmc uses a 
function, called initialise, to order them.  If you use samples.dmc or 
h.samples.dmc to initialise a samples.  The parameter order in models' p.vector
might not be the same as the order in p.prior and/or pp.prior (depending on
how you set up p.prior and pp.prior). 

Specifically, most of the time DMC uses p.vector attribute in the model as 
p.names.  You should see a vector with usually NA values. The vector's name is  
often used as a reference for the parameter sequence. 

```
rm(list=ls()); 
library(ggdmc); 
source ("tutorial/file_utils.R")
load_model ("ddm.R") 

model <- model.dmc(
  p.map     = list(a="1",v="F",z="1",d="1",sz="1",sv="F",t0="1",st0="1"),
  match.map = list(M=list(s1="r1",s2="r2")),
  factors=list(S=c("s1","s2"), F=c("f1","f2")),
  constants = c(st0=0,d=0),
  responses = c("r1","r2"),
  type = "rd")
 
# Population distribution
pop.mean <- c(a=1,  v.f1=1, v.f2=.2, z=0.5, sz=0.3, 
  sv.f1=0.25, sv.f2=0.23, t0=0.3)
pop.scale <-c(a=0.2,v.f1=.2, v.f2=.2, z=0.1, sz=0.05, 
  sv.f1=.05, sv.f2=.05,t0=0.05)
pop.prior <- prior.p.dmc(
  dists = rep("tnorm", 8),
  p1=pop.mean,                           
  p2=pop.scale,
  lower=c(0,-5, -5, 0, 0, 0, 0,0),
  upper=c(2, 5, 5, 1, 2, 2, 1, 1)
)

## Simulate some data
## This example simulates 10 participants. In the following random-effect fit
## more participants will help parameter recovery. See example4. 
raw.data   <- h.simulate.dmc(model, p.prior = pop.prior, n = 150, ns = 10)
data.model <- data.model.dmc(raw.data, model)
ps <- round( attr(raw.data, "parameters"), 2); 

## Set up also hyper level to demo random-effect fits. 
mu.prior <- prior.p.dmc(
  dists = rep("tnorm", 8),
  p1=pop.mean,
  p2=pop.scale*6,
  lower=c(0,-5, -5, 0, 0, 0, 0,0),
  upper=c(2, 5, 5, 1, 2, 2, 1, 1)
)

sigma.prior <- prior.p.dmc(
  dists = rep("beta", length(pop.prior)),
  p1=c(a=1, v.f1=1,v.f2=1, z=1, sz=1, sv.f1=1, sv.f2=1,  t0=1),
  p2=c(1,1,1,1,1,1,1,1)
)

pp.prior <- list(mu.prior, sigma.prior)

## Use DMC's routines
samples0 <- h.samples.dmc(nmc=500, pop.prior, data.model, thin = 3)
samples1 <- add.hyper.dmc(samples0, pp.prior = pp.prior, 
  start.from = 1, thin=3)

```

The setting argument is optional. If not set, run function will choose a default
set, as listed below. When fitting fixed-effect, ggdmc use the parallel pacakge 
to allocate each individual fit to all the available cores.  Setting cores to 
more than one will trigger the parallel package.

```
## Use ggdmc's initialise to replace h.samples.dmc and add.hyper.dmc at one go
settingi0 <- list(add=FALSE, rp=0.001, thin=1, start.from=1, restart=TRUE,
                   verbose=TRUE, remove=FALSE)
hsamples0 <- initialise(nmc=1000, data=data.model, p.prior=p.prior, 
             pp.prior=pp.prior, setting=settingi0)
                         
```

Perform a fixed-effect model fit for one participant. Take ~60 seconds

```
settingr0      <- list(p.migrate=0.5, gamma.mult=2.38, h.gamma.mult=2.38,  
                       report=200, cores=1)
samples_one    <- run(samples0[[1]])
samples_one_cp <- run(samples0[[1]], setting=settingr0)  ## identical

# The reporting number has been corrected.
# Run fixed-effect fit ...
# Only one subject to be fit.
# 100 200 300 400 500    
#  user    system elapsed 
#  60.344   0.012  60.399 
 
```

Perform a fixed-effect model fit for 6 participants. This takes quite a
while about 9.1 mins.

```
samples_multiple <- run(samples0)

# Run fixed-effect fit ...
# Fitting 6 participants
# Participant 1: 100 200 300 400 
# Participant 2: 100 200 300 400 
# Participant 3: 100 200 300 400 
# Participant 4: 100 200 300 400 
# Participant 5: 100 200 300 400 
# Participant 6: 100 200 300 400 
#    user  system elapsed 
# 545.872   0.088 545.924 


```

### Fitting random-effect model
Please see example3, 4, and 5 in examples folder.


### Diagnosis
Using traceplot.dmc to examine the results. You can use lapply to plot all 
participants at once, but RStudio may hang. These plotting routines take up 
a lot of memory, becase data are huge.  Load required packages. The dependency
should be resolved while installing them.

```
library(ggmcmc); library(ggthemes) ;library(gridExtra)
traceplot.dmc(samples_one)

# use lapply to plot all at once
lapply(samples_multiple, traceplot.dmc)

# or plot each participant one by one
traceplot.dmc(samples_multiple[[1]])
traceplot.dmc(samples_multiple[[2]])
traceplot.dmc(samples_multiple[[3]])
traceplot.dmc(samples_multiple[[4]])
traceplot.dmc(samples_multiple[[5]])
traceplot.dmc(samples_multiple[[6]])

```

### Check against DMC
ggdmc took ~13 mins to fit 10 participants with 500 iterations and thinning 
length = 5. DMC took ~ 75 mins with 4 cores.  This test was conducted on Intel
core i7.   


```
system.time(run1.run <- run(samples0))
# user  system elapsed 
# 782.096   0.000 781.192 

samples0 <- h.samples.dmc(nmc=500, pop.prior, data.model, thin = 5)
system.time(run1.run.dmc <- h.run.dmc(samples0, p.migrate = 0, cores=4))
# user   system  elapsed 
# 11.340    6.372 4507.278 
```

### Example for using _initialise_ funciton

```
rm(list = ls())
library(ggdmc)
load("data/dmc_4_7_fixed.RData")

##### Case 1 #####
samples_s1_0 <- initialise(nmc=100, pList = samples[[1]]$p.prior, data = samples[[1]]$data)

##### Case 2 #####
removeVector <- c(1:50)
setting1 <- list(rp=0.001, thin=3, start.from=1, n.chains=samples[[1]]$n.chains, 
                 remove=removeVector, add=FALSE, restart=FALSE, verbose=TRUE)

##### Case 3 #####
## Nothing been removed, because the user has not entered a samples
samples_s1_1_1 <- initialise(nmc=100, pList = samples[[1]]$p.prior, 
                  data = samples[[1]]$data, setting=setting1)

##### Case 4 #####
## You must enter a samples to remove
samples_s1_1_2 <- initialise(nmc=100, pList = samples[[1]]$p.prior, 
                  data = samples[[1]]$data, samples=samples[[1]],  
                  setting=setting1)

##### Case 5 #####
## When adding new nmc, remove must be FALSE 
setting2 <- list(rp=0.001, thin=3, start.from=1, n.chains=samples[[1]]$n.chains, 
                 remove=FALSE, add=TRUE, restart=FALSE, verbose=TRUE)

samples_s1_2 <- initialise(nmc=50, pList = samples[[1]]$p.prior, 
                data = samples[[1]]$data, samples=samples[[1]],  
                setting=setting2)


##### Case 6 #####
## add and remove must be FALSE, if you want to restart
setting3 <- list(rp=0.001, thin=3, start.from=51, n.chains=samples[[1]]$n.chains, 
                 remove=FALSE, add=FALSE, restart=TRUE, verbose=TRUE)

samples_s1_3 <- initialise(nmc=50, pList = samples[[1]]$p.prior, 
                data = samples[[1]]$data, samples=samples[[1]],  
                setting=setting3)
                
```

### Parallel Processing 
See examples folder.

### LBA C++ density working in progress 
Show the speed up of two LBA-related function in rtdists while I am waiting 
the snowfall and parallel memory tests.

```
library(Rcpp)
sourceCpp("tests/dlba_norm.cpp")

## rtdists pnormP and dnormP functions
pnormP2 <- function(x,mean=0,sd=1,lower.tail=TRUE) {
  ifelse(abs(x)<7, pnorm(x, mean=mean, sd=sd,lower.tail=lower.tail),
    ifelse(x<0, 0, 1))
}

dnormP2 <- function(x,mean=0,sd=1) {
  ifelse(abs(x)<7, dnorm(x,mean=mean,sd=sd), 0)
}

rt <- -100:100
all(pnormP(rt) == pnormP2(rt))
all(dnormP(rt) == dnormP2(rt))

library(microbenchmark)
microbenchmark(pnormP(rt), pnormP2(rt))
# Unit: microseconds
# expr           min      lq     mean median      uq     max neval cld
# pnormP(rt)   3.423  4.2965  5.55685  5.309  6.1120  23.259   100   a 
# pnormP2(rt) 68.656 81.4020 82.63640 82.764 83.8465 135.844   100   b

microbenchmark(dnormP(rt), dnormP2(rt))
# Unit: microseconds
# expr           min      lq     mean  median      uq     max neval cld
# dnormP(rt)   2.794  3.4580  4.16349  3.7375  4.1565  23.258   100   a 
# dnormP2(rt) 36.598 37.5405 40.12995 38.1350 39.0420 114.681   100   b

```


### Prerequisities
This release uses Rcpp, RcppArmadillo and C++11 libraries. So you will need
Rcpp (>= 0.12.3) and RcppArmadillo.  RcppArmadillo requires g++ >= 4.6 to 
compile. cl uses g++/gcc 4.4.7, which is earlier than 4.6.  So until cl upgrades
to a g++ later than 4.4.7, RcppArmadillo will not compile. 

~~Note that you need to install RcppArmadillo before adding the c++11 flag. 
Installing RcppArmadillo will fail, if you flag g++ compiler to use c++11.~~
(Resolved. A local Makevars has been added in src folder.)

traceplot.dmc requires a large collection ggplot2-related packages, so ggdmc
will not install them for you.  Major ones are ggmcmc, ggthemes and
gridExtra. ggmcmc depends on GGally, which depends on a defunct package, grid.
Therefore, to sort out GGally's problem, you need to install
GGally from its github <https://github.com/ggobi/ggally>, instead of CRAN. 
This will requires devtools package, which depends on another huge collection 
of C and R packages to facilitate installation github-based packages. For the 
C parts, you will need libcurl4-openssl-dev and libssl-dev.

ggdmc has been tested on Ubuntu 14.04.4 LTS, Ubuntu 15.04, Ubuntu 15.10.  All 
are 64-bit systems.  It has also been tested on Windows 10, and so far there 
are still some hurdles waiting to be resolved on Windows system. OS X should 
work, because it is also Unix-like system. 

If you want to use DMC's loo facility, you need to install loo too. ggdmc will 
not do that for you.

### Installing

The usual way to install an R package from a tar.gz source. 

```
install.packages("ggdmc_0.1.2.0.tar.gz", repos = NULL, type="source")

```

## Others
C++ rtnorm pits against R's rtnorm. 1.381 vs 1.795 s

```
## simulate1 uses msm::rtnorm()
## simulate2 uses ggdmc::rtnorm()
rbenchmark::benchmark(replications=rep(1000, 3),
  data1 <- simulate1.dmc(p.vector,model,n=1),
  data2 <- simulate2.dmc(p.vector,model,n=1),
  columns=c('test', 'elapsed', 'replications'))


##                                             test elapsed
## 1 data1 <- simulate1.dmc(p.vector, model, n = 1)   1.795
## 3 data1 <- simulate1.dmc(p.vector, model, n = 1)   1.634
## 5 data1 <- simulate1.dmc(p.vector, model, n = 1)   1.743
## 2 data2 <- simulate2.dmc(p.vector, model, n = 1)   1.381
## 4 data2 <- simulate2.dmc(p.vector, model, n = 1)   1.399
## 6 data2 <- simulate2.dmc(p.vector, model, n = 1)   1.428
##   replications
## 1         1000
## 3         1000
## 5         1000
## 2         1000
## 4         1000
## 6         1000

```

## Contributors

The C++ codes are developed by [Yi-Shin Lin](http://www.tascl.org/yi-shin-lin.html) 
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
* DEMCMC sampler is based on Turner and Sederberg (2012), Turnder et al (2013) 
and Ter Braak (2006).
* Thanks to Matthew Gretton's consulation about rtdists's intenal mechanism
