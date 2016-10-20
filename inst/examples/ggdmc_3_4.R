##################  DMC Lesson 3: Sampling


### Lesson 3.4 Extra: Use C++ function to run sampling a single DDM subject


## Speed Difference Demo ----------------------------------------------
rm(list=ls()) 
## Current working directory must be set to the top-level folder  
## containing the dmc and tutorial subfolders
setwd("/media/yslin/OWL/Documents/ggdmc/tests/")
source ("tutorial/file_utils.R")
load_model ("ddm.R")

## Use the simulation data and settings
load("/media/yslin/OWL/Documents/ggdmc/examples/ggdmc_3_4.RData")

samples0 <- samples.dmc(nmc=100, p.prior=p.prior, data=mdi1)
system.time(samples0 <- run.dmc(samples0, p.migrate=.05, report=100))
##    user  system elapsed 
##   6.540   0.004   6.544 

samples0.cpp <- ggdmc::samples.dmc(nmc=100, p.prior=p.prior, data=mdi1)
system.time(samples0.cpp <- ggdmc::run.dmc(samples0.cpp, p.migrate=.05))

##    user  system elapsed
##   1.100   0.000   1.101 
## It is about 6 times quicker in one-subject model fit


rm(list=ls())
## Set up model and import/simulate data ------------------------------------
## To avoid name conflict, this lesson does not load DMC functions.
## See speed demo at the end of the script
library(ggdmc)
load("ggdmc_3_4.RData")

## Set up a Ratcliff's Diffusion (rd) model
m1 <- model.dmc(
  p.map     = list(a="1",v="1",z="1",d="1",sz="1",sv="1", t0="1",st0="1"),
  constants = c(st0=0, d=0),
  match.map = list(M=list(s1="r1", s2="r2")),
  factors   = list(S=c("s1", "s2")),
  responses = c("r1", "r2"),
  type      = "rd")

## Set up priror distributions for the 6 parameters
p.prior <- prior.p.dmc(
  dists = rep("tnorm", 6),
  p1    = c(a=2,   v=2.5, z=0.5,  sz=1,   sv=1,  t0=0.3),
  p2    = c(a=0.5, v=.5,  z=0.25, sz=0.1, sv=.3, t0=0.05),
  lower = c(0,-5, 0, 0, 0, 0),
  upper = c(5, 7, 2, 2, 2, 2))



## Check if all prior distributions are set appropriately
## 1. print out their content
view(p.prior)
##    mean   sd lower upper log  dist  untrans
## a     2  0.5     0     5   1 tnorm identity
## v   2.5  0.5    -5     7   1 tnorm identity
## z   0.5 0.25     0     2   1 tnorm identity
## sz    1  0.1     0     2   1 tnorm identity
## sv    1  0.3     0     2   1 tnorm identity
## t0  0.3 0.05     0     2   1 tnorm identity

## 2. plot all prior distributions
plot_priors(p.prior)

## Arbitrarily decide a set of true parameter values
p.vector <- c(a=1, v=1, z=0.5, sz=.25, sv=.32, t0=.15)


## Simulate data based on m1 model and true parameters set in p.vector
dat1 <- simulate(m1, nsim=1e2, p.vector=p.vector)


## Make a data model instance. Because this lesson fit only one subject,
## mdi1 is a data frame
mdi1 <- data.model.dmc(dat1, m1)

head(mdi1)
##    S  R        RT
## 1 s1 r2 0.6878078
## 2 s1 r1 0.1832346
## 3 s1 r2 0.1601428
## 4 s1 r1 1.2150433
## 5 s1 r1 0.2002556
## 6 s1 r1 0.3593219

# Accuracy around 66% and 76%
par(mfrow=c(1,2))
plot_cell_density(data.cell=mdi1[mdi1$S=="s1", ], C="r1", xlim=c(0,2))
plot_cell_density(data.cell=mdi1[mdi1$S=="s2", ], C="r2", xlim=c(0,2))
par(mfrow=c(1,1))

par(mfrow=c(2,3));
profile(mdi1, "a",   .1,  2, p.vector)
profile(mdi1, "v",   .1,  2, p.vector)
profile(mdi1, "z",   .2, .8, p.vector)
profile(mdi1, "sz",  .1, .9, p.vector)
profile(mdi1, "sv",  .1,  2, p.vector)
profile(mdi1, "t0", .01, .5, p.vector)
par(mfrow=c(1,1));


## Initialise a DMC sample
## (1) niter=100; (2) nthin=1 (default)
samples0 <- samples.dmc(nmc=100, p.prior=p.prior, data=mdi1)
names(samples0)
samples0$n.chains  ## n.chains == 18 (3 * n.pars)
samples0$n.pars    ## 6 DDM parameters
# [1] "theta"            "summed_log_prior" "log_likelihoods"
# [4] "data"             "p.prior"          "start"
# [7] "n.pars"           "p.names"          "rp"
# [10] "nmc"              "thin"             "n.chains"


## Demonstrate a few run.dmc arguments
## 1. Run model fit without migration. Default 100 iteration once
samples0 <- run.dmc(samples0)

## 2. If the user wants to use migration sampler, p.migrate=.05 indicates
## 5% of chance run.dmc will use migration sampler. Otherwise, it uses
## crossover sampler
samples0 <- run.dmc(samples0, p.migrate=.05)

## 3. Report more frequently
samples0 <- run.dmc(samples0, report=20)

## For more details about samples.dmc and run.dmc
?samples.dmc
?run.dmc

## Run Model Fit --------------------------------------------
samples0 <- samples.dmc(nmc=500, p.prior=p.prior, data=mdi1)

## Use 5% migration and burnin the first 500 iterations
## This will take about 10 s
samples0 <- run.dmc(samples0, p.migrate=.05)

## Turn off migration and get 1000 new iterations
samples1 <- run.dmc(samples.dmc(nmc=1000, p.prior=p.prior, samples=samples0))

## If the user wants to add more iteration onto an existing one
samples2 <- samples.dmc(nmc=1000, p.prior=p.prior, samples=samples1, add=TRUE)
samples2 <- run.dmc(samples2)



## Diagnosis --------------------------------------------
plot(samples2) ## Traceplot
plot(samples2, density=TRUE)  ## Trace- and density-plot
plot(samples2, density=TRUE, start=1001)      ## starting from 1001st iteraton
plot(samples2, pll.chain=TRUE) ## check posterior density


## Check Gelman and Rubin Convergence Diagnostic
gelman.diag.dmc(samples2)
## Potential scale reduction factors:
## 
##    Point est. Upper C.I.
## a        1.00       1.01
## v        1.00       1.01
## z        1.01       1.01
## sz       1.00       1.01
## sv       1.01       1.01
## t0       1.01       1.02
## 
## Multivariate psrf
## 
## 1.01


## Effective sample size
effectiveSize.dmc(samples2)
##    a    v    z   sz   sv   t0 
## 1653 1630 1725 1636 1709 1479

est <- summary(samples2); est
## Iterations = 1:2000
## Thinning interval = 1 
## Number of chains = 18 
## Sample size per chain = 2000 
## 
## 1. Empirical mean and standard deviation for each variable,
##    plus standard error of the mean:
## 
##      Mean       SD  Naive SE Time-series SE
## a  1.2464 0.062941 3.317e-04      1.591e-03
## v  1.4019 0.212905 1.122e-03      5.317e-03
## z  0.5309 0.019281 1.016e-04      4.677e-04
## sz 0.6550 0.056684 2.987e-04      1.433e-03
## sv 1.3355 0.271658 1.432e-03      6.706e-03
## t0 0.1541 0.002762 1.456e-05      7.855e-05
## 
## 2. Quantiles for each variable:
## 
##      2.5%    25%    50%    75%  97.5%
## a  1.1341 1.2014 1.2431 1.2873 1.3818
## v  0.9921 1.2574 1.3974 1.5456 1.8290
## z  0.4927 0.5181 0.5311 0.5440 0.5681
## sz 0.5385 0.6187 0.6567 0.6938 0.7602
## sv 0.7895 1.1476 1.3371 1.5312 1.8478
## t0 0.1472 0.1527 0.1546 0.1561 0.1580


## Estimations are OK with only 100 trials. 
p.vector
##    a    v    z   sz   sv   t0 
## 1.00 1.00 0.50 0.25 0.32 0.15 
round(est$quantiles[,3], 2)
##    a    v    z   sz   sv   t0 
## 1.24 1.40 0.53 0.66 1.34 0.15



# No signs of autocorrelation
acf.dmc(samples2)

# Confirm that no chains are stuck
pick.stuck.dmc(samples2, cut=10, verbose=TRUE)
## Deviation of mean chain log-likelihood from median of means
##     5     6    10    13    16     1     4     2    11     9    14     7    18 
##  0.34  0.33  0.26  0.18  0.12  0.04  0.01  0.01  0.00  0.00 -0.01 -0.03 -0.04 
##    15     8    12     3    17 
## -0.09 -0.15 -0.17 -0.26 -0.41 
## Bad chains:numeric(0)


## DDM parameters are correlated with one another
## esp. a-sv and a-v; v-sv
pairs(samples2)

# Sample posterior predictive to check fit
pp <- post.predict.ggdmc(samples2)


## PDF is good, CDF is poor
plot(pp, style="both")

save(pp, samples2, samples1, samples0, mdi1, dat1, p.vector, p.prior, m1,
     file="ggdmc_3_4.RData")


