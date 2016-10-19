##################  DMC Lesson 4: Multiple Subjects

### Lesson 4.7 Extra: Use C++ function to run hierarchical-Bayesian mdoel


rm(list=ls())
## Set up model and import/simulate data ------------------------------------
## To avoid name conflict, this lesson does not load DMC functions.
## See speed demo at the end of the script
setwd("/media/yslin/OWL/Documents/ggdmc/examples/")
library(ggdmc)
load("/media/yslin/OWL/Documents/ggdmc/examples/ggdmc_4_7_random.RData")

## Set up a DDM Model, rate effect of factor F
m1 <- model.dmc(
  p.map     = list(a="1",v="F",z="1",d="1",sz="1",sv="1",t0="1",st0="1"),
  match.map = list(M=list(s1="r1",s2="r2")),
  factors   = list(S=c("s1","s2"), F=c("f1","f2")),
  constants = c(st0=0,d=0),
  responses = c("r1","r2"),
  type      = "rd")

# Population distribution
pop.mean  <- c(a=1.6, v.f1=1.35, v.f2=0.75,  z=.5, sz=.3, sv=1,  t0=.3)
pop.scale <- c(a=.15, v.f1=.25,   v.f2=.25,  z=.1, sz=.1, sv=.3, t0=.05)
pop.prior <- prior.p.dmc(
  dists = rep("tnorm",7),
  p1    = pop.mean,
  p2    = pop.scale,
  lower = c(0,-5, -5, 0, 0, 0, 0),
  upper = c(5, 7,  7, 2, 2, 2, 2))


##  Check population distributions
plot_priors(pop.prior)

# Simulate a moderate size of data
dat1 <- h.simulate.dmc(m1, nsim=100, ns=30, p.prior=pop.prior)
mdi1 <- data.model.dmc(dat1, m1)


## Take a look at the first 10 subjects data
## The user should adjust canvas size, otherwise, the error of figure margins too large
## may come up
## The accuracies are around 60-80%
par(mfrow=c(6,5)) 
for (i in 1:length(mdi1))
{
  plot_cell_density(data.cell =mdi1[[i]][mdi1[[i]]$S=="s1", ], C="r1")
}
par(mfrow=c(1,1))



## Take a look at parameters
ps1 <- round( attr(dat1, "parameters"), 2); ps1
##       a v.f1 v.f2    z   sz   sv   t0
## 1  1.61 1.68 0.63 0.46 0.09 0.85 0.29
## 2  1.80 2.12 0.57 0.62 0.26 0.79 0.33
## 3  1.16 1.76 0.97 0.33 0.19 1.09 0.26
## 4  1.60 1.81 0.98 0.58 0.26 0.88 0.33
## 5  1.64 1.04 0.46 0.56 0.31 0.71 0.30
## 6  1.69 1.59 0.86 0.63 0.09 0.84 0.39
## 7  1.60 1.11 0.24 0.64 0.40 1.30 0.27
## 8  1.78 1.52 0.71 0.34 0.16 1.16 0.29
## 9  1.81 1.54 0.85 0.53 0.24 1.15 0.42
## 10 1.75 1.28 0.79 0.48 0.13 0.69 0.31
## 11 1.59 1.01 0.86 0.45 0.22 0.75 0.31
## 12 1.77 1.62 0.87 0.64 0.29 1.02 0.28
## 13 1.78 1.37 0.58 0.37 0.27 1.34 0.22
## 14 1.64 1.10 0.79 0.38 0.37 1.19 0.34
## 15 1.32 1.28 0.89 0.59 0.30 0.97 0.21
## 16 1.65 1.55 1.05 0.51 0.34 1.90 0.25
## 17 1.64 1.53 0.88 0.55 0.31 0.83 0.26
## 18 1.67 1.40 0.18 0.54 0.35 0.85 0.29
## 19 1.71 0.99 0.45 0.52 0.26 1.25 0.33
## 20 1.78 0.92 0.78 0.54 0.13 1.27 0.41
## 21 1.92 1.38 0.82 0.71 0.39 1.23 0.27
## 22 1.61 1.02 0.86 0.54 0.35 0.31 0.30
## 23 1.70 1.52 0.81 0.38 0.29 0.67 0.39
## 24 1.52 1.41 0.99 0.52 0.38 1.07 0.24
## 25 1.61 1.58 0.88 0.50 0.25 1.41 0.24
## 26 1.69 1.45 0.59 0.55 0.29 0.91 0.33
## 27 1.56 1.37 0.82 0.33 0.39 0.70 0.38
## 28 1.41 1.21 0.59 0.47 0.13 1.59 0.35
## 29 1.68 0.99 0.65 0.42 0.35 0.38 0.28
## 30 1.53 1.61 0.71 0.60 0.40 1.10 0.29

## FIT RANDOM EFFECTS ------------------------------------
## specify a data-level prior distributions
p.prior <- prior.p.dmc(
    dists = rep("tnorm", length(pop.mean)),
    p1    = pop.mean,
    p2    = pop.scale*5,
    lower = c(0,-5, -5, 0, 0, 0, 0),
    upper = c(5, 7,  7, 2, 2, 2, 2))

## specifiy a hyper-level prior distributions
mu.prior <- prior.p.dmc(
    dists = rep("tnorm", length(pop.mean)),
    p1    = pop.mean,
    p2    = pop.scale*5,
    lower = c(0,-5, -5, 0, 0, 0, 0),
    upper = c(5, 7,  7, 2, 2, 2, 2))

sigma.prior <- prior.p.dmc(
    dists = rep("beta", length(p.prior)),
    p1    = c(a=1, v.f1=1,v.f2 = 1, z=1, sz=1, sv=1, t0=1),
    p2    = c(1,1,1,1,1,1,1),
    upper = c(2,2,2,2,2,2,2))

pp.prior <- list(mu.prior, sigma.prior) # Make a hyper-prior list

plot_priors(p.prior)
plot_priors(mu.prior)
plot_priors(sigma.prior)

## Run model fits
hsamples0 <- h.samples.dmc(nmc=50, p.prior=p.prior, pp.prior=pp.prior,
                           data=mdi1, thin=5)
hsamples0 <- h.run.dmc(hsamples0, report=10, p.migrate=.05, h.p.migrate=.05)

## Looking well, as usual t0 and st0 have a stuck chain
plot(hsamples0, hyper=TRUE, density=TRUE)

## Add 250 iterations on top of the original one
hsamples1 <- h.samples.dmc(nmc=250, p.prior=p.prior, pp.prior=pp.prior,
                           samples=hsamples0, thin=5, add=TRUE)
hsamples1 <- h.run.dmc(hsamples1, report=50, p.migrate=.05, h.p.migrate=.05)


## There is one posterior chain still get stuck
plot(hsamples1, hyper=TRUE, pll.chain=TRUE, start=101)

## This takes almost 30 mins
hsamples2 <- h.samples.dmc(nmc=500, p.prior=p.prior, pp.prior=pp.prior,
                           samples=hsamples1, thin=5)
system.time(hsamples2 <- h.run.dmc(hsamples2, report=100, p.migrate=.05, h.p.migrate=.05))


## This almost take 1hrs
hsamples3 <- h.samples.dmc(nmc=1000, p.prior=p.prior, pp.prior=pp.prior,
                           samples=hsamples2, thin=5)
system.time(hsamples3 <- h.run.dmc(hsamples3, report=200))


## t0.h1 is almost converged
## t0.h2 still has one stuck chain
## sv.h1 and sv.h2 are as usual very difficult to recover
plot(hsamples3, hyper=TRUE)
plot(hsamples3, hyper=TRUE, pll.chain=TRUE)  ## hyper posterior has one stuck chain
plot(hsamples3, pll.chain=TRUE)

## This only problem is t0. 
gelman.diag.dmc(hsamples3, hyper=TRUE)
## Potential scale reduction factors:
## 
##         Point est. Upper C.I.
## a.h1          1.01       1.01
## v.f1.h1       1.01       1.01
## v.f2.h1       1.01       1.02
## z.h1          1.01       1.02
## sz.h1         1.01       1.02
## sv.h1         1.05       1.06
## t0.h1         3.49       7.44
## a.h2          1.01       1.02
## v.f1.h2       1.01       1.01
## v.f2.h2       1.02       1.03
## z.h2          1.01       1.02
## sz.h2         1.01       1.02
## sv.h2         1.04       1.07
## t0.h2         5.17       6.74
## 
## Multivariate psrf
## 
## 5.81


## Reasonable large effective sample sizes
effectiveSize.dmc(hsamples3, hyper=TRUE)
##    a.h1 v.f1.h1 v.f2.h1    z.h1   sz.h1   sv.h1   t0.h1    a.h2 v.f1.h2 v.f2.h2 
##    1495    1322    1434    1435     460    1265    1318    1333    1339    1165 
##    z.h2   sz.h2   sv.h2   t0.h2 
##    1379     771     798    1230

# Less so at individual subject level
es <- effectiveSize.dmc(hsamples1)

round(apply(data.frame(es),1,mean))
round(apply(data.frame(es),1,min))
##    a v.f1 v.f2    z   sz   sv   t0 
## 1068 1239 1446 1682 1010  930 1586 
##  627  628  751  853  738  403  557 


# Parameter recovery
summary(hsamples3, hyper.means=TRUE)
##           a      v.f1      v.f2         z        sz        sv        t0
## h1 1.638787 0.7619536 0.2265866 0.3275682 0.2715193 0.1066117 0.3758328
## h2 1.404375 0.5058083 0.8962963 0.1801540 0.2143202 0.1842407 0.1483826

hest1 <- summary(hsamples3)

assess.bias <- function(ps=NULL, hest=NULL)
{
  tru.h1   <- colMeans(ps) 
  est.h1   <- hest$quantiles[1:7, 3] 
  bias.h1  <- colMeans(ps) - hest$quantiles[1:7, 3]
  abias.h1 <- abs(colMeans(ps) - hest$quantiles[1:7, 3])
  location <- rbind(tru.h1, est.h1, bias.h1, abias.h1)
  print(apply(location, 2, round, 3))

  tru.h2   <- apply(ps, 2, sd)
  est.h2   <- hest$quantiles[8:14, 3]
  bias.h2  <- apply(ps, 2, sd) - hest$quantiles[8:14, 3]
  abias.h2 <- abs(apply(ps, 2, sd) - hest$quantiles[8:14, 3])
  scale <- rbind(tru.h2, est.h2, bias.h2, abias.h2)
  print(apply(scale, 2, round, 3))
}


## 1. This stuck chain at t0 does not affect overall estimation,
## because there are 21 chains and t0 is relatively easy to recover.
## 2. The scale parameter for sv is estiamted well. This
## may be due to the good data quality (error rate around 60-70)
assess.bias(ps1, hest1)
##              a   v.f1   v.f2     z    sz    sv     t0
## tru.h1   1.641  1.392  0.737 0.509 0.273 1.007  0.305
## est.h1   1.639  1.403  0.760 0.506 0.245 0.898  0.306
## bias.h1  0.002 -0.011 -0.023 0.004 0.028 0.109 -0.001
## abias.h1 0.002  0.011  0.023 0.004 0.028 0.109  0.001
##               a  v.f1 v.f2      z     sz     sv     t0
## tru.h2    0.151 0.286 0.21  0.100  0.095  0.337  0.055
## est.h2    0.177 0.265 0.21  0.104  0.177  0.342  0.062
## bias.h2  -0.026 0.021 0.00 -0.004 -0.083 -0.005 -0.007
## abias.h2  0.026 0.021 0.00  0.004  0.083  0.005  0.007


## save(hsamples3, hsamples2, hsamples1, hsamples0, m1, mdi1, dat1, pp.prior, ps1,
##       p.prior, pop.mean, pop.scale, pop.prior, file="ggdmc_4_7_random.RData")



## Speed Difference Demo ----------------------------------------------
rm(list=ls()) 
## Current working directory must be set to the top-level folder  
## containing the dmc and tutorial subfolders
setwd("/media/yslin/OWL/Documents/ggdmc/tests/")
source ("tutorial/file_utils.R")
load_model ("ddm.R")
## Use the same simulation data
load("/media/yslin/OWL/Documents/ggdmc/examples/ggdmc_4_7_random.RData")

## One CPU core competition
hsamples0 <- h.samples.dmc(nmc=50, p.prior=p.prior, pp.prior=pp.prior,
                           data=mdi1, thin=5)
system.time(hsamples0 <- h.run.dmc(hsamples0, report=10, p.migrate=.05, h.p.migrate=.05))


##     user   system  elapsed
## 1767.632    0.060 1767.907 

hsamples0.cpp <- ggdmc::h.samples.dmc(nmc=50, p.prior=p.prior, pp.prior=pp.prior,
                                      data=mdi1, thin=5)
system.time(hsamples0.cpp <- ggdmc::h.run.dmc(hsamples0.cpp, report=10, p.migrate=.05,
                                              h.p.migrate=.05 ))

dev.new(); plot.dmc(hsamples0, hyper=TRUE)
dev.new(); plot(hsamples0.cpp, hyper=TRUE)
##    user  system elapsed
## 145.360   0.016 145.382 
## It is about 12.16 times quicker in hierarchical-Bayesian multiple-subject model 

## Multiple cores competition with large data ---------------------------
## This data set simulates 15 participants, 1e4 trials in each conditon
## Unfortunately, if you get a lot of CPUs (12 in this case), you would be better
## of by using DMC's R codes, comparing to 1 CPU C++ codes (397 s vs 466 s).
## Fortunately, multi-code C++ still won (71s)
## 12 cores C++ is about 5.6 times and 6.6 times faster than R codes and 1-core
## C++ cores, respectively.
rm(list=ls())
setwd("/media/yslin/OWL/Documents/ggdmc/tests/")
source ("tutorial/file_utils.R")
load_model ("ddm.R")
load("/media/yslin/OWL/Documents/ggdmc/examples/ggdmc_4_7_random_speed.RData")

mdi11 <- data.model.dmc(raw.data11, model)
head(raw.data11)
## d <- dplyr::tbl_dt(raw.data11)
## d[, .N, .(s,S,F)]
length(unique(raw.data11$s))  ## [1] 15 subjects
## Source: local data table [60 x 4]
## 
##         s      S      F     N
##    (fctr) (fctr) (fctr) (int)
## 1       1     s1     f1 10000
## 2       1     s2     f1 10000
## 3       1     s1     f2 10000
## 4       1     s2     f2 10000
## 5       2     s1     f1 10000
## 6       2     s2     f1 10000
## 7       2     s1     f2 10000
## 8       2     s2     f2 10000
## 9       3     s1     f1 10000
## 10      3     s2     f1 10000
## ..    ...    ...    ...   ...


##      user  system elapsed
## 2398.540   90.132  396.867 
hsamples1 <- h.samples.dmc(nmc=50, p.prior=p.prior, pp.prior=pp.prior,
                           data=mdi11, thin=1)
system.time(hsamples1 <- h.run.dmc(hsamples1, report=10, cores=12,
                                   p.migrate=.05, h.p.migrate=.05))


## Note that because Open MPI allows you to control your PC.
## It will be only a fair comparsion, if you leave you PC alone when running
## the following codes.
##     user  system elapsed
## 787.972   1.860  70.662 
hsamples1.cpp <- ggdmc::h.samples.dmc(nmc=50, p.prior=p.prior, pp.prior=pp.prior,
                           data=mdi11, thin=1)
system.time(hsamples1.cpp <- ggdmc::h.run.dmc(hsamples1.cpp, report=10, cores=2,
                                    p.migrate=.05, h.p.migrate=.05 ))

##    user  system elapsed
## 437.364   8.968 446.355 
hsamples2.cpp <- ggdmc::h.samples.dmc(nmc=50, p.prior=p.prior, pp.prior=pp.prior,
                           data=mdi11, thin=1)
system.time(hsamples2.cpp <- ggdmc::h.run.dmc(hsamples2.cpp, report=10, cores=1,
                                    p.migrate=.05, h.p.migrate=.05 ))

