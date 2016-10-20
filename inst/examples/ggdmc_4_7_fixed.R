##################  DMC Lesson 4: Multiple Subjects

### Lesson 4.7 Extra: Use C++ function to run hierarchical-Bayesian mdoel


rm(list=ls())
## Set up model and import/simulate data ------------------------------------
## To avoid name conflict, this lesson does not load DMC functions.
## See speed demo at the end of the script
library(ggdmc)
load("ggdmc_4_7_fixed.RData")

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
dat1 <- h.simulate.dmc(m1, nsim=50, ns=12, p.prior=pop.prior)
mdi1 <- data.model.dmc(dat1, m1)


## Take a look at the first 10 subjects data
par(mfrow=c(3,4)) # Row 1 = subjects 1..5, row2 = 6..10
for (i in 1:length(mdi1))
{
  plot_cell_density(data.cell =mdi1[[i]][mdi1[[i]]$S=="s1", ], C="r1")
}
par(mfrow=c(1,1))



## Take a look at parameters
ps1 <- round( attr(dat1, "parameters"), 2); ps1
##       a v.f1 v.f2    z   sz   sv   t0
## 1  1.61 1.35 1.04 0.52 0.34 1.02 0.36
## 2  1.32 1.28 0.58 0.49 0.24 1.34 0.25
## 3  1.60 1.24 0.92 0.50 0.28 1.03 0.25
## 4  1.49 1.12 0.62 0.55 0.35 0.63 0.27
## 5  1.55 1.33 0.81 0.48 0.40 0.89 0.28
## 6  1.44 1.49 0.51 0.48 0.17 0.73 0.37
## 7  1.98 1.41 0.83 0.61 0.33 1.02 0.33
## 8  1.59 1.74 1.04 0.58 0.24 0.98 0.25
## 9  1.58 1.41 0.56 0.44 0.21 1.26 0.28
## 10 1.39 1.54 0.87 0.49 0.29 0.84 0.31
## 11 1.49 1.53 0.88 0.38 0.35 1.42 0.27
## 12 1.49 0.98 1.01 0.44 0.14 1.14 0.29


## FIT FIXED EFFECTS --------------------------------------
## specify a broader prior than the true population distribution
p.prior <- prior.p.dmc(
    dists= rep("tnorm", length(pop.mean)),
    p1=pop.mean,
    p2=pop.scale*5,
    lower=c(0,-5, -5, 0, 0, 0, 0),
    upper=c(5, 7,  7, 2, 2, 2, 2))


samples0 <- h.samples.dmc(nmc=500, p.prior=p.prior, data=mdi1, thin=2)


## Use 12 CPU cores to fit 12 participants in parllel
samples0 <- h.run.dmc(samples0, cores=12)

plot(samples0) ## Default plot 1st subject
plot(samples0, subject=2)
plot(samples0, subject=3)
plot(samples0, subject=4)
plot(samples0, subject=5)
plot(samples0, subject=6)
plot(samples0, subject=7)
plot(samples0, subject=8)
plot(samples0, subject=9)
plot(samples0, subject=10)
plot(samples0, subject=11)
plot(samples0, subject=12)  ## Converge around 500th


## Get 500 new iterations
samples1 <- h.run.dmc(h.samples.dmc(nmc=500, p.prior=p.prior, samples=samples0, thin=2),
                      cores=12)

## Add 500 more 
samples2 <- h.run.dmc(h.samples.dmc(nmc=500, p.prior=p.prior, samples=samples1,
                                    thin=2, add=TRUE), cores=12)


plot(samples2, pll.chain=TRUE)
plot(samples2, subject=2, density=TRUE)
plot(samples2, subject=3, density=TRUE)

## For more arguments about plot.dmc.list 
?plot.dmc.list


# All 12 participants show chains  converged.
gelman.diag.dmc(samples2)
##   11    1    6    2    4    9    3    8   10    5   12    7 
## 1.01 1.02 1.02 1.02 1.02 1.02 1.02 1.02 1.02 1.02 1.02 1.03 
## Mean
## [1] 1.02

## Effective sample sizes are reasonable large
es <- effectiveSize.dmc(samples2)
round(apply(do.call(rbind,es),2,median))
apply(do.call(rbind,es),2,min)
##    a v.f1 v.f2    z   sz   sv   t0 
## 1356 1477 1461 1536 1394 1447 1444 
## 1206 1311 1333 1353 1148 1223 1184


## Check against true values: fairly good parameter recovery.
est <- summary(samples2)

## More information about summary function
?summary.dmc.list


est.mean <- t(data.frame(lapply(est,function(x){
  x$statistics[,"Mean"]})))[,dimnames(ps)[[2]]]

## The data quality is good, so sz and sv recover very
## well. v.f2 is poorly recovery.

round(apply(est.mean,2,mean),3)
round(apply(ps,2,mean), 3)
round(apply(est.mean,2,mean)-apply(ps,2,mean),3)
##     a   v.f1  v.f2       z     sz     sv    t0 
## 1.598  1.432  0.912  0.492  0.334  1.053  0.300 
## 1.544  1.368  0.806  0.497  0.278  1.025  0.292 
## 0.054  0.064  0.106 -0.005  0.055  0.028  0.007


## Scale parameters recover well, too.
round(apply(est.mean,2,sd),3)
round(apply(ps, 2, sd), 3)
round(apply(est.mean,2,sd)-apply(ps,2,sd),3)
##     a   v.f1   v.f2      z     sz     sv    t0 
## 0.188  0.206  0.288  0.055  0.077  0.205  0.043 
## 0.164  0.202  0.193  0.063  0.080  0.238  0.042 
## 0.024  0.004  0.096 -0.009 -0.003 -0.032  0.001


## Good fits? Error responses are fit poorly
pps1  <- post.predict.ggdmc(samples2[[1]])
pps2  <- post.predict.ggdmc(samples2[[2]])
pps3  <- post.predict.ggdmc(samples2[[3]])
pps4  <- post.predict.ggdmc(samples2[[4]])
pps5  <- post.predict.ggdmc(samples2[[5]])
pps6  <- post.predict.ggdmc(samples2[[6]])
pps7  <- post.predict.ggdmc(samples2[[7]])
pps8  <- post.predict.ggdmc(samples2[[8]])
pps9  <- post.predict.ggdmc(samples2[[9]])
pps10 <- post.predict.ggdmc(samples2[[10]])
## pps11 <- post.predict.ggdmc(samples2[[11]]) ## Too fews errors
pps12 <- post.predict.ggdmc(samples2[[12]])

plot(pps1,  gpvar1="F", style="both")
plot(pps2,  gpvar1="F", style="both")
plot(pps3,  gpvar1="F", style="both")
plot(pps4,  gpvar1="F", style="both")
plot(pps5,  gpvar1="F", style="both")
plot(pps6,  gpvar1="F", style="both")
plot(pps7,  gpvar1="F", style="both")
plot(pps8,  gpvar1="F", style="both")
plot(pps9,  gpvar1="F", style="both")
plot(pps10, gpvar1="F", style="both")
## plot(pps11, gpvar1="F", style="both")
plot(pps12, gpvar1="F", style="both")


## For more details
?plot.pp.ggdmc

## save(pps1, pps2, pps3, pps4, pps5, pps6, pps7, pps8, pps9, pps10, pps12
##      samples2, samples1, samples0, m1, mdi1, dat1,
##      p.prior, pop.mean, 
##      pop.scale, pop.prior, file="ggdmc_4_7_fixed.RData")

