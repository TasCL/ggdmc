# System functions for the DMC (Dynamic Models of Choice)
#    Functions specific to extract and analyze samples
#    Usually user does not need to edit

#' Find Stuck Chains
#'
#' This function picks out the chains that get stuck.
#'
#' @param samples a DMC object/sample
#' @param hyper whether to work on hyper parameters, namely phi
#' @param cut cutoff point. Default 10
#' @param start starting from which iteration. Default 1
#' @param end calculating until which iteration. Default is the last one,
#' extracted from the input DMC sample.
#' @param verbose whether print out debugging information
#' @param digits print how many digits. Default 2
#' @keywords pick.stuck
#' @return index of stuck chains (deviaton of ll from median < -cut)
#' @export
pick.stuck.dmc <- function(samples, hyper=FALSE, cut=10, start=1,end=NA,
                           verbose=FALSE, digits=2)
{
  if (hyper) { # MG: TODO: Haven't split up hyper$pll yet
    hyper <- attr(samples,"hyper")
    if ( is.na(end) ) end <- dim(hyper$pll)[1]
    if ( end <= start )
      stop("End must be greater than start")
    mean.ll <- apply(hyper$pll[start:end,],2,mean)
    names(mean.ll) <- 1:length(mean.ll)
  } else {
    if ( is.na(end) ) end <- samples$nmc
    if ( end <= start ) stop("End must be greater than start")

    mean.ll <- apply(samples$log_likelihoods[start:end,] +
                     samples$summed_log_prior[start:end,],2,mean)
    names(mean.ll) <- 1:length(mean.ll)
  }
  dev <- -(sort(mean.ll)-median(mean.ll))
  bad <- as.numeric(names(dev)[dev>cut])
  if (verbose) {
    cat("Deviation of mean chain log-likelihood from median of means\n")
    print(round(dev,digits))
    cat("Bad chains:")
  }
  bad
}


pick.stucks.dmc <- function(samples,cut=10,start=1,end=NA,
                           verbose=FALSE,digits=2)
# List named by subject with chain numbers that are stuck.
{
  tmp <- lapply(lapply(samples,
    pick.stuck.dmc,start=start,cut=cut,digits=digits),
    function(x){if (length(x)>0) x})
  empty <- unlist(lapply(tmp,is.null))
  if (any(empty)) tmp[-c(1:length(tmp))[empty]] else tmp
}


unstick.dmc <- function(samples, bad)
  # remove stuck chains from a samples object
{
  if (length(bad)>0) {
    if (!all(bad %in% 1:samples$n.chains))
      stop(paste("Index of bad chains must be in 1 :",samples$n.chains))
    samples$theta <- samples$theta[-bad,,]
    samples$summed_log_prior <- samples$summed_log_prior[,-bad]
    samples$log_likelihoods  <- samples$log_likelihoods[,-bad]
    samples$n.chains <- samples$n.chains-length(bad)
  }
  samples
}


#' Convert Theta to a mcmc List
#'
#' \code{theta.as.mcmc.list} extract theta from a DMC sample and convert it to
#' a \pkg{coda} mcmc.list object.
#'
#' @param samples a DMC sample
#' @param start start iteration
#' @param end end iteraton
#' @importFrom coda mcmc mcmc.list
#' @export
theta.as.mcmc.list <- function(samples, start=1, end=NA)
{
  if ( is.na(end) ) end <- dim(samples$theta)[3]
  n.chains <- dim(samples$theta)[1]
  lst <- vector(mode="list",length=n.chains)
  for (i in 1:n.chains)
    lst[[i]] <- mcmc(t(samples$theta[i,,start:end]))
  mcmc.list(lst)
}

#' Convert Phi to a Theta Vector
#'
#' \code{phi.as.mcmc.list} extract phi from a sample's hyper attribute and
#' convert it to theta vector to a mcmc.list \pkg{coda} object
#'
#' @param hyper a hyper attribute extracted from a hierarchical sample
#' @param start start iteration
#' @param end end iteration
#' @importFrom coda mcmc mcmc.list
#' @export
phi.as.mcmc.list <- function(hyper, start=1, end=NA)
{
  # Parameters that are not constants
  ok1 <- lapply(hyper$pp.prior, function(x) {
    lapply(x,function(y){attr(y,"dist")!="constant"})}
  )
  ok2 <- paste(names(ok1[[2]])[unlist(ok1[[2]])], "h2", sep=".")
  ok1 <- paste(names(ok1[[1]])[unlist(ok1[[1]])], "h1", sep=".")
  if ( is.na(end) ) end <- dim(hyper$phi[[1]])[3]
  n.chains <- dim(hyper$h_log_likelihoods)[2]
  lst <- vector(mode="list",length=n.chains)
  for (i in 1:n.chains) {
    tmp1 <- t(hyper$phi[[1]][i,,start:end])
    dimnames(tmp1)[[2]] <- paste(dimnames(tmp1)[[2]],"h1",sep=".")
    tmp1 <- tmp1[,ok1]
    tmp2 <- t(hyper$phi[[2]][i,,start:end])
    dimnames(tmp2)[[2]] <- paste(dimnames(tmp2)[[2]],"h2",sep=".")
    tmp2 <- tmp2[,ok2]
    # Remove cases with !has.sigma
    tmp2 <- tmp2[,!apply(tmp2,2,function(x){all(is.na(x))})]
    lst[[i]] <- mcmc(cbind(tmp1,tmp2))
  }
  mcmc.list(lst)
}

#' Gelman and Rubin Convergence Diagnostic
#'
#' \code{gelman.diag.dmc}  calls \pkg{coda} gelman.diag to get R hats for one
#' or a list of subjects. It can calculate at the data or hyper level.
#'
#' @param x a DMC sample
#' @param hyper a switch to extract hyper attribute and calculate it
#' @param digits print out how many digits
#' @param start start iteration
#' @param autoburnin turn on auto burnin
#' @param transform turn on transform
#' @param end end iteraton
#' @param ... arguments passing to \code{coda} gelman.diag.
#' @importFrom coda gelman.diag
#' @export
#' @examples
#' m1 <- model.dmc(
#'      p.map     = list(a="1",v="F",z="1",d="1",sz="1",sv="1",t0="1",st0="1"),
#'      match.map = list(M=list(s1="r1",s2="r2")),
#'      factors   = list(S=c("s1","s2"),F=c("f1","f2")),
#'      constants = c(st0=0,d=0),
#'      responses = c("r1","r2"),
#'      type      = "rd")
#'
#' pop.mean  <- c(a=1.15, v.f1=1.25, v.f2=1.85, z=.55,  sz=.15, sv=.32, t0=.25)
#' pop.scale <- c(a=.10,  v.f1=.8,   v.f2=.5,   z=0.1,  sz=.05, sv=.05, t0=.05)
#' pop.prior <- prior.p.dmc(
#'   dists = rep("tnorm", length(pop.mean)),
#'   p1    = pop.mean,
#'   p2    = pop.scale,
#'   lower = c(0,-5, -5, 0, 0,   0, 0),
#'   upper = c(5, 7,  7, 1, 0.5, 2, 2))
#'
#' dat  <- h.simulate.dmc(m1, nsim=30, ns=8, p.prior=pop.prior)
#' mdi1 <- data.model.dmc(dat, m1)
#' ps   <- attr(dat,  "parameters")
#' ### FIT RANDOM EFFECTS
#' p.prior <- prior.p.dmc(
#'   dists = c("tnorm","tnorm","tnorm","tnorm","tnorm", "tnorm", "tnorm"),
#'   p1=pop.mean,
#'   p2=pop.scale*5,
#'   lower=c(0,-5, -5, 0, 0, 0, 0),
#'   upper=c(5, 7,  7, 2, 2, 2, 2))
#'
#' mu.prior <- prior.p.dmc(
#'   dists = c("tnorm","tnorm","tnorm","tnorm","tnorm", "tnorm", "tnorm"),
#'   p1=pop.mean,
#'   p2=pop.scale*5,
#'   lower=c(0,-5, -5, 0, 0, 0, 0),
#'   upper=c(5, 7,  7, 2, 2, 2, 2))
#'
#' sigma.prior <- prior.p.dmc(
#'   dists = rep("beta", length(p.prior)),
#'   p1=c(a=1, v.f1=1,v.f2 = 1, z=1, sz=1, sv=1, t0=1),p2=c(1,1,1,1,1,1,1),
#'   upper=c(2,2,2,2,2, 2, 2))
#'
#' pp.prior <- list(mu.prior, sigma.prior)
#'
#' hsamples0 <- h.samples.dmc(nmc=200, p.prior=p.prior, pp.prior=pp.prior,
#'   data=mdi1, thin=1)
#' hsamples0 <- h.run.dmc(hsamples0, p.migrate=.05, h.p.migrate=.05)
#' gelman.diag.dmc(hsamples0, hyper=TRUE)
#' ## Potential scale reduction factors:
#' ##
#' ## Point est. Upper C.I.
#' ## a.h1          1.08       1.14
#' ## v.f1.h1       1.31       1.53
#' ## v.f2.h1       1.18       1.30
#' ## z.h1          1.16       1.29
#' ## sz.h1         1.05       1.09
#' ## sv.h1         1.13       1.22
#' ## t0.h1         1.11       1.19
#' ## a.h2          1.10       1.18
#' ## v.f1.h2       1.06       1.10
#' ## v.f2.h2       1.07       1.12
#' ## z.h2          1.10       1.17
#' ## sz.h2         1.08       1.14
#' ## sv.h2         1.12       1.20
#' ## t0.h2         1.10       1.18
#' ## Multivariate psrf
#'
#' gelman.diag.dmc(hsamples0)
#' ##    2    8    4    3    1    5    6    7
#' ## 1.21 1.21 1.27 1.36 1.37 1.55 1.57 1.77
#' ## Mean
#' ## [1] 1.41
gelman.diag.dmc <- function(x, hyper=FALSE, digits=2, start=1,
                            autoburnin=FALSE, transform=TRUE, end=NA, ...)
{
  if (hyper) {
    hyper0 <- attr(x, "hyper")
    if (is.na(end)) end <- hyper0$nmc
    gelman.diag(phi.as.mcmc.list(hyper0, start=start,end=end),
                autoburnin=autoburnin, transform=transform,...)
  } else {
    if ( !is.null(x$theta) ) {
      if (is.na(end)) end <- x$nmc
      gelman.diag(theta.as.mcmc.list(x, start=start,end=end),
                  autoburnin=autoburnin,transform=transform,...)
    } else {
      out <- lapply(x, function(xx){
        if (is.na(end)) end <- xx$nmc
        gelman.diag(theta.as.mcmc.list(xx, start=start),
                    autoburnin=autoburnin,transform=transform)
      })
      tmp <- unlist(lapply(out,function(xx){xx$mpsrf}))
      print(round(sort(tmp), digits))
      cat("Mean\n")
      print(round(mean(tmp), digits))
      invisible(out)
    }
  }
}

#' Effective Sample Size for Estimating the Mean
#'
#' \code{effectiveSize.dmc} calls \pkg{coda} effectiveSize to effective size
#' for either single or multiple subjects. It can calculate at the data or
#' hyper level, too.
#'
#' @param x a DMC sample
#' @param hyper a switch to extract hyper attribute and calculate it
#' @param digits print out how many digits
#' @param start start iteration
#' @param end end iteraton
#' @importFrom coda effectiveSize
#' @export
#' @examples
#' m1 <- model.dmc(
#'      p.map     = list(a="1",v="F",z="1",d="1",sz="1",sv="1",t0="1",st0="1"),
#'      match.map = list(M=list(s1="r1",s2="r2")),
#'      factors   = list(S=c("s1","s2"),F=c("f1","f2")),
#'      constants = c(st0=0,d=0),
#'      responses = c("r1","r2"),
#'      type      = "rd")
#'
#' pop.mean  <- c(a=1.15, v.f1=1.25, v.f2=1.85, z=.55,  sz=.15, sv=.32, t0=.25)
#' pop.scale <- c(a=.10,  v.f1=.8,   v.f2=.5,   z=0.1,  sz=.05, sv=.05, t0=.05)
#' pop.prior <- prior.p.dmc(
#'   dists = rep("tnorm", length(pop.mean)),
#'   p1    = pop.mean,
#'   p2    = pop.scale,
#'   lower = c(0,-5, -5, 0, 0,   0, 0),
#'   upper = c(5, 7,  7, 1, 0.5, 2, 2))
#'
#' dat  <- h.simulate.dmc(m1, nsim=30, ns=8, p.prior=pop.prior)
#' mdi1 <- data.model.dmc(dat, m1)
#' ps   <- attr(dat,  "parameters")
#' ### FIT RANDOM EFFECTS
#' p.prior <- prior.p.dmc(
#'   dists = c("tnorm","tnorm","tnorm","tnorm","tnorm", "tnorm", "tnorm"),
#'   p1=pop.mean,
#'   p2=pop.scale*5,
#'   lower=c(0,-5, -5, 0, 0, 0, 0),
#'   upper=c(5, 7,  7, 2, 2, 2, 2))
#'
#' mu.prior <- prior.p.dmc(
#'   dists = c("tnorm","tnorm","tnorm","tnorm","tnorm", "tnorm", "tnorm"),
#'   p1=pop.mean,
#'   p2=pop.scale*5,
#'   lower=c(0,-5, -5, 0, 0, 0, 0),
#'   upper=c(5, 7,  7, 2, 2, 2, 2))
#'
#' sigma.prior <- prior.p.dmc(
#'   dists = rep("beta", length(p.prior)),
#'   p1=c(a=1, v.f1=1,v.f2 = 1, z=1, sz=1, sv=1, t0=1),p2=c(1,1,1,1,1,1,1),
#'   upper=c(2,2,2,2,2, 2, 2))
#'
#' pp.prior <- list(mu.prior, sigma.prior)
#'
#' hsamples0 <- h.samples.dmc(nmc=200, p.prior=p.prior, pp.prior=pp.prior,
#'   data=mdi1, thin=1)
#' hsamples0 <- h.run.dmc(hsamples0, p.migrate=.05, h.p.migrate=.05)
#' es <- effectiveSize.dmc(hsamples0)
#' ## $`1`
#' ##   a v.f1 v.f2    z   sz   sv   t0
#' ## 253  193   93  292  134  184  370
#' ##
#' ## $`2`
#' ##   a v.f1 v.f2    z   sz   sv   t0
#' ## 169   78  142  254  193  187  343
#' ##
#' ## $`3`
#' ##   a v.f1 v.f2    z   sz   sv   t0
#' ## 118   56  105  121  157  129  210
#' ##
#' ## $`4`
#' ##   a v.f1 v.f2    z   sz   sv   t0
#' ## 179   77   95  236  121  175  124
#' ##
#' ## $`5`
#' ##   a v.f1 v.f2    z   sz   sv   t0
#' ## 199   49   64  205  153  124  339
#' ##
#' ## $`6`
#' ##   a v.f1 v.f2    z   sz   sv   t0
#' ## 220   65  103  208  128  123  231
#' ##
#' ## $`7`
#' ##   a v.f1 v.f2    z   sz   sv   t0
#' ## 234   58  108  219  191  190  203
#' ##
#' ## $`8`
#' ##   a v.f1 v.f2    z   sz   sv   t0
#' ## 390   81   90  182  249  191  257
#'
#' hes <- effectiveSize.dmc(hsamples0, hyper=TRUE)
#' ##   a.h1  v.f1.h1 v.f2.h1    z.h1   sz.h1   sv.h1   t0.h1    a.h2 v.f1.h2
#' ##     130     127     134     154      81     108     172     126     100
#' ## v.f2.h2    z.h2   sz.h2   sv.h2   t0.h2
#' ##     118     100     153      92     123
effectiveSize.dmc <- function(x, hyper=FALSE, digits=0,start=1,end=NA)
{
  if (hyper) {
    hyper0 <- attr(x, "hyper")
    if (is.na(end)) end <- hyper0$nmc
    round(effectiveSize(phi.as.mcmc.list(attr(x, "hyper"), start=start,
      end=end)), digits)
  } else {
      if (!is.null(x$theta)) {
        if (is.na(end)) end <- x$nmc
        round(effectiveSize(theta.as.mcmc.list(x, start=start,end=end)),digits)
      } else
          lapply(x, function(xx){
            if (is.na(end)) end <- xx$nmc
            round(effectiveSize(theta.as.mcmc.list(xx, start=start, end=end)),
              digits)
          })
  }
}

#' Summarise a DMC Sample with One Participant
#'
#' Call coda package to summarise the model parameters in a DMC samples
#'
#' @param object a model samples
#' @param start summarise from which MCMC iteration. Default uses the first
#' iteration.
#' @param end summarise to the end of MCMC iteration. For example, set
#' \code{start=101} and \code{end=1000}, instructs the function to calculate
#' from 101 to 1000 iteration. Default uses the last iteration.
#' @param ... other arguments
#' @keywords summary
#' @seealso \code{\link{summary.dmc.list}, \link{summary.hyper}}
#' @export
#' @examples
#' m1 <- model.dmc(
#' p.map=list(a="1",v="1",z="1",d="1",sz="1",sv="1", t0="1",st0="1"),
#' constants=c(st0=0,d=0),
#' match.map=list(M=list(s1="r1",s2="r2")),
#' factors=list(S=c("s1","s2")),
#' responses=c("r1","r2"),
#' type="rd")
#'
#' p.prior <- prior.p.dmc(
#'   dists = rep("tnorm", 6),
#'   p1=c(a=2,   v=2.5, z=0.5, sz=0.3, sv=1,  t0=0.3),
#'   p2=c(a=0.5, v=.5,  z=0.1, sz=0.1, sv=.3, t0=0.05),
#'   lower=c(0,-5, 0, 0, 0, 0),
#'   upper=c(5, 7, 2, 2, 2, 2))
#'
#' pVec <- c(a=1,v=1, z=0.5, sz=0.25, sv=0.2,t0=.15)
#'
#' dat1 <- simulate(m1, nsim=1e2, p.vector=pVec)
#' mdi1 <- data.model.dmc(dat1, m1)
#'
#' samples0 <- samples.dmc(nmc=500, p.prior=p.prior, data=mdi1)
#' samples0 <- run.dmc(samples0, p.migrate=.05)
#' gelman.diag.dmc(samples0)
#' class(samples0)
#' ## [1] "dmc"
#' summary(samples0)
summary.dmc <- function(object, start=1, end=NA, ... )
{
  if (is.na(end)) end <- object$nmc
  summary(theta.as.mcmc.list(object, start=start, end=end))
}

#' Summarise a DMC Sample with Multiple Participants
#'
#' Call coda package to summarise the model parameters in a DMC samples with
#' multiple participants
#'
#' @param object a model samples
#' @param digits how many digits to print
#' @param start summarise from which MCMC iteration. Default uses the first
#' iteration.
#' @param end summarise to the end of MCMC iteration. For example, set
#' \code{start=101} and \code{end=1000}, instructs the function to calculate
#' from 101 to 1000 iteration.  Default uses the last iteration.
#' @param ... other aruguments
#' @keywords summary
#' @seealso \code{\link{summary.dmc}, \link{summary.hyper}}
#' @export
#' @examples
#'
#' m1 <- model.dmc(
#'        p.map     = list(a="1",v="F",z="1",d="1",sz="1",sv="1",t0="1",st0="1"),
#'        match.map = list(M=list(s1="r1",s2="r2")),
#'        factors   = list(S=c("s1","s2"),F=c("f1","f2")),
#'        constants = c(st0=0,d=0),
#'        responses = c("r1","r2"),
#'        type      = "rd")
#'
#' pop.mean  <- c(a=1.15, v.f1=1.25, v.f2=1.85, z=0.55, sz=0.15, sv=0.32, t0=0.25)
#' pop.scale <-c(a=0.10, v.f1=.8,   v.f2=.5,   z=0.1,  sz=0.05, sv=0.05, t0=0.05)
#' pop.prior <- prior.p.dmc(
#'   dists = rep("tnorm", length(pop.mean)),
#'   p1    = pop.mean,
#'   p2    = pop.scale,
#'   lower = c(0,-5, -5, 0, 0,   0, 0),
#'   upper = c(5, 7,  7, 1, 0.5, 2, 2))
#'
#' dat  <- h.simulate.dmc(m1, nsim=30, ns=8, p.prior=pop.prior)
#' mdi1 <- data.model.dmc(dat, m1)
#' ps   <- attr(dat,  "parameters")
#'
#' p.prior <- prior.p.dmc(
#'   dists= rep("tnorm", length(pop.mean)),
#'   p1=pop.mean,
#'   p2=pop.scale*5,
#'   lower=c(0,-5, -5, 0, 0, 0, 0),
#'   upper=c(5, 7,  7, 2, 2, 2, 2))
#'
#' samples0 <- h.samples.dmc(nmc=100, p.prior=p.prior, data=mdi1, thin=1)
#' samples0 <- h.run.dmc(samples0, report=20)
#' class(samples0)
#' ## [1] "dmc.list"
#' gelman.diag.dmc(samples0)
#' ##    8    7    1    3    6    4    5    2
#' ## 1.68 1.74 1.86 1.92 1.96 1.98 2.02 2.43
#' ## Mean
#' ## [1] 1.95
#' summary(samples0)
#' ##       a v.f1 v.f2    z   sz   sv   t0
#' ## 1    1.59 1.12 1.94 0.67 0.33 0.67 0.23
#' ## 2    1.46 2.62 1.65 0.57 0.32 0.69 0.28
#' ## 3    1.38 1.02 0.92 0.56 0.35 0.64 0.22
#' ## 4    1.43 1.24 1.62 0.57 0.38 0.70 0.28
#' ## 5    1.35 1.68 1.24 0.48 0.31 0.54 0.24
#' ## 6    1.47 1.07 2.71 0.55 0.34 0.57 0.19
#' ## 7    1.40 1.13 2.24 0.52 0.43 0.67 0.23
#' ## 8    1.68 1.79 2.05 0.65 0.38 0.68 0.16
#' ## Mean 1.47 1.46 1.79 0.57 0.35 0.64 0.23
#'
#' summary.samples0 <- summary(samples0)
summary.dmc.list <- function(object, digits=2, start=1, end=NA, ...)
{
  out <- lapply(object, function(x) {
    if (is.na(end)) end <- x$nmc
    summary(theta.as.mcmc.list(x, start=start, end=end))
  })
  tmp <- t(data.frame(lapply(out,function(x){x[[1]][,1]})))
  tmp <- rbind(tmp, apply(tmp,2,mean))
  row.names(tmp) <- c(names(object), "Mean")
  print(round(tmp, digits))
  invisible(out)
}

#' Summarise a DMC Sample with Multiple Participant at the Hyper-level
#'
#' Call coda package to summarise the model parameters in a DMC samples with
#' multiple participants at the hyper level.
#'
#' @param object a model samples
#' @param start summarise from which MCMC iteration.
#' @param end summarise to the end of MCMC iteration. For example, set
#' \code{start=101} and \code{end=1000}, instructs the function to calculate
#' from 101 to 1000 iteration.
#' @param hyper.means default as FALSE
#' @param ... other arguments
#' @keywords summary
#' @seealso \code{\link{summary.dmc}, \link{summary.dmc.list}}
#' @export
#' @examples
#' m1 <- model.dmc(
#'      p.map     = list(a="1",v="F",z="1",d="1",sz="1",sv="1",t0="1",st0="1"),
#'      match.map = list(M=list(s1="r1",s2="r2")),
#'      factors   = list(S=c("s1","s2"),F=c("f1","f2")),
#'      constants = c(st0=0,d=0),
#'      responses = c("r1","r2"),
#'      type      = "rd")
#'
#' pop.mean  <- c(a=1.15, v.f1=1.25, v.f2=1.85, z=.55,  sz=.15, sv=.32, t0=.25)
#' pop.scale <- c(a=.10,  v.f1=.8,   v.f2=.5,   z=0.1,  sz=.05, sv=.05, t0=.05)
#' pop.prior <- prior.p.dmc(
#'   dists = rep("tnorm", length(pop.mean)),
#'   p1    = pop.mean,
#'   p2    = pop.scale,
#'   lower = c(0,-5, -5, 0, 0,   0, 0),
#'   upper = c(5, 7,  7, 1, 0.5, 2, 2))
#'
#' dat  <- h.simulate.dmc(m1, nsim=30, ns=8, p.prior=pop.prior)
#' mdi1 <- data.model.dmc(dat, m1)
#' ps   <- attr(dat,  "parameters")
#' ### FIT RANDOM EFFECTS
#' p.prior <- prior.p.dmc(
#'   dists = c("tnorm","tnorm","tnorm","tnorm","tnorm", "tnorm", "tnorm"),
#'   p1=pop.mean,
#'   p2=pop.scale*5,
#'   lower=c(0,-5, -5, 0, 0, 0, 0),
#'   upper=c(5, 7,  7, 2, 2, 2, 2))
#'
#' mu.prior <- prior.p.dmc(
#'   dists = c("tnorm","tnorm","tnorm","tnorm","tnorm", "tnorm", "tnorm"),
#'   p1=pop.mean,
#'   p2=pop.scale*5,
#'   lower=c(0,-5, -5, 0, 0, 0, 0),
#'   upper=c(5, 7,  7, 2, 2, 2, 2))
#'
#' sigma.prior <- prior.p.dmc(
#'   dists = rep("beta", length(p.prior)),
#'   p1=c(a=1, v.f1=1,v.f2 = 1, z=1, sz=1, sv=1, t0=1),p2=c(1,1,1,1,1,1,1),
#'   upper=c(2,2,2,2,2, 2, 2))
#'
#' pp.prior <- list(mu.prior, sigma.prior)
#'
#' hsamples0 <- h.samples.dmc(nmc=200, p.prior=p.prior, pp.prior=pp.prior,
#'   data=mdi1, thin=1)
#' hsamples0 <- h.run.dmc(hsamples0, p.migrate=.05, h.p.migrate=.05)
#'
#' class(hsamples0)
#' ## [1] "hyper"
#'
#' summary(hsamples0)
#' ## Iterations = 1:200
#' ## Thinning interval = 1
#' ## Number of chains = 21
#' ## Sample size per chain = 200
#' ##
#' ## 1. Empirical mean and standard deviation for each variable,
#' ## plus standard error of the mean:
#' ##
#' ##           Mean     SD Naive SE Time-series SE
#' ## a.h1    1.1728 0.4781 0.007377        0.04151
#' ## v.f1.h1 2.2737 1.7060 0.026324        0.15071
#' ## v.f2.h1 1.9654 1.0858 0.016754        0.09888
#' ## z.h1    0.4848 0.3275 0.005054        0.03275
#' ## sz.h1   0.3164 0.1952 0.003013        0.02378
#' ## sv.h1   0.4785 0.1814 0.002798        0.01769
#' ## t0.h1   0.2534 0.1671 0.002578        0.01489
#' ## a.h2    0.7267 0.4055 0.006258        0.04262
#' ## v.f1.h2 1.0990 0.3420 0.005277        0.03714
#' ## v.f2.h2 0.8329 0.2819 0.004350        0.02816
#' ## z.h2    0.7673 0.5092 0.007857        0.06046
#' ## sz.h2   0.4201 0.1735 0.002677        0.01727
#' ## sv.h2   1.0427 0.3516 0.005426        0.03977
#' ## t0.h2   0.2912 0.1921 0.002964        0.02008
#' ##
#' ## 2. Quantiles for each variable:
#' ##
#' ##             2.5%    25%    50%    75%  97.5%
#' ## a.h1     0.20575 0.8726 1.1892 1.4226 2.1530
#' ## v.f1.h1 -3.28390 1.6922 2.4203 3.2849 4.7110
#' ## v.f2.h1 -0.66095 1.5270 2.0492 2.5424 4.3922
#' ## z.h1     0.03883 0.2603 0.4218 0.5918 1.2787
#' ## sz.h1    0.03091 0.1552 0.3271 0.4464 0.7030
#' ## sv.h1    0.08784 0.3679 0.4696 0.6010 0.8265
#' ## t0.h1    0.02915 0.1352 0.2287 0.3214 0.7357
#' ## a.h2     0.13628 0.4108 0.6441 1.0149 1.5645
#' ## v.f1.h2  0.33562 0.8834 1.1464 1.3299 1.7378
#' ## v.f2.h2  0.24369 0.6600 0.8265 1.0120 1.3514
#' ## z.h2     0.10954 0.3481 0.6890 1.1156 1.9091
#' ## sz.h2    0.13985 0.3030 0.4053 0.5025 0.8777
#' ## sv.h2    0.28573 0.8155 1.0124 1.3448 1.6827
#' ## t0.h2    0.05304 0.1473 0.2531 0.3969 0.7472
#'
#' summary(hsamples0, hyper.means=TRUE)
#' ##           a      v.f1      v.f2         z        sz        sv        t0
#' ## h1 1.172753 1.9653973 0.3164101 0.2533810 1.0990106 0.7673155 1.0427178
#' ## h2 2.273730 0.4848155 0.4784983 0.7266657 0.8329073 0.4201312 0.2912233
summary.hyper <- function(object, start=1, end=NA, hyper.means=FALSE, ...)
{
  hyper <- attr(object, "hyper")
  if (is.na(end)) end <- hyper$nmc
  if (hyper.means) {
    tmp <- summary(phi.as.mcmc.list(hyper,start=start,end=end))
    matrix(tmp$statistics[,"Mean"], nrow=2,
      dimnames=list(c("h1","h2"), object[[1]]$p.names))
  } else {
    summary(phi.as.mcmc.list(hyper,start=start,end=end))
  }
}



### Priors, likelihoods and model selection

pll.dmc <- function(samples,hyper=FALSE,start=1,end=NA,prior=FALSE,like=FALSE,subject=NA)
  # extracts posterior log-likelihoods
{
  if ( hyper ) {
    hyper <- attr(samples,"hyper")
    if (is.null(hyper))
      stop("There are no hyper-parameters to plot.")
    if ( is.na(end) ) end <- hyper$nmc
    if ( end <= start )
      stop("End must be greater than start")
    if (prior) out <- hyper$h_summed_log_prior[start:end,] else
    if (like) out <-  hyper$h_log_likelihoods[start:end,] else
    out <- hyper$h_summed_log_prior[start:end,] +
           hyper$h_log_likelihoods[start:end,]
    dimnames(out) <- list(start:end,1:dim(out)[2])
  } else {
    if ( is.null(samples$theta) & !is.na(subject) )
      samples <- samples[[subject]]
    if ( is.null(samples$theta) ) { # multiple subjects
      if ( is.na(end) ) end <- samples[[1]]$nmc
      if ( end <= start ) stop("End must be greater than start")
      out <- lapply(samples,function(x){
      if (prior) out <- x$summed_log_prior[start:end,] else
        if (like) out <- x$log_likelihoods[start:end,] else
          out <- x$summed_log_prior[start:end,] + x$log_likelihoods[start:end,]
        dimnames(out) <- list(start:end,1:dim(out)[2])
        out
     })
    } else { # single subejct
      if ( is.na(end) ) end <- samples$nmc
      if ( end <= start ) stop("End must be greater than start")
      if (prior) out <- samples$summed_log_prior[start:end,] else
      if (like) out <-  samples$log_likelihoods[start:end,] else
      out <- samples$summed_log_prior[start:end,] +
             samples$log_likelihoods[start:end,]
      dimnames(out) <- list(start:end,1:dim(out)[2])
    }
  }
  out
}

#' Calculate Dstats of DDM Density
#'
#' Calculates and returns list with mean (meanD), variance (varD) and min
#' (minD) of Deviance, deviance of mean theta (Dmean) and (optionally) matrix
#' of deviances for each theta.
#'
#' This function is still under tweaking. Please use DMC's \code{Dstat.dmc},
#' instead.
#'
#' @param samples a DMC sample/object
#' @param save a save switch
#' @param fast choose different calculation routine
#' @keywords Dstats
#' @export
Dstats.ddm <- function(samples, save=FALSE, fast=TRUE) {
  if (fast) {
    summed_ll <- apply(samples$log_likelihoods, c(1,2), sum) ## This done nothing?
    D <- -2*samples$log_likelihoods
  } else {
    D <- apply(samples$theta, c(3,1), function(x) {
    -2*sum(log(likelihood.ddm(x, samples$data)))})
  }

  mtheta <- apply(samples$theta, 2, mean)
  Dmean <- -2*sum(log(likelihood.ddm(mtheta, samples$data)))
  minD <- min(D)

  if (save)
    list(np=dim(samples$theta)[2],meanD=mean(D),varD=var(as.vector(D)),
         minD=min(D),Dmean=Dmean,D=D) else
    list(np=dim(samples$theta)[2],meanD=mean(D),varD=var(as.vector(D)),
         minD=min(D),Dmean=Dmean,D=NA)
}


pd.dmc <- function(ds)
  # effective number of parameters calculated by mean, min and var methods
{
  list(Pmean=ds$meanD-ds$Dmean,Pmin=ds$meanD-ds$minD,Pvar=ds$varD/2)
}


posterior.lr.dmc <- function(D1,D2,main="",
                             plot=FALSE,plot.density=TRUE)
  # Aitkins posterior deviance likelihood ratio test
{
  if (is.list(D1)) D1 <- D1$D
  if (is.list(D2)) D2 <- D2$D
  n <- pmin(length(D1),length(D2))
  dD <- D1[1:n]-D2[1:n]
  if (plot) if (plot.density)
    plot(density(dD),xlab="D1-D2",main=main) else
      hist(dD,breaks="fd",xlab="D1-D2",main=main)
  if (min(dD)<0 & max(dD)>0) abline(v=0)
  c(pD1=mean(dD<0))
}


IC.ddm <- function(ds=NULL,samples=NULL,DIC=FALSE,fast=TRUE,use.pd=NA)
  # Calcualte IC1 (i.e., BPIC, default) or DIC
{
  if (is.null(ds)) if (is.null(samples))
    stop("Must supply samples or deviance statistics") else
      ds <- Dstats.ddm(samples,TRUE,fast)
    pds <- pd.dmc(ds)
    if (is.na(use.pd)) {
      if (ds$minD < ds$Dmean) pd <- pds$Pmin else pd <- pds$Pmean
    } else {
      if (use.pd %in% names(pds))
        pd <- pds[[use.pd]] else
          stop(paste("use.pd must be one of:",paste(names(pds),collapse=",")))
    }
    if (DIC) ds$meanD+pd else ds$meanD+2*pd
}


wIC.dmc <- function(ds=NULL,samples=NULL,
                    DIC=FALSE,fast=TRUE,use.pd=NA,...)
  # Calculate weights for a set of models
{

  ICs <- function(samples,DIC=FALSE,fast=TRUE,use.pd=NA)
    IC.dmc(samples=samples,DIC=DIC,fast=fast,use.pd=use.pd)

  if (is.null(ds)) if (is.null(samples))
    stop("Must supply samples or deviance statistics")
  if (is.null(samples))
    ics <- unlist(lapply(ds,IC.dmc,DIC=DIC,fast=fast,use.pd=use.pd)) else
      ics <- unlist(lapply(ds,ICs,DIC=DIC,fast=fast,use.pd=use.pd))
    d = ics-min(ics)
    w = exp(-d/2)/sum(exp(-d/2))
    if (!is.null(names(ics))) mnams=names(ics) else
      mnams <- paste("M",1:length(d),sep="")
    print(t(matrix(c(d,w),nrow=length(d),
             dimnames=list(mnams,c("IC-min","w")))),...)
}


waic.dmc <- function(trial_log_likes,mc_se=FALSE,save=FALSE,...)
  # Calclaute WAIC
{
  cat("From ggdmc\n")
  out <- waic(matrix(trial_log_likes,ncol=dim(trial_log_likes)[3]))
  if (mc_se) {
    n.chains <- dim(trial_log_likes)[2]
    waics <- vector(mode="list",length=n.chains)
    for (i in 1:n.chains) waics[[i]] <- waic(trial_log_likes[,i,])
    mc <- sd(unlist(lapply(waics,function(x){x$waic})))/sqrt(n.chains)
    print(c(p=out$p_waic,se_p=out$se_p,waic=out$waic,se_waic=out$se_waic,
      mc_se_waic=mc),...)
  } else
    print(c(p=out$p_waic,se_p=out$se_p,waic=out$waic,se_waic=out$se_waic),...)
  if (save) out
}

#' @importFrom loo loo
looic.dmc <- function(trial_log_likes,mc_se=FALSE,save=FALSE,...)
  # Calcuate looic
{
  out <- loo(matrix(trial_log_likes,ncol=dim(trial_log_likes)[3]))
  if (all(out$pareto_k<.5))
    cat("All Pareto k estimates OK (k < 0.5)\n") else {
    msg <- "See PSIS-LOO description (?'loo-package') for more information"
    if (any(out$pareto_k>1))
      msg1 <- paste(sum(out$pareto_k>1)," (",
                    round(100*sum(out$pareto_k>1)/length(out$pareto_k),1),
                    "%) Pareto k estimates greater than 1\n",sep="") else
      msg1 <- paste(sum(out$pareto_k>0.5)," (",
                    round(100*sum(out$pareto_k>1)/length(out$pareto_k),1),
                    "%) Pareto k estimates estimates between 0.5 and 1\n",sep="")
     warning(msg1,msg)
  }
  if (mc_se) {
    n.chains <- dim(trial_log_likes)[2]
    loos <- vector(mode="list",length=n.chains)
    for (i in 1:n.chains) loos[[i]] <- loo(trial_log_likes[,i,])
    mc <- sd(unlist(lapply(loos,function(x){x$looic})))/sqrt(n.chains)
    print(c(p=out$p_loo,se_p=out$se_p,looic=out$looic,se_looic=out$se_looic,
      mc_se_loo=mc),...)
  } else
    print(c(p=out$p_loo,se_p=out$se_p,looic=out$looic,se_looic=out$se_looic),...)
  if (save) out
}

#' @importFrom loo compare
loocompare.dmc <- function(loo1,loo2=NULL,...)
  # Model comparision of objects produced by waic.dmc or looic.dmc
{
  if ( !is.null(loo2) ) {
    tmp <- compare(loo1,loo2)
    out <- c(waic_diff=-2*tmp[[1]],se=2*tmp[[2]])
    print(out,...)
  } else {
    if ( !is.list(loo1) )
      stop("loo1 argument must be a list if loo2 not given")
    if ( any(names(loo1[[1]])=="looic") ) icnam <- "looic" else
                                          icnam <- "waic"
    ics <- unlist(lapply(loo1,function(x){x[[icnam]]}))
    d = ics-min(ics)
    w = exp(-d/2)/sum(exp(-d/2))
    if (!is.null(names(ics))) mnams=names(ics) else
      mnams <- paste("M",1:length(d),sep="")
    print(t(matrix(c(d,w),nrow=length(d),
      dimnames=list(mnams,c("IC-min","w")))),...)
  }
}

p.fun.dmc <- function(samples,fun,hyper=FALSE,ptype=1)
  # applies fun to samples
  {
   if (!hyper) as.vector(apply(samples$theta,c(1,3),fun)) else
     as.vector(apply(attr(samples,"hyper")$phi[[ptype]],c(1,3),fun))
}

### Plausible values

#' @importFrom hypergeo genhypergeo
posteriorRho <- function(r, n, spacing=.01, kappa=1)
  # Code provided by Dora Matzke, March 2016, from Alexander Ly
  # Reformatted into a single funciton
  # Defaults provide a grid with .01 resolution, kappa=1 implies uniform prior
{

  .bf10Exact <- function(n, r, kappa=1) {
	  # Ly et al 2015
	  # This is the exact result with symmetric beta prior on rho
	  # with parameter alpha. If kappa = 1 then uniform prior on rho
	  #
	  if (n <= 2){
		  return(1)
	  } else if (any(is.na(r))){
		  return(NaN)
	  }
	  # TODO: use which
	  check.r <- abs(r) >= 1 # check whether |r| >= 1
	  if (kappa >= 1 && n > 2 && check.r) {
		  return(Inf)
	  }

	  log.hyper.term <- log(hypergeo::genhypergeo(U=c((n-1)/2, (n-1)/2),
	                                              L=((n+2/kappa)/2), z=r^2))
	  log.result <- log(2^(1-2/kappa))+0.5*log(pi)-lbeta(1/kappa, 1/kappa)+
		lgamma((n+2/kappa-1)/2)-lgamma((n+2/kappa)/2)+log.hyper.term
	  real.result <- exp(Re(log.result))
	  return(real.result)
  }

  .jeffreysApproxH <- function(n, r, rho) {
	  result <- ((1 - rho^(2))^(0.5*(n - 1)))/((1 - rho*r)^(n - 1 - 0.5))
	  return(result)
  }

  .bf10JeffreysIntegrate <- function(n, r, kappa=1) {
	# Jeffreys' test for whether a correlation is zero or not
	# Jeffreys (1961), pp. 289-292
	# This is the exact result, see EJ
	##
	if (n <= 2){
		return(1)
	} else if ( any(is.na(r)) ){
		return(NaN)
	}

	# TODO: use which
	if (n > 2 && abs(r)==1) {
		return(Inf)
	}
	hyper.term <- Re(hypergeo::genhypergeo(U=c((2*n-3)/4, (2*n-1)/4), L=(n+2/kappa)/2, z=r^2))
	log.term <- lgamma((n+2/kappa-1)/2)-lgamma((n+2/kappa)/2)-lbeta(1/kappa, 1/kappa)
	result <- sqrt(pi)*2^(1-2/kappa)*exp(log.term)*hyper.term
	return(result)
}


  # 1.0. Built-up for likelihood functions
  .aFunction <- function(n, r, rho) {
	  #hyper.term <- Re(hypergeo::hypergeo(((n-1)/2), ((n-1)/2), (1/2), (r*rho)^2))
	  hyper.term <- Re(hypergeo::genhypergeo(U=c((n-1)/2, (n-1)/2), L=(1/2), z=(r*rho)^2))
	  result <- (1-rho^2)^((n-1)/2)*hyper.term
	  return(result)
  }

  .bFunction <- function(n, r, rho) {
	  #hyper.term.1 <- Re(hypergeo::hypergeo((n/2), (n/2), (1/2), (r*rho)^2))
	  #hyper.term.2 <- Re(hypergeo::hypergeo((n/2), (n/2), (-1/2), (r*rho)^2))
	  #hyper.term.1 <- Re(hypergeo::genhypergeo(U=c(n/2, n/2), L=(1/2), z=(r*rho)^2))
	  #hyper.term.2 <- Re(hypergeo::genhypergeo(U=c(n/2, n/2), L=(-1/2), z=(r*rho)^2))
	  #result <- 2^(-1)*(1-rho^2)^((n-1)/2)*exp(log.term)*
	  #	((1-2*n*(r*rho)^2)/(r*rho)*hyper.term.1-(1-(r*rho)^2)/(r*rho)*hyper.term.2)
	  #
	  hyper.term <- Re(hypergeo::genhypergeo(U=c(n/2, n/2), L=(3/2), z=(r*rho)^2))
	  log.term <- 2*(lgamma(n/2)-lgamma((n-1)/2))+((n-1)/2)*log(1-rho^2)
	  result <- 2*r*rho*exp(log.term)*hyper.term
	  return(result)
  }

  .hFunction <- function(n, r, rho) {
	  result <- .aFunction(n, r, rho) + .bFunction(n, r, rho)
	  return(result)
  }

  .scaledBeta <- function(rho, alpha, beta){
	  result <- 1/2*dbeta((rho+1)/2, alpha, beta)
	  return(result)
  }


  .priorRho <- function(rho, kappa=1) {
	  .scaledBeta(rho, 1/kappa, 1/kappa)
  }

  # Main body
  rho <- seq(-1,1,spacing)
  if (!is.na(r) && !r==0){
		return(spacing/.bf10Exact(n, r, kappa)*.hFunction(n, r, rho)*.priorRho(rho, kappa))
	} else if (!is.na(r) && r==0){
		return(spacing/.bf10JeffreysIntegrate(n, r, kappa)*
		         .jeffreysApproxH(n, r, rho)*.priorRho(rho, kappa))
	}

}

postRav <- function(r, n, spacing=.01, kappa=1)
  # r is a vector, returns average density.
{
  result <- apply(sapply(r,posteriorRho,n=n,spacing=spacing,kappa=kappa),1,mean)
  names(result) <- seq(-1,1,spacing)
  attr(result,"n") <- n
  attr(result,"kappa") <- kappa
  result
}

postRav.Density <- function(result)
  # Produces desnity class object
{
  x.vals <- as.numeric(names(result))
  result <- result/(diff(range(x.vals))/length(x.vals))
  out <- list(x=x.vals,y=result,has.na=FALSE,
    data.name="postRav",call=call("postRav"),
    bw=mean(diff(x.vals)),n=attr(result,"n"))
  class(out) <- "density"
  out
}

postRav.mean <- function(pra) {
  # Average value of object produced by posteriorRhoAverage
  sum(pra*as.numeric(names(pra)))
}

postRav.p <- function(pra,lower=-1,upper=1) {
  # probability in an (inclusive) range of posteriorRhoAverage object
  x.vals <- as.numeric(names(pra))
  sum(pra[x.vals <= upper & x.vals >= lower])
}

