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
#' @importFrom stats median
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
#' hsamples0 <- h.samples.dmc(nmc=10, p.prior=p.prior, pp.prior=pp.prior,
#'   data=mdi1, thin=1)
#' hsamples0 <- h.run.dmc(hsamples0)
#' gelman.diag.dmc(hsamples0, hyper=TRUE)
#' gelman.diag.dmc(hsamples0)
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
#' dat  <- h.simulate.dmc(m1, nsim=30, ns=4, p.prior=pop.prior)
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
#' hsamples0 <- h.samples.dmc(nmc=10, p.prior=p.prior, pp.prior=pp.prior,
#'   data=mdi1, thin=1)
#' hsamples0 <- h.run.dmc(hsamples0)
#' es <- effectiveSize.dmc(hsamples0)
#' hes <- effectiveSize.dmc(hsamples0, hyper=TRUE)
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
#' dat1 <- simulate(m1, nsim=30, p.vector=pVec)
#' mdi1 <- data.model.dmc(dat1, m1)
#'
#' samples0 <- samples.dmc(nmc=100, p.prior=p.prior, data=mdi1)
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
#' m1 <- model.dmc(
#'       p.map     = list(a="1",v="F",z="1",d="1",sz="1",sv="1",t0="1",
#'                        st0="1"),
#'       match.map = list(M=list(s1="r1",s2="r2")),
#'       factors   = list(S=c("s1","s2"),F=c("f1","f2")),
#'       constants = c(st0=0,d=0),
#'       responses = c("r1","r2"),
#'       type      = "rd")
#'
#' pop.mean  <- c(a=1.15, v.f1=1.25, v.f2=1.85, z=0.55, sz=0.15, sv=0.32,
#'                t0=0.25)
#' pop.scale <- c(a=0.10, v.f1=.8,   v.f2=.5,   z=0.1,  sz=0.05, sv=0.05,
#'                t0=0.05)
#' pop.prior <- prior.p.dmc(
#'   dists = rep("tnorm", length(pop.mean)),
#'   p1    = pop.mean,
#'   p2    = pop.scale,
#'   lower = c(0,-5, -5, 0, 0,   0, 0),
#'   upper = c(5, 7,  7, 1, 0.5, 2, 2))
#'
#' dat  <- h.simulate.dmc(m1, nsim=30, ns=4, p.prior=pop.prior)
#' mdi1 <- data.model.dmc(dat, m1)
#' ps   <- attr(dat,  "parameters")
#'
#' p.prior <- prior.p.dmc(
#'   dists= rep("tnorm", length(pop.mean)),
#'   p1=pop.mean,
#'   p2=pop.scale*5,
#'   lower=c(0,-5, -5, 0, 0, 0, 0),
#'   upper=c(5, 7,  7, 2, 2, 2, 2))
#' samples0 <- h.samples.dmc(nmc=30, p.prior=p.prior, data=mdi1, thin=1)
#' samples0 <- h.run.dmc(samples0)
#' class(samples0)
#' ## [1] "dmc.list"
#' gelman.diag.dmc(samples0)
#'
#' ## summary calls theta.as.mcmc.list, which is very slow.
#' ## summary(samples0)
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
#' dat  <- h.simulate.dmc(m1, nsim=30, ns=4, p.prior=pop.prior)
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
#' hsamples0 <- h.samples.dmc(nmc=30, p.prior=p.prior, pp.prior=pp.prior,
#'   data=mdi1, thin=1)
#' hsamples0 <- h.run.dmc(hsamples0, p.migrate=.05, h.p.migrate=.05)
#'
#' class(hsamples0)
#' ## [1] "hyper"
#'
#' ## summary.hyper calls phi.as.mcmc.list, which is very slow.
#' ## summary(hsamples0)
#' ## summary(hsamples0, hyper.means=TRUE)
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

pll.dmc <- function(samples,hyper=FALSE,start=1,end=NA,prior=FALSE,
  like=FALSE,subject=NA)
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
#' @importFrom stats var
#' @export
Dstats.ddm <- function(samples, save=FALSE, fast=TRUE) {
  if (fast) {
    summed_ll <- apply(samples$log_likelihoods, c(1,2), sum) ## This done nothing?
    D <- -2*samples$log_likelihoods
  } else {
    D <- apply(samples$theta, c(3,1), function(x) {
    -2*sum(log(likelihood.rd(x, samples$data)))})
  }

  mtheta <- apply(samples$theta, 2, mean)
  Dmean <- -2*sum(log(likelihood.rd(mtheta, samples$data)))
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

#' @importFrom graphics plot abline
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


# wIC.dmc <- function(ds=NULL,samples=NULL,
#                     DIC=FALSE,fast=TRUE,use.pd=NA,...)
#   # Calculate weights for a set of models
# {
#
#   ICs <- function(samples,DIC=FALSE,fast=TRUE,use.pd=NA)
#     IC.dmc(samples=samples,DIC=DIC,fast=fast,use.pd=use.pd)
#
#   if (is.null(ds)) if (is.null(samples))
#     stop("Must supply samples or deviance statistics")
#   if (is.null(samples))
#     ics <- unlist(lapply(ds,IC.dmc,DIC=DIC,fast=fast,use.pd=use.pd)) else
#       ics <- unlist(lapply(ds,ICs,DIC=DIC,fast=fast,use.pd=use.pd))
#     d = ics-min(ics)
#     w = exp(-d/2)/sum(exp(-d/2))
#     if (!is.null(names(ics))) mnams=names(ics) else
#       mnams <- paste("M",1:length(d),sep="")
#     print(t(matrix(c(d,w),nrow=length(d),
#              dimnames=list(mnams,c("IC-min","w")))),...)
# }


# waic.dmc <- function(trial_log_likes,mc_se=FALSE,save=FALSE,...)
#   # Calclaute WAIC
# {
#   cat("From ggdmc\n")
#   out <- waic(matrix(trial_log_likes,ncol=dim(trial_log_likes)[3]))
#   if (mc_se) {
#     n.chains <- dim(trial_log_likes)[2]
#     waics <- vector(mode="list",length=n.chains)
#     for (i in 1:n.chains) waics[[i]] <- waic(trial_log_likes[,i,])
#     mc <- sd(unlist(lapply(waics,function(x){x$waic})))/sqrt(n.chains)
#     print(c(p=out$p_waic,se_p=out$se_p,waic=out$waic,se_waic=out$se_waic,
#       mc_se_waic=mc),...)
#   } else
#     print(c(p=out$p_waic,se_p=out$se_p,waic=out$waic,se_waic=out$se_waic),...)
#   if (save) out
# }

#' @importFrom loo loo
#' @importFrom stats sd
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

