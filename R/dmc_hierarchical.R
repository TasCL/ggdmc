# System functions for the DMC (Dynamic Models of Choice)
#    Functions for hierarchical models (versions of non-hierarchical functions
#    and hierarchical model-specific functions)
#    Usually user does not need to edit


#' Simulate Choice-RT Data for Multiple Participants
#'
#' \code{h.simulate.dmc} generates random choice-RT responses based on an EAM
#' model, which usually set up by \code{model.dmc} and a prior distribution
#' settings, which usually created by \code{prior.p.dmc}.
#'
#' The simulation engine of the diffusoin model calls \code{rdiffusion} in
#' \pkg{rtdists}, which is adapted Voss & Voss' \code{construct-sample} C
#' function in their \pkg{fast-dm} C software
#'
#' \code{h.simulate.dmc} behaves like \code{simulate.dmc}, but creates data
#' for a set of \code{ns} subjects. It simulates based either on a
#' fixed-effect model or on a random-effect model. For the latter,
#' (i.e., hierarchical model), the user must supply a \code{p.prior}, in which
#' case \code{ps}, (used like \code{p.vector} in \code{simulate.dmc}) is
#' ignored.  That is, \code{h.simulate.dmc} will generate a set of true EAM
#' parameters stochastically based on the supplied prior distributions (stored
#' in the user supplied \code{p.prior}).
#'
#' \code{ps} stands for "true" model \emph{p}arameter\emph{s}. It can be a
#' vector exactly like \code{p.vector}, in which case each subject has
#' identical parameters. All equal to the values in \code{p.vector}. \code{ps}
#' can also be a matrix with one row per subject, in which case must have
#' \code{ns} rows. It is saved (in expanded form) as "parameters" attribute in
#' the generated data.
#'
#' \code{p.prior} is a list of distributions from which subject parameters
#' are sampled. It is usually created by \code{prior.p.dmc} and will be saved
#' as \code{p.prior} attribute.
#'
#' \code{nsim} can be a single number for a balanced design or matrix for an
#' unbalanced design, where rows are subjects and columns are design cells. If
#' the matrix has one row then all subjects have the same n in each cell, if it
#' has one column then all cells have the same n, otherwise each entry
#' specifies the n for a particular design subject x design cell combination.
#'
#' @param object a DMC model
#' @param nsim number of trials. Default is 2
#' @param seed for comparible reason. Default is NULL.
#' @param ns the number of subjects to be simulated.  Default is 1
#' @param ps a parameter x subject matrix
#' @param SSD stop-signal parameter. is for use only with stop-signal designs,
#' specified like n except a vector form is also allowed that is the same
#' length as the data. It must have Inf in all go cells
#' @param p.prior prior parameter list
#' @param staircase staircase modifies SSD (see simulate.dmc).
#' @param subject.cv a data frame, column names are in names(p.prior),
#' @keywords h.simulate.dmc
#' @export
#' @examples
#' ## Set up a DDM Model, rate effect of factor F
#' m1 <- model.dmc(
#'  p.map     = list(a="1", v="F", z="1", d="1", sz="1", sv="1", t0="1",
#'  st0="1"),
#'  match.map = list(M=list(s1="r1", s2="r2")),
#'  factors   = list(S=c("s1","s2"), F=c("f1","f2")),
#'  constants = c(st0=0,d=0),
#'  responses = c("r1","r2"),
#'  type      = "rd")
#'
#' ## Population distribution
#' p.mean  <- c(a=2,   v.f1=2.5, v.f2=1.5, z=0.5, sz=0.3, sv=1,  t0=0.3)
#' p.scale <- c(a=0.5, v.f1=.5,  v.f2=.5,  z=0.1, sz=0.1, sv=.3, t0=0.05)
#' pop.prior <- prior.p.dmc(
#'   dists = rep("tnorm",7),
#'   p1    = p.mean,
#'   p2    = p.scale,
#'   lower = c(0,-5, -5, 0, 0, 0, 0),
#'   upper = c(5, 7,  7, 2, 2, 2, 2))
#'
#' ##  Check population distributions
#' plot_priors(pop.prior)
#'
#' ## Simulate some data
#' dat1 <- h.simulate.dmc(m1, nsim=20, ns=5, p.prior=pop.prior)
#' head(dat1)
#' ##    S  F  R        RT
#' ## 1 s1 f1 r1 0.9227881
#' ## 2 s1 f1 r1 0.7878554
#' ## 3 s1 f1 r1 0.4814711
#' ## 4 s1 f1 r1 0.6864110
#' ## 5 s1 f1 r1 0.5068179
#' ## 6 s1 f1 r1 0.6356547
#'
#' ## Use LBA as an example
#' m2 <- model.dmc(
#'       p.map     = list(A="1", B="1", mean_v="M", sd_v="M", t0="1",
#'       st0="1"),
#'       match.map = list(M=list(s1=1, s2=2)),
#'       factors   = list(S=c("s1", "s2")),
#'       constants = c(st0= 0, sd_v.false=1),
#'       responses = c("r1", "r2"),
#'       type      = "norm")
#'
#' ## Population distribution
#' p.mean  <- c(A=.4,B=.6,mean_v.true=1,mean_v.false=0,sd_v.true = .5,t0=.3)
#' p.scale <- c(A=.1,B=.1,mean_v.true=.2,mean_v.false=.2,sd_v.true = .1,
#' t0=.05)
#' pop.prior <- prior.p.dmc(
#'     dists = c("tnorm","tnorm","tnorm","tnorm","tnorm","tnorm"),
#'     p1=p.mean, p2=p.scale,
#'     lower=c(0,0,NA,NA,0,.1),upper=c(NA,NA,NA,NA,NA,1))
#'
#' ## A data frame
#' dat2 <- h.simulate.dmc(m2, nsim=30, p.prior=pop.prior, ns=10) ##
#'
#' ## A 10-element list; each element is a participant
#' mdi2 <- data.model.dmc(dat2, m2)
#' head(mdi2[[1]])     ## Check the first participant (s1)
#' ##    S  R        RT
#' ## 1 s1 r1 0.7793124
#' ## 2 s1 r2 0.5736594
#' ## 3 s1 r2 0.6900489
#' ## 4 s1 r2 2.3713993
#' ## 5 s1 r1 1.5139890
#' ## 6 s1 r2 0.5649499
h.simulate.dmc <- function(object, nsim=2, seed=NULL, ns=1, ps=NA, SSD=Inf,
  p.prior=NA, staircase=NA, subject.cv=NULL)
{
  model <- object
  # Check nsim
  ncell <- prod(unlist(lapply(attr(model, "factors"), length)))
  if (!is.matrix(nsim) & (is.vector(nsim) & length(nsim) != 1))
    stop("nsim must be a scalar or matrix")
  if (is.matrix(nsim))
  {
    if (dim(nsim)[1]==1) # only cells differ
      nsim <- matrix(rep(nsim, each=ns), nrow=ns)
    if (dim(nsim)[2]==1) # only subjects differ
      nsim <- matrix(rep(nsim,times=ncell),nrow=ns)
  } else nsim <- matrix(rep(nsim,ns*ncell),nrow=ns)
  if ( ns != dim(nsim)[1] )
    stop(paste("The nsim matrix must have",ns,"rows"))
  if ( ncell != dim(nsim)[2] )
    stop(paste("The nsim matrix must have",ncell,"columns"))

  # Check SSD
  ndata <- sum(nsim)
  if ( is.matrix(SSD) )
  {
    if ( dim(SSD)[1]==1 ) # only cells differ
      SSD <- matrix(rep(SSD,each=ns),nrow=ns)
    if ( dim(SSD)[2]==1 ) # only subjects differ
      SSD <- matrix(rep(SSD,times=ncell),nrow=ns)
    if ( ns != dim(SSD)[1] )
      stop(paste("The SSD matrix must have",ns,"rows"))
    if ( ncell != dim(SSD)[2] )
      stop(paste("The SSD matrix must have",ncell,"columns"))
    SSD <- rep(t(SSD),as.vector(nsim))
  } else if ( length(SSD)==1 ) SSD <- rep(SSD,ndata)
  if ( length(SSD) != ndata )
    stop(paste("SSD vector is not the same length as the data (",ndata,")"))

  # check ps/p.prior
  if ( any(is.na(p.prior)) )
  {
    if ( !is.matrix(ps) )
    {
      if (check.p.vector(ps,model)) stop()
      ps <- matrix(rep(ps, each=ns),nrow=ns,dimnames=list(NULL,names(ps)))
    }
    if ((ns != dim(ps)[1]))
      stop("ps matrix must have ns rows")
    if (check.p.vector(ps[1,],model)) stop()
  } else { # Check p.prior
    if (!all(sort(names(attr(model,"p.vector")))==sort(names(p.prior))))
      stop("p.prior incompatible with model")
  }

  # random effects
  if ( !any(is.na(p.prior)) )
    ps <- as.matrix(rprior.dmc(p.prior,ns))

  # Subject covariates
  if ( !is.null(subject.cv) ) {
    if ( !is.data.frame(subject.cv) || (dim(subject.cv)[1]!=ns) ||
      !all(names(subject.cv %in% dimnames(ps)[[2]])) )
      stop("subject.cv must be a data frame with names from p.prior and ns rows")
    for (i in names(subject.cv)) ps[,i] <- ps[,i] + subject.cv[[i]]
  }


  row.names(ps) <- 1:ns
  ndatai <- cumsum(c(0,apply(nsim,1,sum)))
  datr <- (ndatai[1]+1):(ndatai[2])
  data <- cbind(s=rep(1,length(datr)),
    simulate.dmc(model, nsim=nsim[1,], p.vector=ps[1,],
      SSD=SSD[datr], staircase = staircase))
  if (ns>1) for (i in 2:ns)
  {
    datr <- (ndatai[i]+1):(ndatai[i+1])
    data <- rbind(data,cbind(s=rep(i,length(datr)),
      simulate.dmc(model, nsim=nsim[i,], p.vector=ps[i,],
        SSD=SSD[datr],staircase = staircase)))
  }
  data$s <- factor(data$s)
  attr(data,"parameters") <- ps
  if (!any(is.na(p.prior))) attr(data,"p.prior") <- p.prior
  data
}


# Hierarchical prior, likelihoods and posterior

h.summed.log.prior=function (pp, pp.prior)
{
  suppressWarnings( sum(log.prior.dmc(pp[[1]],pp.prior[[1]])) +
                    sum(log.prior.dmc(pp[[2]],pp.prior[[2]])) )
}


h.log.likelihood.dmc <- function(ps,pp,p.prior)
  # log-likelihood of subject parameters ps under population distribuiton p.prior
{
  suppressWarnings( apply(apply(ps,1,log.prior.dmc,
                                p.prior=assign.pp(pp,p.prior)),2,sum) )
}


h.log.posterior.dmc <- function(ps,p.prior,pp,pp.prior)
  # log-likelihood of subject parameters ps under population distribution
  # p.prior,given population parameters pp (phi) with priors pp.prior
{
  # sum over subjects of likelihood of pp (pop pars) given ps (subjects pars)
  sum(h.log.likelihood.dmc(ps,pp,p.prior)) +
  # prior probability of pp
  h.summed.log.prior(pp,pp.prior)
}


make.hstart <- function(fsamples,
  mu.sd=1,sigma.sd=1,
  lower.mu=NULL,upper.mu=NULL,
  lower.sigma=NULL,upper.sigma=NULL,
  do.plot=FALSE,layout=c(2,5),digits=2)
  # Uses the results of fixed effects fitting to get a start prior for hierarchcial
  # By default sets prior sd to 1% of mean (tight), can specify as a vector
  # Prints prior parameters and can plot priors
  ## Change to ggdmc version
{

  mns <- lapply(fsamples,function(x){apply(x$theta,2,mean)})
  mns <- matrix(unlist(mns),nrow=length(mns[[1]]),dimnames=list(names(mns[[1]]),NULL))

  mn <- apply(mns,1,mean)
  if (length(mu.sd)==1) mu.sd <- abs(mn)*mu.sd/100
  if (length(mu.sd)!=length(mn))
    stop(paste("mu.sd must have length 1 or length of parameters:",length(mn)))
  if (is.null(lower.mu)) lower.mu <- rep(-Inf,length(mn))
  if (length(lower.mu)!=length(mn))
    stop(paste("lower.mu must have length of parameters:",length(mn)))
  if (is.null(upper.mu)) upper.mu <- rep(Inf,length(mn))
  if (length(upper.mu)!=length(mn))
    stop(paste("upper.mu must have length of parameters:",length(mn)))

  sd <- apply(mns,1,sd)
  if (length(sigma.sd)==1) sigma.sd <- sd*sigma.sd/100
  if (length(sigma.sd)!=length(sd))
    stop(paste("sigma.sd must have length 1 or length of parameters:",length(sd)))
  if (is.null(lower.sigma)) lower.sigma <- rep(-Inf,length(mn))
  if (length(lower.sigma)!=length(mn))
    stop(paste("lower.sigma must have length of parameters:",length(mn)))
  if (is.null(upper.sigma)) upper.sigma <- rep(Inf,length(mn))
  if (length(upper.sigma)!=length(mn))
    stop(paste("upper.sigma must have length of parameters:",length(mn)))

  mu.prior <- prior.p.dmc(p1=mn,p2=mu.sd,lower=lower.mu,upper=upper.mu)
  sigma.prior <- prior.p.dmc(p1=sd,p2=sigma.sd,lower=lower.sigma,upper=upper.sigma)

  cat("Mu prior\n")
  cat("Mean\n"); print(round(mn,digits))
  cat("SD\n"); print(round(mu.sd,digits))
  cat("Sigma prior\n")
  cat("Mean\n"); print(round(sd,digits))
  cat("SD\n"); print(round(sigma.sd,digits))

  if (do.plot) { plot_priors(mu.prior); plot_priors(sigma.prior) }

  list(mu=mu.prior,sigma=sigma.prior)
}


# samples=hsamples1;nmc=NA;thin=NA;max.try=100
# hstart.prior=NULL;phi1=NULL;start.from=NA

add.hyper.dmc <- function(samples, pp.prior, nmc=NA, thin=NA, max.try=100,
                          hstart.prior=NULL, phi1=NULL, p.prior=NULL)
  # Adds hyper level to fixed subject effect samples object
{

  update_constants <- function(samplesi,hyper,mci=1)
  # update summed_log_prior and log_likelihoods with no.sigma hypers
  {
    consts <- hyper$phi[[1]][,,mci][,!hyper$has.sigma,drop=FALSE]
    p.names <- dimnames(consts)[[2]]
    pp.prior <- hyper$pp.prior[[1]][p.names]
    datai <- samplesi$data
    for (i in 1:dim(samplesi$summed_log_prior)[2]) { # chains
      samplesi$summed_log_prior[mci,i] <- samplesi$summed_log_prior[mci,i] +
        summed.log.prior(consts[i,],pp.prior,p.names)
      attr(attr(datai,"model"),"all.par")[p.names] <- consts[i,]
      samplesi$log_likelihoods[mci,i] <-
        sum(log.likelihood(samplesi$theta[i,,mci],datai))
    }
    if (any(is.na(samplesi$summed_log_prior[mci,])))
      stop("prior=NA for hyper parameter with no sigma, narror prior?")
    if (any(is.na(samplesi$log_likelihoods[mci,])))
      stop("likelihood=NA for hyper parameter with no sigma, narror prior?")
    samplesi
  }


  if (names(samples[[1]])[1]!="theta")
    stop("samples must be a list of subject sample objects")
  if (is.null(p.prior)) p.prior <- samples[[1]]$p.prior else {
    for (i in 1:length(samples)) samples[[i]]$p.prior <- p.prior
  }
  # check pp.prior
  if (!is.list(pp.prior)) stop("pp.prior must be a list")
  if ( length(pp.prior[[1]])<length(pp.prior[[2]]) )
    stop("Hyper-prior must have as many or more location than scale parameters")
  has.sigma <- names(pp.prior[[1]]) %in% names(pp.prior[[2]])
  fixed <- names(pp.prior[[1]])[ !has.sigma ]
  bad <- fixed[!(fixed %in% names(attr(attr(samples[[1]]$data,"model"),"constants")))]
  if ( length(bad)>0 )
    stop(paste("Fixed hyper parameters not constants in data level:",bad))
  pp.names <- lapply(pp.prior,names)
  has.hyper <- names(p.prior) %in% pp.names[[2]]
  if ( !all(names(p.prior)[has.hyper] %in% pp.names[[1]][has.sigma]) )
      stop("pp.prior location parameters not compatible with p.prior")
  if ( !all(names(p.prior)[has.hyper]==pp.names[[2]]) )
      stop("pp.prior scale parameters not compatible with p.prior")
  if (!is.null(hstart.prior)) { # check hstart.prior
    if (!is.list(hstart.prior)) stop("hstart.prior must be a list")
    if (!all(sort(names(pp.prior[[1]]))==sort(names(hstart.prior[[1]]))))
        stop("Incompatible location hstart.prior")
    if (!all(sort(names(pp.prior[[2]]))==sort(names(hstart.prior[[2]]))))
        stop("Incompatible scale hstart.prior")
  }
  if (is.na(nmc)) nmc <- samples[[1]]$nmc
  if (is.na(thin)) thin <- samples[[1]]$thin
  samples <- lapply(samples,function(x){
    samples.dmc(nmc, samples=x, thin=thin)
  })
  n.chains <- samples[[1]]$n.chains
  p.names <- names(pp.prior[[1]])
  n.pars <- length(p.names)
  phi <- array(NA,c(n.chains,n.pars,nmc))
  dimnames(phi)[[2]] <- p.names
  phi <- list(phi,phi)  # pairs of sampled hyper-parameters
  if (!is.null(phi1)) { # check phi1
    if (!all(dim(phi[[1]][,,1])==dim(phi1[[1]])))
      stop("phi1 location not compatible")
    if (!all(dim(phi[[1]][,,1])==dim(phi1[[1]])))
      stop("phi1 scale not compatible")
  }
  h_summed_log_prior <- array(-Inf,c(nmc,n.chains)) # hyper log-likelihoods
  h_log_likelihoods <- array(-Inf,c(nmc,n.chains))  # hyper log-likelihoods
  ntry <- 1
  repeat {
    if ( is.null(phi1) ) {
      cat("Generating hyper-start points for each chain: ")
      for( i in 1:n.chains ) {
        cat(".")
        if ( is.null(hstart.prior) ) { # sample from prior
          phi[[1]][i,,1] <- rprior.dmc(pp.prior[[1]])[,p.names]
          phi[[2]][i,has.sigma,1] <-
            rprior.dmc(pp.prior[[2]])[,p.names[has.sigma]]
        } else {                        # sample from hstart.prior
          phi[[1]][i,,1] <- rprior.dmc(hstart.prior[[1]])[,p.names]
          phi[[2]][i,has.sigma,1] <-
            rprior.dmc(hstart.prior[[2]])[,p.names[has.sigma]]
        }
      }
      cat("\n")
    } else {
      phi[[1]][,,1] <- phi1[[1]]
      phi[[2]][,,1] <- phi1[[2]]
    }
    for (i in 1:n.chains) {
      h_summed_log_prior[1,i] <-
        sum(log.prior.dmc(phi[[1]][i,,1],pp.prior[[1]])) +
        sum(log.prior.dmc(phi[[2]][i,has.sigma,1],pp.prior[[2]]))
    }

    # Fill in hyper-likelihoods
    cps <- lapply(samples,function(x){x$theta[,,1]}) # subject list: chains x pars
    for( i in 1:n.chains ) {
      h_log_likelihoods[1,i] <- sum(h.log.likelihood.dmc(p.prior=p.prior[has.hyper],
        ps=t(data.frame(lapply(cps,function(x){x[i,has.hyper]}))),
        pp=list(phi[[1]][i,has.sigma,1],phi[[2]][i,has.sigma,1])))
    }
    if ( any(!is.finite(h_log_likelihoods[1,])) ) {
      ntry <- ntry + 1
      if (ntry > max.try)
        stop(paste("For",max.try,"attempts sampling start points was valid for data but not hyper level\n
          Try a tighter pp.prior or hstart.prior"))
    } else break
  }
  cat("\n")

  # Update subject level priors to fit with hyper ...
  for (j in 1:length(samples)) for ( i in 1:n.chains) {
    data <- samples[[j]]$data
    samples[[j]]$p.prior <- assign.pp(list(phi[[1]][i,has.sigma,1],
      phi[[2]][i,has.sigma,1]),p.prior=samples[[i]]$p.prior[has.hyper])
    samples[[j]]$summed_log_prior[1,i]  <- summed.log.prior (samples[[j]]$theta[i,,1],
      samples[[j]]$p.prior,p.names=names(attr(attributes(data)$model,"p.vector")))
  }

  hyper <- list(phi=phi,h_log_likelihoods=h_log_likelihoods,
    h_summed_log_prior=h_summed_log_prior,pp.prior=pp.prior,start=1,thin=thin,
    n.pars=n.pars,p.names=p.names,rp=samples[[1]]$rp,nmc=nmc,n.chains=n.chains,
    has.sigma=has.sigma,has.hyper=has.hyper)
  if (any(!has.sigma)) # Update data level for !has.sigma hyper parameters
    samples <- lapply(samples,update_constants,hyper=hyper)
  attr(samples,"hyper") <- hyper
  samples
}

#' Set up a DMC Sample with Multiple Participants
#'
#' \code{h.samples.dmc} initialise a DMC object with each particpant as a list
#' element in a list. A participant is himself/herself a list element, carrying
#' the DMC setting, such as theta, summed_log_likelihood, etc.
#'
#' @param nmc number of Markov Chain Monte Carlo iteration
#' @param p.prior prior distribution setting
#' @param data a model data instance created by \code{data.model.dmc}
#' @param pp.prior prior distribution setting, hyper level
#' @param samples a DMC posterior sample
#' @param thin thinning length. Default 1
#' @param theta1 A user supplied initial theta cube
#' @param phi1 A user supplied initial phi cube
#' @param start.prior A user supplied (different) prior distribution setting
#' @param hstart.prior A user supplied (different) hyper-prior
#' distribution setting
#' @param add whether add new MCMC iteration on top of previous DMC sample
#' @param rp a DEMCMC tuning parameter
#' @param setting a list container to store all DMC setting
#' @param verbose whether to print debugging information
#' @keywords h.samples.dmc
#' @export
#' @examples
#' ## Set up a DDM Model, rate effect of factor F
#' m1 <- model.dmc(
#'   p.map     = list(a="1",v="F",z="1",d="1",sz="1",sv="1",t0="1",st0="1"),
#'   match.map = list(M=list(s1="r1",s2="r2")),
#'   factors=list(S=c("s1","s2"), F=c("f1","f2")),
#'   constants = c(st0=0,d=0),
#'   responses = c("r1","r2"),
#'   type = "rd")
#'
#' ## Population distribution
#' pop.mean  <- c(a=2,   v.f1=2.5, v.f2=1.5, z=0.5, sz=0.3, sv=1,  t0=0.3)
#' pop.scale <- c(a=0.5, v.f1=.5,  v.f2=.5,  z=0.1, sz=0.1, sv=.3, t0=0.05)
#' pop.prior <- prior.p.dmc(
#'   dists = rep("tnorm",7),
#'   p1    = pop.mean,
#'   p2    = pop.scale,
#'   lower = c(0,-5, -5, 0, 0, 0, 0),
#'   upper = c(5, 7,  7, 2, 2, 2, 2))
#'
#' ##  Check population distributions
#' plot_priors(pop.prior)
#'
#' ## Simulate some data
#' dat <- h.simulate.dmc(m1, nsim=20, ns=4, p.prior=pop.prior)
#' mdi <- data.model.dmc(dat, m1)
#' head(dat)
#' ##    S  F  R        RT
#' ## 1 s1 f1 r1 0.9227881
#' ## 2 s1 f1 r1 0.7878554
#' ## 3 s1 f1 r1 0.4814711
#' ## 4 s1 f1 r1 0.6864110
#' ## 5 s1 f1 r1 0.5068179
#' ## 6 s1 f1 r1 0.6356547
#'
#' ## Take a look at true parameters
#' ps <- round( attr(dat, "parameters"), 2)
#' ps
#' ##      a v.f1 v.f2    z   sz   sv   t0
#' ## 1 2.83 2.91 1.41 0.66 0.30 0.65 0.30
#' ## 2 2.37 2.42 2.24 0.48 0.28 1.14 0.31
#' ## 3 1.91 2.49 0.98 0.74 0.33 1.20 0.18
#' ## 4 2.14 2.67 2.34 0.65 0.31 1.74 0.27
#'
#' ## FIT FIXED EFFECTS
#' ## specify a broader prior than the true population distribution
#' p.prior <- prior.p.dmc(
#'   dists= rep("tnorm", length(pop.mean)),
#'   p1=pop.mean,
#'   p2=pop.scale*5,
#'   lower=c(0,-5, -5, 0, 0, 0, 0),
#'   upper=c(5, 7,  7, 2, 2, 2, 2))
#'
#' ## Set up multiple participant DMC sample
#' samples0 <- h.samples.dmc(nmc=100, p.prior=p.prior, data=mdi, thin=1)
#'
#' ## Add 400 more iterations and change thinning length to 2
#' samples1 <- h.samples.dmc(nmc=100, p.prior=p.prior, samples=samples0,
#' thin=2, add=TRUE)
#'
#' samples1[[1]]$nmc
#' ## [1] 200
h.samples.dmc <- function(nmc, p.prior=NULL, data=NULL, pp.prior=NULL,
                samples=NULL, thin=1, theta1=NULL, phi1=NULL,
                start.prior=NULL, hstart.prior=NULL, add=FALSE, rp=.001,
                setting=NULL, verbose=FALSE)
{
  chainFamily <- 3

  if(!is.null(samples) & is.list(samples))
  {
    if( names(samples)[1] != "theta")
    {
      if(verbose)
      {
        cat("Take p.prior and pp.prior from 1st participant's p.prior\n")
        cat("and from samples's hyper\n")
      }
      pp.prior <- attr(samples, "hyper")$pp.prior
      p.prior  <- samples[[1]]$p.prior
      n.chains <- samples[[1]]$n.chains
    }
  }

  if(is.null(samples) & is.null(data))
    stop("Neither samples nor data was found!\n")
  model <- attr(data, "model")

  if( is.null(setting) & is.null(samples) )
  {
    if(is.null(model)) ## Correct for multiple participants
    {
      if(verbose)
        cat("Reset model using the model from the first participants.\n")
      model    <- attr(data[[1]], "model")
      n.chains <- length(attr(attr(data[[1]], "model"),"p.vector"))*3
    }
    p.names <- names(attr(model, "p.vector"))

    n.chains <- chainFamily*length(p.names)
    if (verbose) {
      cat("Apply default setting: add=FALSE, verbose=TRUE, ")
      cat("restart=TRUE, rp=0.001, thin=1, n.chains=", n.chains)
      cat(", start from = 1\n")
    }
    setting_in <- list(add=add, verbose=verbose, rp=rp, restart=TRUE,
      thin=thin, start.from=1, n.chains=n.chains, remove=FALSE)
  } else if (is.null(setting) & !is.null(samples) & names(samples)[1] != "theta")
  {
    if (verbose) {
       cat("Get setting from the first subject in a multi-subject samples and restart\n")
    }
    thin       <- ifelse(is.null(thin), samples[[1]]$thin, thin)
    new.start  <- dim(samples[[1]]$theta)[3]
    setting_in <- list(add=add, verbose=verbose, rp=samples[[1]]$rp,
      restart=TRUE, thin=thin, start.from=new.start,
      n.chains=samples[[1]]$n.chains, remove=FALSE)
  } else if (is.null(setting) & !is.null(samples) & names(samples)[1] == "theta")
  {
    n.chains <- samples$n.chains
    if (verbose) cat("Data: Get setting from samples and restart\n")
    setting_in <- list(add=add, verbose=verbose,
      rp=samples$rp, thin=samples$thin, restart=TRUE,
      start.from=dim(samples$theta)[3], n.chains=samples$n.chains,
      remove=FALSE)
  } else {
    setting$add     <- ifelse(is.null(setting$add), FALSE, setting$add)
    setting$verbose <- ifelse(is.null(setting$verbose), verbose, setting$verbose)
    setting$rp      <- ifelse(is.null(setting$rp), 0.001, setting$rp)
    setting$thin       <- ifelse(is.null(setting$thin), 1, setting$thin)
    setting$start.from <- ifelse(is.null(setting$start.from), 1, setting$start.from)
    setting$remove     <- ifelse(is.null(setting$remove), FALSE, setting$remove)
    setting$restart    <- ifelse(is.null(setting$restart), TRUE, setting$restart)

    if (verbose) {
      cat("Get setting from setting argument: add = "); cat(setting$add)
      cat(", verbose = "); cat(setting$verbose)
      cat(", rp =  "); cat(setting$rp);
      cat(", thin =  "); cat(setting$thin);
      cat(", start from iteration "); cat(setting$start.from)
      cat(" remove these iterations "); cat(setting$remove); cat("\n")
    }
    setting_in <- setting
  }

  ## Choose initialise_data or initialise_hyper
  if (is.null(samples) & is.list(data) & is.null(pp.prior)) {
    if (verbose) cat("Initialise non-hierarchical multiple-subject samples using data\n")
    out <- vector(mode="list", length=length(data))
    for(i in 1:length(data))
    {
      out[[i]] <- initialise_data(
        nmc         = nmc,            pList      = p.prior,
        data        = data[[i]],      samples    = NULL,
        theta1      = theta1,         startPList = start.prior,
        setting     = setting_in)
      attr(out[[i]]$theta, "dimnames") <- list(NULL, out[[i]]$p.names, NULL)
    }
    names(out) <- 1:length(data)
  } else if (is.null(samples) & is.list(data) & !is.null(pp.prior)) {
    if (verbose) cat("Initialise hierarchical samples for multiple participants using data\n")

    ## If data was found, check the setting of data level
    if ( !is.null(data) && !is.list(data) ) { stop("data must be a list") }

    ## If pp.prior was found, check the setting of hierarchical level
    ## [[1]] and [[2]] are location and scale, respectively
    if (!is.null(pp.prior)) {
      if (!is.list(pp.prior))    { stop("pp.prior must be a list") }
      if (length(pp.prior) != 2) { stop("pp.prior size must be 2") }
      if ( length(p.prior) < length(pp.prior[[1]]) ) {
        stop("Location hyperprior list cannot have more elements than p.prior")
      }
      if ( length(pp.prior[[1]]) < length(pp.prior[[2]]) ) {
        cat("Location hyperprior must have as many or more elements "); cat("\n")
        stop("than scale hyperprior.\n")
      }

      ## Check if a location prior has a corresponding scale prior
      ## If it does not, assume that the user wants to fixed the scale prior
      has.hyper <- names(p.prior) %in% names(pp.prior[[2]])
      has.sigma <- names(pp.prior[[1]]) %in% names(pp.prior[[2]])
      fixed     <- names(pp.prior[[1]])[ !has.sigma ]
      pp.names  <- lapply(pp.prior, names)

      ## Check Subject parameter has a hyper sigma
      ## Presuming has.hyper is for location, so correct Andrew's pp.names[[2]] to
      ## [[1]]
      if ( !all(names(p.prior)[has.hyper] %in% pp.names[[1]][has.sigma]) )
        stop("pp.prior location parameters not compatible with p.prior")
      if ( !all(names(p.prior)[has.hyper]==pp.names[[2]]) )
        stop("pp.prior scale parameters not compatible with p.prior")

      ## check hstart.prior
      if (!is.null(hstart.prior)) {
        if (!is.list(hstart.prior)) { stop("hstart.prior must be a list") }
        if (!all(sort(names(pp.prior[[1]]))==sort(names(hstart.prior[[1]]))))
          stop("Incompatible location hstart.prior")
        if (!all(sort(names(pp.prior[[2]]))==sort(names(hstart.prior[[2]]))))
          stop("Incompatible scale hstart.prior")
      }

      ## Check if constant parameters in the model match "fixed", stop
      bad <- fixed[!(fixed %in% names(attr(attr(data[[1]], "model"), "constants")))]
      if ( length(bad) > 0 ) {
        cat("Fixed hyper parameters not constants at data level: ");
        cat(bad); cat("\n"); stop("Check pp.prior")
      }

      ## Remake a hyperprior for scale parameter, when some are missing
      ## This guarantees location and scale prior are even and same length
      ## as p.prior; enable initialise.cpp to run faster
      newSca <- list()
      for(i in 1:length(p.prior)) {
        if (has.sigma[i]) {
          for (j in 1:length(pp.prior[[2]]))
          {
            if(names(p.prior[i]) == names(pp.prior[[2]][j]))
            {
              newSca[i] <- pp.prior[[2]][j]
            }
          }
        } else {
          p        <- c(0, 0, -Inf, Inf, 1)
          names(p) <- c("mean","sd","lower","upper", "log")
          p        <- as.list(p)
          attr(p,"dist")    <- "constant"
          attr(p,"untrans") <- "identity"
          newSca[[i]]       <- p
        }
      }
      names(newSca) <- names(p.prior)
      pp.prior <- list(pp.prior[[1]], newSca)

    }
    if(!is.null(phi1) & !is.list(theta1)) {
      stop("If phi1 specified, theta1 must be a list of thetas for each subject")
    }
    if (!is.null(hstart.prior)) {
      newSca <- list()
      for(i in 1:length(p.prior)) {
        if (has.sigma[i]) {
          for (j in 1:length(hstart.prior[[2]]))
          {
            if(names(p.prior[i]) == names(hstart.prior[[2]][j]))
            {
              newSca[i] <- hstart.prior[[2]][j]
            }
          }
        } else {
          p <- c(0, 0, -Inf, Inf, 1)
          names(p) <- c("mean","sd","lower","upper", "log")
          p <- as.list(p)
          attr(p,"dist") <- "constant"
          attr(p,"untrans") <- "identity"
          newSca[[i]] <- p
        }
      }
      names(newSca) <- names(p.prior)
      hstart.prior <- list(hstart.prior[[1]], newSca)

    }

    out <- initialise_hyper(
      nmc         = nmc,    pList      = p.prior,
      data        = data,   samples    = samples,
      theta1      = theta1, startPList = start.prior,
      phi1        = phi1,   ppList     = pp.prior,
      hStartPList = hstart.prior,
      setting     = setting_in )

    for(i in 1:length(out))
    {
      attr(out[[i]]$theta, "dimnames") <- list(NULL, out[[i]]$p.names, NULL)
    }
    names(out) <- 1:length(samples)

  } else if ( is.null(data) & is.null(pp.prior) & names(samples)[1]!="theta" ) {
    if (verbose) cat("Initialise non-hierarchical multiple-subject using samples \n")
    out <- vector(mode="list", length=length(samples))
    for(i in 1:length(samples))
    {
      out[[i]] <- initialise_data(
        nmc         = nmc,            pList      = samples[[i]]$p.prior,
        data        = NULL,           samples    = samples[[i]],
        theta1      = theta1,         startPList = start.prior,
        setting     = setting_in)
      attr(out[[i]]$theta, "dimnames") <- list(NULL, out[[i]]$p.names, NULL)
    }
    names(out) <- 1:length(samples)
  } else  {
    if (verbose) cat("Initialise hierarchical samples for multiple participants using samples\n")
    out <- initialise_hyper(
      nmc         = nmc,    pList      = samples[[1]]$p.prior,
      data        = NULL,   samples    = samples,
      theta1      = theta1, startPList = start.prior,
      phi1        = phi1,   ppList     = attr(samples, "hyper")$pp.prior,
      hStartPList = hstart.prior,
      setting     = setting_in )
  }

  cat("\n")
  return(out)

}


#' @importFrom stats runif
h.migrate <- function(use.phi,use.logprior,use.loglike,p.prior,ps,rp,pp.prior,
                      has.sigma,has.hyper,is.constant)
  # DEMCMC migrate set, all chains, hyper level
{

  # Which pars are not constants?
  de=list(!unlist(is.constant[[1]]),!unlist(is.constant[[2]]))

  n.chains <- dim(use.phi[[1]])[1]
  pars <- dim(use.phi[[1]])[2]
  lnum1 <- sample(c(1:n.chains),1)  # determine how many groups to work with
  lnum2 <- sort(sample(c(1:n.chains),lnum1,replace=F))  # get groups
  phiset <- list( # initialize
    matrix(NA,lnum1,pars,dimnames=list(NULL,dimnames(use.phi[[1]])[[2]])),
    matrix(NA,lnum1,pars,dimnames=list(NULL,dimnames(use.phi[[1]])[[2]]))
  )

  propset.logprior <- propset.loglike <- numeric(lnum1)
  currentset <- propset <- propw <- currw <- numeric(lnum1)

  index=numeric(lnum1)
  for (i in 1:lnum1) {
    index[i] <- sample(1:n.chains,1,replace=F)
    # create a set of these particles to swap
    phiset[[1]][i,] <- use.phi[[1]][lnum2[i],]
    phiset[[2]][i,] <- use.phi[[2]][lnum2[i],]

    # perturb non-constant parameters
    phiset[[1]][i,de[[1]]] <- phiset[[1]][i,de[[1]]] + runif(1,-rp,rp)
    phiset[[2]][i,de[[2]]] <- phiset[[2]][i,de[[2]]] + runif(1,-rp,rp)

    propset.logprior[i] <- h.summed.log.prior (
      pp=list(phiset[[1]][i,has.sigma],phiset[[2]][i,has.sigma]),pp.prior=pp.prior)
    propset.loglike[i]  <- sum(h.log.likelihood.dmc(ps=ps[i,,has.hyper],
      pp=list(phiset[[1]][i,has.sigma],phiset[[2]][i,has.sigma]),p.prior=p.prior[has.hyper]))

    propset[i] <- propset.logprior[i] + propset.loglike[i]
    propw[i] <- propset[i]
    currentset[i] <- use.logprior[lnum2[i]] + use.loglike[lnum2[i]]

  }
  currw <- currentset

  mh <- exp(propw[lnum1] - currw[1])
  if ( !is.na(mh) && (runif(1) < mh) ) {
    use.phi[[1]][lnum2[1],] <- phiset[[1]][lnum1,]	# swap the 1st with last
    use.phi[[2]][lnum2[1],] <- phiset[[2]][lnum1,]  #  (creating a circle)
    use.logprior[lnum2[1]] <- propset.logprior[lnum1]
    use.loglike[lnum2[1]] <- propset.loglike[lnum1]
  }
  if ( lnum1!=1 ) {										# make sure we are not done yet
    for(i in 1:(lnum1-1)){
      mh <- exp(propw[i] - currw[i+1])
      if ( !is.na(mh) && (runif(1) < mh) ) {
        use.phi[[1]][lnum2[i+1],] <- phiset[[1]][i,]
        use.phi[[2]][lnum2[i+1],] <- phiset[[2]][i,]
        use.logprior[lnum2[i+1]] <- propset.logprior[i]
        use.loglike[lnum2[i+1]] <- propset.loglike[i]
      }
    }
  }
  cbind(use.logprior,use.loglike,use.phi[[1]],use.phi[[2]])
}

#' @importFrom stats runif
blocked.h.crossover <- function(k,blocks,n.pars,use.phi,use.logprior,
  use.loglike, p.prior,ps,pp.prior,rp,has.sigma,has.hyper,is.constant,
  force=FALSE,gamma.mult=2.38,h.gamma.mult=2.38)
  # hyper level crossover for hypers with a sigma, as a series of blocks
{

  h.crossover <- function(k,pars,use.phi,use.logprior,use.loglike,p.prior,ps,
    is.constant,pp.prior,rp,has.sigma,has.hyper,h.gamma.mult=2.38,force=FALSE)
    # DEMCMC crossover update of one chain, phi level
  {

    # Which pars are not constants?
    de=list(!unlist(is.constant[[1]][names(is.constant[[1]])[pars]]),
            !unlist(is.constant[[2]][names(is.constant[[1]])[pars]]))
    if ( all(!unlist(de)) ) # all are constants
      return(c(use.logprior[k],use.loglike[k],use.phi[[1]][k,],use.phi[[2]][k,]))

    # step size
    if ( is.na(h.gamma.mult) ) hgamma <- runif(1,0.5,1) else
      hgamma <-  h.gamma.mult/sqrt(2*length(pars)*2) # extra *2 as p1 and p2

    # Update use.loglike for new ps
    phi <- list(use.phi[[1]][k,],use.phi[[2]][k,])
    use.loglike[k] <- sum(h.log.likelihood.dmc(ps=ps[k,,has.hyper],
      pp=list(phi[[1]][has.sigma],phi[[2]][has.sigma]),p.prior=p.prior[has.hyper]))

    # Calcualte use.post
    use.post <-  use.logprior[k] + use.loglike[k]

    # DE step

    # pick two other chains
    index <- sample(c(1:dim(use.phi[[1]])[1])[-k],2,replace=F)

    # update mu
    phi[[1]][pars[de[[1]]]] <- use.phi[[1]][k,pars[de[[1]]]] + runif(1,-rp,rp) +
      hgamma*(use.phi[[1]][index[1],pars[de[[1]]]]-use.phi[[1]][index[2],pars[de[[1]]]])

    # update sigma
    phi[[2]][pars[de[[2]]]] <- use.phi[[2]][k,pars[de[[2]]]] + runif(1,-rp,rp) +
      hgamma*(use.phi[[2]][index[1],pars[de[[2]]]]-use.phi[[2]][index[2],pars[de[[2]]]])

    # Get new post
    logprior <- h.summed.log.prior(pp.prior=pp.prior,
      pp=list(phi[[1]][has.sigma],phi[[2]][has.sigma])) # only for has.sigma
    loglike <- sum(h.log.likelihood.dmc(ps=ps[k,,has.hyper],p.prior=p.prior[has.hyper],
                   pp=list(phi[[1]][has.sigma],phi[[2]][has.sigma])))
    post <- logprior + loglike

    # Metropolis step
    epup <- exp(post-use.post)
    if ( force || (!is.na(epup)  && (runif(1) < epup)) ) {
      use.phi[[1]][k,pars] <- phi[[1]][pars]
      use.phi[[2]][k,pars] <- phi[[2]][pars]
      use.logprior[k] <- logprior
      use.loglike[k] <- loglike
    }

    c(use.logprior[k],use.loglike[k],use.phi[[1]][k,],use.phi[[2]][k,])
  }


  temp <- h.crossover(k,pars=blocks[[1]],use.phi=use.phi,force=force,
                      use.logprior=use.logprior,use.loglike=use.loglike,
                      p.prior=p.prior,ps=ps,pp.prior=pp.prior,rp=rp,
                      h.gamma.mult=h.gamma.mult,is.constant=is.constant,
                      has.sigma=has.sigma,has.hyper=has.hyper)
  if ( length(blocks)>1 ) for ( b in 2:length(blocks) ) {
    use.logprior[k]  <- temp[1]
    use.loglike[k]   <- temp[2]
    use.phi[[1]][k,] <- temp[3:(n.pars+2)]
    use.phi[[2]][k,] <- temp[(n.pars+3):(2*n.pars+2)]
    temp <- h.crossover(k,pars=blocks[[b]],use.phi=use.phi,force=force,
                        use.logprior=use.logprior,use.loglike=use.loglike,
                        p.prior=p.prior,ps=ps,pp.prior=pp.prior,rp=rp,
                        h.gamma.mult=h.gamma.mult,is.constant=is.constant,
                        has.sigma=has.sigma,has.hyper=has.hyper)
  }
  temp
}

#' @importFrom stats runif
crossover.h <- function(k,pars,use.theta,use.logprior,use.loglike,p.priors,data,
                        rp,gamma.mult=2.38,consts=NULL,pp.priors=NULL,
  force=FALSE)
  # Data level crossover wity different priors for each chain, p.priors is
  # a list
{
  # step size
  if (is.na(gamma.mult)) gamma <- runif(1,0.5,1) else
    gamma <-  gamma.mult/sqrt(2*length(pars))

  # DE step
  # pick two other chains
  index <- sample(c(1:dim(use.theta)[1])[-k],2,replace=F)

  # Update theta
  theta <- use.theta[k,]
  names(theta)=names(attr(attributes(data)$model,"p.vector"))
  theta[pars] <-
    use.theta[k,pars] +
    gamma*(use.theta[index[1],pars]-use.theta[index[2],pars]) +
    runif(1,-rp,rp)

  # Combine use to get old post
  summed.use.post <-  use.loglike[k] + use.logprior[k]

  # Get prior for new point
  log.prior <- summed.log.prior (p.vector=theta,p.prior=p.priors[[k]],
    p.names=names(attr(attributes(data)$model,"p.vector")))

  if ( !is.null(consts) ) { # !has.sigma hyper to incorporte
    log.prior <- log.prior + summed.log.prior(
      consts[k,],pp.priors,p.names=dimnames(consts)[[2]])
    attr(attr(data,"model"),"all.par")[dimnames(consts)[[2]]] <- consts[k,]
  }
  loglike  <- sum(log.likelihood (p.vector=theta,data=data))

  # post for new point
  post <- log.prior + loglike

  # Metropolis step
  epup <- exp(post-summed.use.post)
  if (force || (!is.na(epup) && (runif(1) < epup)) ) {
    use.theta[k,]   <- theta
    use.logprior[k] <- log.prior
    use.loglike[k]  <- loglike
  }

  c(use.logprior[k], use.loglike[k], use.theta[k,])
}


migrate.h <- function(use.theta,use.logprior,use.loglike,
                      p.priors,data,rp,consts=NULL,pp.priors=NULL)
  # Data level migrate with different priors for each chain, p.priors is a list
{
  n.chains <- dim(use.theta)[1]
  pars  <- dim(use.theta)[2]
  n.groups <- sample(c(1:n.chains),1)  # determine how many groups to work with
  groups <- sort(sample(c(1:n.chains),n.groups,replace=F))	# get groups

  thetaset <- matrix(NA,n.groups,pars)									   # initialize
  currentset.logprior <- propset.logprior <- numeric(n.groups)
  propw.logprior      <- currw.logprior   <- numeric(n.groups)
  currentset.loglike  <- propset.loglike  <- numeric(n.groups)
  propw.loglike       <- currw.loglike    <- numeric(n.groups)

  for (i in 1:n.groups) {
    # create a set of these particles to swap
    thetaset[i,] <- use.theta[groups[i],] + runif(1,-rp,rp)

    currentset.logprior[i] <- use.logprior[groups[i]]
    currentset.loglike[i]  <- use.loglike[groups[i]]

    # Proposed new prior
    propset.logprior[i] <- summed.log.prior (p.vector=thetaset[i,],
      p.prior=p.priors[[i]],p.names=names(attr(attributes(data)$model,"p.vector")))

    if ( !is.null(consts) ) { # !has.sigma hyper to incorporte
      propset.logprior[i] <- propset.logprior[i] + summed.log.prior(
        consts[i,],pp.priors,p.names=dimnames(consts)[[2]])
      attr(attr(data,"model"),"all.par")[dimnames(consts)[[2]]] <- consts[i,]
    }

    propset.loglike[i]  <- sum(log.likelihood (p.vector=thetaset[i,],data=data))

    propw.logprior[i]      <- propset.logprior[i]
    propw.loglike[i]       <- propset.loglike[i]

    currw.logprior[i]      <- currentset.logprior[i]
    currw.loglike[i]       <- currentset.loglike[i]

  }

  epup <- exp( (propset.loglike[n.groups] + propset.logprior[n.groups]) -
                      (currentset.loglike[1] + currentset.logprior[1]))
  if (!is.na(epup) && (runif(1) < epup)) {
    use.theta[groups[1],]   <- thetaset[n.groups,]	# swap the 1st with last
    use.logprior[groups[1]] <- propset.logprior[n.groups]
    use.loglike[groups[1]]  <- propset.loglike[n.groups]
  }

  if ( n.groups!=1 ) {										# make sure we are not done yet
    for(i in 1:(n.groups-1)) {
      epup <- exp((propset.loglike[i] + propset.logprior[i]) -
                  (currentset.loglike[i+1] + currentset.logprior[i+1]))
      if(!is.na(epup) && (runif(1) < epup)) {
        use.theta[groups[i+1],] <- thetaset[i,]
        use.logprior[groups[i+1]] <- propset.logprior[i]
        use.loglike[groups[i+1]] <- propset.loglike[i]
      }
    }
  }

  cbind (use.logprior, use.loglike, use.theta)
}


#' Fit an EAM with Multiple Participants
#'
#' \code{h.run.dmc} calls \code{run_data} in the case of fixed-effect model, or
#' \code{run_hyper} in the case of random-effect model. At the moment, it fits
#' either fixed-effect or random-effect model by looking for the
#' \code{hyper} attribute in a DMC sample/object.
#'
#' @param samples a DMC samples, usually generated by calling
#' \code{h.samples.dmc}.
#' @param report the progress report. Default is every 100 iterations report
#' once
#' @param cores a switch for parallel computation.  In the case of fixed-effect
#' model, \code{h.run.dmc} calls \pkg{snow} (via \pkg{snowfall}) or
#' \pkg{parallel}. When setting for example \code{cores=4}, \code{h.run.dmc}
#' requests 4 parallel R processes. In the case of random-effect model,
#' \code{h.run.dmc} calls OpenMP library, which will automatically decide the
#' number of cores to use, so setting core number higher than 1 has the same
#' effect.  Use OpenMP only when the data set is large (e.g., more than 1000
#' trial per condition).
#' @param blocks a DEMCMC parameter for blocking. Not implemented yet. Default
#' is NA
#' @param p.migrate the probability of applying migration sampler. Set it
#' greater than will trigger the migration sampler. For example
#' \code{p.migrate=0.05} will use migration sampler about 5\% chance.
#' @param h.p.migrate similar to p.migrate for hyper level
#' @param force.hyper This argument has no funciton. It is here only for
#' compatibility (with DMC) reason.
#' @param force.data This argument has no funciton. It is here only for
#' compatibility (with DMC) reason.
#' @param gamma.mult a DE-MCMC tuning parameter, affecting the size of jump.
#' Default is 2.38.
#' @param h.gamma.mult Similar with gamm.mult.  This argument has no funciton.
#' It is here only for compatibility reason.
#' @param setting a list to store setting arguments, such as p.migrate,
#' gamma.mult, etc, sending to C++ function. The user should not use this
#' argument
#' @param verbose a swtich to see debugging information
#' @keywords h.run.dmc
#' @import snowfall
#' @import rtdists
#' @import pracma
#' @import statmod
#' @importFrom parallel mclapply
#' @export
#' @examples
#' m1 <- model.dmc(
#' p.map     = list(a="1",v="1",z="1",d="1",sz="1",sv="1",t0="1",st0="1"),
#' match.map = list(M=list(s1="r1", s2="r2")),
#' factors   = list(S=c("s1", "s2")),
#' constants = c(st0=0, d=0),
#' responses = c("r1","r2"),
#' type      = "rd")
#'
#' ## Population distribution
#' pop.prior <- prior.p.dmc(
#'   dists = rep("tnorm", 6),
#'   p1    = c(a=2,  v=2.5, z=.5, sz=.3, sv=1,  t0=.3),
#'   p2    = c(a=.5, v=.5,  z=.1, sz=.1, sv=.3, t0=.05),
#'   lower = c(0,-5, 0, 0, 0, 0),
#'   upper = c(5, 7, 2, 2, 2, 2))
#'
#' dat <- h.simulate.dmc(m1, nsim=30, p.prior=pop.prior, ns=8)
#' mdi <- data.model.dmc(dat, m1)
#' p.prior  <- prior.p.dmc(
#'   dists = rep("tnorm", 6),
#'   p1    = c(a=2,  v=2.5, z=.5, sz=.3, sv=1,  t0=.3),
#'   p2    = c(a=.5, v=.5,  z=.1, sz=.1, sv=.3, t0=.05) * 5,
#'   lower = c(0,-5, 0, 0, 0, 0),
#'   upper = c(5, 7, 2, 2, 2, 2))
#'
#' ## Fixed-effect model
#' samplesInit <- h.samples.dmc(nmc=20, p.prior=p.prior, data=mdi, thin=1)
#' samples0    <- h.run.dmc(samples=samplesInit, report=10)
#'
#' ## Make a hyper-prior list
#' mu.prior <- prior.p.dmc(
#'   dists = rep("tnorm", 6),
#'   p1    = c(a=2,  v=2.5, z=.5, sz=.3, sv=1,  t0=.3),
#'   p2    = c(a=.5, v=.5,  z=.1, sz=.1, sv=.3, t0=.05) * 5,
#'   lower = c(0,-5, 0, 0, 0, 0),
#'   upper = c(5, 7, 2, 2, 2, 2))
#'
#' sigma.prior <- prior.p.dmc(
#'   dists = rep("beta", 6),
#'   p1    = c(a=1, v=1, z=1, sz=1, sv=1, t0=1),
#'   p2    = c(1,1,1,1,1,1),
#'   upper = c(2,2,2,2,2,2))
#'
#' pp.prior <- list(mu.prior, sigma.prior)
#'
#' ## Random-effect model
#' hsamplesInit <- h.samples.dmc(nmc=20, p.prior=p.prior, pp.prior=pp.prior,
#'   data=mdi, thin=1)
#' hsamples0 <- h.run.dmc(samples=hsamplesInit, report=10, p.migrate=.05,
#'   h.p.migrate=.05)
h.run.dmc <- function(samples, report=100, cores=1, blocks=NA, p.migrate=0,
  h.p.migrate=0, force.hyper=NA, force.data=NA, gamma.mult=2.38,
  h.gamma.mult=NA, setting=NULL, verbose=FALSE)
{
  os <- get_os()

  ## Check 2
  if( is.null(setting) )
  {
    setting_in <- list(p.migrate=p.migrate, h.p.migrate=h.p.migrate,
      gamma.mult=gamma.mult, h.gamma.mult=2.38, report=report, cores=cores)
    if(verbose)
    {
      cat("R uses default setting: p.migrate="); cat(p.migrate);
      cat(", h.p.migrate=0, gamma.mult="); cat(gamma.mult);
      cat(", h.gamma.mult=", h.p.migrate, "\n"); cat("run will use "); cat(cores);
      cat(" CPU and report progress every "); cat(report); cat(" step \n")
    }

  } else {
    isSet <- names(setting) %in% c("p.migrate","h.p.migrate","gamma.mult",
      "h.gamma.mult","cores","report")
    setting$p.migrate    <- ifelse(isSet[1], setting$p.migrate,    "0.05")
    setting$h.p.migrate  <- ifelse(isSet[2], setting$h.p.migrate,  "0.05")
    setting$gamma.mult   <- ifelse(isSet[3], setting$gamma.mult,   "2.38")
    setting$h.gamma.mult <- ifelse(isSet[4], setting$h.gamma.mult, "2.38")
    setting$cores        <- ifelse(isSet[5], setting$cores, "1")
    setting$report       <- ifelse(isSet[6], setting$report, "100")
    if(verbose)
    {
      cat("R: Using user: p.migrate = "); cat(setting$p.migrate)
      cat(" h.p.migrate = ");  cat(setting$h.p.migrate)
      cat(" gamma.mult = ");   cat(setting$gamma.mult)
      cat(" h.gamma.mult = "); cat(setting$h.gamma.mult)
      cat(".\nrun will use "); cat(setting$cores)
      cat(" CPU(s) and report every "); cat(setting$report);
      cat(" step.\n");
    }
    setting_in <- setting
  }

  ## ---------------------------------------------------------------------- */
  ## Gyakuzuki: not a very good way to overload function. However Rcpp
  ## has yet had a very conventient overloading method. Either overload
  ## inside C++ or overload with S3/S4
  ## ---------------------------------------------------------------------- */
  ## If hyper attribute is not set, it must be data level
  if ( is.null( attr(samples, ("hyper"))) ) ## fixed-effect
  {
    ## If the 1st name is theta, assume we get only one participant
    ## Otherwise, mulitple participants
    if (verbose) { cat("Run fixed-effect: ") }
    sNames <- names(samples)    ## samples or subjects' names
    if (sNames[1] == "theta" )  ## fixed-effect for one subject
    {
      if (debug==TRUE) {
        out <- run_data(samples, setting_in, debug=TRUE) ;
      } else if (cores > 1) {
        cat("Run in parallel, using Open MP. Set cores=1 to turn off parallel\n")
        out <- run_data_parallel(samples, setting_in) ;
      } else {
        out <- run_data(samples, setting_in) ;
      }
      dimnames(out$theta) <- list(NULL, out$p.names, NULL)
    } else {                    ## fixed-effect for multiple subjects
      if (os == "windows" & cores > 1) {
        snowfall::sfInit(parallel=TRUE, cpus=cores, type="SOCK")
        snowfall::sfClusterSetupRNG()
        snowfall::sfLibrary(ggdmc)
        snowfall::sfLibrary(rtdists)
        snowfall::sfLibrary(pracma)
        snowfall::sfLibrary(statmod)
        snowfall::sfExportAll()

        out <- snowfall::sfLapply(samples, run_data, setting_in);
        snowfall::sfStop()
        for(i in 1:length(out)) {
          dimnames(out[[i]]$theta) <- list(NULL, out[[i]]$p.names, NULL)
        }

      } else if (cores > 1) {
        ## require(parallel, quietly=TRUE)
        if (verbose) { cat("Use parallel package\n") }

        out <- parallel::mclapply(samples, run_data, setting_in, mc.cores=cores)
        for(i in 1:length(out)) {
          dimnames(out[[i]]$theta) <- list(NULL, out[[i]]$p.names, NULL)
        }
      } else {
        if (verbose) { cat("serially using lapply\n") }
        out <- lapply(samples, run_data, setting_in)
        for(i in 1:length(out)) {
          dimnames(out[[i]]$theta) <- list(NULL, out[[i]]$p.names, NULL)
        }
      }
    }
    out <- dmc.list(out)
    ## random-effect
  } else if (!is.null(attr(samples, "hyper"))) {
    if (verbose) { cat("Run random-effect fit\n") }
    if (cores > 1) {
      out <- run_hyper_parallel(samples, setting_in)
    } else {
      out <- run_hyper(samples, setting_in)
    }

    for(i in 1:length(out)) {
      dimnames(out[[i]]$theta) <- list(NULL, out[[i]]$p.names, NULL)
    }

    hyper0 <- attr(out,"hyper")
    for(j in 1:length(hyper0$phi)) {
      dimnames(hyper0$phi[[j]]) <- list(NULL, out[[1]]$p.names, NULL)
    }
    attr(out, "hyper") <- hyper0
    names(out) <- 1:length(out)
    out <- hyper(out)

  } else {
    stop("Unexpected condition occurred")
  }

  cat("\n")
  return(out)
}


h.run.unstuck.dmc <- function(samples,nmc=NA, report=10,cores=1,
                              cut=10,nbad=2,max.try=10,gamma.mult=2.38,
                              h.gamma.mult=NA)
  # repeats sampling until <= nbad stuck chains as defined by cut
{
  if (!is.null(samples$theta))
    stop("For a single subject use h.run.unstuck.dmc")
  if (any(is.na(samples[[1]]$theta[,,2])))
    samples <- h.run.dmc(samples=samples, report=report,cores=cores,
                         gamma.mult=gamma.mult,h.gamma.mult=h.gamma.mult)
  n.chains <- unlist(lapply(samples,function(x){x$n.chains}))
  if ( is.na(nmc) ) nmc <- samples[[1]]$nmc
  try.num <- 1
  repeat {
    stucks <- lapply(samples,pick.stuck.dmc)
    ok <- unlist(lapply(stucks,function(x){length(x)<=nbad}))
    if (all(ok)) break
    samples[!ok] <- h.run.dmc(h.samples.dmc(nmc=nmc,samples=samples[!ok]),
                              report=report,cores=cores,
                              gamma.mult=gamma.mult,h.gamma.mult=h.gamma.mult)
    if (try.num >= max.try) break
    try.num <- try.num + 1
  }
  for (i in 1:length(samples)) {
    samples[[i]] <-
      unstick.dmc(samples[[i]],pick.stuck.dmc(samples[[i]],cut=cut))
    cat(paste("Subject",i,": "))
    print(pick.stuck.dmc(samples[[i]],verbose=TRUE))
  }
  tmp <- rbind(Originally=n.chains,
               Removed=n.chains-unlist(lapply(samples,function(x){x$n.chains})))
  dimnames(tmp)[[2]] <- names(samples)
  cat("Removed chains\n"); print(tmp)
  samples
}


h.run.converge.dmc <- function(samples,hyper=FALSE,nmc=NA,report=10,cores=1,
                cut=1.1,max.try=10,transform=TRUE,autoburnin=FALSE,blocks=NA,
                gamma.mult=2.38,h.gamma.mult=NA)
  # repeats sampling until gelman.diag Multivariate psrf < cut for all subjects
{

  do.gd <- function(samples,transform)
    gelman.diag(theta.as.mcmc.list(samples),
                transform=transform,autoburnin=autoburnin)$mpsrf

  if (!is.null(samples$theta))
    stop("For a single subject use h.run.unstuck.dmc")
  if (any(is.na(samples[[1]]$theta[,,2])))
    samples <- h.run.dmc(samples=samples,report=report,cores=cores,
                         blocks=blocks,
                         gamma.mult=gamma.mult,h.gamma.mult=h.gamma.mult)
  if ( is.na(nmc) ) nmc <- samples[[1]]$nmc
  try.num <- 1
  thin <- samples[[1]]$thin
  if ( hyper ) {
    gd <- gelman.diag(phi.as.mcmc.list(attr(samples,"hyper")),
                      transform=transform,autoburnin=autoburnin)$mpsrf
    repeat {
      ok <- gd<cut
      if ( ok) break
      print(paste("Multivariate psrf achieved =",round(gd,3)))
      samples <- h.run.dmc(h.samples.dmc(nmc=nmc,samples=samples,thin=thin),
                           report=report,cores=cores,blocks=blocks,
                           gamma.mult=gamma.mult,h.gamma.mult=h.gamma.mult)
      if (try.num >= max.try) break
      gd <- gelman.diag(phi.as.mcmc.list(attr(samples,"hyper")),
                        transform=transform,autoburnin=autoburnin)$mpsrf
      try.num <- try.num + 1
    }
  } else {
    gd <- unlist(lapply(samples,do.gd,transform=transform))
    repeat {
      ok <- gd<cut
      if ( all(ok)) break
      print(paste("Multivariate psrf above cutoff =",
                  paste(round(sort(gd[!ok]),3),collapse=" ")))
      new.samples <- h.run.dmc(h.samples.dmc(nmc=nmc,samples=samples[!ok],thin=thin),
                               report=report,cores=cores,blocks=blocks,
                               gamma.mult=gamma.mult,h.gamma.mult=h.gamma.mult)
      gd.new <- unlist(lapply(new.samples,do.gd,transform=transform))
      for (i in 1:length(gd.new)) if (gd.new[i]<gd[!ok][i])
        samples[!ok][i] <- new.samples[i]
      gd[!ok] <- unlist(lapply(samples[!ok],do.gd,transform=transform))
      if (try.num >= max.try) break
      try.num <- try.num + 1
    }
  }
  print(paste("Final multivariate psrf =",
              paste(round(sort(gd),3),collapse=" ")))
  samples
}




