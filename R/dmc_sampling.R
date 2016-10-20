# System functions for the DMC (Dynamic Models of Choice)
#    Functions specific to sampling and related operations
#    Usually user does not need to edit

######### PRIOR and POSTERIOR

dbeta_lu <- function(x,shape1,shape2,lower,upper,log=FALSE)
  # Used with beta prior
{
  if (log) {dbeta((x-lower)/(upper-lower),shape1,shape2,log=log)}
  else {dbeta((x-lower)/(upper-lower),shape1,shape2,log=log)/(upper-lower)}
}


rbeta_lu <- function(n,shape1,shape2,lower,upper)
  # Used with beta prior
{
  lower + rbeta(n,shape1,shape2)*(upper-lower)
}

dgamma_l <- function(x,shape,scale,lower,log=FALSE)
  # Used with gamma prior
{
  dgamma(x-lower,shape=shape,scale=scale,log=log)
}

rgamma_l <- function(n,shape,scale,lower)
  # Used with gamma prior
{
  lower + rgamma(n,shape=shape,scale=scale)
}

dlnorm_l <- function(x,meanlog,sdlog,lower,log=FALSE)
  # Used with lognormal prior
{
  dlnorm(x-lower,meanlog,sdlog,log=log)
}

rlnorm_l <- function(n,meanlog,sdlog,lower)
  # Used with lognormal prior
{
  lower + rlnorm(n,meanlog,sdlog)
}

dconstant <- function(x,constant,log=FALSE)
  # Used with constant prior
{
  if (log) rep(0,length(constant)) else
           rep(1,length(constant))
}


rconstant <- function(n,constant)
  # Used with constant prior
{
  rep(constant,n)
}

#' Makes a list of prior distribution parameters.
#'
#' \code{prior.p.dmc} creates a list of prior distribution an array object ("model") with a set of attributes
#' specifying a particular model and parameterization. Call \pkg{coda} to
#' summarise the model parameters in a DMC samples with multiple participants
#' at the hyper level.
#'
#' @param p1 the values of location parameters for each prior distribution, set
#' as a double vector
#' @param p2 ditto for scale parameter vector
#' @param lower lower support boundary
#' @param upper upper support boundary
#' @param dists indicate which prior distribution, e.g., uniform, beta etc.
#' @param untrans whether do log transformation or not. Default is identity,
#' namely not to transform
#' @param dist.types allowed prior distributions in current version of DMC
#' @keywords prior.p.dmc
#' @export
prior.p.dmc <- function(p1,p2,
                        lower=rep(NA,length(p1)),
                        upper=rep(NA,length(p1)),
                        dists=rep("tnorm",length(p1)),
                        untrans=rep("identity",length(p1)),
                        dist.types=c("tnorm","beta","gamma","lnorm","constant"))
{

  if (length(p2)==1)
    p2 <- rep(p2,length(p1))
  if ( length(p1)!=length(p2) )
    stop("p1 and p2 must have the same length")
  if ( length(p1)!=length(lower) )
    stop("p1 and lower must have the same length")
  if ( length(p1)!=length(upper) )
    stop("p1 and upper must have the same length")
  both.not.na <- !is.na(upper) & !is.na(lower)
  if ( any(upper[both.not.na]<=lower[both.not.na]) )
    stop("All elements of upper must be greater than lower")
  if ( length(p1)!=length(dists) )
    stop("p1 and dists must have the same length")
  if ( !all(dists %in% dist.types) )
    stop(paste("Unsupported distribution, allowable types are:",
               paste(dist.types,collapse=", ")))
  name.untrans <- length(untrans) != length(p1)
  if (name.untrans & (is.null(names(untrans)) | is.null(names(p1))))
    stop("If untrans vector is not the same length as p1 it must have p1 names")
  if (!(all(names(untrans) %in% names(p1))))
    stop("untrans vector has names not in p1 names")

  prior <- vector(mode="list",length=length(p1))
  names(prior) <- names(p1)
  for ( i in 1:length(p1) ) {
    prior[[i]] <- switch(dists[i],
                         tnorm={
                           if (is.na(lower[i])) lower[i] <- -Inf
                           if (is.na(upper[i])) upper[i] <- Inf
                           p <- c(p1[i],p2[i],lower[i],upper[i])
                           names(p) <- c("mean","sd","lower","upper")
                           p <- as.list(p)
                           attr(p,"dist") <- "tnorm"
                           p
                         },
                         beta={
                           if (is.na(lower[i])) lower[i] <- 0
                           if (is.na(upper[i])) upper[i] <- 1
                           p <- c(p1[i],p2[i],lower[i],upper[i])
                           names(p) <- c("shape1","shape2","lower","upper")
                           p <- as.list(p)
                           attr(p,"dist") <- "beta_lu"
                           p
                         },
                         gamma={
                           if (is.na(lower[i])) lower[i] <- 0
                           p <- c(p1[i],p2[i],lower[i])
                           names(p) <- c("shape","scale","lower")
                           p <- as.list(p)
                           attr(p,"dist") <- "gamma_l"
                           p
                         },
                         lnorm={
                           if (is.na(lower[i])) lower[i] <- 0
                           p <- c(p1[i],p2[i],lower[i])
                           names(p) <- c("meanlog","sdlog","lower")
                           p <- as.list(p)
                           attr(p,"dist") <- "lnorm_l"
                           p
                         },
                         {
                           p <- p1[i]
                           names(p) <- c("constant")
                           p <- as.list(p)
                           attr(p,"dist") <- "constant"
                           p
                         }
    )
    prior[[i]]$log <- TRUE
    if (!name.untrans) attr(prior[[i]],"untrans") <- untrans[i] else
      if (is.na(untrans[names(p1)[i]]))
        attr(prior[[i]],"untrans") <- "identity" else
        attr(prior[[i]],"untrans") <- untrans[names(p1)[i]]
  }
  prior
}


log.prior.dmc <- function(p.vector,p.prior)
  # log of prior density for p.vector
{
  prior <- numeric(length(p.vector))
  names(prior) <- names(p.vector)
  for ( i in names(p.vector) ) {
    p.prior[[i]]$x <- p.vector[i]
    prior[i] <- do.call(paste("d",attr(p.prior[[i]],"dist"),sep=""),
                        p.prior[[i]])
  }
  prior
}


rprior.dmc <- function(p.prior,n=1)
  # sample from prior
{
  prior <- matrix(nrow=n,ncol=length(p.prior))
  dimnames(prior) <- list(NULL,names(p.prior))
  for ( i in 1:length(names(p.prior)) ) {
    p.prior[[i]]$n <- n
    prior[,i] <- do.call(paste("r",attr(p.prior[[i]],"dist"),sep=""),
                         p.prior[[i]][names(p.prior[[i]])!="log"])
  }
  prior
}


log.posterior.dmc <- function(p.vector,p.prior,data,min.like=1e-10)
  # Summed log posterior likelihood
{
  sum (log.likelihood(p.vector,data,min.like=min.like)) +
    summed.log.prior (p.vector, p.prior,
      p.names=names(attr(attributes(data)$model,"p.vector")))

}


log.likelihood <- function(p.vector, data, min.like=1e-10)
  ## Get log likelihood
  ## TODO likelihood.dmc to likelihood for ddm lba_B classes
{
  names(p.vector)=names(attr(attributes(data)$model,"p.vector"))
  suppressWarnings( log(likelihood.rd(p.vector,data,min.like=min.like) ) )
}


summed.log.prior <- function (p.vector, p.prior, p.names)
{
  names(p.vector)=p.names
  suppressWarnings( sum(log.prior.dmc(p.vector,p.prior)) )
}


assign.pp <- function(pp,p.prior)
  # Slot pp values into p.prior
{
  for (i in 1:length(p.prior))
    p.prior[[i]][1:2] <- c(pp[[1]][i],pp[[2]][i])
  p.prior
}

######### Sampling

#' Initialising a DMC samples
#'
#' A wrapper function for C++ \code{initialise_data}. \code{samples.dmc}
#' initialise a DMC sample for one participant. The user needs to minimally
#' enter three arguments: \code{nmc}, \code{p.prior} and \code{data} or
#' \code{nmc}, \code{p.prior} and \code{samples}. Note the current version
#' initialise only drift-diffusion model samples.
#'
#' \emph{initialise_data} generates random samples based on the probability
#' functions listed in prior parameter list (ie p.prior). This sets up a MCMC
#' sample (ie samples). \emph{initialise_hyper} generates a samples
#' hierarchically based on hyper-prior and prior parameter lists (ie p.prior
#' and pp.prior). For initialising multiple participants, please use
#' \code{h.samples.dmc}.
#'
#' @param nmc the number of MCMC iterations.
#' @param p.prior prior parameter list
#' @param data a model data instance.
#' @param thin thinning length. Default is 1.
#' @param samples a DMC posterior sample
#' @param theta1 a nChains x npar matrix
#' @param restart whether to restart a new samples
#' @param add whether to add more iterations on top of previous DMC sample
#' @param remove whether to remove some of the samples
#' @param start.prior a user indicated alternative prior parameter list
#' @param rp DE-MCMC tuning parameter. The user usually needs not set this
#' argument. Default 0.001
#' @param setting a list carrying thin, restart, add, remove, start.from,
#' start.prior, rp, and verbose setting. The user usually needs not set this
#' argument.
#' @param verbose whether to print debugging information
#' @export
#' @examples
#' m1 <- model.dmc(
#'     p.map     = list(a="1",v="F",z="1",d="1",sz="1",sv="1", t0="1",st0="1"),
#'     constants = c(st0=0,d=0),
#'     match.map = list(M=list(s1="r1",s2="r2")),
#'     factors   = list(S=c("s1","s2"), F=c("f1", "f2")),
#'     responses = c("r1","r2"),
#'     type      = "rd")
#'
#' ## m1 is "dmc" class
#' class(m1)
#' ## [1] "dmc"
#'
#' pVec <- c(a=1, v.f1=1, v.f2=1.5, z=0.5, sz=0.25, sv=0.2,t0=.15)
#' dat  <- simulate(m1, nsim=1e2, p.vector=pVec)
#' str(dat)
#' ## 'data.frame':	400 obs. of  4 variables:
#' ## $ S : Factor w/ 2 levels "s1","s2": 1 1 1 1 1 1 1 1 1 1 ...
#' ## $ F : Factor w/ 2 levels "f1","f2": 1 1 1 1 1 1 1 1 1 1 ...
#' ## $ R : Factor w/ 2 levels "r1","r2": 1 1 1 2 1 1 1 1 2 1 ...
#' ## $ RT: num  0.26 0.255 0.572 0.25 0.518 ...
#'
#' ## mdi1 is dmc as well as data.frame class
#' mdi1 <- data.model.dmc(dat, m1)
#' class(mdi1)
#' ## [1] "dmc"        "data.frame"
#'
#' p.prior <- prior.p.dmc(
#'    dists = rep("tnorm", 7),
#'    p1    = c(a=2,  v.f1=2.5, v.f2=1.25, z=.5, sz=.3, sv=1,  t0=.3),
#'    p2    = c(a=.5, v.f1=.5,  v.f2=.35,  z=.1, sz=.1, sv=.3, t0=.05),
#'    lower = c(0,-5, -5, 0, 0, 0, 0),
#'    upper = c(5, 7,  7, 2, 2, 2, 2))
#'
#' ## Set up a new DMC sample with 200 iteration. The default thinning length
#' ## is 1
#' samples0 <- samples.dmc(nmc=200, p.prior=p.prior, data=mdi1)
#' samples0$nmc
#' ## [1] 200
#'
#' ## Run a fixed-effect model with 5% chance of using migration sampler
#' samples0 <- run.dmc(samples0, p.migrate=.05)
#'
#' ## Add 200 more iteration on to sample0
#' samples1 <- samples.dmc(nmc=200, p.prior=p.prior, samples=samples0, add=TRUE)
#' ## samples1$nmc
#' ## [1] 400
samples.dmc <- function(nmc,
  p.prior      = NULL,  data         = NULL,
  thin         = 1,     samples      = NULL,
  theta1       = NULL,  restart      = TRUE,
  add          = FALSE, remove       = FALSE,
  start.prior  = NULL,  rp           = .001,
  setting      = NULL,  verbose      = FALSE)
{
  chainFamily <- 3
  if(is.null(samples) & is.null(data))
    { stop("Neither samples nor data was found!\n") }

  if( is.null(setting) & is.null(samples) ) {
    model    <- attr(data, "model")
    p.names  <- names(attr(model, "p.vector"))
    n.chains <- chainFamily*length(p.names) ## length(p.names) == npars
    if (verbose) {
      cat("Apply default setting: add=FALSE, verbose=FALSE, ")
      cat("restart=TRUE, rp=0.001, thin=1, n.chains=", n.chains)
      cat(", start from = 1\n")
    }
    setting_in <- list(add=add, verbose=verbose, rp=rp, restart=restart,
      thin=thin, start.from=1, n.chains=n.chains, remove=remove)
  } else if (is.null(setting) & !is.null(samples) ) {
    n.chains <- samples$n.chains
    sNames <- names(samples)
    if(sNames[1] != "theta") { stop("Use h.run.dmc for multiple participants")}

    if(verbose) {cat("Data: Get setting from samples and restart\n")}
    setting_in <- list(add=add, verbose=verbose,
      rp=samples$rp, thin=samples$thin, restart=restart,
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
  if (is.null(samples) & is.data.frame(data)) {
    if(verbose){cat("Initialise one-subject samples using data\n")}
    out <- initialise_data(
      nmc         = nmc,    pList      = p.prior,
      data        = data,   samples    = samples,
      theta1      = theta1, startPList = start.prior,
      setting     = setting_in)
    attr(out$theta, "dimnames") <- list(NULL, out$p.names, NULL)
  } else {
    if(verbose){cat("Initialise one-subject samples using samples\n")}
    out <- initialise_data(
      nmc         = nmc,    pList      = samples$p.prior,
      data        = NULL,   samples    = samples,
      theta1      = theta1, startPList = start.prior,
      setting     = setting_in)
    attr(out$theta, "dimnames") <- list(NULL, out$p.names, NULL)
  }

  cat("\n")
  return(out)
}


crossover <- function(k,pars,use.theta,use.logprior,use.loglike,p.prior,data,
                      rp,gamma.mult=2.38,force=FALSE)
  # DEMCMC crossover update of one chain, data level
{
  # step size
  if (is.na(gamma.mult)) gamma <- runif(1,0.5,1) else
    gamma <-  gamma.mult/sqrt(2*length(pars))
  # pick two other chains
  index <- sample(c(1:dim(use.theta)[1])[-k],2,replace=F)

  # DE step
  theta <- use.theta[k,]
  names(theta)=names(attr(attributes(data)$model,"p.vector"))

  theta[pars] <-
    use.theta[k,pars] +
    gamma*(use.theta[index[1],pars]-use.theta[index[2],pars]) +
    runif(1,-rp,rp)

  # Old post
  if (force) {
    use.logprior[k]<- summed.log.prior (p.vector=use.theta[k,],p.prior=p.prior,
      p.names=names(attr(attributes(data)$model,"p.vector")))
    use.loglike[k] <- sum(log.likelihood (p.vector=use.theta[k,],data=data))
  }
  summed.use.post <-  use.loglike[k] + use.logprior[k]

  # new post
  log.prior <- summed.log.prior (p.vector=theta,p.prior=p.prior,
    p.names=names(attr(attributes(data)$model,"p.vector")))
  loglike  <- sum(log.likelihood (p.vector=theta,data=data))
  post <- log.prior + loglike
  if ( is.na(post) ) post <- -Inf

  # Metropolis step
  if ( runif(1) < exp(post-summed.use.post) ) {
    use.theta[k,]   <- theta
    use.logprior[k] <- log.prior
    use.loglike[k]  <- loglike
  }

  c(use.logprior[k], use.loglike[k], use.theta[k,])
}


migrate <- function(use.theta,use.logprior,use.loglike,
                    p.prior,data,rp)
  # DEMCMC migrate set, all chains, data level
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

    propset.logprior[i] <- summed.log.prior (p.vector=thetaset[i,],
      p.prior=p.prior,p.names=names(attr(attributes(data)$model,"p.vector")))
    propset.loglike[i]  <- sum(log.likelihood (p.vector=thetaset[i,],
                                                 data=data))
    propw.logprior[i]      <- propset.logprior[i]
    propw.loglike[i]       <- propset.loglike[i]

    currw.logprior[i]      <- currentset.logprior[i]
    currw.loglike[i]       <- currentset.loglike[i]

  }

  mh <- exp( (propset.loglike[n.groups] + propset.logprior[n.groups]) -
                      (currentset.loglike[1] + currentset.logprior[1]))
  if ( !is.na(mh) && (runif(1) < mh) ) {
    use.theta[groups[1],]   <- thetaset[n.groups,]	# swap the 1st with last
    use.logprior[groups[1]] <- propset.logprior[n.groups]
    use.loglike[groups[1]]  <- propset.loglike[n.groups]
  }

  if ( n.groups!=1 ) {										# make sure we are not done yet
    for(i in 1:(n.groups-1)) {
      mh <- exp((propset.loglike[i] + propset.logprior[i]) -
                        (currentset.loglike[i+1] + currentset.logprior[i+1]))
      if( !is.na(mh) && (runif(1) < mh) ) {
        use.theta[groups[i+1],] <- thetaset[i,]
        use.logprior[groups[i+1]] <- propset.logprior[i]
        use.loglike[groups[i+1]] <- propset.loglike[i]
      }
    }
  }

  cbind (use.logprior, use.loglike, use.theta)
}

#' run function
#'
#' a wrapper function calls run_data C++ function to do hierarchical Bayesian
#' sampling. It selects data-level or hyper-level sampling by looking for
#' hyper attribute and the numbers of participants in a sample list.
#'
#' @param samples a sample list generated by calling DMC's samples.dmc.
#' You can also use ggdmc's own initialise to generate a samples.
#' @param report how many iterations to return a report
#' @param cores a switch for computing the prob density for each trial in
#' parallel. Turn it on by setting any number > 1.
#' @param p.migrate set it greater than 0 to use migration samplers. For example
#' p.migrate=0.05 will use migration in 5\% chance.
#' @param gamma.mult a DEMC tuning parameter, affecting the size of jump
#' @param farjump No funciton for compatibility reason
#' @param force No funciton for compatibility reason
#' @param setting a list run setting arguments, such as p.migrate, gamma.mult,
#' report and cores, send to C++ function.
#' @param verbose Turn it on to print loads of debugging informatoin
#' @param debug default as FALSE to use random chain sequence. TRUE uses ordered
#' chain sequence.
#' @export
#' @return a DMC sample
#' @examples
#' m1 <- model.dmc(
#'     p.map     = list(a="1",v="1",z="1",d="1",sz="1",sv="1", t0="1",st0="1"),
#'     constants = c(st0=0, d=0),
#'     match.map = list(M=list(s1="r1", s2="r2")),
#'     factors   = list(S=c("s1", "s2")),
#'     responses = c("r1", "r2"),
#'     type      = "rd")
#'
#' ## Use 6 prior truncated normal distributions
#' p.prior <- prior.p.dmc(
#'   dists = rep("tnorm", 6),
#'   p1    = c(a=2,  v=2.5, z=.5, sz=.3, sv=1,  t0=.3),
#'   p2    = c(a=.5, v=.5,  z=.1, sz=.1, sv=.3, t0=.05),
#'   lower = c(0,-5, 0, 0, 0, 0),
#'   upper = c(5, 7, 2, 2, 2, 2))
#'
#' ## parameter vector. These are the trun values simulating data
#' p.vector <- c(a=1,v=1, z=.5, sz=.25, sv=.2,t0=.15)
#' dat1 <- simulate(m1, nsim=1e2, p.vector=p.vector)
#' mdi1 <- data.model.dmc(dat1, m1)
#'
#' ## Use DMC's plot_cell_density to examine distributions
#' ## Accuracy around 70%
#' par(mfrow=c(1,2))
#' plot_cell_density(data.cell=mdi1[mdi1$S=="s1", ], C="r1", xlim=c(0,2))
#' plot_cell_density(data.cell=mdi1[mdi1$S=="s2", ], C="r2", xlim=c(0,2))
#' par(mfrow=c(1,1))
#'
#' ## ---------------------------
#' ## Profiles all 6 parameters
#' par(mfrow=c(2,3));
#' profile(mdi1, "a",  .1,  2, p.vector)
#' profile(mdi1, "v",  .1,  2, p.vector)
#' profile(mdi1, "z",  .2, .8, p.vector)
#' profile(mdi1, "sz", .1, .9, p.vector)
#' profile(mdi1, "sv", .1,  2, p.vector)
#' profile(mdi1, "t0", .01, .5, p.vector)
#' par(mfrow=c(1,1));
#'
#' ## Initialse a DMC sample
#' ## nthin == 1 (default)
#' ## niter == 100
#' ## prior distributions as listed in p.prior
#' ## data  == model data instance 1
#' ## do not use migrate sampler (default p.migrate=0)
#' samples0 <- samples.dmc(nmc=50, p.prior=p.prior, data=mdi1)
#' samples0 <- run.dmc(samples0)
run.dmc <- function(samples, report=100, cores=1, p.migrate=0, gamma.mult=2.38,
  farjump=NA, force=FALSE, setting=NULL, verbose=FALSE, debug=FALSE)
{
  ## Check 1
  if (!is.na(farjump) && !is.numeric(farjump) && (farjump<1)) {
    stop("farjump must be a number greater than or equal to 1")
  } else {
    farjump <- round (farjump)
  }

  ## Check 2
  if( is.null(setting) )
  {
    setting_in <- list(p.migrate=p.migrate, h.p.migrate=0,
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
  ## Gyakuzuki:
  ## If hyper attribute is not set, it must be data level
  ## If the 1st name is theta, assume we get only one participant
  ## Otherwise, mulitple participants
  ## ---------------------------------------------------------------------- */
  sNames <- names(samples)
  if(sNames[1] != "theta") { stop("Use h.run.dmc for multiple participants")}
  if (debug==TRUE) {
        out <- run_data(samples, setting_in, debug=TRUE) ;
  } else if (cores > 1) {
        out <- run_data_parallel(samples, setting_in) ;
  } else {
        out <- run_data(samples, setting_in) ;
  }
  dimnames(out$theta) <- list(NULL, out$p.names, NULL)

  cat("\n")
  return(dmc(out))
}


trial_log_likes <- function(samples,thin_pointwise=1,
                            chain_dot=TRUE,subject_dot=FALSE)
  # Get pointwise log-likelihoods
{

  n.trials <- dim(samples$data)[1]

  nmc_thin <- seq(thin_pointwise,samples$nmc,by=thin_pointwise)
  trial_log_likes  <-
    array(-Inf,c(length(nmc_thin),samples$n.chains, dim(samples$data)[1]))
  dimnames(trial_log_likes) <- list(nmc_thin,NULL,NULL)
  if (chain_dot) cat("Processing chains: ")
  for (j in 1:samples$n.chains) {
    for (i in nmc_thin)
      trial_log_likes[as.character(i),j,]  <-
          log.likelihood(samples$theta[j,,i],samples$data)
      if (chain_dot) cat(".")
  }
  if (chain_dot) cat("\n")
  if (subject_dot) cat(".")
  trial_log_likes
}



group_trial_log_likes <- function(samples,thin_pointwise=1,max_size=1e8)
  # extracts trial_log_likes from a list of subject fits and concatanates
  {

  cat("Processing subjects: ")
  tll <- trial_log_likes(samples[[1]],thin_pointwise=thin_pointwise,
                         chain_dot=FALSE,subject_dot=TRUE)
  sdim <- dim(tll)
  size <- sum(length(samples)*prod(sdim))
  if (size>max_size) stop(paste("Output size",size,
                                "too large,adjust max_size (",max_size,")"))
  tlls <- lapply(samples[-1],trial_log_likes,thin_pointwise=thin_pointwise,
                 chain_dot=FALSE,subject_dot=TRUE)
  tlls[[length(samples)]] <- tll
  sdims <- cbind(matrix(unlist(lapply(tlls,dim)),nrow=3))
  if ( !all(sdims[1,1]==sdims[1,-1]) || !all(sdims[2,1]==sdims[2,-1]) )
    stop("Subjects must have the same number of interations and chains")
  out <- array(dim=c(dim(tlls[[1]])[-3],sum(sdims[3,])))
  start <- 1; end <- sdims[3,1]
  for (i in 1:length(samples)) {
    out[,,start:end] <- tlls[[i]]
    if (i<length(samples)) {
      start <- end+1
      end <- start - 1 + sdims[3,i+1]
    }
  }
  out
}


group_subject_log_likes <- function(samples)
  # extracts summed_log_likes from a list of subject fits and concatanates
  {

  tlls <- lapply(samples,function(x){x$log_likelihoods})
  sdims <- matrix(unlist(lapply(tlls,dim)),nrow=2)
  if ( !all(sdims[1,1]==sdims[1,-1]) || !all(sdims[2,1]==sdims[2,-1]) )
    stop("Subjects must have the same number of interations and chains")
  out <- array(dim=c(dim(tlls[[1]]),length(samples)))
  for (i in 1:length(samples))
    out[,,i] <- tlls[[i]]
  out
}


run.unstuck.dmc <- function(samples,nmc=NA,rp=.001,report=10,cores=1,
                            cut=10,nbad=0,max.try=10,
                            gamma.mult=2.38,h.gamma.mult=NA)
  # repeats sampling until <= nbad stuck chains as defined by cut
{
  if (is.null(samples$theta))
    stop("For multiple subjects use h.run.unstuck.dmc")
  n.chain <- dim(samples$theta)[1]
  if (any(is.na(samples$theta[,,2])))
    samples <- run.dmc(samples=samples, report=report,cores=cores,
                       gamma.mult=gamma.mult)
  if ( is.na(nmc) ) nmc <- samples$nmc
  try.num <- 1
  repeat {
    if ( length(pick.stuck.dmc(samples)) <= nbad ) break
    samples <- h.run.dmc(h.samples.dmc(nmc=nmc,samples=samples),
                         report=report,cores=cores,
                         gamma.mult=gamma.mult,h.gamma.mult=h.gamma.mult)
    if (try.num >= max.try) break
    try.num <- try.num + 1
  }
  samples <- unstick.dmc(samples,pick.stuck.dmc(samples,cut=cut))
  cat(paste("Removed",n.chain-dim(samples$theta)[1]),"chains\n")
  print(pick.stuck.dmc(samples,verbose=TRUE))
  samples
}


run.converge.dmc <- function(samples,nmc=NA,report=10,cores=1,gamma.mult=2.38,
                             cut=1.1,max.try=10,transform=TRUE,autoburnin=FALSE)
  # repeats sampling until gelman.diag Multivariate psrf < cut
{
  if (is.null(samples$theta))
    stop("For multiple subjects use h.run.converge.dmc")
  if (any(is.na(samples$theta[,,2])))
    samples <- run.dmc(samples=samples, report=report,cores=cores,
                       gamma.mult=gamma.mult)
  if ( is.na(nmc) ) nmc <- samples$nmc
  try.num <- 0
  best.gd <- Inf
  repeat {
    gd <- gelman.diag(theta.as.mcmc.list(samples),
                      autoburnin=autoburnin,transform=transform)$mpsrf
    if (gd < best.gd) {
      best <- samples
      best.gd <- gd
    }
    if ( gd <= cut ) break
    print(paste("Multivariate psrf achieved =",gd))
    samples <- run.dmc(samples.dmc(nmc=nmc,samples=samples),
                       report=report,cores=cores,gamma.mult=gamma.mult)
    if (try.num >= max.try) break
    try.num <- try.num + 1
  }
  print(paste("Final multivariate psrf =",best.gd))
  best
}


post.predict.dmc <- function(samples,n.post=100,probs=c(1:99)/100,
                             bw="nrd0",report=10,save.simulation=FALSE)
  # make list of posterior preditive density, quantiles and response p(robability)
{

  get.dqp <- function(sim,facs,probs) {

    quantile.names <- function(x,probs=seq(0, 1, 0.25),na.rm=FALSE,type=7, ...)
    {
      out <- quantile(x,probs=probs,na.rm=na.rm,type=type,names=FALSE,...)
      names(out) <- probs*100
      out
    }

    qs <- tapply(sim$RT,sim[,c(facs,"R")],quantile.names,probs=probs,na.rm=TRUE)
#     qs <- apply(qs,1:length(dim(qs)),function(x){
#       if ( is.null(x[[1]]) || all(is.na(x[[1]]))) NULL else x[[1]]})
    p <- tapply(sim$RT,sim[,c(facs,"R")],length)
    p[is.na(p)] <- 0 # In case some cells are empty
    p <- p/rep(apply(p,1:length(facs),sum),times=length(levels(sim$R)))
    # cell names
    cell.names <- dimnames(qs)[[1]]
    n.cell <- length(facs)+1
    for ( i in 2:n.cell )
      cell.names <- outer(cell.names,dimnames(qs)[[i]],"paste",sep=".")
    cell.names <- as.vector(cell.names)
    # Get density and make defective
    dens <- tapply(sim$RT,sim[,c(facs,"R")],function(x){
      if (all(is.na(x))) NULL else {
        x <- x[x>=quantile(x,.01) & x<=quantile(x,.99)]
        if (length(x)<2) NULL else density(x,bw=bw)
      }
    })
    for (i in 1:length(p)) if ( !(p[i]==0) )
    {
      if (!is.null(qs[i][[1]])) {
        names(qs[i][[1]]) <- as.numeric(names(qs[i][[1]]))*p[i]/100
        attr(qs[i][[1]],"cell.name") <- cell.names[i]
      }
      if (!is.null(dens[i][[1]]) ) {
        dens[i][[1]]$y <- dens[i][[1]]$y*p[i]
        attr(dens[i][[1]],"cell.name") <- cell.names[i]
      }
    }
    dnd <- dimnames(dens)
    dnq <- dimnames(qs)
    dens <- apply(dens,1:length(facs),function(x){x})
    qs <- apply(qs,1:length(facs),function(x){x})
    if ( is.null(dim(dens)) ) {
      dens <- array(dens,dim=c(length(dens)))
      dimnames(dens) <- dnd[-length(dnd)]
      qs <- array(qs,dim=c(length(qs)))
      dimnames(qs) <- dnq[-length(dnq)]
    }
    list(pdf=dens,cdf=qs,p=p)
  }

  model <- attributes(samples$data)$model
  facs <- names(attr(model,"factors"))
  resp <- names(attr(model,"responses"))
  ns <- table(samples$data[,facs])
  n.par <- dim(samples$theta)[2]
  thetas <- matrix(aperm(samples$theta,c(1,3,2)),ncol=n.par)
  colnames(thetas) <- dimnames(samples$theta)[[2]]
  posts <- thetas[sample(c(1:dim(thetas)[1]),n.post,replace=F),]
  n.rep <- sum(ns)
  sim <- data.frame(matrix(nrow=n.post*n.rep,ncol=dim(samples$data)[2]))
  names(sim) <- names(samples$data)

  # Tweaks for Stop Signal
  if (!any(names(samples$data)=="SSD")) {
    SSD <- rep(Inf,sum(ns))
    leave.out <- -c(1:dim(samples$data)[2])[names(samples$data)=="RT"]
  } else {
# Assumes last two are SSD and RT! FIX ME.
    SSD <- unlist(tapply(samples$data$SSD,samples$data[,facs],identity))
    leave.out <- -c((dim(samples$data)[2]-1):dim(samples$data)[2])
  }
  cat("\n")
  cat(paste("Simulating (\'.\'=",report,"): ",sep=""))
  for (i in names(samples$data)[leave.out])
    sim[[i]] <- factor(rep(NA,n.post*n.rep),levels=levels(samples$data[[i]]))
  for (i in 1:n.post) {
    tmp <- simulate.dmc(model, nsim=ns, p.vector=posts[i,], SSD=SSD)
    sim[(1+(i-1)*n.rep):(i*n.rep),names(tmp)] <- tmp
    if ( (i %% report) == 0) cat(".")
  }

  if (save.simulation) {
    sim <- cbind(reps=rep(1:n.post,each=dim(samples$data)[1]),sim)
    attr(sim,"data") <- samples$data
    sim
  } else {
    sim.dqp <- get.dqp(sim,facs,probs)
    dat.dqp <- get.dqp(sim=samples$data,facs,probs)
    names(dat.dqp) <- paste("data",names(dat.dqp),sep=".")
    out <- c(sim.dqp,dat.dqp)
    return(out)
  }
}


ppp.dmc <- function(samples,fun=function(x){mean(x$RT)},n.post=500,
                    plot.density=TRUE,main="",bw="nrd0",report=10)
  # posterior predictive (pp) p value for function fun of data (p(observed)>pp)
{

  model <- attributes(samples$data)$model
  facs <- names(attr(model,"factors"))
  resp <- names(attr(model,"responses"))
  ns <- table(samples$data[,facs])
  n.par <- dim(samples$theta)[2]
  thetas <- matrix(aperm(samples$theta,c(1,3,2)),ncol=n.par)
  colnames(thetas) <- dimnames(samples$theta)[[2]]
  posts <- thetas[sample(c(1:dim(thetas)[1]),n.post,replace=F),]
  sim <- vector(mode="list",length=n.post)
  cat("\n")
  cat(paste("Simulating (\'.\'=",report,"): ",sep=""))
  for (i in 1:n.post) {
    sim[[i]] <- simulate.dmc(model, nsim=ns, p.vector=posts[i,])
    if ( (i %% report) == 0) cat(".")
  }
  cat("\n")
  pp <- unlist(lapply(sim,fun))
  obs <- fun(samples$data)
  ppp <- mean(obs>pp)
  if (plot.density) {
    plot(density(pp,bw=bw),main=main)
    abline(v=obs,lty=2)
  }
  ppp
}


#' Truncated Normal Distribution
#'
#' The two wrapper functions call Rcpp functions to speed up computation. The
#' codes are based on Christopher Jackson's \pkg{msm} (1.5) package and Jonathan
#' Olmsted's \pkg{RcppTN} (0.1-8) package.
#'
#' \code{dtnorm} calculates probability density for a truncated normal
#' distribution with mean equal to mean and standard deviation equal to
#' sd before truncation, and truncated on the interval [lower, upper].
#'
#' \code{rtnorm} generates a random number from a truncated Normal
#' distribution with mean equal to mean and standard deviation equal
#' to sd before truncation, and truncated on the interval [lower, upper].
#'
#' @param x,n x in \code{dtnorm} is a vector of quantiles; n in
#' \code{rtnorm} indicates how many random number to generate
#' @param mean vector of means
#' @param sd vector of standard deviations
#' @param lower lower truncation point. Default as -Inf
#' @param upper upper truncation point. Default as Inf.
#' @param log default 0 (FALSE). Enter 1 to make it TRUE. Whether to
#' calculate density. Use only in dtnorm.
#' @param checks a switch to turn on check functions to check inputs and ouput.
#' Default to FALSE for the sake of speed.
#' @export
#' @examples
#' ## Use dtnorm and rtnorm with their default values
#' dtnorm()
#' rtnorm()
#'
#' ## A similar curve plotting example extracted from dnorm functions
#' plot(function(x) dnorm(x, log = FALSE), -2.5, 2.5,
#'      main = "Normal Distribution", ylim=c(0,0.45), ylab="Density")
#' curve(dtnorm(x, lower=-2, upper=2), add=TRUE, col="tomato", lwd=2)
#' mtext("dnorm(x)", adj = 0)
#' mtext("dtnorm(x)", col = "tomato", adj = 1)
dtnorm <- function(x=0, mean=0, sd=1, lower=-Inf, upper=Inf, log=0, checks=FALSE)
{
  mean  <- rep(mean,  length=length(x))
  sd    <- rep(sd,    length=length(x))
  lower <- rep(lower, length=length(x))
  upper <- rep(upper, length=length(x))
  log   <- rep(log, length=length(x))

  if (checks) { checkInputs(mean, sd, lower, upper) }

  ## inside C++, I use islog as the variable name for the boolean log, so
  ## I can use cmath's log function
  out <- .Call("dtn_wrapper", x_ = x, mean_ = mean, sd_ = sd, lower_ = lower,
    upper_ = upper, islog_=log)

  if (checks) { checkOutputs(out) }
  return(out)
}

#' @export
#' @rdname dtnorm
rtnorm <- function(n=1, mean=0, sd=1, lower=-Inf, upper=Inf, checks=TRUE)
{
  ## The size of mean will determine the draw numbers
  if (length(n) > 1) n <- length(n)
  mean  <- rep(mean,  length=n)
  sd    <- rep(sd,    length=n)
  lower <- rep(lower, length=n)
  upper <- rep(upper, length=n)

  if (checks) { checkInputs(mean, sd, lower, upper) }

  out <- .Call("rtn_wrapper", mean_ = mean,   sd_ = sd,
    lower_ = lower, upper_ = upper)

  if (checks) { checkOutputs(out) }
  return(out)
}

checkInputs <- function(m, s, l, u)
{
  if (!(length(m) == length(s) & length(m) == length(l) & length(m) == length(u)))
  {
    stop("Input vectors are not all same length. Nothing done.")
  }
}

checkOutputs <- function(out)
{
  if (any(is.na(out)))
  {
    warning("NAs returned. Check for invalid parameters.")
  }
}



