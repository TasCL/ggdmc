#' Censor missing values and RT outliers
#'
#' \code{censor} requests a data frame with at lease three columns, stimulus
#' type (S factor), response/accumulator type (R factor) and response time
#' (RT, double) as input \code{data}. It adding a boolean
#' column by scoring R column against S column. This produces a C column.
#' Otherwise plots each response. When NAs occur, \code{censor} reports a
#' summary of p(NA).
#'
#' @param x a data frame for a design cell.
#' @param xlim the lower and upper boundaries for censoring
#' @keywords censor
#' @export
#' @examples
#' model <- model.dmc(
#'    p.map     = list(a="1",v="1",z="1",d="1",sz="1",sv="1", t0="1",st0="1"),
#'    constants = c(st0=0,d=0),
#'    match.map = list(M=list(s1="r1",s2="r2")),
#'    factors   = list(S=c("s1","s2")),
#'    responses = c("r1","r2"),
#'    type      = "rd")
#'
#' pVec <- c(a=1,v=1, z=0.5, sz=0.25, sv=0.2,t0=.15)
#' dat1 <- simulate(model, nsim=1e2, p.vector=pVec)
#' mdi1 <- data.model.dmc(dat1, model)
#'
#' ## 1. xlim=c(.2, 3) defines the censored lower and upper bounds at .2 s and
#' ## 3 s
#' ## 2. mdi1[mdi1$S=="s1", ] extracts the data frame with rows/trials equal to
#' ## "s1".
#' mdi1.s1 <- mdi1[mdi1$S=="s1", ]
#' censored.s1 <- censor(mdi1.s1, xlim=c(.2, 3))
#' head(censored.s1)
#' ## Below is printed by dplyr::tbl_df(censored.s1)
#' ## Source: local data frame [94 x 3]
#' ##         S      R        RT
#' ##    (fctr) (fctr)     (dbl)
#' ## 1      s1     r2 0.3370070
#' ## 2      s1     r1 0.2956996
#' ## 3      s1     r2 0.2128732
#' ## 4      s1     r1 0.3710657
#' ## 5      s1     r2 0.3083122
#' ## 6      s1     r2 0.2997714
#' ## 7      s1     r1 0.4123926
#' ## 8      s1     r1 0.2356586
#' ## 9      s1     r2 0.4721079
#' ## 10     s1     r2 0.3524906
#' ## ..    ...    ...       ...
censor <- function(x, xlim=c(0, Inf)) {
  if ( !all(names(x) %in% c("S","R","RT")) )
    stop("A minimial three columns, S, R and RT are required.")

  ## Make R column as.factor
  if (!is.factor(x$R)) x$R <- factor(x$R)
  p.na <- mean(is.na(x$RT))   ## proportion of na response
  is.in <- !is.na(x$RT) ## TRUE/FALSE index which row is not NA
  is.in[is.in] <- x$RT[is.in]>xlim[1] & x$RT[is.in]<xlim[2]
  return(x[is.in,])
}


#' Calculate Probability Density for an Experimental Condition
#'
#' A S3 method derived from \pkg{stats} density to get probability density
#' for a specified experimental condition. The user has to indicate
#' a correctness column.
#'
#' @param x a data.cell as a data frame as welle dmc class, containing
#' only one experimental condition.
#' @param C a string, e.g., "r1" or a logical vector to represent the
#' correctness column in a data frame
#' @param ... other arguments
#' @keywords density.dmc
#' @export
#' @importFrom stats density
#' @examples
#' model <- model.dmc(
#'  p.map     = list(a="1",v="1",z="1",d="1",sz="1",sv="1", t0="1",st0="1"),
#'  constants = c(st0=0,d=0),
#'  match.map = list(M=list(s1="r1",s2="r2")),
#'  factors   = list(S=c("s1","s2")),
#'  responses = c("r1","r2"),
#'  type      = "rd")
#'
#' pVec <- c(a=1,v=1, z=0.5, sz=0.25, sv=0.2,t0=.15)
#' dat1 <- simulate(model, nsim=1e2, p.vector=pVec)
#' mdi1 <- data.model.dmc(dat1, model)
#' censored.s1  <- censor(mdi1[mdi1$S=="s1", ], xlim=c(.2, 3))
#' s1r1.density <- density(censored.s1, C="r1")
#' head(s1r1.density)
#' ## Source: local data frame [1,024 x 3]
#' ##             x           y     C
#' ##         (dbl)       (dbl) (lgl)
#' ## 1  0.02077313 0.004432689  TRUE
#' ## 2  0.02276278 0.004923698  TRUE
#' ## 3  0.02475243 0.005480113  TRUE
#' ## 4  0.02674208 0.006079448  TRUE
#' ## 5  0.02873173 0.006731685  TRUE
#' ## 6  0.03072138 0.007471079  TRUE
#' ## 7  0.03271103 0.008265728  TRUE
#' ## 8  0.03470068 0.009124208  TRUE
#' ## 9  0.03669033 0.010097473  TRUE
#' ## 10 0.03867998 0.011141091  TRUE
#' ## ..        ...         ...   ...
density.dmc <- function(x, C=NA, ...) {
  if (length(C)==1) C <- rep(C, nrow(x));
  is.in <- !is.na(x$RT);
  p.na <- mean(is.na(x$RT));

  if ( !any(is.na(C)) ) {

    if (is.logical(C) & length(C)==nrow(x)) {
      x$C <- C[is.in]
    } else {
      x$C <- x$R==C[is.in];
    }

    if (length(x$RT[x$C]) > 2) {
      dO <- density(x$RT[x$C]) ## correct density
    } else {
      dO <- NULL
    }

    if (length(x$RT[!x$C]) > 2) {
      dX <- density(x$RT[!x$C])
    } else {
      dX <- NULL
    }

    if (is.null(dX) & is.null(dO)) stop("No density to plot")
    acc <- mean(x$C)         ## Calculate the accuracy rate

    ## Adjust pdf's with acc and na rates; error and na rates
    if (!is.null(dO)) dO$y <- dO$y*acc*(1-p.na)
    if (!is.null(dX)) dX$y <- dX$y*(1-acc)*(1-p.na)
  }

  ## Transform pdf data as a data frame
  df <- rbind(data.frame(x=dO$x, y=dO$y, C=TRUE),
    data.frame(x=dX$x, y=dX$y, C=FALSE))
  attr(df, "accuracy") <- acc
  return(df)
}

#' Convert factor levels to a data frame
#'
#' \code{fac2df} takes a model object created by model.dmc and returns a
#' data frame with all combination of factor levels.
#'
#' @param model a model object created by \code{model.dmc}
#' @return a data frame
#' @keywords fac2df
#' @export
#' @examples
#' m1 <- model.dmc(
#'       p.map     = list(a="1",v="F",z="1",d="1",sz="1",sv="1",t0="1",st0="1"),
#'       match.map = list(M=list(s1="r1",s2="r2")),
#'       factors   = list(S=c("s1","s2"), F=c("f1","f2")),
#'       constants = c(st0=0,d=0),
#'       responses = c("r1","r2"),
#'       type      = "rd")
#' df1 <- fac2df(m1)
#' str(df1)
#' ## 'data.frame':	4 obs. of  2 variables:
#' ## $ S: Factor w/ 2 levels "s1","s2": 1 2 1 2
#' ## $ F: Factor w/ 2 levels "f1","f2": 1 1 2 2
#'
#' m2 <- model.dmc(
#'             p.map=list(A="1",B="1",v="M",sv="M",t0="1",st0="1"),
#'             constants=c(st0=0,sv.false=1),
#'             match.map=list(M=list(s1=1,s2=2)),
#'             factors=list(S=c("s1","s2")),
#'             responses=c("r1","r2"),
#'             type="lnorm")
#' fac2df(m2)
#'
#' m3 <- model.dmc(
#'   p.map     = list(A="1",B="R",t0="1",mean_v=c("F","M"),sd_v="M",st0="1"),
#'   match.map = list(M=list(s1=1,s2=2)),
#'   factors   = list(S=c("s1","s2"),F=c("f1","f2")),
#'   constants = c(sd_v.false=1,st0=0),
#'   responses = c("r1","r2"),
#'   type="norm")
#'
#' fac2df(m3)
#' ##    S  F
#' ## 1 s1 f1
#' ## 2 s2 f1
#' ## 3 s1 f2
#' ## 4 s2 f2
fac2df  <- function(model) {
  facs  <- lapply(strsplit(dimnames(model)[[1]],"\\."), function(x){x[-length(x)]})
  facs  <- facs[ 1: ( length(facs) / length(attr(model,"responses")) ) ]
  fnams <- names(attr(model, "factors"))
  facs  <- data.frame(t(matrix(unlist(facs), nrow=length(fnams))))
  names(facs) <- fnams
  return(facs)
}

makePrior.dmc <- function(p1, p2, p.names=rep(NA,length(p1)),
  lower=rep(NA,length(p1)), upper=rep(NA,length(p1)),
  dists=rep("tnorm",length(p1)),
  untrans=rep("identity",length(p1)),
  dist.types=c("tnorm","beta","gamma","lnorm","constant"),log=FALSE) {
  ## From prior.p.dmc.R
  ## Makes a list of prior distribution parameters.
  lenp1 <- length(p1);
  lenp2 <- length(p2)
  len.lower <- length(lower)
  len.upper <- length(upper)

  if (lenp2==1) p2 <- rep(p2,lenp1)
  if (lenp1 != lenp2)     stop("p1 and p2 must have the same length")
  if (lenp1 != len.lower) stop("p1 and lower must have the same length")
  if (lenp1 != len.upper) stop("p1 and upper must have the same length")

  both.not.na <- !is.na(upper) & !is.na(lower)

  if (any(upper[both.not.na]<=lower[both.not.na]) )
    stop("All elements of upper must be greater than lower")
  if (lenp1 != length(dists)) stop("p1 and dists must have the same length")
  if (!all(dists %in% dist.types)) stop(paste("Unsupported distribution, allowable types are:",
                                              paste(dist.types,collapse=", ")))
  name.untrans <- length(untrans) != lenp1

  if (name.untrans & (is.null(names(untrans)) | is.null(names(p1))))
    stop("If untrans vector is not the same length as p1 it must have p1 names")
  if (!(all(names(untrans) %in% names(p1)))) stop("untrans vector has names not in p1 names")

  prior <- vector(mode="list", length=lenp1)
  names(p1) <- p.names
  names(prior) <- p.names
  for (i in 1:lenp1) {
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
    prior[[i]]$log <- log
    if (!name.untrans) attr(prior[[i]],"untrans") <- untrans[i] else
      if (is.na(untrans[names(p1)[i]]))
        attr(prior[[i]],"untrans") <- "identity" else
          attr(prior[[i]],"untrans") <- untrans[names(p1)[i]]
  }
  return(prior)
}

isMatOrVec <- function(n) {
  out <- (!is.matrix(n) & (is.vector(n) & length(n) != 1))
  return(out)
}

#' Create a mcmc.list in DMC format
#'
#' If not hyper, by default changes samples$theta to a mcmc.list object and
#' plots If pll.chain=TRUE changes samples$pll to an mcmc object and plots
#' posterior log-likelihood of chains. If pll.barplot=TRUE plots their means
#' as a barplot. only.prior and only.like pick out one or other component of
#' pll. If hyper same but at the the hyper  (phi) level.
#'
#' @param samples a data frame stored choice-RT data
#' @param hyper a DMC model
#' @param start starting iteration
#' @param end end iteration
#' @param save.ll whether save log likelihood
#' @param ask Default TRUE
#' @param main.pll Default NULL
#' @param pll.chain Default FALSE
#' @param pll.together Default TRUE
#' @param pll.barplot Default FALSE
#' @param only.prior plot only prior distribution
#' @param only.like plot only likelihood
#' @param subject plot which participant. Default the 1st
#' @param layout graphic layout. Default NA
#' @param ... other arguments
#' @export
mcmc.list.dmc <- function(samples,hyper=FALSE,start=1,end=NA, save.ll=FALSE,
  ask=TRUE, main.pll=NULL, pll.chain=FALSE, pll.together=TRUE,
  pll.barplot=FALSE, only.prior=FALSE, only.like=FALSE, subject=1,
  layout=NA, ...) {
  auto.layout <- any(is.na(layout))
  ## hyper-parameter
  if ( hyper ) {
    hyper <- attr(samples,"hyper")
    if (is.null(hyper))
      stop("There are no hyper-parameters to plot.")
    if ( is.na(end) ) end <- hyper$nmc
    if ( end <= start )
      stop("End must be greater than start")
    if ( pll.chain | pll.barplot ) {
      chain.pll <- hyper$h_summed_log_prior[start:end,] +
        hyper$h_log_likelihoods[start:end,]
      colnames(chain.pll) <- 1:dim(chain.pll)[2]
      if (pll.barplot) {
        mean.ll <- apply(chain.pll,2,mean)
        names(mean.ll) <- 1:length(mean.ll)
        barplot(mean.ll,ylab="Mean Post Like",main=main.pll)
        if (save.ll) mean.ll
      } else {
        if (!auto.layout) par(mfrow=layout)
        if (!pll.together)
          plot(mcmc(chain.pll),auto.layout=auto.layout,ask=ask,...) else
            plot(mcmc.list(lapply(data.frame(chain.pll),function(x){
              mcmc(x)})),auto.layout=auto.layout,ask=ask,...)
      }
    } else {
      if (!auto.layout) par(mfrow=layout)
      plot(window(phi.as.mcmc.list(hyper,start=start,end=end),thin=1),
           auto.layout=auto.layout,ask=ask,...)
    }
    ## non-hyper
  } else {
    if (is.null(samples$theta)) samples <- samples[[subject]]
    if ( is.na(end) ) end <- samples$nmc
    if ( end <= start ) stop("End must be greater than start")

    if ( pll.chain | pll.barplot ) {
      if (only.prior) {
        chain.pll <- samples$summed_log_prior[start:end,]
      } else if (only.like) {
        chain.pll <- samples$log_likelihoods[start:end,]
      } else {
        chain.pll <- samples$summed_log_prior[start:end,] + samples$log_likelihoods[start:end,]
      }
      colnames(chain.pll) <- 1:dim(chain.pll)[2]
      d <- chain.pll

      if (pll.barplot) {
        mean.ll <- apply(chain.pll,2,mean)
        names(mean.ll) <- 1:length(mean.ll)
        barplot(mean.ll,ylab="Mean Post Like",main=main.pll)
        if (save.ll) mean.ll
      } else {
        if (!auto.layout) par(mfrow=layout)
        if (!pll.together) {
          ## plot(mcmc(chain.pll),auto.layout=auto.layout,ask=ask,...)
          ## d <- mcmc(chain.pll)
        } else {
          #plot(mcmc.list(lapply(data.frame(chain.pll),function(x){
          #  mcmc(x)})),auto.layout=auto.layout,ask=ask,...)
          lst <- list()
          for(i in 1:samples$n.chains) { lst[[i]] <- mcmc(matrix(chain.pll[,i])) }
          d <- mcmc.list(lst)
        }

      }
    } else {
      if (!auto.layout) par(mfrow=layout)
      d <- window(theta.as.mcmc.list(samples),start=start,end=end,thin=1)
    }
  }
  return(d)
}

## A prototype function preparing for converting to C++
post_predict_dmc <- function(samples,n.post=100,probs=c(1:99)/100,
  bw="nrd0", report=10, save.simulation=FALSE) {
  # make list of posterior preditive density, quantiles and response p(robability)
  get.dqp <- function(sim, facs, probs) {

    quantile.names <- function(x,probs=seq(0, 1, 0.25),na.rm=FALSE,type=7, ...)
    {
      out <- quantile(x, probs=probs, na.rm=na.rm, type=type, names=FALSE,...)
      names(out) <- probs*100
      return(out)
    }

    qs <- tapply(sim$RT,sim[,c(facs,"R")],quantile.names,probs=probs,na.rm=TRUE)
    p  <- tapply(sim$RT,sim[,c(facs,"R")],length)
    p[is.na(p)] <- 0 # In case some cells are empty
    p <- p/rep(apply(p,1:length(facs),sum),times=length(levels(sim$R)))

    ## cell names
    cell.names <- dimnames(qs)[[1]]
    n.cell     <- length(facs)+1
    for ( i in 2:n.cell ) {
      cell.names <- outer(cell.names,dimnames(qs)[[i]],"paste",sep=".")
    }

    cell.names <- as.vector(cell.names)
    ## Get density and make defective
    dens <- tapply(sim$RT,sim[,c(facs,"R")],function(x){
      if (all(is.na(x))) NULL else {
        x <- x[x>=quantile(x,.01) & x<=quantile(x,.99)]
        if (length(x)<2) NULL else density(x,bw=bw)
      }
    })
    for (i in 1:length(p)) if ( !(p[i]==0) )
    {
      names(qs[i][[1]]) <- as.numeric(names(qs[i][[1]]))*p[i]/100
      #        as.numeric(unlist(strsplit(names(qs[i][[1]]),"%")))*p[i]/100
      attr(qs[i][[1]],"cell.name") <- cell.names[i]
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
  facs  <- names(attr(model,"factors")); nfacs <- length(facs)
  resp  <- names(attr(model,"responses")) ## sometimes do not have names
  ns <- table(samples$data[,facs])       ## number of level of stimulus factor
  n.rep <- sum(ns)                       ## total number of data points
  n.par <- dim(samples$theta)[2]  ## dim[1]=chains; dim[2]=list of para; dim[3]=samples

  ## accumulate chain1, chain2, chain3, ... all on one column for each para
  thetas <- matrix(aperm(samples$theta,c(1,3,2)),ncol=n.par)
  colnames(thetas) <- dimnames(samples$theta)[[2]]

  ## random select n.post samples from thetas for each para
  posts <- thetas[sample(c(1:nrow(thetas)), n.post, replace=F),]

  ## Construct a NA container with same dimension as data to store simulation data
  ## Should change to data.table or tcl_df
  sim <- data.frame(matrix(nrow=n.post*n.rep, ncol=ncol(samples$data)))
  names(sim) <- names(samples$data)

  # Tweaks for Stop Signal; unclear comments only work for himself
  if (!any(names(samples$data)=="SSD")) {  ## When no SSD column
    ## Create an Inf container (vector) with same number of samples
    SSD <- rep(Inf,sum(ns))
    leave.out <- -ncol(samples$data)
  } else {  ## Usualy not this route
    ## This works only when data contains a SSD column
    SSD <- unlist( tapply(samples$data$SSD, samples$data[,facs], identity) )
    ## From the last column to one column before it
    leave.out <- -c( (ncol(samples$data)-1):ncol(samples$data) )
  }

  for (i in names(samples$data)[leave.out]) {
    sim[[i]] <- factor(rep(NA,n.post*n.rep),levels=levels(samples$data[[i]]))
  }
  ## Every "report" step prints one dot
  message(paste0("\nSimulating ('.'= ", report, " steps): "))

  ## Key step
  ## 1:1*nrep, 20001:2*n.rep, 40001:3*n.rep
  for (i in 1:n.post) {
    sim[(1+(i-1)*n.rep):(i*n.rep),] <- simulate.dmc(model,
      p.vector=posts[i,], nsim=ns, SSD=SSD)
    if ( (i %% report) == 0) message(".", appendLF = FALSE)
  }

  sim <- cbind(reps=rep(1:n.post,each=dim(samples$data)[1]),sim)
  attr(sim, "data") <- samples$data
  return(sim)
  cat("\n")
}

#' Simulate Post-predictive Sample
#'
#' Make a list of probability density, quantiles and resposnes based on a
#' posterior DMC sample. This function is to produce a data frame that can be
#' used by ggplot2 functions.
#'
#' @param samples a DMC sample/object
#' @param n.post how many data point to simulate
#' @param probs the quantile probability. Default c(1:99)/100
#' @param bw bin width method. Default "nrd0", passing to density function
#' @param report report progress
#' @param save.dat whether to save simulation data to further modify data
#' frame for plotting
#' @importFrom ggmcmc ggs ggs_autocorrelation
#' @export
post.predict.ggdmc <- function(samples, n.post=100, probs=c(1:99)/100,
                             bw="nrd0", report=10, save.dat=FALSE)
{
  # make list of posterior preditive density, quantiles and response p(robability)
  get.dqp <- function(sim, facs, probs) {

    quantile.names <- function(x,probs=seq(0, 1, 0.25),na.rm=FALSE,type=7, ...)
    {
      out <- quantile(x, probs=probs, na.rm=na.rm, type=type, names=FALSE,...)
      names(out) <- probs*100
      return(out)
    }

    qs <- tapply(sim$RT,sim[,c(facs,"R")],quantile.names,probs=probs,na.rm=TRUE)
    p <- tapply(sim$RT,sim[,c(facs,"R")],length)
    p[is.na(p)] <- 0 # In case some cells are empty
    p <- p/rep(apply(p,1:length(facs),sum),times=length(levels(sim$R)))

    ## cell names
    cell.names <- dimnames(qs)[[1]]
    n.cell <- length(facs)+1
    for ( i in 2:n.cell ) {
      cell.names <- outer(cell.names,dimnames(qs)[[i]],"paste",sep=".")
    }

    cell.names <- as.vector(cell.names)
    ## Get density and make defective
    dens <- tapply(sim$RT,sim[,c(facs,"R")],function(x){
      if (all(is.na(x))) NULL else {
        x <- x[x>=quantile(x,.01) & x<=quantile(x,.99)]
        if (length(x)<2) NULL else density(x,bw=bw)
      }
    })
    for (i in 1:length(p)) if ( !(p[i]==0) )
    {
      names(qs[i][[1]]) <- as.numeric(names(qs[i][[1]]))*p[i]/100
      #        as.numeric(unlist(strsplit(names(qs[i][[1]]),"%")))*p[i]/100
      attr(qs[i][[1]],"cell.name") <- cell.names[i]
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
  facs  <- names(attr(model,"factors")); nfacs <- length(facs)
  resp  <- names(attr(model,"responses")) ## sometimes do not have names
  ns <- table(samples$data[,facs])       ## number of level of stimulus factor
  n.rep <- sum(ns)                       ## total number of data points
  n.par <- dim(samples$theta)[2]         ## dim[1]=chains; dim[2]=list of para; dim[3]=samples

  ## accumulate chain1, chain2, chain3, ... all on one column for each para
  thetas <- matrix(aperm(samples$theta,c(1,3,2)),ncol=n.par)
  colnames(thetas) <- dimnames(samples$theta)[[2]]

  ## random select n.post samples from thetas for each para
  posts <- thetas[sample(c(1:nrow(thetas)), n.post, replace=F),]

  ## Construct a NA container with same dimension as data to store simulation data
  ## Should change to data.table or tcl_df
  sim <- data.frame(matrix(nrow=n.post*n.rep, ncol=ncol(samples$data)))
  names(sim) <- names(samples$data)

  # Tweaks for Stop Signal; unclear comments only work for himself
  if (!any(names(samples$data)=="SSD")) {  ## When no SSD column
    ## Create an Inf container (vector) with same number of samples
    SSD <- rep(Inf,sum(ns))
    leave.out <- -ncol(samples$data)
  } else {  ## Usualy not this route
    ## This works only when data contains a SSD column
    SSD <- unlist( tapply(samples$data$SSD, samples$data[,facs], identity) )
    ## From the last column to one column before it
    leave.out <- -c( (ncol(samples$data)-1):ncol(samples$data) )
  }

  for (i in names(samples$data)[leave.out]) {
    sim[[i]] <- factor(rep(NA,n.post*n.rep),levels=levels(samples$data[[i]]))
  }
  ## Every "report" step prints one dot
  cat(paste0("Simulating ('.'= ", report, " steps): "))

  # simulate.dmc <- function(object, nsim=2, p.vector=NULL, SSD=Inf, staircase=NA,
  #   TRIALS=NA)

  ## Key step; 1:1*nrep, 20001:2*n.rep, 40001:3*n.rep
  for (i in 1:n.post)
  {
    sim[(1+(i-1)*n.rep):(i*n.rep),] <- simulate.dmc(model, nsim=ns,
      p.vector=posts[i,], SSD=SSD)
    if ( (i %% report) == 0) cat(".")
  }


  if (save.dat) {
    sim <- cbind(reps=rep(1:n.post,each=dim(samples$data)[1]),sim)
    attr(sim,"data") <- samples$data
    return(sim)
  } else {
    sim.dqp <- get.dqp(sim,facs,probs)
    dat.dqp <- get.dqp(sim=samples$data,facs,probs)  ## replace sim with samples$data
    x1 <- x3 <- NULL
    pdfList <- list(dat.dqp[[1]],sim.dqp[[1]])
    cdfList <- list(dat.dqp[[2]],sim.dqp[[2]])

    ### PDF-data
    for(k in 1:length(pdfList)) {
      dat.pdf <- lapply(pdfList[[k]], function(x) {  ## each factor level
        lapply(x, function(s){                       ## each response type
          x <- s$x; y <- s$y
          if(length(x) == length(y)) {
            npost <- length(x)
          } else {
            stop("x, y unequal legnth")
          }
          cell.name <- strsplit( attr(s,"cell.name"), "\\.")[[1]]
          out <- matrix(rep(NA,npost*(nfacs+3)), ncol=(nfacs+3))
          out[,1:2] <- cbind(x, y)
          out[,3:(ncol(out)-1)] <- matrix( rep(cell.name[1:(length(cell.name)-1)],
                                               each=npost), ncol=nfacs )
          out[,ncol(out)] <- rep(cell.name[length(cell.name)], npost)
          return( data.frame(out) )
          })
      })

      x0 <- NULL
      for(i in 1:length(dat.pdf)) {
        for(j in 1:length(dat.pdf[[i]])) {
          x0 <- rbind(x0,dat.pdf[[i]][[j]])
        }
      }
      x0$mode <- ifelse(k == 1,"data", "simulation")
      x1 <- rbind(x1,x0)
    }
    x1$pcdf <- "pdf"

    ## CDF-data
    for(k in 1:length(cdfList)) {
      dat.cdf <- lapply(cdfList[[k]], function(x){
        lapply(x, function(s){
          y <- as.numeric( names( s ))
          x <- as.vector( s )
          if(length(x) == length(y)) { npost <- length(x) } else { stop("x, y unequal
            legnth") }
          cell.name <- strsplit( attr(s,"cell.name"), "\\.")[[1]]
          out <- matrix(rep(NA,npost*(nfacs+3)), ncol=(nfacs+3))
          out[,1:2] <- cbind(x, y)
          out[,3:(ncol(out)-1)] <- matrix( rep(cell.name[1:(length(cell.name)-1)],
                                               each=npost), ncol=nfacs )
          out[,ncol(out)] <- rep(cell.name[length(cell.name)], npost)
          return( data.frame(out) )
        })
      })

      x2 <- NULL
      for(i in 1:length(dat.cdf)) {
        for(j in 1:length(dat.cdf[[i]])) {
          x2 <- rbind(x2,dat.cdf[[i]][[j]])
        }
      }
      x2$mode <- ifelse(k == 1,"data","simulation")
      x3 <- rbind(x3,x2)
    }
    x3$pcdf <- "cdf"

    ## Final tidy up
    out <- data.frame(rbind(x1,x3))
    names(out) <- c("x","y",facs,"R","mode","pcdf")
    out$x <- as.numeric(as.character(out$x))
    out$y <- as.numeric(as.character(out$y))

    ## Add multiple classes
    attr(out, "class") <- c("pp.ggdmc", "data.frame")
    return(out)
  }
  cat("\n")
}

h.post.predict.dmc <- function(samples,n.post=100,probs=c(1:99)/100,bw="nrd0")
  # apply lost.predict to each subject
{
  lapply(samples, ggdmc::post.predict.ggdmc, n.post=n.post, probs=probs, bw=bw)
}


ppp.dmc <- function(samples,fun=function(x){mean(x$RT)}, n.post=500,
  bw="nrd0",report=10)
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
  print(ppp)
  return( list(obs,pp,ppp) )
}

# hyper2db.dmc <- function(d) {
#   ## Convert Andrew's hyperparameter data frame deep inside a mdi's attribute
#   ## to a compact and efficient tbl_df
#   df <- tbl_df( data.frame(attr(d,"parameters") ))
#   tmp1 <- apply(df, 2, function(s) {
#     data.frame(x=density(s)$x,y=density(s)$y)
#   })
#
#   x0 <- NULL
#   for(i in 1:length(tmp1)) {
#     tmp1[[i]]$para <- names(tmp1)[i]
#     x0 <- rbind(x0,tmp1[[i]])
#   }
#   return(list( df, tbl_df(x0)))
# }

#' view function
#'
#' a convenient function to rearrange p.prior or an element in pp.prior as a
#' data frame to easily inspect. That is another function called View
#' (uppercase) in utils packages.
#'
#' @param p.prior a prior list
#' @return a p.prior data frame
#' @export
#' @examples
#' # Use a single value to replace all missing values
#' pop.mean  <- c(a=1,  v.f1=1,  v.f2=.2, z=0.5, sz=0.3,  sv.f1=0.25,
#'     sv.f2=0.23, t0=0.3)
#' pop.scale <-c(a=0.2,v.f1=.2, v.f2=.2, z=0.1, sz=0.05, sv.f1=.05,
#'     sv.f2=.05,  t0=0.05)
#' p.prior <- prior.p.dmc(
#'   dists = rep("tnorm", 8),
#'   p1=pop.mean,
#'   p2=pop.scale,
#'   lower=c(0,-5, -5, 0, 0, 0, 0,0),
#'   upper=c(2, 5, 5, 1, 2, 2, 1, 1))
#' view(p.prior)
view <- function(p.prior)
{
  colLen <- length(p.prior[[1]]) + 2; colLen
  npars  <- length(p.prior); npars
  tmp <- matrix(numeric(  npars*colLen), nrow=npars); tmp
  for(i in 1:npars)
  {
    add1 <- attr(p.prior[[i]],"dist"); add1
    add2 <- attr(p.prior[[i]],"untrans"); add2
    rowObj <- c(unlist(p.prior[[i]]), add1, add2); rowObj
    if(add1=="gamma_l")
    {
      rowObj <- c(rowObj[1:3], Inf, rowObj[4:6])
    }
    tmp[i,] <- rowObj
  }

  out <- data.frame(tmp)
  names(out) <- c(names(unlist(p.prior[[1]])), "dist", "untrans")
  rownames(out) <- names(p.prior)
  return(out)
}


#' get_os Function
#'
#' To probe what OS is hosting R.
#'
#' @export
#' @examples
#' get_os()
get_os <- function()
{
  sysinf <- Sys.info()
  ostype <- .Platform$OS.type
  R.version$os

  ## Probe using Sys.info: Windows, Linux or Darwin
  if (!is.null(sysinf)) {
    os <- sysinf["sysname"]
    if (os == "Darwin") os <- "osx"
  } else {
  ## If something gets wrong with Sys.info, probe using .Platform
    os <- .Platform$OS.type
    if (grepl("^darwin", R.version$os))   os <- "osx"
    if (grepl("unix", R.version$os))      os <- "osx"
    if (grepl("linux-gnu", R.version$os)) os <- "linux"
    if (grepl("mingw32", R.version$os))   os <- "windows"
  }
  tolower(os)
}

checkHyper <- function(samples, n){
  hyper <- attr(samples, "hyper")
  nmc <- hyper$nmc
  str(hyper$h_summed_log_prior)
  str(hyper$h_log_likelihoods)
  cat("location 1 to "); cat(n, ": "); cat("\n")
  for(i in 1:n) {
    str(hyper$phi$location[,,i])
  }
  cat("scale 1 to "); cat(n, ": "); cat("\n")
  for(i in 1:n) {
    str(hyper$phi$scale[,,i])
  }

  cat("last 4 location"); cat("\n")
  for(j in (nmc-3):nmc) {
    str(hyper$phi$location[,,j])
  }

  cat("last 4 scale"); cat("\n")
  for(j in (nmc-3):nmc) {
    str(hyper$phi$scale[,,j])
  }

  thinChains  <- paste0("thins and chains: ", hyper$thin, " ", hyper$n.chains)
  cat(thinChains); cat("\n")
}

checkTheta <- function(samples, n){
  if( !is.list(samples) ) {stop("I handle multiple participants only")}
  s1 <- samples[[1]]
  nmc <- samples[[1]]$nmc
  str(s1$summed_log_prior)
  str(s1$log_likelihoods)
  cat("theta 1 to "); cat(n, ": "); cat("\n")
  for(i in 1:n) {
    str(s1$theta[,,i])
  }

  cat("last 4 theta"); cat("\n")
  for(j in (nmc-3):nmc) {
    str(s1$theta[,,j])
  }
  thinChains  <- paste0("thins and chains: ", s1$thin, " ", s1$n.chains)
  cat(thinChains); cat("\n")
}

