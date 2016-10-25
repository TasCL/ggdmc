# System functions for the DMC (Dynamic Models of Choice)
#    Functions specific to generating graphical output
#    Usually user does not need to edit

#' Plot Distributions for Each Cell
#'
#' If !is.na(C) plots density for correct and error responses for a data
#' frame with columns R (a factor) and RT, adding a boolean score column
#' for R=C. Otherwise plots each response. Can deal with NA in the RT column,
#' in which case it provides a summary of p(NA)
#' @param data.cell a data frame with only onn experimental conditoin
#' @param C a correctness column
#' @param xlim x censor range
#' @param ymax the upper bound for y axis when plotting
#' @param save.density whether to save density data
#' @param digits print how many digits
#' @param main main title for the figure
#' @param show.mean whether to show mean
#' @export
#' @importFrom graphics plot lines legend
#' @examples
#' m1 <- model.dmc(
#'   p.map     = list(a="1",v="1",z="1",d="1",sz="1",sv="1", t0="1",st0="1"),
#'   constants = c(st0=0,d=0),
#'   match.map = list(M=list(s1="r1",s2="r2")),
#'   factors   = list(S=c("s1","s2")),
#'   responses = c("r1","r2"),
#'   type      = "rd")
#'
#' p.prior <- prior.p.dmc(
#'   dists = rep("tnorm", 6),
#'   p1    = c(a=2,   v=2.5, z=0.5, sz=0.3, sv=1,  t0=0.3),
#'   p2    = c(a=0.5, v=.5,  z=0.1, sz=0.1, sv=.3, t0=0.05),
#'   lower = c(0,-5, 0, 0, 0, 0),
#'   upper = c(5, 7, 2, 2, 2, 2))
#'
#' p.vector <- c(a=1,v=1, z=0.5, sz=0.25, sv=0.2,t0=.15)
#'
#' dat1 <- simulate(m1, nsim=1e2, p.vector=p.vector)
#' mdi1 <- data.model.dmc(dat1, m1)
#'
#' ## Accuracy around 70%
#' par(mfrow=c(1,2))
#' plot_cell_density(data.cell=mdi1[mdi1$S=="s1", ], C="r1", xlim=c(0,2))
#' plot_cell_density(data.cell=mdi1[mdi1$S=="s2", ], C="r2", xlim=c(0,2))
#' par(mfrow=c(1,1))
plot_cell_density <- function(data.cell, C=NA, xlim=c(0,Inf), ymax=NA,
  save.density=FALSE, digits=3, main="", show.mean=FALSE)
{
  if (!is.factor(data.cell$R)) data.cell$R <- factor(data.cell$R)
  if (length(C)==1) C <- rep(C,dim(data.cell)[1])
  p.na <- mean(is.na(data.cell$RT))
  is.in <- !is.na(data.cell$RT)
  is.in[is.in] <- data.cell$RT[is.in]>xlim[1] & data.cell$RT[is.in]<xlim[2]
  dat <- data.cell[is.in,]
  if ( !any(is.na(C)) ) {
    if ( is.logical(C) & length(C)==dim(data.cell)[1] )
      dat$C <- C[is.in] else dat$C <- dat$R==C[is.in]
      if (length(dat$RT[dat$C])>2)
        dns.correct <- density(dat$RT[dat$C]) else dns.correct <- NULL
        if (length(dat$RT[!dat$C])>2)
          dns.error <- density(dat$RT[!dat$C]) else dns.error <- NULL
          if (is.null(dns.error) & is.null(dns.correct))
            stop("There are no densities to plot")
          acc <- mean(dat$C)
          if (!is.null(dns.correct))
            dns.correct$y <- dns.correct$y*acc*(1-p.na)
          if (!is.null(dns.error))
            dns.error$y <- dns.error$y*(1-acc)*(1-p.na)
          if (is.na(ymax)) ymax <- max(c(dns.correct$y,dns.error$y))
          if (!is.null(dns.correct)) {
            plot(dns.correct,xlab="RT",ylab="density",ylim=c(0,ymax),main=main)
            if (!is.null(dns.error)) lines(dns.error,col="red")
          } else {
            plot(dns.error,xlab="RT",ylab="density",ylim=c(0,ymax),
              main=main,col="red")
            if (!is.null(dns.correct)) lines(dns.correct)
          }
          nams <- "Accuracy ="
          ps <- round(acc,2)
          if (p.na!=0) {
            nams <- c(nams,"p(NA) =")
            ps <- c(ps,round(p.na,2))
          }
          legend("topright",paste(nams,ps),bty="n")
          legend("topright",xjust=0, inset=c(0,0.1), c("correct","error"),
            bty="n", lty=c(1,1), col=c("black","red"))
          if ( save.density ) list(correct=dns.correct,error=dns.error)
  } else {
    rs <- levels(dat$R)
    dns <- vector(mode="list",length=length(rs))
    names(dns) <- rs
    ps <- table(dat$R)/dim(dat)[1]
    for (i in rs) {
      rt <- dat$RT[dat$R==i]
      if (length(rt)>2) {
        dns[[i]] <- density(rt)
        dns[[i]]$y <- ps[i]*dns[[i]]$y*(1-p.na)
      }
    }
    ymax <- suppressWarnings(max(unlist(lapply(dns,function(x){max(x$y)}))))
    no.dns <- unlist(lapply(dns,is.null))
    if (all(no.dns))
      stop("There are no densities to plot!")
    dns1 <- dns[!no.dns]
    ltys <- c(1:length(dns1))
    plot(dns1[[1]],xlab="RT",ylab="density",ylim=c(0,ymax),lty=ltys[1],
      main=main)
    if (length(dns1)>1) for (i in 2:length(dns1)) lines(dns1[[i]],lty=ltys[i])
    nams <- paste("p(",names(dns1),") =",sep="")
    ps <- round(ps[!no.dns]*(1-p.na),2)
    if ( p.na!=0 ) {
      nams <- c(nams,"p(NA) =")
      ps <- c(ps,round(p.na,2))
      lty <- c(ltys,NA)
    }
    legend("topright",paste(nams,ps),lty=ltys,bty="n")
    if ( save.density ) dns
  }
}



#' Profile a DMC Object
#'
#' Profile a DMC model based a model data instance (ie \code{fitted}). For
#' the parameter \code{p.name}, e.g., the boundary separation \emph{a} in a DDM,
#' extract it our from a \code{p.vector} and draws the profile likelihood for
#' data and returns the maximum (on a grid of resolution \code{n.point})
#'
#' Set a range for the x axis, which is the profiled parameter. Initiate a
#' log-likelihodd vector with 0 everywhere and length n.point. keep the other
#' parmaters value fixed and change only the target parameter value to ps[i],
#' and then calculate its sum of log-likelihood. Finally, store it in the i
#' position of ll vector. This function currently works for DDM density only.
#'
#' @param fitted a DMC model data instance
#' @param p.name indicate which paramter in theta to plot. For example, in a
#' LBA model with a \code{pVec <- c(A="1",B="1",mean_v="M",sd_v="1",t0="1",
#' st0="1")}, one can assign \code{p.name <- "A"} to ask \code{profile.dmc} to
#' profile \emph{A} parameter
#' @param min.p minimal number of points to plot
#' @param max.p maximal number of points to plot
#' @param p.vector a parameter vector. Use Lognromal LBA model as an example,
#' \code{pVec <- c(A=.25, B=.35, meanlog_v.true=1, meanlog_v.false=.25,
#' sdlog_v.true=.5,t0=.2)}
#' @param n.point grid  resolution
#' @param digits print out how many digits
#' @param ... other graphical parameters (see par)
#' @export
#' @importFrom graphics plot
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
#' dat1 <- simulate(m1, nsim=1e2, p.vector=pVec)
#' mdi1 <- data.model.dmc(dat1, m1)
#'
#' ## ---------------------------
#' par(mfrow=c(2,3));
#' profile(mdi1, "a",  .1,   2, pVec)
#' profile(mdi1, "v",  .1,   2, pVec)
#' profile(mdi1, "z",  .2,  .8, pVec)
#' profile(mdi1, "sz", .1,  .9, pVec)
#' profile(mdi1, "sv", .1,   2, pVec)
#' profile(mdi1, "t0", .01, .5, pVec)
#' par(mfrow=c(1,1));
profile.dmc <- function(fitted, p.name, min.p, max.p, p.vector, n.point=100,
  digits=2, ...)
{
  if (!(p.name %in% names(p.vector))) stop("p.name not in p.vector")

  p  <- p.vector
  ps <- seq(min.p, max.p, length.out=n.point)
  ll <- numeric(n.point)
  for (i in 1:n.point)
  {
    p[p.name] <- ps[i]
    ## ll[i] <- sum(log(likelihood.dmc(p, data)))
    ll[i] <- summed_log_likelihood(p, fitted)
  }
  names(ll) <- round(ps,digits)
  plot(ps,ll,type="l",xlab=p.name,ylab="log-likelihood")
  ll[ll==max(ll)]
}

#' Plot Prior Probability Density
#'
#' Plot the i'th member of list created by \code{prior.p.dmc}.  If
#' \code{trans=TRUE} then plot on natural (logarithmic) scale using transform
#' specified in \code{attr(p.prior[[i]],"trans")}. This is akin to DMC's
#' \code{plot.prior} but uses \code{ggplot2}. Use \code{plot.priors} to plot all
#' parameters at once.
#'
#' \code{plot.prior} checks if  any of the elements in trans is \code{NA}. If
#' \code{trans} contains NAs, then it checks if \code{natural} is set
#' \code{TRUE} (by default it's \code{TRUE}). If \code{natural} is set also
#' \code{TRUE} (i.e., the user wants to do transformation), then checks what has
#' been set in the "untrans" attribute for the plotted parameter (where the i
#' argument indicates). Otherwise, we set it abitrarily as \code{identity}
#' (i.e., stay the same).
#'
#' @param i a numeral or name index indicating which parameter to plot
#' @param p.prior a prior parameter list, created usually by prior.p.dmc
#' @param xlim set the range of on x axis. This is usually the range for each
#' EAM parameter. Default is NA, letting plot facilitate to determine the range
#' automatically.
#' @param natural default TRUE
#' @param n.point default to plot 1e2
#' @param trans default NA. trans can be a scalar or vector.
#' @param main No function. A legency argument for compartibility reason
#' @param save.dat whether to save a data frame for re-plotting
#' @param ... other arguments
#' @keywords plot_prior
#' @import ggplot2
#' @export
#' @examples
#' pop.prior <- prior.p.dmc(
#'              dists = rep("tnorm", 7),
#'              p1=c(a=2,  v.f1=4, v.f2=3, z=0.5,sv=1, sz=0.3, t0=0.3),
#'              p2=c(a=0.5,v.f1=.5,v.f2=.5,z=0.1,sv=.3,sz=0.1,t0=0.05),
#'              lower=c(0,-5, -5, 0, 0, 0, 0),
#'              upper=c(5, 7,  7, 1, 2, 1, 1))
#' plot_prior("a", p.prior=pop.prior)
#' plot_prior(2,  p.prior=pop.prior)
plot_prior <- function(i, p.prior, xlim=NA, natural=TRUE, n.point=1e2,
  trans=NA, main="", save.dat=FALSE, ... )
{
  if (any(is.na(trans)))   ### Check 1 ###
  {trans <- ifelse(natural, attr(p.prior[[i]], "untrans"), "identity")}
  if (is.numeric(i)) i <- names(p.prior)[i]   ### Check 2 ###
  ## Stop the function running, if the parameters are out-of-range
  if (!(i %in% names(p.prior))) stop("Parameter not in prior")
  ## Finish checks. Calculate the data frame
  p <- p.prior[[i]]    ## Save all parameter setting to another object called p

  ## Do an educated guess for xlim
  ## xlim can be a scalar or vector. test if any of its elements is NA. If it
  if ( any(is.na(xlim)) ) {
    ## does contain NAs, we set a xlim for the user
    ## check what distribution the user wants to as her prior distribution
    ## 4 options: tnorm (default), beta_lu, gamma_l, and lnorm_l
    ## Basically, we use parameter 1 and 2 to do some (random?) arithmetic
    ## and pick (1) the bigger and smaller number (tnorm); (2) use lower and
    ## upper directly (beta_lu); (3) use lower and a different arithmetic
    ## (gamma_l); (4) use lower and another different arithmetic (lnorm_l)
    xlim <- switch(attr(p, "dist"),
      tn = c(pmax(p$lower, p[[1]]-3*p[[2]]),
        pmin(p$upper, p[[1]]+3*p[[2]])),
      tnorm    = c(pmax(p$lower, p[[1]]-3*p[[2]]),
        pmin(p$upper, p[[1]]+3*p[[2]])),
      beta_lu  = c(p$lower, p$upper),
      gamma_l  = c(p$lower, p[[1]]*p[[2]]+3*sqrt(p[[1]])*p[[2]]),
      lnorm_l  = c(p$lower, exp(p[[1]]+2*p[[2]])),
      constant = c(-10,10)
    )
  }

  ## Now we get xlim. Then we want to set all the x points for plotting.
  ## By default we plot 100 point (n.point = 1e2)
  x <- seq(xlim[1], xlim[2], length.out = n.point)
  p$x <- x
  p$log <- FALSE    ## Turn off logarithmic transformation

  ## 1. call a log/identity transform function to change x value
  ## 2. call a density function for a distribution to set y (density) value
  df <- data.frame(
    xpos  = do.call(trans, list(x = x)),
    ypos  = do.call(paste("d", attr(p,"dist"), sep=""), p),
    gpvar = rep( names(p.prior[i]), length(x)))

  if (save.dat==TRUE) {
    return(df)
  } else {
    print(ggplot(df, aes_string(x="xpos", y="ypos")) + geom_line() +
      xlab("Value")+ ylab("Density") + facet_grid(.~gpvar))
  }
}

#' Plot Prior Probability Density
#'
#' Plot each the prior density distribution for each EAM parameter and then
#' plot them all together in a canvas.
#'
#' @param p.prior a prior parameter list, created usually by prior.p.dmc
#' @param save.dat a boolean switch for save data frame, if the user wants to
#' modifiy the plotting parameter
#' @param ... other arguments
#' @keywords plot_priors
#' @import ggplot2
#' @export
#' @examples
#' pop.prior <- prior.p.dmc(
#'              dists = rep("tnorm", 7),
#'              p1=c(a=2,  v.f1=4, v.f2=3, z=0.5,sv=1, sz=0.3, t0=0.3),
#'              p2=c(a=0.5,v.f1=.5,v.f2=.5,z=0.1,sv=.3,sz=0.1,t0=0.05),
#'              lower=c(0,-5, -5, 0, 0, 0, 0),
#'              upper=c(5, 7,  7, 1, 2, 1, 1))
#'
#' view(pop.prior)
#' plot_priors(pop.prior)
#' d <- plot_priors(pop.prior, save.dat=TRUE)
#'
#' require(ggplot2)
#' p2 <- ggplot(d, aes(x = xpos, y = ypos)) + geom_line() +
#'       facet_wrap(~gpvar, scales="free") + theme_bw(base_size =14)
plot_priors <- function(p.prior, save.dat=FALSE, ...) {
  pD <- NULL
  for(j in names(p.prior))
    invisible(pD <- rbind(pD, plot_prior(i=j, p.prior=p.prior, save.dat=TRUE)))
  if (save.dat) { return(pD) } else {
    p0 <- ggplot(pD, aes_string(x = "xpos", y = "ypos")) +
        geom_line() +
        xlab("")+ ylab("")+
        facet_wrap(~gpvar, scales="free")
      print(p0)
    }
}

#' @importFrom graphics par
ppl.barplots.dmc <- function(samples,start=1,end=NA,layout=c(5,2))
  # Grid of barplots of pll for set of subjects
{
  if (!is.null(samples$theta))
    stop("This function cannot be applied to a single subject.")
  par(mfrow=layout)
  for (i in 1:length(samples))
    plot.dmc(samples,subject=i,pll.barplot=TRUE,start=start,end=end,
             main.pll=paste("Subject",i))
}


# style="pdf"; layout=NULL; pos=NULL; dname="Data";mname="Model"
#' @importFrom graphics par plot lines legend
plot.pp.dmc <- function(pp, style="pdf",layout=NULL,pos=NULL,
                        dname="Data",mname="Model")
  # pdf or cdf of data and posterior predictive fits
{
  cdf.av <- function(pp,cdf.name="cdf")
  {
      out <- outn <- outp <- pp[[1]][[cdf.name]]
      for (j in 1:length(out)) for (k in 1:length(out[[j]])) {
        # Initialize count
        outn[j][[1]][[k]] <- as.numeric(!is.null(outn[j][[1]][[k]]))
        # Intialize p
        if (!is.null(outp[j][[1]][[k]]))
          outp[j][[1]][[k]] <- as.numeric(names(outp[j][[1]][[k]]))
      }
      for (i in 2:length(pp)) {
        for (j in 1:length(out)) for (k in 1:length(out[[j]])) {
          if ( !is.null(pp[[i]]$cdf[j][[1]][[k]]) ) { # update
            if ( is.null(out[j][[1]][[k]]) ) { # first time
                out[j][[1]][[k]] <- pp[[i]]$cdf[j][[1]][[k]]
                outp[j][[1]][[k]] <- as.numeric(names(pp[[i]]$cdf[j][[1]][[k]]))
            } else { # increment values
              out[j][[1]][[k]] <- out[j][[1]][[k]] + pp[[i]]$cdf[j][[1]][[k]]
              outp[j][[1]][[k]] <-
                outp[j][[1]][[k]] + as.numeric(names(pp[[i]]$cdf[j][[1]][[k]]))
            }
            outn[j][[1]][[k]] <- outn[j][[1]][[k]] + 1 # increment count
          }
        }
      }
      for (j in 1:length(out)) for (k in 1:length(out[[j]]))
        if (!outn[j][[1]][[k]]==0) {
          out[j][[1]][[k]] <- out[j][[1]][[k]]/outn[j][[1]][[k]]
          outp[j][[1]][[k]] <- outp[j][[1]][[k]]/outn[j][[1]][[k]]
        }
      av.cdf <- out
      for ( j in 1:length(out) ) for ( k in 1:length(out[[j]]) )
        names(av.cdf[j][[1]][[k]]) <- outp[j][[1]][[k]]
      pp[[1]][[cdf.name]] <- av.cdf
      pp[[1]]
  }

  if ( !all(names(pp)[1:2]==c("pdf","cdf")) ) {
    style <- "cdf"
    dname <- paste("Average",dname)
    mname <- paste("Average",mname)
    av.cdf <- cdf.av(pp,"cdf")
    av.cdf <- cdf.av(pp,"data.cdf")
    pp <- av.cdf
  }
  if (!any(style %in% c("pdf","cdf")))
    stop("style must be either \"pdf\" or \"cdf\"")
  facs <- dimnames(pp[[1]])
  n.facs <- length(facs)
  resp <- names(pp[[1]][1][[1]])
  n.resp <- length(resp)
  n.levs <- lapply(facs,length)
  if ( is.null(layout) ) layout <- unlist(n.levs)
  if (length(layout)==1) layout <- c(1,layout) else
    layout <- layout[1:2]
  par(mfrow=layout,mar=c(2,2,0,0))
  lgnd <- NULL
  data.style <- paste("data",style,sep=".")
  for ( i in 1:length(pp[[1]]) ) if ( !all(unlist(lapply(pp[[1]][[i]],is.null))) ) {
    ok.n <- c(1:n.resp)[!unlist(lapply(pp[[1]][[i]],is.null))]
    lgnd <- unlist(lapply(pp[[style]][[i]],function(x){attr(x,"cell.name")}))
    if ( style=="pdf" ) {
      if (is.null(pos)) pos <- "topright"
      ylim <- c(0,max(
        c(unlist(lapply(pp[[style]][[i]],function(x){
            if (length(x$y)==0) NA else max(x$y)})),
          unlist(lapply(pp[[data.style]][[i]],function(x){
            if (length(x$y)==0) NA else max(x$y)}))),na.rm=TRUE))
      xlim <- c(Inf,-Inf)
      for (j in 1:length(ok.n)) {
        xlim[1] <- ifelse(min(pp[[style]][i][[1]][[ok.n[j]]]$x)<xlim[1],
                          min(pp[[style]][i][[1]][[ok.n[j]]]$x),xlim[1])
        xlim[2] <- ifelse(max(pp[[style]][i][[1]][[ok.n[j]]]$x)>xlim[2],
                          max(pp[[style]][i][[1]][[ok.n[j]]]$x),xlim[2])
      }
      for (j in 1:length(ok.n)) {
        if ( j == 1 )
          plot(pp[[style]][i][[1]][[ok.n[j]]],main="",col="red",ylim=ylim,xlim=xlim) else
            lines(pp[[style]][i][[1]][[ok.n[j]]],lty=j,col="red")
        lines(pp[[data.style]][i][[1]][[ok.n[j]]],lty=j)
      }
    } else {
      if (is.null(pos)) pos <- "topleft"
      ylim <- c(0,1)
      xlim <- c(Inf,-Inf)
      for (j in 1:length(ok.n)) {
        xlim[1] <- ifelse(min(pp[[style]][i][[1]][[ok.n[j]]])<xlim[1],
                          min(pp[[style]][i][[1]][[ok.n[j]]]),xlim[1])
        xlim[2] <- ifelse(max(pp[[style]][i][[1]][[ok.n[j]]])>xlim[2],
                          max(pp[[style]][i][[1]][[ok.n[j]]]),xlim[2])
      }
      for (j in 1:length(ok.n)) {
        y <- as.numeric(names(pp[[style]][i][[1]][[ok.n[j]]]))
        if ( j != 1 ) lines(pp[[style]][i][[1]][[ok.n[j]]],y,col="red",lty=j) else
          plot(pp[[style]][i][[1]][[ok.n[j]]],y,main="",type="l",col="red",
               ylim=ylim,xlim=xlim)
        lines(pp[[data.style]][i][[1]][[ok.n[j]]],
              as.numeric(names(pp[[data.style]][i][[1]][[ok.n[j]]])),lty=j)
      }
    }
    if (!is.na(pos)) legend(pos,c(paste(dname,lgnd),paste(mname,lgnd)),
           lty=c(1:length(ok.n),1:length(ok.n)),bty="n",
           col=c(rep("black", length(ok.n)),rep("red",length(ok.n))) )
  }
}

#' @importFrom graphics hist abline legend
plot.deviance.dmc <- function(ds=NULL,samples=NULL,digits=2,fast=TRUE,
                              main="",xlim=NA)
  ## Posterior deviance histogram
  ## TODO correct Dstats.dmc class problem
{
  if (is.null(ds)) if (is.null(samples))
    stop("Must supply samples or deviance statistics") else
      ds <- Dstats.ddm(samples, TRUE, fast)
    if (any(is.na(xlim)))
      hist(ds$D,breaks="fd",xlab="Posterior Deviance",main=main) else
        hist(ds$D,breaks="fd",xlab="Posterior Deviance",main=main,xlim=xlim)
    abline(v=ds$Dmean,lty=2)
    abline(v=ds$meanD,lty=3)
    pds <- pd.dmc(ds)
    legend("topright",c(paste("NP =",ds$np),
                        paste("NPmean =",round(pds$Pmean,digits)),
                        paste("NPmin =",round(pds$Pmin,digits)),
                        paste("NPvar =",round(pds$Pvar,digits)),
                        "D(mean theta)","mean(D)"),
           bty="n",lty=c(NA,NA,NA,NA,2:3))
}



### Fixed Effects

#' @importFrom graphics arrows
plotSpar.dmc <- function(est,p.name, a.len=.05)
  # plots ordered cis for p.name, est produced by summary.dmc
{
  lo <- unlist(lapply(est,function(x){x$quantiles[p.name,c("2.5%")]}))
  med <- unlist(lapply(est,function(x){x$quantiles[p.name,c("50%")]}))
  hi <- unlist(lapply(est,function(x){x$quantiles[p.name,c("97.5%")]}))
  ord <- order(med)
  n <- length(med)
  plot(1:n,med[ord],ylim=c(min(lo),max(hi)),ylab=p.name,xlab="Subject",pch=16,
    main="Medians and 95% credible intervals")
  arrows(1:n,med[ord],1:n,lo[ord],angle=90, length=a.len)
  arrows(1:n,med[ord],1:n,hi[ord],angle=90, length=a.len)
}


### Hierachical

#' @importFrom graphics plot
h.profile.dmc <- function(p.name,p.num,min.p,max.p,ps,p.prior,n.point=100,
                          digits=3,ylim=NA)
  # for parameter p.name at position p.num (1 or 2) in p.prior given subject
  # parameters ps draws a likelihood profile and returns the maximum (on a grid
  # of resolution n.point)
{
  if (!(p.name %in% dimnames(ps)[[2]]))
    stop("p.name not in ps")
  if (!(p.num %in% 1:2))
    stop("p.num must be either 1 or 2")
  pps <- seq(min.p,max.p,length.out=n.point)
  ll <- numeric(n.point)
  pop <- p.prior[[p.name]]
  pop$x <- ps[,p.name]
  for (i in 1:n.point)
  {
    pop[[p.num]] <- pps[i]
    ll[i] <- sum(do.call(paste("d",attr(pop,"dist"),sep=""),pop))
  }
  names(ll) <- round(pps,digits)
  p.name <- paste(p.name,names(pop)[p.num],sep=".")
  if (any(is.na(ylim)))
    plot(pps,ll,type="l",xlab=paste(p.name,p.num,sep="."),
         ylab="log-likelihood") else
           plot(pps,ll,type="l",xlab=paste(p.name,p.num,sep="."),
                ylab="log-likelihood",ylim=ylim)
  ll[ll==max(ll)]
}

#' @importFrom stats cor
cor.plausible <- function(hsamples,p.name,cv,plot=FALSE,
                          xlab="r",ylab="Density",main=p.name,
                          fun=NULL,...)
  # Correlations for each interation and chain with of p.name with cv$p.name
{
  if ( !is.null(fun) ) subject.pars <- do.call(rbind,lapply(hsamples,function(x){
    pars <- aperm(x$theta,c(1,3,2));
    pars <- matrix(as.vector(pars),ncol=dim(pars)[3]);
    dimnames(pars)[[2]] <- dimnames(x$theta)[[2]];
    apply(pars,1,fun)
  })) else subject.pars <- do.call(rbind,lapply(hsamples,function(x){
    as.vector(x$theta[,p.name,])}))
  out <- apply(subject.pars,2,function(x){cor(x,cv[[p.name]])})
  if (plot) {
    dns <- density(out)
#    dns$y <- dns$y*diff(range(dns$x))/length(dns$x)
    plot(dns,xlab=xlab,ylab=ylab,main=main,...)
  }
  invisible(out)
}

