###############################################################################
##### Functions taken from rtdists to enable accumulator vectorized t0 ########
############  WILL NOT BE NECESSARY AFTER NEXT UPDATE OF RTDISTS ##############
############### NB: Use getAnywhere(name) for rtdists functions ###############
###############################################################################

#######  Standard model ###

make.r <- function (drifts, b, A, n_v, t0, st0 = 0, n)
{
    drifts[drifts < 0] <- 0
    starts <- matrix(runif(min = 0, max = A, n = n * n_v), nrow = n_v)
    ttf <- t0 + (b - starts)/drifts
    resp <- apply(ttf, 2, which.min)
    rt <- ttf[cbind(resp,1:n)]
    if (st0[1]>0) rt <- rt + runif(min = 0, max = st0[1], n = n)
    bad <- !is.finite(rt)
    if (any(bad)) {
        warning(paste(sum(bad), "infinite RTs removed and less than",
            n, "rts returned"))
        resp <- resp[!bad]
        rt <- rt[!bad]
    }
    data.frame(rt = rt, response = resp)
}


rlba.norm <- function (n,A,b,t0,mean_v,sd_v,st0=0,posdrift=TRUE)
{
    if (any(b < A)) stop("b cannot be smaller than A!")
    n_v <- length(mean_v)
    if (posdrift) {
      drifts <- matrix(ggdmc::rtnorm(n = n * n_v, mean = mean_v, sd = sd_v,
        lower = 0, checks=TRUE), nrow = n_v)
    } else {
      drifts <- matrix(ggdmc::rtnorm(n = n * n_v, mean = mean_v, sd = sd_v),
        nrow = n_v)
    }
    make.r(drifts = drifts, b = b, A = A, n_v = n_v, t0 = t0, st0 = st0, n = n)
}


dlba.norm.core <- function (t,A,b,mean_v,sd_v,posdrift=TRUE,robust=FALSE,nn)
# like dlba_norm_core but t0 dealt with outside (removed from t), and
# default added to nn
{

  pnormP <- function (x, mean = 0, sd = 1, lower.tail = TRUE)
    ifelse(abs(x) < 7, pnorm(x, mean = mean, sd = sd, lower.tail = lower.tail),
      ifelse(x < 0, 0, 1))

  dnormP <- function (x, mean = 0, sd = 1)
    ifelse(abs(x) < 7, dnorm(x, mean = mean, sd = sd), 0)

    if (robust) {
        pnorm1 <- pnormP
        dnorm1 <- dnormP
    }
    else {
        pnorm1 <- pnorm
        dnorm1 <- dnorm
    }
    if (posdrift)
        denom <- pmax(pnorm1(mean_v/sd_v), 1e-10)
    else denom <- rep(1, nn)
    if (any(A < 1e-10)) {
        A_small <- A < 1e-10
        out_A <- pmax(0, ((b[A_small]/t[A_small]^2) *
                            dnorm1(b[A_small]/t[A_small],
            mean_v[A_small], sd = sd_v[A_small]))/denom[A_small])
        zs <- t[!A_small] * sd_v[!A_small]
        zu <- t[!A_small] * mean_v[!A_small]
        chiminuszu <- b[!A_small] - zu
        chizu <- chiminuszu/zs
        chizumax <- (chiminuszu - A[!A_small])/zs
        out_o <- pmax(0, (mean_v[!A_small] * (pnorm1(chizu) -
            pnorm1(chizumax)) + sd_v[!A_small] * (dnorm1(chizumax) -
            dnorm1(chizu)))/(A[!A_small] * denom[!A_small]))
        out <- numeric(nn)
        out[!A_small] <- out_o
        out[A_small] <- out_A
        return(out)
    }
    else {
        zs <- t * sd_v
        zu <- t * mean_v
        chiminuszu <- b - zu
        chizu <- chiminuszu/zs
        chizumax <- (chiminuszu - A)/zs
        return(pmax(0, (mean_v * (pnorm1(chizu) - pnorm1(chizumax)) +
            sd_v * (dnorm1(chizumax) - dnorm1(chizu)))/(A * denom)))
    }
}


plba.norm.core <- function (t,A,b,mean_v,sd_v,posdrift=TRUE,robust=FALSE,nn)
# like plba_norm_core but t0 dealt with outside (removed from t), and
# default added to nn
{

    pnormP <- function (x, mean = 0, sd = 1, lower.tail = TRUE)
      ifelse(abs(x) < 7, pnorm(x, mean = mean, sd = sd, lower.tail=lower.tail),
        ifelse(x < 0, 0, 1))

    dnormP <- function (x, mean = 0, sd = 1)
      ifelse(abs(x) < 7, dnorm(x, mean = mean, sd = sd), 0)

    if (robust) {
        pnorm1 <- pnormP
        dnorm1 <- dnormP
    }
    else {
        pnorm1 <- pnorm
        dnorm1 <- dnorm
    }
    if (posdrift)
        denom <- pmax(pnorm1(mean_v/sd_v), 1e-10)
    else denom <- 1
    if (any(A < 1e-10)) {
        A_small <- A < 1e-10
        out_A <- pmin(1, pmax(0, (pnorm1(b[A_small]/t[A_small],
            mean = mean_v[A_small], sd = sd_v[A_small],
            lower.tail = FALSE))/denom[A_small]))
        zs <- t[!A_small] * sd_v[!A_small]
        zu <- t[!A_small] * mean_v[!A_small]
        chiminuszu <- b[!A_small] - zu
        xx <- chiminuszu - A[!A_small]
        chizu <- chiminuszu/zs
        chizumax <- xx/zs
        tmp1 <- zs * (dnorm1(chizumax) - dnorm1(chizu))
        tmp2 <- xx * pnorm1(chizumax) - chiminuszu * pnorm1(chizu)
        out_o <- pmin(pmax(0, (1 + (tmp1 + tmp2)/A[!A_small])/denom[!A_small]),
            1)
        out <- numeric(length(mean_v))
        out[!A_small] <- out_o
        out[A_small] <- out_A
        return(out)
    }
    else {
        zs <- t * sd_v
        zu <- t * mean_v
        chiminuszu <- b - zu
        xx <- chiminuszu - A
        chizu <- chiminuszu/zs
        chizumax <- xx/zs
        tmp1 <- zs * (dnorm1(chizumax) - dnorm1(chizu))
        tmp2 <- xx * pnorm1(chizumax) - chiminuszu * pnorm1(chizu)
        return(pmin(pmax(0, (1 + (tmp1 + tmp2)/A)/denom), 1))
    }
}


dlba.norm <- function (t, A, b, mean_v, sd_v, posdrift = TRUE, robust = FALSE)
{
    nn <- length(t)
    A <- rep(A, length.out = nn)
    b <- rep(b, length.out = nn)
    mean_v <- rep(mean_v, length.out = nn)
    sd_v <- rep(sd_v, length.out = nn)
    if (any(b < A)) # b cannot be smaller than A!
        return(rep(0,nn))
    dlba.norm.core(t = t, A = A, b = b, mean_v = mean_v,
        sd_v = sd_v, posdrift = posdrift, robust = robust, nn = nn)
}


plba.norm <- function (t, A, b, mean_v, sd_v, posdrift = TRUE, robust = FALSE)
{
    nn <- length(t)
    A <- rep(A, length.out = nn)
    b <- rep(b, length.out = nn)
    mean_v <- rep(mean_v, length.out = nn)
    sd_v <- rep(sd_v, length.out = nn)
    if (any(b < A)) # b cannot be smaller than A!
        return(rep(0,nn))
    plba.norm.core(t = t, A = A, b = b, mean_v = mean_v,
        sd_v = sd_v, posdrift = posdrift, robust = robust, nn = nn)
}


n1PDFfixedt0.norm=function(dt,A,b,mean_v,sd_v,
                           posdrift=TRUE,robust = FALSE)
# Generates defective PDF for responses on node= 1
# dt (decison time) is a matrix with length(mean_v) rows, one row for
# each accumulator to allow for different start times
{

  dt[1,] <- dlba.norm(dt[1,],A=A[1],b=b[1],mean_v=mean_v[1],sd_v=sd_v[1],
                      posdrift=posdrift,robust=robust)
  if (length(mean_v)>1) for (i in 2:length(mean_v))
    dt[1,] <- dt[1,]*(1-plba.norm(dt[i,],
      A=A[i],b=b[i],mean_v=mean_v[i],sd_v=sd_v[i],
      posdrift=posdrift,robust=robust))
  dt[1,]
}

################################################################################
############################ ADDITIONAL MODELS #################################
################################################################################

# IMPLEMENTED WITH standard t0
# LBA Gamma non-decision
# LBA Lognormal non-decision

# IMPLEMENTED WITH ACCUMULATOR VECTORIZED t0
# LNR n-choice
# LNR go-nogo
# LBA go-nogo
# LNR stop signal
# ExGaussian stop signal.

#####################  Gamma and Lognormal Non-decison time ## ----

{
n1PDFfixedt0.norm.old <- function(dt,A,b,mean_v,sd_v,
                           posdrift=TRUE,robust = FALSE)
# Generates defective PDF for responses on node= 1
# dt (decison time) is a vector with t0 removed
{
  out <- dlba.norm(t=dt,A=A[1],b=b[1],mean_v=mean_v[1],sd_v=sd_v[1],
                  posdrift=posdrift,robust=robust)
  if (length(mean_v)>1) for (i in 2:length(mean_v))
    out <- out*(1-plba.norm(t=dt,
      A=A[i],b=b[i],mean_v=mean_v[i],sd_v=sd_v[i],
      posdrift=posdrift,robust=robust))
  out
}
}

####################### LBA Gamma non-decision ----
{
rlba.norm.gamma <- function(n,A,b,t0,mean_v,sd_v,shape,scale,posdrift=TRUE)
{
    if (any(b < A))
        stop("b cannot be smaller than A!")
    n_v <- length(mean_v)
    if (posdrift)
      drifts <- matrix(rtnorm(n = n * n_v, mean = mean_v, sd = sd_v,
                       lower = 0), nrow = n_v) else
      drifts <- matrix(rnorm(n = n * n_v, mean = mean_v, sd = sd_v), nrow = n_v)
    out <- make.r(drifts=drifts,b=b,A=A,n_v=n_v,t0=t0,n = n)
    if (shape>0 & scale>0)
      cbind.data.frame(RT=out$rt + rgamma(n,shape=shape[1],scale=scale[1]),
                       R=factor(out$response)) else
      cbind.data.frame(RT=out$rt,R=factor(out$response))
}


n1PDF.norm.gamma <- function (t,A,b,mean_v,sd_v,t0,shape,scale,
                           posdrift=TRUE,robust = FALSE)
{
  dt <- pmax(t-t0[1], 0)
  if (shape[1] == 0 | scale[1] == 0)
    return(n1PDFfixedt0.norm.old(dt=dt,A=A,b=b,mean_v=mean_v,sd_v=sd_v,
                              posdrift=posdrift,robust = robust))
    else {
      tmpf <- function(tau,ti0,A,b,mean_v,sd_v,scale,shape)
        n1PDFfixedt0.norm.old(tau,A,b,mean_v,sd_v)*
          dgamma(ti0-tau,shape=shape,scale=scale)
    outs <- numeric(length(dt))
    for (i in c(1:length(outs))) {
      tmp <- try(integrate(f=tmpf,
        lower=0,upper=dt[i],ti0=dt[i],
        A=A,b=b,mean_v=mean_v,sd_v=sd_v,
        scale=scale[1],shape=shape[1])$value,silent=T)
      if (is.numeric(tmp))
        outs[i] <- tmp else
        outs[i] <- 0
    }
    return(outs)
  }
}


# # Check
# n=1e6
# shape=0; scale=0
# shape=1; scale=1
# t0=.2; A=c(.5,.5); b=c(1,1); mean_v=c(1,.5); sd_v=c(1,1)
# sim <- rlba.norm.gamma(n=n,A=A,b=b,mean_v=mean_v,sd_v=sd_v,
#                        shape=shape,scale=scale,t0=t0)
# dns <- plot.cell.density(sim,C=1,xlim=c(0,7),save.density=TRUE)
# # n1PDF check corret
# d <- n1PDF.norm.gamma(dns$correct$x,A=A,b=b,mean_v=mean_v,sd_v=sd_v,t0=t0,
#                       shape=shape,scale=scale)
# # red=black?
# plot(dns$correct$x,dns$correct$y,type="l",ylab="density",xlab="RT")
# lines(dns$correct$x,d,col="red")
# # n1PDF check, error
# d <- n1PDF.norm.gamma(dns$error$x,A=A[2:1],b=b[2:1],
#                       mean_v=mean_v[2:1],sd_v=sd_v[2:1],t0=t0,
#                       shape=shape,scale=scale)
# lines(dns$error$x,dns$error$y,lty=2)
# lines(dns$error$x,d,col="red",lty=2)
}

####### LBA Gamma non-decision, MR linear angle effects (AS) on shape  ----
{
rlba.norm.gammaMR <- function(n,A,b,t0,mean_v,sd_v,shape,scale,AS,posdrift=TRUE)
{
  shape <- shape*AS
  if (any(b < A))
        stop("b cannot be smaller than A!")
    n_v <- length(mean_v)
    if (posdrift)
      drifts <- matrix(rtnorm(n = n * n_v, mean = mean_v, sd = sd_v,
                       lower = 0), nrow = n_v) else
      drifts <- matrix(rnorm(n = n * n_v, mean = mean_v, sd = sd_v), nrow = n_v)
    out <- make.r(drifts=drifts,b=b,A=A,n_v=n_v,t0=t0,n = n)
    if (shape>0 & scale>0)
      cbind.data.frame(RT=out$rt + rgamma(n,shape=shape[1],scale=scale[1]),
                       R=factor(out$response)) else
      cbind.data.frame(RT=out$rt,R=factor(out$response))
}


n1PDF.norm.gammaMR <- function (t,A,b,mean_v,sd_v,t0,shape,scale,AS,
                           posdrift=TRUE,robust = FALSE)
{
  shape <- shape*AS
  dt <- pmax(t-t0[1], 0)
  if (shape[1] == 0 | scale[1] == 0)
    return(n1PDFfixedt0.norm.old(dt=dt,A=A,b=b,mean_v=mean_v,sd_v=sd_v,
                              posdrift=posdrift,robust = robust))
    else {
      tmpf <- function(tau,ti0,A,b,mean_v,sd_v,scale,shape)
        n1PDFfixedt0.norm.old(tau,A,b,mean_v,sd_v)*
          dgamma(ti0-tau,shape=shape,scale=scale)
    outs <- numeric(length(dt))
    for (i in c(1:length(outs))) {
      tmp <- try(integrate(f=tmpf,
        lower=0,upper=dt[i],ti0=dt[i],
        A=A,b=b,mean_v=mean_v,sd_v=sd_v,
        scale=scale[1],shape=shape[1])$value,silent=T)
      if (is.numeric(tmp))
        outs[i] <- tmp else
        outs[i] <- 0
    }
    return(outs)
  }
}


# # Check
# n=1e6
# shape=1; scale=1
# AS=0
# AS=1
# t0=.2; A=c(.5,.5); b=c(1,1); mean_v=c(1,.5); sd_v=c(1,1)
# sim <- rlba.norm.gammaMR(n=n,A=A,b=b,mean_v=mean_v,sd_v=sd_v,
#                        shape=shape,scale=scale,t0=t0,AS=AS)
# dns <- plot.cell.density(sim,C=1,xlim=c(0,7),save.density=TRUE)
# # n1PDF check corret
# d <- n1PDF.norm.gammaMR(dns$correct$x,A=A,b=b,mean_v=mean_v,sd_v=sd_v,t0=t0,
#                       shape=shape,scale=scale,AS=AS)
# # red=black?
# plot(dns$correct$x,dns$correct$y,type="l",ylab="density",xlab="RT")
# lines(dns$correct$x,d,col="red")
# # n1PDF check, error
# d <- n1PDF.norm.gammaMR(dns$error$x,A=A[2:1],b=b[2:1],
#                       mean_v=mean_v[2:1],sd_v=sd_v[2:1],t0=t0,
#                       shape=shape,scale=scale,AS=AS)
# lines(dns$error$x,dns$error$y,lty=2)
# lines(dns$error$x,d,col="red",lty=2)
}

####################### LBA Lognormal non-decision ----
{
rlba.norm.lnorm <- function(n,A,b,t0,mean_v,sd_v,meanlog,sdlog,posdrift=TRUE)
{
    if (any(b < A))
        stop("b cannot be smaller than A!")
    n_v <- length(mean_v)
    if (posdrift)
        drifts <- matrix(rtnorm(n = n * n_v, mean = mean_v, sd = sd_v,
            lower = 0), nrow = n_v)
    else drifts <- matrix(rnorm(n = n * n_v, mean = mean_v, sd = sd_v),
        nrow = n_v)
    out <- make.r(drifts=drifts,b=b,A=A,n_v=n_v,t0=t0,n = n)
    if (sdlog>0)
      cbind.data.frame(RT=out$rt + rlnorm(n,meanlog=meanlog,sdlog=sdlog),
                       R=factor(out$response)) else
      cbind.data.frame(RT=out$rt,R=factor(out$response))
}

n1PDF.norm.lnorm <- function (t,A,b,mean_v,sd_v,t0,meanlog,sdlog,
                           posdrift=TRUE,robust = FALSE)
{
  dt <- pmax(t-t0[1], 0)
  if (sdlog[1] == 0)
    return(n1PDFfixedt0.norm.old(dt=dt,A=A,b=b,mean_v=mean_v,sd_v=sd_v,
                              posdrift=posdrift,robust = robust))
    else {
      tmpf <- function(tau,ti0,A,b,mean_v,sd_v,meanlog,sdlog)
        n1PDFfixedt0.norm.old(tau,A,b,mean_v,sd_v)*
          dlnorm(ti0-tau,meanlog=meanlog,sdlog=sdlog)
    outs <- numeric(length(dt))
    for (i in c(1:length(outs))) {
      tmp <- try(integrate(f=tmpf,
        lower=0,upper=dt[i],ti0=dt[i],
        A=A,b=b,mean_v=mean_v,sd_v=sd_v,
        meanlog=meanlog[1],sdlog=sdlog[1])$value,silent=T)
      if (is.numeric(tmp))
        outs[i] <- tmp else
        outs[i] <- 0
    }
    return(outs)
  }
}

# # Check
# n=1e6
# sdlog=1; meanlog=0
# t0=.2; A=c(.5,.5); b=c(1,1); mean_v=c(1,.5); sd_v=c(1,1)
# sim <- rlba.norm.lnorm(n=n,A=A,b=b,mean_v=mean_v,sd_v=sd_v,
#                        meanlog=meanlog,sdlog=sdlog,t0=t0)
# dns <- plot.cell.density(sim,C=1,xlim=c(0,7),save.density=TRUE)
# # n1PDF check, correct
# d <- n1PDF.norm.lnorm(dns$correct$x,A=A,b=b,mean_v=mean_v,sd_v=sd_v,t0=t0,
#                       meanlog=meanlog,sdlog=sdlog)
# # red=black?
# plot(dns$correct$x,dns$correct$y,type="l",ylab="density",xlab="RT")
# lines(dns$correct$x,d,col="red")
# # n1PDF check, error
# d <- n1PDF.norm.lnorm(dns$error$x,A=A[2:1],b=b[2:1],
#                       mean_v=mean_v[2:1],sd_v=sd_v[2:1],t0=t0,
#                       meanlog=meanlog,sdlog=sdlog)
# lines(dns$error$x,dns$error$y,lty=2)
# lines(dns$error$x,d,col="red",lty=2)
}

####################### LNR n-choice ----
{
rlnr <- function (n, meanlog, sdlog, t0, st0 = 0)
# Race among n_acc accumulators, mealnlog and sdlog can be n_acc length
# vectors or n_acc x n matrices. t0 can be
# a) a scalar, b) a vector of length number of accumulators or
# c) a matrix with 1 row per accumulator, when start times differ on each trial
# st0, range of non-decison time variability, must be a scalar, as the same
# variability is assumed in a common encoding/production stage

{
    n_acc <- ifelse(is.null(dim(meanlog)),length(meanlog),dim(meanlog)[1])
    dt <- matrix(rlnorm(n = n*n_acc, meanlog = meanlog, sdlog = sdlog),
        nrow = n_acc) + t0
    winner <- apply(dt,2,which.min)
    if (st0[1]==0) data.frame(RT=dt[cbind(winner,1:n)],R=winner) else
      data.frame(RT=dt[cbind(winner,1:n)]+runif(n,0,st0[1]),R=winner)
}

n1PDFfixedt0.lnr=function(dt,meanlog,sdlog)
# Generates defective PDF for responses first among n_acc accumulator at
# dt (decison time), a matrix with one row for each accumulator (allowing for
# different start times per accumulator)
{

  n_acc <- ifelse(is.null(dim(meanlog)),length(meanlog),dim(meanlog)[1])
  if (!is.matrix(meanlog)) meanlog <- matrix(rep(meanlog,dim(dt)[2]),nrow=n_acc)
  if (!is.matrix(sdlog))     sdlog <- matrix(rep(  sdlog,dim(dt)[2]),nrow=n_acc)
  # winner
  dt[1,] <- dlnorm(dt[1,],meanlog[1,],sdlog[1,])
  # loosers
  if (dim(meanlog)[1]>1) for (i in 2:dim(meanlog)[1])
    dt[1,] <- dt[1,]*plnorm(dt[i,],meanlog[i,],sdlog[i,],lower.tail=FALSE)
  dt[1,]
}

n1PDF.lnr <- function(dt,meanlog,sdlog,t0,st0=0)
# dt (decision time) is a vector, meanlog and sdlog have same length = number of
# accumulators, t0 is the lower bound of non-decision time, it can be:
# a) a scalar, b) a vector of length number of accumulators or
# c) a matrix with 1 row per accumulator, when start times differ on each trial
# st0, range of non-decison time variability, must be a scalar, as the same
# variability is assumed in a common encoding/production stage
{

  # NOTE: No need to flag negative dt in dt-t0 as plnorm/dlnorm return 0

  n_acc <- ifelse(is.null(dim(meanlog)),length(meanlog),dim(meanlog)[1])
  if (!is.matrix(meanlog)) meanlog <- matrix(rep(meanlog,length(dt)),nrow=n_acc)
  if (!is.matrix(sdlog))     sdlog <- matrix(rep(  sdlog,length(dt)),nrow=n_acc)
  if ( st0 < 1e-3 ) # smaller values can cause numerical integration problems
    return(n1PDFfixedt0.lnr(meanlog=meanlog,sdlog=sdlog,dt=matrix(
      pmax(rep(dt,each=n_acc)-t0,0),nrow=n_acc))) else
  {

    integrate.f <- function(dt,meanlog,sdlog,t0,st0)
      n1PDFfixedt0.lnr(meanlog=meanlog,sdlog=sdlog,dt=matrix(
        pmax(rep(dt,each=length(meanlog))-t0,0),nrow=length(meanlog)))/st0

    outs <- numeric(length(dt))
    for (i in 1:length(outs)) {
      tmp <- try(integrate(f=integrate.f,
        lower=dt[i]-st0[1],upper=dt[i],
        meanlog=meanlog[,i],sdlog=sdlog[,i],t0=t0,st0=st0[1])$value,silent=TRUE)
      if (is.numeric(tmp))
        outs[i] <- tmp else
        outs[i] <- 0
    }
    return(outs)
  }
}


# # Check
# n=1e5
# meanlog=c(.5,.75,1); sdlog=c(1,1,1)
# n_acc <- length(meanlog)
#
# # check scalar t0
# t0=.2
# t0=c(.2,1,1)
#
# # check t0 noise
# st0=1
# st0=0
#
# sim <- rlnr(n=n,meanlog,sdlog,t0,st0)
# dns <- plot.cell.density(sim,xlim=c(0,7),save.density=TRUE)
#
# # Check matrix forms work
# meanlog <- matrix(rep(meanlog,length(dns[[ichar]]$x)),nrow=n_acc)
# sdlog <- matrix(rep(sdlog,length(dns[[ichar]]$x)),nrow=n_acc)
# t0 <- matrix(rep(t0,length.out=n_acc*length(dns[[ichar]]$x)),nrow=n_acc)
#
# check_n1 <- TRUE   # n1PDF check
# check_n1 <- FALSE  # n1PDFfixedt0.lnr check (ignores st0)
#
# for (i in 1:n_acc) {
#   ichar <- as.character(i)
#   indx <- c(i,c(1:n_acc)[-i])
#   if (!is.matrix(meanlog)) meanlogi <- meanlog[indx] else meanlogi <- meanlog[indx,]
#   if (!is.matrix(sdlog))     sdlogi <- sdlog[indx] else     sdlogi <- sdlog[indx,]
#   if (!is.matrix(t0)) {
#     if (length(t0)!=1) t0i <- t0[indx] else t0i <- t0
#   } else t0i <- t0[indx,]
#   if (check_n1) d <- n1PDF.lnr(dt=dns[[ichar]]$x,meanlogi,sdlogi,t0i,st0=st0) else
#     d <- n1PDFfixedt0.lnr(meanlog=meanlogi,sdlog=sdlogi,
#       dt=matrix(pmax(rep(dns[[ichar]]$x,each=n_acc)-t0i,0),nrow=n_acc))
#   if (i!=1) lines(dns[[ichar]]$x,dns[[ichar]]$y,lty=i) else
#     plot(dns[[ichar]]$x,dns[[ichar]]$y,type="l",xlab="RT",ylab="Density")
#   lines(dns[[ichar]]$x,d,col="red",lty=i)
# }
}

######################## PNRWald ----

# n-choice postivie (zero truncated) normal rate trial to trial ~N(v,sv) Wald
#    race, with t0, v, sv, a (boundary) parameterixaiton

# require("statmod") invgauss random, pdf, cdf etc.
# This has parameterizaiton mean (expected value) and dispersion = 1/shape. In
# the accumulator interpritaiton dispersion is scaled relative to threshold (a)
# so without loss of generality we can set diffusion coefficeint to 1 and so
# shape = a^2. Using the notation v = mean rate, expected value (mean) = a/v
# and varaince = mean^3 / a^2 = a / v^3

### Single accumulator model

rPNRWald <- function(n,a,v,sv,t0)
  # random function for single acumulator
{
    drifts <- matrix(rtnorm(n = n , mean = v, sd = sv, lower = 0),
                     nrow = length(v))
    t0+rinvgauss(n,mean=a/drifts,shape=a^2)
}


dPNRWald <- function(t,a,v,sv,posdrift=TRUE)
  # density for single accumulator
{
  # Convert to Desmond and Yang's notation
  d <- v/a # mean of delta
  l <- a^2 # lambda assuming diffusive variance is 1
  v <- sv^2 # variance of delta
  sqrt(l/(2*pi*t^3*(v*t+1)))*exp(-(l*(d*t-1)^2)/(2*t*(v*t+1)))*
  pnorm((v+d)*sqrt(l/(t*v^2+v)))/pnorm(d*sqrt(l/v)) # normalize
}


# # Check density and normalization
# n=1e6; a=3; v=2; sv=1
# integrate(dPNRWald,lower=0,upper=Inf,a=a,v=v,sv=sv)
# sim <- rPNRWald(n=n,a=a,v=v,sv=sv,t0=0)
# bad <- sim>10; mean(bad) # if this is large the check isnt valid
# dns <- density(sim[!bad])
# x <- dns$x[dns$x>0]
# d <- dPNRWald(x,a=a,v=v,sv=sv)
# plot(dns)
# lines(x,d,col="red")

pPNRWald <- function(t,a,v,sv)
  # cumulative density for single accumulator, NB: direct integration of
  # dPNRWald fails! This method using the analytic for the un-mixed
  # Wald CDF works much better
  ## @importFrom statmod pinvgauss
{
  if ( length(c(a,v,sv))>3 )
    stop("pPNRWald can only handle scalar parameters")
  ifun <- function(x,a,v,sv,t)
    dtnorm(x,v,sv,lower=0)*pinvgauss(t,mean=a/x,shape=a^2)
  for (i in 1:length(t)) t[i] <-
      integrate(ifun,lower=0,upper=Inf,a=a,v=v,sv=sv,t=t[i])$value
#   t[t<0] <- 0; t[t>1] <- 1 # doesnt seem necessary
  t
}


# # Check cumulative density
# n=1e6; a=.1; v=.2; sv=1
# sim <- rPNRWald(n=n,a=a,v=v,sv=sv,t0=0)
# probs=1:990/1000
# qs <- quantile(sim,probs=probs)
# cd <- pPNRWald(qs,a=a,v=v,sv=sv)
# plot(qs,probs,type="l")
# lines(qs,cd,col="red")


### Race model

rPNRWaldRace <- function(n,a,v,sv,t0)
  # random function for PNRWald race, if all accumualtors have non-poistive
  # rates race deleted and warning returned.
{
    n_v <- length(v)
    drifts <- matrix(rtnorm(n=n*n_v,mean=v,sd=sv,lower=0),nrow=n_v)
    ttf <- matrix(t0 + rinvgauss(n*n_v,mean=a/drifts,shape=a^2),nrow=n_v)
    resp <- apply(ttf, 2, which.min)
    data.frame(RT = ttf[cbind(resp,1:n)], R = factor(apply(ttf, 2, which.min)))
}


n1PNRWald <- function(dt,a,v,sv,t0,posdrift=TRUE)
# Generates defective PDF for responses on node=1dt (decison time) is a vector of times
{
  # Some protection for negative dt values, remove if not needed
  dt <- dt-t0
  bad <- dt <= 0
  d <- numeric(length(dt))
  d[!bad] <- dPNRWald(dt[!bad],a=a[1],v=v[1],sv=sv[1])
  for (i in 2:length(v))
    d[!bad] <- d[!bad]*(1-pPNRWald(dt[!bad],a=a[i],v=v[i],sv=sv[i]))
  d
}


# # Check
# n=1e5
# v=c(4,0.5); sv=c(1,1); a=c(2,2); t0=0.5
# sim <- rPNRWaldRace(n=n,a=a,v=v,sv=sv,t0=t0)
# par(mfrow=c(1,2))
# dns <- plot.cell.density(sim,C=1,xlim=c(0,4),save.density=TRUE)
# # NB: Direclty integrating again appears very inaccurate
# integrate(n1PNRWald,lower=t0,upper=Inf,a=a,v=v,sv=sv,t0=t0)$value
# dt=dns$correct$x
# d <- n1PNRWald(dt,a=a,v=v,sv=sv,t0=t0)
# # red=black?
# plot(dns$correct$x,dns$correct$y,type="l")
# lines(dns$correct$x,d,col="red")

####################### Robust integrator used in GNG and SS code ----

# fix me
my.integrate <- function(...,big=10)
# Avoids but in integrate upper=Inf that uses only 1  subdivision
# Use of  big=10 is arbitary ...
{
  out <- try(integrate(...,upper=Inf),silent=TRUE)
  if (class(out)=="try-error") 0 else
  {
    cat("Subdivisions\n")
    print(out$subdivisions)

    if (out$subdivisions==1)
    {
      out <- try(integrate(...,upper=big),silent=TRUE)
      if (class(out)=="try-error") 0 else
      {
         if (out$subdivisions==1) 0 else out$value
      }
    } else out$value
  }
}


####################### LNR go-nogo ----
{
rlnrgng <- function (n, meanlog, sdlog, t0, st0 = 0)
# Race among length(meanlog) accumulators, first of which is a no-go
# accumulator. For trials with winning first accumulator RT set to NA.
{
   out <- rlnr(n,meanlog,sdlog,t0,st0)
   out[out$R==1,"RT"] <- NA
   out$R <- factor(as.numeric(out$R))
   data.frame(out)
}


n1PDFfixedt0.lnrgng=function(dt,meanlog,sdlog)
# Same as n1PDFfixedt0.lnr except dt=NA done by integration
{

  stopfn <- function(t,meanlog,sdlog)
  {
    n1PDFfixedt0.lnr(
      matrix(rep(t,each=length(meanlog)),nrow=length(meanlog)),
      meanlog,sdlog
    )
  }

  n.trials <- dim(dt)[2]
  out <- numeric(n.trials)
  is.stop <- is.na(dt[1,])
  dt[1,!is.stop] <- n1PDFfixedt0.lnr(dt[,!is.stop,drop=F],meanlog,sdlog)
  # tmp <- try(integrate(f=stopfn,lower=0,upper=Inf,
  #     meanlog=meanlog,sdlog=sdlog)$value,silent=TRUE)
  # if (!is.numeric(tmp)) tmp <- 0
  tmp <- my.integrate(f=stopfn,lower=0,meanlog=meanlog,sdlog=sdlog)
  dt[1,is.na(dt[1,])] <- tmp
  dt[1,]
}

n1PDF.lnrgng <- function(dt,meanlog,sdlog,t0,st0=0)
# Same as n1PDF.lnr except NAs dealt with by integration
{

  if ( st0 < 1e-3 ) # smaller values can cause numerical integration problems
    return(n1PDFfixedt0.lnrgng(meanlog=meanlog,sdlog=sdlog,dt=matrix(
      pmax(rep(dt,each=length(meanlog))-t0,0),nrow=length(meanlog)))
    ) else
  {

    integrate.f <- function(dt,meanlog,sdlog,t0,st0)
      n1PDFfixedt0.lnrgng(meanlog=meanlog,sdlog=sdlog,dt=matrix(pmax(
        rep(dt,each=length(meanlog))-t0,0),nrow=length(meanlog)))/st0

    outs <- numeric(length(dt))
    for (i in 1:length(outs)) {
      tmp <- try(integrate(f=integrate.f,
        lower=dt[i]-st0[1],upper=dt[i],
        meanlog=meanlog,sdlog=sdlog,t0=t0,st0=st0[1])$value,silent=TRUE)
      if (is.numeric(tmp))
        outs[i] <- tmp else
        outs[i] <- 0
    }
    return(outs)
  }
}


# # Check
# n=1e5
# meanlog=c(.5,.75); sdlog=c(1,1); t0=.2
# sim <- rlnrgng(n=n,meanlog,sdlog,t0,st0=0)
# dns <- plot.cell.density(sim,xlim=c(0,7),save.density=TRUE)
# d <- n1PDF.lnrgng(dns$'2'$x,meanlog[c(2,1)],sdlog[c(2,1)],t0,st0=0)
# lines(dns$'2'$x,d,col="red")
#
# # p(Stop check)
# mean(is.na(sim$RT))
# n1PDFfixedt0.lnrgng(dt=matrix(rep(NA,2),ncol=1),meanlog,sdlog)
}

####################### LBA go-nogo ----
{

rlba.normgng <- function (n, A, b, t0, mean_v, sd_v, st0 = 0, posdrift = TRUE)
# Race among length(mean_v) accumulators, first of which is a no-go
# accumulator. For trials with winning first accumulator RT set to NA.
{
   out <- rlba.norm(n,A,b,t0,mean_v,sd_v,st0,posdrift)
   out[out$response==1,"rt"] <- NA
   out$response <- factor(as.numeric(out$response))
   data.frame(out)
}



n1PDFfixedt0.normgng=function(dt,A,b,mean_v,sd_v,st0=0,posdrift=TRUE)
# Same as n1PDFfixedt0 except dt=NA done by integration
{

  stopfn <- function(t,A,b,mean_v,sd_v,st0=0,posdrift=TRUE)
  {
    n1PDFfixedt0.norm(
      matrix(rep(t,each=length(mean_v)),nrow=length(mean_v)),
      A=A,b=b,mean_v=mean_v,sd_v=sd_v,posdrift=posdrift
    )
  }

  n.trials <- dim(dt)[2]
  out <- numeric(n.trials)
  is.stop <- is.na(dt[1,])
  if (any(!is.stop)) dt[1,!is.stop] <- n1PDFfixedt0.norm(dt[,!is.stop,drop=F],
    A=A,b=b,mean_v=mean_v,sd_v=sd_v,posdrift=posdrift)
  # tmp <- try(integrate(f=stopfn,lower=0,upper=Inf,
  #     A=A,b=b,mean_v=mean_v,sd_v=sd_v,posdrift=posdrift)$value,silent=TRUE)
  # if (!is.numeric(tmp)) tmp <- 0
  tmp <- my.integrate(f=stopfn,lower=0,
      A=A,b=b,mean_v=mean_v,sd_v=sd_v,posdrift=posdrift)
  dt[1,is.na(dt[1,])] <- tmp
  dt[1,]
}

n1PDF.normgng <- function(dt,A,b,t0,mean_v,sd_v,st0=0,posdrift=TRUE)
# Same as n1PDF except NAs dealt with by integration
{

  if ( st0 < 1e-3 ) # smaller values can cause numerical integration problems
    return(n1PDFfixedt0.normgng(
      A=A,b=b,mean_v=mean_v,sd_v=sd_v,posdrift=posdrift,
      dt=matrix(pmax(rep(dt,each=length(mean_v))-t0,0),nrow=length(mean_v))
    )) else
  {

    integrate.f <- function(A,b,t0,mean_v,sd_v,st0,posdrift)
      n1PDFfixedt0.normgng(
        A=A,b=b,mean_v=mean_v,sd_v=sd_v,posdrift=posdrift,
        dt=matrix(pmax(rep(dt,each=length(mean_v))-t0,0),
                  nrow=length(mean_v)))/st0

    outs <- numeric(length(dt))
    for (i in 1:length(outs)) {
      tmp <- try(integrate(f=integrate.f,
        lower=dt[i]-st0[1],upper=dt[i],A=A,b=b,t0=t0,mean_v=mean_v,sd_v=sd_v,
        st0=st0[1],posdrift=posdrift)$value,silent=TRUE)
      if (is.numeric(tmp))
        outs[i] <- tmp else
        outs[i] <- 0
    }
    return(outs)
  }
}


# # Check
# n=1e5
# A=c(1,1);b=c(2,2);t0=.2;mean_v=c(1,0);sd_v=c(1,1);st0=0;posdrift=TRUE
# sim <- rlba.normgng(n,A,b,t0,mean_v,sd_v,st0,posdrift)
# names(sim) <- c("RT","R")
# dns <- plot.cell.density(sim,xlim=c(0,7),save.density=TRUE)
# d <- n1PDF.normgng(dt=dns$'2'$x,A=A[2:1],b=b[2:1],mean_v=mean_v[2:1],
#   sd_v=sd_v[2:1],t0=t0,st0=st0,posdrift=posdrift)
# # d <- n1PDF.normgng(dt=dns$'2'$x,A=A,b=b,mean_v=mean_v,
# #   sd_v=sd_v,t0=t0,st0=st0,posdrift=posdrift)
#
# lines(dns$'2'$x,d,col="red")
#
# # p(Stop check)
# mean(is.na(sim$RT))
# n1PDFfixedt0.normgng(dt=matrix(rep(NA,2),ncol=1),A=A,b=b,
#                      mean_v=mean_v,sd_v=sd_v)

}


####################### LNR stop signal ----
{
# NB: no st0

rlnrss <- function (n, meanlog, sdlog, t0, t0sg, tf=0, gf=0, ts = 0,
                    SSD=Inf, TRIALS = NA, staircase=NA)
# Race among length(meanlog) accumulators (if meanlog a vector),
# or dim(meanlog)[1] (if a matrix), first of which is a stop accumulator.
# Acts the same as rlnr except NA returned for RT when winner = 1.
# Optional SSD argument can be used to adjust start time for first
# accumulator. SSD can be a scalar or vector length n; output has an SSD column
# For trials with winning first accumulator RT and R set to NA.
# tf = trigger failure probability, gf = go failure probability
# If any !is.na in staircase runs a staircase
# t0 is a scalar with the standard interpritaiton for GO accumulators
# sgt0 = stop encoding - go encoding time + t0 (non-negative scalar)
# ts = slope of slowing (speeding if negative) over TRIALS, meanlog - ts*TRIALS
# This has a linear effect on mean and sd

{
  if ( length(SSD)==1 ) SSD <- rep(SSD,n)
  if ( any(is.na(SSD)) || length(SSD) != n )
    stop("SSD cannot have NAs and must be a scalar or same length as n")

  n_acc <- ifelse(is.null(dim(meanlog)),length(meanlog),dim(meanlog)[1])

  # t0sg -> sg for stop accumulator, relative to 0 for go accumulators
  t0S <- matrix(rep(c(t0sg-t0[1],rep(0,n_acc-1)),length.out=n*n_acc),nrow=n_acc)

  if (!is.matrix(meanlog)) meanlog <- matrix(rep(meanlog,n),nrow=n_acc)
  if (!is.matrix(sdlog)) sdlog <- matrix(rep(sdlog,n),nrow=n_acc)

  if ( !any(is.na(TRIALS)) ) {
    if (length(TRIALS)!=n)
      stop("TRIALS must have length n")
    meanlog[-1,] <- meanlog[-1,] + rep(ts*TRIALS,each=n_acc-1)
  }

  if ( gf > 0 ) # Setup for GO failure
    is.gf <- as.logical(rbinom(length(SSD),1,gf)) else
    is.gf <- logical(length(SSD))

  if ( all(!is.finite(SSD)) ) {              # ALL GO
    out <- rlnr(n,meanlog[-1,,drop=FALSE],sdlog[-1,,drop=FALSE],t0)
    out$R <- out$R+1
  } else {                                   # SOME STOP
    if ( any(is.na(staircase)) ) {           # STOP fixed SSD
      # add SSD to stop accumulator
      t0S[1,] <- t0S[1,] + SSD
      out <- rlnr(n,meanlog,sdlog,t0S)
      if ( tf>0 ) {
        is.tf <- logical(length(SSD))
        is.tf[is.finite(SSD)][as.logical(rbinom(sum(is.finite(SSD)),1,tf))] <- TRUE
        if ( any(is.tf) ) {
          out[is.tf,] <- rlnr(sum(is.tf),meanlog[-1,is.tf,drop=FALSE],
            sdlog[-1,is.tf,drop=FALSE],t0=0)
          out[is.tf,"R"] <- out[is.tf,"R"]+1
        }
      }
    } else {                                 # STOP, staircase
      if ( !is.numeric(staircase) | length(staircase)!=1 )
        stop("Staircase must be a numeric vector of length 1 specifying the step.")
      SSDi <- SSD[is.finite(SSD)][1] # begining SSD
      dt <- matrix(rlnorm(n = n*n_acc, meanlog = meanlog, sdlog = sdlog),
                 nrow = n_acc) + t0S
      # Setup
      winner <- numeric(n)
      for ( i in c(1:n) ) {
        if ( !is.finite(SSD[i]) )   # not staircase
          dt[1,i] <- dt[1,i] + SSD[i] else
          dt[1,i] <- dt[1,i] + SSDi # staircase
        if ( runif(1)<tf ) # Trigger failure
          winner[i] <- which.min(dt[2:n_acc,i])+1 else
          winner[i] <- which.min(dt[,i])
        if (is.gf[i]) winner[i] <- 1
        if ( is.finite(SSD[i]) ) { # update staircase
          SSD[i] <- SSDi
          if ( winner[i]==1 )
            SSDi <- SSDi + staircase else
            SSDi <- SSDi - staircase
            if (SSDi<1e-10) SSDi <- 0
        }
      }
      out <- data.frame(RT=dt[cbind(winner,1:n)],R=winner)
    }
    out$RT <- out$RT + t0 # Add t0 for go responses
  }
  out[out$R==1,"RT"] <- NA
  if (gf > 0) {
    out$RT[is.gf] <- NA
    out$R[is.gf] <- 1
  }
  if ( any(is.na(TRIALS)) ) cbind.data.frame(out,SSD=SSD) else
                            cbind.data.frame(out,SSD=SSD,TRIALS=TRIALS)
}




n1PDF.lnrss <- function(rt,meanlog,sdlog,t0,t0sg,tf=0,gf=0,ts=0,
                        SSD=Inf,TRIALS=NA,Si)
# Same as n1PDF.lnr except SSD is either a scalar or vector of length(rt)
# stop accumulator must have name "NR". SSD is subtracted stop accumulator time
# and dt=NA done by integration.
#
# tf= probabiliy of trigger failure, where
# L = trigger fail & respond + trigger and respond + trigger and no-response
#   = tf*L(N-1)+(1-tf)[L(N)+p(S)],
# L(N-1) = choice race likelihood (no stop accumulator),
# L(N) = full N unit race likelihood given they did respond,
# p(S) probability of stop winning
#
# gf = probabiliy of go failure.
# On go trials:   L = go fail (so no response) + go and L above
# L = gf + (1-gf)*[tf*L(N-1)+(1-tf)[L(N)+p(S)]]  or similarly
# L =    [ p(non-response) ]    +           [ p(response) ]
#   = [ gf + (1-gf)(1-tf)p(S) ] + [ (1-gf){(tf*Ln(n-1) + (1-tf)*L(N))} ]
#
# NB:rt is NOT decision time, but rather full RT as t0 has to be passed
#    in order to include properly in cases where RT is NA (i.e., sucessful stop)

{

  stopfn <- function(t,meanlogj,sdlogj,t0,t0sg,SSD,Si)
  {
    # d = s-g = tD-t0(GO), then add SSDstart time to get
    # start time go - start time stop, i.e., finishing time advantage for GO
    t0S <- rep(0,length(meanlogj)) # Set relative to stop accumulator finish time
    t0S[-Si] <- t0S[-Si]+t0sg-t0+SSD
    # min(t0) keeps GO time positive
    dt <- matrix(rep(t,each=length(meanlogj)),nrow=length(meanlogj))+t0S-min(t0S)
    i <- c(Si,c(1:length(meanlogj))[-Si])
    n1PDFfixedt0.lnr(dt[i,,drop=FALSE],meanlogj[i],sdlogj[i])
  }

  # NOTE: t0 is not subtracted when making dt but passed to handle RT=NA case

  if ( length(SSD)==1 ) SSD <- rep(SSD,length(rt))
  if (length(SSD) != length(rt))
    stop("SSD must be a scalar or same length as rt")
  n_acc <- ifelse(is.null(dim(meanlog)),length(meanlog),dim(meanlog)[1])

  rt <- matrix(rep(rt,each=n_acc),nrow=n_acc)
  is.stop <- is.na(rt[1,])

  if (!is.matrix(meanlog)) meanlog <- matrix(rep(meanlog,dim(rt)[2]),nrow=n_acc)
  if (!is.matrix(sdlog)) sdlog <- matrix(rep(sdlog,dim(rt)[2]),nrow=n_acc)

  if ( any(is.na(TRIALS)) | ts == 0 ) {
    p <- SSD[is.stop]
    pj <- c(1:sum(is.stop))[!duplicated(p)] # index of unique SSD
  } else {
    meanlog[-Si,] <- meanlog[-Si,] + ts*TRIALS
    p <- apply(
        rbind(meanlog[,is.stop,drop=FALSE],SSD[is.stop]),
      2,paste,collapse="")
    pj <- c(1:sum(is.stop))[!duplicated(p)] # index of unique p and SSD
  }

  if ( any(!is.stop) )
  {
    rt[Si,!is.stop] <- rt[Si,!is.stop] - t0sg - SSD[!is.stop]
    rt[-Si,!is.stop] <- rt[-Si,!is.stop]-t0
    if ( tf > 0 )
    {
      rt[1,!is.stop] <- (1-gf)*(
        tf*n1PDFfixedt0.lnr(rt[-Si,!is.stop,drop=FALSE],
          meanlog[-Si,!is.stop,drop=FALSE],sdlog[-Si,!is.stop,drop=FALSE]) +
        (1-tf)*n1PDFfixedt0.lnr(rt[,!is.stop,drop=FALSE],
          meanlog[,!is.stop,drop=FALSE],sdlog[,!is.stop,drop=FALSE])
      )
    } else
      rt[1,!is.stop] <- (1-gf)*n1PDFfixedt0.lnr(rt[,!is.stop,drop=FALSE],
        meanlog[,!is.stop,drop=FALSE],sdlog[,!is.stop,drop=FALSE])
  }

  if ( any(is.stop) ) for (j in pj) {
# fix me
print (j)
    #     tmp <- ifelse(!is.finite(SSD[j]),0,
#       try(integrate(f=stopfn,lower=0,upper=Inf,meanlogj=meanlog[,j],
#         sdlogj=sdlog[,j],t0=t0,t0sg=t0sg,SSD=SSD[j],Si=Si)$value,silent=TRUE))
#     if (!is.numeric(tmp)) tmp <- 0
    if (!is.finite(SSD[j])) tmp <- 0 else
      tmp <- my.integrate(f=stopfn,lower=0,meanlogj=meanlog[,j],
        sdlogj=sdlog[,j],t0=t0,t0sg=t0sg,SSD=SSD[j],Si=Si)
     rt[1,is.stop][p %in% p[j]] <- gf +(1-gf)*(1-tf)*tmp
  }
  rt[1,]
}


# # VERY EXTENSIVE TESTING WITH Two different SSDs
#
# # ########### TWO ACCUMULATOR CASE
#
# n=1e5
# meanlog=c(.75,.75); sdlog=c(.5,1)
# SSD = rep(c(1,10)/10,each=n/2)
# do.trials=FALSE
# do.trials = TRUE # requires very differnet plotting check, can be SLOW!
#
# #### RUN ONE OF THE FOLLOWING THREE LINES, all assume .2s go Ter
# t0=.2; t0sg= 0  # minimum possible value of t0sg, stop encoding .2 less than go
# t0=.2; t0sg=.2 # equal stop and go enconding times
# t0=.2; t0sg=.4 # stop .2 slower than go
#
# ### RUN ONE OF THE FOLLOWING FOUR LINES
# # Without trigger failure or go failure
# tf=0; gf=0
# # With trigger failure, no go failure
# tf=.1;gf=0
# # Without trigger failure, with go failure
# tf=0; gf=.1
# # With trigger failure and go failure
# tf=.1;gf=.1
#
# if (do.trials) {
#   ts=.5; TRIALS=log10(seq(1,10,length.out=n)) # 1..10 natural so 0-1 on log
#   TRIALS <- as.vector(t(matrix(TRIALS,nrow=2))) # interleave SSDs
#   # Plot slowing in GO (usually nice and linear, up to smooting overfitting)
#   sim.go <- rlnrss(n=n,meanlog,sdlog,t0,t0sg,tf=tf,gf=gf,TRIALS=TRIALS,ts=ts)
#   is.in <- !is.na(sim.go$RT) # in case go failure
#   plot(smooth.spline(TRIALS[is.in],sim.go$RT[is.in]),ylab="Smooth",xlab="TRIALS",type="l")
# } else {TRIALS=NA;ts=0}
#
# # Simulate stop trials
# sim <- rlnrss(n=n,meanlog,sdlog,t0,t0sg,SSD=SSD,tf=tf,gf=gf,TRIALS=TRIALS,ts=ts)
#
# # Plot densities
# par(mfrow=c(1,2))
# dns1 <- plot.cell.density(sim[sim$SSD==.1,],xlim=c(0,7),save.density=TRUE,main="SSD=.1")
# dns2 <- plot.cell.density(sim[sim$SSD!=.1,],xlim=c(0,7),save.density=TRUE,main="SSD=1")
# x1c <- dns1$'2'$x; x2c <- dns2$'2'$x
#
# # Signal respond RT
# dat <- sim; dat <- dat[!is.na(dat$RT),]; dat$R <- factor(as.character(dat$R))
# round(tapply(dat$RT,dat[,c("R","SSD")],mean),2)
#
# if (do.trials) {
#   tmp <- n1PDF.lnrss(sim$RT[!is.na(sim$RT)],meanlog[2:1],sdlog[2:1],t0,t0sg,
#     ts=ts,TRIALS=TRIALS[!is.na(sim$RT)],SSD=SSD[!is.na(sim$RT)],Si=2,tf=tf,gf=gf)
#   par(mfrow=c(1,2))
#   # red=black?
#   plot(x1c,dns1$'2'$y,pch=".",main="SSD=.1",ylab="Density",xlab="RT")
#   lines(smooth.spline(sim$RT[!is.na(sim$RT) & SSD==.1],
#    tmp[c(SSD==.1)[!is.na(sim$RT)]]),col="red")
#   # red=black?
#   plot(x2c,dns2$'2'$y,pch=".",main="SSD=1",ylab="Density",xlab="RT")
#   lines(smooth.spline(sim$RT[!is.na(sim$RT) & SSD==1],
#     tmp[c(SSD==1)[!is.na(sim$RT)]]),col="red")
#   print(tapply(is.na(sim$RT),sim$SSD,mean)) # empirical
#   tmp <- n1PDF.lnrss(rep(NA,n),meanlog,sdlog,t0,t0sg,SSD=SSD,Si=1,tf=tf,gf=gf,ts=ts,TRIALS=TRIALS)
#   print(mean(tmp[SSD==.1]))
#   print(mean(tmp[SSD==1]))
#   plot(TRIALS,tmp,pch=".",xlab="TRIALS",ylab="p(NA)",ylim=c(0,1))
#   lines(smooth.spline(TRIALS[SSD==.1],as.numeric(is.na(sim$RT)[SSD==.1])),col="red")
#   lines(smooth.spline(TRIALS[SSD==1],as.numeric(is.na(sim$RT)[SSD==1])),col="red")
# } else {
#   # Save simulated densities
#   r1 <- c(2,1)
#   d.r1 <- n1PDF.lnrss(rt=c(x1c,x2c),meanlog[r1],sdlog[r1],t0,t0sg,
#     SSD=c(rep(.1,length(x1c)),rep(1,length(x2c))),Si=2,tf=tf,gf=gf)
#   # Plot simulated (black) and theoretical (red) densities
#   par(mfrow=c(1,2))
#   # red=black?
#   plot(x1c,dns1$'2'$y,type="l",main="SSD=.1",ylab="Density",xlab="RT",
#      ylim=c(0,max(dns1$'2'$y)))
#   lines(x1c,d.r1[1:length(x1c)],col="red")
#   # red=black?
#   plot(x2c,dns2$'2'$y,type="l",main="SSD=1",ylab="Density",xlab="RT",
#      ylim=c(0,max(dns2$'2'$y)))
#   lines(x2c,d.r1[(length(x2c)+1):(2*length(x2c))],col="red")
#
#   # p(Stop check)
#   print(tapply(is.na(sim$RT),sim$SSD,mean)) # empirical
#   print(n1PDF.lnrss(NA,meanlog,sdlog,t0,t0sg,SSD=.1,Si=1,tf=tf,gf=gf))
#   print(n1PDF.lnrss(NA,meanlog,sdlog,t0,t0sg,SSD=1,Si=1,tf=tf,gf=gf))
# }
#
#
#
# ########### THREE ACCUMULATOR CASE
#
# n=1e5
# meanlog=c(.75,.75,1); sdlog=c(.5,1,1)
# SSD = rep(c(1,10)/10,each=n/2)
#
# do.trials=FALSE
# do.trials = TRUE # requires very differnet plotting check, can be SLOW!
#
# #### RUN ONE OF THE FOLLOWING THREE LINES, all assume .2s go Ter
# t0=.2; t0sg= 0  # minimum possible value of t0sg, stop encoding .2 less than go
# t0=.2; t0sg=.2 # equal stop and go enconding times
# t0=.2; t0sg=.4 # stop .2 slower than go
#
# ### RUN ONE OF THE FOLLOWING FOUR LINES
# # Without trigger failure or go failure
# tf=0; gf=0
# # With trigger failure, no go failure
# tf=.1;gf=0
# # Without trigger failure, with go failure
# tf=0; gf=.1
# # With trigger failure and go failure
# tf=.1;gf=.1
#
# if (do.trials) {
#   ts=.5; TRIALS=log10(seq(1,10,length.out=n)) # 1..10 natural so 0-1 on log
#   TRIALS <- as.vector(t(matrix(TRIALS,nrow=2))) # interleave SSDs
#   # Plot slowing in GO (usually nice and linear, up to smooting overfitting)
#   sim.go <- rlnrss(n=n,meanlog,sdlog,t0,t0sg,tf=tf,gf=gf,TRIALS=TRIALS,ts=ts)
#   is.in <- !is.na(sim.go$RT) # in case go failure
#   plot(smooth.spline(TRIALS[is.in],sim.go$RT[is.in]),ylab="Smooth",xlab="TRIALS",type="l")
# } else {TRIALS=NA;ts=0}
#
# # Simulate stop trials
# sim <- rlnrss(n=n,meanlog,sdlog,t0,t0sg,SSD=SSD,tf=tf,gf=gf,TRIALS=TRIALS,ts=ts)
#
# par(mfrow=c(1,2))
# dns1 <- plot.cell.density(sim[sim$SSD==.1,],xlim=c(0,7),save.density=TRUE,main="SSD=.1")
# dns2 <- plot.cell.density(sim[sim$SSD!=.1,],xlim=c(0,7),save.density=TRUE,main="SSD=1")
# x1c <- dns1$'2'$x; x2c <- dns2$'2'$x
# x1e <- dns1$'3'$x; x2e <- dns2$'3'$x
#
# # Signal respond RT
# dat <- sim; dat <- dat[!is.na(dat$RT),]; dat$R <- factor(as.character(dat$R))
# round(tapply(dat$RT,dat[,c("R","SSD")],mean),2)
#
# if (do.trials) {
#   r1 <- c(2,1,3)
#   is.in1 <- !is.na(sim$RT) & sim$R==2
#   d.r1 <- n1PDF.lnrss(sim$RT[is.in1],meanlog[r1],ts=ts,TRIALS=TRIALS[is.in1],
#                    sdlog[r1],t0,t0sg,SSD=SSD[is.in1],Si=2,tf=tf,gf=gf)
#   r2 <- c(3,1,2)
#   is.in2 <- !is.na(sim$RT) & sim$R==3
#   d.r2 <- n1PDF.lnrss(sim$RT[is.in2],meanlog[r2],ts=ts,TRIALS=TRIALS[is.in2],
#                    sdlog[r2],t0,t0sg,SSD=SSD[is.in2],Si=2,tf=tf,gf=gf)
#   par(mfrow=c(1,3))
#   # red=black?
#   plot(x1c,dns1$'2'$y,pch=".",main="SSD=.1",ylab="Density",xlab="RT",type="l")
#   lines(x1e,dns1$'3'$y,lty=2)
#   lines(smooth.spline(sim$RT[is.in1 & sim$SSD==.1],
#                       d.r1[c(sim$SSD==.1)[is.in1]]),col="red")
#   lines(smooth.spline(sim$RT[is.in2 & sim$SSD==.1],d.r2[c(sim$SSD==.1)[is.in2]]),
#         lty=2,col="red")
#   # red=black?
#   plot(x2c,dns2$'2'$y,pch=".",main="SSD=1",ylab="Density",xlab="RT",type="l")
#   lines(x2e,dns2$'3'$y,lty=2)
#   lines(smooth.spline(sim$RT[is.in1 & sim$SSD==1],
#                       d.r1[c(sim$SSD==1)[is.in1]]),col="red")
#   lines(smooth.spline(sim$RT[is.in2 & sim$SSD==1],
#                       d.r2[c(sim$SSD==1)[is.in2]]),col="red",lty=2)
#
#   print(tapply(is.na(sim$RT),sim$SSD,mean)) # empirical
#   tmp <- n1PDF.lnrss(rep(NA,n),meanlog,sdlog,t0,t0sg,SSD=SSD,Si=1,tf=tf,gf=gf,ts=ts,TRIALS=TRIALS)
#   print(mean(tmp[SSD==.1]))
#   print(mean(tmp[SSD==1]))
#   plot(TRIALS,tmp,pch=".",xlab="TRIALS",ylab="p(NA)",ylim=c(0,1))
#   lines(smooth.spline(TRIALS[SSD==.1],as.numeric(is.na(sim$RT)[SSD==.1])),col="red")
#   lines(smooth.spline(TRIALS[SSD==1],as.numeric(is.na(sim$RT)[SSD==1])),col="red")
# } else {
#   # Save simulated densities
#   r1 <- c(2,1,3)
#   d.r1 <- n1PDF.lnrss(rt=c(x1c,x2c),meanlog[r1],sdlog[r1],t0,t0sg,
#     SSD=c(rep(.1,length(x1c)),rep(1,length(x2c))),Si=2,tf=tf,gf=gf)
#   r2 <- c(3,1,2)
#   d.r2 <- n1PDF.lnrss(rt=c(x1e,x2e),meanlog[r2],sdlog[r2],t0,t0sg,
#     SSD=c(rep(.1,length(x1e)),rep(1,length(x2e))),Si=2,tf=tf,gf=gf)
#   # Plot simulated (black) and theoretical (red) densities
#   par(mfrow=c(1,2))
#   # red=black?
#   plot(x1c,dns1$'2'$y,type="l",main="SSD=.1",ylab="Density",xlab="RT",
#      ylim=c(0,max(dns1$'2'$y)))
#   lines(x1c,d.r1[1:length(x1c)],col="red")
#   lines(x1e,dns1$'3'$y,lty=2)
#   lines(x1e,d.r2[1:length(x1e)],col="red",lty=2)
#   # red=black?
#   plot(x2c,dns2$'2'$y,type="l",main="SSD=1",ylab="Density",xlab="RT",
#      ylim=c(0,max(dns2$'2'$y)))
#   lines(x2c,d.r1[(length(x2c)+1):(2*length(x2c))],col="red")
#   lines(x2e,dns2$'3'$y,lty=2)
#   lines(x2e,d.r2[(length(x2e)+1):(2*length(x2e))],col="red",lty=2)
#
#   # p(Stop check)
#   print(tapply(is.na(sim$RT),sim$SSD,mean)) # empirical
#   print(n1PDF.lnrss(NA,meanlog,sdlog,t0,t0sg,SSD=.1,Si=1,tf=tf,gf=gf))
#   print(n1PDF.lnrss(NA,meanlog,sdlog,t0,t0sg,SSD=1,Si=1,tf=tf,gf=gf))
# }

}


####################### ExGaussian stop signal----
{
# WHY TRIAL SCALING IS BAD IN EXG, mean not proportional to sd
# TRIALS=c(0:10)/10
# mu=.5; sigma=.025; tau=.05; ts=.25
# sigma=sigma+ts*TRIALS
# tau=tau + ts*TRIALS
# mn= mu + tau
# sd=sqrt( (sigma)^2 + (tau)^2 )
# ratio=mn/sd
# round(rbind(sigma,tau,mn,sd,ratio),3)
# ratio[1]/ratio[10]


# Modified from gamlss.dist to make cdf in nu > 0.05 * sigma case robust,
# and robust to -Inf and Inf inputs, returns NA for bad sigma or tau
pexGAUS <- function (q, mu = 5, sigma = 1, nu = 1, lower.tail = TRUE, log.p = FALSE)
{
  if (any(sigma <= 0)) return(rep(NA,length(q)))
  if (any(nu <= 0)) return(rep(NA,length(q)))
  ly <- length(q)
  sigma <- rep(sigma, length = ly)
  mu <- rep(mu, length = ly)
  nu <- rep(nu, length = ly)
  index <- seq(along = q)
  z <- q - mu - ((sigma^2)/nu)
  cdf <- ifelse(is.finite(q),
      ifelse(nu > 0.05 * sigma,
        pnorm((q - mu)/sigma) -
          exp(log(pnorm(z/sigma)) + ((mu + (sigma^2/nu))^2 - (mu^2) -
            2 * q * ((sigma^2)/nu))/(2 * sigma^2)),
        pnorm(q, mean = mu, sd = sigma)),
      ifelse(q<0,0,1)
    )
  if (lower.tail == TRUE)
      cdf <- cdf
    else cdf <- 1 - cdf
  if (log.p == FALSE)
      cdf <- cdf
    else cdf <- log(cdf)
  cdf
}

# gamlss.dist function, but returns NA for bad sigma or tau
dexGAUS <- function (x, mu = 5, sigma = 1, nu = 1, log = FALSE)
{
    if (any(sigma <= 0)) return(rep(NA,length(x)))
    if (any(nu <= 0)) return(rep(NA,length(x)))
    ly <- length(x)
    sigma <- rep(sigma, length = ly)
    mu <- rep(mu, length = ly)
    nu <- rep(nu, length = ly)
    z <- x - mu - ((sigma^2)/nu)
    logfy <- ifelse(nu > 0.05 * sigma, -log(nu) - (z + (sigma^2/(2 *
        nu)))/nu + log(pnorm(z/sigma)), dnorm(x, mean = mu, sd = sigma,
        log = TRUE))
    if (log == FALSE)
        fy <- exp(logfy)
    else fy <- logfy
    fy
}

rexg <- function (mu,sigma,tau)
# NOTE: gamlss.dist rexGAUS is very slow so didnt use.
# Ex-gaussian race, mu is a matrix with 1 row per accumulator, columns
# are trials, ouput is a two column (RT, R) data frame, on row for each trial.
# Sigma and tau can be vectors with one value, or one value for each accumulator
# or a matrix in the same form as mu. NB: CAN PRODUCE NEGATIVE RTS!
{
    dt <- matrix(rnorm(n=length(mu), mean=mu,sd=sigma) +
                 rexp(n=length(mu),rate=1/tau),nrow = dim(mu)[1])
    winner <- apply(dt,2,which.min)
    data.frame(RT=dt[cbind(winner,1:dim(mu)[2])],R=winner)
}

n1PDF.exg <- function(dt,mu,sigma,tau)
# Generates defective PDF for responses on node= 1
# dt (decison time) is a matrix with length(mu) rows, one row for
# each accumulator to allow for different start times
{

  dt[1,] <- dexGAUS(dt[1,],mu[1],sigma[1],tau[1])
  if (length(mu)>1) for (i in 2:length(mu))
    dt[1,] <- dt[1,]*pexGAUS(q=dt[i,],mu=mu[i],sigma=sigma[i],nu=tau[i],lower.tail=FALSE)
  dt[1,]
}


rexgss <- function (n, mu, sigma, tau, tf=0, gf=0,
                    SSD=Inf, staircase=NA)
# Race among n accumulators, first of which is a stop accumulator.
# NA returned for RT when winner = 1. Optional SSD arguement can be used to
# adjust mu for first accumulator to mu+SSD. SSD can be a scalar or vector
# length n. For trials with winning first accumulator RT and R set to NA.
# Adds SSD column to output.
# tf = trigger failure probability, gf = go failure probability
{
  if ( length(SSD)==1 ) SSD <- rep(SSD,n)
  if ( any(is.na(SSD)) || length(SSD) != n )
    stop("SSD cannot have NAs and must be a scalar or same length as n")
  if ( !is.matrix(mu) )
    mu <- matrix(rep(mu,times=n),nrow=length(mu))
  if (gf > 0) # Setup for GO failure
    is.gf <- as.logical(rbinom(length(SSD),1,gf)) else
    is.gf <- logical(length(SSD))

  if ( all(!is.finite(SSD)) ) {              # ALL GO
    mu[1,] <- mu[1,] + SSD
    out <- rexg(mu,sigma,tau)
  } else {                                   # SOME STOP
    if ( any(is.na(staircase)) ) {           # STOP fixed SSD
      mu[1,] <- mu[1,] + SSD
      out <- rexg(mu,sigma,tau)
      if ( tf>0 ) {
        is.tf <- logical(length(SSD))
        is.tf[is.finite(SSD)][as.logical(rbinom(sum(is.finite(SSD)),1,tf))] <- TRUE
        if ( any(is.tf) ) {
          out[is.tf,] <-
            rexg(mu[-1,is.tf,drop=FALSE],sigma[-1],tau[-1])
          out[is.tf,"R"] <- out[is.tf,"R"]+1
        }
      }
    } else {                                 # STOP, staircase
      if ( !is.numeric(staircase) | length(staircase)!=1 )
        stop("Staircase must be a numeric vector of length 1 specifying the step.")
      n_acc <- dim(mu)[1]
      dt <- matrix(rnorm(n=length(mu), mean=mu,sd=sigma) +
                 rexp(n=length(mu),rate=1/tau),nrow = dim(mu)[1])
      winner <- numeric(n)
      for (i in c(1:n)) {
        dt[1,i] <- dt[1,i] + SSD[i]
        if ( runif(1)<tf ) # Trigger failure
          winner[i] <- which.min(dt[2:n_acc,i])+1 else
          winner[i] <- which.min(dt[,i])
        if (is.gf[i]) winner[i] <- 1
        if (i!=n) if (winner[i]==1)
          SSD[i+1] <- SSD[i] + staircase else
          SSD[i+1] <- pmax(SSD[i] - staircase,0)
      }
      out <- data.frame(RT=dt[cbind(winner,1:n)],R=winner)
    }
  }
  out[out$R==1,"RT"] <- NA
  if (gf > 0) {
    out$RT[is.gf] <- NA
    out$R[is.gf] <- 1
  }
  cbind.data.frame(out,SSD=SSD)
}


n1PDF.exgss=function(dt,mu,sigma,tau,tf=0,gf=0,
                     SSD,Si)
# SSD is either a scalar or vector of length(dt)
# stop accumulator must have name "NR"
# SSD is subtracted from dt[Si,]
# (i.e., stop accumulator RT) and dt=NA done by integration.
# ts = slope of slowing (speeding if negative) over TRIALS, linear affect on
# mean and sd:  sigma + ts*TRIALS, tau + ts*trials (truncated to be positive)
#
# tf= probabiliy of trigger failure, where
# L = trigger fail & respond + trigger and respond + trigger and no-response
#   = tf*L(N-1)+(1-tf)[L(N)+p(S)],
# L(N-1) = choice race likelihood (no stop accumulator),
# L(N) = full N unit race likelihood given they did respond,
# p(S) probability of stop winning
#
# gf = probabiliy of go failure.
# On go trials:   L = go fail (so no response) + go and L above
# L = gf + (1-gf)*[tf*L(N-1)+(1-tf)[L(N)+p(S)]]  or similarly
# L =    [ p(non-response) ]    +           [ p(response) ]
#   = [ pf + (1-gf)(1-tf)p(S) ] + [ (1-gf){(tf*Ln(n-1) + (1-tf)*L(N))} ]
{

  stopfn <- function(t,mu,sigma,tau,SSD,Si)
  {
    dt <- matrix(rep(t+SSD,each=length(mu)),nrow=length(mu))
    dt[Si,] <- dt[Si,]-SSD
    i <- c(Si,c(1:length(mu))[-Si])
    n1PDF.exg(dt[i,],mu[i],sigma[i],tau[i])
  }

  dt=matrix(rep(dt,each=length(mu)),nrow=length(mu))

  is.stop <- is.na(dt[1,])
  dt[Si,] <- dt[Si,] - SSD
  if ( any(!is.stop) )
  {
    if ( tf > 0 )
    {
      dt[1,!is.stop] <- (1-gf)*(
        tf*n1PDF.exg(dt[-Si,!is.stop,drop=FALSE],mu[-Si],sigma[-Si],tau[-Si]) +
        (1-tf)*n1PDF.exg(dt[,!is.stop,drop=FALSE],mu,sigma,tau))
    } else
      dt[1,!is.stop] <- (1-gf)*n1PDF.exg(dt[,!is.stop,drop=FALSE],mu,sigma,tau)
  }
  if ( any(is.stop) ) for ( i in unique(SSD[is.stop]) ) {
    # tmp <- ifelse(!is.finite(i),0,
    #   try(integrate(f=stopfn,lower=-Inf,upper=Inf,
    #   mu=mu,sigma=sigma,tau=tau,SSD=i,Si=Si)$value,silent=TRUE))
    # if (!is.numeric(tmp)) tmp <- 0
    tmp <- ifelse(!is.finite(i),0,
      my.integrate(f=stopfn,lower=-Inf,
                   mu=mu,sigma=sigma,tau=tau,SSD=i,Si=Si))
    dt[1,is.stop & (SSD==i)] <- gf + (1-gf)*(1-tf)*tmp
  }
  dt[1,]
}


# ########### TWO ACCUMULATOR CASE
#
# # Check, two different SSDs
# n <- 1e5
# SSD <- rep(c(1,10)/10,each=n/2)
# TRIALS <- NA
#
# # Stop and one go accumulator
# mu=c(.75,1); sigma=c(.25,.25); tau=c(.25,.5)
#
# ### RUN ONE OF THE FOLLOWING FOUR LINES
# # Without trigger failure or go failure
# tf=0; gf=0
# # With trigger failure, no go failure
# tf=.1;gf=0
# # Without trigger failure, with go failure
# tf=0; gf=.1
# # With trigger failure and go failure
# tf=.1;gf=.1
#
# # Simulate data and get simlate stop rate
# sim <- rexgss(n=n,mu,sigma,tau,SSD=SSD,tf=tf,gf=gf,TRIALS=TRIALS,ts=ts)
# # Plot data
# par(mfrow=c(1,2))
# dns1 <- plot.cell.density(sim[sim$SSD==.1,],xlim=c(0,7),save.density=TRUE,main="SSD=.1")
# dns2 <- plot.cell.density(sim[sim$SSD!=.1,],xlim=c(0,7),save.density=TRUE,main="SSD=1")
#
# # Signal respond RT
# dat <- sim; dat <- dat[!is.na(dat$RT),]; dat$R <- factor(as.character(dat$R))
# round(tapply(dat$RT,dat[,c("R","SSD")],mean),2)
#
# # Save simulated densities
# x1 <- dns1$'2'$x; x2 <- dns2$'2'$x
# SSD <- c(rep(.1,length(x1)),rep(1,length(x2)))
# r1 <- c(2,1)
# d.r1 <- n1PDF.exgss(dt=c(x1,x2),mu=mu[r1],sigma=sigma[r1],tau=tau[r1],
#                     SSD=SSD,Si=2,tf=tf,gf=gf)
#
# # Plot simulated (black) and theoretical (red) densities
# par(mfrow=c(1,2))
# # red=black?
# plot(x1,dns1$'2'$y,type="l",main="SSD=.1",ylab="Density",xlab="RT",
#      ylim=c(0,max(dns1$'2'$y)))
# lines(x1,d.r1[1:length(x1)],col="red")
# # red=black?
# plot(x2,dns2$'2'$y,type="l",main="SSD=1",ylab="Density",xlab="RT",
#      ylim=c(0,max(dns2$'2'$y)))
# lines(x2,d.r1[(length(x2)+1):(2*length(x2))],col="red")
#
# # p(Stop check)
# tapply(is.na(sim$RT),sim$SSD,mean) # empirical
# # Theoretical
# n1PDF.exgss(NA,mu,sigma,tau,SSD=.1,Si=1,tf=tf,gf=gf)
# n1PDF.exgss(NA,mu,sigma,tau,SSD=1,Si=1,tf=tf,gf=gf)
#
# ########### THREE ACCUMULATOR CASE
#
# # Check, two different SSDs
# n=1e5
# SSD = rep(c(1,10)/10,each=n/2)
#
# # Stop, first go accumulator correct, second error
# mu=c(.75,1,1.25); sigma=c(.25,.25,.25); tau=c(.25,.5,.5)
#
# ### RUN ONE OF THE FOLLOWING FOUR LINES
# # Without trigger failure or go failure
# tf=0; gf=0
# # With trigger failure, no go failure
# tf=.1;gf=0
# # Without trigger failure, with go failure
# tf=0; gf=.1
# # With trigger failure and go failure
# tf=.1;gf=.1
#
# # Simulate data and get simlate stop rate
# sim <- rexgss(n=n,mu,sigma,tau,SSD=SSD,tf=tf,gf=gf)
# # Plot data
# par(mfrow=c(1,2))
# dns1 <- plot.cell.density(sim[sim$SSD==.1,],xlim=c(0,7),save.density=TRUE,main="SSD=.1")
# dns2 <- plot.cell.density(sim[sim$SSD!=.1,],xlim=c(0,7),save.density=TRUE,main="SSD=1")
#
# # Signal respond RT
# dat <- sim; dat <- dat[!is.na(dat$RT),]; dat$R <- factor(as.character(dat$R))
# round(tapply(dat$RT,dat[,c("R","SSD")],mean),2)
#
# # Save simulated densities
# x1c <- dns1$'2'$x; x2c <- dns2$'2'$x
# x1e <- dns1$'3'$x; x2e <- dns2$'3'$x
#
# SSD <- c(rep(.1,length(x1c)),rep(1,length(x2c)))
# r1 <- c(2,1,3)
# d.r1 <- n1PDF.exgss(dt=c(x1c,x2c),mu=mu[r1],sigma=sigma[r1],tau=tau[r1],
#                     SSD=SSD,Si=2,tf=tf,gf=gf)
# SSD <- c(rep(.1,length(x1e)),rep(1,length(x2e)))
# r2 <- c(3,1,2)
# d.r2 <- n1PDF.exgss(dt=c(x1e,x2e),mu=mu[r2],sigma=sigma[r2],tau=tau[r2],
#                     SSD=SSD,Si=2,tf=tf,gf=gf)
#
# # Plot simulated (black) and theoretical (red) densities
# par(mfrow=c(1,2))
# # red=black?
# plot(x1c,dns1$'2'$y,type="l",main="SSD=.1",ylab="Density",xlab="RT",
#      ylim=c(0,max(dns1$'2'$y)))
# lines(x1c,d.r1[1:length(x1c)],col="red")
# lines(x1e,dns1$'3'$y,lty=2)
# lines(x1e,d.r2[1:length(x1e)],col="red",lty=2)
# # red=black?
# plot(x2c,dns2$'2'$y,type="l",main="SSD=1",ylab="Density",xlab="RT",
#      ylim=c(0,max(dns2$'2'$y)))
# lines(x2c,d.r1[(length(x2c)+1):(2*length(x2c))],col="red")
# lines(x2e,dns2$'3'$y,lty=2)
# lines(x2e,d.r2[(length(x2e)+1):(2*length(x2e))],col="red",lty=2)
#
# # p(Stop check)
# tapply(is.na(sim$RT),sim$SSD,mean) # empirical
# # Theoretical
# n1PDF.exgss(NA,mu,sigma,tau,SSD=.1,Si=1,tf=tf,gf=gf)
# n1PDF.exgss(NA,mu,sigma,tau,SSD=1,Si=1,tf=tf,gf=gf)
}


