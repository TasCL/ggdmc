#######  Standard model ###

#' @importFrom stats runif
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
      drifts <- matrix(rtnorm(n = n * n_v, mean = mean_v, sd = sd_v,
        lower = 0, checks=TRUE), nrow = n_v)
    } else {
      drifts <- matrix(rtnorm(n = n * n_v, mean = mean_v, sd = sd_v),
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

#' @importFrom stats pnorm dnorm
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


n1PDFfixedt0.norm <- function(dt,A,b,mean_v,sd_v,
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


#' @importFrom stats rlnorm
rlnr <- function (n, meanlog, sdlog, t0, st0 = 0)
{
  # Race among n_acc accumulators, mealnlog and sdlog can be n_acc length
  # vectors or n_acc x n matrices. t0 can be
  # a) a scalar, b) a vector of length number of accumulators or
  # c) a matrix with 1 row per accumulator, when start times differ on each trial
  # st0, range of non-decison time variability, must be a scalar, as the same
  # variability is assumed in a common encoding/production stage
  n_acc <- ifelse(is.null(dim(meanlog)),length(meanlog),dim(meanlog)[1])
  dt <- matrix(rlnorm(n = n*n_acc, meanlog = meanlog, sdlog = sdlog),
               nrow = n_acc) + t0
  winner <- apply(dt,2,which.min)
  if (st0[1]==0) data.frame(RT=dt[cbind(winner,1:n)],R=winner) else
    data.frame(RT=dt[cbind(winner,1:n)]+runif(n,0,st0[1]),R=winner)
}


#' @importFrom stats rlnorm
rlnrgng <- function (n, meanlog, sdlog, t0, st0 = 0)
  # Race among length(meanlog) accumulators, first of which is a no-go
  # accumulator. For trials with winning first accumulator RT set to NA.
{
  out <- rlnr(n,meanlog,sdlog,t0,st0)
  out[out$R==1,"RT"] <- NA
  out$R <- factor(as.numeric(out$R))
  data.frame(out)
}



rlba.normgng <- function (n, A, b, t0, mean_v, sd_v, st0 = 0, posdrift = TRUE)
{
  # Race among length(mean_v) accumulators, first of which is a no-go
  # accumulator. For trials with winning first accumulator RT set to NA.
  out <- rlba.norm(n,A,b,t0,mean_v,sd_v,st0,posdrift)
  out[out$response==1,"rt"] <- NA
  out$response <- factor(as.numeric(out$response))
  data.frame(out)
}


#' @importFrom stats rlnorm rbinom runif
rlnrss <- function (n, meanlog, sdlog, t0, t0sg, tf=0, gf=0, ts = 0,
                    SSD=Inf, TRIALS = NA, staircase=NA)
{
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


#' @importFrom stats rexp rnorm
rexg <- function (mu,sigma,tau)
{
  # NOTE: gamlss.dist rexGAUS is very slow so didnt use.
  # Ex-gaussian race, mu is a matrix with 1 row per accumulator, columns
  # are trials, ouput is a two column (RT, R) data frame, on row for each trial.
  # Sigma and tau can be vectors with one value, or one value for each accumulator
  # or a matrix in the same form as mu. NB: CAN PRODUCE NEGATIVE RTS!
  dt <- matrix(rnorm(n=length(mu), mean=mu,sd=sigma) +
                 rexp(n=length(mu),rate=1/tau),nrow = dim(mu)[1])
  winner <- apply(dt,2,which.min)
  data.frame(RT=dt[cbind(winner,1:dim(mu)[2])], R=winner)
}


#' @importFrom stats rnorm rbinom rexp runif
rexgss <- function (n, mu, sigma, tau, tf=0, gf=0,
                    SSD=Inf, staircase=NA)
{
  # Race among n accumulators, first of which is a stop accumulator.
  # NA returned for RT when winner = 1. Optional SSD arguement can be used to
  # adjust mu for first accumulator to mu+SSD. SSD can be a scalar or vector
  # length n. For trials with winning first accumulator RT and R set to NA.
  # Adds SSD column to output.
  # tf = trigger failure probability, gf = go failure probability
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

#' @importFrom stats rnorm rgamma rnorm
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

#' @importFrom stats rnorm rgamma
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

#' @importFrom stats rnorm rlnorm
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
