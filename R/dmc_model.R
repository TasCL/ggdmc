# System functions for the DMC (Dynamic Models of Choice)
#    Usually user does not need to edit

make.level.array <- function(factors=list(S=c("s1","s2")))
# Create array of all factorial combinations of factor levels
{
  all.levels <- factors[[1]]

  if ( length(factors)>1 ) {
    for (i in 2:length(factors))
      all.levels <- outer(all.levels,factors[[i]],"paste",sep=".")
    dimnames(all.levels) <- factors
  } else
    all.levels <- array(all.levels,dim=length(factors[[1]]),dimnames=factors)

  all.levels
}



# p.map = list(A="1",B=c("R"),t0="1",
#   mean_v=c("A","S","M"),sd_v="M",st0="1")
#   match.map = list(M=list(mirror="MIRROR",normal="NORMAL"))
#   factors=list(A=c("00","45","90","135","180"),S=c("mirror","normal"))
#   constants = c(sd_v.false=1,st0=0,
#     mean_v.00.mirror.false=0,mean_v.45.mirror.false=0,
#     mean_v.90.mirror.false=0,mean_v.135.mirror.false=0,
#     mean_v.180.mirror.false=0,mean_v.00.normal.false=0,
#     mean_v.45.normal.false=0,mean_v.90.normal.false=0,
#     mean_v.135.normal.false=0,mean_v.180.normal.false=0)
#   responses = c("NORMAL","MIRROR")
#   type="norm"

# p.map=list(t0="1",A="1",B=c("E","R"),mean_v=c("E","V"),sd_v="1",st0="1")
# factors=list(E=c("nc", "wc"),S=c("n","w","pn","pw"))
# responses=c("N","W","P")
# match.map=list(M=list(n="N",w="W",pw="P",pn="P"),V=map3)
# constants=c(sd_v=1,st0=0)
# type="norm"; posdrift=TRUE

#' Creating a Model Object
#'
#' \code{model.dmc} Creates an array object ("model") with a set of attributes
#' specifying a particular model and parameterisation. Call \pkg{coda} to
#' summarise the model parameters in a DMC samples with multiple participants
#' at the hyper level.
#'
#' \code{model.dmc} creates a matrix used by \code{get.par.mat} to arrange
#' elements of a parameter vector appropriately. The attributes of output used by
#' \code{get.par.mat} to add in constants and check \code{transform.par}
#' creates the right parameters.
#'
#' @param p.map list factors and constants for parameters
#' @param factors Factor names and levels
#' @param responses Response (accumulator) names
#' @param match.map Scores responses
#' @param constants Parameters set to constant value
#' @param type model type. Should go to class in the future
#' @param posdrift only used by norm (ie LBA model)
#' @param verbose Print p.vector, constants and type
#' @keywords model.dmc
#' @export
#' @examples
#' m1 <- model.dmc(
#'    p.map     = list(a="1",v="F",z="1",d="1",sz="1",sv="1",t0="1",st0="1"),
#'    match.map = list(M=list(s1="r1",s2="r2")),
#'    factors=list(S=c("s1","s2"), F=c("f1","f2")),
#'    constants = c(st0=0,d=0),
#'    responses = c("r1","r2"),
#'    type = "rd")
model.dmc <- function(p.map, factors=list(S=c("s1","s2")),
    responses=c("r1","r2"), match.map=list(M=list(s1=1,s2=2)),
    constants=numeric(0), type="norm", posdrift=TRUE, verbose=TRUE)
{

  # Check responnses
  if ( type =="rd"  & ( length(responses)!=2 ) )
     stop("DDM only applicable for two responses")

  # Check match.map contains at least name M
  if ( length(match.map)<1 || class(match.map[[1]]) != "list" )
    stop("match.map must be a list of lists")

  if ( !any(names(match.map) %in% "M") )
    stop("match.map has no list named M")

  map.names <- names(match.map)[names(match.map)!="M"]
  map.levels <- unlist(lapply(match.map[names(match.map)!="M"],levels))

  # convert match.map$M to responses and check
  if ( is.numeric(unlist(match.map$M)) )
    match.map$M <- lapply(match.map$M,function(x){responses[x]})
  if ( !all(unlist(match.map$M) %in% responses) )
    stop("match.map$M has index or name not in response names")
  if ( !(all(sort(responses)==sort(unique(unlist(match.map$M))))) )
    stop("Not all response names are scored by match.map$M")
  if ( length(match.map)>1 &&
    !all(lapply(match.map[names(match.map)!="M"],class)=="factor") )
    stop("Entries in match.map besides M must be factors")

  # Check factors and add resposne
  if ( any(names(factors) %in% c("1","R","M","s")) )
    stop("Do not use s, M, R or 1 as a factor name")
  if ( any(names(factors) %in% names(match.map)) )
    stop(paste(map.names,"used in match.map, can not use as a factor name"))
  if ( length(unlist(factors)) != length(unique(unlist(factors))) )
    stop("All factors levels must be unqiue")
  if ( length(unlist(map.levels)) != length(unique(unlist(map.levels))) )
    stop("All match.map levels must be unqiue")
  if ( any(unlist(factors) %in% c("true","false")) )
    stop("\"true\" and \"false\" cannot be used as factor levels")
  if ( any(map.levels %in% c("true","false")) )
    stop("\"true\" and \"false\" cannot be used as match.map levels")
  if ( length(unlist(c(factors,map.levels))) !=
         length(unique(unlist(c(factors,map.levels)))) )
    stop("Factor levels cannot overlap match.map levels")

  # Check no parameter names have a dot
  has.dot <- unlist(lapply(strsplit(names(p.map),".",fixed=TRUE),length))>1
  if ( any(has.dot) )
    stop(paste("Dots not allowed in p.map names, fix:",paste(names(p.map)[has.dot])))

  # Check M and R are last
  if (any(unlist(lapply(p.map,function(x){any(x=="M") && x[length(x)]!="M"}))))
      stop("M factors must always be last")
  if (any(unlist(lapply(p.map,function(x){any(x=="R") && x[length(x)]!="R"}))))
      stop("R factors must always be last")

  factors.short <- factors
  factors$R <- responses
  factors$M <- c("true","false")
  for (i in map.names)
    factors[[i]] <- levels(match.map[[i]])

  # protect against grep problems
  # for (i in unlist(factors)) if ( length(grep(i,unlist(factors)))!=1 )
  #   stop("Factor, response or map level is not unique or is substring of another level!" )

  # Make parameter names
  names.par <- character(0)
  for ( i in names(p.map) )
  {
    if ( length(p.map[[i]])==1 && p.map[[i]] == "1" ) new.names <- i else
    {
      new.names <- paste(i,factors[[p.map[[i]][1]]],sep=".")
      if ( length(p.map[[i]])>1 ) for ( j in 2:length(p.map[[i]]) )
        new.names <- as.vector(outer(
          new.names,factors[[p.map[[i]][j]]],"paste",sep="."
        ))
    }
    names.par <- c(names.par,new.names)
  }

  # Make level array for manifest design and accumulators
  level.array <- make.level.array(factors[1:(length(factors.short)+1)])

  # Check match.map and expand
  for ( i in names(match.map$M) )
  {
    match.position <- grep(i,level.array)
    if ( length(match.position)==0 )
      stop(paste(i,"in match.map is not in the design"))
  }

  # Does the par use the match factor?
  is.M <- unlist(lapply(p.map,function(x){
    any(unlist(strsplit(x,".",fixed=T)=="M"))
  }))

  # Does the par use a response factor
  is.R <- unlist(lapply(p.map,function(x){
    any(unlist(strsplit(x,".",fixed=T)=="R"))
  }))

  if ( type =="rd"  & ( any(is.M) | any(is.R) ) )
     stop("Cannot use M or R in DDM p.map")

  # Does the par use a map factor
  is.map <- unlist(lapply(p.map,function(x){
    any(unlist(strsplit(x,".",fixed=T) %in% map.names))
  }))

  if ( any(is.map) ) {
    p.map.name <- lapply(p.map,function(x){
      unlist(strsplit(x,".",fixed=T))[
        unlist(strsplit(x,".",fixed=T)) %in% map.names]
    })
    nr <- length(responses)
    n <- length(level.array)
    map.shuffle <- matrix(aperm(array(1:n,dim=c(n/nr,nr,nr)),c(1,3,2)),ncol=nr)
  }

  if ( any(apply(cbind(is.M,is.R,is.map),1,function(x){sum(x)>1})) )
    stop("Parameters cannot have more than one of match.map and R factors")

  # use.par = boolean matrix for parameter use, cells x pars x resposnes
  use.par <- array(NA,
    dim=c(length(level.array),length(names.par),length(responses)))
  dimnames(use.par) <-
    list(as.vector(level.array),names.par,responses)

  # col.par = column parameter type (1st name)
  col.par <- strsplit(dimnames(use.par)[[2]],".",fixed=T)
  col.fac <- lapply(col.par,function(x){x[-1]})
  col.par <- unlist(lapply(col.par,function(x){x[1]}))
  # split into fac and resp
  col.fac <- lapply(col.fac,function(x){
    if ( length(x)==0 ) out <- c(NA,NA)
    if ( length(x)==1 ) {
      if ( x %in% c(responses,"true","false",map.levels) )
        out <- c(NA,x) else out <- c(x,NA)
    }
    if ( length(x)>1 )
      if ( x[length(x)] %in% c(responses,"true","false",map.levels) )
        out <- c(paste(x[-length(x)],collapse="."),x[length(x)]) else
        out <- paste(x,collapse=".")
    out
  })
  col.resp <- unlist(lapply(col.fac,function(x){x[2]}))
  col.fac <- unlist(lapply(col.fac,function(x){x[1]}))

  row.fac <- strsplit(dimnames(use.par)[[1]],".",fixed=T)
#  row.resp <- unlist(lapply(row.fac,function(x){x[length(x)]}))
  row.fac <- unlist(lapply(row.fac,function(x){
    paste(x[-length(x)],collapse=".")}))

  # Fill use.par array
  for ( p in unique(col.par) )
  { # parameters
    is.col <- p==col.par
    ncols <- sum(is.col)
    if ( ncols==1 ) use.par[,is.col,] <- TRUE else
    { # there are parameter subtypes
      for ( i in 1:ncols )
      { # each parameter subtype
        # select rows based on factors
        tmp <- col.fac[is.col][i]
        is.fac.row <- rep(TRUE,dim(use.par)[1])
        if ( !is.na(tmp) ) is.fac.row[!grepl(tmp,row.fac)] <- FALSE
        # set not applicable rows to false
        use.par[!is.fac.row,is.col,][,i,] <- FALSE
        if ( is.M[p] )
        { # has a match factor
          for ( j in names(match.map$M) )
          { # response cell
            correct.response <- match.map$M[[j]]
            is.rcell <- is.fac.row & grepl(j,row.fac)
            for ( k in responses )
            { # responses
               if ( k==correct.response )
               {
                 if ( grepl("true",col.resp[is.col][i]) )
                   use.par[,is.col,][is.rcell,i,k] <- TRUE else
                   use.par[,is.col,][is.rcell,i,k] <- FALSE
               } else {
                 if ( grepl("false",col.resp[is.col][i]) )
                   use.par[,is.col,][is.rcell,i,k] <- TRUE else
                   use.par[,is.col,][is.rcell,i,k] <- FALSE
               }
            }
          }
        } else if ( is.R[p] ) {
          for ( k in responses )
            use.par[is.fac.row,is.col,k][,i] <- k==col.resp[is.col][i]
        }  else if ( is.map[p] ) {
           use.par[is.fac.row,is.col,][,i,] <-
             match.map[[ p.map.name[[p]] ]] [map.shuffle[is.fac.row,]]==col.resp[is.col][i]
        } else use.par[is.fac.row,is.col,][,i,] <- TRUE
      }
    }
  }

  if ( any(is.na(use.par)) )
    stop("Some cells of the map were not assigned!")

  # add in constants
  all.par <- use.par[1,,1]
  all.par[1:length(all.par)] <- NA
  if (length(constants)>0) {
    if ( !any(names(constants) %in% names(all.par)) )
      stop("Name(s) in constants not in p.map")
    all.par[names(constants)] <- constants
  }
  attr(use.par,"all.par") <- all.par
  attr(use.par,"p.vector") <- all.par[is.na(all.par)]

  if (length(attr(use.par,"p.vector"))<2)
    stop("DMC cant handle models with less than two parameters")

  if (verbose) {
    cat("\nParameter vector names are: ( see attr(,\"p.vector\") )\n")
    print(names(all.par[is.na(all.par)]))
    cat("\nConstants are (see attr(,\"constants\") ):\n")
    print(constants)
    mod <- paste("\nModel type =",type)
    if (type=="norm") mod <- paste(mod,"(posdrift=",posdrift,")")
    cat(paste(mod,"\n\n"))
  }

  # Names of parameter types (cannot have a period)
  attr(use.par,"par.names") <- unique(col.par)
  attr(use.par,"type.par.names") <- switch(type,
    norm   =c("A","b","t0","mean_v","sd_v","st0"),
    normN  =c("A","b","t0","mean_v","sd_v","st0","N"),
    normgng=c("A","b","t0","mean_v","sd_v","st0"),
    normgamma=c("A","b","t0","mean_v","sd_v","shape","scale"),
    normgammaMR=c("A","b","t0","mean_v","sd_v","shape","scale","AS"),
    normlnorm=c("A","b","t0","mean_v","sd_v","meanlog","sdlog"),
    lnr   =c("meanlog","sdlog","t0","st0"),
    lnrdr =c("meanlog","sdlog","t0","st0"),
    lnrgng=c("meanlog","sdlog","t0","st0"),
    lnrss =c("meanlog","sdlog","meanlogS","sdlogS","t0","t0sg","tf","gf","ts"),
    exgss=c("mu","sigma","tau","muS","sigmaS","tauS","tf","gf"),
    rd=c("a","v","t0","z","d","sz","sv","st0")
  )
  attr(use.par,"type") <- type
  attr(use.par,"factors") <- factors.short
  attr(use.par,"responses") <- responses
  attr(use.par,"constants") <- constants
  attr(use.par,"posdrift") <- posdrift

  # save "n1" orders
  resp <- unlist(lapply(strsplit(level.array,".",fixed=TRUE),function(x){
    x[length(x)]}))
  nr <- length(responses)
  n1.order <- matrix(nrow=length(resp),ncol=nr)
  for (i in 1:length(resp))
    n1.order[i,] <- c(c(1:nr)[resp[i]==responses],c(1:nr)[resp[i]!=responses])
  row.names(n1.order) <- row.names(use.par)

  # Boolean for matching cells
  match.cell <- logical(length(row.names(n1.order)))
  names(match.cell) <- row.names(n1.order)
  for (i in 1:length(match.map$M)) {
    match.num <- grep(match.map$M[i],names(match.cell))
    match.cell[match.num[match.num %in%
        grep(names(match.map$M[i]),names(match.cell))]] <- TRUE
  }

  attr(use.par,"n1.order") <- n1.order
  attr(use.par,"match.cell") <- match.cell
  attr(use.par,"match.map") <- match.map

  if (type=="rd") # r1 cells have z flipped
  {
    is.r1 <- rep(FALSE,length(row.names(use.par)))
    names(is.r1) <- row.names(use.par)
    is.r1[unlist(lapply(
      lapply(
        as.list(names(match.map$M)[match.map$M==responses[1]]),
        function(x)grepl(x,row.names(use.par))
      ),
      function(x)which(x==TRUE)))] <- TRUE
    attr(use.par,"is.r1") <- is.r1
  }

  ## Add class ##
  attr(use.par, "class") <- switch(type,
    norm        = c("array", "dmc", "norm"),
    normN       = c("dmc", "normN"),
    normgng     = c("dmc", "normgng"),
    normgamma   = c("dmc", "normgamma"),
    normgammaMR = c("dmc", "normgammaMR"),
    normlnorm   = c("dmc", "normlnorm"),
    lnr    = c("dmc", "lnr"),
    lnrdr  = c("dmc", "lnrdr"),
    lnrgng = c("dmc", "lnrgng"),
    lnrss  = c("dmc", "lnrss"),
    exgss  = c("dmc", "exgss"),
    rd     = c("array", "dmc", "rd"))

  return(use.par)
}


update.model.dmc <- function(model,p.map=NULL,constants=NULL,
  match.map=NULL,factors=NULL,responses=NULL,type=NULL,posdrift=NULL)
# updates model with changed arguements
{

  if (is.null(p.map)) p.map <- attr(model,"p.map")
  if (is.null(constants)) constants <- attr(model,"constants")
  if (is.null(match.map)) match.map <-attr(model,"match.map")
  if (is.null(factors)) factors <- attr(model,"factors")
  if (is.null(responses)) responses <- attr(model,"responses")
  if (is.null(type)) type <- attr(model,"type")
  if (is.null(posdirft)) posdrift <- attr(model,"posdrift")
  model.dmc(p.map=p.map,
            constants=constants,
            match.map=match.map,
            factors=factors,
            responses=responses,
            type=type
  )
}


# Is p.vector compatible with model?
check.p.vector <- function(p.vector,model)
{
  is.match <- names(attr(model,"p.vector")) %in% names(p.vector)
  bad <- any(!is.match)
  if ( any(!is.match) ) warning(paste("Parameter",
   names(attr(model,"p.vector"))[!is.match],"in model not present in p.vector\n"))
  is.match <- names(p.vector) %in% names(attr(model,"p.vector"))
  bad <- bad | any(!is.match)
  if ( any(!is.match) ) warning(paste("Parameter",
    names(p.vector)[!is.match],"in p.vector not present in model\n"))
  invisible(bad)
}


#' Print accumulator x internal parameter type matrix for each cell
#'
#' A DMC model array, created by \code{model.dmc} is computationally efficient
#' but it is not easy to see how parameters map to design cells.
#' \code{print.cell.p} makes it easy to loop through cells and print the
#' mapping.
#'
#' @param p.vector a parameter vector,
#' @param model a DMC model
#' @export
#' @examples
#' model <- model.dmc(
#' p.map     = list(A = "1", B = "1", mean_v = "M", sd_v = "M", t0 = "1",
#'                  st0 = "1"),
#' match.map = list(M = list(s1 = 1, s2 = 2)),
#' factors   = list(S = c("s1", "s2")),
#' constants = c(st0 = 0, sd_v.false=1),
#' responses = c("r1","r2"),
#' type      = "norm")
#' p.vector <- c(A=3, B=4, mean_v.true=2, mean_v.false=-1, sd_v.true=0.5, t0=.2)
#' print_cell_p(p.vector, model)
print_cell_p <- function(p.vector, model)
{
  for (i in 1:dim(model)[1])
  {
    print(row.names(model)[i])
    print(p.df.dmc(model, p.vector, i))
    cat("\n")
  }
}

#' Gets Parameter Data Frame
#'
#' Gets parameter data frame (one row for each accumulator) for a design cell
#' (specified by name or index) with model picking out the appropriate
#' elements of par, and function transform.par transforms them appropriately
#' for model specified in model.
#'
#' @param model a model, created by model.dmc
#' @param p.vector a parameter vector
#' @param cell a string indicating a design cell, e.g., \code{s1.f1.r1}
#' @param n1order a boolean switch to use specific LBA ordering for its
#' parameters
#' @return  rows in natural (r1, r2 etc., used by simulate.dmc) or "n1" order
#' (used by likelihood.dmc)
#' @export
#' @examples
#' m1 <- model.dmc(
#'   p.map     = list(a="1",v="F",z="1",d="1",sz="1",sv="1", t0="1",st0="1"),
#'   constants = c(st0=0,d=0),
#'   match.map = list(M=list(s1="r1",s2="r2")),
#'   factors   = list(S=c("s1","s2"), F=c("f1", "f2")),
#'   responses = c("r1","r2"),
#'   type      = "rd")
#' pVec <- c(a=1, v.f1=1, v.f2=1.5, z=0.5, sz=0.25, sv=0.2,t0=.15)
#' p.df.dmc(m1, pVec, "s1.f1.r1", n1order=TRUE)
p.df.dmc <- function(model, p.vector, cell, n1order=TRUE)
{
  # Fill in non-constants
  attr(model,"all.par")[is.na(attr(model,"all.par"))] <-
    p.vector[names(attr(model,"p.vector"))]

  # Make parameter matrix
  par.mat <- matrix(rep(attr(model,"all.par"), times=dim(model)[3])[model[cell,,]],
         byrow=TRUE,nrow=dim(model)[3])
  dimnames(par.mat) <- list(dimnames(model)[[3]],attr(model,"par.names"))

  # For rd flip z
  if ( attributes(model)$type=="rd" && attributes(model)$is.r1[cell] )
    par.mat[,dimnames(par.mat)[[2]]=="z"] <-
      1-par.mat[,dimnames(par.mat)[[2]]=="z"]

  if (n1order) { par.mat <- par.mat[attr(model,"n1.order")[cell,],] }

  ## For allowing the user to do data.frame operation
  ## e.g., par.df$b <- par.df$B + par.df$A
  ## We have to make par.mat as data.frame and out as data.frame
  ## For simulation.dmc to aslo do data .frame operation
  ## Note most data.frame operations are slow and consume lots of memory
  ## One solution is to ask the user to learn data.table operations
  out <- transform(model, data.frame(par.mat))
  return(data.frame(out))

}


#' Simulate Responses from an EAM
#'
#' Simulate one or more responses from an EAM model object. This creates a data
#' frame using the model.
#'
#' p.vector=ps[1,,drop=FALSE]; n=n[1,]; SSD=SSD[datr]
#' !!! TO DO !!!
#' !!! Add an update ability where n=data and RT is updated,
#' also on first creation add index to speed update
#'
#' @param object is a DMC model object.
#' @param nsim nsim can be a single number for a balanced design or a vector
#' (one number per cell) to create an unbalanced design.
#' @param seed for compatible to \pkg{stats} simulate. Default NA
#' @param p.vector a parameter vector. Use Lognromal LBA model as an example,
#' \code{pVec <- c(A=.25, B=.35, meanlog_v.true=1, meanlog_v.false=.25,
#' sdlog_v.true=.5,t0=.2)}
#' @param SSD is for use only with stop-signal designs, it must be a scalar, or
#' a vector the same length as the number of cells or the same length as the
#' data and have Inf in all go cells
#' @param staircase is used only in stop-signal designs, overrides SSD, and
#' simulates a tracking algorighm setting SSD=start initially then moving up
#' step for each stop fail and down step for each stop-success.
#' @param TRIALS is a trial number covariate used in stop-signal paradimgs to
#' account for slowing or speeding over the course of the experiment.
#' @param ... other arguments
#' @export
#' @importFrom stats simulate
#' @examples
#' model <- model.dmc(
#'    p.map     = list(a="1",v="1",z="1",d="1",sz="1",sv="1", t0="1",st0="1"),
#'    constants = c(st0=0,d=0),
#'    match.map = list(M=list(s1="r1",s2="r2")),
#'    factors   = list(S=c("s1","s2")),
#'    responses = c("r1","r2"),
#'    type      = "rd")
#'
#' pop.prior <- prior.p.dmc(
#'   dists = rep("tnorm", 6),
#'   p1    = c(a=2,   v=2.5, z=0.5, sz=0.3, sv=1,  t0=0.3),
#'   p2    = c(a=0.5, v=.5,  z=0.1, sz=0.1, sv=.3, t0=0.05),
#'   lower = c(0,-5, 0, 0, 0, 0),
#'   upper = c(5, 7, 2, 2, 2, 2))
#'
#' p.vector <- c(a=1,v=1, z=0.5, sz=0.25, sv=0.2,t0=.15)
#' d <- simulate(model, nsim=1e2, p.vector=p.vector)
#' head(d)
#' ## dplyr::tbl_dt(d)
#' ## Source: local data table [200 x 3]
#' ##        S      R        RT
#' ##    (fctr) (fctr)     (dbl)
#' ## 1      s1     r1 0.6339412
#' ## 2      s1     r2 0.5783174
#' ## 3      s1     r2 0.2005078
#' ## 4      s1     r1 0.2973437
#' ## 5      s1     r2 0.4195281
#' ## 6      s1     r2 0.1946740
#' ## 7      s1     r1 0.2696773
#' ## 8      s1     r1 0.3917966
#' ## 9      s1     r1 0.8721902
#' ## 10     s1     r1 0.2786452
#' ## ..    ...    ...       ...
simulate.dmc <- function(object, nsim=2, seed=NULL, p.vector=NULL, SSD=Inf,
  staircase=NA, TRIALS=NA, ...)
{
  if (is.null(p.vector)) stop("Must supply a p.vector\n")

  if (check.p.vector(p.vector, object)) stop()
  # create factor data frame
  facs <- lapply(strsplit(dimnames(object)[[1]],".",fixed=TRUE),
         function(x){x[-length(x)]})
  facs <- facs[1:(length(facs)/length(attr(object,"responses")))]
  fnams <- names(attr(object,"factors"))
  facs <- data.frame(t(matrix(unlist(facs),nrow=length(fnams))))
  names(facs) <- fnams

  # check nsim
  if ( length(nsim) ==1 ) nsim <- rep(nsim, dim(facs)[1])
  if ( length(nsim)!= dim(facs)[1] )
    stop(paste("nsim must either be length 1 or", dim(facs)[1],"for this model."))

  # create data data frame
  data <- data.frame(lapply(facs, rep, times=nsim))
  data$R <- NA
  data$RT <- NA
  row1 <- 1

  if ( attr(object,"type") %in% c("lnrss","exgss") ) # stop signal model
  {
    if (any(is.na(SSD))) stop("SSD cannot contain NAs!")
    if ( !(length(SSD) %in% c(1,length(nsim), dim(data))) )
      stop(paste("SSD must have length =",dim(facs)[1],
        "(number of cells) or length =",dim(data)[1],"(number of data points)"))
    if ( length(SSD)==1 ) {
      if ( any(is.na(staircase)) ) warning(
        paste("You have specified SSD as a scalar, unless it is Inf that\n",
        " means all cells have a stop signal, if it is Inf it means none do!")
      )
      SSD <- rep(SSD,dim(data)[1])
    }
    if ( length(SSD)==length(nsim) ) SSD <- rep(SSD, nsim)
    data$SSD <- NA
    SSD.name <- "SSD"
  } else SSD.name <- NULL
  if ( attr(object,"type") %in% c("lnrss") ) # stop signal with slowing/speeding
  {
    if (length(TRIALS)==1) {
      if (!is.na(TRIALS)) stop("TRIALS must be NA if length = 1 (default)")
      TRIALS.name <- NULL
    } else {
      if ( length(TRIALS) != dim(data)[1] | any(is.na(TRIALS)) )
      stop(paste("TRIALS cannot contain NAs and must have length =",
                 dim(data)[1],"(number of data points)"))
      TRIALS.name <- "TRIALS"
      data$TRIALS <- NA
    }
  } else TRIALS.name <- NULL


  for (i in 1:dim(facs)[1]) if ( nsim[i]>0 ) {
    ## object is a model
    p.df <- p.df.dmc(object, p.vector, i, n1order=FALSE)
    rown <- row1 + nsim[i]-1
    data[row1:rown,c("RT","R",SSD.name,TRIALS.name)] <- switch(attr(object,"type"),

      # rtdists
      norm=rlba.norm(nsim[i], A=p.df$A,b=p.df$b,t0=p.df$t0,
        mean_v=p.df$mean_v,sd_v=p.df$sd_v,st0=p.df$st0[1],
        posdrift = attr(object,"posdrift")),
      normN=rlba.norm(nsim[i], A=p.df$A[1:p.df$N[1]],b=p.df$b[1:p.df$N[1]],
        t0=p.df$t0[1:p.df$N[1]],mean_v=p.df$mean_v[1:p.df$N[1]],
        sd_v=p.df$sd_v[1:p.df$N[1]],st0=p.df$st0[1],
        posdrift = attr(object,"posdrift")),
      rd=rdiffusion(nsim[i], a=p.df$a[1],v=p.df$v[1],t0=p.df$t0[1],
        sz=p.df$sz[1]*p.df$a[1],z=p.df$z[1]*p.df$a[1], # convert to absolute
        d=p.df$d[1],sv=p.df$sv[1],st0=p.df$st0[1]),


      # rtdists_extras
      lnr=      rlnr(nsim[i],meanlog=p.df$meanlog,sdlog=p.df$sdlog,
                     t0=p.df$t0,st0=p.df$st0[1]),

      lnrdr=    rlnr(nsim[i],meanlog=p.df$meanlog,sdlog=p.df$sdlog,
                     t0=p.df$t0,st0=p.df$st0[1]),
      # GNG models
      lnrgng=rlnrgng(nsim[i],meanlog=p.df$meanlog,sdlog=p.df$sdlog,
                     t0=p.df$t0, st0=p.df$st0[1]),
      normgng=rlba.normgng(nsim[i],A=p.df$A,b=p.df$b,t0=p.df$t0,
        mean_v=p.df$mean_v,sd_v=p.df$sd_v,st0=p.df$st0[1],
        posdrift = attr(object,"posdrift")),

      # Stop-signal models
      lnrss=rlnrss(nsim[i],
                   meanlog=p.df$meanlog,sdlog=p.df$sdlog,
                   t0=p.df$t0[1],t0sg=p.df$t0sg[1],
                   tf=p.df$tf[1],gf=p.df$gf[1],
                   ts=p.df$ts[1],TRIALS=TRIALS[row1:rown],
                   SSD=SSD[row1:rown],staircase=staircase),
      exgss=rexgss(nsim[i],
                   mu=p.df$mu,sigma=p.df$sigma,tau=p.df$tau,
                   tf=p.df$tf[1],gf=p.df$gf[1],
                   SSD=SSD[row1:rown],staircase=staircase),

      # LBA with gamma non-decision time
      normgamma=rlba.norm.gamma(nsim[i],A=p.df$A,b=p.df$b,t0=p.df$t0,
                                mean_v=p.df$mean_v,sd_v=p.df$sd_v,
                                shape=p.df$shape[1],scale=p.df$scale[1]),

      # LBA with gamma non-decision time and linear AS effect on shape
      normgammaMR=rlba.norm.gammaMR(nsim[i],A=p.df$A,b=p.df$b,t0=p.df$t0,
                                mean_v=p.df$mean_v,sd_v=p.df$sd_v,
                                shape=p.df$shape[1],scale=p.df$scale[1],
                                AS=p.df$AS[1]),

      # LBA with lnorm non-decision time
      normlnorm=rlba.norm.lnorm(nsim[i],A=p.df$A,b=p.df$b,t0=p.df$t0,
                                mean_v=p.df$mean_v,sd_v=p.df$sd_v,
                                meanlog=p.df$meanlog[1],sdlog=p.df$sdlog[1])
    )
    row1 <- rown+1
  }
  data$R <- factor(data$R,levels=1:length(attr(object,"responses")),
                   labels=attr(object,"responses"))

  if ( attr(object,"type") == "rd" )
  # Flip responses for cells mapped to response 1
  {
    cell.names <- apply(data[,1:length(names(facs)),drop=F],1,paste,collapse=".")
    M <- attr(object,"match.map")$M
    R <- attr(object,"responses")
    for ( i in names(M) )
      if ( M[[i]] == R[1] )
        data[grep(i,cell.names),"R"] <- as.character(
          factor(as.character(data[grep(i,cell.names),"R"]),
            levels=R,labels=R[2:1]))
  }
  return(data)
}

#' Bind Data and Models
#'
#' Binding a data frame and a DMC model description. This function also check
#' if they are compatible and adding a \code{cell.index} and many other
#' attributes to the data frame in order to speed likelihood computation.
#'
#' \code{data.model.dmc} adds a \code{dmc} class.
#'
#' @param data a data frame stored choice-RT data
#' @param model a DMC model
#' @export
#' @examples
#' m1 <- model.dmc(
#'    p.map     = list(a="1",v="1",z="1",d="1",sz="1",sv="1", t0="1",st0="1"),
#'    constants = c(st0=0,d=0),
#'    match.map = list(M=list(s1="r1",s2="r2")),
#'    factors   = list(S=c("s1","s2")),
#'    responses = c("r1","r2"),
#'    type      = "rd")
#'
#' pVec <- c(a=1,v=1, z=0.5, sz=0.25, sv=0.2,t0=.15)
#' dat1 <- simulate(m1, nsim=1e2, p.vector=pVec)
#' mdi1 <- data.model.dmc(dat1, m1)
#' class(mdi1)
#' ## [1] "data.frame"        "dmc"
data.model.dmc <- function(data, model)
{
  # check data
  if ( !is.data.frame(data) ) stop("data must be a data frame")
  fnams <- names(attr(model,"factors"))
  if ( !all(c(fnams,"R","RT") %in% names(data)) )
    stop(paste("data must have columns named:",
               paste(fnams,collapse=" "),"R","RT"))
  for ( i in names(attr(model,"factors")) )
    if ( !all(sort(attr(model,"factors")[[i]]) == sort(levels(data[,i]))) )
      stop(paste("Factor",i,"must have levels:",
                 paste(attr(model,"factors")[[i]],collapse=" ")))
  if ( !all(sort(attr(model,"responses")) == sort(levels(data[,"R"]))) )
      stop(paste("data$R must have levels:",
                 paste(attr(model,"responses"),collapse=" ")))
  if ( !is.numeric(data$RT) )
    stop("data$RT must be of type numeric")

  if ( any(names(data)=="s") ) # more than one subject
  {
    dat <- data
    data <- vector(mode="list",length=length(levels(dat$s)))
    names(data) <- levels(dat$s)
    is.sim <- !is.null(attr(dat,"parameters"))
    for (s in names(data))
    {
      data[[s]] <- dat[dat$s==s,names(dat)!="s"]
      # add model and index attribute to data
      cells <- apply(data[[s]][,c(fnams,"R")],1,paste,collapse=".")
      cell.index <- vector(mode="list",length=dim(model)[1])
      names(cell.index) <- row.names(model)
      for ( i in names(cell.index) )
        cell.index[[i]] <- cells %in% i
      attr(data[[s]],"cell.index") <- cell.index
      attr(data[[s]],"cell.empty") <-
        unlist(lapply(cell.index,function(x){sum(x)}))==0
      attr(data[[s]],"model") <- model
      if (is.sim) attr(data[[s]],"parameters") <- attr(dat,"parameters")[s,]
    }

    ## Add multiple classes
    attr(data, "class") <- c("list", "dmc")

  } else { # one subject
    # add model and index attribute to data
    cells <- apply(data[,c(fnams,"R")],1,paste,collapse=".")
    cell.index <- vector(mode="list",length=dim(model)[1])
    names(cell.index) <- row.names(model)
    for ( i in names(cell.index) )
      cell.index[[i]] <- cells %in% i
    attr(data,"cell.index") <- cell.index
    attr(data,"cell.empty") <-
      unlist(lapply(cell.index,function(x){sum(x)}))==0
    attr(data,"model") <- model

    ## Add classes
    attr(data, "class") <- c("data.frame", "dmc")
  }

  return(data)
}



######### CUSTOM MAP MAKING

empty.map <- function(FR,levels)
{
  level.array <- make.level.array(FR)
  map <- rep("",length=length(level.array))
  names(map) <- level.array
  factor(map,levels=levels)
}


assign.map <- function(map,value="",eq.list=list(),
                       funs=NULL,include=NA,match.values=NA)
{

  if (any(is.na(include))) ok <- rep(TRUE,length(map)) else
  {
    ok <- grepl(include[1],names(map))
    if (length(include)>1) for (i in 2:length(include))
      ok <- ok | grepl(include[i],names(map))
  }
  if ( length(eq.list)==0 ) # Fill in empty elements if included
    map[is.na(map) & ok] <- value else if (all(unlist(lapply(eq.list,length))==1))
  {
    if (all(is.na(match.values)) || length(match.values)!=length(eq.list))
      stop("If eq.list only has length one entries match.value must contain target for each entry")
    match.mat <- matrix(unlist(lapply(eq.list,function(x){
      (unlist(lapply(strsplit(names(map),".",fixed=T),function(y){
        y[x]})))})),ncol=length(eq.list))
    map[apply(match.mat,1,function(x){all(x==match.values)})] <- value
  } else {
    if ( is.null(funs) ) funs <- lapply(eq.list,function(x){
      list("identity","identity")
    })
    if (length(funs)!=length(eq.list))
      stop("Must specify one function pair in funs for each entry in eq.list")
    if ( class(funs)!="list" || !all(unlist(lapply(funs,class))=="list") )
      stop("funs must be  list of lists")
    if ( !all(unlist(lapply(funs,length))==2) )
      stop("Each entry in funs must be length 2")
    funs <- lapply(funs,function(x){lapply(x,function(y){
      if ( is.character(y) && y=="" ) "identity" else y
    })})
    pair.list <- lapply(eq.list,function(x){
      matrix(unlist(lapply(strsplit(names(map),".",fixed=T),function(y){
        y[x]})),nrow=2)})
    map.mat <- matrix(nrow=length(map),ncol=length(eq.list))
    for ( i in 1:length(eq.list) )
      map.mat[,i] <- apply(pair.list[[i]],2,function(x){
        do.call(funs[[i]][[1]],list(x[1])) ==
        do.call(funs[[i]][[2]],list(x[2]))
      })
      map[apply(map.mat,1,all) & ok] <- value
  }
  map
}

# map <- empty.map()
# map <- assign.map( map,value="true",eq.list=list(c(1,2)) )
# map <- assign.map(map,value="false")
