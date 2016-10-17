#' Transform Parameter Data Frame
#'
#' \code{transfrom.ddm} transforms parameters to a form suitbale for the model
#' being used. Called inside of get.par.mat. "par.df" is a data frame of
#' parameters types , some of which may need to be transformed, or new
#' columns created, so that the full set of internal parameter types, specified
#' in "type.par.names", required by the type of evidence accumulation model
#' being used ("norm" etc.) is present.
#'
#' \code{transform.lba_B} template setup for n-choice LBA,
#' B=b-A parameterization
#' \itemize{
#'   \item External parameters types: A, B, t0, mean_v, sd_v, st0 = 0 (optional)
#'   \item Internal parameters types: A, b, t0, mean_v, sd_v, st0 = 0 (optional)
#' }
#'
#' @param _data original called \code{par.df}, which is a data frame.
#' @param type.par.names default is NULL
#' @param ... other arguments
#' @keywords transform
#' @export
transform.ddm <- function( `_data`, type.par.names=NULL, ...)
{
  # User supplied tranforms go here
  #   par.df$b <- par.df$B+par.df$A

  #   # COMMENT OUT this check for speed after debugging
  #   if ( !all(type.par.names %in% names(par.df)) )
  #     stop("Trasform has not created parameter(s) required by the model.")

  ## _data <- par.df
  `_data`[, type.par.names]
}


#' Calculate Log-Likelihood
#'
#' \code{likelihood.ddm} is a legacy function for computing diffusion
#' model density. If the speed of computation is a concern, the user may want to
#' use \code{ddmc}, which implemnts a streamlined C++ diffusion model density.
#'
#' The precision for computing DDM density is usually set in range of 2
#' (fast but less accruate) to 3 (slower but more accurate). The default is
#' 2.5. Use \code{ddmc_parallel} if a high precision is needed or fitting l
#' arge samples. Set \coce{precision) (for DDM) or \code{ok.types} (for LBA) by
#' sending them via ... arguments.
#'
#' @param p.vector a parameter vector
#' @param data a model data instance
#' @param min.like minimal likelihood, 1e-10
#' @param ... other parameters, such as \emph{precision} for \code{ddm},
#' \emph{ok.types} for \code{lba_B}.
#' @keywords likelihood
#' @return a vector of likelihoods for each RT in data (in same order)
#' @export
likelihood.ddm <- function(p.vector, data, min.like=1e-10, ...)
{
  precision <- 2.5
  bad <- function(p)
    ## Stops ddiffusion crashing if given bad values, may need to add extra
    ## protection for v and upper bounds for a, sv and t0.
  {
    (p$a[1]<0)      | (p$z[1] <1e-6) | (p$z[1] >.999999) | (p$t0[1]<1e-6)  |
      (p$sz[1]<0) | (p$st0[1]<0)    | (p$sv[1]<0) |
      (p$sz[1]>1.999999*min(c(p$z[1],1-p$z[1]))) |
      (p$t0[1]-p$st0[1]/2<0)
  }

  bound <- rep("lower",dim(attr(data,"model"))[1])
  bound[as.vector(sapply(
    paste(names(attr(attributes(data)$model,"match.map")$M),
      attr(attributes(data)$model,"match.map")$M,sep="*"),
    function(x){grep(glob2rx(x),row.names(attr(data,"model")))}))] <- "upper"
  names(bound) <- row.names(attr(data,"model"))

  likelihood <- numeric(dim(data)[1])
  for ( i in row.names(attr(data,"model")) ) if ( !attr(data,"cell.empty")[i] )
  {
    p.df <- p.df.dmc(p.vector,i,attributes(data)$model,n1order=TRUE)
    # Avoid numerical errors
    if ( bad(p.df) ) likelihood[ attr(data,"cell.index")[[i]] ] <- 0 else
    {
      likelihood[ attr(data,"cell.index")[[i]] ] <-
        abs(ddiffusion(data$RT[attr(data,"cell.index")[[i]]],
          a=p.df$a[1],
          v=p.df$v[1],
          t0=p.df$t0[1],
          z=p.df$z[1]*p.df$a[1],   # convert to absolute z.
          #z=p.df$z[1],                # rtdists 0.6-6 will do converting
          d=p.df$d[1],
          sz=p.df$sz[1]*p.df$a[1], # convert to absolute sz.
          #sz=p.df$sz[1],              # rtdists 0.6-6 will do converting
          sv=p.df$sv[1],
          st0=p.df$st0[1],
          response=bound[i],
          precision=precision
        ))
    }
  }
  pmax(likelihood, min.like)
}

#' @export
#' @rdname transform.ddm
transform.lba_B <- function(`_data`, type.par.names=NULL, ...)
{
  par.df <- `_data`
  # User supplied tranforms go here
  par.df$b <- par.df$B+par.df$A

  # COMMENT OUT this check for speed after debugging
  # if ( !all(type.par.names %in% names(par.df)) )
  # stop("Trasform has not created parameter(s) required by the model.")
  par.df[,type.par.names]
}


#' @export
#' @rdname likelihood.ddm
likelihood.lba_B <- function(p.vector, data, min.like=1e-10, ... )

  # !!! TO DO: types other than norm
{
  ok.types <- c("norm")
  #   # COMMENT OUT this check for speed after debugging
  #   if ( !all(attr(model,"type") %in% ok.types) )
  #     stop("Distribution funciton type not supported by likelihood.dmc")

  likelihood <- numeric(dim(data)[1])
  for ( i in row.names(attr(data,"model")) ) if ( !attr(data,"cell.empty")[i] )
  {
    p.df <- p.df.dmc(p.vector,i,attributes(data)$model,n1order=TRUE)
    likelihood[ attr(data,"cell.index")[[i]] ] <-
      switch(attr(attributes(data)$model,"type"),
        norm=n1PDF(data$RT[attr(data,"cell.index")[[i]]],
          A=as.list(p.df$A),
          b=as.list(p.df$b),
          t0=p.df$t0[1],
          mean_v=p.df$mean_v,
          sd_v=p.df$sd_v,
          st0=p.df$st0[1],
          distribution=attr(attr(data,"model"),"type"),
          args.list = list(posdrift=attr(attr(data,"model"),"posdrift")),
          silent=TRUE)
      )
  }
  pmax(likelihood, min.like)
}


