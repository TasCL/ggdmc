#' Transform Parameter Data Frame
#'
#' Transforms parameters to a form suitbale for the model been used
#' The function is called inside of \code{get.par.mat}.
#' \code{par.df} is a data frame of parameters types, some of which may need
#' to be transformed, or new columns created, so that the full set of
#' internal parameter types, specified in \code{type.par.names} attribute in
#' a model data instance, required by the type of evidence accumulation model
#' being used, such as the \code{norm} class in one of the LBA models, is
#' present.
#'
#' \code{transform.norm} template sets up for n-choice LBA, \eqn{B=b-A}
#' parameterization.
#'
#' External parameters types: \code{A}, \code{B}, \code{t0}, \code{mean_v},
#' \code{sd_v}, \eqn{st0 = 0} (optional).
#' Internal parameters types: \code{A}, \code{b}, \code{t0}, \code{mean_v},
#' \code{sd_v}, \eqn{st0 = 0} (optional).
#'
#' @param mdi a model data instance
#' @param par.df a data frame
#' @param ... other arguments
#' @keywords transform
#' @export
transform <- function(mdi, par.df, ...)
{
  UseMethod("transform", mdi)
}


#' @export
#' @rdname transform
transform.rd <- function(mdi, par.df, ...)
{

  type.par.names <- attr(mdi, "type.par.names")
  ## print("Use rd transform.dmc method.")

  #   # COMMENT OUT this check for speed after debugging
  #   if ( !all(type.par.names %in% names(par.df)) )
  #     stop("Trasform has not created parameter(s) required by the model.")

  ## mdi <- par.df
  par.df[, type.par.names]
}

#' @export
#' @rdname transform
transform.norm <- function(mdi, par.df, ...)
{
  type.par.names <- attr(mdi, "type.par.names")

  # User supplied tranforms go here
  par.df$b <- par.df$B + par.df$A

  # COMMENT OUT this check for speed after debugging
  # if ( !all(type.par.names %in% names(par.df)) )
  # stop("Trasform has not created parameter(s) required by the model.")
  par.df[, type.par.names]
}

#' Calculate Log-Likelihood
#'
#' \code{likelihood.rd} is a legacy function for computing diffusion
#' model density. If the speed of computation is a concern, the user may want to
#' use \code{ddmc}, which implemnts a streamlined C++ diffusion model density.
#'
#' \code{likelihood.norm} calculates LBA probability density.
#'
#' The precision for computing DDM density is usually set in range of 2
#' (fast but less accruate) to 3 (slower but more accurate). The default is
#' 2.5. Use \code{ddmc_parallel} if a high precision is needed or fitting l
#' arge samples. Set \coce{precision) (for DDM) or \code{ok.types} (for LBA) by
#' sending them via ... arguments.
#'
#' @param data a model data instance
#' @param p.vector a parameter vector
#' @param min.like minimal likelihood, 1e-10
#' @param ... other parameters, such as \emph{precision} for \code{ddm},
#' \emph{ok.types} for \code{lba_B}.
#' @keywords likelihood
#' @return a vector of likelihoods for each RT in data (in same order)
#' @export
likelihood <- function(data, p.vector, min.like=1e-10, ...)
{
  UseMethod("likelihood", data)
}

#' @rdname likelihood
#' @export
likelihood.default <- function(data, p.vector, min.like=1e-10, ...)
{
  print("Call likelihood default.")
}


#' @rdname likelihood
#' @importFrom utils glob2rx
#' @importFrom rtdists ddiffusion
#' @export
likelihood.rd <- function(data,  p.vector, min.like=1e-10, ...)
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

#' @importFrom rtdists n1PDF
#' @rdname likelihood
#' @export
likelihood.norm <- function(data, p.vector, min.like=1e-10, ... )
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


