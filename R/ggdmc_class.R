dmc <- function (samples=NA)
{
  attr(samples, "class") <- "dmc"
  return(samples)
}

dmc.list <- function (samples=NA)
{
  class(samples) <- "dmc.list"
  return(samples)
}

hyper <- function (samples=NA)
{
  class(samples) <- "hyper"
  return(samples)
}

pp.ggdmc <- function (samples=NA)
{
  class(samples) <- "pp.ggdmc"
  return(samples)
}

ddm <- function(object=NA)
{
  class(object) <- "ddm"
  return(object)
}

lba_B <- function(object=NA)
{
  class(object) <- "lba_B"
  return(object)
}

