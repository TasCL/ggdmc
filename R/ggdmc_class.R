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

rd <- function(object=NA)
{
  class(object) <- "rd"
  return(object)
}

lba_B <- function(object=NA)
{
  class(object) <- "norm"
  return(object)
}




