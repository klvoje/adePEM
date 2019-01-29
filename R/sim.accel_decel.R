#' @title Simulate time-series where the rate of evolution increase or decrease with time
#'
#' @description Simulate the evolution of a trait according to a random walk where where the rate of evolution increase or decrease with time.
#'
#' @param ns number of samples in time-series
#' 
#' @param r parameter describing the decreasing/increasing rate of change through time.
#'
#' @param vs variance of the step distribution.
#' 
#' @param vp within-population trait variance.
#'
#' @param nn vector of the number of individuals in each sample
#' 
#' @param tt vector of sample ages, increases from oldest to youngest
#'
#' @details The function simulate the evolution of a trait where the rate of evolution increase or decrease with time:
#' The trait evolves accoridng to a random walk where the rate of evolution show either an exponential
#' increase or decrease with time. See reference below for details on the model. 
#'
#' @return First part of the output summarizes the number of iterations in the parametric bootstrap and the
#' confidence level for judging whether a model is an adequate statistical description of the data. The last
#' part of the output is a data frame with the adequacy tests as columns and the following rows:
#'
#' @return A paleoTS object.
#'
#'@author Kjetil L. Voje
#'
#'@references Harmon, L. J., J. B. Losos, T. J. Davies, R. G. Gillespie, J. L. Gittleman, W. B. Jennings, K. H. Kozak, et al. 2010. Early bursts of body size and shape evolution are rare in comparative data. \emph{Evolution} 64:2385–2396.
#'
#'@export
#'
#'@examples
#'## generate a paleoTS objects by simulating a trait that shows an increasing rate of evolution with time. 
#'x <- sim.accel_decel(ns=20, r=-1, vs=0.1)
#'
#'## Investigate if the time series pass all thee adequacy tests
#'fit3adequacy.accel(x)
#'

sim.accel_decel<-function (ns = 20, r=-5, vs = 0.1, vp = 0.01, nn = rep(20, ns), tt = 0:(ns - 1)/ns) 
{
  MM <- array(dim = ns)
  mm <- array(dim = ns)
  vv <- array(dim = ns)
  dt <- diff(tt)
  #inc <- rnorm(ns - 1, 0, sqrt(vs*(exp(r*tt[2:ns]))))
  #inc <- rnorm(ns - 1, 0, sqrt(vs*(exp(r*dt))))
  inc <- rnorm(ns - 1, 0, sqrt(vs*(exp(r*tt))))
  MM <- cumsum(c(0, inc))
  mm <- MM + rnorm(ns, 0, sqrt(vp/nn))
  vv <- rep(vp, ns)
  gp <- c(vs, r)
  names(gp) <- c("vstep", "r")
  res <- as.paleoTS(mm = mm, vv = vv, nn = nn, tt = tt, MM = MM, 
                    genpars = gp, label = "Created by sim.accel_decel()", reset.time = FALSE)
  return(res)
}

