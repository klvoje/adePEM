#' @title Simulate early burst time-series
#'
#' @description Simulate the evolution of a trait according to the early burst model.
#'
#' @param ns number of samples in time-series
#' 
#' @param alpha parameter describing the decreasing rate of change through time.
#'
#' @param vs variance of the step distribution.
#' 
#' @param vp within-population trait variance.
#'
#' @param nn vector of the number of individuals in each sample
#' 
#' @param tt vector of sample ages, increases from oldest to youngest
#'
#' @details The functin simulate the evolution of a trait according to the early burst model: A model with 
#' initial higher rates of evolution (early burst) followed by an exponential drop in the rate of evolution. 
#' See reference below for details on the model. 
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
#'## generate a paleoTS objects by simulating early burst
#'x <- sim.EB(ns=20, alpha=-1, vs=0.1)
#'
#'## Investigate if the time series pass all thee adequacy tests
#'fit3adequacy.EB(x)
#'

sim.EB<-function (ns = 20, alpha=-1, vs = 0.1, vp = 0.01, nn = rep(20, ns), tt = 0:(ns - 1)) 
{
  MM <- array(dim = ns)
  mm <- array(dim = ns)
  vv <- array(dim = ns)
  dt <- diff(tt)
  inc <- rnorm(ns - 1, 0, sqrt(vs*(exp(alpha*tt[2:ns]))))
  MM <- cumsum(c(0, inc))
  mm <- MM + rnorm(ns, 0, sqrt(vp/nn))
  vv <- rep(vp, ns)
  gp <- c(vs, alpha)
  names(gp) <- c("vstep", "alpha")
  res <- as.paleoTS(mm = mm, vv = vv, nn = nn, tt = tt, MM = MM, 
                    genpars = gp, label = "Created by sim.EB()", reset.time = FALSE)
  return(res)
}

