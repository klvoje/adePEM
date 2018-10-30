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

