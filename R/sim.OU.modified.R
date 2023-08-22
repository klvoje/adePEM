sim.OU.modified<-function (ns = 20, anc = 0, theta = 10, alpha = 0.3, vstep = 0.1, 
                           vp = 1, nn = rep(20, ns), tt = 0:(ns - 1)) 
{
  MM <- array(dim = ns)
  mm <- array(dim = ns)
  vv <- array(dim = ns)
  dt <- diff(tt)
  MM[1] <- anc
  x <- rnorm(nn[1], mean = MM[1], sd = sqrt(vp))
  mm[1] <- mean(x)
  if (length(x) == 1) vv[1] <- 0 else vv[1] <- var(x)
  for (i in 2:ns) {
    ex <- paleoTS:::ou.M(MM[i - 1], theta, alpha, dt[i - 1])
    vx <- paleoTS:::ou.V(vstep, alpha, dt[i - 1])
    MM[i] <- rnorm(1, ex, sqrt(vx))
    x <- rnorm(nn[i], mean = MM[i], sd = sqrt(vp))
    mm[i] <- mean(x)
    if (length(x) == 1) vv[i] <- 0 else vv[i] <- var(x)
  }
  gp <- c(anc, theta, alpha, vstep)
  names(gp) <- c("anc", "theta", "alpha", "vstep")
  res <- as.paleoTS(mm = as.vector(mm), vv = as.vector(vv), 
                    nn = nn, tt = tt, MM = MM, genpars = gp, label = "Created by sim.OU()", 
                    reset.time = FALSE)
  return(res)
}