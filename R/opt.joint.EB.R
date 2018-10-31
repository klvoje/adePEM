opt.joint.EB<-function (y, pool = TRUE, cl = list(fnscale = -1), meth = "L-BFGS-B", hess = FALSE) 
{
  if (pool) 
    y <- pool.var(y, ret.paleoTS = TRUE)
  if (y$tt[1] != 0) 
    stop("Initial time must be 0.  Use as.paleoTS() or read.paleoTS() to correctly process ages.")
  p0 <- array(dim = 3)
  p0[1] <- y$mm[1]
  p0[2] <- min(c(mle.URW(y), 1e-07))
  p0[3] <- -1
  names(p0) <- c("anc", "vstep", "alpha")
  if (is.null(cl$ndeps)) 
    cl$ndeps <- abs(p0/10000)
  cl$ndeps[cl$ndeps == 0] <- 1e-08
  if (meth == "L-BFGS-B") 
    w <- optim(p0, fn = logL.joint.EB, control = cl, method = meth, 
               lower = c(NA, 0, NA), upper = (c(NA, NA,0)), hessian = hess, y = y)
  else w <- optim(p0, fn = logL.joint.EB, control = cl, method = meth, lower = c(NA, 0, NA), 
                  upper = (c(NA, NA,0)), hessian = hess, y = y)
  
  if (hess) 
    w$se <- sqrt(diag(-1 * solve(w$hessian)))
  else w$se <- NULL
  wc <- as.paleoTSfit(logL = w$value, parameters = w$par, modelName = "EB", 
                      method = "Joint", K = 3, n = length(y$mm), se = w$se)
  return(wc)
}