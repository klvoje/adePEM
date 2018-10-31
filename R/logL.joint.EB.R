logL.joint.EB<-function (p, y) 
{
  anc <- p[1]
  vs <- p[2]
  alpha <- p[3]
  n <- length(y$mm)
  VV <- vs * outer((exp(alpha*y$tt) -1)/alpha, (exp(alpha*y$tt) -1)/alpha, FUN = pmin)
  diag(VV) <- diag(VV) + y$vv/y$nn
  M <- rep(anc, n)
  S <- dmnorm(y$mm, mean = M, varcov = VV, log = TRUE)
  return(S)
}