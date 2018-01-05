net.change.test.stasis<-function(y, theta, omega, nrep=1000, conf=0.95, plot=TRUE, save.replicates=TRUE){

  x<-y$mm
  v<-y$vv
  n<-y$nn
  time<-y$tt

  lower<-(1-conf)/2
  upper<-(1+conf)/2

  obs.net.change.test<-net.change.test(x, model="stasis")

  ### Parametric bootstrap routine ###

  #Matrix that will contain the test statistic for each simluated data set (time series)
  bootstrap.matrix<-matrix(data = NA, nrow = nrep, ncol = 1)


  # parametric boostrap
  for (i in 1:nrep){

    x.sim<-sim.Stasis(ns = length(x), theta = theta, omega = omega, vp = v, nn = n, tt = time)

    bootstrap.matrix[i,1]<-net.change.test(x.sim$mm, model="stasis")

  }

  # Estimating the ratio of how often the observed net change is smaller than the net change in the simulated data
  bootstrap.net.test<-length(bootstrap.matrix[,1][bootstrap.matrix[,1]>obs.net.change.test])/nrep

  # Calculating the "p-value" and whether the observed data passed the test statistic
  if (bootstrap.net.test<(round(lower,3)*2)) pass.net.change.test<-"FAILED" else pass.net.change.test<-"PASSED"

  # Plot the test statistics estimated from the simulated data
  if (plot==TRUE) {
    layout(1:1)
    plot.distributions(bootstrap.matrix[,1],obs.net.change.test, test="net.change.test", xlab="Simulated data", main="Net evolution");
  }

  #Prepearing the outout
  output<-as.data.frame(cbind(round(obs.net.change.test,5), round(min(bootstrap.matrix),5), round(max(bootstrap.matrix),5), bootstrap.net.test, pass.net.change.test), nrow=5, byrow=TRUE)
  rownames(output)<-"net.change.test"
  colnames(output)<-c("estimate", "min.sim" ,"max.sim","p-value", "Result")


  summary.out<-as.data.frame(c(nrep, conf))
  rownames(summary.out)<-c("replications", "confidense level")
  colnames(summary.out)<-("Value")
  if (save.replicates==FALSE)
  {
    out<- list("info" = summary.out, "summary" = output)
    return(out)
  }
  else
  {
    out<- list("replicates" = bootstrap.matrix, "info" = summary.out, "summary" = output)
    return(out)
  }
}
