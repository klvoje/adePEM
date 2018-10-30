#' @title Applying the constant variance test to the random walk model
#'
#' @description Investigates if the random walk model is an adequate statistical description of an evolutionary
#' time series by applying the constant variance test.
#'
#' @param y a paleoTS object
#'
#' @param vstep the variance of the step distribution estimated from the observed data.
#'
#' @param nrep number of iterations in the parametric bootstrap (number of simulated time series); default is 1000.
#'
#' @param conf confidence level for judging whether a model is an adequate statistical description of the data.
#' Number must be between 0 and 1. A higher number means less strict judgment of whether a model is adequate; default
#' is 0.95. Tests are two-tailed (except for the net evolution test), which means a model is judged adequate if the observed test statistic is within the 2.5
#' percent of the extreme values of the calculated test statistics on the simulated data given the default confidence
#' value of 0.95.
#'
#' @param plot logical; if TRUE, the value of the test statistic calculated based on the observed fossil
#' time series is plotted on the distribution of test statistics calculated on the simulated time series;
#' default is TRUE.
#'
#' @param save.replicates logical; if TRUE, the values of the test statistic calculated on the simulated time
#' series is saved and can be accessed later for plotting purposes; default is TRUE.
#'
#' @details Estimates the slope of the least square regression of the size of the detrended data (their absolute value) from the average
#' as a function of time.as a function of time.
#'
#' @return First part of the output summarizes the number of iterations in the parametric bootstrap and the
#' confidence level for judging whether a model is an adequate statistical description of the data. The last
#' part of the output is:
#'
#' @return
#'  \item{estimate}{The calculated test statistic on the observed data.}
#'  \item{min.sim}{The smallest test statistic calculated on the simulated data.}
#'  \item{max.sim}{The largest test statistic calculated on the simulated data.}
#'  \item{p-value}{Not a real p-value, but is calculated as the fraction of simulated test statistics
#'  that is larger (or smaller) than the calculated test statistic on the observed data divided by 0.5.
#'  A value of 1 means 50 percent of the test statistics on the simulated data are larger and smaller
#'  than the calculated statistic on the observed data. A value of 0.10 means 90 percent of the test
#'  statistics on the simulated data are larger or smaller than the test statistic on the observed time
#'  series.}
#'  \item{result}{Whether the model PASSED or FAILED the adequacy test. The outcome depends on the
#'  confidence level.}
#'
#'@author Kjetil L. Voje
#'
#'@references Voje, K.L., Starrfelt, J., and Liow, L.H. Model adequacy and microevolutionary explanations for stasis in the fossil record. \emph{The American Naturalist}. In press.
#'
#'@seealso \code{\link{fit3adequasy.RW}}, \code{\link{slope.test.stasis}}, \code{\link{slope.test.trend}}
#' @export
#'@examples
#'## generate a paleoTS objects by simulating a directional trend
#'x <- sim.GRW(ns=40, ms=0, vs=0.1)
#'
#'## estimate the variance of the step distribution
#'vstep <- mle.URW(x)[1]
#'
#'## investigate if the time series pass the adequacy test
#'slope.test.RW(x,vstep)
#'


slope.test.EB<-function(y, alpha, vstep, nrep=1000, conf=0.95, plot=TRUE, save.replicates=TRUE){

  x<-y$mm
  v<-y$vv
  n<-y$nn
  time<-y$tt

  lower<-(1-conf)/2
  upper<-(1+conf)/2
  lower.2<-(1-conf)
  upper.2<-conf

  obs.slope.test<-slope.test_EB(x,time, model="EB")

  ### Parametric bootstrap routine ###

  #Matrix that will contain the test statistic for each simulated data set (time series)
  bootstrap.matrix<-matrix(data = NA, nrow = nrep, ncol = 1)

  # parametric boostrap
  for (i in 1:nrep){

    x.sim<-sim.EB(ns=length(x), alpha=alpha, vs=vstep, vp=mean(v), nn=n, tt=time)

    bootstrap.matrix[i,1]<-slope.test_EB(x.sim$mm,time, model="EB")

  }

  # Estimating the ratio of how often the observed slope statistic is smaller than the slope tests in the simulated data
  bootstrap.slope.test<-length(bootstrap.matrix[,1][bootstrap.matrix[,1]>obs.slope.test])/nrep
  bootstrap.slope.test.2<-length(bootstrap.matrix[,1][bootstrap.matrix[,1]>0])/nrep


  # Calculating the "p-value" and whether the observed data passed the test statistic
  if (bootstrap.slope.test>round(upper,3) | bootstrap.slope.test<round(lower,3) | bootstrap.slope.test.2>round(lower.2,3)) pass.slope.test<-"FAILED" else pass.slope.test<-"PASSED"
  if(bootstrap.slope.test>0.5) bootstrap.slope.test<-1-bootstrap.slope.test

  # Plot the test statistics estimated from the simulated data
  if (plot==TRUE) {
    layout(1:1)
    plot.distributions(bootstrap.matrix[,1],obs.slope.test, test="slope.test", xlab="Simulated data", main="Reduced variance");
  }

  #Preparing the output
  output<-as.data.frame(cbind(round(obs.slope.test,5), round(min(bootstrap.matrix),5), round(max(bootstrap.matrix),5), bootstrap.slope.test/0.5, bootstrap.slope.test.2, pass.slope.test), nrow=6, byrow=TRUE)
  rownames(output)<-"slope.test"
  colnames(output)<-c("estimate","min.sim" ,"max.sim","p-value", "% > 0","result")


  summary.out<-as.data.frame(c(nrep, conf))
  rownames(summary.out)<-c("replications", "confidence level")
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
