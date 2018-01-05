#' @title Applying 4 adequasy tests to the stasis model
#'
#' @description Investigating if the stasis model is an adequate statistical description of an evolutionary
#' time series by applying the following tests (1) autocorrelation, (2) runs test, (3) constant variation, and (4)
#' net evolution.
#'
#' @param y a paleoTS object
#'
#' @param theta evolutionary optimum estmated from the observed data.
#'
#' @param omega evolutonary variance estmated from the observed data.
#'
#' @param nrep number of iterations in the parametric boostrap (number of simulated time series); default is 1000.
#'
#' @param conf confidence level for judging whether a model is an adequate statistical description of the data.
#' Number must be between 0 and 1. A higher number means less strict judgment of whether a model is adequate; default
#' is 0.95. Tests are two-tailed (except for the net evolution test), which means a model is judged adequate if the observed test statistic is within the 2.5
#' percent of the extreeme values of the calculated test statistics on the simulated data given the default confidence
#' value of 0.95.
#'
#' @param plot logical; if TRUE, the value of the test statistic calculated based on the observed fossil
#' time series is plotted on the distribution of test statistics calculated on the simulated time series;
#' default is TRUE.
#'
#' @details A wrapper function for investigating adequasy of the directional trend model
#' applying all three tests at the same time.
#'
#'
#' @return First part of the output summarizes the number of iterations in the parametric boostrap and the
#' confidence level for judging whether a model is an adequate statistical description of the data. The last
#' part of the output is a data frame with the adequasy tests as columns and the following rows:
#'
#' @return
#'  \item{estimate}{The calculated test statistic on the observed data.}
#'  \item{min.sim}{The smallest test statistic calculated on the simulated data.}
#'  \item{max.sim}{The largest test statistic calculated on the simulated data.}
#'  \item{p-value}{Not a real p-value, but is calculated as the fraction of simulated test statistics
#'  that is larger (or smaller) than the calculated test statistic on the observed data divided by 0.5.
#'  A value of 1 means 50 percent of the test statistics on the simulated data are largen and smaller
#'  than the calculated statistic on the observed data. A value of 0.10 means 90 percent of the test
#'  statistics on the simulated data are larger or smaller than the test statistic on the observed time
#'  series.}
#'  \item{result}{Whether the model PASSED or FAILED the adequasy test. The outcome depends on the
#'  confidence level.}
#'
#'@author Kjetil L. Voje
#'
#'@references Voje, K.L., Starrfelt, J., and Liow, L.H. Model adequacy and microevolutionary explanations for stasis in the fossil record. \emph{The American Naturalist}. In press.
#'
#'@seealso \code{\link{fit3adequasy.BM}}, \code{\link{fit3adequasy.trend}}
#' @export
#'@examples
#'## generate a paleoTS objects by simulating a stasis time series
#'x <- sim.Stasis(ns = 40, theta = 0, omega = 0.1)
#'
#'## estimate the evolutionary optimum
#'theta <- mle.Stasis(x)[1]
#'
#'## estimate the evolutionary variance
#'omega <- mle.Stasis(x)[2]
#'
#'## Investigate if the time series pass all four adequasy tests
#'fit4adequasy.stasis(x,theta,omega)
#'
fit4adequasy.stasis<-function(y, theta, omega, nrep=1000, conf=0.95, plot=TRUE){

  x<-y$mm
  v<-y$vv
  n<-y$nn
  time<-y$tt

  lower<-(1-conf)/2
  upper<-(1+conf)/2


  # Compute the test statistics for the observed time series
  obs.auto.corr<-auto.corr(x, model="stasis")
  obs.runs.test<-runs.test(x, model="stasis", theta=theta)
  obs.slope.test<-slope.test(x,time, model="stasis")
  obs.net.change.test<-net.change.test(x, model="stasis")

  #Run parametric bootstrap
  out.auto<-auto.corr.test.stasis(y,theta,omega, nrep, conf, plot=FALSE)
  out.runs<-runs.test.stasis(y,theta,omega, nrep, conf, plot=FALSE)
  out.slope<-slope.test.stasis(y,theta,omega, nrep, conf, plot=FALSE)
  out.net<-net.change.test.stasis(y,theta,omega, nrep, conf, plot=FALSE)

  #Prepearing the output
  output<-c(as.vector(matrix(unlist(out.auto[[3]]),ncol=5,byrow=FALSE)),
                        as.vector(matrix(unlist(out.runs[[3]]),ncol=5,byrow=FALSE)),
                        as.vector(matrix(unlist(out.slope[[3]]),ncol=5,byrow=FALSE)),
                        as.vector(matrix(unlist(out.net[[3]]),ncol=5,byrow=FALSE)))

  output<-as.data.frame(cbind(c(output[c(1,6,11,16)]), c(output[c(2,7,12,17)]), c(output[c(3,8,13,18)]), c(output[c(4,9,14,19)]), c(output[c(5,10,15,20)])), ncol=5)

  rownames(output)<-c("auto.corr", "runs.test", "slope.test", "net.change.test")
  colnames(output)<-c("estimate", "min.sim" ,"max.sim","p-value", "Result")

  if (plot==TRUE) {
    par(mfrow=c(2,2));
    model.names<-c("auto.corr", "runs.test", "slope.test", "net.change.test")
    plot.distributions(out.auto$replicates,obs.auto.corr, model.names[1], xlab="Simulated data", main="Autocorrelation");
    plot.distributions(out.runs$replicates,obs.runs.test, model.names[2], xlab="Simulated data", main="Runs");
    plot.distributions(out.slope$replicates,obs.slope.test, model.names[3], xlab="Simulated data", main="Fixed variance");
    plot.distributions(out.net$replicates,obs.net.change.test, model.names[4], xlab="Simulated data", main="Net evolution");
  }
  summary.out<-as.data.frame(c(nrep, conf))
  rownames(summary.out)<-c("replications", "confidense level")
  colnames(summary.out)<-("Value")
  out<- list("info" = summary.out, "summary" = output)
  return(out)
}
