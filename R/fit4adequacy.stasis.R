#' @title Applying 4 adequacy tests to the stasis model
#'
#' @description Investigating if the stasis model is an adequate statistical description of an evolutionary
#' time series by applying the following tests (1) autocorrelation, (2) runs test, (3) constant variation, and (4)
#' net evolution.
#'
#' @param y a paleoTS object
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
#' @param omega evolutonary variance estimated from the observed data. This parameter is automatically estimated from the data, if not set 
#' by the user (usually not recommended).
#'
#' @details A wrapper function for investigating adequacy of the directional trend model
#' applying all three tests at the same time.
#'
#'
#' @return First part of the output summarizes the number of iterations in the parametric bootstrap and the
#' confidence level for judging whether a model is an adequate statistical description of the data. The last
#' part of the output is a data frame with the adequacy tests as columns and the following rows:
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
#'@references Voje, K.L. 2018. Assessing adequacy of models of phyletic evolution in the fossil record. \emph{Methods in Ecology and Evoluton}. (in press).
#'@references Voje, K.L., Starrfelt, J., and Liow, L.H. 2018. Model adequacy and microevolutionary explanations for stasis in the fossil record. \emph{The American Naturalist}. 191:509-523.
#'
#'@seealso \code{\link{fit3adequacy.RW}}, \code{\link{fit3adequacy.trend}}
#' @export
#'@examples
#'## Generate a paleoTS objects by simulating a stasis time series
#'x <- sim.Stasis(ns = 40, theta = 0, omega = 0.1)
#'
#'## Investigate if the stasis model is an adequate description of the data
#'fit4adequacy.stasis(x)
#'
fit4adequacy.stasis<-function(y, nrep=1000, conf=0.95, plot=TRUE, omega=NULL){

  x<-y$mm
  v<-y$vv
  n<-y$nn
  tt<-y$tt
  
  theta<-opt.joint.Stasis(y)$parameters[1]
  if (is.null(omega)) omega<-opt.joint.Stasis(y)$parameters[2]

  lower<-(1-conf)/2
  upper<-(1+conf)/2


  # Compute the test statistics for the observed time series
  obs.auto.corr<-auto.corr(x, model="stasis")
  obs.runs.test<-runs.test(x, model="stasis", theta=theta)
  obs.slope.test<-slope.test(x,tt, model="stasis", theta=theta)
  obs.net.change.test<-net.change.test(x, model="stasis")

  #Run parametric bootstrap
  out.auto<-auto.corr.test.stasis(y, nrep, conf, plot=FALSE, theta, omega)
  out.runs<-runs.test.stasis(y, nrep, conf, plot=FALSE, theta, omega)
  out.slope<-slope.test.stasis(y, nrep, conf, plot=FALSE, theta, omega)
  out.net<-net.change.test.stasis(y, nrep, conf, plot=FALSE, theta, omega)

  #Preparing the output
  output<-c(as.vector(matrix(unlist(out.auto[[3]]),ncol=5,byrow=FALSE)),
                        as.vector(matrix(unlist(out.runs[[3]]),ncol=5,byrow=FALSE)),
                        as.vector(matrix(unlist(out.slope[[3]]),ncol=5,byrow=FALSE)),
                        as.vector(matrix(unlist(out.net[[3]]),ncol=5,byrow=FALSE)))

  output<-as.data.frame(cbind(c(output[c(1,6,11,16)]), c(output[c(2,7,12,17)]), c(output[c(3,8,13,18)]), c(output[c(4,9,14,19)]), c(output[c(5,10,15,20)])), ncol=5)

  rownames(output)<-c("auto.corr", "runs.test", "slope.test", "net.change.test")
  colnames(output)<-c("estimate", "min.sim" ,"max.sim","p-value", "result")

  if (plot==TRUE) {
    par(mfrow=c(2,2));
    model.names<-c("auto.corr", "runs.test", "slope.test", "net.change.test")
    plotting.distributions(out.auto$replicates,obs.auto.corr, model.names[1], xlab="Simulated data", main="Autocorrelation");
    plotting.distributions(out.runs$replicates,obs.runs.test, model.names[2], xlab="Simulated data", main="Runs");
    plotting.distributions(out.slope$replicates,obs.slope.test, model.names[3], xlab="Simulated data", main="Fixed variance");
    plotting.distributions(out.net$replicates,obs.net.change.test, model.names[4], xlab="Simulated data", main="Net evolution");
  }
  summary.out<-as.data.frame(c(nrep, conf))
  rownames(summary.out)<-c("replications", "confidence level")
  colnames(summary.out)<-("Value")
  out<- list("info" = summary.out, "summary" = output)
  return(out)
}
