#' @title Applying 3 adequacy tests to the Random walk model
#'
#' @description Investigating if the Random walk model is an adequate statistical description of an evolutionary
#' time series by applying the following tests (1) autocorrelation (2) runs test, and (3) constant variation.
#'
#' @param y a paleoTS object
#'
#' @param nrep number of iterations in the parametric bootstrap (number of simulated time series); default is 1000.
#'
#' @param conf confidence level for judging whether a model is an adequate statistical description of the data.
#' Number must be between 0 and 1. A higher number means less strict judgment of whether a model is adequate; default
#' is 0.95. Tests are two-tailed, which means a model is judged adequate if the observed test statistic is within the 2.5
#' percent of the extreme values of the calculated test statistics on the simulated data given the default confidence
#' value of 0.95.
#' 
#' @param cutoff confidence level for how often the distance traveled in morphospace as a function of time on average should be 
#' larger compared to a linear model describing a constant rate of directional change from the ancestral trait state to the 
#' last population trait mean. Number must be between 0 and 1. Default is 0.80.
#'
#' @param plot logical; if TRUE, the value of the test statistic calculated based on the observed fossil
#' time series is plotted on the distribution of test statistics calculated on the simulated time series;
#' default is TRUE.
#'
#' @param vstep the variance of the step distribution. This parameter is automatically estimated from the data, if not set
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
#'@seealso \code{\link{fit3adequacy.trend}}, \code{\link{fit4adequacy.stasis}}
#' @export
#'@examples
#'## generate a paleoTS objects by simulating random walk
#'x <- sim.GRW(ns=40, ms=0, vs=0.1)
#'
#'## Investigate if the random walk model is an adequate description of the data
#'fit3adequacy.RW(x)
#'

fit3adequacy.OU<-function(y, nrep=1000, conf=0.95, cutoff=0.80, plot=TRUE){

  x<-y$mm
  v<-y$vv
  n<-y$nn
  tt<-y$tt

  anc<-opt.joint.OU(y)$parameters[1]
  vstep<-opt.joint.OU(y)$parameters[2]
  theta<-opt.joint.OU(y)$parameters[3]
  alpha<-opt.joint.OU(y)$parameters[4]
  
  tmp_OU<-opt.joint.OU(y)
  pred_OU<-est.OU(y, tmp_OU, tt=tt)
  detrended_OU<-x-pred_OU$ee
  
  lower<-(1-conf)/2
  upper<-(1+conf)/2

  # Compute the test statistics for the observed time series
  obs.auto.corr<-auto.corr(detrended_OU, model="OU", tt)
  obs.runs.test<-runs.test(detrended_OU, model="OU", tt)
  obs_sum_of_residuals<-0
  
  #Run parametric bootstrap
    out.auto<-auto.corr.test.OU(y, nrep, conf, plot=FALSE, save.replicates = TRUE)
    out.runs<-runs.test.OU(y,nrep, conf, plot=FALSE, save.replicates = TRUE)
    out.var<-variance.test.OU(y,nrep, cutoff, plot=FALSE, save.replicates = TRUE)

  #Preparing the output
output<-c(as.vector(matrix(unlist(out.auto[[3]]),ncol=5,byrow=FALSE)),
          as.vector(matrix(unlist(out.runs[[3]]),ncol=5,byrow=FALSE)),
          as.vector(matrix(unlist(out.var[[3]]),ncol=5,byrow=FALSE)))

  output<-as.data.frame(cbind(c(output[c(1,6,11)]), c(output[c(2,7,12)]),
                              c(output[c(3,8,13)]), c(output[c(4,9,14)]),
                              c(output[c(5,10,15)])), ncol=5)

  rownames(output)<-c("auto.corr", "runs.test", "slope.test")
  colnames(output)<-c("estimate", "min.sim" ,"max.sim","p-value", "result")

  if (plot==TRUE) {
    par(mfrow=c(1,3))
    model.names<-c("auto.corr", "runs.test", "slope.test")
    plotting.distributions(out.auto$replicates,obs.auto.corr, model.names[1], xlab="Simulated data", main="Autocorrelation");
    plotting.distributions(out.runs$replicates,obs.runs.test, model.names[2], xlab="Simulated data", main="Runs");
    plotting.distributions(out.var$replicates,obs_sum_of_residuals, model.names[3], xlab="Simulated data", main="Initial rapid change");

  }
  summary.out<-as.data.frame(c(nrep, conf, cutoff))
  rownames(summary.out)<-c("replications", "confidence level", "cut-off faster evolution")
  colnames(summary.out)<-("Value")
  out<- list("info" = summary.out, "summary" = output)
  return(out)
}
