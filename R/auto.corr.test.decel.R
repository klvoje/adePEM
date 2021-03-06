#' @title Applying the autocorrelation test to the decelerated evolution model
#'
#' @description Investigates if the decelerated evolution model is an adequate statistical description of an evolutionary
#' time series by applying the autocorrelation test.
#'
#' @param y a paleoTS object
#'
#' @param r parameter describing the decrease in rate of evolution through time. r is restricted to values smaller than zero 
#' (the model reduces to the BM model when r = 0).
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
#' @details This function calculates the autocorrelation in a vector of sample means
#' defined as the correlation of the first n-1 observations with the last n-1. The
#' autocorrelation is calculated directly on the sample means if the evaluated model is stasis.
#' If a different model is evaluated (e.g. random walk or directional trend), the data is
#' detrended prior to the calculation of autocorrelation.
#'
#' @return First part of the output summarizes the number of iterations in the parametric boostrap and the
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
#'@references Voje, K.L. 2018. Assessing adequacy of models of phyletic evolution in the fossil record. \emph{Methods in Ecology and Evoluton}. (in press).
#'@references Voje, K.L., Starrfelt, J., and Liow, L.H. 2018. Model adequacy and microevolutionary explanations for stasis in the fossil record. \emph{The American Naturalist}. 191:509-523.
#'
#'@seealso \code{\link{fit3adequasy.EB}}, \code{\link{auto.corr.test.trend}}, \code{\link{auto.corr.test.stasis}}
#' @export
#'@examples
#'## generate a paleoTS objects by simulating early burst
#'x <- sim.accel_decel(ns=40, r=-1, vs=0.1)
#'
#'## investigate if the time series pass the adequasy test
#'auto.corr.test.decel(x)
#'

auto.corr.test.decel<-function(y, r=NULL, vstep=NULL, nrep=1000, conf=0.95, plot=TRUE, save.replicates=TRUE){

  x<-y$mm
  v<-y$vv
  n<-y$nn
  time<-y$tt

  if (is.null(vstep)) vstep<-opt.joint.decel(y)$parameters[2]
  if (is.null(r)) r<-opt.joint.decel(y)$parameters[3]
  
  lower<-(1-conf)/2
  upper<-(1+conf)/2

  obs.auto.corr<-auto.corr(x, model="accel_decel")

  ### Parametric bootstrap routine ###

  #Matrix that will contain the test statistic for each simuluated data set (time series)
  bootstrap.matrix<-matrix(data = NA, nrow = nrep, ncol = 1)


  # parametric boostrap
  for (i in 1:nrep){

    x.sim<-sim.accel_decel(ns=length(x), r=r, vs=vstep, vp=mean(v), nn=n, tt=time)

    bootstrap.matrix[i,1]<-auto.corr(x.sim$mm, model="accel_decel")

  }

  # Estimating the ratio of how often the observed autocorrelation is smaller than the autocorrelation in the simulated data
  bootstrap.auto.corr<-length(bootstrap.matrix[,1][bootstrap.matrix[,1]>obs.auto.corr])/nrep

  # Calculating the "p-value" and whether the observed data passed the test statistic
  if (bootstrap.auto.corr>round(upper,3) | bootstrap.auto.corr<round(lower,3)) pass.auto.corr.test<-"FAILED" else pass.auto.corr.test<-"PASSED"
  if(bootstrap.auto.corr>0.5) bootstrap.auto.corr<-1-bootstrap.auto.corr

  # Plot the test statistics estimated from the simulated data
  if (plot==TRUE){
    layout(1:1)
    plotting.distributions(bootstrap.matrix[,1],obs.auto.corr, test="auto.corr", xlab="Simulated data", main="Autocorrelation");
  }

  #Preparing the output
  output<-as.data.frame(cbind(round(obs.auto.corr,5), round(min(bootstrap.matrix),5), round(max(bootstrap.matrix),5), bootstrap.auto.corr/0.5, pass.auto.corr.test), nrow=5, byrow=TRUE)
  rownames(output)<-"auto.corr"
  colnames(output)<-c("estimate", "min.sim" ,"max.sim", "p-value", "result")

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
