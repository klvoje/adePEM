#' @title Applying the autocorrelation test to the stasis model
#'
#' @description Investigates if the stasis model is an adequate statistical description of an evolutionary
#' time series by applying the autocorrelation test.
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
#' @param save.replicates logical; if TRUE, the values of the test statistic calculated on the simulated time
#' series is saved and can be accessed later for plotting purposes; default is TRUE.
#' 
#' @param omega evolutonary variance estimated from the observed data. This parameter is automatically estimated from the data, if not set 
#' by the user (usually not recommended).
#'
#' @details This function calculates the autocorrelation in a vector of sample means
#' defined as the correlation of the first n-1 observations with the last n-1. The
#' autocorrelation is calculated directly on the sample means if the evaluated model is stasis.
#' If a different model is evaluated (random walk or directional trend), the data is
#' detrended prior to the calculation of autocorrelation.
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
#'@references Voje, K.L. 2018. Assessing adequacy of models of phyletic evolution in the fossil record. \emph{Methods in Ecology and Evoluton}. (in press).
#'@references Voje, K.L., Starrfelt, J., and Liow, L.H. 2018. Model adequacy and microevolutionary explanations for stasis in the fossil record. \emph{The American Naturalist}. 191:509-523.
#'
#'@seealso \code{\link{fit4adequacy.stasis}}, \code{\link{auto.corr.test.RW}}, \code{\link{auto.corr.test.trend}}
#' @export
#'@examples
#'## generate a paleoTS objects by simulating a stasis time series
#'x <- sim.Stasis(ns = 40, theta = 0, omega = 0.1)
#'
#'## investigate if the time series pass the adequacy test
#'auto.corr.test.stasis(x)
#'

auto.corr.test.stasis<-function(y, nrep=1000, conf=0.95, plot=TRUE, save.replicates=TRUE, omega=NULL){

  x<-y$mm
  v<-y$vv
  n<-y$nn
  tt<-y$tt
  
  lower<-(1-conf)/2
  upper<-(1+conf)/2

  theta<-opt.joint.Stasis(y)$parameters[1]
  if (is.null(omega)) omega<-opt.joint.Stasis(y)$parameters[2]
  
  obs.auto.corr<-auto.corr(x, model="stasis")

  ### Parametric bootstrap routine ###

  #Matrix that will contain the test statistic for each simuluated data set (time series)
  bootstrap.matrix<-matrix(data = NA, nrow = nrep, ncol = 1)


  # parametric boostrap
  for (i in 1:nrep){

    x.sim<-sim.Stasis(ns = length(x), theta = theta, omega = omega, vp = mean(v), nn = n, tt = tt)

    bootstrap.matrix[i,1]<-auto.corr(x.sim$mm, model="stasis")

  }

  # Estimating the ratio of how often the observed autocorrelation is smaller than the autocorrelation in the simulated data
  bootstrap.auto.corr<-length(bootstrap.matrix[,1][bootstrap.matrix[,1]>obs.auto.corr])/nrep

  # Calculating the "p-value" and whether the observed data passed the test statistic
  if (bootstrap.auto.corr>round(upper,3) | bootstrap.auto.corr<round(lower,3)) pass.auto.corr.test<-"FAILED" else pass.auto.corr.test<-"PASSED"
  if(bootstrap.auto.corr>0.5) bootstrap.auto.corr<-1-bootstrap.auto.corr

  # Plot the test statistics estimated from the simulated data
  if (plot==TRUE){
    par(mfrow=c(1,1))
    plotting.distributions(bootstrap.matrix[,1],obs.auto.corr, test="auto.corr", xlab="Simulated data", main="Autocorrelation");
  }

  #Preparing the outout
  output<-as.data.frame(cbind(round(obs.auto.corr,5), round(min(bootstrap.matrix),5), round(max(bootstrap.matrix),5), bootstrap.auto.corr/0.5, pass.auto.corr.test), nrow=5, byrow=TRUE)
  rownames(output)<-"auto.corr"
  colnames(output)<-c("estimate", "min.sim" ,"max.sim", "'p-value'", "result")

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
