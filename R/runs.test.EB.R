#' @title Applying the runs test to the Early burst model
#'
#' @description Investigates if the Early burst model is an adequate statistical description of an evolutionary
#' time series by applying the runs test.
#'
#' @param y a paleoTS object
#' 
#' @param alpha parameter describing the decreasing rate change through time. alpha is restricted to values below zero 
#' (the model reduces to the BM model when alpha = 0).
#'
#' @param vstep the variance of the step distribution estimated from the observed data.
#'
#' @param nrep number of iterations in the parametric bootstrap (number of simulated time series); default is 1000.
#'
#' @param conf confidence level for judging whether a model is an adequate statistical description of the data.
#' Number must be between 0 and 1. A higher number means less strict judgment of whether a model is adequate; default
#' is 0.95. Tests are two-tailed, which means a model is judged adequate if the observed test statistic is within the 2.5
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
#' @details This function applies a runs test in order to investigate if the random walk model can be judged an
#' adequate statistical description of the data. After detrending, there should be no tendency in the data to successively deviate
#' from the average in the same direction and the runs test is applied to the sign of the residuals (i.e. θ – trait value)
#' to identify series that have non-random patterns in the sign of deviations. For a time series of length n,
#' the number of runs (one run is a sequence of consecutive numbers with same sign), is approximately normal
#' with mean μ=(2(n_+ n_-))/n+1 and variance (μ-1)(μ-2)/(n-1), where n+ and n- are the number of residuals
#' above and below the optimum respectively. The mean and variance are used to calculate the standard/Z-score
#' implemented as the test statistic.
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
#'@seealso \code{\link{runs.test.stasis}}, \code{\link{runs.test.RW}}, \code{\link{fit3adequasy.trend}}
#' @export
#'@examples
#'## generate a paleoTS objects by simulating early burst
#'x <- sim.EB(ns=40, alpha=-1, vs=0.1)
#'
#'## investigate if the time series pass the adequacy test
#'runs.test.EB(x,vstep)
#'

runs.test.EB<-function(y, alpha, vstep, nrep=1000, conf=0.95, plot=TRUE, save.replicates=TRUE){

  x<-y$mm
  v<-y$vv
  n<-y$nn
  time<-y$tt

  lower<-(1-conf)/2
  upper<-(1+conf)/2

  obs.runs.test<-runs.test.EB(x, model="EB")

  ### Parametric bootstrap routine ###

  #Matrix that will contain the test statistic for each simulated data set (time series)
  bootstrap.matrix<-matrix(data = NA, nrow = nrep, ncol = 1)

  # parametric boostrap
  for (i in 1:nrep){

    x.sim<-sim.EB(ns=length(x), alpha=alpha, vs=vstep, vp=mean(v), nn=n, tt=time)

    bootstrap.matrix[i,1]<-runs.test.EB(x.sim$mm, model="EB")

  }

  # Estimating the ratio of how often the observed runs test is smaller than the runs tests in the simulated data
  bootstrap.runs.test<-length(bootstrap.matrix[,1][bootstrap.matrix[,1]>obs.runs.test])/nrep

  # Calculating the "p-value" and whether the observed data passed the test statistic
  if (bootstrap.runs.test>round(upper,3) | bootstrap.runs.test<round(lower,3)) pass.runs.test<-"FAILED" else pass.runs.test<-"PASSED"
  if(bootstrap.runs.test>0.5) bootstrap.runs.test<-1-bootstrap.runs.test

  # Plot the test statistics estimated from the simulated data
  if (plot==TRUE) {
    layout(1:1)
    plot.distributions(bootstrap.matrix[,1],obs.runs.test, test="runs.test", xlab="Simulated data", main="Runs");
  }

  #Preparing the output
  output<-as.data.frame(cbind(round(obs.runs.test,5), round(min(bootstrap.matrix),5), round(max(bootstrap.matrix),5), bootstrap.runs.test/0.5, pass.runs.test), nrow=5, byrow=TRUE)
  rownames(output)<-"runs.test"
  colnames(output)<-c("estimate","min.sim" ,"max.sim","p-value", "result")

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
