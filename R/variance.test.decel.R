#' @title Applying the constant variance test to the decelerated evolution model
#'
#' @description Investigates if the decelerated evolution model is an adequate statistical description of an evolutionary
#' time series by applying the constant variance test.
#'
#' @param y a paleoTS object
#' 
#' @param r parameter describing the decreasing rate of evolution through time. r is restricted to values smaller than zero 
#' (the model reduces to the BM model when r = 0).
#'
#' @param vstep the initital variance of the step distribution estimated from the observed data.
#'
#' @param nrep number of iterations in the parametric bootstrap (number of simulated time series); default is 1000.
#'
#' @param cutoff confidence level for judging whether a model is an adequate statistical description of the data.
#' Number must be between 0 and 1. Default is 0.80.
#'
#' @param plot logical; if TRUE, the value of the test statistic calculated based on the observed fossil
#' time series is plotted on the distribution of test statistics calculated on the simulated time series;
#' default is TRUE.
#'
#' @param save.replicates logical; if TRUE, the values of the test statistic calculated on the simulated time
#' series is saved and can be accessed later for plotting purposes; default is TRUE.
#'
#' @details Tests if the distances traveled in morphospace as a function of time on average will be larger 
#' compared to a linear model describing a constant rate of directional change from the ancestral trait state 
#' to the last population trait mean.
#'
#' @return First part of the output summarizes the number of iterations in the parametric bootstrap and the
#' confidence level for judging whether a model is an adequate statistical description of the data. The last
#' part of the output is:
#'
#' @return
#'  \item{estimate}{Set to 0.}
#'  \item{min.sim}{The smallest test statistic calculated on the simulated data.}
#'  \item{max.sim}{The largest test statistic calculated on the simulated data.}
#'  \item{p-value}{Not a real p-value, but is calculated as the fraction of simulated test statistics
#'  that is larger than 0. A value of 0.90 means 90 percent of the test
#'  statistics on the simulated data are larger than 0.}
#'  \item{result}{Whether the model PASSED or FAILED the adequacy test. The outcome depends on the
#'  confidence level.}
#'
#'@author Kjetil L. Voje
#'
#'@seealso \code{\link{fit3adequasy.decel}}, \code{\link{slope.test.stasis}}, \code{\link{slope.test.trend}}
#'@export
#'@examples
#'## generate a paleoTS objects by simulating early burst
#'x <- sim.accel_decel(ns=20)
#'
#'## investigate if the time series pass the adequacy test
#'variance.test.decel(x)
#'


variance.test.decel<-function(y, r=NULL, vstep=NULL, nrep=1000, cutoff=0.80, plot=TRUE, save.replicates=TRUE){

  x<-y$mm
  v<-y$vv
  n<-y$nn
  time<-y$tt
  
  if (is.null(vstep)) vstep<-opt.joint.decel(y)$parameters[2]
  if (is.null(r)) r<-opt.joint.decel(y)$parameters[3]

  ### Parametric bootstrap routine ###

  #Matrix that will contain the test statistic for each simulated data set (time series)
  bootstrap.matrix<-matrix(data = NA, nrow = nrep, ncol = 1)
  
  obs_sum_of_residuals<-0
  # parametric boostrap
  for (i in 1:nrep){

    x.sim<-sim.accel_decel(ns=length(x), r=r, vs=vstep, vp=mean(v), nn=n, tt=time)
    dist_trav_morphospace.sim<-dist.in.morphospace(x.sim, correct= FALSE,iter = 10000)$observed.accumulated.change.not.bias.cor
    slope_linear_model_sim<-max(dist_trav_morphospace.sim)/max(x.sim$tt)
    bootstrap.matrix[i,1]<-sum(c(0,dist_trav_morphospace.sim)-(slope_linear_model_sim*x.sim$tt))
    
  }

  # Estimating the ratio of how often the observed slope statistic is smaller than the slope tests in the simulated data
  bootstrap.var.test<-length(bootstrap.matrix[,1][bootstrap.matrix[,1]>obs_sum_of_residuals])/nrep

  if (bootstrap.var.test<round(cutoff,3)) pass.var.test<-"FAILED" else pass.var.test<-"PASSED"
  
  
  
  # Plot the test statistics estimated from the simulated data
  if (plot==TRUE) {
    layout(1:1)
    plotting.distributions(bootstrap.matrix[,1],obs_sum_of_residuals, test="slope.test", xlab="Simulated data", main="Initial rapid evolution");
  }

  #Preparing the output
  output<-as.data.frame(cbind(round(obs_sum_of_residuals,5), round(min(bootstrap.matrix),5), round(max(bootstrap.matrix),5), bootstrap.var.test, pass.var.test), nrow=5, byrow=TRUE)
  rownames(output)<-"var.test"
  colnames(output)<-c("cut-off","min.sim" ,"max.sim","p-value", "result")

  summary.out<-as.data.frame(c(nrep, cutoff))
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

