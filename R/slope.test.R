#' @title Estimating the linear time-dependency of the trait variance
#'
#' @description This function estimates the linear dependency of the trait variance on time.
#' This function is used by other functions and will generally not be used directly by users.
#'
#' @param x vector of sample means
#'
#' @param tt vector of sample ages
#'
#' @param model the model being evaluated
#'
#' @param theta if the tested model is stasis, theta is the estimated theta from the data.
#'
#' @details If model = stasis: Estimates the slope of the least square regression of the size of deviations
#' (their absolute value) from the optimal phenotype as a function of time. If model = random walk or directional
#' trend: Estimates the slope of the least square regression of the size of the detrended data as a function of time.
#' This function is used by other functions and will generally not be used directly by users.
#' @export
#' @return least-square slope estimate
#'

slope.test <- function(x, model, tt, theta=NULL, int=NULL, mstep=NULL){
  if (model=="RW")
  {
    x<-x-x[1]
    x<-diff(x,1)
    tt<-tt[-(length(tt))]
    slope.est<-(lm((abs(x))~tt)$coeff[2])
  }
  
  if (model=="EB")
  {
    x<-x-x[1]
    x<-diff(x,1)
    time<-time[-(length(time))]
    slope.est<-(lm(log(abs(x+0.00000001))~time)$coeff[2])
  }
  
  if (model =="trend")
  {
    x<-x-(int+mstep*tt)
    slope.est<-(lm((abs(x))~tt)$coeff[2])
  }

   if (model =="stasis")
  {
    resid_stasis<-abs(x-theta)
    slope.est<-lm(resid_stasis~tt)$coeff[2]
  }

  return(slope.est)
}

