#' @title Estimating the linear time-dependency of the trait variance
#'
#' @description This function estimates the linear dependency of the trait variance on time.
#' This function is used by other functions and will generally not be used directly by users.
#'
#' @param x vector of sample means
#'
#' @param time vector of sample ages
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

slope.test <- function(x, time, model, theta=NULL, anc=NULL, mstep=NULL){
  if (model=="RW")
  {
    x<-x-x[1]
    x<-diff(x,1)
    time<-time[-(length(time))]
    slope.est<-(lm((abs(x))~time)$coeff[2])
  }

  if (model =="trend")
  {
    x<-x-(anc+mstep*time)
    slope.est<-(lm((abs(x))~time)$coeff[2])
  }

   if (model =="stasis")
  {
    resid_stasis<-abs(x-theta)
    slope.est<-lm(resid_stasis~time)$coeff[2]
  }

  return(slope.est)
}

