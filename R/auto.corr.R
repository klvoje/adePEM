#' @title Calculating autocorrelation
#'
#' @description Estimates the autocorrelation in the data.
#' This function is used by other functions and will generally not be used directly by users.
#'
#' @param x vector of sample means
#'
#' @param model the model being evaluated
#'
#' @details This function calculates the autocorrelation in a vector of sample means
#' defined as the correlation of the first n-1 observations with the last n-1. The
#' autocorrelation is calculated directly on the sample means if the evaluated model is stasis.
#' If a different model is evaluated (random walk or directional trend), the data is
#' detrended prior to the calculation of autocorrelation.
#' @export
#'
#' @return the autocorrelation
#'
#'


auto.corr <- function(x, model, tt=NULL, int=NULL, mstep=NULL){
  if (model=="RW" | model=="EB")
    {
    x<-x-x[1]
    x<-diff(x,1)
    x<-c(0,x)
  }
  if (model=="trend")
    {
    x<-x-(int+mstep*tt)
  }
  return(cor(x[1:length(x)-1],x[2:length(x)]))
}
