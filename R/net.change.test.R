#' @title Calculating net evolution
#'
#' @description Calculates the net evolution in the data, i.e. the absolute difference
#' in trait value between the first and last sample mean in the time series.
#' This function is used by other functions and will generally not be used directly by users.
#'
#' @param x vector of sample means
#'
#' @param model the model being evaluated
#'
#' @details This function calculates the net evolution in a vector of sample means
#' defined as the the absolute difference in trait value between the first and last
#' sample mean in the time series.
#' @export
#' @return net evolution
#'
#'

net.change.test<-function(x, model){
  if(model=="RW" | model=="trend" | model=="stasis"){
    net.change<-abs(head(x,1)-tail(x,1))
    return(net.change)
  }

}

