#' @title Calculating the number of runs
#'
#' @description This function calculates the number of runs of successive positive or
#' negative deviations from theta (if model = stasis) or the mean of the data.
#' This function is used by other functions and will generally not be used directly by users.
#'
#' @param x vector of sample means
#'
#' @param model the model being evaluated
#'
#' @param theta the value of theta if model = stasis
#'
#' @details The runs test is applied to the sign of the residuals (i.e. θ – trait value)
#' to identify series that have non-random patterns in the sign of deviations. For a time
#' series of length n, the number of runs (one run is a sequence of consecutive numbers with
#' same sign), is approximately normal with mean μ=(2(n_+*n_-))/n+1 and variance (μ-1)(μ-2)/(n-1),
#' where n+ and n- are the number of residuals above and below the optimum respectively.
#' @export
#' @return net evolution
#'



runs.test <- function(x, model, tt=NULL, theta=NULL, int=NULL, mstep=NULL){
  if (model=="RW")
  {
    x<-x-x[1]
    x<-diff(x,1)
    x<-c(0,x)

    mu = 2*(sum(x>mean(x))*sum(x<mean(x)))/length(x) +1;
    # with variance
    vr = (mu-1)*(mu-2)/(length(x)-1);
    z = (sum(diff(sign(x-mean(x)))!=0)+1 - mu)/sqrt(vr)
  }
  if (model=="trend")
  {
    x<-x-(int+mstep*tt)
    
    mu = 2*(sum(x>mean(x))*sum(x<mean(x)))/length(x) +1;
    # with variance
    vr = (mu-1)*(mu-2)/(length(x)-1);
    z = (sum(diff(sign(x-mean(x)))!=0)+1 - mu)/sqrt(vr)
  }

  if (model=="stasis")
  {
  mu = 2*(sum(x>theta)*sum(x<theta))/length(x) +1;
  # with variance
  vr = (mu-1)*(mu-2)/(length(x)-1);
  z = (sum(diff(sign(x-theta))!=0)+1 - mu)/sqrt(vr)
  }

  return(z)
}
