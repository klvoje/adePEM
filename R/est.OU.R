#' @title Predicting trait values for a given OU model
#'
#' @description 
#' This function calculates the predicted trait values for a given OU model.
#'This function is used by other functions and will generally not be used directly by users.
#'
#' @param x a paleoTS object
#'
#' @param w an object of class paleoTSfit
#'
#' @param tt vector of sample ages, increases from oldest to youngest
#'
#' @export
#' 
#' @return Predicted trait values for an OU model with given parameters.
#' 
#' @author Kjetil L. Voje
#' 
#'

est.OU<-function (x, w, tt) 
{
  ee <- ii <- array(dim = length(tt))
  mn <- w$modelName
  mp <- w$par
  x0 <- ifelse(w$method == "AD", x$mm[1], w$par["anc"])
  ttp <- tt
  
  
  
  ee <- mp["theta"] * (1 - exp(-mp["alpha"] * ttp)) + mp["anc"] * exp(-mp["alpha"] * ttp)
  vv <- (mp["vstep"]/(2 * mp["alpha"])) * (1 - exp(-2 * mp["alpha"] * ttp))
  
  
  res <- list(tt = ttp, ee = ee, ll = ee - 1.96 * sqrt(vv), 
              uu = ee + 1.96 * sqrt(vv))
  return(res)
}