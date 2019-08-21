#' @title Calculating distances in morphospace.
#'
#' @description 
#' This function calculates the distances in morphospace between sample means in a time series.
#'This function is used by other functions and will generally not be used directly by users.
#'
#' @param y a paleoTS object
#'
#' @param correct logical, if TRUE, distances between population means are corrected for sampling variance.  
#'
#' @param iter number of iterations in the parametric bootstrap; default is 1000.
#'
#' @export
#' 
#' @return distances in morphospace between sample means in time series.
#' 
#' @author Kjetil L. Voje
#' 
#' @references Voje, K.L. 2016. Tempo does not correlate with mode in the fossil record. \emph{Evoluton} 70:2678–2689.
#'


### Function calculating distances traveled between sample means in a time series.

# mm = vector of sample means.
# vv = vector of sample variances.
# nn = vector of sample sizes.
# tt = vector of sample ages.
# correct = If TRUE, correcting for estimation error in sample means when estimating distances travelled in phenotype space.
# iter = number of iterations in bootstrap routine when correcting for estimation error in sample means.

dist.in.morphospace<-function (y, correct= TRUE, iter = 10000)
{
  
  mm<-y$mm
  vv<-y$vv
  nn<-y$nn
  tt<-y$tt
  
  Number.of.samples<-length(mm)
  
  obs.abs.change<-rep(NA, Number.of.samples-1)
  for (i in 1:Number.of.samples)
      {
    obs.abs.change[i]<-abs(mm[i]- mm[i+1])
      }

  
  obs.abs.change<-obs.abs.change[1:(Number.of.samples-1)]
  if (correct == TRUE) 
  {
  trait<-matrix(NA, ncol=Number.of.samples, nrow=iter) 
    
    for (j in 1:iter)
    for (i in 1:(length(mm)))
      {
        trait[j,i]<- mean(rnorm(nn[i], mean = 0, sd = sqrt(vv[i])))
      }

  abs.sim.data<-matrix(NA, ncol=length(obs.abs.change), nrow=iter) 
  for (j in 1:iter)
    for (i in 1:length(obs.abs.change))
      {
        abs.sim.data[j,i]<-abs(trait[j,i]- trait[j,i+1])
      }
  
  data.corrected<-matrix(NA, ncol=length(mm)-1, nrow=iter) 
  for (i in 1:iter)
    {
      data.corrected[i,]<-obs.abs.change-abs.sim.data[i,]
    }
  
  
  obs.abs.change.sim<-rep(NA,length(obs.abs.change)) 
  obs.sd.sim<-rep(NA,length(obs.abs.change)) 
  for (i in 1:length(data.corrected[1,]))
    {
      obs.abs.change.sim[i]<-mean(data.corrected[,i])
      obs.sd.sim[i]<-sd(data.corrected[,i])
    }

  
  for (i in 1:(length(obs.abs.change.sim)))
    {
      if (obs.abs.change.sim[i]<0) obs.abs.change.sim[i]<-0
    }
  
  
  accu.evo.after.correction<-rep(NA, length(obs.abs.change))
  for (i in 1:length(obs.abs.change))
    {
    accu.evo.after.correction[i]<-sum(obs.abs.change.sim[1:i])
    }
  accu.evo.after.correction<-c(0,accu.evo.after.correction)
  obs.sd.sim<-c(0, obs.sd.sim)
  
  
  }


  if (correct == TRUE) output<-list(length.of.time.series = Number.of.samples, observed.absolute.change.not.bias.cor = obs.abs.change, observed.accumulated.change.not.bias.cor=cumsum(obs.abs.change),  net.evolution = as.numeric(abs(tail(mm, n=1)-head(mm, n=1))), observed.absolute.change.bias.cor=obs.abs.change.sim, observed.accumulated.change.bias.cor = accu.evo.after.correction,SE.observed.accumulated.change.bias.cor = obs.sd.sim, average.distance.traveled.between.samples = sum(obs.abs.change.sim)/length(obs.abs.change))
  else output<-list(length.of.time.series = Number.of.samples, observed.absolute.change.not.bias.cor = obs.abs.change, observed.accumulated.change.not.bias.cor=cumsum(obs.abs.change), net.evolution = as.numeric(abs(tail(mm, n=1)-head(mm, n=1))))

  return(output)
  
  }