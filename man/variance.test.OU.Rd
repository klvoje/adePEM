% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/variance.test.OU.R
\name{variance.test.OU}
\alias{variance.test.OU}
\title{Applying the constant variance test to the OU evolution model}
\usage{
variance.test.OU(y, nrep = 1000, conf = 0.95, plot = TRUE,
  save.replicates = TRUE)
}
\arguments{
\item{y}{a paleoTS object}

\item{nrep}{number of iterations in the parametric bootstrap (number of simulated time series); default is 1000.}

\item{conf}{confidence level for judging whether a model is an adequate statistical description of the data.
Number must be between 0 and 1. A higher number means less strict judgment of whether a model is adequate; default
is 0.95. Tests are two-tailed (except for the net evolution test), which means a model is judged adequate if the observed test statistic is within the 2.5
percent of the extreme values of the calculated test statistics on the simulated data given the default confidence
value of 0.95.}

\item{plot}{logical; if TRUE, the value of the test statistic calculated based on the observed fossil
time series is plotted on the distribution of test statistics calculated on the simulated time series;
default is TRUE.}

\item{save.replicates}{logical; if TRUE, the values of the test statistic calculated on the simulated time
series is saved and can be accessed later for plotting purposes; default is TRUE.}
}
\value{
First part of the output summarizes the number of iterations in the parametric bootstrap and the
confidence level for judging whether a model is an adequate statistical description of the data. The last
part of the output is:

\item{estimate}{The calculated test statistic on the observed data.}
 \item{min.sim}{The smallest test statistic calculated on the simulated data.}
 \item{max.sim}{The largest test statistic calculated on the simulated data.}
 \item{p-value}{Not a real p-value, but is calculated as the fraction of simulated test statistics
 that is larger (or smaller) than the calculated test statistic on the observed data divided by 0.5.
 A value of 1 means 50 percent of the test statistics on the simulated data are larger and smaller
 than the calculated statistic on the observed data. A value of 0.10 means 90 percent of the test
 statistics on the simulated data are larger or smaller than the test statistic on the observed time
 series.}
 \item{result}{Whether the model PASSED or FAILED the adequacy test. The outcome depends on the
 confidence level.}
}
\description{
Investigates if the OU is an adequate statistical description of an evolutionary
time series by applying the constant variance test.
}
\details{
Estimates the slope of the least square regression of the size of the detrended data (their absolute value) from the average
as a function of time.as a function of time.
}
\examples{
## generate a paleoTS objects by simulating early burst
x <- sim.accel_decel(ns=20)

## investigate if the time series pass the adequacy test
variance.test.decel(x)

}
\references{
Voje, K.L. 2018. Assessing adequacy of models of phyletic evolution in the fossil record. \emph{Methods in Ecology and Evoluton}. (in press).

Voje, K.L., Starrfelt, J., and Liow, L.H. 2018. Model adequacy and microevolutionary explanations for stasis in the fossil record. \emph{The American Naturalist}. 191:509-523.
}
\seealso{
\code{\link{fit3adequasy.RW}}, \code{\link{slope.test.stasis}}, \code{\link{slope.test.trend}}
}
\author{
Kjetil L. Voje
}
