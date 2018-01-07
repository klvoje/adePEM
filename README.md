
This package allows assessing the adequasy of the three canonical models of evolution in the fossil record, stasis, directional trends and random walk developed by Hunt (2006) and implemented in the R package paleoTS.  
Passing the adequacy tests suggests the model provides a good statistical explanation for the trait dynamics observed in the data and that meaningful inferences can be drawn from the fitted model parameters. 
The package includes functions to simulate datasets, calculate summary statistics and plot results. 
The method for assessing the adequasy of the stasis model was first described in the paper: Model adequacy and microevolutionary explanations for stasis in the fossil record. Voje, K.L., Starrfelt, J., and Liow, L.H. The American Naturalist. In press. The manuscript where this r package is presented is currently under review.


## Installation

Install the package from github using devtools:

```
install.packages("devtools")
library("devtools")

devtools::install_github("klvoje/adePEM")
```

The package depends on the `paleoTS` and `pracma` libraries, which should load automatically when installing adePEM from github.

## Example

We are interested in analyzing the evolution of element length in the conodont Pterospathodus. The data is available as part of the adePEM package and was originally published by Jones (2009).
The data (Jones.2009.element.length) is already a paleoTS object. We first plot the data and run the fot3modells function from the PaleoTS package to check the relatice fit of the stasis, random walk and directional trend model to the data.

```
plot(element.length)
```
![adequate.DT](https://github.com/klvoje/adePEM/blob/master/extra/phenetic.evolution.element.length.png)

```
fit3models(element.length, pool=TRUE)

Comparing 3 models [n = 31, method = Joint]

           logL K      AICc Akaike.wt
GRW    25.38445 3 -43.88002     0.262
URW    25.12370 2 -45.81882     0.690
Stasis 22.47400 2 -40.51943     0.049
```

The plot shows a trend in the data, but it is the random walk (URW) model that has the best fit to the data according to AICc. However, the difference is small (<2 AICc units) relative to he direcional trend model (GRW). 

Let's investigate if the random walk represents an adequate statistical description of the trait dynamics in te data. To do that, we run the function fit3adequasy.BM from the adePEM package. This is a wrapper function that runs 3 statistical tests at the same time and prints the output.  

The first thing we need to do is to estimate the step variance of the random walk model from the real data. This 'observed' step variance is needed to simulate a large number of time series where the trait evolves according to a random walk. Test statistics calculated on these simulate data will then be compared to the test statistics calculated on the real data.   
```
# Estimate the vstep parameter from the data:
vstep<-mle.URW(element.length)[1]

# Run adequasy test for the random walk model:
fit3adequasy.BM(element.length, vstep=vstep)


$info
                   Value
replications     1000.00
confidense level    0.95

$summary
           estimate  min.sim max.sim p-value Result
auto.corr    -0.318 -0.69919 0.38265   0.436 PASSED
runs.test   1.09003 -2.73719 3.34277   0.526 PASSED
slope.test  0.01199  -0.0217 0.04293   0.618 PASSED
```

![adequate.DT](https://github.com/klvoje/adePEM/blob/master/extra/adequasy.BM.png)


```
# Estimate the mstep and vstep parameter from the data:
mstep<-mle.GRW(element.length)[1]
vstep<-mle.GRW(element.length)[2]

# Run adequasy test for the directional trend model:
fit3adequasy.BM(element.length, vstep=vstep)

$info
                   Value
replications     1000.00
confidense level    0.95

$summary
           estimate  min.sim max.sim p-value Result
auto.corr   0.15903 -0.07639 0.93098   0.034 FAILED
runs.test  -0.41179 -4.93166 1.71372   0.042 FAILED
slope.test  0.00355 -0.07232 0.06886   0.862 PASSED
```

![adequate.DT](https://github.com/klvoje/adePEM/blob/master/extra/adequasy.trend.png)
