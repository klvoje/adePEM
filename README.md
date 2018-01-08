## About adePEM (Assessing adequasy of phyletic evolution models) 

This package allows assessing the adequasy of the three canonical models of evolution in the fossil record, stasis, directional trends and random walk developed by Hunt (2006) and implemented in the R package `paleoTS`.  

Model fit in `paleoTS` is evaluated using AICc. However, the best model among a list of candidates according to an information criterion may not describe the data particularly well. This is true because any set of candidate models will only reflect a subset of ways of portraying evolutionary dynamics in a lineage. 

Passing adequacy tests suggests the model provides an adequate statistical description of the trait dynamics observed in the data and that meaningful inferences can be drawn from the fitted model parameters. 

The package includes functions to simulate datasets, calculate summary statistics and plot results. 

The methods for assessing adequasy of the stasis model were first described in the paper: Model adequacy and microevolutionary explanations for stasis in the fossil record. Voje, K.L., Starrfelt, J., and Liow, L.H. The American Naturalist. [In press] (http://www.amnat.org/an/newpapers/AprVoje-A.html).

The manuscript where the `adePEM` package is presented is currently under review.


## Installation

Install the package from github using devtools:

```
install.packages("devtools")
library("devtools")

devtools::install_github("klvoje/adePEM")
```

The package depends on the `paleoTS` and `pracma` libraries, which should load automatically when installing `adePEM` from github.


## Example

We are interested in analyzing the evolution of element length (measured in mm) in the conodont Pterospathodus. The data is available as part of the `adePEM` package and was originally published by Jones (2009). The data (`element.length`) is already a `paleoTS` object. We first plot the data. 

```
plot.paleoTS(element.length)
```
![time seires](https://github.com/klvoje/adePEM/blob/master/extra/time.series.png)

Time (the x-axis) is in millions of years and the trait is measured in millimeters. Error bars represent one standard error.

We then run the `fit3models` function from the `paleoTS` package to check the relative fit of the stasis, random walk and directional trend models to the data.
```
fit3models(element.length, pool=TRUE)

Comparing 3 models [n = 31, method = Joint]

           logL K      AICc Akaike.wt
GRW    25.38445 3 -43.88002     0.262
URW    25.12370 2 -45.81882     0.690
Stasis 22.47400 2 -40.51943     0.049
```

The random walk (URW) model has the best fit to the data according to the AICc scores. However, the difference in the AICc score is small (<2 units) relative to the directional trend model (GRW). 

Let's investigate if the random walk represents an adequate statistical description of the trait dynamics in the data. To do that, we run the function `fit3adequasy.RW` from the `adePEM` package. This is a wrapper function that runs 3 adequasy tests at the same time. 

Before we run the adequasy tests, we need to estimate the step variance of the random walk model from the real data. This is done using the `mle.URW` function from the `PaleoTS` package. The estimated step variance is used to simulate a large number of time series where the trait evolves according to a random walk. Test statistics calculated on these simulated data will then be compared to the test statistics calculated on the real data.  
```
# Estimate the vstep parameter from the data:
vstep<-mle.URW(element.length)[1]

# Run adequasy test for the random walk model:
fit3adequasy.RW(element.length, vstep=vstep)


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

The first part of the output summarizes the number of bootstrap replications (the number of simulated data sets) used for assessing model adequasy and the confidence level. 1000 replications and a confidence level of 0.95 are the default values.

The second part of the output prints information on the adequasy tests. The first column contains the adequasy tests and the second column gives the test statistic computed on the real data. From the second column, we see that the autocorrelation is calculated to be negative and the slope test (which is the least-square slope estimate of how deviations from the mean trait value changes with time) is positive. 

The third and fourth columns reports the smallest and largest test statistics calculated on the simulated data sets. As we can see, all the three test statistics computed on the real data (second column) are not close to the extreme values reported in columns three and four. 

The fifth column is not a real p-value, but is calculated as the fraction of simulated test statistics that is larger (or smaller) than the calculated test statistic on the observed data divided by 0.5. A value of 1 means 50 percent of the test statistics on the simulated data are larger and smaller than the calculated statistic on the observed data, respectively. A value of 0.10 means 90 percent of the test statistics on the simulated data are larger (or smaller) than the test statistic on the observed time series. 

The sixth column indicates whether our model passed the adequasy tests. Since we set our confidence level to 0.95 and all values in the fifth column is larger than 0.05, this means the random walk passed all tests for our data set. 

That the random walk model passed all tests can also be seen in the visual representation of the distributions of test statistics, where the test statistics computed for the real data is indicated with a broken (red) line. These plots are generated automatically if `plot = TRUE` (which is the default setting) when we run the `fit3adequasy.RW` function.   

![RW distributions](https://github.com/klvoje/adePEM/blob/master/extra/adequasy.bm.png)


To summarize: Among the three candidate models stasis, random walk and directional change, random walk has the best relative model fit to the data based on AICc. However, a relative better fit for a model (in this case, the random walk model) to a phyletic fossil time series is no guarantee that the model represents a sufficiently good statistical explanation for the trait dynamics observed in the data. We therefore assessed to what extent the random walk model also fitted the data in an absolute sense by running adequasy models. The random walk model passed all adequasy tests, which suggest the random walk model represents an adequate statistical description of the phyletic time series.

If we take a look at the plot of how the trait changes over 6 million years, it seems to suggest a trend towards becoming bigger. Therefore, let's assess the adequasy for the directional trend model on the data. This model did indeed show a quite similar fit to the data based on their AICc scores. 
     
Again, the first thing we need to do is to estimate the model parameters for the real data. The directional change model has two parameters: the `mstep` is the mean of the step distribution while the `vstep` is the variance of the step distribution. We estimate these using the `mle.GRW` function in the paleoTS package. We then run the wrapper function to run all three adequasy tests simultaneously.
```
# Estimate the mstep and vstep parameter from the data:
mstep<-mle.GRW(element.length)[1]
vstep<-mle.GRW(element.length)[2]

# Run adequasy test for the directional trend model:
fit3adequasy.trend(element.length, mstep=mstep, vstep=vstep)

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
The directional trend model fails both the autocorrelation test and the runs test, and passes the slope test. This suggests that the directional trend model is not an adequate statistical description of the data.

![trend distributions](https://github.com/klvoje/adePEM/blob/master/extra/adequasy.trend.png)

Functions for running each test alone are provided in the package (e.g. `auto.corr.test.stasis`, `runs.test.RW`, `slope.test.trend`). The wrapper function for investigating the adequasy for the stasis model is `fit4adequasy.stasis`. One of the tests run by this function is only implemented for stasis (test of to large levers of net evolution), as low amounts of net evolution is part of the verbal definition of stasis, but not for random walk and directional trend.


## Author

Kjetil L. Voje <k.l.voje@gmail.com>


## References

Hunt, G. 2006. Fitting and Comparing Models of Phyletic Evolution: Random Walks and beyond. Paleobiology 32:578–601. 

Jones, D. 2009. Directional evolution in the conodont Pterospathodus. Paleobiology 35: 413-431.

Voje, K.L., Starrfelt, J., and Liow, L.H. 2018. Model adequacy and microevolutionary explanations for stasis in the fossil record. The American Naturalist [in press].

