## About adePEM (Assessing adequacy of phyletic-evolution models) 

This package allows assessing the adequacy of the three canonical models of evolution in the fossil record stasis, trend and random walk developed by Hunt (2006) and implemented in the R package `paleoTS`.  

Model fit in `paleoTS` is evaluated using AICc. However, the best model among a list of candidates according to an information criterion may not describe the data particularly well. This is true because any set of candidate models will only reflect a subset of ways of portraying evolutionary dynamics in a lineage. 

Passing adequacy tests suggests the model provides an adequate statistical description of the trait dynamics observed in the data and that meaningful inferences can be drawn from the fitted model parameters. 

The package includes functions to simulate datasets, calculate summary statistics and plot results. 

The methods for assessing adequacy of the stasis model were first described in the paper: Model adequacy and microevolutionary explanations for stasis in the fossil record. Voje, K.L., Starrfelt, J., and Liow, L.H. The American Naturalist. [In press] (http://www.journals.uchicago.edu/doi/10.1086/696265).

The manuscript where the `adePEM` package is presented is currently under review.


## Installation

Install the package from github using devtools:

```
install.packages("devtools")

devtools::install_github("klvoje/adePEM")

require(adePEM)
```

The `adePEM` package depends on the `paleoTS` package, which can be installed from CRAN. 

```
install.packages("paleoTS")

require(paleoTS)
```


## Example

We are interested in analyzing the evolution of element length (measured in mm) in the conodont Pterospathodus. The data is available as part of the `adePEM` package and was originally published by Jones (2009). The data (`element.length`) is already a `paleoTS` object. We first plot the data. 

```
plot.paleoTS(element.length)
```
![time seires](https://github.com/klvoje/adePEM/blob/master/extra/time.series.png)

Time (the x-axis) is in millions of years and the trait is measured in millimeters. Error bars represent one standard error.

We then run the `fit3models` function from the `paleoTS` package to check the relative fit of the stasis, random walk and trend models to the data.
```
fit3models(element.length, pool=TRUE)

Comparing 3 models [n = 31, method = Joint]

           logL K      AICc Akaike.wt
GRW    25.38445 3 -43.88002     0.262
URW    25.12370 2 -45.81882     0.690
Stasis 22.47400 2 -40.51943     0.049
```

The random walk (URW) model has the best fit to the data according to the AICc scores. However, the difference in the AICc score is small (<2 units) relative to the trend model (GRW). 

Let's investigate if the random walk represents an adequate statistical description of the trait dynamics in the data. To do that, we run the function `fit3adequacy.RW` from the `adePEM` package. This is a wrapper function that runs 3 adequacy tests at the same time. 
  
```
# Run adequacy tests for the random walk model:
fit3adequacy.RW(element.length)


$info
                   Value
replications     1000.00
confidence level    0.95

$summary
           estimate  min.sim max.sim p-value result
auto.corr    -0.318 -0.66474 0.36931   0.494 PASSED
runs.test   1.09003  -2.0057 4.09368    0.68 PASSED
slope.test  0.01199 -0.02635 0.04576   0.498 PASSED
```

The first part of the output summarizes the number of bootstrap replications (the number of simulated data sets) used for assessing model adequacy and the confidence level. 1000 replications and a confidence level of 0.95 are the default settings, but both can be defined by the user when running the `fit3adequacy.RW` function.

The second part of the output contains information on the results of the adequacy tests. The first column names the adequacy tests. Please check out [Voje et al. 2018 (AmNat)](http://www.journals.uchicago.edu/doi/pdfplus/10.1086/696265) for detailed info on each adequacy test. The second column gives the test statistic computed on the real data. From the second column, we see that the autocorrelation is calculated to be negative and the slope test (which is the least-squares slope of how the (detrended) data changes with time) is positive. 

The third and fourth columns reports the smallest and largest test statistics calculated on the simulated data sets. As we can see, all the three test statistics computed on the real data (second column) are not close to the extreme values reported in columns three and four. 

The fifth column is not a real p-value, but is the fraction of simulated test statistics that is larger (or smaller) than the calculated test statistic on the observed data, divided by 0.5. A value of 1 means 50 percent of the test statistics on the simulated data are larger and smaller than the calculated statistic on the observed data, respectively. A value of 0.10 means 90 percent of the test statistics on the simulated data are larger (or smaller) than the test statistic on the observed time series. 

The sixth column indicates whether our model passed the adequacy tests. Since we set our confidence level to 0.95 and all values in the fifth column is larger than 0.05, this means the random walk model passed all tests for our data set. 

That the random walk model passed all tests can also be seen in the visual representation of the distributions of test statistics, where the test statistics computed for the real data is indicated with a broken (red) line. These plots are generated automatically if `plot = TRUE` (which is the default setting) when we run the `fit3adequacy.RW` function.   

![RW distributions](https://github.com/klvoje/adePEM/blob/master/extra/adequacy.bm.png)


To summarize: Among the three candidate models stasis, random walk and directional change, random walk has the best relative model fit to the data based on AICc. However, a relative better fit for a model (in this case, the random walk model) to a phyletic fossil time series is no guarantee that the model represents a sufficiently good statistical explanation for the trait dynamics. We therefore assessed to what extent the random walk model also fitted the data in an absolute sense by running adequacy models. The random walk model passed all adequacy tests, which suggest the random walk model represents an adequate statistical description of the phyletic time series.

If we take a look at the plot of how the trait changes over 6 million years, it seems to suggest a trend towards becoming bigger. Therefore, let's assess the adequacy for the trend model on the data. This model did indeed show a quite similar fit to the data based on their AICc scores. We run the wrapper function `fit3adequacy.trend` to run all three adequacy tests simultaneously.

```
# Run adequacy tests for the trend model:
fit3adequacy.trend(element.length)

$info
                   Value
replications     1000.00
confidence level    0.95

$summary
           estimate  min.sim max.sim p-value result
auto.corr    0.0464 -0.02127 0.97608   0.004 FAILED
runs.test  -0.54272 -5.29741 0.23802   0.028 FAILED
slope.test  0.00637 -0.13183  0.2477   0.934 PASSED
```
The trend model fails the autocorrelation test and the runs test, and passes the slope test. This suggests that the trend model is not an adequate statistical description of the data.

![trend distributions](https://github.com/klvoje/adePEM/blob/master/extra/adequacy.trend.png)

Functions for running each adequacy test alone are provided in the package (e.g. `auto.corr.test.stasis`, `runs.test.RW`, `slope.test.trend`). The wrapper function for investigating the adequacy of the stasis model is `fit4adequacy.stasis`. This function runs a fourth adequacy test that is only implemented for testing the adequacy of the stasis model. A low amount of net evolution is part of the verbal definition of stasis, but not for random walk and directional trend. The `fit4adequacy.stasis` function therefore runs a test to check if the amount of net evolution is larger than expected given the  model parameters the stasis model. The function `net.change.test.stasis` runs this test alone on the data.



## Author

Kjetil L. Voje <k.l.voje@gmail.com>


## References

Hunt, G. 2006. Fitting and Comparing Models of Phyletic Evolution: Random Walks and beyond. Paleobiology 32:578–601. [link](http://www.bioone.org/doi/abs/10.1666/05070.1)

Jones, D. 2009. Directional evolution in the conodont Pterospathodus. Paleobiology 35: 413-431. [link](http://www.bioone.org/doi/abs/10.1666/0094-8373-35.3.413)

Voje, K.L., Starrfelt, J., and Liow, L.H. 2018. Model adequacy and microevolutionary explanations for stasis in the fossil record. The American Naturalist [in press]. [PDF](http://www.journals.uchicago.edu/doi/pdfplus/10.1086/696265)

