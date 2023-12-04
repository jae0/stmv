# STMV

Space-Time Models of Variability [(STMV)](https://github.com/jae0/stmv) is a simple spatiotemporal modelling and prediction framework that addresses *non-stationarity* of spatiotemporal processes, that is where mean and variance varies across space and time. STMV comes from the tradition of *Kriging with external drift*, where external drift represents some externally imposed gradient that is first modelled at a *global* (domain-wide) level and then the residuals are modelled at a *local* (moving spatial window) level as some spatial autocorrelation process. This is, therefore, a divide and conquer approach to large space-time modelling and prediction problems. It is, therefore, an ad hoc hierarchical model. It is primarily a prediction-focussed approach, where predictions of overlapping locations are represented as variance-weighted spatial(-temporal) averages.

More specifically, external forcing ("drift") is treated separately (*global* modelled via MGCV, GLM or any other framework) from the random spatial processes. The random processes are subjected to a moving window that incrementally grow depending upon local estimates of spatial structure via autocorrelation functions. More refined results would be expected if Expectation-Maximization approaches or Gibbs sampling were to be used. However, due to computational load, this was not available when initially constructed (2010). This may change soon and I might update the approach, time permitting. 

This package provides the core functions required for prediction and inference of STMV-type spatiotemporal models. 

[For a short self-contained example, see the test example here modelling ocean bottom depths.](inst/scripts/01_bathymetry_stmv_example.md)

Other examples can be found in the various projects in https://github.com/jae0/aegis.*, but these are more involved as they use much larger data sets and spatial/temporal domains.

## Installation


To install, run the following:

```r
  remotes::install_github( "jae0/aegis")  # helper functions
  remotes::install_github( "jae0/stmv")
``` 

You probably will want to have an Rprofile set up properly such as:

```r
homedir = path.expand("~")
code_root = file.path( homedir, "bio" )   ### replace with correct path to the parent directory of your git-projects
data_root = file.path( homedir, "bio.data" )   ### replace with correct path to your data

require( aegis )
require( stmv )

```
 