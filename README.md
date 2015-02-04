mrds - Mark-Recapture Distance Sampling
=======================================

[![Build Status](https://travis-ci.org/DistanceDevelopment/mrds.svg?        branch=master)](https://travis-ci.org/DistanceDevelopment/mrds)

# What is `mrds`?

This package for R analyzes single or double observer distance sampling data for line or point sampling.  It is used in program [DISTANCE](http://distancesampling.org/) as one of the analysis engines. Supported double observer configurations include independent, trial and removal. All options in mrds are not yet fully supported via DISTANCE.

# Getting `mrds`

The easiest way to ensure you have the latest version of `mrds`, is to install Hadley Wickham's `devtools` package:

      install.packages("devtools")

then install `mrds` from github:

      library(devtools)
      install_github("DistanceDevelopment/mrds")

Otherwise:

  * One can download a Windows package binary using the ["Releases" tab in github](https://github.com/DistanceDevelopment/mrds/releases). To install in R, from the R menu, use "Packages\Install from Local Zip file" and browse to location of downloaded zip. 
  * Or, download [package source files](https://github.com/jlaake/mrds/archive/master.zip).
  * Finally the current stable version of `mrds` is [available on CRAN](http://cran.r-project.org/web/packages/mrds/index.html), though this may be up to a month out of date due to CRAN policy.


# References

The following are references for the methods used in the package.

Buckland, S. T., J. Laake, et al. (2010). "Double observer line transect methods: levels of independence." Biometrics 66: 169-177.

Borchers, D. L., J. L. Laake, et al. (2006). "Accommodating unmodeled heterogeneity in double-observer distance sampling surveys." Biometrics 62(2): 372-378.

Buckland, S. T., D. R. Anderson, et al., Eds. (2004). Advanced distance sampling: estimating abundance of biological populations. Oxford, UK; New York, Oxford University Press. (see chapter 6).
