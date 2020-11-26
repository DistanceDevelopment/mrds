mrds - Mark-Recapture Distance Sampling
=======================================


[![Build Status](https://travis-ci.org/DistanceDevelopment/mrds.svg?branch=master)](https://travis-ci.org/DistanceDevelopment/mrds)
[![CRAN (RStudio Mirror) Downloads](http://cranlogs.r-pkg.org/badges/mrds)](https://www.r-pkg.org/pkg/mrds)
[![CRAN Version](http://www.r-pkg.org/badges/version/mrds)](https://cran.r-project.org/package=mrds)
[![Codecov test coverage](https://codecov.io/gh/DistanceDevelopment/mrds/branch/master/graph/badge.svg)](https://codecov.io/gh/DistanceDevelopment/mrds?branch=master)

# What is `mrds`?

This package for R analyzes single or double observer distance sampling data for line or point sampling.  It is used in program [DISTANCE](http://distancesampling.org/) as one of the analysis engines. Supported double observer configurations include independent, trial and removal. Not all options in mrds are fully supported via DISTANCE.

If you only wish to perform a conventional or multiple covariate distance sampling analysis (CDS/MCDS) (as opposed to a double observer analysis), you may want to try the `Distance` R package, which has a simplified interface and is available from [https://github.com/DistanceDevelopment/Distance](https://github.com/DistanceDevelopment/Distance).

# Getting `mrds`

The easiest way to ensure you have the latest version of `mrds`, is to install using the `remotes` package:

      install.packages("remotes")

then install `mrds` from github:

      library(remotes)
      install_github("DistanceDevelopment/mrds")

Otherwise:

  * One can download a Windows package binary using the ["Releases" tab in github](https://github.com/DistanceDevelopment/mrds/releases). To install in R, from the R menu, use "Packages\Install from Local Zip file" and browse to location of downloaded zip.
  * Or, download [package source files](https://github.com/DistanceDevelopment/mrds/archive/master.zip).
  * Finally the current stable version of `mrds` is [available on CRAN](https://cran.r-project.org/package=mrds), though this may be up to a month out of date due to CRAN policy.


# References

The following are references for the methods used in the package.

Burt, M. L., D. L. Borchers, K. J. Jenkins and T. A. Marques. (2014). "[Using mark-recapture distance sampling methods on line transect surveys](https://besjournals.onlinelibrary.wiley.com/doi/abs/10.1111/2041-210X.12294)." Methods in Ecology and Evolution 5: 1180-1191.

Buckland, S. T., J. Laake, et al. (2010). "Double observer line transect methods: levels of independence." Biometrics 66: 169-177.

Borchers, D. L., J. L. Laake, et al. (2006). "Accommodating unmodeled heterogeneity in double-observer distance sampling surveys." Biometrics 62(2): 372-378.

Buckland, S. T., D. R. Anderson, et al., Eds. (2004). Advanced distance sampling: estimating abundance of biological populations. Oxford, UK; New York, Oxford University Press. (see chapter 6).
