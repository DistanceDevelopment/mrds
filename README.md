mrds - Mark-Recapture Distance Sampling
=======================================

This package for R analyzes single or double observer distance sampling data for line or
point sampling.  It is used in program [DISTANCE](http://www.ruwpa.st-and.ac.uk/distance/) as one of the analysis engines. 
Supported double observer configurations include independent, trial and removal. All options in mrds are not yet fully supported via DISTANCE.

Download [Windows package binary](https://github.com/downloads/jlaake/mrds/mrds_2.0.9.zip). From R menu, use Packages\Install from Local Zip file and browse to location of downloaded zip. As an
alternative or for other operating systems, install the [devtools package](http://cran.r-project.org/web/packages/devtools/index.html) and use the 
following (may require [Rtools](cran.r-project.org/bin/windows/Rtools/):

```
library(devtools)
install_github("mrds",user="jlaake",subdir="mrds")
```

Download [package source files](https://github.com/jlaake/mrds/archive/master.zip)

The following are references for the methods used in the package.

Buckland, S. T., J. Laake, et al. (2010). "Double observer line transect methods: levels of independence." Biometrics 66: 169-177.

Borchers, D. L., J. L. Laake, et al. (2006). "Accommodating unmodeled heterogeneity in double-observer distance sampling surveys." Biometrics 62(2): 372-378.

Buckland, S. T., D. R. Anderson, et al., Eds. (2004). Advanced distance sampling: estimating abundance of biological populations. Oxford, UK; New York, Oxford University Press. (see chapter 6).