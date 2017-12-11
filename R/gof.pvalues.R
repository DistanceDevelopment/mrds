# Goodness of fit testing for detection functions
#
# See ddf.gof
#
# @param model a fitted model object
# @param nboot number of bootstraps to do
# @param ks logical whether to calculate Kolmogorov-Smirnov test statistic
# @param progress logical whether to display progress
#
# @return list with p-values for the two tests (\code{ks.p}, \code{cramer.p}) an dtest statistics (\code{Dn} for K-S and \code{W} for C-vM).
#' @importFrom utils setTxtProgressBar txtProgressBar
# @author David L Miller
#
gof.pvalues <- function(model, nboot=500, ks=FALSE, progress=FALSE){
  
  # storage
  ks.boot <- rep(0, nboot)
  cramer.boot <- rep(0, nboot)
  
  if(progress) pb <- txtProgressBar(0, nboot)
  
  boot_success <- 0
  
  # now calculate test statistics for the model
  model_edf_cdf <- get_edf_cdf(model)
  model_Dn <- max(c(abs(model_edf_cdf$lower.edf - model_edf_cdf$cdfvalues),
                    abs(model_edf_cdf$upper.edf - model_edf_cdf$cdfvalues)))
  model_W <- 1/(12*model_edf_cdf$n) + sum((model_edf_cdf$cdfvalues - ((1:model_edf_cdf$n)-.5)/model_edf_cdf$n)^2)
  
  # can skip the bootstraps if we want
  if(ks && nboot>0){
    for (i in 1:nboot) {
      # simulate data from the model
      # fit a model
      refit <- sample_ddf(model)
      
      if(!inherits(refit, "try-error")){
        ## calculate test statistics and store them
        edf_cdf <- get_edf_cdf(refit)
        # Kolmogorov-Smirnov test statistic
        ks.boot[i] <- max(c(abs(model_edf_cdf$lower.edf - edf_cdf$cdfvalues),
                            abs(model_edf_cdf$upper.edf - edf_cdf$cdfvalues)))
        boot_success <- boot_success + 1
      }
      if(progress) setTxtProgressBar(pb, i)
    }
    
    # do the boostrap test
    ks.p <- mean(model_Dn <= ks.boot, na.rm=TRUE)
  }else{
    ks.p <- NA
    model_Dn <- NA
  }
  
  # do regular Cramer-von Mises thingo
  cramer.p <- 1-pcramer(model_W)
  
  # put everything together
  res <- list(ks     = ks.p,
              cramer = cramer.p,
              Dn     = model_Dn,
              W      = model_W)
  
  attr(res, "boot_success") <- boot_success
  attr(res, "nboot") <- nboot
  
  return(res)
}
