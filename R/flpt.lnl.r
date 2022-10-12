flpt.lnl <- function(fpar, ddfobj, misc.options){

  # flpt.lnl - computes negative log-likelihood values of line or point
  # transect binned/unbinned distances
  #
  # Arguments: see flnl for description of arguments
  #
  # value: vector of negative log-likelihood values for observations at current
  #        set of parameters (fpar)
  #
  # Functions Used: assign.par, tablecgf, predict(predict.smooth.spline),
  # detfct, integratedetfct, integratedetfct.logistic

  # Assign parameter values into ddfobj
  ddfobj <- assign.par(ddfobj, fpar)

  # Setup integration ranges
  int.range <- misc.options$int.range
  if(is.vector(int.range)){
    int.range <- matrix(int.range, nrow=1)
    samelimits <- TRUE
  }else{
    #int.range <- int.range[2:nrow(int.range), , drop=FALSE]
    samelimits <- FALSE
  }
  left <- int.range[, 1]
  right <- int.range[, 2]

  x <- ddfobj$xmat
  z <- ddfobj$scale$dm
  if(is.null(z)){
    z <- matrix(1, nrow=nrow(x), ncol=1)
  }

  width <- misc.options$width

  # set "standardize" to FALSE, as this cancels in the likelihood
  # evaluation. Note that it has to be the same value (TRUE/FALSE)
  # through all computations (for both detection function and integral)
  # or the below makes no sense
  standardize <- FALSE

  # Compute log-likelihood for any binned data
  lnl <- rep(0, nrow(x))
  if(any(x$binned)){
    # we only need evaluate the likelihood at unique bin-covariate combinations

    # Get bins and create unique set of bins/covariates and indices
    # (int.index) in that set
    bins <- as.matrix(x[x$binned, c("distbegin", "distend")])
    allbins <- apply(cbind(bins, z[x$binned, , drop=FALSE]), 1,
                     paste, collapse="")

    uniquevals <- !duplicated(allbins)
    uniquebins <- bins[uniquevals, , drop=FALSE]

    int.index <- match(allbins, allbins[uniquevals])

    # which were those observations?
    which.obs <- x$binned
    which.obs[!uniquevals] <- FALSE
    # evaluate those observations only!
    int.bin <- integratepdf(ddfobj, select=which.obs, width=width,
                            int.range=uniquebins,
                            standardize=standardize, point=misc.options$point,
                            left=misc.options$left)

    if(any(int.bin<0)){
      int.bin <- integratepdf(ddfobj, select=which.obs, width=width,
                              int.range=uniquebins,
                              standardize=standardize, point=misc.options$point,
                              left=misc.options$left)
    }

    if(any(int.bin<=0) && misc.options$showit > 0){
      warning("\nDetection function integral <=0. Setting integral to 1E-25\n")
      int.bin[int.bin<=0] <- 1E-25
    }

    int.bin <- int.bin[int.index]

    # Compute integral from left to right
    if(ddfobj$intercept.only & samelimits){
      int.all <- integratepdf(ddfobj, select=c(TRUE, rep(FALSE, nrow(x)-1)),
                              width=width,int.range=int.range,
                              standardize=standardize,
                              point=misc.options$point, left=misc.options$left)
    }else{
      if(nrow(int.range)==1){
        int.range <- int.range[rep(1, nrow(x)), ]
      }

      bins <- int.range[x$binned,,drop=FALSE]
      allbins <- apply(cbind(bins,z[x$binned,,drop=FALSE]),1,paste,collapse="")

      uniquevals <- !duplicated(allbins)
      uniquebins <- bins[uniquevals, , drop=FALSE]

      int.index <- match(allbins, allbins[uniquevals])

      which.obs <- x$binned
      which.obs[!uniquevals] <- FALSE

      int.all <- integratepdf(ddfobj, select=which.obs, width=width,
                              int.range=uniquebins,
                              standardize=standardize, point=misc.options$point,
                              left=misc.options$left)
      int.all <- int.all[int.index]
    }

    # Replace infinite integral values
    if(is.vector(left)){
      int.all[is.infinite(int.all)] <- right - left
    }else{
      int.all[is.infinite(int.all)] <- right[is.infinite(int.all)] -
                                        left[is.infinite(int.all)]
    }

    # Negative log-likelihood values for binned data
    lnl[x$binned] <- -(log(int.bin) -log(int.all))
  }

  # Compute log-likelihood for any unbinned data
  if(!all(x$binned)){
    p1 <- distpdf(x$distance[!x$binned], ddfobj=ddfobj, select=!x$binned,
                  width=right, standardize=standardize,
                  point=misc.options$point, left=left)
    p1[p1<1.0e-15] <- 1.0e-15
    p1[is.nan(p1)] <- 1.0e-15

    # Compute integrals - repeat with doeachint if not set and any 
    # integral values < 0
    int1 <- -1
    i <- 0
    doeachint <- FALSE
    while(any(int1 < 0) & (i < 2)){
      if(ddfobj$intercept.only & samelimits){
        int1 <- integratepdf(ddfobj,
                             select=c(TRUE, rep(FALSE, nrow(ddfobj$xmat))),
                             width=width, int.range=int.range,
                             doeachint=doeachint,
                             point=misc.options$point, standardize=standardize,
                             left=left)
      }else{
        if(nrow(int.range)>1){
          intrange <- int.range
        }else{
          intrange <- int.range[rep(1, nrow(x)), ]
        }

        int1 <- integratepdf(ddfobj, select=!x$binned, width=right,
                             int.range=intrange[!x$binned, ],
                             doeachint=doeachint,
                             point=misc.options$point, standardize=standardize,
                             left=left)
      }

      doeachint <- TRUE
      # set to -1 to flag above if we get NaNs
      int1[is.nan(int1)] <- -1
      i <- i + 1
    }
    # if the integral is <=0 something Very Bad has happened
    # but usually we recover as this is just a bad par combination
    if(any(int1<=0) && misc.options$showit > 0){
      int1[int1<=0] <- NaN
      warning("\n Problems with integration (integral <=0).\n")
    }

    if(is.matrix(int1) && dim(int1)==c(1, 1)){
      int1 <- as.vector(int1)
    }
    # Negative log-likelihood values for unbinned data
    lnl[!x$binned] <- -(log(p1) -log(int1))
  }

  if(any(is.nan(lnl))){
    lnl <- NA
  }

  return(lnl)
}
