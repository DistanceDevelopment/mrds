flpt.lnl <- function(fpar,ddfobj,misc.options){

  # flpt.lnl - computes negative log-likelihood values of line or point 
  # transect binned/unbinned distances
  #
  # Arguments: see flnl for description of arguments
  # 
  # value: vector of negative log-likelihood values for observations at current
  #        set of parameters (fpar)
  #
  # Functions Used: assign.par, tablecgf, predict(predict.smooth.spline), 
  # detfct, integratedetfct, integratedetfct.logistic

  # Assign parameter values into ddfobj
  ddfobj <- assign.par(ddfobj,fpar)

  # Setup integration ranges
  int.range <- misc.options$int.range
  if(is.vector(int.range)){
    int.range <- matrix(int.range,nrow=1)
    samelimits <- TRUE
  }else{
    int.range <- int.range[2:nrow(int.range),,drop=FALSE]
    samelimits <- FALSE
  }
  left <- int.range[,1]
  right <- int.range[,2]

  # If detection function has a shape parameter or adjustments are used
  # and !doeachint, compute the integration spline with tablecgf with this
  # set of parameter values
  doeachint <- misc.options$doeachint
  x <- ddfobj$xmat
  z <- ddfobj$scale$dm
  if(is.null(z)){
    z <- matrix(1,nrow=nrow(x),ncol=1)
  }

  width <- misc.options$width
  if((!is.null(ddfobj$shape) | !is.null(ddfobj$adjustment)) & !doeachint){
    ddfobj$cgftab <- tablecgf(ddfobj,width=width,standardize=FALSE,
                              point=misc.options$point)
  }

  # Compute log-likelihood for any binned data
  lnl <- rep(0,dim(x)[1])
  if(any(x$binned)){
    # Get bins and create unique set of bins/covariates and indices 
    # (int.index) in that set
    bins <- as.matrix(x[x$binned,c("distbegin","distend")])
    allbins <- apply(cbind(bins,z[x$binned,,drop=FALSE]),1,paste,collapse="")

    uniquevals <- !duplicated(allbins)
    uniquebins <- bins[uniquevals,,drop=FALSE]

    int.index <- match(allbins,
                       apply(cbind(bins,
                              z[x$binned,,drop=FALSE])[uniquevals,,drop=FALSE],
                              1,paste,collapse=""))

    which.obs <- x$binned
    which.obs[!uniquevals] <- FALSE
    int.bin <- integratepdf(ddfobj,select=which.obs,width=width,
                            int.range=uniquebins,doeachint=doeachint,
                            standardize=FALSE,point=misc.options$point)

    if(any(int.bin<0)){
      int.bin <- integratepdf(ddfobj,select=which.obs,width=width,
                              int.range=uniquebins,doeachint=TRUE,
                              standardize=FALSE,point=misc.options$point)
    }

    if(any(int.bin<0)){
      warning("\nProblems with integration. integral <0. Setting prob=0\n")
      int.bin[int.bin<0] <- 0
    }

    int.bin <- int.bin[int.index]

    # Compute integral from left to right
    if(ddfobj$intercept.only & samelimits){
      int.all <- integratepdf(ddfobj,select=c(TRUE,rep(FALSE,nrow(x)-1)),
                              width=width,int.range=int.range,
                              doeachint=doeachint,standardize=FALSE,
                              point=misc.options$point)
    }else{
      if(nrow(int.range)==1){
        int.range <- int.range[rep(1,nrow(x)),]
      }

      bins <- int.range[x$binned,,drop=FALSE]
      allbins <- apply(cbind(bins,z[x$binned,,drop=FALSE]),1,paste,collapse="")

      uniquevals <- !duplicated(allbins)
      uniquebins <- bins[uniquevals,,drop=FALSE]

      int.index <- match(allbins,
                         apply(cbind(bins,
                               z[x$binned,,drop=FALSE])[uniquevals,,drop=FALSE],
                               1,paste,collapse=""))

      which.obs=x$binned
      which.obs[!uniquevals]=FALSE

      int.all <- integratepdf(ddfobj,select=which.obs,width=width,
                              int.range=uniquebins,doeachint=doeachint,
                              standardize=FALSE,point=misc.options$point)
      int.all <- int.all[int.index]
    }

    # Replace infinite integral values
    if(is.vector(left)){
      int.all[is.infinite(int.all)]<- right - left
    }else{
      int.all[is.infinite(int.all)]<- right[is.infinite(int.all)] -
                                       left[is.infinite(int.all)]
    }

    # Negative log-likelihood values for binned data
    lnl[x$binned] <- -log(int.bin/int.all)
  }

  # Compute log-likelihood for any unbinned data
  if(!all(x$binned)){
    p1 <- distpdf(x$distance[!x$binned],ddfobj=ddfobj,select=!x$binned,
                  width=width,standardize=FALSE,point=misc.options$point)
    p1[p1<1.0e-15] <- 1.0e-15
    p1[is.nan(p1)] <- 1.0e-15
    # Compute integrals - repeat with doeachint if not set and any 
    # integral values < 0
    int1 <- -1
    i <- 0
    while (any(int1 < 0) & (i < 2)){
      if(ddfobj$intercept.only & samelimits){
        int1 <- integratepdf(ddfobj,select=c(TRUE,rep(FALSE,nrow(ddfobj$xmat))),
                           width=width,int.range=int.range,doeachint=doeachint,
                           point=misc.options$point,standardize=FALSE)
      }else{
        if(nrow(int.range)>1){
          intrange <- int.range
        }else{
          intrange <- int.range[rep(1,nrow(x)),]
        }

        int1 <- integratepdf(ddfobj,select=!x$binned,width=width,
                             int.range=intrange[!x$binned,],doeachint=doeachint,
                             point=misc.options$point,standardize=FALSE)
        }
      doeachint <- TRUE
      i <- i + 1
    }
    if(any(int1<0)){
      stop("\n Problems with integration. One or more integrals <0")
    }

    # Handle some special cases that can occur 
    if(is.vector(left)){
      int1[is.infinite(int1)] <-  right - left
      int1[is.nan(int1)] <- right - left
    }else{ 
      int1[is.infinite(int1)] <-  right[is.infinite(int1)] -
                                  left[is.infinite(int1)]
      int1[is.nan(int1)] <- right[is.nan(int1)] - left[is.nan(int1)]
    }
    # Negative log-likelihood values for unbinned data
    lnl[!x$binned] <- -log(p1/int1)
  }

  if(any(is.nan(lnl))){
    lnl <- NA
  }

  return(lnl)
}
