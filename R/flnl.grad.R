#' @title Gradient of the negative log likelihood function
#' 
#' This function derives the gradients of the negative log likelihood function, 
#' with respect to all parameters. It is based on the theory presented in 
#' Introduction to Distance Sampling (2001) and Distance Sampling: Methods and
#' Applications (2015). It is not meant to be called by users of the \code{mrds}
#' and \code{Distance} packages directly but rather by the gradient-based 
#' solver. This solver is used when our distance sampling model is for 
#' single-observer data coming from either line or point transect and only when  
#' the detection function contains an adjustment series but no covariates. It is 
#' implement for the following key + adjustment series combinations for the 
#' detections function: the key function can be half-normal, hazard-rate or 
#' uniform, and the adjustment series can be cosine, simple polynomial or 
#' Hermite polynomial. Data can be either binned or exact, but a combination
#' of the two has not been implemented yet. 
#' 
#' @param pars vector of parameter values for the detection function at which 
#' the gradients of the negative log-likelihood should be evaluated
#' @param ddfobj distance sampling object 
#' @param misc.options a list object containing all additional information such 
#' as the type of optimiser or the truncation width, and is created by 
#' \code{\link{ddf.ds}}
#' @param fitting character string with values "all", "key", "adjust" to
#'   determine which parameters are allowed to vary in the fitting. Not actually
#'   used. Defaults to "all". 
#'   
#' @return The gradients of the negative log-likelihood w.r.t. the parameters
#'
#' @author Felix Petersma
flnl.grad <- function(pars, ddfobj, misc.options, fitting = "all") {
  
  # Assign parameter values into ddfobj
  ddfobj <- assign.par(ddfobj, pars)
  
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
  
  distance <- x$distance
  width <- misc.options$width
  point <- misc.options$point
  binned <- misc.options$binned
  standardize <- FALSE

  gradients <- list()
  
  ## If the distances are binned
  if (binned) {
    if (all(x$binned)) {
      # Get bins and create unique set of bins/covariates and indices
      # (int.index) in that set
      bins <- x[x$binned, c("distbegin", "distend")]
      allbins <- apply(cbind(bins, z[x$binned, , drop=FALSE]), 1,
                       paste, collapse="")
      bins$id <- allbins
      
      ## Values per bin
      freq.bins <- table(allbins)
      
      uniquevals <- !duplicated(allbins)
      uniquebins <- bins[uniquevals, , drop=FALSE]
      
      int.index <- match(allbins, allbins[uniquevals])
      
      ## not sure if this is necessary, but I think it might be  a little faster
      ddfobj2 <- ddfobj
      ddfobj2$xmat$binned <- FALSE
      
      # which were those observations?
      which.obs <- x$binned
      which.obs[!uniquevals] <- FALSE
      
      ## Deriva part1a for every bin i
      part1a.i <- sapply(1:length(freq.bins), function(i) {
        m <- freq.bins[i]
        id <- names(m)
        bin.ir <- unlist(uniquebins[uniquebins$id == id, 
                                    c("distbegin", "distend")])
        
        ## Used integrate(distpdf) instead of integratepdf. Both work.
        ## ddfobj$xmat$binned is TRUE, which is not a porlbem since I use
        ## select = c(TRUE), I think. To be sure, create second ddfobj outside of 
        ## this loop
        return(integratepdf(ddfobj2, select = c(TRUE),
                            point = point,
                            width = width, left = left,
                            int.range = bin.ir,
                            standardize = standardize))
      })      

      
      ## Sum all part1a.i to get part2a
      part2a <- sum(part1a.i)
      
      ## Now loop over the parameter indices
      for (par.index in seq(along = pars)) {
        ## Derive part 1b for every bin i
        part1b.i <- sapply(1:length(freq.bins), function(i) {
          m <- freq.bins[i]
          id <- names(m)
          bin.ir <- unlist(uniquebins[uniquebins$id == id, 
                                      c("distbegin", "distend")])
          
          return(integratepdf.grad(int.range = bin.ir, par.index = par.index,
                                   ddfobj = ddfobj, left = left, point = point,
                                   standardize = standardize, width = width))
        })
        ## Derive part 2b by summing all part1b.i values
        part2b <- sum(part1b.i)
        
        ## Combine to get part 1 and part 2
        part1 <- part1b.i / part1a.i
        part2 <- part2b / part2a 
        
        ## Subtract those, multiply by the bin size and sum the outcomes to get 
        ## the gradient
        gradients[par.index] <- sum(freq.bins * (part2 - part1))
      }
    } else {
    #   ## This part covers the case when some data are binned and some are not. 
    #   ## Separate the data into binned data and exact data.
    #   
    #   ## FTP update: A first setup has been added but it has not been properly 
    #   ## implemented yet
    #   
    #   ## BINNED DATA PART <<<<<<<<<<<<<<<<<
    #   ddfobj.binned <- ddfobj
    #   ddfobj.binned$xmat <- x[x$binned, ]
    #   # misc.options.binned <- misc.options
    #   x.binned <- x[x$binned, ]
    #   
    #   bins <- x.binned[, c("distbegin", "distend")]
    #   allbins <- apply(cbind(bins, z[x$binned, , drop=FALSE]), 1,
    #                    paste, collapse="")
    #   bins$id <- allbins
    #   
    #   ## Values per bin
    #   freq.bins <- table(allbins)
    #   
    #   uniquevals <- !duplicated(allbins)
    #   uniquebins <- bins[uniquevals, , drop=FALSE]
    #   
    #   int.index <- match(allbins, allbins[uniquevals])
    #   
    #   ## not sure if this is necessary, but I think it might be  a little faster
    #   ddfobj.binned.2 <- ddfobj.binned
    #   ddfobj.binned.2$xmat$binned <- FALSE
    #   
    #   # which were those observations?
    #   which.obs <- x$binned
    #   which.obs[!uniquevals] <- FALSE
    #   
    #   ## Deriva part1a for every bin i
    #   part1a.i.binned <- sapply(1:length(freq.bins), function(i) {
    #     m <- freq.bins[i]
    #     id <- names(m)
    #     bin.ir <- unlist(uniquebins[uniquebins$id == id, 
    #                                 c("distbegin", "distend")])
    #     
    #     ## Used integrate(distpdf) instead of integratepdf. Both work.
    #     ## ddfobj$xmat$binned is TRUE, which is not a porlbem since I use
    #     ## select = c(TRUE), I think. To be sure, create second ddfobj outside of 
    #     ## this loop
    #     return(integratepdf(ddfobj.binned.2, select = c(TRUE),
    #                         point = point,
    #                         width = width, left = left,
    #                         int.range = bin.ir,
    #                         standardize = standardize))
    #   })      
    #   
    #   ## Set values in part1a.i that are zero to something really small
    #   # part1a.i[part1a.i < 1e-16] <- 1e-16
    #   
    #   ## Sum all part1a.i to get part2a
    #   part2a.binned <- sum(part1a.i.binned)
    #   
    #   
    #   gradients.binned <- sapply(seq(along = pars), function(par.index){
    #     
    #     ## Derive part 1b for every bin i
    #     part1b.i.binned <- sapply(1:length(freq.bins), function(i) {
    #       m <- freq.bins[i]
    #       id <- names(m)
    #       bin.ir <- unlist(uniquebins[uniquebins$id == id, 
    #                                   c("distbegin", "distend")])
    #       
    #       return(integratepdf.grad(int.range = bin.ir, par.index = par.index,
    #                                ddfobj = ddfobj.binned, left = left, point = point,
    #                                standardize = standardize, width = width))
    #     })
    #     ## Derive part 2b by summing all part1b.i values
    #     part2b.binned <- sum(part1b.i.binned)
    #     
    #     ## Combine to get part 1 and part 2
    #     part1.binned <- part1b.i.binned / part1a.i.binned
    #     part2.binned <- part2b.binned / part2a.binned
    #     
    #     ## Subtract those, multiply by the bin size and sum the outcomes to get 
    #     ## the gradient
    #     return(sum(freq.bins * (part2.binned - part1.binned)))
    #     
    #   })
    # 
    #   ## EXACT DATA PART <<<<<<<<<<<<<<<<<<<<
    #   ddfobj.exact <- ddfobj
    #   ddfobj.exact$xmat <- x[!x$binned, ]
    #   # misc.options.exact <- misc.options
    #   # misc.options.exact$binned <- FALSE
    #   distance.exact <- ddfobj.exact$xmat$distance
    #   
    #   part1a.exact <- distpdf(distance.exact, ddfobj.exact, left = left, 
    #                     standardize = standardize,
    #                     point = point, width = width)
    #   part2a.exact <- length(distance.exact) / (integratepdf(ddfobj.exact, 
    #                                                    select = c(TRUE), 
    #                                                    point = point,
    #                                                    width = width, left = left,
    #                                                    int.range = int.range,
    #                                                    standardize = standardize))
    #   gradients.exact <- sapply(seq(along = pars), function(par.index) {
    #     ## Part 1
    #     # New version as proposed on 24/05, which means we build around distpdf
    #     # and the defintiona g(x)/w (line) and g(x)*2x/w^2 (point). Also, took part1a
    #     # and part2a outside the loop
    #     part1b.exact <- distpdf.grad(distance = distance.exact, par.index = par.index, 
    #                            ddfobj = ddfobj.exact, width = width, left = left, 
    #                            point = point, standardize = standardize)
    #     part1.exact <- sum(part1b.exact / part1a.exact)
    #     
    #     ## Part 2
    #     ## Also here, there is a new version to keep things in line with the current
    #     ## implementation of code [24/05/24]
    #     part2b.exact <- integratepdf.grad(par.index = par.index, 
    #                                 ddfobj = ddfobj.exact, int.range = int.range, 
    #                                 width = width, standardize = standardize, 
    #                                 point = point, left = left)    
    #     part2.exact <- part2a.exact * part2b.exact
    #     
    #     ## Make sure we are getting the gradient of the NEGATIVE log likelihood
    #     ## I should not subtract log(distance) for point transect, but I thought I 
    #     ## had to? look into the theory again.
    #     ## UPDATE: when doing a polynomial or hermite, it DOES require the 
    #     ## subtraction of sum(log(distance)). I am so confused.
    #     if (point) { #& ddfobj$adjustment$series != "cos") {
    #       return(part2.exact - part1.exact) #- sum(log(distance)) 
    #     } else {
    #       return(part2.exact - part1.exact)
    #     }
    #   })
    #   
    #   gradients <- gradients.binned + gradients.exact
    #   return(gradients)
    }
  } else { ## If the distances are exact
    part1a <- distpdf(distance, ddfobj, left = left, standardize = standardize,
                      point = point, width = width)
    ## Set 0 values to something really small
    part1a[part1a < 1e-16] <- 1e-16
    
    part2a <- length(distance) / integratepdf(ddfobj, select = c(TRUE), 
                                              point = point,
                                              width = width, left = left,
                                              int.range = int.range,
                                              standardize = standardize)
    
    for (par.index in seq(along = pars)) {
      ## Part 1
      # New version as proposed on 24/05, which means we build around distpdf
      # and the definitions g(x)/w (line) and g(x)*2x/w^2 (point). Also, took part1a
      # and part2a outside the loop
      part1b <- distpdf.grad(distance = distance, par.index = par.index,
                             ddfobj = ddfobj, width = width, left = left,
                             point = point, standardize = standardize)
      part1 <- sum(part1b / part1a)

      ## Part 2
      ## Also here, there is a new version to keep things in line with the current
      ## implementation of code [24/05/24]
      part2b <- integratepdf.grad(par.index = par.index,
                                  ddfobj = ddfobj, int.range = int.range,
                                  width = width, standardize = standardize,
                                  point = point, left = left)
      part2 <- part2a * part2b

      ## Make sure we are getting the gradient of the NEGATIVE log likelihood
      ## I should not subtract log(distance) for point transect, but I thought I
      ## had to? look into the theory again.
      gradients[par.index] <- part2 - part1 #- sum(log(distance))
    }
  }
  return(unlist(gradients))
}

## Negative of the gradients. Not in use now, but could be useful for future
## implementations. 
flnl.grad.neg <- function(pars, ddfobj, misc.options, fitting = "all") {
  return(-flnl.grad(pars, ddfobj, misc.options))
}
