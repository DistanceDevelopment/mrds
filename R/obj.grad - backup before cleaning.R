#' Gradient of the objective function, which is the negative log likelihood
#' 
#' The current implementation is for unbinned line transect data. 
#'
#' @param j the index of the parameter of interest
#' @param distance vector of distances
#' @param ddfobj the ddf object
#' @param width the truncation width
#' @param left the left truncation (defaults to zero)
#'
#' @author Felix Petersma
obj.grad <- function(pars, ddfobj, misc.options, fitting="all") {
  
  ## THE NEXT PART COMES FROM flpt.lnl.R >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
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

  ### <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< END HERE
  
  gradients <- list()
  
  ## If the distances are binned
  if (binned) {
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
      # return(integrate(distpdf, lower = bin.ir[1], upper = bin.ir[2],
      #                  ddfobj = ddfobj, left = left, point = point,
      #                  standardize = standardize, width = width,
      #                  select = c(TRUE, rep(FALSE, nrow(x) - 1)))$value)
    })      
    
    ## Set values in part1a.i that are zero to something really small
    # part1a.i[part1a.i < 1e-16] <- 1e-16
    
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
  } 
  ## If the distances are exact
  else {
    part1a <- distpdf(distance, ddfobj, left = left, standardize = standardize,
                      point = point, width = width)
    part2a <- length(distance) / (integratepdf(ddfobj, select = c(TRUE), 
                                               point = point,
                                               width = width, left = left,
                                               int.range = int.range,
                                               standardize = standardize))
    for (par.index in seq(along = pars)) {
      ## Part 1
      # New version as proposed on 24/05, which means we build around distpdf
      # and the defintiona g(x)/w (line) and g(x)*2x/w^2 (point). Also, took part1a
      # and part2a outside the loop
      part1b <- distpdf.grad(distance = distance, par.index = par.index, 
                             ddfobj = ddfobj, width = width, left = left, 
                             point = point, standardize = standardize)
      part1 <- sum(part1b / part1a)
      
      ## older version
      # part1a <- detfct(distance, ddfobj, left = left, standardize = standardize)
      # 
      # part1b <- detfct.grad(par.index = par.index, distance = distance, ddfobj = ddfobj,
      #                       width = width, left = left, standardize = standardize)
      # part1 <- sum(part1b / part1a)
      
      ## Alternative for Part 1
      # part1a <- distpdf(distance, ddfobj, left = left, standardize = standardize,
      #                   point = point, width = width)
      # if (point) {
      #   part1a <- part1a * width ^ 2 / 2
      # } else {
      #   part1a <- part1a * width
      # }
      # part1b <- nonnormpdf.grad(par.index = par.index, distance = distance, ddfobj = ddfobj,
      #                       width = width, left = left, standardize = standardize)
      # part1 <- sum(part1b / part1a)
      # 
      ## Initial checking gave the same results for both derivations of part 1.
      ## The first version was chosen as it is slightly quicker
      ## However, does this mean that the problem is in Part 2?
      
      ## Part 2
      ## Also here, there is a new version to keep things in line with the current
      ## implementation of code [24/05/24]
      part2b <- integratepdf.grad(par.index = par.index, 
                                  ddfobj = ddfobj, int.range = int.range, 
                                  width = width, standardize = standardize, 
                                  point = point, left = left)    
      part2 <- part2a * part2b
      # if (!point) {
      #   part2a <- nrow(x) / int.nonnormpdf(int.range = int.range, ddfobj = ddfobj,
      #                                        width = width, standardize = standardize,
      #                                        point = point, left = left)
      #   ## this one integrates g(x)/w, so multiply by w to correct
      #   # part2a <- nrow(x) / (integratepdf(ddfobj, select = c(TRUE), width = width,
      #   #                                   int.range = int.range,
      #   #                                   standardize = standardize) * width) 
      # 
      # } else {
      #   # part2a <- nrow(x) / (integratepdf(ddfobj, select = c(TRUE), width = width, 
      #   #                                   int.range = int.range,
      #   #                               standardize = standardize) / 2 * width ^ 2)
      # 
      # }
      # 
      # part2b <- int.nonnormpdf.grad(par.index = par.index, 
      #                               ddfobj = ddfobj, int.range = int.range, 
      #                               width = width, standardize = standardize, 
      #                               point = point, left = left)                      
      # 
      # 
      # part2 <- part2a * part2b
      
      ## Make sure we are getting the gradient of the NEGATIVE log likelihood
      ## I should not subtract log(distance) for point transect, but I thought I 
      ## had to? look into the theory again.
      ## UPDATE: when doing a polynomial or hermite, it DOES require the 
      ## subtraction of sum(log(distance)). I am so confused.
      if (point) { #& ddfobj$adjustment$series != "cos") {
        gradients[par.index] <- part2 - part1 #- sum(log(distance)) 
      } else {
        gradients[par.index] <- part2 - part1
      }
    }; return(gradients)
  }
  
  return(unlist(gradients))
}

neg.obj.grad <- function(pars, ddfobj, misc.options, fitting="all") {
  return(-obj.grad(pars, ddfobj, misc.options, fitting))
}
