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
obj.grad <- function(fpar, ddfobj, misc.options, fitting="all") {
  
  ## THE NEXT PART COMES FROM flpt.lnl.R >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
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
  
  distance <- x$distance
  width <- misc.options$width
  point <- misc.options$point
  standardize <- FALSE

  ### <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< END HERE
  
  gradients <- list()
  part1a <- distpdf(distance, ddfobj, left = left, standardize = standardize,
                    point = point, width = width)
  part2a <- nrow(x) / (integratepdf(ddfobj, select = c(TRUE), point = point,
                                    width = width, left = 0,
                                    int.range = int.range,
                                    standardize = standardize))
  for (par.index in seq(along = fpar)) {
    ## Part 1
    # New version as proposed on 24/05, which means we build around distpdf
    # and the defintiona g(x)/2 (line) and g(x)*2/w^2 (point). also, took part1a
    # and part2a outside the loop
    part1b <- distpdf.grad(par.index = par.index, distance = distance, ddfobj = ddfobj,
                          width = width, left = left, standardize = standardize)
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
    
    ## Make sure we are getting the gradient of teh NEGATIVE log likelihood
    if (point) gradients[par.index] <- part2 - part1 - sum(log(distance))
    else gradients[par.index] <- part2 - part1
  }; gradients
  
  return(unlist(gradients))
}

neg.obj.grad <- function(fpar, ddfobj, misc.options, fitting="all") {
  return(-obj.grad(fpar, ddfobj, misc.options, fitting))
}