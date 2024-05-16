#' Gradient of the objective function. 
#' 
#' The current implementation is for unbinned line transect data. 
#'
#' @param j the index of the parameter of interest
#' @param distance vector of distances
#' @param ddfobj the ddf object
#' @param width the truncation width
#' @param left the left truncation (defaults to zero)
#' 
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
  
  for (par.index in seq(along = fpar)) {
    # Part 1: d g(x)/ d(theta_j) * (1 / g(x))
    part1a <- detfct(distance, ddfobj, left = left, standardize = standardize)
    part1b <- detfct.grad(par.index = par.index, distance = distance, ddfobj = ddfobj,
                          width = width, left = left, standardize = standardize)
    part1 <- sum(part1a / part1b)
    
    part2a <- nrow(x) / integrate.detfct(int.range = int.range, ddfobj = ddfobj, 
                                         width = width, standardize = standardize, 
                                         point = point, left = left)
    part2b <- integrate.detfct.grad(distance = distance, par.index = par.index, 
                                    ddfobj = ddfobj, int.range = int.range, 
                                    width = width,bstandardize = standardize, 
                                    point = point, left = left)
    part2 <- part2a * part2b
    
    gradients[par.index] <- part1 - part2
  }
  
  return(unlist(gradients))
}