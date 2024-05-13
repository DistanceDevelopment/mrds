#' Gradient of the objective function
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
  
  width <- misc.options$width
  point <- misc.options$point
  standardize = TRUE

  ### <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< END HERE
  
  gradients <- list()
  
  for (par.index in seq.along(fpar)) {
    # Part A: d g(x)/ d(theta_j) * (1 / g(x))
    A1 <- detfct(x$distance, ddfobj, left = left, standardize = standardize)
    A2 <- detfct.grad(par.index = par.index, distance = x$distance, ddfobj = ddfobj,
                width = width, left = left)
    A <- sum(A1 / A2)
    
    B1 <- nrow(x) / beta(int.range = int.range, ddfobj = ddfobj, width = width, 
                         standardize = standardize, point = point, left = left)
    B2 <- grad.beta()
    
    B <- B1 * B2
    
    gradients[par.index] <- A - B
  }
}