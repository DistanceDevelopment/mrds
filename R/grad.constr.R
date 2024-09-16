#' Gradients of constraint function
#'
#' The function derives the gradients of the constraint function for 
#' all model parameters, in the following order:
#' 1. Scale parameter (if part of key function)
#' 2. Shape parameter (if part of key function)
#' 3. Adjustment parameter 1
#' 4. Adjustment parameter 2
#' 5. Etc. 
#' 
#' The constraint function itself is formed of a specified number of non-linear
#' constraints, which defaults to 20 and is specified through 
#' \code{misc.options$mono.points}. The constraint function checks whether the 
#' standardised detection function is 1) weakly/strictly monotonic at the 
#' points and 2) non-negative at all the points. \code{grad.constr} returns 
#' the gradients of those constraints w.r.t. all parameters of the detection
#' function, i.e., 2 times \code{mono.points} gradients for every parameter. 
#'
#' This function mostly follows the same structure as \code{\link{flnl.constr}} in 
#' \code{detfct.fit.mono.R}. 
#' 
#' @param pars vector of parameter values for the detection function at which 
#' the gradients of the negative log-likelihood should be evaluated
#' @param ddfobj distance sampling object 
#' @param misc.options a list object containing all additional information such 
#' as the type of optimiser or the truncation width, and is created by 
#' \code{\link{ddf.ds}}
#' 
#' @returns a matrix of gradients for all constraints (rows) w.r.t to every
#' parameters (columns)
grad.constr.neg <- function(pars, ddfobj, misc.options, fitting = "all") {
  
  if (is.null(ddfobj$adjustment)) {
    # This never gets called from ddf()
    stop("This should not happen")
  } else {
    ddfobj <- assign.par(ddfobj, pars)
    
    ## Apply the constraints only if mono == TRUE
    constr <- misc.options$mono
    ## Apply strict monotonicity only if mono.strict == TRUE
    strict <- misc.options$mono.strict
    
    ## Start the constraint function evaluation
    ## ========================================
    ## Number of points (distances) at which the DF is evaluated
    no.points <- misc.options$mono.points
    
    ## Extract the reference points
    ref.points <- getRefPoints(no.points, misc.options$int.range)
    
    ## Reference point associated with distance=0
    ref.point0 <- 0
    
    ## Extract some additional information
    no.pars <- length(pars)      # number of parameters
    key <- ddfobj$type           # the key function type
    width <- misc.options$width  # the truncation width
    left <- misc.options$left    # the left truncation distance
    point <- misc.options$point  # whether the data comes from point transects
    no.data <- nrow(ddfobj$xmat)
    
    # ## To get detfct to play nice need to mudge ddfobj a bit...
    ## MADE REDUNDANT BY USING select
    # if(!is.null(ddfobj$scale)){
    #   ddfobj$scale$dm <- rep(1, no.points)
    # }
    # if(!is.null(ddfobj$shape)){
    #   ddfobj$shape$dm <- rep(1, no.points)
    # }
    
    # g.refvals <- distpdf(distance = ref.points, ddfobj = ddfobj,
    #                      width = width, left = left, point = point,
    #                      select = c(rep(TRUE, no.points), 
    #                                 rep(FALSE, no.data - no.points)))
    # if (!point) {
    #   g.refvals <- g.refvals * (width - left)
    # } else {
    #   g.refvals <- g.refvals / (2 * ref.points / (width ^ 2))
    # }
    
    g.refvals <- detfct(distance = ref.points, ddfobj = ddfobj,
                        width = width, left = left,
                        select = c(rep(TRUE, no.points), 
                                   rep(FALSE, no.data - no.points)))
    
    # ## Again, to get detfct to play nice need to mudge ddfobj a bit...
    ## MADE REDUNDANT BY USING select
    # if(!is.null(ddfobj$scale)){
    #   ddfobj$scale$dm <- 1
    # }
    # if(!is.null(ddfobj$shape)){
    #   ddfobj$shape$dm <- 1
    # }
    
    g.0 <- detfct(distance = ref.point0, ddfobj = ddfobj, width = width,
                  left = left, select = c(TRUE, rep(FALSE, no.data - 1)))
    
    # g.0 <- distpdf(distance = ref.point0, ddfobj = ddfobj, width = width,
    #                left = left, point = point, 
    #                select = c(TRUE, rep(FALSE, no.data - 1)))
    
    ## Add a very small number if g.0 = 0 to avoid dividing by zero, however,
    ## I am not sure if this is correct since dividing by a very small number is 
    ## kinda sketchy. 
    # if (g.0 == 0) g.0 <- 1e-10
    # g.0 <- 1
    
    
    ## Create matrix to start the constraint values, which are
    ## no.points values for the monotonicity constraint and no.points values
    ## for the positivity constraint, for no.pars parameters.
    grad.constraints <- sapply(1:no.pars, function(par.index) {
      ## Evaluate the gradient of the detection function at the reference points
      ## Note that we must standardize so 0<=g(x)<=1. this is done at line 124
      dg.0 <- distpdf.grad(distance = ref.point0, par.index = par.index,
                           ddfobj = ddfobj, standardize = FALSE,
                           point = point, left = left, width = width,
                           pdf.based = FALSE) #,
      # select = c(TRUE, rep(FALSE, no.data - 1)))
      # if (!point) {
      #   dg.0 <- dg.0 * (width - left)
      # } else {
      #   dg.0 <- dg.0 / (2 * ref.point0 / (width ^ 2))
      # }
      
      # dg.0 <- 0
      ## Standardised gradient at distance 0 is always 0 for line (but not point?)
      # std.dg.0 <- (dg.0 - g.0 * dg.0 / g.0) / g.0 ## Always 0
      std.dg.0 <- 0
      
      dg.refvals <- distpdf.grad(distance = ref.points, par.index = par.index, 
                                 ddfobj = ddfobj, standardize = FALSE, 
                                 point = point, left = left, width = width,
                                 pdf.based = FALSE) #,
      # select = c(rep(TRUE, no.points), 
      # rep(FALSE, no.data - no.points)))
      # if (!point) {
      #   dg.refvals <- dg.refvals * (width - left)
      # } else {
      #   dg.refvals <- dg.refvals / (2 * ref.points / (width ^ 2))
      # }
      ## Derive the gradients of the standardised detection function [see notes 17-06-2024]
      std.dg.refvals <- (dg.refvals - g.refvals * dg.0 / g.0) / g.0 # * dlink ## Can probably remove dlink
      # std.dg.refvals <- dg.refvals
      # inequality constraints ensuring the
      # (weak or strict) monotonicity
      grads.mono.par <- NULL
      if(constr){
        # set the reference point to be the detection function
        # value at 0
        std.dg.prev <- std.dg.0
        grads.mono.par <- double()
        for(i in 1:no.points){
          grads.mono.par[i] <- (std.dg.prev - std.dg.refvals[i])
          if(strict){
            # if we have strict monotonicity, then change the ref
            # point to be the last point
            std.dg.prev <- std.dg.refvals[i]
          }
        }
      }
      
      # inequality constraints ensuring that the detection function is 
      # always >=0
      # grads.pos.par <- double(no.points) # FTP: why create ic_p and then overwrite it? Commented out for now
      grads.pos.par <- std.dg.refvals
      
      #  set of inequality constraints
      return(c(grads.mono.par, grads.pos.par))
    })
    # })
    
    # ## Below is version that uses loops and is way slower. 
    # grad.constraints <- matrix(NA, nrow = no.points * 2, ncol = no.pars)
    # 
    # ## Now loop through the parameters
    # for (par.index in 1:no.pars) {
    #   ## To get detfct to play nice need to mudge ddfobj a bit...
    #   ## MADE REDUNDANT BY USING select
    #   # if(!is.null(ddfobj$scale)){
    #   #   ddfobj$scale$dm <- 1
    #   # }
    #   # if(!is.null(ddfobj$shape)){
    #   #   ddfobj$shape$dm <- 1
    #   # }
    #   # 
    #   ## Evaluate the gradient of the detection function at the reference points
    #   ## Note that we must standardize so 0<=g(x)<=1. this is done at line 124
    #   dg.0 <- distpdf.grad(distance = ref.point0, par.index = par.index, 
    #                        ddfobj = ddfobj, standardize = FALSE, 
    #                        left = left, width = width, 
    #                        select = c(TRUE, rep(FALSE, no.data - 1)))
    #   
    #   # std.dg.0 <- (dg.0 - g.0 * dg.0 / g.0) / g.0 ## Always 0
    #   std.dg.0 <- 0
    #   
    #   ## To get detfct to play nice need to mudge ddfobj a bit... 
    #   ## MADE REDUNDANT BY USING select
    #   # if(!is.null(ddfobj$scale)){
    #   #   ddfobj$scale$dm <- rep(1, no.points)
    #   # }
    #   # if(!is.null(ddfobj$shape)){
    #   #   ddfobj$shape$dm <- rep(1, no.points)
    #   # }
    #   dg.refvals <- distpdf.grad(distance = ref.points, par.index = par.index, 
    #                              ddfobj = ddfobj, standardize = FALSE, 
    #                              left = left, width = width,
    #                              select = c(rep(TRUE, no.points), 
    #                                         rep(FALSE, no.data - no.points)))
    #   
    #   ## Derive the gradients of the standardised detection function [see notes 17-06-2024]
    #   std.dg.refvals <- (dg.refvals - g.refvals * dg.0 / g.0) / g.0 # * dlink ## Can probably remove dlink
    #   
    #   # grad.ref.points <- as.vector(detfct.grad(distance = ref.points, 
    #   #                                          par.index = par.index,
    #   #                                          ddfobj = ddfobj, 
    #   #                                          width = misc.options$width,
    #   #                                          standardize = TRUE))
    #   
    #   # ## Evaluate the gradient of the detection function at 0
    #   # grad.ref.point0 <- as.vector(detfct.grad(distance = ref.point0,
    #   #                                          par.index = par.index,
    #   #                                          ddfobj = ddfobj, 
    #   #                                          width = misc.options$width,
    #   #                                          standardize = TRUE))
    #   
    #   # inequality constraints ensuring the
    #   # (weak or strict) monotonicity
    #   grads.mono.par <- NULL
    #   if(constr){
    #     # set the reference point to be the detection function
    #     # value at 0
    #     std.dg.prev <- std.dg.0
    #     grads.mono.par <- double()
    #     for(i in 1:no.points){
    #       grads.mono.par[i] <- (std.dg.prev - std.dg.refvals[i])
    #       if(strict){
    #         # if we have strict monotonicity, then change the ref
    #         # point to be the last point
    #         std.dg.prev <- std.dg.refvals[i]
    #       }
    #     }
    #   }
    #   
    #   # inequality constraints ensuring that the detection function is 
    #   # always >=0
    #   # grads.pos.par <- double(no.points) # FTP: why create ic_p and then overwrite it? Commented out for now
    #   grads.pos.par <- std.dg.refvals
    #   
    #   #  set of inequality constraints
    #   grad.constraints[, par.index] <- c(grads.mono.par, grads.pos.par)
    # }
  }
  return(-1 * grad.constraints)
}

## Negative of the constraint gradient
grad.constr <- function(pars, ddfobj, misc.options, fitting = "all") {
  return(-1 * grad.constr.neg(pars = pars, 
                              ddfobj = ddfobj, 
                              misc.options = misc.options))
}
