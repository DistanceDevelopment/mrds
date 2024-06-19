#' Gradient of constrained function
#'
#' The function consists of the gradients of the constrained function for 
#' all model parameters, in the following order:
#' 1. Scale parameter (if part of key function)
#' 2. Shape parameter (if part of key function)
#' 3. Adjustment parameter 1
#' 4. Adjustment parameter 2
#' 5. Etc. 
#' 
#' The constrained function itself is formed of between 10 and 20 non-linear
#' constraints. The function there return a gradient value for each parameter
#' and each constraint.
#'
#' This function mostly follows the same structure as \code{flnl.constr()} in 
#' \code{detfct.fit.mono.R}. 
#' 
#' @param pars
#' @param ddfobj
#' @param misc.options
#' @param fitting
#' 
#' @value a matrix of gradients for all constraints (rows) w.r.t to all
#' parameters (columns)
grad.constr <- function(pars, ddfobj, misc.options, fitting="all") {
  
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
    
    g.refvals <- distpdf(distance = ref.points, ddfobj = ddfobj,
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
    g.0 <- distpdf(distance = ref.point0, ddfobj = ddfobj, width = width,
                   left = left, point = point, 
                   select = c(TRUE, rep(FALSE, no.data - 1)))
    
    ## Create matrix to start the constraint values, which are
    ## no.points values for the monotonicity constraint and no.points values
    ## for the positivity constraint, for no.pars parameters.
    grad.constraints <- sapply(1:no.pars, function(par.index) {
      ## Evaluate the gradient of the detection function at the reference points
      ## Note that we must standardize so 0<=g(x)<=1. this is done at line 124
      dg.0 <- distpdf.grad(distance = ref.point0, par.index = par.index, 
                           ddfobj = ddfobj, standardize = FALSE, 
                           left = left, width = width) #, 
                           # select = c(TRUE, rep(FALSE, no.data - 1)))
      
      ## Standardised gradient at distance 0 is always 0
      # std.dg.0 <- (dg.0 - g.0 * dg.0 / g.0) / g.0 ## Always 0
      std.dg.0 <- 0

      dg.refvals <- distpdf.grad(distance = ref.points, par.index = par.index, 
                                 ddfobj = ddfobj, standardize = FALSE, 
                                 left = left, width = width) #,
                                 # select = c(rep(TRUE, no.points), 
                                            # rep(FALSE, no.data - no.points)))
      
      ## Derive the gradients of the standardised detection function [see notes 17-06-2024]
      std.dg.refvals <- (dg.refvals - g.refvals * dg.0 / g.0) / g.0 # * dlink ## Can probably remove dlink
      
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
  return(grad.constraints)
}

## Negative of the constraint gradient
grad.constr.neg <- function(pars, ddfobj, misc.options, fitting = "all") {
  return(-grad.constr(pars = pars, 
                      ddfobj = ddfobj, 
                      misc.options = misc.options,
                      fitting = fitting))
}
