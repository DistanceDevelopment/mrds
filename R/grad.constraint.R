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
flnl.constr <- function(pars, ddfobj, misc.options, fitting) {
  
  if (is.null(ddfobj$adjustment)) {
    # This never gets called from ddf()
    ineq_constr <- rep(10, 2 * misc.options$mono.points)
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
    
    no.pars <- length(pars)
    # par.indices <- getpar(ddfobj, index = TRUE)
    # k <- sum(par.indices[-3])
    # m <- par.indices[3] - k
    # key.scale <- ifelse(par.indices[2] != 0, pars[par.indices[2]], NULL)
    # key.shape <- ifelse(par.indices[1] != 0, pars[par.indices[1]], NULL)
    
    ## Create matrix to start the constraint values, which are
    ## no.points values for the monotonicity constraint and no.points values
    ## for the positivity constraint, for no.pars parameters.
    grad.constraints <- matrix(NA, nrow = no.points * 2, ncol = no.pars)

    for (par.index in 1:no.pars) {
      ## To get detfct to play nice need to mudge ddfobj a bit...
      if(!is.null(ddfobj$scale)){
        ddfobj$scale$dm <- rep(1, no.points)
      }
      if(!is.null(ddfobj$shape)){
        ddfobj$shape$dm <- rep(1, no.points)
      }
      
      ## Evaluate the gradient of the detection function at the reference points
      ## Note that we must standardize so 0<=g(x)<=1
      grad.ref.points <- as.vector(detfct.grad(distance = ref.points, 
                                               par.index = par.index,
                                               ddfobj = ddfobj, 
                                               width = misc.options$width,
                                               standardize = TRUE))
      
      
      ## Again, to get detfct to play nice need to mudge ddfobj a bit...
      if(!is.null(ddfobj$scale)){
        ddfobj$scale$dm <- 1
      }
      if(!is.null(ddfobj$shape)){
        ddfobj$shape$dm <- 1
      }
      
      ## Evaluate the gradient of the detection function at 0
      grad.ref.point0 <- as.vector(detfct.grad(distance = ref.point0,
                                               par.index = par.index,
                                               ddfobj = ddfobj, 
                                               width = misc.options$width,
                                               standardize = TRUE))
      

      ## Gradients of the constraints w.r.t. current parameter
      ## Version de Felix
      ## =====================================================
      ## Gradients of inequality constraints at the reference points
      if (constr) {
        if (strict) {
          grads.mono.par <- grad.ref.points - c(grad.ref.point0, 
                                           grad.ref.points[-no.points])
        } else {
          grads.mono.par <- grad.ref.points - grad.ref.point0
        }
      }
      
      ## Gradients of positivity constraint at the reference points
      # grads.pos.par <- double(no.points) # FTP: why create ic_p and then overwrite it?
      grads.pos.par <- grad.ref.points
      
      #  set of inequality constraints
      grad.constraints[, par.index] <- c(grads.mono.par, grads.pos.par)
    }
  }
  return(grad.constraints)
}