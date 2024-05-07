#' Wrapper for the auglag function from nloptr to overcome a scoping issue
#'
#' Some description.
#'
#' @param x0 the starting values
#' @param fn the objective function
#' @param ddfobj detection function object
#' @param misc.options miscellaneous options
#' @param ... Add additional arguments for \code{\link{nloptr::auglag}}
#' @return 
#' @author Dave Miller; Jeff Laake; Lorenzo Milazzo
#' @importFrom nloptr auglag
auglagWrapper <- function (x0, 
                           fn, 
                           ddfobj,
                           misc.options,
                           ...) {
  ## FTP: When trying to run nloptr::auglag(), I get the following error:
  # Error in .hin(x) : argument "ddfobj" is missing, with no default
  # I think this is because of the following lines in nloptr::auglag()
  #           hin <- function(x) (-1) * .hin(x)
  # What happens here is that only x is passed on, but not ... . Therefore, the 
  # additional information that flnl.contr requires, such as ddfobj and 
  # misc.options, are not passed to the hin() function, 
  # and thus it cannot calculate the constraint.
  # Cal and I came up with the following solution:
  # Write a wrapper function for the optimiser which takes ddfobj ect as 
  # arguments, and within that wrapper define the constraint function 
  # flnl.contr(pars) the same as in detfct.fit.mono.R but only give it one 
  # argument. As it also needs ddfobj and misc.options, but cannot find it, it 
  # will look up one higher scope, which is the wrapper function, in which it 
  # will be able to find these variables. Now the optimiser is happy as it only 
  # needs to pass the current parameters as an argument to the constraint 
  # function, and the constraint function is happy as it can find the other 
  # necessary arguments by looking in a higher scope. 
  
  ## Add the ... to the enviroement -- not necessary if naming is correct
  # list2env(list(...), envir = environment())    # add these to the environment

  # 1. Define the constraint function similarly to detfct.fit.mono.R
  # It will only get parameters as an argument and look elsewhere for other 
  # information. 
  #
  # Input:
  #  pars           - parameters
  #
  constr <- function(pars){
    
    if(is.null(ddfobj$adjustment)){
      # this never gets called from ddf()
      ineq_constr <- rep(10,2*misc.options$mono.points)
    }else{
      ddfobj <- assign.par(ddfobj, pars)
      
      # apply the constraints?
      constr <- misc.options$mono
      # apply strict monotonicity?
      strict <- misc.options$mono.strict
      
      ### Constraint stuff here:
      # number of points (distances)
      # at which the DF is evaluated
      no_d <- misc.options$mono.points
      # reference points (distances)
      ref_p <- getRefPoints(no_d, misc.options$int.range)
      # to get detfct to play nice need to mudge ddfobj a bit...
      if(!is.null(ddfobj$scale)){
        ddfobj$scale$dm <- rep(1, no_d)
      }
      if(!is.null(ddfobj$shape)){
        ddfobj$shape$dm <- rep(1, no_d)
      }
      
      # evaluate the detection function at the reference points
      # note that we must standardize so 0<=g(x)<=1
      df_v_rp <- as.vector(detfct(ref_p,ddfobj,width=misc.options$width,
                                  standardize=TRUE))
      
      # reference point associated with distance=0
      ref_p0 <- 0
      # again, to get detfct to play nice need to mudge ddfobj a bit...
      if(!is.null(ddfobj$scale)){
        ddfobj$scale$dm <- 1
      }
      if(!is.null(ddfobj$shape)){
        ddfobj$shape$dm <- 1
      }
      
      # evaluate the detection function at 0
      df_v_rp0 <- as.vector(detfct(ref_p0,ddfobj,width=misc.options$width,
                                   standardize=TRUE))
      
      # START OLD CODE >>>>>
      # # inequality constraints ensuring the
      # # (weak or strict) monotonicity
      # ic_m <- NULL
      # if(constr){
      #   # set the reference point to be the detection function
      #   # value at 0
      #   df_v_rp_p <- df_v_rp0
      #   ic_m <- double(no_d)
      #   for(i in 1:no_d){
      #     ic_m[i] <- (df_v_rp_p - df_v_rp[i])
      #     if(strict){
      #       # if we have strict monotonicity, then change the ref
      #       # point to be the last point
      #       df_v_rp_p <- df_v_rp[i]
      #     }
      #   }
      # }
      # <<<<< END OLD CODE
      # FTP: I think, rather than what is done above, simply create extract
      # the differences with the point before. Then, when enforcing weak 
      # monotonicity, simple keep the numbers as they are, as the inequality 
      # constraint checks for >=0, not only >0; when enforcing strong 
      # monononicity, subtract a very small number (e., 1e-6) from the differences
      # (ic_m) to make the zero differences negative, thereby failing the 
      # monotonicity constraint. 
      # Find the (in my opinion) correct code below:
      # START NEW CODE >>>>>>
      ic_m <- NULL
      if(constr){
        # set the reference point to be the detection function
        # value at 0
        df_v_rp_p <- df_v_rp0
        ic_m <- double(no_d)
        for(i in 1:no_d){
          ic_m[i] <- (df_v_rp_p - df_v_rp[i])
          # update the reference point to the most recently evaluated point
          df_v_rp_p <- df_v_rp[i]
        }
      }
      if (!strict) {
        # subtract small number to make the zero-differences negative
        ic_m <- ic_m - 1e-6
      }
      # <<<<< END NEW CODE
      
      # inequality constraints ensuring that
      # the detection function is always >=0
      # ic_p <- double(no_d) # FTP: why create ic_p and then overwrite it? Commented it out
      ic_p <- df_v_rp
      
      #  set of inequality constraints
      ineq_constr <- c(ic_m, ic_p)
    }
    return(ineq_constr)
  }
  
  # 2. Now run the auglag solver
  auglag(x0 = x0, 
         fn = fn, 
         hin = constr, 
         # ineqLB = lowerbounds.ic, ineqUB = upperbounds.ic,
         # lower = lowerbounds, upper = upperbounds,
         ddfobj = ddfobj, misc.options = misc.options,
         control = list(xtol_rel = misc.options$mono.tol), ...)
}
