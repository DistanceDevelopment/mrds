#' Check parameters bounds during optimsations
#'
#' Simple internal function to check that the optimisation didn't hit bounds. Based on code that used to live in \code{detfct.fit.opt}.
#'
#' @param pars estimated parameters
#' @param lowerbounds current lower bounds
#' @param upperbounds current upper bounds
#' @param ddfobj ddf object
#' @param showit debug level
#' @param setlower were lower bounds set by the user
#' @param setupper were upper bounds set by the user
#' @return \code{TRUE} if bounded (ie parameters close to bound), else \code{FALSE}
#'
#' @author Dave Miller; Jeff Laake
check.bounds <- function(pars, lowerbounds, upperbounds, ddfobj, showit,
                         setlower, setupper){

  tol <- 1e-6

  # function to check upper/lower bounds
  chk.bnds <- function(par, bounds, bound.label, set, tol){

    if(set) return(FALSE)

    if(bound.label=="lower" && (any(par<bounds) | any(abs((bounds-par))<tol))){
      bounded <- TRUE
    }else if(bound.label=="upper" && (any(par>bounds) | any(abs(((bounds-par)<tol))))){
      bounded <- TRUE
    }else{
      bounded <-FALSE
    }

    # Issue message if any of the parameters are at their bounds
    if(bounded & showit>=1){
      message(paste("One or more parameters was at a", bound.label, "bound\n",
                    "Parameters:", paste(par, collapse=", "), "\n",
                    bound.label, "bounds:", paste(bounds, collapse=", ")))
    }
    return(bounded)
  }

  ## check lower bounds
  # handle hazard rate power par
  if(ddfobj$type=="hr"){
    bounded <- chk.bnds(pars[2:length(pars)],
                        lowerbounds[2:length(pars)], "lower", setlower, tol)
  }else{
    bounded <- chk.bnds(pars, lowerbounds, "lower", setlower, tol)
  }

  ## check upper bounds
  bounded <- bounded | chk.bnds(pars, upperbounds, "upper", setupper, tol)

  return(bounded)
}
