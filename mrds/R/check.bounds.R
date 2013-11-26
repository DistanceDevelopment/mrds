#' Check bounds during optimsations
#'
#' Simple internal function to check that the optimisation didn't hit bounds.
#' Based on code that used to live in \code{detfct.fit.opt}.
#'
#' @param lt optimisation object
#' @param lowerbounds current lower bounds
#' @param upperbounds current upper bounds
#' @param ddfobj ddf object
#' @param showit debug level
#' @param setlower were lower bounds set by the user
#' @param setupper were upper bounds set by the user
#'
#' @author Dave Miller; Jeff Laake
check.bounds <- function(lt,lowerbounds,upperbounds,ddfobj,showit,
                         setlower,setupper){

  # function to check upper/lower bounds
  chk.bnds <- function(par,bounds,bound.label,set){
    if(any(abs(par-bounds)<0.000001)){
      if(showit>=1){
        # Issue warning if any of the parameters are at their bounds
        message(paste("One or more parameters was at a",bound.label,"bound\n",
                      "Parameters:",paste(par,collapse=", "),"\n",
                      bound.label,"bounds:",paste(bounds,collapse=", ")))
      }
      if(!set){
        return(TRUE)
      }
    }
    return(FALSE)
  }

  ## check lower bounds
  # handle hazard rate power par
  if(ddfobj$type=="hr"){
    bounded <- chk.bnds(lt$par[2:length(lt$par)],
                        lowerbounds[2:length(lt$par)],"lower",setlower)
  }else{
    bounded <- chk.bnds(lt$par,lowerbounds,"lower",setlower)
  }

  ## check upper bounds
  bounded <- bounded | chk.bnds(lt$par,upperbounds,"upper",setupper)

  return(bounded)
}
