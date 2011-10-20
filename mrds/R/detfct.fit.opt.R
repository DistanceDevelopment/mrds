#' Fit detection function using key-adjustment functions
#' 
#' Fit detection function to observed distances
#' using the key-adjustment
#' function approach.  If adjustment functions are included it will alternate
#' between fitting parameters of key and adjustment functions and then all
#' parameters much like the approach in the CDS and MCDS Distance FORTRAN code.
#' This function is called by the driver functioin \code{detfct.fit}.  This function does
#' the calls the R optim function
#' which does not allow non-linear constraints so inclusion of adjustments does
#' allow the detection function to be non-monotone.
#' 
#' 
#' @import optimx
#' @aliases detfct.fit.opt
#' @param ddfobj detection function object
#' @param optim.options control options for optim
#' @param bounds bounds for the parameters
#' @param misc.options miscellaneous options
#' @param fitting character string with values "all","key","adjust" to
#'   determine which parameters are allowed to vary in the fitting
#' @return fitted detection function model object with the following list
#'   structure \item{par}{final parameter vector} \item{value}{final negative
#'   log likelihood value} \item{counts}{number of function evaluations}
#'   \item{convergence}{see codes in optim} \item{message}{string about
#'   convergence} \item{hessian}{hessian evaluated at final parameter values}
#'   \item{aux}{ a list with 21 elements \itemize{ \item maxit: maximum number
#'   of iterations allowed for optimization \item lower: lower bound values for
#'   parameters \item upper: upper bound values for parameters \item setlower:
#'   TRUE if they are user set bounds \item setupper: TRUE if they are user set
#'   bounds \item point: TRUE if point counts and FALSE if line transect \item
#'   int.range: integration range values \item showit: integer value that
#'   determines information printed during iteration \item doeachint: if TRUE
#'   each integral is computed numerically rather using tabled values \item
#'   integral.numeric if TRUE compute logistic integrals numerically \item
#'   breaks: breaks in distance for defined fixed bins for analysis \item
#'   maxiter: maximum iterations used \item refit: if TRUE, detection function
#'   will be fitted more than once if parameters are at a boundary or when
#'   convergence is not achieved \item nrefits: number of refittings \item
#'   parscale: parameter scale values \item fdebug: not sure what this means
#'   \item nonmono: logical as to whether montonicity should be enforced \item
#'   strict: if TRUE, then strict monotonicity is enforced; otherwise weak
#'   \item width: radius of point count or half-width of strip \item
#'   standardize: if TRUE, detection function is scaled so g(0)=1 \item ddfobj:
#'   distance detection function object; see \code{\link{create.ddfobj}} \item
#'   bounded: TRUE if parameters ended up a boundary (I think) \item model:
#'   list of formulas for detection function model (probably can remove this)
#'   }}
#' @author Dave Miller; Jeff Laake
detfct.fit.opt <- function(ddfobj,optim.options,bounds,misc.options,fitting="all")
{
#
# detfct.fit.opt
#
# A function to actually do the optimisation.
#
# Arguments:
#
#  ddfobj      		    - distance sampling object for the detection function
#  control.options 		- control options for optim
#  bounds				- bounds for the function
#  misc.options			- things that wouldn't fit in elsewhere
#  fitting				- "all","key","adjust"
#
# Value:
#
#  lt				- lt object (something more descriptive)
#
# Functions Used:
#
#  assign.par, detfct.fit.opt, errors, get.par 
	
  initialvalues=getpar(ddfobj)
# dlm 31-May-2006 Initial work started.
# dlm 05-June-2006 Added an extra level of showit for this level.
  bounded<-TRUE
  refit.count<-0
# Set some shortcuts
  lowerbounds<-bounds$lower
  upperbounds<-bounds$upper
  showit<-misc.options$showit
  refit<-misc.options$refit
  nrefits<-misc.options$nrefits
#
# jll 18-sept-2006; added this code to get the logicals that indicate whether
# lower/upper bound settings were specified by the user
  setlower<-bounds$setlower
  setupper<-bounds$setupper
#
# 30 Jan 06; jll- modified default parscale values
  if(any(is.na(misc.options$parscale)))
     misc.options$parscale<-abs(initialvalues)
  else{
     if(length(misc.options$parscale)!=length(initialvalues)){
        errors("Incorrect length of parameter scale vector; using default values\n")
        misc.options$parscale<-abs(initialvalues)
     }
  }


#    while parameters are bounded continue refitting and adjusting bounds
  while(bounded){
    itconverged<-FALSE

#    Continue fitting until convergence occurs as long as refitting is requested
    while(!itconverged){
#
#    Call optimization routine to find constrained mles; upon completion add the user specified 
#    models and return the list.
#
#    change 17-Aug-05 added parscale and maxit controls in call to optim
#                     added conditions here and below that control printing based on showit

      lt <- optimx(initialvalues, flnl, method="nlminb", control=c(optim.options),
		hessian=TRUE, lower = lowerbounds, upper = upperbounds,  
		ddfobj=ddfobj, fitting=fitting,misc.options=misc.options,TCI= FALSE)
      lt <-attr(lt,"details")[[1]]
	  lt$hessian<-lt$nhatend
#      lt$shapemodel=misc.options$shapemodel
      if(showit==3)
        errors(paste("Converge = ",lt$conv,"\nlnl = ",lt$value,"\nFinal values = ",paste(lt$par,collapse=", ")))
# If we do have convergence what do we do
      if(lt$conv==0|!refit){
        itconverged<-TRUE 

        lt$aux <- c(optim.options,bounds,misc.options)
		ddfobj=assign.par(ddfobj,lt$par)
		lt$aux$ddfobj=ddfobj
      }else{
# If we don't have converenge what do we do
        refit.count<-refit.count+1
        if(is.null(nrefits)|refit.count<=nrefits){
         
          if(showit>=1)
            errors("No convergence. Refitting ...")
          
          if(lt$conv==1)
            initialvalues<-lt$par
          else
            initialvalues<-lt$par*(runif(length(initialvalues))+.5)
		    initialvalues[is.na(initialvalues)]=0
        }else{
          itconverged<-TRUE
        }
      }
    }   
#
#  Issue warning if any of the parameters are at their bounds
    bounded <- FALSE
	if(any(is.na(lt$par)))
	{
		errors("Problems with fitting data. Did not converge")
		stop()
	}
# fix jll 17-Aug-05 allows for constrained power parameter in hazard rate
    if(ddfobj$type=="hr"){
      if(any(abs(lt$par[2:length(lt$par)]-lowerbounds[2:length(lt$par)])<0.000001)){
        if(!setlower) 
          bounded<-TRUE
  
        if(showit>=1)
          errors(paste("One or more parameters was at a lower bound.\nParameters:",paste(lt$par,collapse=", "),
                       "\nLower bounds:",paste(lowerbounds,collapse=", ")))
        
      }
    }else{
      if(any(abs(lt$par-lowerbounds)<0.000001)){
        if(!setlower)
          bounded<-TRUE
        
        if(showit>=1)
          errors(paste("One or more parameters was at a lower bound\nParameters: ",
						  paste(lt$par,collapse=", "),"\nLower bounds: ", paste(lowerbounds,collapse=", ")))

      }
    }
    if(any(abs(lt$par-upperbounds)<0.000001)){ 
      if(!setupper)
        bounded<-TRUE

      if(showit>=1)
        errors(paste("One or more parameters was at an upper bound\nParameters: ",
						paste(lt$par,collapse=", "),"\nUpper bounds: ", paste(upperbounds,collapse=", ")))
      
    }
    lt$bounded <- bounded
    if(!refit)
        bounded<-FALSE
    else{
        refit.count=refit.count+1
        if(!is.null(nrefits))
          if(refit.count>nrefits){ 
            bounded<- FALSE
            lt$message<-"FALSE CONVERGENCE"
            lt$conv<-1
        }
    }
#
# fix: jll 18-Nov-04; previous code would get stuck at 0 if sign of initial value was opposite
#      of the mle.  Additional statement skips over 0.
#
    if(bounded){
        bound.low<-abs(lt$par-lowerbounds)<0.000001
        bound.hi<-abs(lt$par-upperbounds)<0.000001
        if(!setlower) {
          lowerbounds[bound.low] <- lowerbounds[bound.low] - 0.5*abs(lowerbounds[bound.low])
          lowerbounds[bound.low & lowerbounds>0 &lowerbounds < 0.5]=-0.5
        }
        if(!setupper){
          upperbounds[bound.hi] <- upperbounds[bound.hi] + 0.5*abs(upperbounds[bound.hi])
          upperbounds[bound.hi & upperbounds<0 &upperbounds > -0.5]= 0.5
        }
      
        if(showit>=1)
          errors("Refitting ...")

# change jll 17-Aug-05 added to constrain haz exponent to be > 1
#        if(ddfobj$type=="hr")
#          lowerbounds[1]=1
      
    }
  }
  lt$model <- list(scalemodel=misc.options$scalemodel) 
  lt$converge<-lt$conv
  lt$conv<-NULL

  return(lt)


}

