#' Alternate detection function fitting code
#' 
#' A wrapper function for foptim.dll -- the Fortran optimisation routine. Most
#' of the code juggles the data before making the call to ensure Fortran
#' doesn't fail.
#' 
#' foptim.dll is a re-implementation of the CDS engine found in Distance. It
#' uses the DNCONG/DNCONF routines from the IMSL Fortran library in order to
#' perform constrained likelihood optimisation. The likelihood is maximised
#' ensuring that the detection function remains monotonic. This is done by
#' evaluating the detection function at 20 equally-spaced distances, ensuring
#' that the probability of detection is always decreasing with increasing
#' distance.
#' 
#' NB: this function is ONLY available under Windows. Covariates are not
#' handled at the moment.
#' 
#' This routine should only be called from higher level functions. YMMV.
#' 
#' 
#' @aliases foptim foptim.dll
#' @param initialvalues vector of initial values of the parameters
#' @param bounds list with two elements, "upper" and "lower" containing the
#'   upper and lower bounds for parameters.
#' @param detfct.options options defining detection function; list with
#'   elements: ftype - key function, "unif", "hn" or "hr"
#' 
#' adj.series - adjustment terms to use: "NULL", "poly", "cos" or "herm"
#' 
#' adj.order - order of adjustment
#' 
#' adj.scale - scaling used for the adjustment terms (scale or width)
#' 
#' width - truncation distance
#' @param misc.options miscellaneous options; list with elements: showit - as
#'   in other functions, shows debugging information
#' 
#' nonmono - should monotonicity be enforced (default no)
#' 
#' strict - should strict monotonicity (default no)
#' 
#' maxiter - maximum number of iterations for the Fortran code (default 50)
#' 
#' x - data frame
#' @return fitted model object \item{par}{parameter estimates}
#'   \item{fvalue}{value of the log-likelihood at the parameter estimates}
#'   \item{message}{if the procedure has converged the message "CONVERGENCE"}
#'   \item{aux}{list of other auxiliary information: z - covariate matrix (not
#'   used)
#' 
#' zdim - dimension of z (not used)
#' 
#' itercept.only - (not used)
#' 
#' int.range - integration range (not used)
#' 
#' width - truncation distance
#' 
#' x - data frame
#' 
#' ftype - key function type
#' 
#' adj.series - adjustment series type ("NULL" for none)
#' 
#' adj.order - adjustment order
#' 
#' adj.scale - scaling for adjustments}
#' 
#' \item{misc.options}{contents of misc.options, above} \item{bounded}{TRUE if
#'   the parameters hit their upper or lower bounds}
#' @author David L. Miller \email{mrds@@ninepointeightone.net}
foptim<-function(initialvalues,bounds,detfct.options,misc.options)
{
#
# detfct.fit.opt is called with: initialvalues,optim.options,bounds,
#				detfct.options,misc.options,opt.options 
# foptim
#   A wrapper function for foptim.dll. It also does a little data
#   juggling to make sure Fortran doesn't barf.

# foptim is called in Fortran as follows:
# Note that all of the arrays are horrible sizes due to COMMON.
# Arguments:
#   n			number of data points of the above.
#   x			vector distances 
#   xguess		first guess at parameter values.
#   nfpar		number of parameters to estimate
#			NB. This is worked out below, don't panic.
#   parub		Upper bound on the parameters.
#   parlb		Lower bound on the parameters.
#   maxiter		Maximum number of iterations to go through
#   maxnlcs		
#   stricts
#   nonmonos
#   detkey		Key function to use for the detection function.
#			 0 - uniform
#			 1 - half-normal
#			 2 - hazard-rate
#   detadj	
#			 0 - none
#			 1 - simple polynomial
#			 2 - cosine
#			 3 - Hermite polynomial
#   detadjorder
#   detnadjorder
#   truncwidth
#   debuglevel		IMSL DNCONG debug level
#			  0 - No output
#			  1 - Only final convergence analysis
#			  2 - One line of output per iter
#			  3 - Detailed information
# Returns:
#   fvalue		ln(likelhood) value for the solution
#   fpar		parameter values at solution point
#
# dlm 09-June-2006 Initial work started.
# dlm 03-July-2007 Started doing Real Work(TM)
#  NAG website is useful: http://www.nag.co.uk/numeric/RunderWindows.asp
#  Acadia U also has useful information: 
#    http://ace.acadiau.ca/math/ACMMaC/howtos/Fortran_R.html


# Check that we can load the dll before we do anything

# This needs to be done properly
  library.dynam("foptim2",package=c("mrds"))

# Set the data
# Things not passed (yet)
# covariates
#  detfct.options$zdim
# structure of misc.options
#  z=z, int.range=int.range,showit=showit,
#  intercept.only=intercept.only,cgftab=cgftab,
#  doeachint=doeachint,scalemodel=scalemodel,
#  integral.numeric=integral.numeric, breaks=breaks,
#  refit=refit,nrefits=nrefits,
#  parscale=parscale)

  showit<-misc.options$showit

  # Convert the key and adjustments to integer codes
  if(is.null(detfct.options$ftype)){
    detkey<-0
  }else{
    detkey<-switch(detfct.options$ftype,unif=0,hn=1,hr=2)
  }

  if(is.null(detfct.options$adj.series)){
    detadj<-0
  }else{
    detadj<-switch(detfct.options$adj.series,poly=1,cos=2,herm=3)
  }

  # Need to make sure that if we use hazard rate, we have the parameters
  # the right way around: scale,shape not shape,scale
  if(detfct.options$ftype=="hr"){
    tmp<-initialvalues[1]
    initialvalues[1]<-initialvalues[2]
    initialvalues[2]<-tmp
    parub<-c(bounds$upper[2],bounds$upper[1],bounds$upper[-c(1,2)])
    parlb<-c(bounds$lower[2],bounds$lower[1],bounds$lower[-c(1,2)])
    rm(tmp)
  }else{
    parub<-bounds$upper
    parlb<-bounds$lower
  }

# jll 18-sept-2006; added this code to get the logicals that indicate whether
# lower/upper bound settings were specified by the user
  setlower<-bounds$setlower
  setupper<-bounds$setupper
  
  # Just some hacky bits here
  x<-misc.options$x
  if(is.null(detfct.options$adj.order)){
    adj.order<-0
  }else{
    adj.order<-detfct.options$adj.order
  }

  # By default no strict monotonicity
  if(is.null(misc.options$strict)){
    misc.options$strict<-0
  }else{
    if(misc.options$strict){
      misc.options$strict<-1
    }else{
      misc.options$strict<-0
    }
  }

  # By default no monotonicity
  if(is.null(misc.options$nonmono)){
    misc.options$nonmono<-0
  }else{
    if(misc.options$nonmono){
      misc.options$nonmono<-1
    }else{
      # If we already selected strict monotonicity then choose weak
      if(misc.options$strict){
	misc.options$nonmono<-1
      }else{
        misc.options$nonmono<-0
      }
    }
  }

  # set the scaling
  # 0 = key scale parameter
  # 1 = width
  if(is.null(detfct.options$adj.scale)){
    adj.scaling<-0
  }else{
    adj.scaling<-switch(detfct.options$adj.scale,scale=0,width=1)
  }

  # Set some default values
  if(is.null(misc.options$fdebug))
    misc.options$fdebug<-0

  # set the default number of iterations for the Fortran code
  # to be 50. It's usually 12, this is too few.
  if(misc.options$maxiter<50)
    misc.options$maxiter<-50

  #    while parameters are bounded continue refitting and adjusting bounds
  bounded<-TRUE
  while(bounded){

    # do the call to Fortran
    foptim_call <- .Fortran("foptim",
		     vn=length(x$distance),
		     vx=as.double(x$distance),
	 	     xguess=as.double(initialvalues),
		     nfpar=as.integer(length(initialvalues)),
		     parub=as.double(parub),
  		   parlb=as.double(parlb),
	  	   maxiter=as.integer(misc.options$maxiter),
	 	     stricts=as.integer(misc.options$strict),
    	   nonmonos=as.integer(misc.options$nonmono),
		     vkey=as.integer(detkey),
		     vadj=as.integer(detadj),
		     vadjorder=as.integer(adj.order),
		     vnadjorder=as.integer(length(adj.order)),
		     vadjscaling=as.integer(adj.scaling),
		     vwidth=as.double(detfct.options$width),
		     vshowit=as.integer(misc.options$showit),
		     fvalue=double(1),
		     fpar=double(length(initialvalues)), PACKAGE="mrds"
       )
               
    #  Issue warning if any of the parameters are at their bounds
    bounded <- FALSE
    # allows for constrained power parameter in hazard rate
    if(detfct.options$ftype=="hr"){
      if(any(abs(foptim_call$fpar[2:length(foptim_call$fpar)]-parlb[2:length(foptim_call$fpar)])<0.000001)){
        if(!setlower) 
          bounded<-TRUE
  
        if(showit>=1)
          errors(paste("One or more parameters was at a lower bound.\nParameters:",foptim_call$fpar,
                       "\nLower bounds:",parlb))
        
      }
    }else{
      if(any(abs(foptim_call$fpar-parlb)<0.000001)){
        if(!setlower)
          bounded<-TRUE
        
        if(showit>=1)
          errors(paste("One or more parameters was at a lower bound\nParameters: ",
                       foptim_call$fpar,"\nLower bounds: ", parlb))

      }
    }
    if(any(abs(foptim_call$fpar-parub)<0.000001)){ 
      if(!setupper)
        bounded<-TRUE

      if(showit>=1)
        errors(paste("One or more parameters was at an upper bound\nParameters: ",
                     foptim_call$fpar,"\nUpper bounds: ", parub))
      
    }

    if(bounded){
      bound.low<-abs(foptim_call$fpar-parlb)<0.000001
      bound.hi<-abs(foptim_call$fpar-parub)<0.000001
      if(!setlower) {
        parlb[bound.low] <- parlb[bound.low] - 0.5*abs(parlb[bound.low])
        parlb[bound.low & parlb>0 & parlb < 0.5]=-0.5
      }
      if(!setupper){
        parub[bound.hi] <- parub[bound.hi] + 0.5*abs(parub[bound.hi])
        parub[bound.hi & parub<0 &parub > -0.5]= 0.5
      }
      
      if(showit>=1)
        errors("Refitting ...")

      # change jll 17-Aug-05 added to constrain haz exponent to be > 1
      if(detfct.options$ftype=="hr")
        parlb[1]=1
    }
  }

  # Format a nice object to send back!
  res<-list()

  # Parameters
  res$par<-foptim_call$fpar
  # log-likelihood value
  res$value<-foptim_call$fvalue

  # just tell mrds that everything is fine
  res$message<-"CONVERGENCE"
 

  # set a whole bunch of variables 
  res$aux$z<-misc.options$z
  res$aux$zdim<-detfct.options$zdim
  res$aux$intercept.only<-misc.options$intercept.only
  res$aux$int.range<-misc.options$int.range
  res$aux$width<-detfct.options$width
  res$aux$x<-misc.options$x
  res$aux$ftype<-detfct.options$ftype
  res$aux$adj.series<-detfct.options$adj.series
  res$aux$adj.order<-detfct.options$adj.order
  res$aux$adj.scale<-detfct.options$adj.scale

  # shove in the misc.options that got sent
  res$misc.options<-misc.options

  res$bounded<-bounded

#  showit<-showit
#  ,cgftab=cgftab,
#  doeachint=doeachint,scalemodel=scalemodel,
#  integral.numeric=integral.numeric, breaks=breaks,
#  refit=refit,nrefits=nrefits,
#  parscale=parscale)

  # unload the fortran module
  library.dynam.unload("foptim2",paste(.libPaths()[1],"/mrds",sep="")) 

  return(res)

}
