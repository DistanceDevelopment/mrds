#' Extract coefficients
#'
#' Extract coefficients and provide a summary of parameters and estimates
#' from the output of \code{\link{ddf}} model objects.
#'
#' @aliases coefficients coef.ds coef.io coef.io.fi coef.trial coef.trial.fi
#'   coef.rem coef.rem.fi
#' @param object ddf model object of class \code{ds}, \code{io}, \code{io.fi},
#'    \code{trial}, \code{trial.fi}, \code{rem}, or \code{rem.fi}.
#' @param \dots unspecified arguments that are unused at present
#' @usage \method{coef}{ds}(object,...)
#'        \method{coef}{io}(object,...)
#'        \method{coef}{io.fi}(object,...)
#'        \method{coef}{trial}(object,...)
#'        \method{coef}{trial.fi}(object,...)
#'        \method{coef}{rem}(object,...)
#'        \method{coef}{rem.fi}(object,...)
#' @S3method coef ds
#' @S3method coef io
#' @S3method coef io.fi
#' @S3method coef trial
#' @S3method coef trial.fi
#' @S3method coef rem
#' @S3method coef rem.fi
#' @return
#'
#' For \code{coef.ds} List of data frames for coefficients (scale and exponent
#'   (if hazard)) \item{scale}{dataframe of scale coefficent estimates and
#'   standard errors} \item{exponent}{dataframe with exponent estimate and
#'   standard error if hazard detection function}
#'
#' For all others Data frame containing each coefficient and standard error
#' @note These functions are called by the generic function \code{coef} for any
#'   \code{ddf} model object.  It can be called directly by the user, but it is
#'   typically safest to use \code{coef} which calls the appropriate function
#'   based on the type of model.
#' @author Jeff Laake;
coef.ds <- function(object,...){
  ltmodel <- object$ds
  coeff <- NULL
  if(is.null(object$par)){
    vcov <- NULL
  }else{
    vcov <- solvecov(object$hessian)$inv
  }

  ddfobj <- ltmodel$aux$ddfobj
  se <- sqrt(diag(vcov))

  # make a fake ddf object to easily extract the SEs
  se.dd <- assign.par(ddfobj,se)
  key.shape.se <- se.dd$pars$shape
  key.scale.se <- se.dd$pars$scale
  adj.parm.se <- se.dd$pars$adjustment

  # Always get the scale parameter
  if(!is.null(ddfobj$scale)){
    coeff <- as.data.frame(cbind(estimate=ddfobj$pars$scale,
                                 se=key.scale.se))
    row.names(coeff) <- colnames(ddfobj$scale$dm)
  }

  # If we have a non-null shape parameter, get that too
  if(!is.null(ddfobj$shape)){
    exp.coeff <- as.data.frame(cbind(estimate=ddfobj$pars$shape,
                                     se=key.shape.se))
    row.names(exp.coeff) <- colnames(ddfobj$shape$dm)
  }

  # Get the adjustment parameters if necessary
  if(!is.null(ddfobj$adjustment)){
    adj.coeff <- as.data.frame(cbind(estimate=ddfobj$pars$adjustment,
                                     se=adj.parm.se))
    adj.names <- NULL
    for(i in 1:nrow(adj.coeff)){
      adj.names[i] <- paste(ddfobj$adjustment$series,
                            ", order ",ddfobj$adjustment$order[i],sep="")
    }
    row.names(adj.coeff) <- adj.names

    # Return values if we have adjustment terms
    if(!is.null(ddfobj$shape)){
      return(list(scale=coeff,exponent=exp.coeff,adjustment=adj.coeff))
    }else if(is.null(ddfobj$scale)){
      return(list(adjustment=adj.coeff))
    }else{
      return(list(scale=coeff,adjustment=adj.coeff))
    }

  }

  # Return values if no adjustment terms
  if(!is.null(ddfobj$shape)){
    return(list(scale=coeff,exponent=exp.coeff))
  }else{
    return(list(scale=coeff))
  }
}
