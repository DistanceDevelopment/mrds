#' log-likelihood value for a fitted detection function
#'
#' Extract the log-likelihood from a fitted detection function.
#'
#' @param object a fitted detection function model object
#' @param \dots included for S3 completeness, but ignored
#' @return a numeric value giving the log-likelihood with two attributes: \code{"df"} the "degrees of freedom" for the model (number of parameters) and \code{"nobs"} the number of observations used to fit the model
#' @export
#' @author David L Miller
#' @aliases logLik.ds logLik.io logLik.io.fi logLik.rem logLik.rem.fi logLik.trial logLik.trial.fi
logLik.ddf <- function(object, ...){

  # see ?logLik for information on why

  ret <- object$lnl

  attr(ret, "df") <- length(object$par)
  attr(ret, "nobs") <- nrow(object$data)

  class(ret) <- "logLik"
  return(ret)
}
