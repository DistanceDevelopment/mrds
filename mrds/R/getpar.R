#' Extraction and assignment of parameters to vector
#'
#' Extracts parameters of a particular type (scale,
#' shape, adjustments or g0 (p(0))) from the vector of parameters in
#' \code{ddfobj}. All of the parameters are kept in a single vector for
#' optimization even though they have very different uses.  \code{assign.par}
#' parses them from the vector based on a known structure and assigns them into
#' \code{ddfobj}.  \code{getpar} extracts the requested types to be extracted
#' from \code{ddfobj}.
#'
#' @aliases getpar
#' @param ddfobj distance sampling object (see \code{\link{create.ddfobj}})
#' @param fitting character string which is either "all","key","adjust" which
#'   determines which parameters are retrieved.
#' @note Internal functions not intended to be called by user.
#' @author Jeff Laake, David L Miller
#' @seealso assign.par
#' @keywords utility
getpar <- function(ddfobj,fitting="all"){

  # this needs to be generalized to non K+A setting
  if(fitting=="adjust"){
    ddfobj$pars$scale <- rep(NA,length(ddfobj$pars$scale))
    ddfobj$pars$shape <- rep(NA,length(ddfobj$pars$shape))
  }else if(fitting=="key"){
    ddfobj$pars$adjustment <- rep(NA,length(ddfobj$pars$adjustment))
  }

  return(as.vector(unlist(ddfobj$pars)))

}
