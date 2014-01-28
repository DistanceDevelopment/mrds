#' Check order of adjustment terms
#'
#' 'adj.check.order' checks that the Cosine, Hermite or simple polynomials are
#' of the correct order.
#'
#'
#' Only even functions are allowed as adjustment terms. Also Hermite
#' polynomials must be of degree at least 4 and Cosine of order at least 3. If
#' incorrect terms are supplied then \code{stop} is called.
#'
#' @param adj.series Adjustment series used
#'   ('\code{cos}','\code{herm}','\code{poly}')
#' @param adj.order Integer to check
#' @return Logical value, true if there have been errors, false otherwise.
#'  Usually doesn't matter as we already called \code{stop()} if something went
#'  wrong.
#'
#' @author David L Miller
#' @seealso \code{\link{adjfct.cos}}, \code{\link{adjfct.poly}},
#'   \code{\link{adjfct.herm}}, \code{\link{detfct}}, \code{\link{mcds}},
#'   \code{\link{cds}}
#' @references S.T.Buckland, D.R.Anderson, K.P. Burnham, J.L. Laake. 1993.
#'   Robust Models. In: Distance Sampling, eds. S.T.Buckland, D.R.Anderson,
#'   K.P. Burnham, J.L. Laake. Chapman & Hall.
#' @keywords methods
adj.check.order <- function(adj.series,adj.order){

  # Nothing has gone wrong yet!
  err <- FALSE

  if(adj.series == "poly"){
    # check even polynomial
    if(any(as.integer(adj.order/2) != (adj.order/2))){
      stop("Only even polynomial adjustment terms may be used")
      err <- TRUE
    }

  }else if(adj.series == "herm"){
    # Hermite, check even and greater than (or equal to) order 4
    if(any(adj.order < 4)){
      stop("Hermite polynomial adjustment terms must be of order >= 4")
      err <- TRUE
    }

    if(any(as.integer(adj.order/2) != (adj.order/2))){
      stop("Hermite polynomial adjustment terms must be even")
      err <- TRUE
    }
  }

  invisible(err)
}
