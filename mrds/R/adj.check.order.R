#' Check order of adjustment terms
#'
#' 'adj.check.order' checks that the Cosine, Hermite or simple polynomials are
#' of the correct order.
#'
#'
#' Only even functions are allowed as adjustment terms. Also Hermite
#' polynomials must be of degree at least 4 and Cosine of order at least 3. If
#' incorrect terms are supplied then \code{errors} is called.
#'
#' @param adj.series Adjustment series used
#'   ('\code{cos}','\code{herm}','\code{poly}')
#' @param adj.order Integer to check
#' @param key key function to be used with this adjustment series
#' @return Nothing! Just calls \code{stop} if something goes wrong.
#'
#' @author David Miller
#' @seealso \code{\link{adjfct.cos}}, \code{\link{adjfct.poly}},
#'   \code{\link{adjfct.herm}}, \code{\link{detfct}}, \code{\link{mcds}},
#'   \code{\link{cds}}
#' @references S.T.Buckland, D.R.Anderson, K.P. Burnham, J.L. Laake. 1993.
#'   Robust Models. In: Distance Sampling, eds. S.T.Buckland, D.R.Anderson,
#'   K.P. Burnham, J.L. Laake. Chapman & Hall.
#' @keywords methods
adj.check.order <- function(adj.series,adj.order,key){

  if(adj.series == "poly"){
    # If polynomial, check even in a very crude way
    if(any(as.integer(adj.order/2) != (adj.order/2))){
      stop("Odd polynomial adjustment terms selected")
    }

  }else if(adj.series == "herm"){
    # If hermite, check even and greater than (or equal to) order 4
    if(any(adj.order < 4)){
      stop("Hermite polynomial adjustment terms of order < 4 selected")
    }

    if(any(as.integer(adj.order/2) != (adj.order/2))){
      errors("Odd Hermite polynomial adjustment terms selected")
    }

  }else if(adj.series=="cos"){
    if((key %in% c("hn","hr")) & any(adj.order<2)){
      stop("Cosine adjustments must be of order >2 for half-normal and hazard-rate key functions")
    }
  }

  invisible()
}
