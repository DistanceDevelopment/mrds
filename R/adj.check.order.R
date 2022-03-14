#' Check order of adjustment terms
#'
#' 'adj.check.order' checks that the Cosine, Hermite or simple polynomials are
#' of the correct order.
#'
#' Only even functions are allowed as adjustment terms, per p.47 of Buckland et
#' al (2001). If incorrect terms are supplied then an error is throw via
#' \code{stop}.
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
adj.check.order <- function(adj.series, adj.order, key){

  if(adj.series %in% "poly"){
    # If polynomial, check even in a very crude way
    if(any(as.integer(adj.order/2) != (adj.order/2))){
      stop("Odd polynomial adjustment terms selected")
    }
    if(key != "unif" & any(adj.order < 4)){
      stop("Polynomial adjustment terms of order < 4 selected")
    }

  }else if(adj.series == "herm"){
    # If hermite, check even and greater than (or equal to) order 4
    if(any(as.integer(adj.order/2) != (adj.order/2))){
      stop("Odd Hermite polynomial adjustment terms selected")
    }

    if(key != "unif" & any(adj.order < 4)){
      stop("Hermite polynomial adjustment terms of order < 4 selected")
    }

  }

  adj.name <- switch(adj.series,
                     cos = "Cosine",
                     herm = "Hermite",
                     poly = "Simple polynomial",
                     NULL)

  key.name <- switch(key,
                     hn = "half-normal",
                     hr = "hazard-rate",
                     NULL)

  if((key %in% c("hn","hr")) & any(adj.order<2)){
    stop(paste0(adj.name," adjustments must be of order >2 for ",
                key.name," key functions"))
  }

  invisible()
}
