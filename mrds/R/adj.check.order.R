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
#' @return
#' 
#' Logical value, true if there have been errors, false otherwise.
#' @author David Miller
#' @seealso \code{\link{adjfct.cos}}, \code{\link{adjfct.poly}},
#'   \code{\link{adjfct.herm}}, \code{\link{detfct}}, \code{\link{mcds}},
#'   \code{\link{cds}}
#' @references S.T.Buckland, D.R.Anderson, K.P. Burnham, J.L. Laake. 1993.
#'   Robust Models. In: Distance Sampling, eds. S.T.Buckland, D.R.Anderson,
#'   K.P. Burnham, J.L. Laake. Chapman & Hall.
#' @keywords methods
adj.check.order <-
		function(adj.series,adj.order)
{
# dlm 30-May-2006 Added any() so we can handle multiple adjustment orders.
	
# Nothing has gone wrong yet!
	err <- FALSE
	
	if(adj.series == "poly"){
# If polynomial, check even in a very crude way
		
		if(any(as.integer(adj.order/2) != (adj.order/2))){
			errors("Odd polynomial adjustment terms selected")
			err <- TRUE
		}
		
	}else if(adj.series == "herm"){
# If hermite, check even and greater than (or equal to) order 4
		
		if(any(adj.order < 4)){
			errors("Hermite polynomial adjustment terms of order < 4 selected")
			err <- TRUE
		}
		
		if(any(as.integer(adj.order/2) != (adj.order/2))){
			errors("Odd Hermite polynomial adjustment terms selected")
			err <- TRUE
		}
		
	}else{
		
# dlm 17-July-2007  Not sure that this is really necessary...
# Check Cosine is at least 2(this is safe to assume since we already
# checked that we had cos/herm/poly/null
#    if(any(adj.order < 2)){
#      errors("Cosine adjustment terms of order < 2 selected")
#      err <- TRUE
#    }
	}
	
	return(err)
	
}
