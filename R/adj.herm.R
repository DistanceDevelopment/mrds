#' Hermite polynomial adjustment term, not the series.
#' 
#' For internal use only -- not to be called by 'mrds' or 'Distance' users 
#' directly.
#'
#' @param distance perpendicular distance vector/scalar
#' @param scaling scale parameter
#' @param adj.order the adjustment order
#' 
#' value in \param{distance}
#' @returns scalar or vector containing the Hermite adjustment term for every
#'
#' @author Felix Petersma
adj.herm <- function(distance, scaling, adj.order) {
  return(hermite.poly((distance / scaling), adj.order))
}
