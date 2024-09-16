#' Simple polynomial adjustment term, not the series.
#' 
#' For internal use only -- not to be called by 'mrds' or 'Distance' users 
#' directly.
#'
#' @param distance perpendicular distance vector/scalar
#' @param scaling scale parameter
#' @param adj.order the adjustment order
#' 
#' value in \param{distance}
#' @returns scalar or vector containing the polynomial adjustment term for every
#'
#' @author Felix Petersma
adj.poly <- function(distance, scaling, adj.order) {
  return((distance / scaling) ^ adj.order)
}