#' Cosine adjustment term, not the series.
#' 
#' For internal use only -- not to be called by 'mrds' or 'Distance' users 
#' directly.
#'
#' @param distance perpendicular distance vector/scalar
#' @param scaling scale parameter
#' @param adj.order the adjustment order
#' 
#' @returns scalar or vector containing the cosine adjustment term for every
#' value in \code{distance} argument
#' 
#' @author Felix Petersma
adj.cos <- function(distance, scaling, adj.order) {
  return(cos(adj.order * pi * distance / scaling))
}