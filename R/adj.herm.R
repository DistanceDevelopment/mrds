# Herminte polynomial adjustment term, not the series
#
# distance  - perpendicular distance vector
# scaling     - scale parameter
# adj.order  - scalar of adjustment order
#
#' @aliases  hermite.poly
adj.herm <- function(distance, scaling, adj.order) {
  return(hermite.poly((distance / scaling), adj.order))
}