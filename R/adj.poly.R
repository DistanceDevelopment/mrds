# Simple polynomial adjustment term, not the series
#
# distance  - perpendicular distance vector
# scaling     - scale parameter
# adj.order  - scalar of adjustment order
adj.poly <- function(distance, scaling, adj.order) {
  return((distance / scaling) ^ adj.order)
}