# Cosine adjustment term, not the series
#
# distance  - perpendicular distance vector
# scaling     - scale parameter
# adj.order  - scalar of adjustment order
adj.cos <- function(distance, scaling, adj.order) {
  return(cos(adj.order * pi * distance) / scaling)
}