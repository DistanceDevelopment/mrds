#' Create bins from a set of binned distances and a set of cutpoints.
#'
#' This is an internal routine and shouldn't be necessary in normal analyses.
#'
#' @param data `data.frame` with at least the column `distance`.
#' @param cutpoints vector of cutpoints for the bins
#'
#' @return argument `data` with two extra columns `distbegin` and
#'        `distend`.
#'
#' @author David L. Miller
#' @export
create.bins <- function(data, cutpoints){

  # remove distances outside bins
  in.cp.ind <- data$distance >= cutpoints[1] & data$distance<= cutpoints[length(cutpoints)]
  if(!all(in.cp.ind)){
    warning("Some distances were outside bins and have been removed.")
  }
  data <- data[in.cp.ind, , drop=FALSE]

  # use cut() to create bins
  chopped <- cut(data$distance, 
                 breaks=cutpoints, 
                 include.lowest=TRUE, 
                 labels = FALSE)
  
  distbegin <- cutpoints[1:(length(cutpoints)-1)]
  distend <- cutpoints[2:length(cutpoints)]
  
  # put all that together and make a data.frame
  data <- cbind(data,
                # process to get bin beginnings/endings
                distbegin = distbegin[chopped],
                distend = distend[chopped])
  data <- data.frame(data)

  return(data)
}
