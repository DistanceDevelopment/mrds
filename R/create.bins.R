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

  # lazy typist
  cp <- cutpoints

  # remove distances outside bins
  in.cp.ind <- data$distance>=cp[1] & data$distance<=cp[length(cp)]
  if(!all(in.cp.ind)){
    warning("Some distances were outside bins and have been removed.")
  }
  data <- data[in.cp.ind, , drop=FALSE]

  # use cut() to create bins
  chopped <- as.character(cut(data$distance, breaks=cp, include.lowest=TRUE))

  # use cut() to create bins
  chopped <- cut(data$distance, breaks=cp, include.lowest=TRUE, labels = FALSE)
  
  distbegin <- cp[1:(length(cp)-1)]
  distend <- cp[2:length(cp)]
  
  # put all that together and make a data.frame
  data <- cbind(data,
                # process to get bin beginnings/endings
                distbegin = distbegin[chopped],
                distend = distend[chopped])
  data <- data.frame(data)

  return(data)
}
