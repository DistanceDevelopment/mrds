#' Distribution of probabilities of detection
#'
#' Generate a table of frequencies of probability of detection from a detection function model. This is particularly useful when employing covariates, as it can indicate if there are detections with very small \hat{P}s that can be unduly influential when calculating abundance estimates.
#'
#' @param object fitted detection function
#' @param bins how the results should be binned
#' @param proportion should proportions be returned as well as counts?
#' @return a \code{data.frame} with probability bins, counts and (optionally) proportions. The object has an attribute \code{p_range} which contains the range of estimated detection probabilities
#' @export
#' @author David L Miller
p_dist_table <- function(object, bins=seq(0, 1, by=0.1), proportion=FALSE){

  # get the probabilities from the fitted object
  ps <- fitted(object)

  # do the binning, make a data.frame
  tab <- as.data.frame(table(cut(ps, bins, include.lowest=TRUE)))

  # format column names
  names(tab) <- c("p", "count")

  # format the labels
  tab$p <- sub("\\(", "", tab$p)
  tab$p <- sub("]", "", tab$p)
  tab$p <- sub("\\[", "", tab$p)
  tab$p <- sub(",", " - ", tab$p)

  if(proportion){
    tab$proportion <- tab$count/sum(tab$count)
  }

  attr(tab, "p_range") <- range(ps)

  class(tab) <- c("p_dist_table", "data.frame")
  return(tab)
}


#' Print distribution of probabilities of detection
#'
#' Just a pretty printer for the table of probabilities of detection.

#' @param x output from \code{\link{p_dist_table}}
#' @param \dots other arguments to be passed to \code{\link{print.data.frame}}
#' @param digits number of significant digits to print
#' @return just prints the table and the range of ps
#' @author David L Miller
#' @export
print.p_dist_table <- function(x, digits=2, ...){
  class(x) <- c("data.frame")
  print(x, row.names=FALSE, digits=digits, ...)

  cat("Range of probabilities: ",
      paste0(signif(attr(x, "p_range"), digits), collapse=" - "), "\n")
}
