#' Distribution of probabilities of detection
#'
#' Generate a table of frequencies of probability of detection from a detection
#' function model. This is particularly useful when employing covariates, as it
#' can indicate if there are detections with very small detection probabilities
#' that can be unduly influential when calculating abundance estimates.
#'
#' Because \code{\link{dht}} uses a Horvitz-Thompson-like estimator, abundance
#' estimates can be sensitive to errors in the estimated probabilities. The
#' estimator is based on \eqn{\sum 1/ \hat{P}_a(z_i)}, which means that the
#' sensitivity is greater for smaller detection probabilities. As a rough
#' guide, we recommend that the method be not used if more than say 5\% of the
#' \eqn{\hat{P}_a(z_i)} are less than 0.2, or if any are less than 0.1. If
#' these conditions are violated, the truncation distance w can be reduced.
#' This causes some loss of precision relative to standard distance sampling
#' without covariates.
#'
#' @param object fitted detection function
#' @param bins how the results should be binned
#' @param proportion should proportions be returned as well as counts?
#' @return a \code{data.frame} with probability bins, counts and (optionally)
#' proportions. The object has an attribute \code{p_range} which contains the
#' range of estimated detection probabilities
#' @references Marques, F.F.C. and S.T. Buckland. 2004. Covariate models for
#' the detection function.
#'   In: Advanced Distance Sampling, eds. S.T. Buckland, D.R. Anderson, K.P.
#'   Burnham, J.L. Laake, D.L. Borchers, and L. Thomas. Oxford University
#'   Press.
#'
#' @export
#' @author David L Miller
#' @examples
#' \dontrun{
#' # try out the tee data
#' data(book.tee.data)
#' egdata <- book.tee.data$book.tee.dataframe
#' # fit model with covariates
#' result <- ddf(dsmodel = ~mcds(key = "hn", formula = ~sex+size),
#'               data = egdata[egdata$observer==1, ], method = "ds",
#'               meta.data = list(width = 4))
#' # print table
#' p.dist.table(result)
#' # with proportions
#' p.dist.table(result, proportion=TRUE)
#' }
# Update Distance::p_dist_table when changing pars here!!
p.dist.table <- function(object, bins=seq(0, 1, by=0.1), proportion=FALSE){

  # if we have a ds object from Distance
  if("dsmodel" %in%class(object)){
    object <- object$ddf
  }

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

#' @rdname p.dist.table
#' @export
p_dist_table <- p.dist.table
