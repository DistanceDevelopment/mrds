#' Prints density and abundance estimates
#'
#' Outputs summary statistics, abundance and density by region (if any) and
#' optionally a correlation matrix if more than one region.
#'
#' @method print dht
#' @param x dht object that results from call to dht for a specific ddf object
#' @param cor if TRUE outputs correlation matrix of estimates
#' @param bysample if TRUE, prints results for each sample
#' @param vcmatrices if TRUE, prints variance-covariance matrices
#' @param clusters if TRUE, prints the abundance and density of clusters, as
#' well as individuals
#' @param cluster.summary if TRUE prints the expected and observed group sizes
#' @param \dots unspecified and unused arguments for S3 consistency
#' @export
#' @return None
#' @author Jeff Laake
#' @seealso \code{\link{dht}}
#' @keywords utility
print.dht <- function(x, cor=FALSE, bysample=FALSE, vcmatrices=FALSE,
                      clusters=FALSE, cluster.summary=FALSE, ...){

  print.tables <- function(x, cor, bysample, vcmatrices){
    if("N" %in% names(x)){
      cat("\nAbundance:\n")
      # Change Label to Region for output purposes only (won't break dependent packages or code!)
      xN <- x$N
      index <- which(names(x$N) == "Label")
      if(length(index) > 0){
        names(xN)[index] <- "Region"
      }
      print(xN)
    }
    cat("\nDensity:\n")
    # Change Label to Region for output purposes only (won't break dependent packages or code!)
    xD <- x$D
    index <- which(names(x$D) == "Label")
    if(length(index) > 0){
      names(xD)[index] <- "Region"
    }
    print(xD)
    if(vcmatrices){
      cat("\nAbundance variance-covariance matrix component from estimating detection function\n" )
      print(x$vc$detection$variance)
      cat("\n")
      cat("Abundance variance-covariance matrix component from sample selection")
      cat("\n")
      print(x$vc$er)
    }
    if(cor){
      cat("\nCorrelation matrix:\n")
      print(x$cormat)
    }
    if(bysample){
      cat("\nEstimates by sample:\n")
      print(x$bysample)
    }
  }

  # general information
  cat("Abundance and density estimates from distance sampling\n")
  cat("Variance       :", paste0(attr(x, "ER_var")[1], ","),
      ifelse(attr(x, "ER_var")[3], "binomial",
             ifelse(attr(x, "ER_var")[2], "N/L", "n/L")), "\n")

  # now print
  if(is.null(x$clusters)){
    # summary statistics
    cat("\nSummary statistics\n\n")
    print(x$individuals$summary)
    print.tables(x$individuals, cor, bysample, vcmatrices)
    cat("\nVariance contributions\n")
    print(attr(x, "variance_components")$individuals)
  }else{
    # summary statistics
    cat("\nSummary statistics\n\n")
    print(x$clusters$summary)

    if(clusters){
      cat("\nSummary for clusters\n")
      print.tables(x$clusters, cor, bysample, vcmatrices)
      cat("\nVariance contributions\n")
      print(attr(x, "variance_components")$clusters)
    }

    cat("\nSummary for individuals\n")
    print.tables(x$individuals, cor, bysample, vcmatrices)

    if(cluster.summary){
      cat("\nCluster size summary\n")
      print(as.data.frame(x$Expected.S))
    }

    cat("\nVariance contributions\n")
    print(attr(x, "variance_components")$individuals)
  }

  invisible()
}
