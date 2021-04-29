#' Observation detection tables
#'
#' Plot the tables created by \code{\link{det.tables}}. Produces a series of
#' tables for dual observer data that shows the number missed and detected for
#' each observer within defined distance classes.
#'
#' Plots that are produced are as follows (controlled by the \code{which}
#' argument):
#' \describe{
#'   \item{1}{Detected by either observer/Detected by observer 1}
#'   \item{2}{Detected by either observer/Detected by observer 2}
#'   \item{3}{Seen by both observers}
#'   \item{4}{Seen by either observer}
#'   \item{5}{Detected by observer 2/Detected by observer 1 | 2}
#'   \item{6}{Detected by observer 1/Detected by observer 2 | 1}
#' }
#'
#'
#' @aliases plot.det.tables
#' @export
#' @param x object returned by \code{\link{det.tables}}
#' @param which items in x to plot (vector with values in 1:6)
#' @param angle shading angle for hatching
#' @param density shading density for hatching
#' @param col1 plotting colour for total histogram bars.
#' @param col2 plotting colour for subset histogram bars.
#' @param new if \code{TRUE} new plotting window for each plot
#' @param \dots other graphical parameters, passed to plotting functions
#' @return Just plots.
#' @author Jeff Laake, David L Miller
#' @importFrom grDevices dev.new
#' @importFrom graphics legend
#' @examples
#' \donttest{
#' data(book.tee.data)
#' region <- book.tee.data$book.tee.region
#' egdata <- book.tee.data$book.tee.dataframe
#' samples <- book.tee.data$book.tee.samples
#' obs <- book.tee.data$book.tee.obs
#' xx <- ddf(mrmodel=~glm(formula=~distance*observer),
#'           dsmodel = ~mcds(key = "hn", formula = ~sex),
#'           data = egdata, method = "io", meta.data = list(width = 4))
#' tabs <- det.tables(xx,breaks=c(0,.5,1,2,3,4))
#' par(mfrow=c(2,3))
#' plot(tabs,which=1:6,new=FALSE)
#' }
plot.det.tables <- function(x, which=1:6, angle=NULL, density=NULL,
                            col1="white", col2="lightgrey", new=TRUE,...){

  # plotting function that actually does the work
  plot_seen <- function(x, col1, col2, leg.title, ...){
    if(new & .Platform$GUI=="Rgui") dev.new()
    missed <- x[, "Missed"]
    detected <- x[, "Detected"]
    ymax <- max(missed+detected)

    # plot "background" total
    histline(missed+detected, breaks=breaks, lineonly=TRUE, ylim=c(0, ymax),
             xlab="Distance", ylab="Frequency", angle=angle,
             density=density, col=col1, add=FALSE, ...)

    # plot "foreground" detected
    histline(detected, breaks, lineonly=TRUE, col=col2, add=TRUE,
             density=density, det.plot=TRUE, ...)

    legend("topright", legend=leg.title, fill=c(col1, col2))
  }

  breaks <- x$breaks

  if(is.element(1, which) & !is.null(x$Observer1)){
    plot_seen(x$Observer1, col1, col2, c("Detected by either observer",
                                         "Detected by observer 1"), ...)
  }

  if(is.element(2, which) & !is.null(x$Observer2)){
    plot_seen(x$Observer2, col1, col2, c("Detected by either observer",
                                         "Detected by observer 2"), ...)
  }

  if(is.element(3, which) & !is.null(x$Duplicates)){
    ymax <- max(x$Duplicates)
    histline(x$Duplicates, breaks=breaks, lineonly=FALSE, xlab="Distance",
             ylab="Frequency", angle=angle, density=density, col=col1,
             ylim=c(0, ymax), ...)
    legend("topright", legend=c("Seen by both observers"), fill=c(col1))
  }

  if(is.element(4, which) & !is.null(x$Pooled)){
    ymax <- max(x$Pooled)
    histline(x$Pooled, breaks=breaks, lineonly=FALSE, xlab="Distance",
             ylab="Frequency", angle=angle, density=density, col=col1,
             ylim=c(0, ymax), ...)
    legend("topright", legend=c("Seen by either observer"), fill=c(col1))
  }

  if(is.element(5, which) & !is.null(x$Obs1_2)){
    plot_seen(x$Obs1_2, col1, col2, c("Detected by observer 2",
                                      "Detected by observer 1 | 2"), ...)
  }

  if(is.element(6, which) & !is.null(x$Obs2_1)){
    plot_seen(x$Obs2_1, col1, col2, c("Detected by observer 1",
                                      "Detected by observer 2 | 1"), ...)
  }

  invisible()
}
