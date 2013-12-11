#' Plot fit of detection functions and histograms of data from distance
#' sampling independent observer model with full independence (\code{io.fi})
#'
#' Plots the fitted detection functions for a distance sampling model and
#' histograms of the distances (for unconditional detection functions) or
#' proportion of observations detected within distance intervals (for
#' conditional detection functions) to compare visually the fitted model and
#' data.
#'
#' The structure of the histogram can be controlled by the user-defined
#' arguments \code{nc} or \code{breaks}.  The observation specific detection
#' probabilities along with the line representing the fitted average detection
#' probability.
#'
#' It is not intended for the user to call any of \code{plot.io.fi} but the
#' arguments are documented here. The generic \code{plot} command should be used#' and will call the appropriate function based on the \code{ddf} object.
#'
#' @aliases plot.io.fi
#' @S3method plot io.fi
#' @method plot io.fi
#' @export
#' @param x fitted model from \code{ddf}
#' @param which index to specify which plots should be produced.
#'  \tabular{ll}{1 \tab Plot histogram of observed distances for observer 1\cr
#'               2 \tab Plot histogram of observed distances for observer 2\cr
#'               3 \tab Plot primary unconditional detection function \cr
#'               4 \tab Plot secondary unconditional detection function \cr
#'               5 \tab Plot pooled unconditional detection function \cr
#'               6 \tab Plot duplicate unconditional detection function \cr
#'               7 \tab Plot primary conditional detection function \cr
#'               8 \tab Plot secondary conditional detection function \cr}
#' @param breaks user define breakpoints
#' @param nc number of equal-width bins for histogram
#' @param maintitle main title line for each plot
#' @param showpoints logical variable; if TRUE plots predicted value for each
#'   observation
#' @param showlines logical variable; if TRUE a line representing the average
#'   detection probability is plotted
#' @param ylim range of y axis; defaults to (0,1)
#' @param angle shading angle for hatching
#' @param density shading density for hatching
#' @param col plotting colour
#' @param jitter scaling option for plotting points.  Jitter is applied to
#'   points by multiplying the fitted value by a random draw from a normal
#'   distribution with mean 1 and sd jitter.
#' @param divisions number of divisions for averaging line values; default = 25
#' @param pages the number of pages over which to spread the plots. For
#'  example, if \code{pages=1} then all plots will be displayed on one page.
#'  Default is 0, which prompts the user for the next plot to be displayed.
#' @param xlab label for x-axis
#' @param ylab label for y-axis
#' @param subtitle if TRUE, shows plot type as sub-title
#' @param \dots other graphical parameters, passed to the plotting functions
#'   (\code{plot}, \code{hist}, \code{lines}, \code{points}, etc)
#' @return Just plots.
#' @author Jeff Laake, Jon Bishop, David Borchers, David L Miller
#' @keywords plot
#' @examples
#' \donttest{
#' library(mrds)
#' data(book.tee.data)
#' egdata <- book.tee.data$book.tee.dataframe
#' result.io.fi <- ddf(mrmodel=~glm(~distance), data = egdata, method = "io.fi",
#'               meta.data = list(width = 4))
#'
#' # just plot everything
#' plot(result.io.fi)
#'
#' # Plot primary and secondary unconditional detection functions on one page
#' # and  primary and secondary conditional detection functions on another
#' plot(result.io.fi,which=c(1,2,5,6),pages=2)
#' }
plot.io.fi <- function(x, which=1:8, breaks=NULL, nc=NULL, maintitle="",
                       showlines=TRUE, showpoints=TRUE,ylim=c(0,1),angle=-45,
                       density=20,col="black",jitter=NULL,divisions=25,pages=0,
                       xlab="Distance",ylab="Detection probability",
                       subtitle=TRUE,...){
  # Functions used: process.data, predict.io.fi, plot_uncond,plot_cond


  # Retrieve values from model object
  model <- x
  xmat <- model$mr$data
  xmat$offsetvalue <- 0
  cond.det <- predict(model,newdata=xmat,integrate=FALSE)
  fitted <-cond.det$fitted
  p1 <- cond.det$p1
  p2 <- cond.det$p2
  width <- model$meta.data$width
  left <- model$meta.data$left

  # If number of classes for histogram intervals was not set
  # compute a reasonable default
  if(is.null(nc)){
    nc<-round(sqrt(min(length(xmat$distance[xmat$observer==1&xmat$detected==1]),
                       length(xmat$distance[xmat$observer==2&xmat$detected==1]),
                       length(xmat$distance[xmat$observer==1&
                                            xmat$timesdetected==2]))),0)
  }

  # Set up default break points unless specified
  if(model$meta.data$binned){
    breaks <- model$meta.data$breaks
    nc <- length(breaks)-1
  }else{
    if(is.null(breaks)){
      breaks <- left + ((width-left)/nc)*(0:nc)
    }else{
      nc <- length(breaks)-1
    }
  }

  # make a list of the gxvalues we would want
  gxlist <- list(p1,
                 p2,
                 p1+p2-p1*p2,
                 p1*p2,
                 p1[xmat$detected[xmat$observer==2]==1],
                 p2[xmat$detected[xmat$observer==1]==1])

  # do the plotting layout
  oask <- plot.layout(which,pages)
  on.exit(devAskNewPage(oask))

  # plot histograms
  for(wh in which[which < 3]){
    if(maintitle!="") maintitle <- paste(maintitle,"\n",sep="")
    mt <- paste(maintitle, "Observer = ",wh, " detections")

    hist(xmat$distance[xmat$observer==wh & xmat$detected==1],breaks=breaks,
         main=mt,angle=angle,density=density,col=col,xlab=xlab,ylab=ylab)
  }


  # 3 - Plot primary unconditional detection function
  # 4 - Plot secondary unconditional detection function
  # 5 - Plot pooled unconditional detection function
  # 6 - Plot duplicate unconditional detection function
  for(wh in which[which > 2 & which < 7]){
    observer <- wh-2
    plot_uncond(model,observer,xmat,gxvalues=gxlist[[observer]],nc,
                finebr=(width/divisions)*(0:divisions),
                breaks,showpoints,showlines,maintitle,ylim,
                angle=angle,density=density,
                col=col,jitter=jitter,xlab=xlab,ylab=ylab,...)
  }

  # 7 - Plot conditional detection functions
  data <- model$mr$data
  data$offsetvalue <- 0
  if(is.element(7,which)){
    gxvalues <- p1[xmat$detected[xmat$observer==2]==1]
    plot_cond(1,data,gxvalues,model,nc,breaks,
              finebr=(width/divisions)*(0:divisions),showpoints,showlines,
              maintitle,ylim,angle=angle,density=density,col=col,jitter=jitter,
              xlab=xlab,ylab=ylab,subtitle=subtitle,...)
  }

  # 8 - Plot secondary conditional detection function
  if(is.element(8,which)){
    gxvalues <-p2[xmat$detected[xmat$observer==1]==1]
    plot_cond(2,data,gxvalues,model,nc,breaks,
              finebr=(width/divisions)*(0:divisions),showpoints,showlines,
              maintitle,ylim,angle=angle,density=density,col=col,jitter=jitter,
              xlab=xlab,ylab=ylab,subtitle=subtitle,...)
  }
  invisible(NULL)
}
