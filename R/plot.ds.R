#' Plot fit of detection functions and histograms of data from distance sampling model
#'
#' Plots the fitted detection function(s) with a histogram of the observed distances to compare visually the fitted model and data.
#'
#' The structure of the histogram can be controlled by the user-defined arguments \code{nc} or \code{breaks}. The observation specific detection probabilities along with the line representing the fitted average detection probability.
#'
#' It is not intended for the user to call \code{plot.ds} but its arguments are documented here. Instead the generic \code{plot} command should be used and it will call the appropriate function based on the class of the \code{ddf} object.
#'
#' @aliases plot.ds
#' @export
#' @param x fitted model from \code{ddf}.
#' @param which index to specify which plots should be produced:
#'  \tabular{ll}{1 \tab histogram of observed distances\cr
#'               2 \tab histogram of observed distanes with fitted line and points (default)\cr}
#' @param breaks user defined breakpoints
#' @param nc number of equal width bins for histogram
#' @param jitter.v scaling option for plotting points.  Jitter is applied to points by multiplying the fitted value by a random draw from a normal distribution with mean 1 and sd \code{jitter.v[j]}.  Where \code{j=1,2} corresponds to observer \code{j} and \code{j=3} corresponds to pooled/duplicate detections.
#' @param showpoints logical variable; if \code{TRUE} plots predicted value for each observation.
#' @param subset subset of data to plot.
#' @param pl.col colours plotting colours for obs 1, obs 2 detections.
#' @param bw.col grayscale plotting colours for obs 1, obs 2 detections.
#' @param black.white logical variable; if \code{TRUE} plots are grayscale.
#' @param pl.den shading density for plots of obs 1, obs 2 detections.
#' @param pl.ang shading angle for plots of obs 1, obs 2 detections.
#' @param main user-specfied plot title.
#' @param ylim user-specified y axis limits.
#' @param pdf plot the histogram of distances with the PDF of the probability of detection overlaid. Ignored (with warning) for line transect models.
#' @param pages the number of pages over which to spread the plots. For example, if \code{pages=1} then all plots will be displayed on one page. Default is 0, which prompts the user for the next plot to be displayed.
#' @param xlab label for the x axis.
#' @param \dots other graphical parameters, passed to the plotting functions (\code{\link{plot}}, \code{\link{hist}}, \code{\link{lines}}, \code{\link{points}}, etc).
#' @return Just plots.
#' @author Jeff Laake, Jon Bishop, David Borchers, David L Miller
#' @keywords plot
#' @importFrom graphics hist par lines points title
#' @importFrom grDevices devAskNewPage grey
#' @importFrom stats rnorm
#' @examples
#' \donttest{
#' # fit a model to the tee data
#' data(book.tee.data)
#' egdata <- book.tee.data$book.tee.dataframe
#' xx <- ddf(dsmodel=~mcds(key="hn", formula=~sex),
#'           data=egdata[egdata$observer==1, ],
#'           method="ds", meta.data=list(width=4))
#'
#' # not showing predicted probabilities
#' plot(xx, breaks=c(0, 0.5, 1, 2, 3, 4), showpoints=FALSE)
#'
#' # two subsets
#' plot(xx, breaks=c(0, 0.5, 1, 2, 3, 4), subset=sex==0)
#' plot(xx, breaks=c(0, 0.5, 1, 2, 3, 4), subset=sex==1)
#'
#' # put both plots on one page
#' plot(xx, breaks=c(0, 0.5, 1, 2, 3, 4), pages=1, which=1:2)
#' }
plot.ds <- function(x, which=2, breaks=NULL, nc=NULL,
                    jitter.v=rep(0,3), showpoints=TRUE, subset=NULL,
                    pl.col='black', bw.col=grey(0), black.white=FALSE,
                    pl.den=rep(20,1), pl.ang=rep(-45,1), main=NULL, pages=0,
                    pdf=FALSE, ylim=NULL, xlab="Distance", ...){

  model<-x
  lower <- 0
  vname <- "distance"
  dat <- model$data

  # ignore pdf=TRUE with line transect data
  if(pdf & !model$meta.data$point){
    warning("Ignoring pdf=TRUE for line transect data")
    pdf <- FALSE
  }

  # decide which plots to show
  show <- rep(FALSE, 2)
  show[which] <- TRUE

  ## Set printing options for plots:
  # By default  pl.col=c('black'),
  #             bw.col=c(grey(0))
  # If greyscale plots are required use the following colours:
  if(black.white){
    byval1 <- bw.col[1]
  }else{
    # If colour plots are required use the following:
    byval1 <- pl.col[1]
  }
  # Density of shaded lines - default is set all to 20
  denval1 <- pl.den[1]
  # Angle of shaded lines - default is set all to -45
  angval1 <- pl.ang[1]

  # code from dpexplot:
  width <- model$meta.data$width
  left <- model$meta.data$left
  ddfobj <- model$ds$aux$ddfobj
  point <- model$ds$aux$point
  if(is.null(model$ds$aux$int.range)){
    int.range <- c(0,width)
  }else{
    int.range <- model$ds$aux$int.range
  }
  if(is.matrix(int.range)){
    max.range <- as.vector(int.range[1,])
    int.range <- int.range[2:dim(int.range)[1],]
    range.varies <- TRUE
  }else{
    max.range <- int.range
    normalize <- FALSE
    range.varies <- FALSE
  }

  if(range.varies & showpoints){
    warning("Point values can be misleading for g(x) when the range varies")
  }

  ## make selection of the data subset to plot
  if(!is.null(substitute(subset))){
    selected <- eval(substitute(subset),ddfobj$xmat)
  }else{
    selected <- rep(TRUE,nrow(ddfobj$xmat))
  }
  # die if there was nothing selected
  if(all(!selected)){
    stop("Specified subset is empty.")
  }

  if(is.matrix(int.range)){
    int.range <- int.range[selected,]
  }

  xmat <- ddfobj$xmat[selected,]
  if(!is.null(ddfobj$scale)){
    z <- ddfobj$scale$dm[selected,,drop=FALSE]
  }else{
    z <- matrix(1,nrow=1,ncol=1)
  }

  if(length(model$fitted)==1){
    pdot <- rep(model$fitted, sum(as.numeric(selected)))
  }else{
    pdot <- model$fitted[selected]
    Nhat <- sum(1/pdot)
  }

  zdim <- dim(z)[2]
  n <- length(xmat$distance)

  if(!is.null(breaks)){
    nc <- length(breaks)-1
  }

  if(is.null(nc)){
    nc <- round(sqrt(n), 0)
  }

  # Set logical hascov=TRUE when detection function has
  #  covariates other than distance and observer
  hascov <- FALSE
  if(!ddfobj$intercept.only){
    hascov <- TRUE
  }

  # Compute a grid for distance (xgrid), and covariates zgrid for
  # plotting of detection functions.
  if(!hascov){
    xgrid <- seq(0, width, length.out=101)
    zgrid <- matrix(rep(z[1,],length(xgrid)), byrow=TRUE, ncol=sum(zdim))
  }

  # create intervals of distance (breaks) for the chosen number of classes (nc).
  if(is.null(breaks)){
    if(is.null(model$meta.data$binned)){
      binned <- FALSE
    }else{
      binned <- model$meta.data$binned
    }
    if(binned){
      breaks <- model$ds$aux$breaks
      nc <- length(breaks)-1
    }else{
      breaks <- c(max(0, (max.range[1])),
                  max.range[1]+((max.range[2]-max.range[1])/nc)*(1:nc))
      if(breaks[1]>left){
        breaks <- c(left,breaks)
        nc <- nc+1
      }
    }
  }

  # test breaks for validity and reset as needed
  breaks <- test.breaks(breaks,model$meta.data$left,width)
  nc <- length(breaks)-1
  lower <- min(breaks)
  upper <- max(breaks)
  dat <- dat[selected,]
  keep <- dat[ ,vname]>=lower & dat[ ,vname]<=upper

  # get the histogram object
  hist.obj <- hist(dat[ ,vname][keep], breaks=breaks, plot=FALSE)
  # what's the top of the largest bar?
  ymax <- max(hist.obj$counts)

  # Rescaling for the histogram
  if(normalize & !point){
    bindata <- function(x, r, breaks){
      return(hist(r[r>=x[1] & r<=x[2]], breaks=breaks, plot=FALSE)$counts)
    }
    sumit<-function(x,n,wt){
      return(sum(x/(wt*n)))
    }
    expected.counts <- apply(int.range, 1, bindata,
                             r=(0:1000)*width/1001, breaks=breaks)
    expected.counts <- apply(expected.counts, 1, sumit, n=1001, wt=pdot)
  }else{
    if(!point){
      expected.counts <- (breaks[2:(nc+1)]-breaks[1:nc])*(Nhat/breaks[nc+1])
    }else{
      if(!pdf){
        expected.counts <- -apply(matrix(c(breaks[2:(nc+1)]^2, breaks[1:nc]^2),
                                         ncol=2, nrow=nc),
                                  1, diff)*(Nhat/breaks[nc+1]^2)
      }else{
        expected.counts <- sum(hist.obj$counts)
      }
    }
  }

  # rescale the histogram object by the expected counts
  # but only if we don't have point/pdf plots
  if(!(pdf & point)){
    hist.obj$density <- hist.obj$counts/expected.counts
    hist.obj$density[expected.counts==0] <- 0
  }
  hist.obj$equidist <- FALSE

  ### Actual plotting starts here

  # do the paging, using devAskNewPage() means we can just call plots and
  # R will make the prompt for us
  if(pages!=1 & sum(show)!=1){
    oask <- devAskNewPage(TRUE)
    on.exit(devAskNewPage(oask))
  }else if(sum(show)!=1){
    opar <- par(mfrow=c(1, sum(show)))
    on.exit(par(opar))
  }

  # Data summary plot
  if(show[1]){
    if(is.null(ylim)) ylim <- c(0, ymax)
    histline(hist.obj$counts, breaks=breaks, lineonly=FALSE, ylim=ylim,
             xlab=xlab, ylab="Frequency", angle=angval1,
             density=denval1, col=byval1, ...)
  }

  ## Detection function plot overlaid on histogram of observed distances
  if(show[2]){

    # area under the histogram
    hist_area <- sum(hist.obj$density*diff(breaks))
    # Detection function/pdf values for points to be plotted
    if(point & pdf){
      point_vals <- distpdf(xmat$distance, ddfobj, width=width, point=TRUE,
                                standardize=TRUE)/
                    integratepdf(ddfobj, select=selected, width=width,
                                 int.range=int.range, standardize=TRUE,
                                 point=TRUE)
    }else{
      point_vals <- detfct(xmat$distance, ddfobj, select=selected, width=width)
    }

    # set y labels, limits and tick marks (det.plot) depending on if we
    # are plotting PDF or df
    if(is.null(ylim)) ylim<-c(0, max(hist.obj$density, max(point_vals)))
    if(pdf){
      ylab <- "Probability density"
      det.plot <- FALSE
    }else{
      ylab <- "Detection probability"
      det.plot <- TRUE
    }

    # plot the histogram
    histline(hist.obj$density, breaks=breaks, lineonly=FALSE,
             xlab=xlab, ylab=ylab, ylim=ylim,
             angle=angval1, density=denval1, col=byval1,
             det.plot=det.plot, ...)

    # when we have covariates
    if(hascov){
      finebr <- seq(0, width, length.out=101)
      xgrid <- NULL
      linevalues <- NULL
      newdat <- xmat
      for(i in 1:(length(finebr)-1)){
        x <- (finebr[i]+finebr[i+1])/2
        xgrid <- c(xgrid, x)
        newdat$distance <- rep(x, nrow(newdat))

        detfct.values <- detfct(newdat$distance, ddfobj, select=selected,
                                width=width)

        if(!normalize&range.varies){
          detfct.values[x<int.range[, 1] | x>int.range[, 2]] <- 0
        }

        if(point & pdf){
          ## calculate the pdf of distances
          # want r g(r) / int r g(r) dr

          # this is 2 r g(r)/w^2
          r_gr <- distpdf(newdat$distance, ddfobj, width=width, point=TRUE,
                          standardize=TRUE)
          # this is the value of int [2 r g(r) /w^2] dr
          int_r_gr <- integratepdf(ddfobj, select=selected, width=width,
                                   int.range=int.range, standardize=TRUE,
                                   point=TRUE)
          # so the pdf values are:
          pdf_vals <- r_gr/int_r_gr

          # now rescale such that area under pdf == area under histogram
          vals <- pdf_vals * hist_area
          linevalues <- c(linevalues, sum(vals/pdot)/sum(1/pdot))
        }else{
          linevalues <- c(linevalues, sum(detfct.values/pdot)/sum(1/pdot))
        }
      }
    ## without covariates
    }else{
      if(!is.null(ddfobj$scale)){
        ddfobj$scale$dm <- ddfobj$scale$dm[rep(1, length(xgrid)), ,drop=FALSE]
      }
      if(!is.null(ddfobj$shape)){
        ddfobj$shape$dm <- ddfobj$shape$dm[rep(1, length(xgrid)), ,drop=FALSE]
      }

      # goofy workaround -- gamma is 0 at x=0, so add a little to the grid
      #  value so we don't have a weird drop
      if(ddfobj$type=="gamma"){
        xgrid[1] <- xgrid[1]+1e-6
      }

      if(point & pdf){
        ## calculate the pdf of distances
        # want r g(r) / int r g(r) dr

        # this is 2 r g(r)/w^2
        r_gr <- distpdf(xgrid, ddfobj, width=width, point=TRUE,
                        standardize=TRUE)
        # this is the value of int [2 r g(r) /w^2] dr
        int_r_gr <- integratepdf(ddfobj, select=TRUE, width=width,
                                 int.range=int.range, standardize=TRUE,
                                 point=TRUE)[1]
        # so the pdf values are:
        pdf_vals <- r_gr/int_r_gr

        # now rescale such that area under pdf == area under histogram
        linevalues <- pdf_vals * hist_area
      }else{
        linevalues <- detfct(xgrid, ddfobj, width=width)
      }
    }

    # actually plot the lines
    lines(xgrid, linevalues, col=byval1, ...)

    if(showpoints){
      jitter.p <- rnorm(length(point_vals), 1, jitter.v[1])
      points(xmat$distance, point_vals*jitter.p, col=byval1, ...)
    }

    # use the user-supplied main= ...
    if(!is.null(main)){
       title(main, cex.main=0.8)
    }
  }
  invisible()
}
