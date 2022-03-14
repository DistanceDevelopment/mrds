#' Plot fit of detection functions and histograms of data from distance
#' sampling model
#'
#' Plots the fitted detection function(s) with a histogram of the observed
#' distances to compare visually the fitted model and data.
#'
#' The structure of the histogram can be controlled by the user-defined
#' arguments \code{nc} or \code{breaks}. The observation specific detection
#' probabilities along with the line representing the fitted average detection
#' probability.
#'
#' It is not intended for the user to call \code{plot.ds} but its arguments are
#' documented here. Instead the generic \code{plot} command should be used and
#' it will call the appropriate function based on the class of the \code{ddf}
#' object.
#'
#' @aliases plot.ds
#' @export
#' @param x fitted model from \code{ddf}.
#' @param which index to specify which plots should be produced:
#'  \tabular{ll}{1 \tab histogram of observed distances\cr
#'               2 \tab histogram of observed distances with fitted line and
#'               points (default)\cr}
#' @param breaks user defined breakpoints
#' @param nc number of equal width bins for histogram
#' @param jitter.v apply jitter to points by multiplying the fitted value by a
#' random draw from a normal distribution with mean 1 and sd \code{jitter.v}.
#' @param showpoints logical variable; if \code{TRUE} plots predicted value for
#' each observation (conditional on its observed distance).
#' @param subset subset of data to plot.
#' @param pl.col colour for histogram bars.
#' @param pl.den shading density for histogram bars.
#' @param pl.ang shading angle for histogram bars.
#' @param main plot title.
#' @param ylim vertical axis limits.
#' @param pdf plot the histogram of distances with the PDF of the probability
#' of detection overlaid. Ignored (with warning) for line transect models.
#' @param pages the number of pages over which to spread the plots. For
#' example, if \code{pages=1} then all plots will be displayed on one page.
#' Default is 0, which prompts the user for the next plot to be displayed.
#' @param xlab horizontal axis label (defaults to "Distance").
#' @param ylab vertical axis label (default automatically set depending on plot
#' type).
#' @param \dots other graphical parameters, passed to the plotting functions
#' (\code{\link{plot}}, \code{\link{hist}}, \code{\link{lines}},
#' \code{\link{points}}, etc).
#' @return Just plots.
#' @author Jeff Laake, Jon Bishop, David Borchers, David L Miller
#' @keywords plot
#' @importFrom graphics hist par lines points title
#' @importFrom grDevices devAskNewPage grey
#' @importFrom stats rnorm
#' @seealso add_df_covar_line
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
                    pl.col="lightgrey",
                    pl.den=NULL, pl.ang=NULL, main=NULL, pages=0,
                    pdf=FALSE, ylim=NULL, xlab="Distance", ylab=NULL, ...){

  model <- x
  lower <- 0
  vname <- "distance"
  dat <- model$data

  # ignore pdf=TRUE with line transect data or with gamma df
  if((pdf & !model$meta.data$point) |
     (pdf & model$ds$aux$ddfobj$type=="gamma")){
    warning("Ignoring pdf=TRUE for line transect data")
    pdf <- FALSE
  }

  # decide which plots to show
  show <- rep(FALSE, 2)
  show[which] <- TRUE

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
    selected <- eval(substitute(subset), ddfobj$xmat)
  }else{
    selected <- rep(TRUE, nrow(ddfobj$xmat))
  }
  # die if there was nothing selected
  if(all(!selected)){
    stop("Specified subset is empty.")
  }

  if(is.matrix(int.range)){
    int.range <- int.range[selected, ]
  }

  xmat <- ddfobj$xmat[selected,]
  if(!is.null(ddfobj$scale)){
    z <- ddfobj$scale$dm[selected, , drop=FALSE]
  }else{
    z <- matrix(1, nrow=1, ncol=1)
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
    zgrid <- matrix(rep(z[1,], length(xgrid)), byrow=TRUE, ncol=sum(zdim))
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
        breaks <- c(left, breaks)
        nc <- nc+1
      }
    }
  }

  # test breaks for validity and reset as needed
  breaks <- test.breaks(breaks, model$meta.data$left, width)
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
    if(is.null(ylab)) ylab <- "Frequency"
    if(is.null(ylim)) ylim <- c(0, ymax)
    histline(hist.obj$counts, breaks=breaks, lineonly=FALSE, ylim=ylim,
             xlab=xlab, ylab=ylab, angle=angval1,
             density=denval1, col=pl.col, ...)
  }

  ## Detection function plot overlaid on histogram of observed distances
  if(show[2]){

    # area under the histogram
    hist_area <- sum(hist.obj$density*diff(breaks))
    # Detection function/pdf values for points to be plotted
    if(point & pdf){
      point_vals <- distpdf(xmat$distance, ddfobj, width=width, point=point,
                                standardize=TRUE)/
                    integratepdf(ddfobj, select=selected, width=width,
                                 int.range=int.range, standardize=TRUE,
                                 point=point)
    }else{
      point_vals <- detfct(xmat$distance, ddfobj, select=selected, width=width)
    }

    # set y labels, limits and tick marks (det.plot) depending on if we
    # are plotting PDF or df
    if(is.null(ylim)) ylim<-c(0, max(hist.obj$density, max(point_vals)))
    if(pdf){
      if(is.null(ylab)) ylab <- "Probability density"
      det.plot <- FALSE
    }else{
      if(is.null(ylab)) ylab <- "Detection probability"
      det.plot <- TRUE
    }

    # plot the histogram
    histline(hist.obj$density, breaks=breaks, lineonly=FALSE,
             xlab=xlab, ylab=ylab, ylim=ylim,
             angle=angval1, density=denval1, col=pl.col,
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

        if((point & pdf) | ddfobj$type=="gamma"){
          ## calculate the pdf of distances
          # want r g(r) / int r g(r) dr

          # this is 2 r g(r)/w^2
          r_gr <- distpdf(newdat$distance, ddfobj, width=width, point=point,
                          standardize=TRUE)
          # this is the value of int [2 r g(r) /w^2] dr
          int_r_gr <- integratepdf(ddfobj, select=selected, width=width,
                                   int.range=int.range, standardize=TRUE,
                                   point=point)
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
      if(ddfobj$type=="gamma" & left==0){
        xgrid[1] <- xgrid[1]+1e-6
      }

      if((point & pdf) | ddfobj$type=="gamma"){
        ## calculate the pdf of distances
        # want r g(r) / int r g(r) dr

        # this is 2 r g(r)/w^2
        r_gr <- distpdf(xgrid, ddfobj, width=width, point=point,
                        standardize=TRUE, left=left)
        # this is the value of int [2 r g(r) /w^2] dr
        int_r_gr <- integratepdf(ddfobj, select=TRUE, width=width,
                                 int.range=int.range, standardize=TRUE,
                                 point=point, left=left)[1]
        # so the pdf values are:
        pdf_vals <- r_gr/int_r_gr

        # now rescale such that area under pdf == area under histogram
        linevalues <- pdf_vals * hist_area
      }else{
        linevalues <- detfct(xgrid, ddfobj, width=width, left=left)
      }
    }

    # remove elements of ... which are used in the histogram but don't
    # make sense for the lines or points next
    dots <- list(...)
    dots$border <- NULL

    # actually plot the lines
    ldots <- dots
    ldots$x <- xgrid
    ldots$y <- linevalues
    do.call(lines, ldots)

    if(showpoints){
      jitter.p <- rnorm(length(point_vals), 1, jitter.v[1])
      pdots <- dots
      pdots$x <- xmat$distance
      pdots$y <- point_vals*jitter.p
      #points(xmat$distance, point_vals*jitter.p, ...)
      do.call(points, pdots)
    }

    # use the user-supplied plot title
    if(!is.null(main)){
       title(main, cex.main=0.8)
    }
  }
  invisible()
}
