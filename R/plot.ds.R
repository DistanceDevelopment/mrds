#' Plot fit of detection functions and histograms of data from distance
#' sampling model
#'
#' Plots the fitted detection function(s) with a histogram of the observed distances to compare visually the fitted model and data.
#'
#' The structure of the histogram can be controlled by the user-defined
#' arguments \code{nc} or \code{breaks}. The observation specific detection
#' probabilities along with the line representing the fitted average detection
#' probability.
#'
#' It is not intended for the user to call \code{plot.ds} but its arguments
#' are documented here. Instead the generic \code{plot} command should be used
#' and it will call the appropriate function based on the class of the
#' \code{ddf} object.
#'
#' @aliases plot.ds
#' @export
#' @param x fitted model from \code{ddf}.
#' @param which index to specify which plots should be produced:
#'  \tabular{ll}{1 \tab histogram of observed distances\cr
#'               2 \tab histogram of observed distanes with fitted line and points (default)\cr}
#' @param breaks user define breakpoints
#' @param nc number of equal-width bins for histogram
#' @param jitter.v scaling option for plotting points.  Jitter is applied to points by multiplying the fitted value by a random draw from a normal distribution with mean 1 and sd \code{jitter.v[j]}.  Where \code{j=1,2} corresponds to observer \code{j} and \code{j=3} corresponds to pooled/duplicate detections.
#' @param showpoints logical variable; if \code{TRUE} plots predicted value for each observation.
#' @param subset subset of data to plot.
#' @param pl.col colours plotting colours for obs 1, obs 2 detections.
#' @param bw.col grayscale plotting colours for obs 1, obs 2 detections.
#' @param black.white logical variable; if \code{TRUE} plots are grayscale.
# @param pl.den shading density for plots of obs 1, obs 2 detections.
# @param pl.ang shading angle for plots of obs 1, obs 2 detections.
#' @param main user-specfied plot title.
#' @param pages the number of pages over which to spread the plots. For
#'  example, if \code{pages=1} then all plots will be displayed on one page.
#'  Default is 0, which prompts the user for the next plot to be displayed.
#'   (\code{\link{plot}}, \code{\link{hist}}, \code{\link{lines}},
#'   \code{\link{points}}, etc).
#' @param covlevels plot covariate level plots, in which case arguments \code{pages} and \code{which} are ignored.
#' @return Just plots.
#' @author Jeff Laake, Jon Bishop, David Borchers, David L Miller
#' @keywords plot
#' @importFrom graphics hist par lines points title
#' @importFrom grDevices devAskNewPage grey
#' @importFrom stats rnorm
#' @examples
#' \donttest{
#' data(book.tee.data)
#' egdata <- book.tee.data$book.tee.dataframe
#' xx <- ddf(dsmodel = ~mcds(key = "hn", formula = ~sex),
#'           data = egdata[egdata$observer==1, ],
#'           method = "ds", meta.data = list(width = 4))
#'
#' # not showing predicted probabilities
#' plot(xx,breaks=c(0,.5,1,2,3,4),showpoints=FALSE)
#'
#' # two subsets
#' plot(xx,breaks=c(0,.5,1,2,3,4),subset=sex==0)
#' plot(xx,breaks=c(0,.5,1,2,3,4),subset=sex==1)
#'
#' # put both plots on one page
#' plot(xx,breaks=c(0,.5,1,2,3,4),pages=1,which=1:2)
#' }
plot.ds <- function(x, which=2, breaks=NULL, nc=NULL,
                    jitter.v=rep(0,3), showpoints=TRUE, subset=NULL,
                    pl.col='black', bw.col=grey(0), black.white=FALSE,
                    #pl.den=0, pl.ang=rep(-45,1),
                    main=NULL, covlevels=TRUE,
                    pages=0, ...){

  #  Uses: setcov, detfct, histline, test.breaks
  model <- x
  show <- rep(FALSE, 2)
  show[which] <- TRUE
  #lower <- 0
  vname <- "distance"
  dpdata <- model$data

  xlab <- "Distance"

  dat<-dpdata

  # code from dpexplot:
  ltmodel <- model$ds
  width <- model$meta.data$width
  left <- model$meta.data$left
  ddfobj <- ltmodel$aux$ddfobj
  point <- ltmodel$aux$point
  if(is.null(ltmodel$aux$int.range)){
    int.range <- c(0,width)
  }else{
    int.range <- ltmodel$aux$int.range
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

  if(!is.null(substitute(subset))){
    selected <- eval(substitute(subset),ddfobj$xmat)
  }else{
    selected <- rep(TRUE,nrow(ddfobj$xmat))
  }

  if(all(!selected)){
    stop("Specified subset is empty.")
  }

  if(is.matrix(int.range)){
    int.range <- int.range[selected,]
  }

  xmat <- ddfobj$xmat[selected,]

  if(length(model$fitted)==1){
    pdot <- rep(model$fitted,sum(as.numeric(selected)))
  }else{
    pdot <- model$fitted[selected]
    Nhat <- sum(1/pdot)
  }

  n <- length(xmat$distance)

  if(!is.null(breaks)){
    nc <- length(breaks)-1
  }

  if(is.null(nc)){
    nc <- round(sqrt(n),0)
  }

  # Set logical hascov=TRUE when detection function has
  #  covariates other than distance and observer
  hascov <- !ddfobj$intercept.only

  # Compute a grid for distance (xgrid)
  xgrid <- (width/100)*(seq(0:100)-1)

  # create intervals of distance (breaks) for the chosen number of classes (nc).
  if(is.null(breaks)){
    if(!is.null(model$meta.data$binned) && model$meta.data$binned){
      breaks <- ltmodel$aux$breaks
      nc <- length(breaks)-1
    }else{
      breaks <- c(max(0,(max.range[1])),
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

  dat <- dpdata[selected,]
  keep <- dat[,vname]>=min(breaks) & dat[,vname]<=max(breaks)
  h1 <- hist(dat[,vname][keep],breaks=breaks,plot=FALSE)
  ymax <- max(h1$counts)

  # Set printing options for plots:
  # By default  pl.col=c('black'),
  #             bw.col=c(grey(0))

  # If greyscale plots are required use the following colours:
  if(black.white){
    byval1 <- bw.col[1]
  }else{
    # If colour plots are required use the following:
    byval1 <- pl.col[1]
  }

  # Density of shaded lines - default is set all to 0
  #denval1 <- pl.den[1]
  # Angle of shaded lines - default is set all to -45
  #angval1 <- pl.ang[1]

  # Scaling for grouped data:
  if(normalize & !point){
    bindata <- function(x,r,breaks){
      return(hist(r[r>=x[1]&r<=x[2]],breaks=breaks,plot=FALSE)$counts)
    }
    sumit<-function(x,n,wt){
      return(sum(x/(wt*n)))
    }
    expected.counts <- apply(int.range,1,bindata,
                             r=(0:1000)*width/1001,breaks=breaks)
    expected.counts <- apply(expected.counts,1,sumit,n=1001,wt=pdot)
  }else{
    if(!point){
      expected.counts <- (breaks[2:(nc+1)]-breaks[1:nc])*(Nhat/breaks[nc+1])
    }else{
      expected.counts <- -apply(matrix(c(breaks[2:(nc+1)]^2,breaks[1:nc]^2),
                                       ncol=2,nrow=nc),
                                1,diff)*(Nhat/breaks[nc+1]^2 )
    }
  }

  hist.obj <- hist(dat[,vname][keep], breaks=breaks, plot=FALSE)
  hist.obj$density <- hist.obj$counts/expected.counts
  hist.obj$density[expected.counts==0] <- 0
  hist.obj$equidist <- FALSE


  ### Actual plotting starts here

  if(!hascov) covlevels <- FALSE

  if(!covlevels){
    # do the paging, using devAskNewPage() means we can just call plots and
    # R will make the prompt for us
    if(pages!=1 & sum(show)!=1){
      oask <- devAskNewPage(TRUE)
      on.exit(devAskNewPage(oask))
    }else if(sum(show)!=1){
      opar <- par(mfrow=c(1,sum(show)))
      on.exit(par(opar))
    }

    # plot histogram of distances
    if(show[1]){
      plot(h1, xlab=xlab, ylab="Frequency",# yaxp=c(0,1,5),
           main="Histogram of distances",freq=TRUE,...)
      box()
    }

    # Detection function plot overlaid on histogram of observed
    # distances - much of this code is taken from earlier plot.ds function.
    if(show[2]){
      # Primary detection function
      gxvalues <- detfct(xmat$distance,ddfobj,select=selected,width=width)
      ymax <- max(hist.obj$density,max(gxvalues))
      plot(hist.obj, ylim=c(0,ymax), xlab=xlab, yaxp=c(0,1,5),
           main="Detection function",freq=FALSE,ylab="Detection probability",...)
      box()

      # when we have covariates
      if(hascov){

        linevalues <- plot_avg_df(xgrid, width, xmat, ddfobj, selected,
                                  int.range, pdot, normalize, range.varies)

      # without covariates
      }else{
        if(!is.null(ddfobj$scale)){
          ddfobj$scale$dm <- ddfobj$scale$dm[rep(1,length(xgrid)),]
        }
        if(!is.null(ddfobj$shape)){
           ddfobj$shape$dm <- ddfobj$shape$dm[rep(1,length(xgrid)),]
        }

        # goofy workaround -- gamma is 0 at x=0, so add a little to the grid
        #  value so we don't have a weird drop
        if(ddfobj$type=="gamma"){
          xgrid[1] <- xgrid[1]+1e-6
        }
        linevalues <- detfct(xgrid,ddfobj,width=width)
      }

      # actually plot the lines
      lines(xgrid,linevalues,col=byval1,...)

      if(showpoints){
        jitter.p <- rnorm(length(gxvalues),1,jitter.v[1])
        points(xmat$distance,gxvalues*jitter.p,col=byval1,...)
      }
    }

    # use the user-supplied main= ...
    if(!is.null(main)){
       title(main,cex.main=0.8)
    }

  }else{

    # plot the covariate levels
    cov_levels <- ds_cov_levels(model, xgrid)
    #lapply(cov_levels, lines, x=xgrid)
    df_line <- cov_levels$df_line

    plot_all <- cov_levels$plot_all

    # set graphics pars
    opar <- par(mfrow=c(1,length(plot_all)+1))
    on.exit(par(opar))

    if(!is.null(main)){
      # this seems stupid but it does mean that folks can write
      # main="" and it blanks all titles
      main <- rep_len(main, length(plot_all)+1)
    }else{
      main <- paste0("Average detection function")
    }

    # first make the average detection function plot
    plot(hist.obj, freq=FALSE, xlab=xlab, ylab="Detection probability", main=main[1], ...)
    lines(x=xgrid, y=plot_avg_df(xgrid, width, xmat, ddfobj, selected,
                                 int.range, pdot, normalize,range.varies))
    box()


    ## actually do some plotting
    for(cov_i in seq_along(plot_all)){
      if(is.null(main) || length(main)==1){
        plot_title <- paste(ifelse(is.factor(plot_all[[cov_i]][[1]]),
                                   "Levels of", "Quantiles of"),
                            names(plot_all)[cov_i])
      }else{
        plot_title <- main[cov_i+1]
      }
      # make the plot
      plot(hist.obj, freq=FALSE, xlab=xlab, ylab="Frequency", main=plot_title, ...)

      # colours/styles for the lines
      # as long as there are <=6 levels to plot, use different styles
      # when > then use grey scales
      if(length(plot_all[[cov_i]])>6){
        cols <- grey(seq(0.1,0.9,len=length(plot_all[[cov_i]])))
        ltys <- rep(1,length(plot_all[[cov_i]]))
      }else{
        cols <- rep("black",length(plot_all[[cov_i]]))
        ltys <- 1:length(plot_all[[cov_i]])
      }

      ii <- 1
      # now iterate over their values
      for(val_j in seq_along(plot_all[[cov_i]])){
        lines(x=xgrid, y=df_line[[ii]], col=cols[val_j], lty=ltys[val_j])
        ii <- ii+1
      }

      # get the legend levels labels
      if(is.factor(plot_all[[cov_i]])){
        legend.levels <- plot_all[[cov_i]]
      }else{
        legend.levels <- as.character(round(plot_all[[cov_i]],3))
      }

      # plot legend
      legend("topright", legend=legend.levels,
             title=names(plot_all)[cov_i],
             col=cols, lty=ltys, inset=c(0.01,0.01))

      # put a box around it
      box()

    }
  }


  invisible()
}
