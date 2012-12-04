\name{plot.rem.fi}
\alias{plot.rem.fi}
\title{Plot fit of detection functions and histograms of data from distance
sampling model}
\usage{
  \method{plot}{rem.fi} (x, which = 1:6, breaks = NULL,
    nc = NULL, maintitle = "", showlines = TRUE,
    showpoints = TRUE, ylim = c(0, 1), angle = -45,
    density = 20, col = "black", jitter = NULL,
    divisions = 25, new = TRUE, ...)
}
\arguments{
  \item{x}{fitted model from \code{ddf}}

  \item{which}{index to specify which plots should be
  produced. 1: uncond det fct, 2:cond det fct}

  \item{breaks}{user define breakpoints}

  \item{nc}{number of equal-width bins for histogram}

  \item{maintitle}{main title line for each plot}

  \item{showpoints}{logical variable; if TRUE plots
  predicted value for each observation}

  \item{showlines}{logical variable; if TRUE a line
  representing the average detection probability is
  plotted}

  \item{ylim}{range of y axis; defaults to (0,1)}

  \item{angle}{shading angle for hatching}

  \item{density}{shading density for hatching}

  \item{col}{plotting colour}

  \item{jitter}{scaling option for plotting points.  Jitter
  is applied to points by multiplying the fitted value by a
  random draw from a normal distribution with mean 1 and sd
  jitter.}

  \item{divisions}{number of divisions for averaging line
  values; default = 25}

  \item{new}{if TRUE, opens new device for each plot; set
  new=FALSE if you use par(mfrow=..) or layout}

  \item{\dots}{other graphical parameters, passed to the
  plotting functions (plot, hist, lines, points, etc)}
}
\value{
  NULL
}
\description{
  Plots the fitted detection functions for a distance
  sampling model and histograms of the distances (for
  unconditional detection functions) or proportion of
  observations detected within distance intervals (for
  conditional detection functions) to compare visually the
  fitted model and data.
}
\details{
  The structure of the histogram can be controlled by the
  user-defined arguments \code{nc} or \code{breaks}.  The
  observation specific detection probabilities along with
  the line representing the fitted average detection
  probability.

  In earlier versions of this package the plot routines
  other than \code{plot.ds} (\code{plot.io} etc) used
  \code{plot.cond} and \code{plot.uncond} to construct the
  conditional and unconditional detection functions
  respectively.  This is no longer the case although the
  advanced user can still call \code{plot.cond} and
  \code{plot.uncond} if required.  It is not intended for
  the user to call any of \code{plot.ds},
  \code{plot.trial.fi},
  \code{plot.trial},\code{plot.rem.fi}, \code{plot.rem},
  \code{plot.io.fi} or \code{plot.io} but the arguments are
  documented here. Instead the generic \code{plot} command
  should be used and it will call the appropriate function
  based on the type of \code{ddf} object.

  For plot routine \code{plot.ds} the \code{which} command
  allows the user to select which plots are returned from
  the following options: \tabular{ll}{ \code{which} \tab
  \code{plot}\cr \code{1} \tab data summary plot - a
  histogram of the observed distances \cr \code{2} \tab a
  scaled histogram of detections with a line giving the
  detection function averaged over the estimated population
  levels of the covariate values, and one dot for each
  observation at its estimated detection probability.\cr }

  For plot routines \code{plot.trial.fi} and
  \code{plot.trial} the following plots are available:
  \tabular{ll}{ \code{which} \tab \code{plot} \cr \code{1}
  \tab data summary plot - a histogram of the observed
  distances for observer 1 \cr \code{2} \tab data summary
  plot - a histogram of the observed distances for observer
  2 \cr \code{3} \tab Observer 1 detection function - a
  scaled histogram of detections with fitted DS model
  scaled from the MR estimated g(0).  The line shows the
  population average detection function and the points
  display estimated detection probability \cr \code{4} \tab
  Conditional MR detection function - observer 1 given obs
  2, giving the proportion of duplicates with fitted MR
  model averaged over population covariate vales and dots
  for each estimated detection probability. \cr } For plot
  routines \code{plot.rem.fi} and \code{plot.rem} the
  following plots are available: \tabular{ll}{ \code{which}
  \tab \code{plot} \cr \code{1} \tab data summary plot - a
  histogram of the observed distances for observer 1 \cr
  \code{2} \tab data summary plot - a histogram of the
  observed distances for observer 2 \cr \code{3} \tab
  Observer 1 detection function - a scaled histogram of
  detections with fitted DS model scaled from the MR
  estimated g(0).  The line shows the population average
  detection function and the points display estimated
  detection probability \cr \code{4} \tab Conditional MR
  detection function - observer 1 given obs 2, giving the
  proportion of duplicates with fitted MR model averaged
  over population covariate vales and dots for each
  estimated detection probability. \cr }

  For plot routines \code{plot.io.fi} and \code{plot.io}
  the following plots are available: \tabular{ll}{
  \code{which} \tab \code{plot} \cr \code{1} \tab data
  summary plot - a histogram of the observed distances for
  observer 1 \cr \code{2} \tab data summary plot - a
  histogram of the observed distances for observer 2 \cr
  \code{3} \tab Observer 1 detection function - a scaled
  histogram of detections with fitted DS model scaled from
  the MR estimated g(0).  The line shows the population
  average detection function and the points display
  estimated detection probability \cr \code{4} \tab
  Observer 2 detection function - as for \code{plot} 3 but
  using the detections from Observer 2\cr \code{5} \tab
  Duplicates detection function - as for \code{plot} 3 but
  using the duplicate detections\cr \code{6} \tab Pooled
  detection function - as for \code{plot} 3 but using the
  pooled detections\cr \code{7} \tab Conditional MR
  detection function - observer 1 given obs 2, giving the
  proportion of duplicates with fitted MR model averaged
  over population covariate vales and dots for each
  estimated detection probability. \cr \code{8} \tab
  Conditional MR detection function - observer 2 given obs
  1, giving the proportion of duplicates with fitted MR
  model averaged over population covariate vales and dots
  for each estimated detection probability. \cr }
}
\examples{
\donttest{
data(book.tee.data)
region<<-book.tee.data$book.tee.region
egdata<<-book.tee.data$book.tee.dataframe
samples<<-book.tee.data$book.tee.samples
obs<<-book.tee.data$book.tee.obs
xx=ddf(dsmodel = ~mcds(key = "hn", formula = ~sex), data = egdata[egdata$observer==1, ], method = "ds", meta.data = list(width = 4))
par(mfrow=c(2,2))
plot(xx,breaks=c(0,.5,1,2,3,4),showpoints=FALSE,new=FALSE)
plot(xx,breaks=c(0,.5,1,2,3,4),subset=sex==0,new=FALSE)
plot(xx,breaks=c(0,.5,1,2,3,4),subset=sex==1,new=FALSE)
}
}
\author{
  Jeff Laake, Jon Bishop, David Borchers
}
\keyword{plot}

