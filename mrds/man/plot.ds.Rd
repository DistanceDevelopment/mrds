\name{plot.ds}
\alias{plot.ds}
\title{Plot fit of detection functions and histograms of data from distance
sampling model}
\usage{
  \method{plot}{ds} (x, which = c(2), byvar = "",
    breaks = NULL, nc = NULL, winht = 4, winwd = 6,
    jitter.v = rep(0, 3), showpoints = TRUE, subset = NULL,
    pl.col = c("black"), bw.col = c(grey(0)),
    black.white = FALSE, pl.den = rep(20, 1),
    pl.ang = rep(-45, 1), main = NULL, ...)
}
\arguments{
  \item{x}{fitted model from \code{ddf}}

  \item{which}{index to specify which plots should be
  produced. 1: histogram, 2: histogram with fitted line and
  points.}

  \item{byvar}{name of variable to be used to color points
  - not currently implemented.}

  \item{breaks}{user define breakpoints}

  \item{nc}{number of equal-width bins for histogram}

  \item{winht}{plot window height (not currently
  implemented)}

  \item{winwd}{plot window width (not currently
  implemented)}

  \item{jitter.v}{scaling option for plotting points.
  Jitter is applied to points by multiplying the fitted
  value by a random draw from a normal distribution with
  mean 1 and sd jitter.v[j].  Where j=1,2 corresponds to
  observer j and j=3 corresponds to pooled/duplicate
  detections.}

  \item{showpoints}{logical variable; if TRUE plots
  predicted value for each observation}

  \item{subset}{subset of data to plot}

  \item{pl.col}{colours plotting colours for obs 1, obs 2
  detections}

  \item{bw.col}{grayscale plotting colours for obs 1, obs 2
  detections}

  \item{black.white}{logical variable; if TRUE plots are
  grayscale}

  \item{pl.den}{shading density for plots of obs 1, obs 2
  detections}

  \item{pl.ang}{shading angle for plots of obs 1, obs 2
  detections}

  \item{main}{user-specfied plot title}

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

  It is not intended for the user to call any of
  \code{plot.ds}, \code{plot.trial.fi},
  \code{plot.trial},\code{plot.rem.fi}, \code{plot.rem},
  \code{plot.io.fi} or \code{plot.io} but the arguments are
  documented here. Instead the generic \code{plot} command
  should be used and it will call the appropriate function
  based on the type of \code{ddf} object.

  The \code{which} command allows the user to select which
  plots are returned. See which argument definition.
}
\examples{
\donttest{
data(book.tee.data)
region<<-book.tee.data$book.tee.region
egdata<<-book.tee.data$book.tee.dataframe
samples<<-book.tee.data$book.tee.samples
obs<<-book.tee.data$book.tee.obs
xx=ddf(dsmodel = ~mcds(key = "hn", formula = ~sex), data = egdata[egdata$observer==1, ],
           method = "ds", meta.data = list(width = 4))
plot(xx,breaks=c(0,.5,1,2,3,4),showpoints=FALSE)
plot(xx,breaks=c(0,.5,1,2,3,4),subset=sex==0)
plot(xx,breaks=c(0,.5,1,2,3,4),subset=sex==1)
}
}
\author{
  Jeff Laake, Jon Bishop, David Borchers
}
\keyword{plot}

