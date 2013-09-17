\name{plot.io.fi}
\alias{plot.io.fi}
\title{Plot fit of detection functions and histograms of data from distance
sampling independent observer model}
\usage{
  \method{plot}{io.fi} (x, which = 1:6, breaks = NULL,
    nc = NULL, maintitle = "", showlines = TRUE,
    showpoints = TRUE, ylim = c(0, 1), angle = -45,
    density = 20, col = "black", jitter = NULL,
    divisions = 25, new = TRUE, xlab = "Distance",
    ylab = "Detection probability", subtitle = TRUE, ...)
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

  \item{xlab}{label for x-axis}

  \item{ylab}{label for y-axis}

  \item{subtitle}{if TRUE, shows plot type as sub-title}

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
\author{
  Jeff Laake, Jon Bishop, David Borchers
}
\keyword{plot}

