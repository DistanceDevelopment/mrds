\name{histline}
\alias{histline}
\title{Plot histogram line}
\usage{
histline(height, breaks, lineonly = FALSE, outline = FALSE, fill = FALSE,
  ylim = range(height), xlab = "x", ylab = "y", det.plot = FALSE, ...)
}
\arguments{
  \item{height}{heights of histogram bars}

  \item{breaks}{cutpoints for x}

  \item{lineonly}{if TRUE, drawn with plot; otherwise with
  lines to allow addition of current plot}

  \item{outline}{if TRUE, only outline of histogram is
  plotted}

  \item{fill}{If fill==TRUE, uses polygon() to fill bars}

  \item{ylim}{limits for y axis}

  \item{xlab}{label for x axis}

  \item{ylab}{label for y axis}

  \item{det.plot}{if TRUE, plot is of detection so yaxis
  limited to unit interval}

  \item{\dots}{Additional unspecified arguments for plot
  (fill==TRUE}
}
\value{
None
}
\description{
Takes bar heights (height) and cutpoints (breaks), and
constructs a line-only histogram from them using the
function plot() (if lineonly==FALSE) or lines() (if
lineonly==TRUE).
}
\author{
???
}

