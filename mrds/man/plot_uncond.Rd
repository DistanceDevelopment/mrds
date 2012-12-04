\name{plot_uncond}
\alias{plot_uncond}
\title{Plot unconditional detection function from distance sampling model}
\usage{
  plot_uncond(model, obs, xmat, gxvalues, nc, finebr,
    breaks, showpoints, showlines, maintitle, ylim,
    point = FALSE, return.lines = FALSE, angle = -45,
    density = 20, col = "black", jitter = NULL, ...)
}
\arguments{
  \item{model}{fitted model from \code{ddf}}

  \item{obs}{value of observer for plot}

  \item{xmat}{processed data}

  \item{gxvalues}{detection function values for each
  observation}

  \item{nc}{number of equal-width bins for histogram}

  \item{finebr}{fine break values over which line is
  averaged}

  \item{breaks}{user define breakpoints}

  \item{showpoints}{logical variable; if TRUE plots
  predicted value for each observation}

  \item{showlines}{logical variable; if TRUE plots average
  predicted value line}

  \item{maintitle}{main title line for each plot}

  \item{ylim}{range of y axis; defaults to (0,1)}

  \item{point}{point count data if TRUE}

  \item{return.lines}{if TRUE, returns values for line}

  \item{angle}{shading angle for hatching}

  \item{density}{shading density for hatching}

  \item{col}{plotting colour}

  \item{jitter}{scaling option for plotting points.  Jitter
  is applied to points by multiplying the fitted value by a
  random draw from a normal distribution with mean 1 and sd
  jitter.}

  \item{\dots}{other graphical parameters, passed to the
  plotting functions (plot, hist, lines, points, etc)}
}
\value{
  NULL
}
\description{
  Plots unconditional detection function for observer=obs
  observations overlays histrogram, average detection
  function and values for individual observations data.
}
\examples{
\donttest{
data(book.tee.data)
region<<-book.tee.data$book.tee.region
egdata<<-book.tee.data$book.tee.dataframe
samples<<-book.tee.data$book.tee.samples
obs<<-book.tee.data$book.tee.obs
xx=ddf(dsmodel = ~mcds(key = "hn", formula = ~sex), data = egdata[egdata$observer==1, ], method = "ds", meta.data = list(width = 4))
plot(xx,breaks=c(0,.5,1,2,3,4),showpoints=FALSE)
plot(xx,breaks=c(0,.5,1,2,3,4),subset=sex==0)
plot(xx,breaks=c(0,.5,1,2,3,4),subset=sex==1)
}
}
\author{
  Jeff Laake, Jon Bishop, David Borchers
}
\keyword{plot}

