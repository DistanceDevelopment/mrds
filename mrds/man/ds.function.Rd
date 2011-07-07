\name{ds.function}
\alias{ds.function}
\title{Distance Sampling Functions...}
\usage{ds.function(model, newdata, obs="All", conditional=FALSE,
    pdf=TRUE, finebr)}

\description{
  Distance Sampling Functions Placeholder -- Not functional
  ----
}
\value{List containing \item{xgrid }{grid of distance values} \item{values }{average detection fct values at the xgrid values}}
\author{Jeff Laake}
\arguments{
  \item{model}{model object}
  \item{newdata}{dataframe at which to compute values; if NULL uses fitting data}
  \item{obs}{1 or 2 for observer 1 or 2, 3 for duplicates, "." for combined and "All" to return all of the values}
  \item{conditional}{if FALSE, computes p(x) based on distance detection function and if TRUE based on mr detection function}
  \item{pdf}{if FALSE, returns p(x) and if TRUE, returns p(x)*pi(x)/integral p(x)*pi(x)}
  \item{finebr}{fine break values over which line is averaged}
}
