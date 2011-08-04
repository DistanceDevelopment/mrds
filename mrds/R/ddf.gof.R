#' Goodness of fit tests for distance sampling models
#' 
#' Generic function that computes chi-square goodness of fit test for ddf
#' models
#' 
#' 
#' @aliases ddf.gof gof.ds gof.io gof.io.fi gof.trial gof.trial.fi gof.rem
#'   gof.rem.fi
#' @export
#' @param model ddf model object
#' @param breaks Cutpoints to use for binning data
#' @param nc Number of distance classes
#' @param qq Flag to indicate whether quantile-quantile plot is desired
#' @param \dots Graphics parameters to pass into qqplot function
#' @return List of class 'ddf.gof' containing \item{chi-square }{Goodness of
#'   fit test statistic} \item{df }{Degrees of freedom associated with test
#'   statistic} \item{p-value }{Significance level of test statistic}
#' @author Jeff Laake
#' @seealso \code{\link{qqplot.ddf}}
#' @keywords utility
ddf.gof <-
function(model,breaks=NULL,nc=NULL,qq=TRUE,...)
#
# gof - generic function that computes chi-square gof test for ddf models
#
# Arguments: 
#
# model  - ddf model object
# breaks - distance cut points
# nc     - number of distance classes
# ...    - graphics params to pass into qqplot
#
# Value:
# 
# result - list with chi-square value, df and p-value
#
# Functions Used: gof.ds, gof.io, gof.io.fi, gof.trial, gof.trial.fi, qqplot.df
#
{
#
# call method specific function
#
  if(!is.null(breaks)){
    breaks=test.breaks(breaks,model$meta.data$left,model$meta.data$width)
    nc=length(breaks)-1
  }
  result=switch(model$method,
         ds=gof.ds(model,breaks,nc),
         io=gof.io(model,breaks,nc),
         io.fi=gof.io.fi(model,breaks,nc),
         trial=gof.trial(model,breaks,nc),
         trial.fi=gof.trial.fi(model,breaks,nc),
         rem=gof.rem(model,breaks,nc),
         rem.fi=gof.rem.fi(model,breaks,nc))

  if(qq & !model$meta.data$binned)
     result=list(chisquare=result,dsgof=qqplot.ddf(model,...))
  else
     result=list(chisquare=result)
  
  class(result)=c("ddf.gof")
  return(result)
}
