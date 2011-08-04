#' Print summary of distance detection function model object
#' 
#' Provides a brief summary of data and fitted detection probability model
#' parameters, model selection criterion, and optionally abundance in the
#' covered (sampled) region and its standard error. What is printed depends
#' on the corresponding call to summary.
#' 
#' @aliases print.summary.rem
#' @method print summary.rem
#' @S3method print summary.rem
#' @param x a summary of \code{ddf} model object
#' @param \dots unspecified and unused arguments for S3 consistency
#' @return NULL
#' @author Jeff Laake
#' @seealso \code{\link{summary.rem}}
#' @keywords utility
print.summary.rem <-
function(x,...)
{
#
# print.summary.rem
#
# Provides a summary of parameters and estimates from the output of rem object
#
# Arguments:
#
# x      - list object from summary.rem
#
# Value: Null
#
  print(x$mr.summary)
  cat("\n\n")
  print(x$ds.summary)
  cat("\n\nSummary for rem object\n")
  cat("\nTotal AIC value = ",x$AIC,"\n")

  parameters=data.frame(Estimate=c(x$average.p,x$Nhat))
  row.names(parameters)=c("Average p", "N in covered region")
  if(!is.null(x$average.p.se))
  {
      parameters$SE=c(x$average.p.se,x$Nhat.se)
      parameters$CV=parameters$SE/parameters$Estimate
  }
  print(parameters)
  invisible(NULL)
}
