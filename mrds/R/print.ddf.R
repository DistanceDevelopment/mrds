#' Simple pretty printer for distance sampling analyses 
#'
#' Simply prints out a brief description of the model which was fitted. For more
#' detailed information use \code{\link{summary}}.
#'
#' @param x a \code{ddf} object 
#' @param ... not passed through, just for S3 compatibility.
#' @S3method print ddf
#' @method print ddf 
#' @aliases print.ddf
#'
#' @author David L. Miller
#' @export
print.ddf<-function(x, ...){

  cat("\nDistance sampling analysis object\n")
  cat("\nDetection function:\n",model.description(x),"\n")

  if(!is.null(x$ds$aux$ddfobj$scale$formula)){
    cat("\nModel formula:",x$ds$aux$ddfobj$scale$formula,"\n")
  }

  coeff<-coef(x)

  cat("\nModel coefficients:\n")
  print(coeff$scale)
  if(x$ds$aux$ddfobj$type %in% c("gamma","hr")) {
    cat("\nShape parameters: ", "\n")
    print(coeff$exponent)
  }
  if (!is.null(coeff$adj.parm)) {
     cat("\nAdjustment term parameter(s): ", "\n")
     print(coeff$adjustment)
  }

  # Remind the user that monotonicity constraints were enforced
  if(x$ds$aux$mono & x$ds$aux$mono.strict){
    cat("\nStrict monotonicity constraints were enforced.\n")
  }else if(x$ds$aux$mono){
    cat("\nMonotonicity constraints were enforced.\n")
  }


  if(!is.null(x$Nhat)){
    cat("\nEstimated abundance in covered region:",x$Nhat,"\n")
  }

  cat("\nAIC:",x$criterion,"\n")

  invisible()
}
