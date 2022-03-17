#' @importFrom methods is
# is what we have a formula, make it if you can, else error
fixformula <- function(formula){
  if(!is(formula, "formula")){
    if(!is(try(as.formula(formula)), "formula")){
      stop("Invalid formula")
    }
  }else{
    formula <- paste(as.character(formula), collapse="")
  }

  return(formula)
}
