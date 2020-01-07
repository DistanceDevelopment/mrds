#' Invert of covariance matrices
#'
#' Tries to invert a matrix by \code{solve}. If this fails because of singularity, an eigenvector decomposition is computed, and eigenvalues below \code{1/cmax} are replaced by \code{1/cmax}, i.e., \code{cmax} will be the corresponding eigenvalue of the inverted matrix.
#' @param m a numeric symmetric matrix.
#' @param cmax a positive value, see above.
#' @return A list with the following components: \code{inv} the inverted matrix, \code{coll} \code{TRUE} if \code{solve} failed because of singularity.
#' @author Christian Hennig \url{http://www.homepages.ucl.ac.uk/~ucakche/}
#' @seealso solve, eigen
#' @section Source:
#' \code{solvecov} code was taken from package \code{fpc}: Christian Hennig \url{http://www.homepages.ucl.ac.uk/~ucakche/}
#' @export
solvecov <- function (m, cmax = 1e+10){
  options(show.error.messages = FALSE)
  covinv <- try(solve(m))
  if("try-error" %in% class(covinv)){
    coll <- FALSE
  }else{
    p <- nrow(m)
    cove <- eigen(m, symmetric = TRUE)
    coll <- TRUE
    if(min(cove$values) < 1/cmax) {
      covewi <- diag(p)
      for (i in 1:p){
        if (cove$values[i] < 1/cmax){
          covewi[i, i] <- cmax
        }else{
         covewi[i, i] <- 1/cove$values[i]
        }
      }
    }else{
      covewi <- diag(1/cove$values, nrow = length(cove$values))
    }
    covinv <- cove$vectors %*% covewi %*% t(cove$vectors)
  }
  options(show.error.messages = TRUE)
  out <- list(inv = covinv, coll = coll)
  return(out)
}
