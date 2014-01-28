#' Hermite polynomial generating function
#'
#' @param z points to evaluate
#' @param n order Hermite polynomial
#'
#' Using formula from http://functions.wolfram.com/Polynomials/HermiteH/02/0001/
#' This may not be the fastest way of finding the value.
hermite.poly <- function(z,n){

  end <- floor(n/2)
  herm <- 0

  for(k in 0:end){
    herm <- herm + ((-1)^k * (2*z)^(n-2*k))/(factorial(k)*factorial(n-2*k))
  }

  return(factorial(n)*herm)

}
