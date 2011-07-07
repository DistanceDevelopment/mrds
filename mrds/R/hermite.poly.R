hermite.poly <-function(z,n)
#
# hermite.poly
#
# Hermite polynomial generating function
#
#
# Arguments:
#
# z	- points to evaluate
# n	- order Hermite polynomial
#
# dlm 23-Aug-05
# Using formula found at http://functions.wolfram.com/Polynomials/HermiteH/02/0001/
# This may not be the fastest way of finding the value.
#
{
	
	end <- floor(n/2)
	herm <- 0
	
	for(k in 0:end){
		herm <- herm + ((-1)^k * (2*z)^(n-2*k))/(factorial(k)*factorial(n-2*k))
	}
	
	return(factorial(n)*herm)
	
}
