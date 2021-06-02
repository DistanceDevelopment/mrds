# Hermite polynomial generating function
#
# z points to evaluate
# n order Hermite polynomial
#
# return Hermite polynomials evaluated at z
# documented in ?distpdf
hermite.poly <- function(z, n){

  # note we've done all the input checking prior to calling
  # this function

  # from Distance CDS HRMADJ
  # note that these are "probabalists" Hermite polynomials
  # according to wikipedia https://en.wikipedia.org/wiki/Hermite_polynomials

  # code adapted from https://github.com/DistanceDevelopment/distance-for-windows/blob/bcce3887fa0eea83824732ac540faeea530a9bff/Distance70/Analysis%20Engines/CDS/Engine/mcds/Funcs.for
  switch(as.character(n),
         "0" = 1,
         "1" = z,
         "2" = z*z-1,
         "3" = z^3-3*z,
         "4" = z^4-6*z*z+3,
         "5" = z^5-10*z^3+15*z,
         "6" = z^6-15*z^4+45*z*z-15,
         "7" = z^7-21*z^5+105*z^3-105*z,
         "8" = z^8-28*z^6+210*z^4-420*z*z+105,
         "9" = z^9-36*z^7+378*z^5-1260*z^3+945*z,
         "10" = z^10-45*z^8+630*z^6-3150*z^4+4725*z*z-945)
}
