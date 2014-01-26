integratelogistic <- function (x, models, beta, lower=0, width, transect){
  # computes integral (numerically) over x from 0 to width of a
  # logistic detection function

  integrate(logisticbyx,lower=lower,upper=width,
            subdivisions=10, rel.tol=0.01,abs.tol=0.01,x=x,
            models=models, beta=beta, transect=transect)$value
}
