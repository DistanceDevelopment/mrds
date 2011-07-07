integratelogisticdup <-
function (x1, x2, models, beta,width, point) 
{ 
#
#  integratelogisticdup - numerical integral of product of logistic detection functions
#
integrate(logisticdupbyx,lower=0,upper=width, subdivisions=10, rel.tol=0.01,abs.tol=0.01,
                                x1=x1,x2=x2, models=models, beta=beta, point=point)$value 
}
