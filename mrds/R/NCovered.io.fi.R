NCovered.io.fi <-
		function(par=NULL,model,group=TRUE,...)
#
# NCovered.io.fi - computes abundance in covered region for io.fi model object
#
# Arguments: 
#
# par    - parameter values (used when computing derivatives wrt parameter uncertainty)
# model  - ddf model object
# group  - if TRUE computes group abundance and if FALSE individual abundance
#
# Value:
# 
# result - abundance estimate
#
# Functions Used: predict (predict.io.fi), compute.Nht 
{
	if(!is.null(par))
	{
		model$mr$coefficients=par
		fitted=predict.io.fi(model,integrate=TRUE)$fitted
	}
	else
		fitted=model$fitted
	if(!group)
	{
		size=model$data$size[model$data$observer==1&model$data$object %in% as.numeric(names(model$fitted))]
		Nhat=sum(compute.Nht(fitted,FALSE,size))
	}
	else
		Nhat=sum(compute.Nht(fitted,TRUE,size=NULL))
	return(Nhat)
}
