NCovered.ds <-
		function(par=NULL,model,group=TRUE,...)
#
# NCovered.ds - computes abundance in covered region for ds model object
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
# Functions Used: predict (predict.ds), compute.Nht 
{
	if(!is.null(par) | is.null(model$fitted))
	{
		model$par=par 
		fitted=predict(model,esw=FALSE,compute=TRUE)$fitted
	}
	else
		fitted=model$fitted
	if(!group)
	{
		size=model$data$size[model$data$object %in% as.numeric(names(model$fitted))]
		Nhat=sum(compute.Nht(fitted,FALSE,size))
	}
	else
		Nhat=sum(compute.Nht(fitted,TRUE,size=NULL))
	return(Nhat)
}
