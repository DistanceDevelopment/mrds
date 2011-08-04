#' MCDS function definition
#' 
#' Creates model formula list for multiple covariate distance sampling using
#' values supplied in call to \code{\link{ddf}}
#' 
#' 
#' @param formula formula for scale function
#' @param key string identifying key function (currently either "hn"
#'   (half-normal),"hr" (hazard-rate), "unif" (uniform) or "gamma" (gamma
#'   distribution)
#' @param adj.series string identifying adjustment functions cos (Cosine), herm
#'   (Hermite polynomials), poly (simple polynomials) or NULL
#' @param adj.order vector of order of adjustment terms to include
#' @param adj.scale whether to scale the adjustment terms by "width" or "scale"
#' @param shape.formula formula for shape function
#' @return A formula list used to define the detection function model
#'   \item{fct}{string "mcds"} \item{key}{key function string}
#'   \item{adj.series}{adjustment function string} \item{adj.order}{adjustment
#'   function orders} \item{adj.scale}{adjustment function scale type}
#'   \item{formula}{formula for scale function} \item{shape.formula}{formula
#'   for shape function}
#' @author Jeff Laake; Dave Miller
#' @keywords utility
mcds <-
		function(formula,key=NULL,adj.series=NULL,adj.order=c(NULL),adj.scale="width",shape.formula=~1)
#
#  mcds - creates model formula list for multiple covariate distance sampling
#
#  Arguments:
#
#  formula	- formula for scale function
#  key		- either hn (half-normal) or hr (hazard rate)
#  adj.series	- cos (Cosine), herm (Hermite polynomials), 
#		  poly (simple polynomials) or NULL
#  adj.order	- order of terms to include
#  adj.scale	- whether to use "width" or "scale" to scale the adjustment terms
#  shape.formula- formula for shape function	
#
#  Value:
#
#   model list
#
# dlm 11/07/05	Added in handling for adjustment terms
#
{
	if(class(formula)!="formula")
		if(class(try(as.formula(formula)))=="formula")
			formula=as.formula(formula)
		else
			stop("Invalid formula") 
	if(!is.null(shape.formula))
		if(class(shape.formula)!="formula")
			if(class(try(as.formula(shape.formula)))=="formula")
				shape.formula=as.formula(shape.formula)
			else
				stop("Invalid shape.formula") 
	
	key <- match.arg(key,c("hn","hr","unif","gamma"))
	if(key%in%c("hn","unif"))shape.formula=NULL
# What to do if we have adjustment terms
	
	if(!is.null(adj.series)){
		adj.series <- match.arg(adj.series,c("cos","herm","poly"))
		adj.scale <- match.arg(adj.scale,c("width","scale"))
		if(adj.check.order(adj.series,adj.order))
			stop("Problem with adjustment terms, see above errors")
	}
	
	return(list(fct="mcds",formula=formula,shape.formula=shape.formula,key=substitute(key),adj.series=adj.series,adj.order=adj.order,adj.scale=adj.scale))
	
}
