#' Covered region estimate of abundance from Horvitz-Thompson-like estimator
#' 
#' Computes H-T abundance within covered region by sample.
#' 
#' 
#' @param obs observations table
#' @param samples samples table
#' @param group if TRUE compute abundance of group otherwise abundance of
#'   individuals
#' @return Nhat.by.sample - dataframe of abundance by sample
#' @note Internal function called by \code{\link{dht}} and related functions
#' @author Jeff Laake
#' @keywords utility
covered.region.dht <-
		function(obs, samples, group)
#
#  covered.region.dht - computes H-T abundance within covered region by sample; used by dht
#  
#
#  Arguments:
#  
#  obs     - table of observations
#  samples - table of samples
#  group   - if TRUE compute abundance of groups
#
#  Value:
# 
#  Nhat.by.sample - dataframe of abundance by sample
{
#
#  Compute abundance in covered region depending on value of
#  group = TRUE (do group abundance); F(do individual abundance)
#  if there are observations of this species
	if(nrow(obs) > 0){
		Nhats = compute.Nht(obs$pdot, group, obs$size)
#
#  Sum abundances by sample within region  
#
		Nhats = by(Nhats, obs$Label, sum)
#
#  Sum observations by sample within region
#  1/25/06 jll - added code to add count of observations
#
		if(group){
			sum.obs = by(obs$object, obs$Label, length)
		}
		else{
			sum.obs = by(obs$size, obs$Label, sum)
		}
#
#  Merge with samples
#
		num.obs = data.frame(Label = names(sum.obs), n = as.vector(sum.obs))
		Nhats = data.frame(Label = names(Nhats), Nhat = as.vector(Nhats))
		Nhat.by.sample = merge(samples, Nhats,by.x = "Label", by.y = "Label", all.x = TRUE)
		Nhat.by.sample = merge(Nhat.by.sample, num.obs, by.x = "Label", by.y = "Label", all.x = TRUE)
		Nhat.by.sample$Nhat[is.na(Nhat.by.sample$Nhat)] = 0
		Nhat.by.sample$n[is.na(Nhat.by.sample$n)] = 0
#
#  Create Nhat.by.sample in the case where there were no sightings
#  28/05/12 lhm - added code to prevent error in the case of no observations
#
	}else{
		Nhat.by.sample <- cbind(samples, Nhat = rep(0,nrow(samples)), n = rep(0,nrow(samples)))
	}
	return(Nhat.by.sample)
}
