#' Extrapolate Horvitz-Thompson abundance estimates to entire surveyed region
#'
#' @param Nhat.by.sample dataframe of abundance by sample
#' @param samples samples table
#' @param width truncation width
#' @param left left truncation if any
#' @param point if TRUE point count otherwise line transect
#' @param areas.supplied if \code{TRUE}, covered area is extracted from the
#' \code{CoveredArea} column of \code{Nhat.by.sample}
#' @return Revised Nhat.by.sample dataframe containing estimates extrapolated
#'   to survey region
#' @note Internal function called by \code{\link{dht}} and related functions.
#' @author Jeff Laake and David L Miller
#' @keywords utility
survey.region.dht <- function(Nhat.by.sample, samples, width, left, point,
                              areas.supplied){
  #  Compute effort in each region
  Effort.by.region <- by(samples$Effort, samples$Region.Label, sum)

  # calculate the areas, unless they were given already
  if(areas.supplied){
    CoveredArea <- as.vector(by(samples$CoveredArea, samples$Region.Label, sum))
    Nhat.by.sample$CoveredArea <- NULL
  }else{
    if(point){
      CoveredArea <- pi*as.vector(Effort.by.region)*width^2 -
                      pi*as.vector(Effort.by.region)*left^2
    }else{
      CoveredArea <- 2*as.vector(Effort.by.region)*(width-left)
    }
  }
  # Scale up abundance in covered region to the survey region
  # unless no areas given
  Nhat.by.region <- merge(Nhat.by.sample,
                          data.frame(Region.Label=names(Effort.by.region),
                                     CoveredArea=CoveredArea,
                                     Effort=as.vector(Effort.by.region)),
                          by.x="Region.Label",
                          by.y="Region.Label",
                          all.x=TRUE)
  Nhat.by.region$Nhat <- Nhat.by.region$Nhat*Nhat.by.region$Area/
                          Nhat.by.region$CoveredArea
  return(Nhat.by.region)
}
