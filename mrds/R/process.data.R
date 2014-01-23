#' Process data for fitting distance sampling detection function
#'
#' Sets up dataframe and does some basic error checking. Adds needed fields to
#' dataframe and to \code{meta.data}.
#'
#' The function does a number of error checking tasks, creating fields and
#' adding to \code{meta.data} including:
#'
#' 1) If \code{mr.check=TRUE}, check to make sure the record structure is okay for
#' mrds data. The number of primary records (observer=1) must equal the number
#' of secondary records (observer=2). Also, a field in the dataframe is created
#' \code{timesseen} which counts the number of times an object was detected
#' 0,1,2; if \code{timesseen=0} then the record is tossed from the analysis.
#' Also if there are differences in the data (distance, size, covariates) for
#' observer 1 and 2 a warning is issued that the analysis may fail.  The code
#' assumes these values are the same for both observers.
#'
#' 2) Based on the presence of fields \code{distbegin} and \code{distend}, a
#' determination is made of whether the data analysis should be based on binned
#' distances and a field \code{binned} is created, which is \code{TRUE} if the
#' distance for the observation is binned.  By assigning for each observation
#' this allows an analysis of a mixture of binned and unbinned distances.
#'
#' 4) Data are restricted such that distances are not greater than the right
#' truncation and not less than left truncaton (if values are specified). If
#' they are not specified then the left truncation defaults to 0 and the right
#' truncation defaults to the largest distance measurement.
#'
#' 5) Determine if an integration range (\code{int.begin} and \code{int.end}
#' has been specified for the observations.  If it has, add the structure to
#' \code{meta.data}.  The integration range is typically used for aerial
#' surveys in which the altitude varies such that the strip width (left to
#' width) changes with a change in altitude.
#'
#' 6) Fields defined as factors are cleaned up such that any unused levels are
#' eliminated.
#'
#' 7) If the restrictions placed on the data, eliminated all of the data, the
#' function stops with an error message
#'
#' @param data dataframe object
#' @param truncation either: a single number giving the right truncation for the
#'  distances; a 2-vector giving the left and right truncations
#'  (\code{c(left, right)}); or a list with elements \code{left} and
#'  \code{right} (at least \code{right} must be supplied in list form). Default
#'  is \code{NULL} which uses the largest observed distance (not usually a good
#'  idea with unbinned data).
#' @param meta.data meta.data options; see \code{\link{ddf}} for a description
#' @param control control options; see \code{\link{ddf}} for a description
#' @param mr.check if \code{TRUE} check data for errors in the mark-recapture
#'  part of the model structure; for detection function (\code{method="ds"})
#'  then \code{mr.check=FALSE}.
#' @return \item{xmat}{processed \code{data.frame} with added fields}
#'   \item{meta.data}{meta.data list}
#' @author Jeff Laake
#' @keywords utility
process.data <- function(data,truncation,meta.data=list(),control=list(),
                         mr.check=TRUE){

  ## detection function only checks
  if(!mr.check){
    if(!is.null(data$distance)){
      data <- data[!is.na(data$distance),]
    }else{
      data <- data[!is.na(data$distbegin)&!is.na(data$distend),]
    }
    if(is.null(data$object)){
      stop("\nobject field is missing in the data\n")
    }

  }


  ## MR only checks
  # Check to make sure the record structure is ok. Number of primary
  # records = number of secondary
  if(mr.check){
    if(length(data$detected[data$observer==1]) !=
        length(data$detected[data$observer==2])){
      stop("number of records for primary observer not equal to number for secondary observer")
    }

    # Create field which counts the number of times an object was detected 0,1,2
    timesdetected <- data$detected[data$observer==1] +
                     data$detected[data$observer==2]
    data$timesdetected <- rep(0,dim(data)[1])
    data$timesdetected[data$observer==1] <- timesdetected
    data$timesdetected[data$observer==2] <- timesdetected

    # If any 00 (not detected by either observer), stop and issue error message
    if(any(data$timesdetected==0)){
      stop("following objects were never detected:",
            paste(data$object[data$observer==1&data$timesdetected==0],
                  collapse=","),"\n")
    }
  }

  # Determine if data are binned by presence of distbegin and distend fields
  if(is.null(data$distend) | is.null(data$distbegin)){
    binned <- FALSE
  }else{
    if(all(is.null(data$distend)) | all(is.null(data$distbegin))){
      binned <- FALSE
    }else{
      if(any(is.null(data$distend) & !is.null(data$distbegin)) |
         any(is.null(data$distbegin) & !is.null(data$distend))){
        stop("mismatched distance intervals - one or more endpoints are missing")
      }else{
        binned <- TRUE
      }
    }
  }

  if(meta.data$binned & !binned){
    stop("binned set to TRUE in meta.data but distbegin and distend fields are missing")
  }

  if(!meta.data$binned & binned){
    warning("data contain distbegin and distend fields but binned=FALSE. Analyzing as not binned",immediate.=TRUE)
    binned <- FALSE
  }

  meta.data$binned <- binned

  if(meta.data$binned & is.null(meta.data$breaks)){
    stop("breaks must be set in meta.data for binned data")
  }

  # Fill in distance field for binned observations and create logical variable
  data$binned <- rep(FALSE,dim(data)[1])
  if(binned){
    meta.data$binned <- TRUE
    data$distance[!is.na(data$distend)]<-(data$distbegin[!is.na(data$distend)]+
                                          data$distend[!is.na(data$distend)])/2
    data$binned[!is.na(data$distbegin)] <- TRUE
  }


  ## Restrict data to width interval
  # If no width set, use largest measured distance as width
  if(is.null(truncation)){
    truncation <- list()
    truncation$right <- ifelse(meta.data$binned,
                               max(c(data$distend,data$distance),na.rm=TRUE),
                               max(data$distance))
    xmat <- data
    warning("no truncation distance specified; using largest observed distance",immediate.=TRUE)
  }else{

    ## ensure that the truncation is in the right format for later or
    ## throw an error if the input was garbage
    if(is.list(truncation)){
      if(!all(names(truncation)%in%c("left","right"))){
        stop("truncation must be a number, a 2-vector or a named list")
      }
      # otherwise truncation has been supplied correctly as a list with
      # elements "left" and "right"
    }else if(length(truncation)==2){
      # convert 2 vector into a list
      truncation <- list(left=min(truncation),right=max(truncation))
    }else if(length(truncation)==1){
      truncation <- list(left=0,right=truncation)
    }else{
      stop("truncation must be a number, a 2-vector or a named list")
    }

    # This piece of code makes sure that the set width is as large as the
    # largest bin end point for binned data.
    if(meta.data$binned){
      if(any(data$binned & data$distend > truncation$right)){
        stop("width must exceed largest interval end point")
      }else{
        xmat <- data[data$binned |
                     (!data$binned&data$distance<=truncation$right),]
      }
    }else{
      xmat <- data[data$distance <= truncation$right,]
    }
  }


  # Determine if integration range has been specified
  if(is.null(xmat$int.begin)|is.null(xmat$int.end)){
    if(any(is.na(meta.data$int.range))){
      meta.data$int.range <- c(truncation$left,truncation$right)
    }
  }else{
      meta.data$int.range <- rbind(c(truncation$left,truncation$right),
                                   cbind(xmat$int.begin,xmat$int.end))
  }

  # If left >0 perform left truncation by restricting values
  if(truncation$left >0){
    if(binned){
      if(any(data$binned & (data$distbegin < truncation$left))){
        stop("left truncation must be smaller than the smallest interval begin point")
      }else{
        xmat <- data[data$binned|(!data$binned &
                                  (data$distance >= truncation$left)),]
      }
    }else{
      xmat <- xmat[xmat$distance >= truncation$left,]
    }
  }

  # Clean up factor levels
  b <- dim(xmat)[2]
  for(i in 1:b){
    if(is.factor(xmat[,i])){
      xmat[,i] <- factor(xmat[,i])
    }
  }

  ## further detection function only checks
##  if(!mr.check){
##    # use all unique detections (observer=1) if observer is present
##    if(!is.null(xmat$observer)){
##      if(control$limit){
##        if(length(levels(factor(xmat$observer)))>1){
##          xmat <- xmat[xmat$observer==levels(factor(xmat$observer))[1],]
##          xmat$detected <- rep(1,dim(xmat)[1])
##        }
##      }
##    }
##
##    # If the frame includes column "detected", use only those with detected=1
##    if(!is.null(xmat$detected)){
##      if(control$limit) xmat <- xmat[xmat$detected==1,]
##    }else{
##      xmat$detected <- rep(1,dim(xmat)[1])
##    }
##
##    #  Make sure object #'s are unique
##    if(length(unique(xmat$object))!=length(xmat$object)){
##      stop("\nSome values of object field are duplicates. They must be unique.\n")
##    }
##  }

  # If the exclusion eliminated all of the data, stop with error message
  if(dim(xmat)[1]==0){
    stop("no data to analyze")
  }

  return(list(xmat=xmat,meta.data=meta.data,truncation=truncation))
}
