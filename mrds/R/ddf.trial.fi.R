#' Mark-Recapture Analysis of Trial Configuration - FI
#'
#' Mark-Recapture Analysis of Trial Observer Configuration with Full
#' Independence
#'
#' The mark-recapture data derived from a trial observer distance sampling
#' survey can be used to derive a conditional detection function (p_1(y)) for
#' observer 1 based on trials (observations) from observer 2. It is a
#' conditional detection function because detection probability for observer 1
#' is based on seeing or not seeing observations made by observer 2. Thus,
#' p_1(y) is estimated by p_1|2(y).  If detections by the observers are
#' independent (full independence) then p_1(y)=p_1|2(y) for each distance y.
#' In fitting the detection functions the likelihood given by eq 6.12 or 6.17
#' in Laake and Borchers (2004) is used. That analysis does not require the
#' usual distance sampling assumption that perpendicular distances are
#' uniformly distributed based on line placement that is random relative to
#' animal distribution.  However, that assumption is used in computing
#' predicted detection probability which is averaged based on a uniform
#' distribution (see eq 6.13 of Laake and Borchers 2004).
#'
#' For a complete description of each of the calling arguments, see
#' \code{\link{ddf}}.  The argument \code{model} in this function is the same
#' as \code{mrmodel} in \code{ddf}.  The argument \code{dataname} is the name
#' of the dataframe specified by the argument \code{data} in \code{ddf}. The
#' arguments \code{control},\code{meta.data},and \code{method} are defined the
#' same as in \code{ddf}.
#'
#' @S3method ddf trial.fi
#' @param model mark-recapture model specification
#' @param data analysis dataframe
#' @param truncation either: a single number giving the right truncation for the
#'  distances; a 2-vector giving the left and right truncations
#'  (\code{c(left, right)}); or a list with elements \code{left} and
#'  \code{right} (at least \code{right} must be supplied in list form). Default
#'  is \code{NULL} which uses the largest observed distance (not usually a good
#'  idea with unbinned data).
#' @param transect is the survey \code{"point"} or \code{"line"} transects?
#'  (Default \code{"line"}.)
#' @param meta.data list containing settings controlling data structure
#' @param control list containing settings controlling model fitting
#' @param call original function call used to call \code{ddf}
#' @param method analysis method; only needed if this function called from
#'   \code{ddf.trial}
#' @return result: a trial.fi model object
#' @author Jeff Laake
#' @seealso
#'   \code{\link{ddf.trial}},\code{\link{summary.trial.fi}},\code{\link{coef.trial.fi}},\code{\link{plot.trial.fi}},
#'   \code{\link{gof.trial.fi}}
#' @references Laake, J.L. and D.L. Borchers. 2004. Methods for incomplete
#'   detection at distance zero. In: Advanced Distance Sampling, eds. S.T.
#'   Buckland, D.R.Anderson, K.P. Burnham, J.L. Laake, D.L. Borchers, and L.
#'   Thomas. Oxford University Press.
#' @keywords Statistical Models
ddf.trial.fi <- function(model, data, truncation=NULL, transect="line",
                         meta.data=list(), control=list(), call="", method){
  #  NOTE: gams are only partially implemented

  # The following are dummy glm and gam functions that are defined here to
  # provide the list of arguments for use in the real glm/gam functions. These
  # dummy functions are removed after they are used so the real ones can be
  # used in the model fitting.
  glm <- function(formula,link="logit"){
    if(class(formula)!="formula"){
      if(class(try(as.formula(formula)))!="formula"){
        stop("Invalid formula")
      }
    }else{
      formula <- paste(as.character(formula),collapse="")
    }
    if(class(link)=="function"){
      link <- substitute(link)
    }
    link <- match.arg(link,c("logit"))
    return(list(fct="glm",formula=formula,link=substitute(link)))
  }

  gam <- function(formula,link="logit"){
    if(class(formula)!="formula"){
      if(class(try(as.formula(formula)))!="formula"){
        stop("Invalid formula")
      }
    }else{
      formula <- paste(as.character(formula),collapse="")
    }
    if(class(link)=="function"){
      link <- substitute(link)
    }
    link <- match.arg(link,c("logit"))
    return(list(fct="gam",formula=formula,link=substitute(link)))
  }

  # Test to make sure that observer not used in mrmodel
  if(length(grep("observer",model))!=0){
     stop("observer cannot be included in models for trial configurations\n")
  }

  # Save current user options and then set design contrasts to treatment style
  save.options <- options()
  options(contrasts=c("contr.treatment","contr.poly"))

  # Set up meta data values
  meta.data <- assign.default.values(meta.data, binned=FALSE, int.range=NA,
                                     mono=FALSE, mono.strict=TRUE)

  # Set up control values
  control <- assign.default.values(control, showit=0, doeachint=FALSE,
                                   estimate=TRUE, refit=TRUE, nrefits=25,
                                   initial=NA, lowerbounds=NA, upperbounds=NA,
                                   mono.points=20)

  # Assign model values; this uses temporarily defined functions glm and gam
  modpaste <- paste(model)
  modelvalues <- try(eval(parse(text=modpaste[2:length(modpaste)])))
  if(class(modelvalues)=="try-error"){
    stop("Invalid model specification: ",model)
  }
  rm(glm,gam)

  # Process data if needed
  if(is.data.frame(data)){
    data.list <- process.data(data,truncation,meta.data)
    meta.data <- data.list$meta.data
    xmat <- data.list$xmat
    truncation <- data.list$truncation
  }else{
    xmat <- data
  }

  # Setup default breaks
  if(meta.data$binned){
    meta.data$breaks <- c(max(0,
                           min(as.numeric(levels(as.factor(xmat$distbegin))))),
                          as.numeric(levels(as.factor(xmat$distend))))
  }

  # Create result list with some arguments
  result <- list(call=call,data=data,model=model,truncation=truncation,
                 meta.data=meta.data,control=control,method="trial.fi",
                 truncation=truncation,transect=transect)
  class(result) <- c("trial.fi","ddf")

  # Fit the conditional detection function using glm; need to add gam capability
  GAM <- FALSE
  if(modelvalues$fct=="gam"){
    GAM <- TRUE
  }

  # Fit the conditional detection function for primary using glm

  # Create formula and model frame
  xmat$offsetvalue <- rep(0,dim(xmat)[1])
  GAM <- FALSE
  if(modelvalues$fct=="gam"){
    GAM <- TRUE
  }
  model.formula <- paste("detected",modelvalues$formula)
  data <- xmat[xmat$observer==1,][xmat$detected[xmat$observer==2]==1,]
  data$distance <- xmat$distance[xmat$detected==1&xmat$observer==2]
  data <- create.model.frame(data,as.formula(model.formula),meta.data)

  # Fit model
  result$mr <- glm(as.formula(model.formula),family=binomial,data=data)
  result$par <- coef(result$mr)
  npar <- length(result$par)
  result$lnl <- -result$mr$deviance/2
  if(GAM){
    result$hessian <- solve(result$mr$Vp)
  }else{
    result$hessian <- solve(summary(result$mr)$cov.unscaled)
  }

  # Compute fitted values
  result$fitted <- predict(result,
                           newdata=xmat[xmat$observer==1&xmat$detected==1,],
                           integrate=TRUE)$fitted

  # If this is not TI mode, compute the lnlU2 portion of the likelihood
  # using the conditional detection functions.
  if(method=="trial.fi"){
    distances <- xmat$distance[xmat$observer==1&xmat$detected==1]
    if(!meta.data$binned){
      if(transect=="point"){
        result$lnl <- result$lnl +
                      sum(log(predict(result,
                               newdat=xmat[xmat$observer==1&xmat$detected==1,],
                               integrate=FALSE)$fitted*
                          2*distances/truncation$right^2))-
                      sum(log(result$fitted))
      }else{
        result$lnl <- result$lnl +
                      sum(log(predict(result,
                               newdat=xmat[xmat$observer==1&xmat$detected==1,],
                               integrate=FALSE)$fitted/truncation$right)) -
                      sum(log(result$fitted))
      }
    }else{
      for(i in 1:(nrow(xmat)/2)){
        int.val <- predict(result,newdata=xmat[(2*(i-1)+1):(2*i),],
                           int.range=as.vector(as.matrix(xmat[(2*(i-1)+1),
                                                    c("distbegin","distend")])),
                           integrate=TRUE)$fitted
        result$lnl <- result$lnl + log(int.val)
      }
      result$lnl <- result$lnl - sum(log(result$fitted))
    }
  }

  # Compute AIC and Nhat in covered region
  result$Nhat <- NCovered(result)
  result$criterion <- -2*result$lnl + 2*npar
  result$par <- coef(result$mr)

  # Restore user options
   options(save.options)

  # Return result
   return(result)
}
