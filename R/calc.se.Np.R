#' Find se of average p and N
#'
#' @param model a \code{ddf} model object
#' @param avgp average p function
#' @param n sample size
#' @param average.p the average probability of detection for the model
#'
#' @author David L. Miller
calc.se.Np <- function(model, avgp, n, average.p){

  # if we have any NAs in the hessian then everything is
  # messed up, so we should just return NAs
  if(any(is.na(model$hessian))){
    se.obj <- list()
    se.obj$Nhat.se <- NA
    se.obj$average.p.se <- NA
    model$hessian <- as.matrix(model$hessian)
    se.obj$Nhatvar.list <-list(variance = matrix(NA, nrow(model$hessian),
                                                 ncol(model$hessian)),
                               partial  = matrix(NA, length(model$par),
                                                 length(model$par)))
    se.obj$vcov <- matrix(NA, nrow(model$hessian), ncol(model$hessian))
  }else{
    # calculate the variance-covariate matrix
    vcov <- solvecov(model$hessian)$inv

    # calculate Nhat uncertainty
    Nhatvar.list <- DeltaMethod(model$par, NCovered, vcov, 0.001,
                                model=model, group=TRUE)
    Nhatvar <- Nhatvar.list$variance + sum((1-model$fitted)/model$fitted^2)
    cvN <- sqrt(Nhatvar)/model$Nhat

    # calculate the average p uncertainty
    var.pbar.list <- prob.se(model, avgp, vcov)
    covar <- t(Nhatvar.list$partial) %*% vcov %*% var.pbar.list$partial+
                 var.pbar.list$covar
    var.pbar <- average.p^2*(cvN^2 + var.pbar.list$var/n^2-
                                  2*covar/(n*model$Nhat))

    # what should we return?
    se.obj <- list()
    se.obj$Nhat.se <- sqrt(Nhatvar)
    se.obj$average.p.se <- sqrt(var.pbar)
    se.obj$Nhatvar.list <- Nhatvar.list
    se.obj$vcov <- vcov
  }

  return(se.obj)
}
