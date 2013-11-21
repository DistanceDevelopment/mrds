# see predict.ds for documentation
#' @S3method predict io.fi
predict.io.fi <- function(object,newdata=NULL,compute=FALSE, int.range=NULL,
                          integrate=FALSE,...){
  # Functions Used: pdot.dsr.integrate.logistic, is.linear.logistic,
  #                 predict.glm (could also use predict.gam eventually)
  model <- object
  width <- model$meta.data$width
  point <- model$meta.data$point
  if(is.null(newdata)){
    newdata <- model$mr$data
  }

  newdata$offsetvalue <- 0

  GAM <- FALSE
  if("gam" %in% class(model$mr)){
    GAM <- TRUE
  }


  if(!integrate){
    fitted <- predict(model$mr,newdata,type="response")
    p1 <- fitted[model$mr$data$observer==1]
    p2 <- fitted[model$mr$data$observer==2]
    fitted <- p1+p2-p1*p2
    if(is.null(newdata)){
      names(fitted) <- model$mr$data$object[model$mr$data$observer==1]
    }else{
      names(fitted) <- newdata$object[newdata$observer==1]
    }
    return(list(fitted=fitted,p1=p1,p2=p2))
  }else{
    left <- model$meta.data$left
    formula <- paste("~",as.character(model$mr$formula)[3],collapse="")
    if("gam" %in% class(model$mr)){
      integral.numeric <- TRUE
    }else{
      integral.numeric <- is.linear.logistic(newdata,formula,
                                             length(coef(model$mr)),width)
    }
    models <- list(g0model=formula,scalemodel=NULL,fullscalemodel=NULL)

    # now int.range is a vector with lower and upper bounds
    if(is.null(int.range)){
      pdot.list <- pdot.dsr.integrate.logistic(width,width, model$mr$coef,
                     newdata,integral.numeric, FALSE, models,GAM, point=point)
    }else{
      pdot.list <- pdot.dsr.integrate.logistic(int.range,width, model$mr$coef,
                       newdata,integral.numeric, FALSE, models,GAM, point=point)
    }

    if(left !=0){
      pdot.list$pdot <- pdot.list$pdot -
                    pdot.dsr.integrate.logistic(left, width, model$mr$coef,
                                   newdata,integral.numeric, FALSE, models,GAM,
                                   point=point)$pdot
    }

    fitted <- pdot.list$pdot
    if(is.null(newdata)){
       names(fitted) <- model$mr$data$object[model$mr$data$observer==1]
    }else{
       names(fitted) <- newdata$object[newdata$observer==1]
    }

    return(list(fitted=fitted))
  }
}
