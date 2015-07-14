#' Calculate values of the average detection function when covariates are present
plot_avg_df <- function(xgrid, width, xmat, ddfobj, selected, int.range, pdot, normalize, range.varies){

  finebr <- (width/100)*(0:100)
  xgrid2 <- NULL
  linevalues <- NULL

  newdat <- xmat
  for(i in 1:length(xgrid)){
    xgrid2 <- c(xgrid2,xgrid[i])
    newdat$distance <- rep(xgrid[i],nrow(newdat))

    detfct.values <- detfct(newdat$distance,ddfobj,select=selected,
                           width=width)

    if(!normalize&range.varies){
      detfct.values[xgrid[i]<int.range[,1]|xgrid[i]>int.range[,2]] <- 0
    }
    linevalues <- c(linevalues,sum(detfct.values/pdot)/sum(1/pdot))
  }
  return(linevalues)
}
