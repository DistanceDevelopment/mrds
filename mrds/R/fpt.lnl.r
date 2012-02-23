fpt.lnl <-
function(fpar, ddfobj,TCI,misc.options)
#
# fpt.lnl - computes negative log-likelihood values of point transect grouped/ungrouped distances
#
# Arguments: see flnl for description of arguments
# 
# value: vector of negative log-likelihood values for observations at current set of parameters (fpar)
#
# Functions Used: assign.par, tablecgf, predict(predict.smooth.spline), detfct, integratedetfct, integratedetfct.logistic
#
{
#
   int.range=misc.options$int.range
   width=misc.options$width
   ddfobj=assign.par(ddfobj,fpar)
#
#   Setup integration ranges
#
  if(is.vector(int.range)){
    left=int.range[1]
    right=int.range[2]
  }else if(is.matrix(int.range)){
    left=int.range[2:dim(int.range)[1],1]
    right=int.range[2:dim(int.range)[1],2]
  }       
#
#   Compute log-likelihood for binned data
#
  z=ddfobj$scale$dm
  x=ddfobj$xmat
  lnl=rep(0,dim(x)[1])
  if(any(x$binned))
  {
    if(ddfobj$type=="hr") ddfobj$cgftab <- tablecgf(ddfobj,width=width,standardize=misc.options$standardize,point=TRUE)
    key.scale<-ddfobj$scale$parameters
    intall <- predict(ddfobj$cgftab, as.vector(right/scalevalue(key.scale, z[x$binned,])))$y * 
                        scalevalue(key.scale,z[x$binned,])^2 -
              predict(ddfobj$cgftab, as.vector(left/scalevalue(key.scale, z[x$binned,])))$y *
                          scalevalue(key.scale,z[x$binned,])^2
      
    intbegin <- predict(ddfobj$cgftab, as.vector(x$distbegin[x$binned]/scalevalue(key.scale, z[x$binned,])))$y * scalevalue(key.scale,z[x$binned,])^2    
    intend <- predict(ddfobj$cgftab, as.vector(x$distend[x$binned]/scalevalue(key.scale, z[x$binned,])))$y * scalevalue(key.scale,z[x$binned,])^2
#
#    24-Aug-05; jll; To avoid numerical problems due to approximation of integral the following
#    line of code was added; added this on 19 Jan 06 to mrds-dev
#
    intend[intend<= intbegin]=intbegin[intend<= intbegin]+1e-8
    if(is.vector(left))
      intall[is.infinite(intall)]<- right - left
    else
      intall[is.infinite(intall)]<- right[is.infinite(intall)] - left[is.infinite(intall)]
      
    lnl[x$binned] <-  - log((intend-intbegin)/intall)
  }
#
#   Compute log-likelihood for unbinned data
#
# dlm 13-Oct-11 setting standardize=FALSE here
#               if not then we get divide by zero errors...
  if(!all(x$binned)){
    p1 <- fr(x$distance[!x$binned],ddfobj=ddfobj,select=!x$binned,width=width,standardize=FALSE)
    p1[p1<1.0e-15]<-1.0e-15
    p1[is.nan(p1)]<-1.0e-15
    int1 <- integratedetfct(ddfobj,select=!x$binned,width=width,int.range=int.range,
			                 doeachint=misc.options$doeachint,standardize=misc.options$standardize,point=TRUE)
    if(is.vector(left)){
      int1[is.infinite(int1)]<-  right - left
      int1[is.nan(int1)]<- right - left
    }else{ 
      int1[is.infinite(int1)]<-  right[is.infinite(int1)] - left[is.infinite(int1)]
      int1[is.nan(int1)]<- right[is.nan(int1)] - left[is.nan(int1)]
    }
  
    lnl[!x$binned] <- - log(p1/int1)

   if(any(is.nan(lnl[!x$binned])))
     lnl[!x$binned] <- NA

  }
  return(lnl)
}
