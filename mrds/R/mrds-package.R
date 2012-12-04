#' Mark-Recapture Distance Sampling (mrds)
#' 
#' This package implements mark-recapture distance sampling
#'     methods as described in D.L. Borchers, W. Zucchini and Fewster,
#'     R.M. (1988), "Mark-recapture models for line transect surveys",
#'     Biometrics 54: 1207-1220. and Laake, J.L. (1999) "Distance sampling
#'    with independent observers: Reducing bias from heterogeneity by
#'     weakening the conditional independence assumption." in Amstrup,
#'     G.W., Garner, S.C., Laake, J.L., Manly, B.F.J., McDonald, L.L. and
#'     Robertson, D.G. (eds) "Marine mammal survey and assessment
#'     methods", Balkema, Rotterdam: 137-148 and Borchers, D.L., Laake,
#'     J.L., Southwell, C. and Paxton, C.L.G. "Accommodating unmodelled
#'     heterogeneity in double-observer distance sampling surveys". 2006.
#'     Biometrics 62:372-378.)
#' 
#' \tabular{ll}{ Package: \tab mrds \cr Type: \tab Package\cr Version:
#' \tab 2.0.5\cr Date: \tab 2012-3-23\cr License: \tab GPL (>=2)\cr LazyLoad: \tab
#' yes\cr }
#' 
#' @name mrds-package
#' @aliases mrds-package mrds
#' @docType package
#' @author Jeff Laake <jeff.laake@@noaa.gov>, David Borchers <dlb@@mcs.st-and.ac.uk>, 
#'    Len Thomas <len@@mcs.st-and.ac.uk>, David Miller <dlm22@@st-and.ac.uk> and Jon Bishop <jonb@@mcs.st-and.ac.uk>
#' @keywords package
#' 
NULL




#' Golf tee data used in chapter 6 of Advanced Distance Sampling examples
#' 
#' Double platform data collected in a line transect survey of golf tees by 2
#' observers at St. Andrews. Field sex was actually colour of the golf tee: 0 -
#' green; 1 - yellow. Exposure was either low (0) or high(1) depending on
#' height of tee above the ground. size was the number of tees in an observed
#' cluster.
#' 
#' 
#' @name book.tee.data
#' @docType data
#' @format The format is: List of 4 $ book.tee.dataframe:'data.frame': 324 obs.
#'   of 7 variables: ..$ object : num [1:324] 1 1 2 2 3 3 4 4 5 5 ...  ..$
#'   observer: Factor w/ 2 levels "1","2": 1 2 1 2 1 2 1 2 1 2 ...  ..$
#'   detected: num [1:324] 1 0 1 0 1 0 1 0 1 0 ...  ..$ distance: num [1:324]
#'   2.68 2.68 3.33 3.33 0.34 0.34 2.53 2.53 1.46 1.46 ...  ..$ size : num
#'   [1:324] 2 2 2 2 1 1 2 2 2 2 ...  ..$ sex : num [1:324] 1 1 1 1 0 0 1 1 1 1
#'   ...  ..$ exposure: num [1:324] 1 1 0 0 0 0 1 1 0 0 ...  $ book.tee.region
#'   :'data.frame': 2 obs. of 2 variables: ..$ Region.Label: Factor w/ 2 levels
#'   "1","2": 1 2 ..$ Area : num [1:2] 1040 640 $ book.tee.samples
#'   :'data.frame': 11 obs. of 3 variables: ..$ Sample.Label: num [1:11] 1 2 3
#'   4 5 6 7 8 9 10 ...  ..$ Region.Label: Factor w/ 2 levels "1","2": 1 1 1 1
#'   1 1 2 2 2 2 ...  ..$ Effort : num [1:11] 10 30 30 27 21 12 23 23 15 12 ...
#'   $ book.tee.obs :'data.frame': 162 obs. of 3 variables: ..$ object : int
#'   [1:162] 1 2 3 21 22 23 24 59 60 61 ...  ..$ Region.Label: int [1:162] 1 1
#'   1 1 1 1 1 1 1 1 ...  ..$ Sample.Label: int [1:162] 1 1 1 1 1 1 1 1 1 1 ...
#' @keywords datasets
NULL

#' Pronghorn aerial survey data from Wyoming
#' 
#' Detections of pronghorn from fixed-wing aerial surveys in Southeastern
#' Wyoming using four angular bins defined by strut marks. Illustrates data
#' where altitude above ground level (AGL) varies during the survey.
#' 
#' Each record is an observed cluster of pronghorn.  The data provide the
#' stratum for the observation, the direction of travel, the AGL at the time of
#' the observation, the angular bin which contained the center of the pronghorn
#' cluster(group), and the number of pronghorn in the group. The angular bins
#' were defined by a combination of two window and five wing strut marks to
#' define bin cutpoints for perpendicular ground distances of 0-65, 65-90,
#' 90-115, 115-165 and 165-265 meters when the plane is 300' (91.4 meters)
#' above ground level. The inner band is considered a blind region due to
#' obstruction of view beneath the plane; thus th the line is offset 65 meters
#' from underneath the plane.
#' 
#' @name pronghorn
#' @docType data
#' @format A data frame with 660 observations on the following 5 variables.
#'   \describe{ \item{STRATUM}{a numeric vector}
#'   \item{direction}{a factor with levels \code{N} \code{S}
#'   representing the survey direction} \item{AGL}{height above ground
#'   level} \item{Band}{a factor with levels \code{A} \code{B} \code{C}
#'   \code{D} which represent angular bands between breaks at
#'   35.42,44.56,51.52,61.02,70.97 degrees.  These angles were set based on
#'   selected distance bins based on the target AGL.}
#'   \item{cluster}{number of pronghorn in the observed cluster} }
#' @references Laake, J., R. J. Guenzel, J. L. Bengtson, P. Boveng, M. Cameron,
#'   and M. B. Hanson. 2008.  Coping with variation in aerial survey protocol
#'   for line-transect sampling. Wildlife Research 35:289-298.
#' @source Data provided courtesy of Rich Guenzel of Wyoming Game and Fish.
#' @keywords datasets
NULL


#' Wooden stake data from 1977 survey
#' 
#' Multiple surveys by different observers of a single 1km transect containing
#' 150 wooden stakes placed randomly throughout a 40 m strip (20m on either
#' side).
#' 
#' 
#' @name stake77
#' @docType data
#' @format A data frame with 150 observations on the following 10 variables.
#'   \describe{ \item{StakeNo}{unique number for each stake 1-150}
#'   \item{PD}{perpendicular distance at which the stake was placed
#'   from the line} \item{Obs1}{0/1 whether missed/seen by observer 1}
#'   \item{Obs2}{0/1 whether missed/seen by observer 2}
#'   \item{Obs3}{0/1 whether missed/seen by observer 3}
#'   \item{Obs4}{0/1 whether missed/seen by observer 4}
#'   \item{Obs5}{0/1 whether missed/seen by observer 5}
#'   \item{Obs6}{0/1 whether missed/seen by observer 6}
#'   \item{Obs7}{0/1 whether missed/seen by observer 7}
#'   \item{Obs8}{0/1 whether missed/seen by observer 8} }
#' @references Burnham, K. P., D. R. Anderson, and J. L. Laake. 1980.
#'   Estimation of Density from Line Transect Sampling of Biological
#'   Populations. Wildlife Monographs:7-202.
#' @source Laake, J. 1978. Line transect estimators robust to animal movement.
#'   M.S. Thesis. Utah State University, Logan, Utah. 55p.
#' @keywords datasets
#' @examples
#' \donttest{
#' data(stake77)
#' # Functions to extract stake data and put in the mrds format for model fitting.
#' extract.stake=function(stake,obs)
#' {
#'   extract.obs=function(obs)
#'   {
#'      example=subset(stake,eval(parse(text=paste("Obs",obs,"==1",sep=""))),select="PD")
#'      example$distance=example$PD
#'      example$object=1:nrow(example)
#'      example$PD=NULL
#'      return(example)
#'  }
#'  if(obs!="all")
#' 	 return(extract.obs(obs=obs))
#'  else
#'  {
#' 	 example=NULL
#'      for(i in 1:(ncol(stake)-2))
#' 	 {
#' 		 df=extract.obs(obs=i)
#' 		 df$person=i
#' 		 example=rbind(example,df)
#' 	 }		 
#' 	 example$person=factor(example$person)
#' 	 example$object=1:nrow(example)
#' 	 return(example)
#'  }   
#' }
#' extract.stake.pairs=function(stake,obs1,obs2,removal=FALSE)
#' {
#'   obs1=paste("Obs",obs1,sep="")
#'   obs2=paste("Obs",obs2,sep="")
#'   example=subset(stake,eval(parse(text=paste(obs1,"==1 |",obs2,"==1 ",sep=""))),select=c("PD",obs1,obs2))
#'   names(example)=c("distance","obs1","obs2")
#'   detected=c(example$obs1,example$obs2)
#'   example=data.frame(object=rep(1:nrow(example),2),distance=rep(example$distance,2),detected=detected,observer=c(rep(1,nrow(example)),rep(2,nrow(example))))
#'   if(removal)example$detected[example$observer==2]=1
#'   return(example)
#' }
#' # extract data for observer 1 and fit a single observer model
#' stakes=extract.stake(stake77,1)
#' ds.model=ddf(dsmodel = ~mcds(key = "hn", formula = ~1), data = stakes, method = "ds", meta.data = list(width = 20))
#' plot(ds.model,breaks=seq(0,20,2),showpoints=TRUE)
#' ddf.gof(ds.model)
#' # extract data from observers 1 and 3 and fit an io model
#' stkpairs=extract.stake.pairs(stake77,1,3,removal=FALSE)
#' io.model=ddf(dsmodel = ~mcds(key = "hn", formula=~1), mrmodel=~glm(formula=~distance),data = stkpairs, method = "io")
#' summary(io.model)
#' par(mfrow=c(3,2))
#' plot(io.model,breaks=seq(0,20,2),showpoints=TRUE,new=FALSE)
#' ddf.gof(io.model)
#' }
NULL


#' Wooden stake data from 1978 survey
#' 
#' Multiple surveys by different observers of a single 1km transect containing
#' 150 wooden stakes placed based on expected uniform distribution throughout a
#' 40 m strip (20m on either side).
#' 
#' The 1997 survey was based on a single realization of a uniform distribution.
#' Because it was a single transect and there was no randomization of the
#' distances for each survey, we repeated the experiment and used distances
#' that provided a uniform distribution but randomly sorted the positions along
#' the line so there was no pattern obvious to the observer.
#' 
#' @name stake78
#' @docType data
#' @format A data frame with 150 observations on the following 13 variables.
#'   \describe{ \item{StakeNo}{unique number for each stake 1-150}
#'   \item{PD}{perpendicular distance at which the stake was placed
#'   from the line} \item{Obs1}{0/1 whether missed/seen by observer 1}
#'   \item{Obs2}{0/1 whether missed/seen by observer 2}
#'   \item{Obs3}{0/1 whether missed/seen by observer 3}
#'   \item{Obs4}{0/1 whether missed/seen by observer 4}
#'   \item{Obs5}{0/1 whether missed/seen by observer 5}
#'   \item{Obs6}{0/1 whether missed/seen by observer 6}
#'   \item{Obs7}{0/1 whether missed/seen by observer 7}
#'   \item{Obs8}{0/1 whether missed/seen by observer 8}
#'   \item{Obs9}{0/1 whether missed/seen by observer 9}
#'   \item{Obs10}{0/1 whether missed/seen by observer 10}
#'   \item{Obs11}{0/1 whether missed/seen by observer 11} }
#' @references Burnham, K. P., D. R. Anderson, and J. L. Laake. 1980.
#'   Estimation of Density from Line Transect Sampling of Biological
#'   Populations. Wildlife Monographs:7-202.
#' @source Laake, J. 1978. Line transect estimators robust to animal movement.
#'   M.S. Thesis. Utah State University, Logan, Utah. 55p.
#' @keywords datasets
#' @examples
#' \donttest{
#' data(stake78)
#' data(stake77)
#' # compare distribution of distances for all stakes
#' hist(stake77$PD)
#' if(.Platform$GUI=="Rgui")dev.new()
#' hist(stake78$PD)
#' # Functions to extract stake data and put in the mrds format for model fitting.
#' extract.stake=function(stake,obs)
#' {
#'   extract.obs=function(obs)
#'   {
#'      example=subset(stake,eval(parse(text=paste("Obs",obs,"==1",sep=""))),select="PD")
#'      example$distance=example$PD
#'      example$object=1:nrow(example)
#'      example$PD=NULL
#'      return(example)
#'  }
#'  if(obs!="all")
#' 	 return(extract.obs(obs=obs))
#'  else
#'  {
#' 	 example=NULL
#'      for(i in 1:(ncol(stake)-2))
#' 	 {
#' 		 df=extract.obs(obs=i)
#' 		 df$person=i
#' 		 example=rbind(example,df)
#' 	 }		 
#' 	 example$person=factor(example$person)
#' 	 example$object=1:nrow(example)
#' 	 return(example)
#'  }   
#' }
#' extract.stake.pairs=function(stake,obs1,obs2,removal=FALSE)
#' {
#'   obs1=paste("Obs",obs1,sep="")
#'   obs2=paste("Obs",obs2,sep="")
#'   example=subset(stake,eval(parse(text=paste(obs1,"==1 |",obs2,"==1 ",sep=""))),select=c("PD",obs1,obs2))
#'   names(example)=c("distance","obs1","obs2")
#'   detected=c(example$obs1,example$obs2)
#'   example=data.frame(object=rep(1:nrow(example),2),distance=rep(example$distance,2),detected=detected,observer=c(rep(1,nrow(example)),rep(2,nrow(example))))
#'   if(removal)example$detected[example$observer==2]=1
#'   return(example)
#' }
#' # extract data for observer 10 and fit a single observer model
#' stakes=extract.stake(stake78,10)
#' ds.model=ddf(dsmodel = ~mcds(key = "hn", formula = ~1), data = stakes, method = "ds", meta.data = list(width = 20))
#' plot(ds.model,breaks=seq(0,20,2),showpoints=TRUE)
#' ddf.gof(ds.model)
#' # extract data from observers 5 and 7 and fit an io model
#' stkpairs=extract.stake.pairs(stake78,5,7,removal=FALSE)
#' io.model=ddf(dsmodel = ~mcds(key = "hn", formula=~1), mrmodel=~glm(formula=~distance),data = stkpairs, method = "io")
#' summary(io.model)
#' par(mfrow=c(3,2))
#' plot(io.model,breaks=seq(0,20,2),showpoints=TRUE,new=FALSE)
#' ddf.gof(io.model)
#' }
#' 
NULL

#' Single observer point count data example from Distance
#' 
#' Single observer point count data example from Distance
#' 
#' 
#' @name ptdata.distance
#' @docType data
#' @format The format is 144 obs of 6 variables: ..$ distance: numeric distance from center $ 
#'   observer: Factor w/ 2 levels "1","2": 1 2 1 2 1 2 1 2 1 2 ...  ..$
#'   detected: numeric 0/1  $ object: sequential object number $Sample.Label: point label $ Region.Label: single region label
#' @keywords datasets
#' @examples
#' \donttest{
#' data(ptdata.distance)
#' xx=ddf(dsmodel = ~cds(key="hn", formula = ~1), data = ptdata.distance, method = "ds", meta.data = list(point=TRUE))
#' summary(xx)
#' plot(xx,main="Distance point count data")
#' ddf.gof(xx)
#' Regions=data.frame(Region.Label=1,Area=1)
#' Samples=data.frame(Sample.Label=1:30,Region.Label=rep(1,30),Effort=rep(1,30))
#' print(dht(xx,sample.table=Samples,region.table=Regions))
#' }
NULL



#' Simulated single observer point count data
#' 
#' Simulated single observer point count data with detection p(0)=1; hn sigma=30; w=100
#' 
#' 
#' @name ptdata.single
#' @docType data
#' @format The format is 341 obs of 4 variables: ..$ distance: numeric distance from center $ 
#'   observer: Factor w/ 2 levels "1","2": 1 2 1 2 1 2 1 2 1 2 ...  ..$
#'   detected: numeric 0/1  $ object : sequential object number
#' @keywords datasets
#' @examples
#' \donttest{
#' data(ptdata.single)
#' xx=ddf(dsmodel = ~cds(key="hn", formula = ~1), data = ptdata.single, method = "ds", meta.data = list(point=TRUE))
#' summary(xx)
#' plot(xx,main="Simulated point count data")
#' }
NULL

#' Simulated dual observer point count data
#' 
#' Simulated dual observer point count data with detection p(0)=0.8; hn sigma=30; w=100
#' for both observers with dependency y>0, gamma=0.1
#' 
#' @name ptdata.dual
#' @docType data
#' @format The format is 420 obs of 6 variables: ..$ distance: numeric distance from center $ 
#'   observer: Factor w/ 2 levels "1","2": 1 2 1 2 1 2 1 2 1 2 ...  ..$
#'   detected: numeric 0/1  $ person: Factor with 2 levels A,B $ pair: Factor with 2 levels "AB" BA" $
#'   object : sequential object number
#' @keywords datasets
#' @examples
#' \donttest{
#' data(ptdata.dual)
#' xx=ddf(mrmodel=~glm(formula=~distance), dsmodel = ~cds(key="hn", formula = ~1), data = ptdata.dual, method = "io", meta.data = list(point=TRUE))
#' summary(xx)
#' plot(xx,main="Simulated point count data")
#' }
NULL


#' Simulated removal observer point count data
#' 
#' Simulated removal observer point count data with detection p(0)=0.8; hn sigma=30; w=100
#' for both observers with dependency y>0, gamma=0.1
#' 
#' @name ptdata.removal
#' @docType data
#' @format The format is 408 obs of 6 variables: ..$ distance: numeric distance from center $ 
#'   observer: Factor w/ 2 levels "1","2": 1 2 1 2 1 2 1 2 1 2 ...  ..$
#'   detected: numeric 0/1  $ person: Factor with 2 levels A,B $ pair: Factor with 2 levels "AB" BA" $
#'   object : sequential object number
#' @keywords datasets
#' @examples
#' \donttest{
#' data(ptdata.removal)
#' xx=ddf(mrmodel=~glm(formula=~distance), dsmodel = ~cds(key="hn", formula = ~1), data = ptdata.removal, method = "rem", meta.data = list(point=TRUE))
#' summary(xx)
#' plot(xx,main="Simulated point count data")
#' }
NULL
