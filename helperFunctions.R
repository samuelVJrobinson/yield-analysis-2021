# Helper functions

seqGroup <- function(x,singles=TRUE){ #Turns sequential IDs into numbered groups.
  #Singles = TRUE, sequential values must be exactly 1 greater (gaps matter)
  #Singles = FALSE, sequential values must be greater (gaps don't matter)
  if(singles){
    xl <- (x-lag(x))!=1
  } else {
    xl <- (x-lag(x))<0
  }
  cumsum(ifelse(is.na(xl),FALSE,xl))+1
}

makePolys <- function(dat,width='w',dist='d',angle='a',backwards=FALSE){
  # Make polygons from width, dist, and angle measurements, centered on location from dat
  
  gType <- dat %>% st_geometry_type(FALSE)
  
  if(gType!='POINT') warning(paste('Input data type is',gType,'not POINT',sep=' '))
  
  rectFun <- function(x,y,w,d,a,b){
    #Function to create corners of rotated rectangle from:
    #   starting location (x,y),
    #   width (w), distance (d), and angle (a)
    #   rotate rectangles 180 degrees? (b)
    rotate <- ifelse(b,90,270)
    a <- (rotate-a)*pi/180 #90-angle in radians. 
    v <- c(x,y) #Starting point
    v1 <- c(d*cos(a),d*sin(a)) #Vectors to add together to get corners
    v2 <- c(-(w/2)*sin(a),(w/2)*cos(a))
    return(rbind(v+v2,v+v1+v2,v+v1-v2,v-v2,v+v2)) #Corners of rotated rectangle
  }
  datCRS <- st_crs(dat) #Coordinate system from dat
  #Apply rectFun to all rows in dat. Probably could be done through map or some other purrr related thing.
  xcoord <- st_coordinates(dat)[,1]
  ycoord  <- st_coordinates(dat)[,2]
  polys <- lapply(1:nrow(dat),function(i){
    r <- rectFun(x = xcoord[i], y = ycoord[i],  w = dat[[width]][i],  d = dat[[dist]][i],
                 a = dat[[angle]][i], b = backwards)
    st_polygon(list(r))
  })
  #Combine dat with new polygon geometry, and add original CRS
  dat2 <- st_sf(st_drop_geometry(dat),st_sfc(polys)) %>% st_set_crs(datCRS)
  return(dat2)
}

mergePoly <- function(dat,fList=NULL){
  #Function to merge together completely overlapping yield polygons (cover or perfect overlap)
  # takes a list of functions (fList)
  # columns not mentioned are dropped
  
  vList <- as.list(names(fList))
  # doesContain <- unique(st_contains(dat)) #Deals with identical overlap
  doesContain <- st_covers(dat) 
  #Find polygons that are overlapped by others
  remove <- unique(unlist(lapply(doesContain,function(x) if(length(x)>1) x[2:length(x)] else logical(0)))) 
  doesContain <- doesContain[-remove] #Remove polygons from list
  
  if(length(doesContain)==0) return(dat) #If polygons completely overlap, return original data
  
  # cMat <- table(rep(1:length(doesContain),sapply(doesContain,length)),unlist(doesContain))
  # diag(cMat) <- 0
  # doesContain <- doesContain[colSums(cMat)[1:nrow(cMat)]==0] #Strips out inverse overlap cases
  # any(sapply(doesContain,length)>1) #Do any polygons overlap completely?
  
  newPolys <- lapply(doesContain,function(i){
    if(length(i)>1){ 
      p1 <- dat[i,] #Dataframe with overlapping polygons
      newData <- map2(vList,fList, ~ st_drop_geometry(p1) %>% summarize(across(all_of(.x),.y))) %>%   
        reduce(data.frame)
      rownames(newData) <- paste(i,collapse='.')
      newGeom <- st_union(p1) %>% st_geometry()
      return(list(data=newData,geom=newGeom))
    } else if(i>0) {
      newData <- dat[i,unlist(vList)] %>% st_drop_geometry() 
      newGeom <- dat[i,] %>% st_geometry()
      return(list(data=newData,geom=newGeom))
    }
  }) 
  
  dat2 <- do.call('rbind',lapply(newPolys,function(x) x$data)) %>%
    st_set_geometry(st_sfc(sapply(newPolys,function(x) x$geom),crs=st_crs(dat)))
  
  return(dat2)
}

#Yield conversions 
# Requires starting and ending units
convertYield <- function(x,inUnits=NULL,outUnits=NULL){
  
  if(identical(inUnits,outUnits)) return(x)
  
  #Convert to g/m2
  x <- switch(inUnits,
              tpha = x/100, #Tonnes/ha
              gpm2 = x #g per m2
              )
  
  #Convert to output units
  xOut <- switch(outUnits,
                 tpha = x*100,
                 gpm2 = x)
  
  return(xOut)
} 

#Function to get smooth predictions from a model list
#smoothLabel = label of smooth to choose
#modList = model list to use
#predRange = range to predict over
#lengthOut = number of predictions over
#postSamp = sample posterior?
#noIntercept = set intercept to 0?
#returnDF = return DF? Otherwise returns vector of predictions
getSmooths <- function(smoothLabel,modList,xvals,postSamp=FALSE,noIntercept=TRUE){
  
  # #Debugging
  # smoothLabel <- 's(dist)'
  # postSamp <- FALSE
  # noIntercept <- TRUE
  # xvals <- seq(predRange[1],predRange[2],length.out=lengthOut) #x-values to use for smoother predictions
  
  smoothLabs <- sapply(modList$smooth,function(x) x$label) #Labels for smoothers
  meanSmoothList <- modList$smooth[[which(smoothLabs==smoothLabel)]] #Get smoother info
  
  #Get locations of coefficients
  meanVars <- c(which(grepl('Intercept',names(modList$coefs)))[1],meanSmoothList$first.para:meanSmoothList$last.para)
  
  if(postSamp){ #Sample from posterior
    meanCoefs <- rnorm(rep(1,length(meanVars)),modList$coefs[meanVars],sqrt(diag(modList$vcv)[meanVars]))
  } else { #Get ML estimate
    meanCoefs <- unname(modList$coefs[meanVars])
  }
  
  #Marginalize across intercept (set intercept coef to 0)
  if(noIntercept) meanCoefs[1] <- 0
  
  predDF <- data.frame(x=xvals) #Dataframe for holding predictions
  names(predDF) <- meanSmoothList$term #Name of term
  
  #Get predictions (coef %*% model matrix)
  predDF$pred <- cbind(rep(1,length(xvals)),PredictMat(meanSmoothList,data=predDF)) %*% meanCoefs
  
  return(predDF)
} 

#Function to extract predictions from modList.Rdata at specified path
# l = number of replicates in each smooth prediction (length.out)
# margInt = marginalize across intercepts
# samp = sample from posterior distribution of coefficients rather than using the mean
getPreds <- function(path,l=c(30,100,500),margInt=c(FALSE,FALSE,FALSE),samp=FALSE){
  
  # #Debugging
  # l <- c(30,100,500)
  # # path <- paste0('./Figures/ModelCheck/',datSource$filename[2],' modList.Rdata') #Model type 1
  # path <- paste0('./Figures/ModelCheck2/',datSource$filename[2],' modList.Rdata') #Model type 2
  # margInt <- c(FALSE,TRUE,TRUE)
  # samp <- FALSE
  
  require(mgcv)
  
  load(path) #Load data
  
  #Polygon area (basically speed) regression
  logPArea <- seq(log(modList$pAreaRange[1]),log(modList$pAreaRange[2]),length.out=l[1])
  pArea <- exp(logPArea) 
  
  meanVars <- which(grepl('(^\\(Intercept\\)$|^log\\(pArea\\)$)',names(modList$coefs)))
  sdVars <- which(grepl('(^\\(Intercept\\)\\.1$|^log\\(pArea\\)\\.1$)',names(modList$coefs)))
  
  #Coefficients (intercept and slope) for mean and logSD relationship
  if(samp){ #Sample from posterior
    meanCoefs <- rnorm(rep(1,length(meanVars)),modList$coefs[meanVars],sqrt(diag(modList$vcv)[meanVars]))
    sdCoefs <- rnorm(rep(1,length(meanVars)),modList$coefs[sdVars],sqrt(diag(modList$vcv)[sdVars]))
  } else {
    meanCoefs <- modList$coefs[meanVars] 
    sdCoefs <- modList$coefs[sdVars]
  }
  
  if(margInt[1]){ #Marginalize across intercept (set intercept coef to 0)
    meanCoefs[1] <- 0 
    sdCoefs[1] <- 0  
  }
  
  pAreaDat <- data.frame(pArea, 
                         mean=cbind(rep(1,length(logPArea)),logPArea) %*% meanCoefs,
                         logSD=cbind(rep(1,length(logPArea)),logPArea) %*% sdCoefs)
  
  # Field boundary distance smoothers
  
  #Figure out if model has only 1 type of distance (model type 1) or multiple distances (model type 2)
  #If first type, returns dataframe, otherwise a list of dataframes
  
  type1 <- !any(grepl('dist_',sapply(modList$smooths,function(x) x$label)))
  
  if(type1){
    dRange <- with(modList$smooths[[which(sapply(modList$smooths,function(x) x$label)=='s(dist)')]],
                   range(Xu)+rep(shift,2)) #Get range from smooths
    d <- seq(dRange[1],dRange[2],length.out=l[2]) #x values to sample from
    
    distDat <- data.frame(dist=d,
                          mean=getSmooths(smoothLabel='s(dist)', modList=modList, xvals = d,postSamp=samp, noIntercept=margInt[2])$pred,
                          logSD=getSmooths(smoothLabel='s.1(dist)', modList=modList, xvals = d,postSamp=samp, noIntercept=margInt[2])$pred)
  } else { #If using different types of distances
    
    #Get names of different distances
    distTypes <- unique(gsub('(\\)|s\\(|s.1\\()','',sapply(modList$smooths,function(x) x$label)))
    distTypes <- distTypes[grepl('dist',distTypes)]
    
    distDat <- lapply(distTypes,function(dType){
      meanLab <- paste0('s(',dType,')') #Gets appropriate text labels for smoothers
      sdLab <- paste0('s.1(',dType,')')
      
      dRange <- with(modList$smooths[[which(sapply(modList$smooths,function(x) x$label)==meanLab)]],
                     range(Xu)+rep(shift,2)) #Get range from smooths
      d <- seq(dRange[1],dRange[2],length.out=l[2]) #x values values to sample from
      distDat <- data.frame(dist=d, #
                            mean=getSmooths(smoothLabel=meanLab, modList=modList, xvals = d,postSamp=samp, noIntercept=margInt[2])$pred,
                            logSD=getSmooths(smoothLabel=sdLab, modList=modList, xvals = d,postSamp=samp, noIntercept=margInt[2])$pred)
      return(distDat)
    })
    names(distDat) <- distTypes  
    
  }
  
  #Point order (time of combining) smoothers
  rRange <- with(modList$smooths[[which(sapply(modList$smooths,function(x) x$label)=='s(r)')]],
                 range(Xu)+rep(shift,2)) #Get range from smooths
  r <- seq(rRange[1],rRange[2],length.out=l[3]) #r values to sample from
  
  rDat <- data.frame(r=r,
                        mean=getSmooths(smoothLabel='s(r)', modList=modList, xvals = r, postSamp=samp, noIntercept=margInt[3])$pred,
                        logSD=getSmooths(smoothLabel='s.1(r)', modList=modList, xvals = r, postSamp=samp, noIntercept=margInt[3])$pred)
  
  #Assemble into list
  datList <- list(pAreaDat=pAreaDat,distDat=distDat,rDat=rDat)
  
  return(datList)
}
# debugonce(getPreds)
# (temp <- getPreds(paste0('./Figures/ModelCheck/',datSource$filename[2],' modList.Rdata'),margInt=c(FALSE,TRUE,TRUE),samp=FALSE)) #Test with first type of model
# temp2 <- getSmooths(paste0('./Figures/ModelCheck/',datSource$filename[1],' modList.Rdata'),margInt=c(FALSE,TRUE,TRUE),samp=TRUE) #Test
# 
# ggplot()+geom_line(data=temp$pAreaDat,aes(x=pArea,y=mean))+ #Works for this (I think)
#   geom_line(data=temp2$pAreaDat,aes(x=pArea,y=mean),linetype='dashed',col='red')
# 
# ggplot()+geom_line(data=temp$dist,aes(x=dist,y=mean))+ #Works for this
#   geom_line(data=temp2$dist,aes(x=dist,y=mean),linetype='dashed',col='red')
# 
# ggplot()+geom_line(data=temp$r,aes(x=r,y=mean))+ #Works for this
#   geom_line(data=temp2$r,aes(x=r,y=mean),linetype='dashed',col='red')

# (temp <- getPreds(paste0('./Figures/ModelCheck2/',datSource$filename[2],' modList.Rdata'),margInt=c(FALSE,TRUE,TRUE),samp=FALSE)) #Test with second type of model
