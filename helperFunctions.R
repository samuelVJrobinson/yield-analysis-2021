# Helper functions

backTrans <- function(x){ #Back-transform units - data was sqrt transformed, so logSD is actually (log(sqrt(sd)))
  pow <- function(x,p) x^p
  expow <- function(x,p) log(exp(x)^p)
  # x$pArea$predMean <- pow(x$pArea$predMean,2)
  # x$pArea$predLogSD <- expow(x$pArea$predLogSD,2)
  
  x$coverDist <- lapply(x$coverDist,function(y){
    y$predMean <- pow(y$predMean,2)
    y$predLogSD <- expow(y$predLogSD,2)
    return(y)
  })
  
  # x$r$predMean <- pow(x$r$predMean,2)
  # x$r$predLogSD <- expow(x$r$predLogSD,2)
  
  return(x)
}

#Scales x between 0 and 1
range01 <- function(x) (x-min(x))/(max(x)-min(x))

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
  dat2 <- st_sf(st_drop_geometry(dat),geometry=st_sfc(polys)) %>% st_set_crs(datCRS)
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
getSmooths <- function(smoothLabel,modList,xvals,postSamp=FALSE,noIntercept=TRUE,returnSE=FALSE){
  
  # #Debugging
  # smoothLabel <- 's(dist)'
  # postSamp <- FALSE
  # noIntercept <- TRUE
  # xvals <- seq(0,100,length.out=101) #x-values to use for smoother predictions
  
  require(mgcv)
  
  smoothLabs <- sapply(modList$smooth,function(x) x$label) #Labels for smoothers
  meanSmoothList <- modList$smooth[[which(smoothLabs==smoothLabel)]] #Get smoother info
  
  #Get locations of coefficients
  whichIntercept <- ifelse(grepl('s.1',smoothLabel),2,1) #Is intercept for mean or logSD?
  meanVars <- c(which(grepl('Intercept',names(modList$coefs)))[whichIntercept],meanSmoothList$first.para:meanSmoothList$last.para)
  
  if(postSamp){ #Sample from posterior
    meanCoefs <- rnorm(rep(1,length(meanVars)),modList$coefs[meanVars],sqrt(diag(modList$vcv)[meanVars]))
  } else { #Get ML estimate
    meanCoefs <- unname(modList$coefs[meanVars])
  }
  
  #Marginalize across intercept (set intercept coef to 0)
  if(noIntercept) meanCoefs[1] <- 0
  
  if('matrix' %in% class(xvals)){ #If xvals is a matrix (2D smooth)
    predDF <- data.frame(x=xvals[,1],y=xvals[,2]) #Dataframe for holding predictions
  } else { #If not (1D smooth)
    predDF <- data.frame(x=xvals) #Dataframe for holding predictions
  }
  names(predDF) <- meanSmoothList$term #Name of terms
  
  if(meanSmoothList$by!='NA'){ #If using a by variable
    predDF$by <- factor(meanSmoothList$by.level)
    names(predDF)[2] <- meanSmoothList$by
  }
  
  #Model matrix
  predMat <- PredictMat(meanSmoothList,data=predDF) #Basis function columns
  predMat <- cbind(rep(1,nrow(predMat)),predMat) #Intercept column
  
  predDF$pred <-  predMat %*% meanCoefs  #Get predictions (coef %*% model matrix)
  
  if(returnSE){ #Gets SE of prediction at each xval - from plot.gam code
    predDF$se <- sqrt(pmax(0,rowSums(predMat %*% modList$vcv[meanVars,meanVars] * predMat)))
  }
  
  return(predDF)
} 

#Function to extract predictions from modList.Rdata at specified path
# l = number of replicates in each smooth prediction (length.out)
# margInt = marginalize across intercepts
# samp = sample from posterior distribution of coefficients rather than using the mean

getPreds <- function(path,l=c(30,100,500),margInt=c(FALSE,FALSE,FALSE),samp=FALSE,reps=1){
  
  # #Debugging
  # l <- c(30,100,500)
  # # path <- paste0('./Figures/ModelCheck/',datSource$filename[2],' modList.Rdata') #Model type 1
  # path <- paste0('./Figures/ModelCheck2/',datSource$filename[13],' modList.Rdata') #Model type 2
  # margInt <- c(FALSE,TRUE,TRUE)
  # samp <- FALSE
  
  require(mgcv)
  load(path) #Load data
  
  # #Polygon area (basically speed) regression
  # logPArea <- seq(log(modList$pAreaRange[1]),log(modList$pAreaRange[2]),length.out=l[1])
  # pArea <- exp(logPArea) 
  
  retList <- lapply(1:reps,function(i){
    
    #Indices for mean/sd intercepts
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
    
    # pAreaDat <- data.frame(pArea, 
    #                        mean=cbind(rep(1,length(logPArea)),logPArea) %*% meanCoefs,
    #                        logSD=cbind(rep(1,length(logPArea)),logPArea) %*% sdCoefs)
    pAreaDat <- NA
    
    # Field boundary distance smoothers
    
    #Figure out if model has only 1 type of distance (model type 1) or multiple distances (model type 2)
    #If first type, returns dataframe, otherwise a list of dataframes
    
    type1 <- !any(grepl('(dist_|\\:)',sapply(modList$smooths,function(x) x$label)))
    
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
        if(grepl('\\:',dType)){ #If smoother used "by" (rather than separate distance columns)
          meanLab <- gsub('dist','s(dist)',dType) #Gets appropriate text labels for smoothers
          sdLab <- gsub('dist','s.1(dist)',dType) #Gets appropriate text labels for smoothers
        } else {
          meanLab <- paste0('s(',dType,')') #Gets appropriate text labels for smoothers
          sdLab <- paste0('s.1(',dType,')')
        }
        
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
    
    # #Point order (time of combining) smoothers
    # rRange <- with(modList$smooths[[which(sapply(modList$smooths,function(x) x$label)=='s(r)')]],
    #                range(Xu)+rep(shift,2)) #Get range from smooths
    # r <- seq(rRange[1],rRange[2],length.out=l[3]) #r values to sample from
    # 
    # rDat <- data.frame(r=r,
    #                       mean=getSmooths(smoothLabel='s(r)', modList=modList, xvals = r, postSamp=samp, noIntercept=margInt[3])$pred,
    #                       logSD=getSmooths(smoothLabel='s.1(r)', modList=modList, xvals = r, postSamp=samp, noIntercept=margInt[3])$pred)
    rDat <- NA
    
    #Assemble into list
    datList <- list(pAreaDat=pAreaDat,distDat=distDat,rDat=rDat)
    return(datList)
  })
    
  if(reps==1){
    return(retList[[1]])
  } else {
    return(retList)
  }
  
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


#Function to run the ith model, using datSource ds, and sampling below nSubSamp points if necessary
#kPar = basis dimensions for distance, E/N, and time smoothers + basis dimensions for logSD smoothers
runModI <- function(i,dS,nSubSamp=50000,kPar=c(12,60,60,12,60,60),modelCheckDir='./Figures/ModelCheck1',resultsDir='./Figures/YieldMaps',filterData=TRUE){ 
  if(dS$modelComplete[i]) return('Already completed')
  if(!dS$use[i]) return('Not used')
  
  require(tidyverse)
  theme_set(theme_bw())
  require(mgcv)
  require(sf)
  require(ggpubr)
  
  source('helperFunctions.R')
  
  csvPath <- dS$dataPath[i] #Path to data
  fieldName <- dS$filename[i] #Field/year ID
  boundaryPath <- dS$boundaryPath2[i] #Path to boundary shapefile
  cropType <- dS$crop[i] #Crop type
  
  print('Reading in data')
  dat <- read.csv(csvPath,stringsAsFactors=TRUE,fileEncoding='latin1') 
  
  #Takes 40 seconds with subsamp of 50000 using full 800000 samples from Alvin French's Al Jr Field
  dat <- dat %>% rename_with(.fn = ~gsub('..L.ha.$','_lHa',.x)) %>%
    rename_with(.fn = ~gsub('..tonne.ha.$','_tHa',.x)) %>%
    rename_with(.fn = ~gsub('..m.s.$','_ms',.x)) %>%
    rename_with(.fn = ~gsub('.km.h.$','_kmh',.x)) %>%
    rename_with(.fn = ~gsub('..tonne.h.$','_th',.x)) %>%
    rename_with(.fn = ~gsub('.ha.h.','_hah',.x)) %>%
    rename_with(.fn = ~gsub('.m.','_m',.x)) %>%
    rename_with(.fn = ~gsub('.deg.','Angle',.x)) %>%
    rename_with(.fn = ~gsub('\\.','',.x)) %>%
    rename('ID'='ObjId','DryYield'='YldMassDry_tHa','Lon'='Longitude','Lat'='Latitude','Pass'='PassNum','Speed'='Speed__m') %>%
    st_as_sf(coords=c('Lon','Lat')) %>% #Add spatial feature info
    st_set_crs(4326) %>% st_transform(3401) %>% #Lat-lon -> UTM
    makePolys(width='SwthWdth_m',dist='Distance_m',angle='TrackAngle') %>%    
    mutate(pArea=as.numeric(st_area(.))) %>% #Area of polygon
    st_centroid() %>% #Convert back to point
    mutate(r=1:n()) %>% #row number
    mutate(Pass=factor(seqGroup(ID,FALSE))) %>% 
    group_by(Pass) %>% mutate(rGroup=1:n()) %>% ungroup() %>% 
    mutate(E=st_coordinates(.)[,1],N=st_coordinates(.)[,2]) %>% 
    mutate(E=E-mean(E),N=N-mean(N)) #Center coordinates
  
  print('Filtering data')  
  
  #Filter data
  dat <- dat %>% mutate(vegaFilt = vegaFilter(.,DryYield,nDist = 30)) #Vega et al 2019 spatial "inlier" filter - takes about a minute
  
  dat <- dat %>% 
    mutate(noBS = DryYield<ifelse(cropType=='Wheat',10.75,8)) %>% #JP recommends maximum filters of 160 bu (10.75 T/ha) for wheat, 120 bu (8 T/ha) for peas & canola
    mutate(Qfilt = QuantileFilter(DryYield,q=0.98)) %>% #Trim dry yield outliers
    mutate(bFilt = bearingFilter(TrackAngle,q=0.98)) %>% #Trim extreme bearing changes (turning)
    mutate(speedFilt = QuantileFilter(Speed,q=0.98)) %>% #Trim absolute speed outliers
    mutate(dSpeedFilt = dSpeedFilter(Speed,l=c(-2,-1,1,2),perc = 0.2)) %>% #Trim speed differences (>20% change 2 steps forward and backward, suggested by Lyle et al 2014)
    mutate(posFilt = posFilter(.,q=0.98)) %>% #Trim points that are far away from eachother
    #Combine filter criteria
    mutate(allFilt = noBS & vegaFilt & Qfilt & bFilt & speedFilt & dSpeedFilt & posFilt) %>%  
    filter(allFilt)
  
  if(nrow(dat)>nSubSamp){ #Limit to nSubSamp sequential samples if too large
    dat <- dat %>% slice(round(seq(1,nrow(dat),length.out=nSubSamp)))
  }
  
  print('Calculating distance from boundary')  
  
  fieldEdge <- read_sf(boundaryPath) %>% st_cast('MULTILINESTRING') #Read in boundary file
  # fieldEdge <- dat %>% st_union() %>% st_buffer(dist=10) %>% st_cast('MULTILINESTRING') #Old method
  
  dat <- dat %>%  #Distance from boundary
    mutate(dist=as.numeric(st_distance(.,fieldEdge))[1:nrow(.)]) 
  
  print('Fitting model')
  #Fit non-isotropic GAM:
  f <- sqrt(DryYield) ~ s(dist,k=kPar[1]) + s(E,N,k=kPar[2]) + s(r,k=kPar[3]) + log(pArea) #Mean model
  f2 <- ~ s(dist,k=kPar[4]) + s(E,N,k=kPar[5]) + s(r,k=kPar[6]) + log(pArea) #Variance model
  flist <- list(f,f2) #List of model formulae
  
  #Fit Gaussian location-scale  model
  #NOTE: this model can't be run using bam because of location-scale modeling
  
  a <- Sys.time() #Takes about 8 mins for a 50000 point model
  mod <- gam(flist,data=dat,family=gaulss())
  fitTime <- paste(as.character(round(Sys.time()-a,2)),units(Sys.time()-a)) #Time taken to fit model
  
  #Model with dist(10) E,N(60) r(30); s.1dist(6) s.1E,N(60) s.1r(60) - Has a very large amount of temporal/spatial autocorrelation, but only takes 6 mins or so
  #Trying with dist(12) E,N(80) r(80); s.1dist(12) s.1E,N(80) s.1r(80) - Similar temporal/spatial autocorrelation, but takes 4 times as long, and gives roughly the same answer
  
  print('Saving model results')
  
  #Plot results - FAILS HERE WHEN ALL.TERMS=TRUE. SOME KIND OF GAM NAMESPACE PROBLEM THAT ONLY OCCURS WITHIN FUNCTIONS
  #https://stackoverflow.com/questions/45918662/plot-gam-from-mgcv-with-all-terms-true-within-a-function
  png(paste0(modelCheckDir,'/',fieldName,' Summary.png'),width=12,height=8,units='in',res=200)
  # plot(mod,scheme=2,too.far=0.01,pages=1,all.terms=TRUE)
  plot(mod,scheme=2,too.far=0.01,pages=1,all.terms=FALSE)
  dev.off()
  
  capture.output({ #Model summary
    print("Time taken: "); print(fitTime)
    print(" ")
    print("SUMMARY---------------------------------")
    summary(mod) 
  },file=paste0(modelCheckDir,'/',fieldName,' results.txt'))
  
  capture.output({ #GAM check results
    print(" ")
    print("GAM.CHECK-------------------------------")
    png(paste0(modelCheckDir,'/',fieldName,' gamCheck.png'),width=8,height=8,units='in',res=200)
    par(mfrow=c(2,2)); gam.check(mod); abline(0,1,col='red'); par(mfrow=c(1,1))
    dev.off()
  },file=paste0(modelCheckDir,'/',fieldName,' results.txt'),append=TRUE) 
  
  #Actual, Predicted, and SE maps
  print('Making yield map')
  yieldMap <- dat %>% filter(DryYield<quantile(DryYield,0.95)) %>% #Chop out upper 5% of data (mainly crazy high numbers)
    ggplot() + 
    geom_sf(data=dat,col='grey',size=1)+
    geom_sf(aes(col=DryYield),alpha=0.3) + 
    geom_sf(data=fieldEdge,col='magenta')+ #Field boundary
    labs(title=paste(fieldName,'Data')) + #guides(alpha='none') + 
    scale_colour_distiller(type='div',palette = "Spectral") +
    theme(legend.position='bottom')
  # dat <- dat %>% mutate(yieldPred=mod$fit[,1],yieldVar=mod$fit[,2],resid=resid(mod))
  
  #Marginalize across r (point order) and pArea (speed)
  pred <- dat %>% 
    mutate(pArea=median(pArea)) %>% #Set pArea to its median value
    predict(mod,newdata = .,type='response',exclude='s(r)') #Exclude s(r)
  # mutate(r=1,pArea=median(pArea)) %>% #Sets r to 1, pArea to median
  # predict(mod,newdata = .,type='response')
  
  # pred2 <- dat %>% 
  #   # mutate(r=1,pArea=mean(pArea)) %>%
  #   predict(mod,newdata = .,se.fit=TRUE)
  # 
  # predVals <- cbind(sapply(pred,function(x) x[,1]),sapply(pred2,function(x) x[,1])) %>% as.data.frame() 
  # colnames(predVals) <- c('mean1','se1','mean2','se2')
  # 
  # predVals %>% pivot_longer(everything()) %>% mutate(par=ifelse(grepl('mean',name),'mean','se')) %>% 
  #   mutate(name=ifelse(grepl('1',name),'type1','type2')) %>% 
  #   ggplot()+geom_histogram(aes(x=value))+facet_grid(name~par,scales='free_x')
  
  dat <- dat %>% mutate(yieldMean=pred[,1],yieldSD=1/pred[,2],resid=resid(mod))
  
  #Predicted yield, marginalizing across pArea and r
  p1 <- ggplot(dat)+geom_sf(aes(col=yieldMean^2),alpha=0.3)+
    labs(title='Predicted yield | Combine Speed, Point Order',fill='Yield (T per Ha)',col='Yield (T per Ha)')+
    scale_colour_distiller(type='div',palette = "Spectral",direction=-1)+theme(legend.position='bottom')
  #Predicted SD
  p2 <- ggplot(dat)+geom_sf(aes(col=yieldSD^2),alpha=0.3)+
    labs(title='Yield SD | Combine Speed, Point Order',fill='SD Yield ',col='SD Yield')+
    scale_colour_distiller(type='div',palette = "Spectral",direction=-1)+theme(legend.position='bottom')
  p <- ggarrange(yieldMap,p1,p2,ncol=3)
  ggsave(paste0(resultsDir,'/',fieldName,'.png'),p,height=6,width=15,dpi=100)  
  
  #Residual plots/covariance plots
  p3 <- ggplot(dat)+geom_sf(aes(col=resid),alpha=0.5)+labs(title='Residuals (space)')+
    scale_colour_distiller(type='div',palette = "Spectral")+theme(legend.position='right')
  p4 <- ggplot(dat)+geom_point(aes(x=r,y=resid),alpha=0.5)+labs(title='Residuals (sequence)',x='Sequence',y='Residuals')+
    geom_hline(yintercept=0,col='red',linetype='dashed')
  
  #STFDF does ST empirical covariance plots well, but STIDF doesn't work for some reason
  # Also, there isn't any spatial replication across time, so not sure this matters
  library(gstat)
  library(spacetime) 
  dat_sp <- as_Spatial(select(dat,ID,resid))
  v <- variogram(resid~1, data=dat_sp,width=5,cutoff=300) #Variogram
  # vmod <- fit.variogram(v,vgm('Mat'),fit.kappa=TRUE) #Fit Matern covariance model
  # plot(v,vmod)
  p5 <- v %>% 
    ggplot()+geom_point(aes(x=dist,y=gamma))+
    labs(title='Residual variogram',x='Distance',y='Semivariance')
  # dat_variogram2 <- gstat::variogram(resid~1, #Used to look at non-stationary semivariance
  #                                   data=dat_sp,width=5,cutoff=150,alpha=c(0,45,90,135)) #North, NE, E, SW
  res_acf <- acf(dat$resid,lag.max=300,type='correlation',plot=FALSE)
  p6 <- with(res_acf,data.frame(acf=acf,lag=lag)) %>% 
    ggplot(aes(x=lag,y=acf))+geom_col()+
    labs(x='Time Lag',y='Autocorrelation',title='Residual autocorrelation')
  p <- ggarrange(p3,p4,p5,p6,ncol=2,nrow=2)
  p <- annotate_figure(p,top = text_grob(paste0(fieldName,' Residuals')))
  ggsave(paste0(modelCheckDir,'/',fieldName,' residuals.png'),p,height=8,width=10,dpi=100)  
  
  #Models are huge, so saving only key components (about 25% of the size)
  modList <- list(cropType=cropType, #Crop type
                  coefs=coef(mod), #Coefficients
                  vcv=vcov(mod), #Covariance matrix
                  smooths=mod$smooth, #Smooth specifications
                  distRange=range(dat$dist), #Range of distances
                  pAreaRange=range(dat$pArea), #Range of polygon areas
                  rRange=range(dat$r)) #Range of r
  save(modList,file=paste0(modelCheckDir,'/',fieldName,' modList.Rdata'))
  
  print('Done')
  gc()
  
  # #Example of reconstructing basis function matrix from new data
  # data.frame(dist=seq(0,100,5),
  #              pred=PredictMat(mod$smooth[[1]],data=data.frame(dist=seq(0,100,5))) %*% coef(mod)[c(mod$smooth[[1]]$first.para:mod$smooth[[1]]$last.para)]) %>%
  #   ggplot()+geom_point(aes(x=dist,y=pred))
  # 
  # x <- seq(0,100,5)
  # betas <- coef(mod)[c(mod$smooth[[1]]$first.para:mod$smooth[[1]]$last.para)] #Coefficients
  # p <- PredictMat(mod$smooth[[1]],data=data.frame(dist=seq(0,100,5))) #Matrix of basis functions
  # 
  # par(mfrow=c(2,1))
  # plot(0,0,type='n',xlim=range(x),ylim=range(p),xlab='Dist',ylab='y',main='Unscaled basis functions') #Unscaled basis functions
  # for(i in 1:ncol(p)) lines(x,p[,i],col=i)
  # #Scaled basis functions
  # plot(0,0,type='n',xlim=range(x),ylim=range(p*outer(rep(1,nrow(p)),betas)),xlab='Dist',ylab='y',main='Scaled basis functions')
  # for(i in 1:ncol(p)) lines(x,p[,i]*betas[i],col=i)
  # lines(x,p%*%betas,lwd=3)
  # par(mfrow=c(1,1))
}

#Function to run the ith model, but with boundary type
#useClosest = use only closest boundary? Otherwise, uses distance from all boundaries
runModII <- function(i,dS,nSubSamp=50000,kPar=c(5,80,5,80),useClosest=TRUE,modelCheckDir='./Figures/ModelCheck2',resultsDir='./Figures/YieldMaps2',filterData=TRUE){ 
  # i <- 295 #Debugging - "Trent_Clark W 34 2019"
  # i <- 5 #Alvin French Zoltan - 1.2 million points
  # dS <- datSource
  # nSubSamp <- 50000
  # useClosest <- TRUE
  # kPar <- c(12,60,60,12,60,60)
  # filterData <- TRUE
  
  if(dS$modelComplete2[i]) return('Already completed')
  if(!dS$use[i]) return('Not used')
  
  require(tidyverse)
  theme_set(theme_bw())
  require(mgcv)
  require(sf)
  require(ggpubr)
  
  source('helperFunctions.R')
  
  csvPath <- dS$dataPath[i] #Path to data
  fieldName <- dS$filename[i] #Field/year ID
  boundaryPath <- dS$boundaryPath2[i] #Path to boundary shapefile
  cropType <- dS$crop[i] #Crop type
  
  print('Reading in data')
  dat <- read.csv(csvPath,stringsAsFactors=TRUE,fileEncoding='latin1') 
  
  uprLim <- 150000 #Upper limit of 150,000 points (takes about 6 mins to do spatial inlier filtering)
  if(nrow(dat)>uprLim){
    dat <- dat %>% slice(round(seq(1,nrow(dat),length.out=uprLim)))
  }
  
  #Takes 40 seconds with subsamp of 50000 using full 800000 samples from Alvin French's Al Jr Field
  dat <- dat %>% rename_with(.fn = ~gsub('..L.ha.$','_lHa',.x)) %>%
    rename_with(.fn = ~gsub('..tonne.ha.$','_tHa',.x)) %>%
    rename_with(.fn = ~gsub('..m.s.$','_ms',.x)) %>%
    rename_with(.fn = ~gsub('.km.h.$','_kmh',.x)) %>%
    rename_with(.fn = ~gsub('..tonne.h.$','_th',.x)) %>%
    rename_with(.fn = ~gsub('.ha.h.','_hah',.x)) %>%
    rename_with(.fn = ~gsub('.m.','_m',.x)) %>%
    rename_with(.fn = ~gsub('.deg.','Angle',.x)) %>%
    rename_with(.fn = ~gsub('\\.','',.x)) %>%
    rename('ID'='ObjId','DryYield'='YldMassDry_tHa','Lon'='Longitude','Lat'='Latitude','Pass'='PassNum','Speed'='Speed__m') %>%
    st_as_sf(coords=c('Lon','Lat')) %>% #Add spatial feature info
    st_set_crs(4326) %>% st_transform(3401) %>% #Lat-lon -> UTM
    # makePolys(width='SwthWdth_m',dist='Distance_m',angle='TrackAngle') %>%    
    # mutate(pArea=as.numeric(st_area(.))) %>% #Area of polygon
    # st_centroid() %>% #Convert back to point
    mutate(pArea=SwthWdth_m*Distance_m) %>% #Area of polygon
    mutate(r=1:n()) %>% #row number
    mutate(Pass=factor(seqGroup(ID,FALSE))) %>% 
    group_by(Pass) %>% mutate(rGroup=1:n()) %>% ungroup() %>% 
    mutate(E=st_coordinates(.)[,1],N=st_coordinates(.)[,2]) %>% 
    mutate(E=E-mean(E),N=N-mean(N)) #Center coordinates
  
  print('Filtering data')  
  NrawDat <- nrow(dat)
  
  #Filter data
  dat <- dat %>% 
    #Vega et al 2019 spatial "inlier" filter - takes about a minute
    mutate(vegaFilt = vegaFilter(.,DryYield,nDist = 30)) 
  
  dat <- dat %>% 
    #BS filter - JP recommends maximum filters of 160 bu (10.75 T/ha) for wheat, 120 bu (8 T/ha) for peas & canola
    mutate(noBS = DryYield<ifelse(cropType=='Wheat',10.75,8)) %>% 
    mutate(Qfilt = QuantileFilter(DryYield,q=0.98)) %>% #Trim dry yield outliers
    mutate(bFilt = bearingFilter(TrackAngle,q=0.98)) %>% #Trim extreme bearing changes (turning)
    mutate(speedFilt = QuantileFilter(Speed,q=0.98)) %>% #Trim absolute speed outliers
    #Trim speed differences (>20% change 2 steps forward and backward, suggested by Lyle et al 2014)
    mutate(dSpeedFilt = dSpeedFilter(Speed,l=c(-2,-1,1,2),perc = 0.2)) %>% 
    mutate(posFilt = posFilter(.,q=0.98)) %>% #Trim points that are far away from eachother
    #Combine filter criteria
    mutate(allFilt = noBS & vegaFilt & Qfilt & bFilt & speedFilt & dSpeedFilt & posFilt) %>%  
    filter(allFilt)
  
  NfiltDat <- nrow(dat)
  
  if(NfiltDat>nSubSamp){ #Limit to nSubSamp sequential samples if too large
    dat <- dat %>% slice(round(seq(1,nrow(dat),length.out=nSubSamp)))
    print(paste0('Limiting analysis to ',nSubSamp,' data points'))
  }
  
  print('Calculating distance from boundary')  
  
  fieldEdge <- read_sf(boundaryPath) %>% st_cast('MULTILINESTRING') #Read in boundary file
  boundaryTypes <- unique(fieldEdge$type) #Get boundary types
  
  if(useClosest){
    #Get distance and type of closest boundary
    dat <- dat %>% 
      bind_cols(.,data.frame(t(apply(st_distance(dat,fieldEdge),1,function(x) c(dist=min(x,na.rm=TRUE),boundaryType=fieldEdge$type[which.min(x)]))))) %>% 
      mutate(dist=as.numeric(dist),boundaryType=factor(boundaryType))
  } else {
    #Get distance columns for each boundary type
    dist <- t(apply(st_distance(dat,fieldEdge),1,function(x) tapply(x,fieldEdge$type,min))) %>% 
      data.frame() %>% rename_with(~paste0('dist_',.x))
    dat <- dat %>% bind_cols(.,dist) #Bind distance columns to data
  }
  
  print('Fitting model')
  
  #Make model formula for non-isotropic gam
  if(useClosest){
    f2 <- paste0('~ s(dist,k=',kPar[1],',bs="ts",by=boundaryType) + s(E,N,k=',kPar[2],')')
    f <- paste0('sqrt(DryYield)~ s(dist,k=',kPar[3],',bs="ts",by=boundaryType) + s(E,N,k=',kPar[4],')')
  } else {
    f2 <- paste0('~ ',paste0('s(dist_',boundaryTypes,',k=',kPar[1],',bs="ts")',collapse=' + '),' + s(E,N,k=',kPar[2],')')
    f <- paste0('sqrt(DryYield) ~ ',paste0('s(dist_',boundaryTypes,',k=',kPar[3],',bs="ts")',collapse=' + '),' + s(E,N,k=',kPar[4],')')
  }
  
  # if(useClosest){ #These ones include a time smoother
  #   f2 <- paste0('~ s(dist,k=',kPar[1],',bs="ts",by=boundaryType) + s(E,N,k=',kPar[2],') + s(r,k=',kPar[3],')')
  #   f <- paste0('sqrt(DryYield)~ s(dist,k=',kPar[4],',bs="ts",by=boundaryType) + s(E,N,k=',kPar[5],') + s(r,k=',kPar[6],')')
  # } else {
  #   f2 <- paste0('~ ',paste0('s(dist_',boundaryTypes,',k=',kPar[1],',bs="ts")',collapse=' + '),' + s(E,N,k=',kPar[2],') + s(r,k=',kPar[3],')')
  #   f <- paste0('sqrt(DryYield) ~ ',paste0('s(dist_',boundaryTypes,',k=',kPar[4],',bs="ts")',collapse=' + '),' + s(E,N,k=',kPar[5],') + s(r,k=',kPar[6],')')
  # }

  flist <- list(as.formula(f),as.formula(f2)) #List of model formulae
  
  #Fit Gaussian location-scale  model
  a <- Sys.time() #Takes about 20 mins for a 50000 point model
  mod <- gam(flist,data=dat,family=gaulss())
  fitTime <- paste(as.character(round(Sys.time()-a,2)),units(Sys.time()-a)) #Time taken to fit model
  
  #Other models:
  # mod2 <- gam(as.formula(f),data=dat) #Standard model (constant variance)
  # f2 <- paste0('~ s(dist,k=',kPar[1],',bs="ts",by=boundaryType) + s(E,N,k=',kPar[2],')')
  # f <- paste0('sqrt(DryYield)~ s(dist,k=',kPar[4],',bs="ts",by=boundaryType) + s(E,N,k=',kPar[5],')')
  # flist <- list(as.formula(f),as.formula(f2)) #List of model formulae
  # mod3 <- gam(flist,data=dat,family=gaulss()) #Model without time smoothers
  # mod4 <- gam(as.formula(f),data=dat) #Standard model without time smoother
  
  print('Saving model results')
  
  png(paste0(modelCheckDir,'/',fieldName,' Summary.png'),width=12,height=8,units='in',res=200)
  # plot(mod,scheme=2,too.far=0.01,pages=1,all.terms=TRUE)
  plot(mod,scheme=2,too.far=0.01,pages=1,all.terms=FALSE)
  dev.off()
  
  capture.output({ #Model summary
    print("Time taken: "); print(fitTime)
    print(" ")
    print(paste0('Filtered out ',NrawDat-NfiltDat,' of ',NrawDat,' data points'))
    if(NfiltDat>nSubSamp) print(paste0('Limiting analysis to ',nSubSamp,' data points'))
    print(" ")
    print("SUMMARY---------------------------------")
    summary(mod) 
  },file=paste0(modelCheckDir,'/',fieldName,' results.txt'))
  
  capture.output({ #GAM check results
    print(" ")
    print("GAM.CHECK-------------------------------")
    png(paste0(modelCheckDir,'/',fieldName,' gamCheck.png'),width=8,height=8,units='in',res=200)
    par(mfrow=c(2,2)); gam.check(mod); abline(0,1,col='red'); par(mfrow=c(1,1))
    dev.off()
  },file=paste0(modelCheckDir,'/',fieldName,' results.txt'),append=TRUE) 
  
  #Actual, Predicted, and SE maps
  print('Making yield map')
  yieldMap <- dat %>% filter(DryYield<quantile(DryYield,0.95)) %>% #Chop out upper 5% of data (mainly crazy high numbers)
    ggplot() + 
    geom_sf(data=dat,col='grey',size=1)+
    geom_sf(aes(col=DryYield),alpha=0.3) + 
    geom_sf(data=fieldEdge,col='black')+ #Field boundary
    labs(title=paste(fieldName,'Raw Data'),colour='Yield (T/ha)') + #guides(alpha='none') + 
    scale_colour_distiller(type='div',palette = "Spectral") +
    theme(legend.position='bottom')
  
  #Marginalize across r (point order) 
  pred <- dat %>% 
    predict(mod,newdata = .,type='response',exclude='s(r)') #Exclude s(r)
  
  dat <- dat %>% mutate(yieldMean=pred[,1],yieldSD=1/pred[,2],resid=resid(mod))
  
  #Predicted yield, marginalizing across pArea and r
  p1 <- ggplot(dat)+geom_sf(aes(col=yieldMean^2),alpha=0.3)+
    labs(title='Predicted Mean | Point Order',fill='Yield (T/ha)',col='Yield (T/ha)')+
    scale_colour_distiller(type='div',palette = "Spectral",direction=-1)+
    geom_sf(data=fieldEdge,col='black')+ #Field boundary
    theme(legend.position='bottom')
  #Predicted SD
  p2 <- ggplot(dat)+geom_sf(aes(col=yieldSD^2),alpha=0.3)+
    labs(title='Predicted SD | Point Order',fill='SD Yield ',col='SD Yield')+
    scale_colour_distiller(type='div',palette = "Spectral",direction=-1)+
    geom_sf(data=fieldEdge,col='black')+ #Field boundary
    theme(legend.position='bottom')
  p <- ggarrange(yieldMap,p1,p2,ncol=3)
  ggsave(paste0(resultsDir,'/',fieldName,'.png'),p,height=6,width=15,dpi=100)  
  
  #Residual plots/covariance plots
  p3 <- ggplot(dat)+geom_sf(aes(col=resid),alpha=0.5)+labs(title='Residuals (space)')+
    scale_colour_distiller(type='div',palette = "Spectral")+theme(legend.position='right')
  p4 <- ggplot(dat)+geom_point(aes(x=r,y=resid),alpha=0.5)+labs(title='Residuals (sequence)',x='Sequence',y='Residuals')+
    geom_hline(yintercept=0,col='red',linetype='dashed')
  
  #STFDF does ST empirical covariance plots well, but STIDF doesn't work for some reason
  # Also, there isn't any spatial replication across time, so not sure this matters
  library(gstat)
  library(spacetime) 
  dat_sp <- as_Spatial(select(dat,ID,resid))
  v <- variogram(resid~1, data=dat_sp,width=5,cutoff=100) #Variogram
  vmod <- fit.variogram(v,vgm('Exp'),fit.kappa=TRUE) #Fit Matern covariance model
  # plot(v,vmod)
  p5 <- v %>% 
    ggplot()+geom_point(aes(x=dist,y=gamma))+
    geom_line(dat=variogramLine(vmod,maxdist = 100,min = 1,n=100),aes(x=dist,y=gamma))+
    labs(title='Residual variogram',x='Distance',y='Semivariance')
  # dat_variogram2 <- gstat::variogram(resid~1, #Used to look at non-stationary semivariance
  #                                   data=dat_sp,width=5,cutoff=150,alpha=c(0,45,90,135)) #North, NE, E, SW
  res_acf <- acf(dat$resid,lag.max=100,type='correlation',plot=FALSE)
  p6 <- with(res_acf,data.frame(acf=acf,lag=lag)) %>% 
    ggplot(aes(x=lag,y=acf))+geom_col()+
    labs(x='Time Lag',y='Autocorrelation',title='Residual autocorrelation')
  p <- ggarrange(p3,p4,p5,p6,ncol=2,nrow=2)
  ggsave(paste0(modelCheckDir,'/',fieldName,' residuals.png'),p,height=8,width=10,dpi=100)  
  
  #Models are huge, so saving only key components (about 25% of the size)
  
  if(useClosest){ #Range of distances
    distRange <- with(dat,tapply(Distance_m,boundaryType,range))
  } else {
    distRange <- st_drop_geometry(dat) %>% select(contains('dist_')) %>% summarize(across(everything(),range)) %>% data.frame()
  } 
  
  modList <- list(cropType=cropType, #Crop type
                  boundaryTypes=boundaryTypes, #Boundary types
                  coefs=coef(mod), #Coefficients
                  vcv=vcov(mod), #Covariance matrix
                  smooths=mod$smooth, #Smooth specifications
                  distRange=distRange,
                  pAreaRange=range(dat$pArea), #Range of polygon areas
                  rRange=range(dat$r)) #Range of r
  save(modList,file=paste0(modelCheckDir,'/',fieldName,' modList.Rdata'))
  
  print('Done')
  gc()
}

#Function to run the ith model, but with no boundary included ("null" model)
runMod0 <- function(i,dS,nSubSamp=50000,kPar=c(60,60,60,60),useClosest=TRUE,modelCheckDir='./Figures/ModelCheck0',resultsDir='./Figures/YieldMaps0',startCoefs=NULL,filterData=TRUE){ 
  i <- 1 #Debugging
  dS <- datSource
  nSubSamp <- 50000
  useClosest <- TRUE
  kPar <- c(60,60,60,60)
  startCoefs <- NULL
  startCoefs <- modList$coefs
  
  if(dS$modelComplete0[i]) return('Already completed')
  if(!dS$use[i]) return('Not used')
  
  require(tidyverse)
  theme_set(theme_bw())
  require(mgcv)
  require(sf)
  require(ggpubr)
  
  source('helperFunctions.R')
  
  csvPath <- dS$dataPath[i] #Path to data
  fieldName <- dS$filename[i] #Field/year ID
  boundaryPath <- dS$boundaryPath2[i] #Path to boundary shapefile
  cropType <- dS$crop[i] #Crop type
  
  print('Reading in data')
  dat <- read.csv(csvPath,stringsAsFactors=TRUE,fileEncoding='latin1') 
  
  #Takes 40 seconds with subsamp of 50000 using full 800000 samples from Alvin French's Al Jr Field
  dat <- dat %>% rename_with(.fn = ~gsub('..L.ha.$','_lHa',.x)) %>%
    rename_with(.fn = ~gsub('..tonne.ha.$','_tHa',.x)) %>%
    rename_with(.fn = ~gsub('..m.s.$','_ms',.x)) %>%
    rename_with(.fn = ~gsub('.km.h.$','_kmh',.x)) %>%
    rename_with(.fn = ~gsub('..tonne.h.$','_th',.x)) %>%
    rename_with(.fn = ~gsub('.ha.h.','_hah',.x)) %>%
    rename_with(.fn = ~gsub('.m.','_m',.x)) %>%
    rename_with(.fn = ~gsub('.deg.','Angle',.x)) %>%
    rename_with(.fn = ~gsub('\\.','',.x)) %>%
    rename('ID'='ObjId','DryYield'='YldMassDry_tHa','Lon'='Longitude','Lat'='Latitude','Pass'='PassNum','Speed'='Speed__m') %>%
    st_as_sf(coords=c('Lon','Lat')) %>% #Add spatial feature info
    st_set_crs(4326) %>% st_transform(3401) %>% #Lat-lon -> UTM
    makePolys(width='SwthWdth_m',dist='Distance_m',angle='TrackAngle') %>%    
    mutate(pArea=as.numeric(st_area(.))) %>% #Area of polygon
    st_centroid() %>% #Convert back to point
    mutate(r=1:n()) %>% #row number
    mutate(Pass=factor(seqGroup(ID,FALSE))) %>% 
    group_by(Pass) %>% mutate(rGroup=1:n()) %>% ungroup() %>% 
    mutate(E=st_coordinates(.)[,1],N=st_coordinates(.)[,2]) %>% 
    mutate(E=E-mean(E),N=N-mean(N)) #Center coordinates
  
  print('Filtering data')  
  
  #Filter data
  dat <- dat %>% mutate(vegaFilt = vegaFilter(.,DryYield,nDist = 30)) #Vega et al 2019 spatial "inlier" filter - takes about a minute
  
  dat <- dat %>% 
    mutate(noBS = DryYield<ifelse(cropType=='Wheat',10.75,8)) %>% #JP recommends maximum filters of 160 bu (10.75 T/ha) for wheat, 120 bu (8 T/ha) for peas & canola
    mutate(Qfilt = QuantileFilter(DryYield,q=0.98)) %>% #Trim dry yield outliers
    mutate(bFilt = bearingFilter(TrackAngle,q=0.98)) %>% #Trim extreme bearing changes (turning)
    mutate(speedFilt = QuantileFilter(Speed,q=0.98)) %>% #Trim absolute speed outliers
    mutate(dSpeedFilt = dSpeedFilter(Speed,l=c(-2,-1,1,2),perc = 0.2)) %>% #Trim speed differences (>20% change 2 steps forward and backward, suggested by Lyle et al 2014)
    mutate(posFilt = posFilter(.,q=0.98)) %>% #Trim points that are far away from eachother
    #Combine filter criteria
    mutate(allFilt = noBS & vegaFilt & Qfilt & bFilt & speedFilt & dSpeedFilt & posFilt) %>%  
    filter(allFilt)
  
  if(nrow(dat)>nSubSamp){ #Limit to nSubSamp sequential samples if too large
    dat <- dat %>% slice(round(seq(1,nrow(dat),length.out=nSubSamp)))
  }
  
  print('Fitting model')
  
  #Make model formula for non-isotropic gam
  f2 <- paste0('~ s(E,N,k=',kPar[1],') + s(r,k=',kPar[2],') + log(pArea)')
  f <- paste0('sqrt(DryYield)~ s(E,N,k=',kPar[3],') + s(r,k=',kPar[4],') + log(pArea)')
  
  flist <- list(as.formula(f),as.formula(f2)) #List of model formulae
  
  #Fit Gaussian location-scale  model
  
  
  a <- Sys.time() 
  mod <- gam(flist,data=dat,family=gaulss(),start=startCoefs)
  fitTime <- paste(as.character(round(Sys.time()-a,2)),units(Sys.time()-a)) #Time taken to fit model
  
  print('Saving model results')
  
  png(paste0(modelCheckDir,'/',fieldName,' Summary.png'),width=12,height=8,units='in',res=200)
  # plot(mod,scheme=2,too.far=0.01,pages=1,all.terms=TRUE)
  plot(mod,scheme=2,too.far=0.01,pages=1,all.terms=FALSE)
  dev.off()
  
  capture.output({ #Model summary
    print("Time taken: "); print(fitTime)
    print(" ")
    print("SUMMARY---------------------------------")
    summary(mod) 
  },file=paste0(modelCheckDir,'/',fieldName,' results.txt'))
  
  capture.output({ #GAM check results
    print(" ")
    print("GAM.CHECK-------------------------------")
    png(paste0(modelCheckDir,'/',fieldName,' gamCheck.png'),width=8,height=8,units='in',res=200)
    par(mfrow=c(2,2)); gam.check(mod); abline(0,1,col='red'); par(mfrow=c(1,1))
    dev.off()
  },file=paste0(modelCheckDir,'/',fieldName,' results.txt'),append=TRUE) 
  
  #Actual, Predicted, and SE maps
  print('Making yield map')
  yieldMap <- dat %>% filter(DryYield<quantile(DryYield,0.95)) %>% #Chop out upper 5% of data (mainly crazy high numbers)
    ggplot() + 
    geom_sf(data=dat,col='grey',size=1)+
    geom_sf(aes(col=DryYield),alpha=0.3) + 
    geom_sf(data=fieldEdge,col='magenta')+ #Field boundary
    labs(title=paste(fieldName,'Data')) + #guides(alpha='none') + 
    scale_colour_distiller(type='div',palette = "Spectral") +
    theme(legend.position='bottom')
  
  #Marginalize across r (point order) and pArea (speed)
  pred <- dat %>% 
    mutate(pArea=median(pArea)) %>% #Set pArea to its median value
    predict(mod,newdata = .,type='response',exclude='s(r)') #Exclude s(r)
  
  dat <- dat %>% mutate(yieldMean=pred[,1],yieldSD=1/pred[,2],resid=resid(mod))
  
  #Predicted yield, marginalizing across pArea and r
  p1 <- ggplot(dat)+geom_sf(aes(col=yieldMean^2),alpha=0.3)+
    labs(title='Predicted yield | Combine Speed, Point Order',fill='Yield (T per Ha)',col='Yield (T per Ha)')+
    scale_colour_distiller(type='div',palette = "Spectral",direction=-1)+theme(legend.position='bottom')
  #Predicted SD
  p2 <- ggplot(dat)+geom_sf(aes(col=yieldSD^2),alpha=0.3)+
    labs(title='Yield SD | Combine Speed, Point Order',fill='SD Yield ',col='SD Yield')+
    scale_colour_distiller(type='div',palette = "Spectral",direction=-1)+theme(legend.position='bottom')
  p <- ggarrange(yieldMap,p1,p2,ncol=3)
  ggsave(paste0(resultsDir,'/',fieldName,'.png'),p,height=6,width=15,dpi=100)  
  
  #Residual plots/covariance plots
  p3 <- ggplot(dat)+geom_sf(aes(col=resid),alpha=0.5)+labs(title='Residuals (space)')+
    scale_colour_distiller(type='div',palette = "Spectral")+theme(legend.position='right')
  p4 <- ggplot(dat)+geom_point(aes(x=r,y=resid),alpha=0.5)+labs(title='Residuals (sequence)',x='Sequence',y='Residuals')+
    geom_hline(yintercept=0,col='red',linetype='dashed')
  
  #STFDF does ST empirical covariance plots well, but STIDF doesn't work for some reason
  # Also, there isn't any spatial replication across time, so not sure this matters
  library(gstat)
  library(spacetime) 
  dat_sp <- as_Spatial(select(dat,ID,resid))
  v <- variogram(resid~1, data=dat_sp,width=5,cutoff=300) #Variogram
  # vmod <- fit.variogram(v,vgm('Mat'),fit.kappa=TRUE) #Fit Matern covariance model
  # plot(v,vmod)
  p5 <- v %>% 
    ggplot()+geom_point(aes(x=dist,y=gamma))+
    labs(title='Residual variogram',x='Distance',y='Semivariance')
  # dat_variogram2 <- gstat::variogram(resid~1, #Used to look at non-stationary semivariance
  #                                   data=dat_sp,width=5,cutoff=150,alpha=c(0,45,90,135)) #North, NE, E, SW
  res_acf <- acf(dat$resid,lag.max=300,type='correlation',plot=FALSE)
  p6 <- with(res_acf,data.frame(acf=acf,lag=lag)) %>% 
    ggplot(aes(x=lag,y=acf))+geom_col()+
    labs(x='Time Lag',y='Autocorrelation',title='Residual autocorrelation')
  p <- ggarrange(p3,p4,p5,p6,ncol=2,nrow=2)
  p <- annotate_figure(p,top = text_grob(paste0(fieldName,' Residuals')))
  ggsave(paste0(modelCheckDir,'/',fieldName,' residuals.png'),p,height=8,width=10,dpi=100)  
  
  #Models are huge, so saving only key components (about 25% of the size)
  modList <- list(cropType=cropType, #Crop type
                  coefs=coef(mod), #Coefficients
                  vcv=vcov(mod), #Covariance matrix
                  ll=logLik(mod), #Log-likelihood
                  aic=AIC(mod), #AIC
                  rmse=sqrt(mean(resid(mod)^2)), #Root mean square error
                  mae=mean(abs(resid(mod))), #Mean absolute error
                  smooths=mod$smooth, #Smooth specifications
                  pAreaRange=range(dat$pArea), #Range of polygon areas (pArea)
                  rRange=range(dat$r)) #Range of r
  save(modList,file=paste0(modelCheckDir,'/',fieldName,' modList.Rdata'))
  
  print('Done')
  gc()
}


#Get model info from saved text file at path
getModelInfo <- function(path){
  if(!file.exists(path)){
    warning('File does not exist')
    retList <- list(timeTaken=NA,
                    percDevianceExplained=NA,
                    REML=NA,AIC=NA,
                    hessian=NA,
                    params=NA,
                    smooths=NA,
                    gamCheck=NA)
    return(retList)
  }
  
  txt <- readLines(path)
  
  require(tidyverse)
  #Get time taken to run model
  timeLine <- which(grepl('Time taken',txt))+1
  if(length(timeLine)==0){
    timeTaken <- NA
  } else {
    timeTaken <- strsplit(txt[timeLine],'\"')[[1]][2] #Time taken
    timeTaken <- strsplit(timeTaken,' ')[[1]]
    timeTaken <- as.difftime(as.numeric(timeTaken[1]),units=timeTaken[2])
  }
  
  #Get deviance explained
  devLine <- which(grepl('Deviance explained =',txt))
  if(length(devLine)==0){
    pExpDev <- NA
  } else {
    pExpDev <- txt[devLine]
    pExpDev <- strsplit(pExpDev,c(' '))[[1]]
    pExpDev <- as.numeric(gsub('%','',pExpDev[length(pExpDev)]))
  }
  
  #Get REML
  remlLine <- which(grepl('-REML =',txt))
  if(length(remlLine)==0){
    reml <- NA
  } else {
    reml <- strsplit(txt[remlLine],' ')[[1]]
    reml <- reml[reml!='']
    reml <- as.numeric(reml[3])  
  }
  
  #Check Hessian matrix
  hessLine <- which(grepl('Hessian positive definite, eigenvalue range',txt))
  hessCheck <- txt[hessLine]

  #Get parametric coefs, SEs, p-vals...
  # smooth edfs, p-vals
  # gam.check values
  
  #Starting lines
  parStartLine <- which(grepl('Parametric coefficients:',txt)) 
  smStartLine <- which(grepl('Approximate significance of smooth terms:',txt))
  gcStartLine <- which(grepl('Basis dimension \\(k\\) checking results.',txt))
  threeDash <- which(txt=='---') #Separators
  #Ending lines
  parEndLine <- threeDash[threeDash<smStartLine]
  smEndLine <- threeDash[threeDash>smStartLine&threeDash<gcStartLine] 
  gcEndLine <- threeDash[threeDash>gcStartLine] 
  
  eraseChars <- function(x){
    x <- gsub('(^\\.$|\\*|<)','',x) #Erases single . , *, and < characters
    x <- x[x!=''] #Drops empty vectors
    return(x)
  }
  
  if(is.na(pExpDev)){
    parTxt <- NA
    smTxt <- NA
  } else {
    #Linear terms
    parTxt <- strsplit(txt[(parStartLine+2):(parEndLine-1)],' ') %>% 
      lapply(.,eraseChars) %>% 
      do.call('rbind',.) %>% data.frame() %>% 
      set_names(c('Parameter','Estimate','SE','Chisq','pval')) %>% 
      mutate(across(Estimate:Chisq,as.numeric))
    
    #Smooth terms
    smTxt <- strsplit(txt[(smStartLine+2):(smEndLine-1)],' ') %>% 
      lapply(.,eraseChars) %>% 
      do.call('rbind',.) %>% data.frame() %>% 
      set_names(c('Smoother','edf','Ref.df','Chisq','pval')) %>% 
      mutate(across(edf:Chisq,as.numeric))
  }
  
  #GAM check results
  gcTxt <- strsplit(txt[(gcStartLine+4):(gcEndLine-1)],' ') %>% 
    lapply(.,eraseChars) %>% 
    do.call('rbind',.) %>% data.frame() %>% 
    set_names(c('Smoother','k','edf','index','pval')) %>% 
    mutate(across(k:index,as.numeric))
  
  #Degrees of freedom (# of linear params + sum(edf) from smooth params)
  #Can be viewed as "marginal" AIC
  Ndf <- nrow(parTxt)+sum(smTxt$edf)
  aic <- (-2*reml)+(2*Ndf)
  
  #List of things to return
  retList <- list(timeTaken=timeTaken,
                  percDevianceExplained=pExpDev,
                  REML=reml,AIC=aic,
                  hessian=hessCheck,
                  params=parTxt,
                  smooths=smTxt,
                  gamCheck=gcTxt)
  
  return(retList)
}
# getModelInfo("/home/rsamuel/Documents/yield-analysis-2021/Figures/ModelCheck0/Alvin_French Al_Jr 2020 results.txt")
# getModelInfo("/home/rsamuel/Documents/yield-analysis-2021/Figures/ModelCheck1/Alvin_French Al_Jr 2020 results.txt")$REML
# getModelInfo("/home/rsamuel/Documents/yield-analysis-2021/Figures/ModelCheck2/Alvin_French Al_Jr 2020 results.txt")$REML
# 
# debugonce(getModelInfo)
# getModelInfo("./Figures/ModelCheck/Dean Hubbard 2018 All_27_11_25 results.txt")


#Yield data filters -----------------

#"Inlier" spatial filtering procedure from Vega et al 2019, written to work with sf + dplyr. Returns boolean
vegaFilter <- function(data,ycol,pvalCutoff=0.05,nDist=40){
  require(sf)
  require(sp)
  require(spdep)
  require(tidyverse)
  
  # #Debugging
  # nDist <- 30
  # yield <- data$DryYield
  
  if(!any(class(data) %in% 'sf')) stop('Dataframe must be sf object')
  
  coords <- st_coordinates(data) #Get coordinates in matrix form
  
  #Get neighbourhood weights from 0 to ndist meters
  nWeights <- dnearneigh(coords,0,nDist) 
  
  if(any(sapply(nWeights,length)==1)) warning('Some points had no neighbours and were removed')
  
  #Get neighbourhood indices for each point (which other points are in this point's neighbourhood?)
  nIndices <- nb2listw(nWeights, style = "W",zero.policy = TRUE) 
  
  yield <- pull(data,{{ycol}}) #Get yield data column
  
  #Local Moran's I
  LM <- localmoran(yield,nIndices,p.adjust.method="bonferroni",alternative ="less")
  
  # # Moran Plot of Yield data against spatially lagged values - not needed for filter, only for plotting
  # MP <- moran.plot(yield,nIndices,quiet=TRUE,labels=TRUE,zero.policy=FALSE,
  #                  xlab='Yield', ylab="Spatially Lagged")
  # results <- data.frame(LM,MP) #Joins into single dataframe
    
  results <- data.frame(LM) %>% 
    rename('pval'=contains('Pr.z.')) %>% 
    mutate(keepThese= Ii > 0 | pval > pvalCutoff) %>% #Filter negative Ii and pvals < 0.05
    mutate(keepThese=ifelse(is.na(keepThese),FALSE,keepThese)) #NAs (points had no neighbours)
  
  ret <- pull(results)
  return(ret)
}

ZscoreFilter <- function(x,z=3){
  xmean <- mean(x,na.rm=TRUE)
  xsd <- sd(x,na.rm=TRUE)
  zscore <- abs(((x-xmean)/xsd))
  zscore<z
} #Function to filter anything above a certain z-score

QuantileFilter <- function(x,quant=0.99){ #Function to filter anything above certain quantiles
  l <- c((1-quant)/2,1-(1-quant)/2) #Symmetric quantiles
  x>quantile(x,l[1]) & x<quantile(x,l[2]) 
}

bearingFilter <- function(bearing,q=NULL,z=NULL,returnDiffs=FALSE){
  if(!xor(is.null(q),is.null(z))&!returnDiffs) stop('Input quantiles or Z-score')
  #Difference in compass bearings (in degrees)
  bearingDiff <- function(x1,x2){
    x <- x1-x2
    x <- ifelse(abs(x)>180,x-(360*sign(x)),x) #Angle differences can't be >180
    return(x)
  }

  #Looks 1 point ahead and behind
  bd <- cbind(bearingDiff(lag(bearing),bearing),
              bearingDiff(lead(bearing),bearing))
  bd <- apply(bd,1,function(x) max(abs(x),na.rm=TRUE)*sign(x[which.max(abs(x))])) #Maximum bearing difference ahead and behind
  
  if(returnDiffs) return(bd) #Return bearing differences only, without filtering
  
  if(!is.null(q)){
    ret <- QuantileFilter(bd,q=q)
  } else {
    ret <- ZscoreFilter(bd,z=z)
  }
  return(ret)
} 


posFilter <- function(data,q=NULL,returnDiffs=FALSE){ #Positional difference filter - filters out very distant and very close points
  if(is.null(q)&!returnDiffs) stop('Input upper quantile')
  # #Debugging
  # data <- dat
  # q <- 0.99
  
  if(units(st_distance(data[1,],data[2,]))$numerator!='m') warning('Position differences not in meters')
  coords <- st_coordinates(data) #Get coordinates
  pdiff <- sapply(1:(nrow(coords)-1),function(i) as.numeric(dist(coords[i:(i+1),]))) #Distances between points
  pdiff <- cbind(c(pdiff,NA),c(NA,pdiff)) #Forward and backward lags
  pdiff <- apply(pdiff,1,max,na.rm=TRUE) #Maximum distance ahead and behind
  if(returnDiffs){
    return(pdiff)
  } else {
    return(QuantileFilter(pdiff,q)) #Note: uses 2-sided quantiles
  }
}


dSpeedFilter <- function(speed,l=c(-1,1),perc=0.2){ #Filter for (forward and backward) lagged speed differences.
  # #Debugging
  # speed <- dat$Speed[1:30]
  # l <- c(2,1,-1,-2) #Backward and forward lags
  # perc <- 0.2
  
  llist <- sapply(l,function(x) (lag2(speed,x)-speed)/lag2(speed,x)) #Matrix of % diffs
  
  #Are any lagged speed values > percent change threshold?
  ret <- !apply(llist,1,function(y) any(abs(y)[!is.na(y)]>perc))
  
  return(ret)
}

#Overloaded lag function that takes negative values
lag2 <- function(x,n){
  if(n==0) {
    return(x) #No lag
  } else if(n>0){
    lag(x,n) #Positive lag
    } else {
      lead(x,abs(n)) #Negative lag
    } 
} 

# Function to sample from posterior of each prediction, amalgamate, fit meta-model, return results of meta-model
# My way of getting around the problem of random effects
# a = arguments from parLapply, not used
# useRows = rows of ds to be used, ds = datSource csv, nX = number of samples along range of X
# rCutoff = restrict point number (some fields have hundreds of thousands of points, but all data uses approximately the same size yield rectangles, so cuts off meta-smoother at "about" the maximum field size - i.e. the most common number of points per field)

#Note: this is more memory-intensive, but quicker. Looks like it maxes out around 12 GB
samplePreds <- function(a=NA,useRows=NULL,ds=datSource,nX=c(30,100,500),rCutoff=50000,margInt=c(FALSE,FALSE,FALSE),Nsamp=1,ncore=NA,silent=TRUE,samp=TRUE,kPar=5){ 
  # #Debugging
  # useRows <- which(datSource$crop=='Canola')
  # ds <- datSource
  # nX=c(30,100,500)
  # rCutoff=50000
  # margInt=c(FALSE,FALSE,FALSE)
  # Nsamp <- 10
  # ncore <- 5
  
  if(is.null(useRows)){
    warning('Rows not specified. Using entire dataset.')
    useRows <- 1:nrow(datSource)
  } 
  require(tidyverse); require(mgcv); require(lme4)
  source('helperFunctions.R')
  getFiles <- ds %>% slice(useRows) %>% filter(use,modelComplete2) #Completed models
  
  if(!silent) print('Sampling...')
  #Sample from posterior, get new lines for each field
  if(is.na(ncore)){
    allSmooths <- lapply(getFiles$modelPath2,getPreds,l=nX,margInt=margInt,samp=samp,reps=Nsamp) %>% set_names(getFiles$filename) #Serial version
    if(Nsamp==1) allSmooths <- lapply(allSmooths,function(x){ret <- list(x)})
  } else { #Parallel version
    require(parallel)
    cluster <- makeCluster(ncore) #Memory usage is OK, so could probably max it out
    clusterExport(cluster,c('datSource','getSmooths'))
    allSmooths <- parLapply(cl=cluster,getFiles$modelPath2,getPreds,l=nX,margInt=margInt,samp=samp,reps=Nsamp) %>% 
      set_names(getFiles$filename) 
  }
  gc()
  
  # #log(Area) meta-models
  # pAreaDf <- lapply(allSmooths,function(x) x$pAreaDat) %>% bind_rows(.id = 'field')
  # 
  # #Adds a small amount of noise so model fit isn't singular
  # m1 <- pAreaDf %>% 
  #   mutate(mean=mean+rnorm(n(),0,(max(mean)-min(mean))*0.005)) %>% 
  #   lmer(mean~log(pArea)+(log(pArea)|field),data=.) 
  # m2 <- pAreaDf %>% mutate(logSD=logSD+rnorm(n(),0,(max(logSD)-min(logSD))*0.005)) %>% 
  #   lmer(logSD~log(pArea)+(log(pArea)|field),data=.)
  
  #cover class meta-models
  
  if(!silent) cat('Fitting meta-models...')
  
  #Get cover class names
  coverNames <- sort(unique(unlist(lapply(allSmooths, function(x) names(x[[1]]$distDat)))))
  
  #Which models have which cover classes
  containsCover <- sapply(coverNames,function(z) sapply(lapply(allSmooths, function(x) names(x[[1]]$distDat)), function(y) z %in% y))
  
  #Empty list for results  
  coverPred <- vector(mode='list',length=length(coverNames)) %>% set_names(gsub('dist:boundaryType','',coverNames))
  
  #Helper function to fit smoothers from each cover class
  metaSmoothFun <- function(j,useThese,cname,allSmooths,nX){
    require(dplyr); require(mgcv)
    coverDF <- lapply(allSmooths[useThese],function(x) x[[j]]$distDat[[cname]]) %>% 
      bind_rows(.id = 'field') %>% mutate(field=factor(field)) #Get all field smoothers from Nsamp repetition
    meanMod <- gam(mean~s(dist,k=kPar)+s(field,bs='re'),data=coverDF) #Fit meta-models
    logSDMod <- gam(logSD~s(dist,k=kPar)+s(field,bs='re'),data=coverDF)
    retDF <- data.frame(dist=seq(min(coverDF$dist),max(coverDF$dist),length.out=nX[2]),field=coverDF$field[1]) %>% #Get predictions
      mutate(predMean=predict(meanMod,newdata=.,exclude="s(field)"),predLogSD=predict(logSDMod,newdata=.,exclude="s(field)")) %>% 
      select(-field)
    rm(coverDF,meanMod,logSDMod); gc() #Cleanup
    return(retDF) #Return DF
  }
  
  for(i in 1:ncol(containsCover)){ #For each cover type
    if(!silent) cat(paste0(names(coverPred)[i],' '))
    useThese <- unname(which(containsCover[,i])) #Which fields have this cover type?
    cname <- colnames(containsCover)[i]
    options(dplyr.summarise.inform = FALSE)
    
    if(Nsamp==1){ #Single replicate
      coverPred[[i]] <- metaSmoothFun(1,useThese=useThese,cname=cname,allSmooths=allSmooths,nX=nX)  
    } else if(is.na(ncore)){ #Serial
      coverPred[[i]] <- lapply(1:Nsamp,metaSmoothFun,useThese=useThese,cname=cname,allSmooths=allSmooths,nX=nX)  
    } else { #Parallel
      coverPred[[i]] <- parLapply(cl=cluster,1:Nsamp,metaSmoothFun,useThese=useThese,cname=cname,allSmooths=allSmooths,nX=nX)  
    }
    # set_names(nm=paste0('samp',1:Nsamp)) %>% bind_rows(.id = 'samp') %>%  #Done outside of function
    # pivot_longer(cols=predMean:predLogSD) %>% 
    # group_by(dist,name) %>% summarize(med=median(value),upr=quantile(value,0.9),lwr=quantile(value,0.1)) %>% #Get median, 90th and 10th percentiles
    # pivot_wider(id_cols = dist,values_from = med:lwr) %>% data.frame()
  }
  if(!silent) cat('Done')
  if(!is.na(ncore)) stopCluster(cluster)
  
  # #r (point order) meta-models
  # rDatDf <- lapply(allSmooths,function(x) x$rDat) %>% 
  #   bind_rows(.id = 'field') %>% mutate(field=factor(field)) %>% 
  #   filter(r<=rCutoff) #Filters out point values above cutoff
  # 
  # r1 <- gam(mean~s(r,k=50)+s(field,bs='re'),data=rDatDf)
  # r2 <- gam(logSD~s(r,k=50)+s(field,bs='re'),data=rDatDf)
  
  # #Assemble predictions
  # pAreaPred <- data.frame(pArea=exp(seq(log(min(pAreaDf$pArea)),log(max(pAreaDf$pArea)),length.out=nX[1]))) %>% 
  #   mutate(predMean=predict(m1,newdata=.,re.form=~0),predLogSD=predict(m2,newdata=.,re.form=~0))
  pAreaPred <- NA
  
  # rPred <- data.frame(r=seq(min(rDatDf$r),max(rDatDf$r),length.out=nX[3]),field=rDatDf$field[1]) %>% 
  #   mutate(predMean=predict(r1,newdata=.,exclude="s(field)"),predLogSD=predict(r2,newdata=.,exclude="s(field)"))
  rPred <- NA
  
  retList <- list(pArea=pAreaPred,coverDist=coverPred,r=rPred) #Return list of results
  
  rm(list=ls()[!ls() %in% 'retList']) #Remove everything except retList
  gc()
  return(retList) 
  
  # #Code to produce figures of field-level smoothers  
  # fieldSmooths <- lapply(allSmooths,function(x) x[[1]]$distDat %>% bind_rows(.id='type') %>% mutate(type=gsub('dist:boundaryType','',type))) %>% 
  #   bind_rows(.id='field')
  # 
  # meanSmooth <- bind_rows(coverPred,.id = 'type') #Mean smoothers
  # 
  # ggplot(fieldSmooths)+
  #   geom_line(aes(x=dist,y=mean,group=field))+
  #   facet_wrap(~type,nrow=1)+
  #   coord_cartesian(ylim=c(-1,1),xlim=c(0,300))
}
# debugonce(samplePreds) #Test
# isCanola <- which(datSource$crop=='Canola') #Canola crops only
# samplePreds(useRows=isCanola,Nsamp = 1) #Test
