#Run non-isotropic GAMs on multiple field years

# Load required libraries ---------------------------------------------------------

library(tidyverse)
theme_set(theme_bw())
library(mgcv)
library(sf)
library(beepr)

rootPath <- "/media/rsamuel/Storage/geoData/Rasters/yieldData/csv files"
datSource <- data.frame(path=dir(rootPath,pattern=".csv",recursive=TRUE)) %>% 
  separate(path,c('grower','year','field'),sep="/",remove=FALSE) %>% 
  mutate(field=gsub('\\.csv','',field)) %>% unite(filename,c(grower:field),sep=' ',remove=FALSE) %>% 
  mutate(completed=filename %in% gsub(' modList.Rdata','',dir('./Figures/ModelCheck',pattern=".Rdata",recursive=TRUE))) #Have files already been processed?

datSource %>% count(grower,year) 

set.seed(1)
datSource <- datSource %>%
  # filter(grower!='Alvin French') %>% #These fields are huge, so leave out for now
  slice_sample(n=12) #Smaller set to experiment with

#Function to run the ith model
runModI <- function(i,dS,rP,nSubSamp=50000){ 
  i <- 3 #Debugging
  dS <- datSource
  rP <- rootPath
  nSubSamp <- 50000
  
  if(dS$completed[i]) return('Already completed')
  
  library(tidyverse)
  theme_set(theme_bw())
  library(mgcv)
  library(sf)
  
  source('helperFunctions.R')
  
  csvPath <- paste(rP,dS$path[i],sep='/')
  fieldName <- with(dS[i,],paste(grower,year,field))
  
  print('Reading in data')
  dat <- read.csv(csvPath,stringsAsFactors=TRUE,fileEncoding='latin1') 
  
  if(nrow(dat)>nSubSamp){
    #Limit to nSubSamp samples
    dat <- dat %>% slice(round(seq(1,nrow(dat),length.out=nSubSamp)))
  }  
  
  #NOTE: NEED TO FIGURE OUT HOW TO PROCESS YIELD DATA. PROBABLY SHOULD JUST LEAVE IT IN POINT FORM FOR NOW RATHER THAN USING DATA-HEAVY POLYGONS
  # COULD PROBABLY ADJUST POINT TO BE IN THE CENTROID, BUT THAT'S ABOUT IT
  
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
    # select(Lon:Distance,TrackAngle,Pass.Num,DryYield,Date) %>% 
    st_as_sf(coords=c('Lon','Lat')) %>% #Add spatial feature info
    st_set_crs(4326) %>% st_transform(3401) %>% #Lat-lon -> UTM
    makePolys(width='SwthWdth_m',dist='Distance_m',angle='TrackAngle') %>%  
    mutate(pArea=as.numeric(st_area(.))) %>% 
    mutate(YieldMass=convertYield(DryYield,'tpha','gpm2')*pArea) %>% 
    mergePoly(fList=lst(Date= first, ID = first, Pass= first, Speed= mean, YieldMass = sum)) %>% #Merges completely overlapping polygons. Takes a few seconds
    mutate(pArea=as.numeric(st_area(.))) %>%
    mutate(DryYield=convertYield(YieldMass/pArea,'gpm2','tpha')) %>%
    mutate(r=1:n()) %>% #row number
    mutate(Pass=factor(seqGroup(ID,FALSE))) %>% 
    group_by(Pass) %>% mutate(rGroup=1:n()) %>% ungroup() %>% 
    mutate(E=st_coordinates(st_centroid(.))[,1],N=st_coordinates(st_centroid(.))[,2]) %>% 
    mutate(E=E-mean(E),N=N-mean(N)) #Center coordinates
  
  print('Calculating field edge')
  fieldEdge <- dat %>% st_union() %>% st_buffer(dist=10) %>% st_cast('MULTILINESTRING')
  
  dat <- dat %>%  #Distance from edge of field
    mutate(dist=as.numeric(st_distance(.,fieldEdge))[1:nrow(.)]) %>% 
    mutate(dist=dist-min(dist)) #Shrink to 0
  
  #Save basic yield map, along with boundary
  print('Making yield map')
  yieldMap <- ggplot(dat)+
    geom_sf(aes(fill=log(DryYield),col=log(DryYield)))+
    geom_sf(data=fieldEdge,col='red')+
    labs(title=fieldName)
  ggsave(paste0('./Figures/YieldMaps/',fieldName,'.png'),yieldMap,height=8,width=8,dpi='screen')  
  
  
  print('Fitting yield model')
  #Fit non-isotropic GAM:
  f <- sqrt(DryYield) ~ s(dist,k=10) + s(E,N,k=60) + s(r,k=30) + log(pArea) #Mean model
  f2 <- ~ s(dist,k=6) + s(E,N,k=60) + s(r,k=60) + log(pArea) #Variance model
  flist <- list(f,f2) #List of model formulae
  
  #Fit Gaussian location-scale  model
  #NOTE: this model can't be run using bam because of location-scale modeling
  a <- Sys.time()
  mod <- gam(flist,data=dat,family=gaulss())
  dS$time[i] <- paste(as.character(round(Sys.time()-a,2)),units(Sys.time()-a)) #Time taken to fit model
  
  print('Saving model results')

  #Plot results
  png(paste0('./Figures/ModelCheck/',fieldName,' Summary.png'),width=8,height=8,units='in',res=200)
  plot(mod,scheme=2,too.far=0.01,pages=1,all.terms=TRUE)
  dev.off()
  
  #GAM check results
  sink(paste0('./Figures/ModelCheck/',fieldName,' results.txt'))
  print("SUMMARY---------------------------------")
  summary(mod)
  
  print("GAM.CHECK-------------------------------")
  
  png(paste0('./Figures/ModelCheck/',fieldName,' gamCheck.png'),width=8,height=8,units='in',res=200)
  par(mfrow=c(2,2)); gam.check(mod); abline(0,1,col='red'); par(mfrow=c(1,1))
  dev.off()
  sink()
  
  #Models are huge, so saving only key components (about 25% of the size)
  modList <- list(coefs=coef(mod),vcv=vcov(mod),smooths=mod$smooth,
                  maxDist=max(dat$dist))
  save(modList,file=paste0('./Figures/ModelCheck/',fieldName,' modList.Rdata'))
  
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

# runModI(1,dS=datSource,rP=rootPath) #Test

library(parallel)
cluster <- makeCluster(6) #Number of cores to use
parLapply(cl=cluster,1:nrow(datSource),runModI,dS=datSource,rP=rootPath)
beep(1)
stopCluster(cluster)
