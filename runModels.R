#Run non-isotropic GAMs on multiple field years

# Load required libraries ---------------------------------------------------------

library(tidyverse)
theme_set(theme_bw())
library(ggpubr)
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
  # i <- 3 #Debugging
  # dS <- datSource
  # rP <- rootPath
  # nSubSamp <- 50000
  
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
  
  #Takes 40 seconds with subsamp of 50000 using full 800000 samples from Alvin French's Al Jr Field
  
  if(nrow(dat)>nSubSamp){
    #Limit to nSubSamp sequential samples
    dat <- dat %>% slice(round(seq(1,nrow(dat),length.out=nSubSamp)))
  }
  
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
  
  ##Old method 42.493 seconds
  # dat <- dat %>% rename_with(.fn = ~gsub('..L.ha.$','_lHa',.x)) %>%
  #   rename_with(.fn = ~gsub('..tonne.ha.$','_tHa',.x)) %>%
  #   rename_with(.fn = ~gsub('..m.s.$','_ms',.x)) %>%
  #   rename_with(.fn = ~gsub('.km.h.$','_kmh',.x)) %>%
  #   rename_with(.fn = ~gsub('..tonne.h.$','_th',.x)) %>%
  #   rename_with(.fn = ~gsub('.ha.h.','_hah',.x)) %>%
  #   rename_with(.fn = ~gsub('.m.','_m',.x)) %>%
  #   rename_with(.fn = ~gsub('.deg.','Angle',.x)) %>%
  #   rename_with(.fn = ~gsub('\\.','',.x)) %>%
  #   rename('ID'='ObjId','DryYield'='YldMassDry_tHa','Lon'='Longitude','Lat'='Latitude','Pass'='PassNum','Speed'='Speed__m') %>%
  #   st_as_sf(coords=c('Lon','Lat')) %>% #Add spatial feature info
  #   st_set_crs(4326) %>% st_transform(3401) %>% #Lat-lon -> UTM
  #   makePolys(width='SwthWdth_m',dist='Distance_m',angle='TrackAngle') %>%  
  #   mutate(pArea=as.numeric(st_area(.))) %>% 
  #   mutate(YieldMass=convertYield(DryYield,'tpha','gpm2')*pArea) %>% 
  #   mergePoly(fList=lst(Date= first, ID = first, Pass= first, Speed= mean, YieldMass = sum)) %>% #Merges completely overlapping polygons. Takes a few seconds
  #   mutate(pArea=as.numeric(st_area(.))) %>%
  #   mutate(DryYield=convertYield(YieldMass/pArea,'gpm2','tpha')) %>%
  #   mutate(r=1:n()) %>% #row number
  #   mutate(Pass=factor(seqGroup(ID,FALSE))) %>% 
  #   group_by(Pass) %>% mutate(rGroup=1:n()) %>% ungroup() %>% 
  #   mutate(E=st_coordinates(st_centroid(.))[,1],N=st_coordinates(st_centroid(.))[,2]) %>% 
  #   mutate(E=E-mean(E),N=N-mean(N)) #Center coordinates
  
  print('Calculating field edge')
  fieldEdge <- dat %>% st_union() %>% st_buffer(dist=10) %>% st_cast('MULTILINESTRING')
  
  dat <- dat %>%  #Distance from edge of field
    mutate(dist=as.numeric(st_distance(.,fieldEdge))[1:nrow(.)]) %>% 
    mutate(dist=dist-min(dist)) #Shrink to 0
  
  print('Fitting yield model')
  #Fit non-isotropic GAM:
  f <- sqrt(DryYield) ~ s(dist,k=10) + s(E,N,k=60) + s(r,k=30) + log(pArea) #Mean model
  f2 <- ~ s(dist,k=6) + s(E,N,k=60) + s(r,k=60) + log(pArea) #Variance model
  flist <- list(f,f2) #List of model formulae
  
  #Fit Gaussian location-scale  model
  #NOTE: this model can't be run using bam because of location-scale modeling
  
  a <- Sys.time() #Takes about 8 mins for a 50000 point model
  mod <- gam(flist,data=dat,family=gaulss())
  fitTime <- paste(as.character(round(Sys.time()-a,2)),units(Sys.time()-a)) #Time taken to fit model
  
  print('Saving model results')

  #Plot results - FAILS HERE WHEN ALL.TERMS=TRUE. SOME KIND OF GAM NAMESPACE PROBLEM THAT ONLY OCCURS WITHIN FUNCTIONS
  #https://stackoverflow.com/questions/45918662/plot-gam-from-mgcv-with-all-terms-true-within-a-function
  png(paste0('./Figures/ModelCheck/',fieldName,' Summary.png'),width=8,height=8,units='in',res=200)
  plot(mod,scheme=2,too.far=0.01,pages=1,all.terms=TRUE)
  dev.off()
  
  #GAM check results
  sink(paste0('./Figures/ModelCheck/',fieldName,' results.txt'))
  print("Time taken: "); print(fitTime)
  print(" ")
  print("SUMMARY---------------------------------")
  summary(mod)
  print(" ")
  print("GAM.CHECK-------------------------------")
  png(paste0('./Figures/ModelCheck/',fieldName,' gamCheck.png'),width=8,height=8,units='in',res=200)
  par(mfrow=c(2,2)); gam.check(mod); abline(0,1,col='red'); par(mfrow=c(1,1))
  sink()
  dev.off()
  
  #Actual, Predicted, and SE maps
  print('Making yield map')
  yieldMap <- ggplot(dat) + geom_sf(aes(col=log(DryYield),alpha=0.3)) + geom_sf(data=fieldEdge,col='black')+
    labs(title=paste(fieldName,'Data'))+guides(alpha='none') + scale_colour_distiller(type='div',palette = "Spectral") +
    theme(legend.position='bottom')
  dat <- dat %>% mutate(yieldPred=mod$fit[,1],yieldVar=mod$fit[,2],resid=resid(mod))
  p1 <- ggplot(dat)+geom_sf(aes(col=yieldPred),alpha=0.5)+labs(title='Predicted yield',fill='(sqrt)\nYield',col='(sqrt)\nYield')+scale_colour_distiller(type='div',palette = "Spectral")+theme(legend.position='bottom')
  p2 <- ggplot(dat)+geom_sf(aes(col=yieldVar),alpha=0.5)+labs(title='Yield SD',fill='SD\nYield',col='SD\nYield')+scale_colour_distiller(type='div',palette = "Spectral")+theme(legend.position='bottom')
  p <- ggarrange(yieldMap,p1,p2,ncol=3)
  ggsave(paste0('./Figures/YieldMaps/',fieldName,'.png'),p,height=6,width=15,dpi=100)  
  
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
  p5 <- gstat::variogram(resid~1, data=dat_sp,width=5,cutoff=150) %>% 
    ggplot()+geom_point(aes(x=dist,y=gamma))+
    geom_line(aes(x=dist,y=gamma))+
    labs(title='Residual variogram',x='Distance',y='Semivariance')
  # dat_variogram2 <- gstat::variogram(resid~1, #Used to look at non-stationary semivariance
  #                                   data=dat_sp,width=5,cutoff=150,alpha=c(0,45,90,135)) #North, NE, E, SW
  res_acf <- acf(dat$resid,lag.max=500,type='correlation')
  p6 <- with(res_acf,data.frame(acf=acf,lag=lag)) %>% 
    ggplot(aes(x=lag,y=acf))+geom_col()+
    labs(x='Time Lag',y='Autocorrelation',title='Residual autocorrelation')
  p <- ggarrange(p3,p4,p5,p6,ncol=2,nrow=2)
  p <- annotate_figure(p,top = text_grob(paste0(fieldName,' Residuals')))
  ggsave(paste0('./Figures/ModelCheck/',fieldName,' residuals.png'),p,height=8,width=10,dpi=100)  
  
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

runModI(1,dS=datSource,rP=rootPath) #Test

library(parallel)
cluster <- makeCluster(6) #Number of cores to use
parLapply(cl=cluster,1:nrow(datSource),runModI,dS=datSource,rP=rootPath)
beep(1)
stopCluster(cluster)
