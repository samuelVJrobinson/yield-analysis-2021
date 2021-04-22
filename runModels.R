#Run non-isotropic GAMs on multiple field years

# Load required libraries ---------------------------------------------------------

library(tidyverse)
theme_set(theme_bw())
library(ggpubr)
library(mgcv)
library(sf)
library(beepr)
source('helperFunctions.R')

# Get paths to data -------------------------------------------------------

# #Generate new dataSource dataframe
# rootPath <- "/media/rsamuel/Storage/geoData/Rasters/yieldData/csv files" #Path to csv files
# datSource <- data.frame(dataPath=dir(rootPath,pattern=".csv",recursive=TRUE,full.names = TRUE)) %>% 
#   mutate(path=gsub('/media/rsamuel/Storage/geoData/Rasters/yieldData/csv files/','',dataPath)) %>% 
#   separate(path,c('grower','year','field'),sep="/",remove=FALSE) %>% select(-path) %>% 
#   mutate(field=gsub('\\.csv','',field)) %>% 
#   unite(filename,c(grower:field),sep=' ',remove=FALSE) %>% 
#   mutate(boundaryPath=paste0('./Figures/FieldBoundaries/',filename,' boundary.shp')) %>% #Paths to boundary shapefiles
#   mutate(modelPath=paste0('./Figures/ModelCheck/',filename,' modList.Rdata'))  #Paths to saved model files
# write.csv(datSource,'./Data/datSource_new.csv',row.names = FALSE)

datSource <- read.csv('./Data/datSource.csv') #Read previous datasource file
datSource <- datSource %>% mutate(boundaryComplete=file.exists(boundaryPath)) %>% #Has boundary been made already?
  mutate(modelComplete=file.exists(modelPath)) #Has model already been run?

datSource %>% count(grower,year) 

# set.seed(1) # Run subset of models
# datSource <- datSource %>%
#   # filter(grower!='Alvin French') %>% #These fields are huge, so leave out for now
#   slice_sample(n=10) #Smaller set to experiment with

# Run models --------------------------------------------------------------

#Function to run the ith model
runModI <- function(i,dS,nSubSamp=50000){ 
  # i <- 1 #Debugging
  # dS <- datSource
  # nSubSamp <- 50000

  if(dS$modelComplete[i]) return('Already completed')
  if(!dS$use[i]) return('Not used')
  
  require(tidyverse)
  theme_set(theme_bw())
  require(mgcv)
  require(sf)
  require(ggpubr)
  
  source('helperFunctions.R')
  
  csvPath <- dS$dataPath[i]
  fieldName <- dS$filename[i]
  boundaryPath <- dS$boundaryPath[i]
  
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
  
  print('Calculating distance from boundary')  
  
  fieldEdge <- read_sf(boundaryPath) %>% st_cast('MULTILINESTRING') #Read in boundary file
  # fieldEdge <- dat %>% st_union() %>% st_buffer(dist=10) %>% st_cast('MULTILINESTRING') #Old method
  
  dat <- dat %>%  #Distance from boundary
    mutate(dist=as.numeric(st_distance(.,fieldEdge))[1:nrow(.)]) 
  
  print('Fitting model')
  #Fit non-isotropic GAM:
  f <- sqrt(DryYield) ~ s(dist,k=12) + s(E,N,k=60) + s(r,k=60) + log(pArea) #Mean model
  f2 <- ~ s(dist,k=12) + s(E,N,k=60) + s(r,k=60) + log(pArea) #Variance model
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
  png(paste0('./Figures/ModelCheck/',fieldName,' Summary.png'),width=12,height=8,units='in',res=200)
  # plot(mod,scheme=2,too.far=0.01,pages=1,all.terms=TRUE)
  plot(mod,scheme=2,too.far=0.01,pages=1,all.terms=FALSE)
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
  ggsave(paste0('./Figures/ModelCheck/',fieldName,' residuals.png'),p,height=8,width=10,dpi=100)  
  
  #Models are huge, so saving only key components (about 25% of the size)
  modList <- list(coefs=coef(mod), #Coefficients
                  vcv=vcov(mod), #Covariance matrix
                  smooths=mod$smooth, #Smooth specifications
                  distRange=range(dat$dist), #Range of distances
                  pAreaRange=range(dat$pArea), #Range of polygon areas
                  rRange=range(dat$r)) #Range of r
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
# 
# a <- Sys.time()
# runModI(2,dS=datSource) #Test
# Sys.time()-a
# beep(1)

library(parallel)
cluster <- makeCluster(8) #10 procs max - uses about 90% of memory
parLapply(cl=cluster,1:nrow(datSource),runModI,dS=datSource) #Takes about 42 mins for running 10 procs. Some seem to take longer than others (weird shaped fields?)
beep(1)
stopCluster(cluster)
Sys.time() #Takes ~ 7 hrs

#Need to identify boundary types at each field

# Get smoother info from models -------------------------------------------

# #Estimate at field 1
# est <- getSmooths(paste0('./Figures/ModelCheck/',datSource$filename[1],' modList.Rdata'),margInt=c(FALSE,TRUE,TRUE),samp=FALSE)$distDat
# #Samples around estimate at field 1
# Nsamp <- 100 #Number of samples
# samp <- lapply(1:Nsamp,
#           function(x) getSmooths(paste0('./Figures/ModelCheck/',datSource$filename[1],' modList.Rdata'),
#                                  margInt=c(FALSE,TRUE,TRUE),samp=TRUE)$distDat) %>% 
#   do.call('rbind',.) %>% mutate(N=rep(1:length(unique(dist)),each=Nsamp))
# ggplot()+geom_line(data=test,aes(x=dist,y=mean,group=N),alpha=0.3)+geom_line(data=est,aes(x=dist,y=mean),col='red',size=3)

#Get smoother from each field
#Takes about 10 seconds
getFiles <- datSource$filename[datSource$use]
allSmooths <- lapply(paste0('./Figures/ModelCheck/',getFiles,' modList.Rdata'),
                     getSmooths,margInt=c(FALSE,TRUE,TRUE))
names(allSmooths) <- gsub(' ','-',getFiles)
names(allSmooths[[1]]) #Variables to get from allSmooths

#Get df of predictions for each variable
allEff <- lapply(names(allSmooths[[1]]),function(y){
  lapply(allSmooths,function(x) x[[y]]) %>% 
    do.call('rbind',.) %>% 
    rownames_to_column('field') %>% 
    mutate(field=gsub('\\.\\d{1,3}','',field))
})

#Plot with individual field-level smoothers

tHa2buAc <- 17.0340 #tonnes per hectare to bushels per acre (https://www.agrimoney.com/calculators/calculators)
ylimVals <- list(mean=c(-1.2,1.2)*tHa2buAc,sd=c(0,5)) #Y limits
ylabMean <- 'Mean Yield (bushels/acre)'
ylabSD <- 'log(SD Yield)'

p1 <- allEff[[1]] %>% ggplot(aes(x=pArea,y=mean*tHa2buAc))+geom_line(aes(group=field),alpha=0.3)+
  geom_smooth(method='lm',formula=y~log2(x),col='blue',se=FALSE,n=500)+
  labs(x='Polygon Area',y=ylabMean)+
  coord_cartesian(xlim = c(0,200))
p2 <- allEff[[1]] %>% ggplot(aes(x=pArea,y=logSD*tHa2buAc))+geom_line(aes(group=field),alpha=0.3)+
  geom_smooth(method='lm',formula=y~log(x),col='blue',se=FALSE,n=500)+
  labs(x='Polygon Area',y=ylabSD)+coord_cartesian(xlim = c(0,200))

p3 <- allEff[[2]] %>% ggplot(aes(x=dist,y=mean*tHa2buAc))+geom_line(aes(group=field),alpha=0.3)+
  labs(x='Boundary Distance',y=ylabMean)+
  geom_hline(yintercept = 0,linetype='dashed',col='red')+
  geom_smooth(method='gam',formula=y~s(x,k=10),col='blue',se=TRUE)+
  coord_cartesian(xlim = c(0,400),ylim=ylimVals$mean)
p4 <- allEff[[2]] %>% ggplot(aes(x=dist,y=log(exp(logSD)*tHa2buAc)))+geom_line(aes(group=field),alpha=0.3)+
  labs(x='Boundary Distance',y=ylabSD)+
  # geom_hline(yintercept = 0,linetype='dashed',col='red')+
  geom_smooth(method='gam',formula=y~s(x),col='blue',se=FALSE)+
  coord_cartesian(xlim = c(0,400),ylim=ylimVals$sd)

p5 <- allEff[[3]] %>% ggplot(aes(x=r,y=mean*tHa2buAc))+geom_line(aes(group=field),alpha=0.3)+
  labs(x='Point Order',y=ylabMean)+
  # geom_hline(yintercept = 0,linetype='dashed',col='red')+
  geom_smooth(method='gam',formula=y~s(x),col='blue',se=FALSE)+
  coord_cartesian(ylim=ylimVals$mean)
p6 <- allEff[[3]] %>% ggplot(aes(x=r,y=log(exp(logSD)*tHa2buAc)))+geom_line(aes(group=field),alpha=0.3)+
  labs(x='Point Order',y=ylabSD)+
  # geom_hline(yintercept = 0,linetype='dashed',col='red')+
  geom_smooth(method='gam',formula=y~s(x),col='blue',se=FALSE)+
  coord_cartesian(ylim=ylimVals$sd)

(p <- ggarrange(p1,p3,p5,p2,p4,p6,ncol=3,nrow=2))
ggsave(paste0('./Figures/ModelSummary.png'),p,height=6,width=12,dpi=300)  

#Get smoother from each field
#Takes about 10 seconds

getFiles <- datSource$filename[datSource$use] #Field names
files <- paste0('./Figures/ModelCheck/',getFiles,' modList.Rdata') #File paths

#Function to get predicted smoother values from list of functions, fit models to these, and return results
# Used for getting 
# f = files to draw from
# s = sample from posterior?
sampleSmooth <- function(f,s,...){
  require(mgcv)
  require(tidyverse)
  
  allSmooths <- lapply(f,getSmooths,margInt=c(FALSE,TRUE,TRUE),samp=s)
  names(allSmooths) <- gsub(' ','-',f)
  names(allSmooths[[1]]) #Variables to get from allSmooths
  
  #Get df of predictions for each variable
  allEff <- lapply(names(allSmooths[[1]]),function(y){
    lapply(allSmooths,function(x) x[[y]]) %>% 
      do.call('rbind',.) %>% 
      rownames_to_column('field') %>% 
      mutate(field=gsub('\\.\\d{1,3}','',field))
  })
  
  #Fit models of sampled estimates
  m1mean <- lm(mean~log(pArea),data=allEff[[1]])
  m1sd <- lm(logSD~log(pArea),data=allEff[[1]])
  m2mean <- gam(mean~s(dist,k=30),data=allEff[[2]])
  m2sd <- gam(logSD~s(dist,k=30),data=allEff[[2]])
  m3mean <- gam(mean~s(r,k=30),data=allEff[[3]])
  m3sd <- gam(logSD~s(r,k=30),data=allEff[[3]])
  
  #Predictions from models
  datList <- list(
    #200 locations along log-pArea line
    pArea=data.frame(pArea=exp(seq(min(log(allEff[[1]]$pArea)),max(log(allEff[[1]]$pArea)),length.out=200))) %>%
      mutate(mean=predict(m1mean,newdata=.),logSD=predict(m1sd,newdata=.)),
    #Unique rounded distance measurements
    dist=data.frame(dist=sort(unique(round(allEff[[2]]$dist)))) %>% 
      mutate(mean=predict(m2mean,newdata=.),logSD=predict(m2sd,newdata=.)), 
    #Unique rounded sequence measurements
    r=data.frame(r=round(seq(min(allEff[[3]]$r),max(allEff[[3]]$r),length.out=1000))) %>% 
      mutate(mean=predict(m3mean,newdata=.),logSD=predict(m3sd,newdata=.))
  )
  rm(list=ls()[ls()!='datList']) #Cleanup
  gc()
  return(datList)
}

temp <- sampleSmooth(f=files,s=FALSE) #Mean smoother

Nrep <- 3
# temp2 <- replicate(Nrep,sampleSmooth(f=files,s=TRUE),simplify=FALSE)
debugonce(sampleSmooth)
debugonce(getSmooths)
temp2 <- sampleSmooth(f=files,s=TRUE)


lapply(temp2,function(x) x$r) %>% set_names(paste0('s',1:Nrep)) %>% 
  do.call('rbind',.) %>%  rownames_to_column('rep') %>%   mutate(rep=gsub('\\.\\d{1,3}','',rep)) %>% 
  group_by(rep) %>% slice(1:3)

library(parallel)
Nrep <- 30
cluster <- makeCluster(8) 
clusterExport(cl = cluster,varlist='getSmooths', envir = .GlobalEnv) #Export function to clusters
tempSamp <- parLapply(cl=cluster,1:Nrep,fun=sampleSmooth,f=files,s=TRUE)
beep(1)
stopCluster(cluster)

#Get range of variability for models 

#Low variability in pArea models
p1 <- lapply(tempSamp,function(x) x$pArea) %>% set_names(paste0('s',1:Nrep)) %>% 
  do.call('rbind',.) %>%  rownames_to_column('rep') %>%   mutate(rep=gsub('\\.\\d{1,3}','',rep)) %>% 
  ggplot(aes(x=pArea,y=mean))+
  geom_line(aes(group=rep),alpha=0.3)+
  geom_line(data=temp$pArea,col='red')
p2 <- lapply(tempSamp,function(x) x$pArea) %>% set_names(paste0('s',1:Nrep)) %>% 
  do.call('rbind',.) %>%  rownames_to_column('rep') %>%   mutate(rep=gsub('\\.\\d{1,3}','',rep)) %>% 
  ggplot(aes(x=pArea,y=logSD))+
  geom_line(aes(group=rep),alpha=0.3)+
  geom_line(data=temp$pArea,col='red')

#Higher variability in dist model
p3 <- lapply(tempSamp,function(x) x$dist) %>% set_names(paste0('s',1:Nrep)) %>% 
  do.call('rbind',.) %>%  rownames_to_column('rep') %>%   mutate(rep=gsub('\\.\\d{1,3}','',rep)) %>% 
  ggplot(aes(x=dist,y=mean))+
  geom_line(aes(group=rep),alpha=0.3)+
  geom_line(data=temp$dist,col='red',size=2)
p4 <- lapply(tempSamp,function(x) x$dist) %>% set_names(paste0('s',1:Nrep)) %>% 
  do.call('rbind',.) %>%  rownames_to_column('rep') %>%   mutate(rep=gsub('\\.\\d{1,3}','',rep)) %>% 
  ggplot(aes(x=dist,y=logSD))+
  geom_line(aes(group=rep),alpha=0.3)+
  geom_line(data=temp$dist,col='red',size=2)

#Low variability in r model
p5 <- lapply(tempSamp,function(x) x$r) %>% set_names(paste0('s',1:Nrep)) %>% 
  do.call('rbind',.) %>%  rownames_to_column('rep') %>%   mutate(rep=gsub('\\.\\d{1,3}','',rep)) %>% 
  ggplot(aes(x=r,y=mean))+
  geom_line(aes(group=rep),alpha=0.3)+
  geom_line(data=temp$r,col='red',size=2)
p6 <- lapply(tempSamp,function(x) x$r) %>% set_names(paste0('s',1:Nrep)) %>% 
  do.call('rbind',.) %>%  rownames_to_column('rep') %>%   mutate(rep=gsub('\\.\\d{1,3}','',rep)) %>% 
  ggplot(aes(x=r,y=logSD))+
  geom_line(aes(group=rep),alpha=0.3)+
  geom_line(data=temp$r,col='red',size=2)
