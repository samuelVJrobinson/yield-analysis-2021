#Seems to be a consensus that filter/cleaning is necessary for yield data. Trying to see how much of a difference this makes for our models

library(tidyverse)
theme_set(theme_classic())
library(sf)
library(mgcv)
library(beepr)
library(ggpubr)
setwd("~/Documents/yield-analysis-2021/gam_TMB")
source('../helperFunctions.R')

nSubSamp <- 50000 #Number of points to start with
csvPath <- "/media/rsamuel/Storage/geoData/Rasters/yieldData/csv files/Alvin_French C41 2020.csv"
boundaryPath <- "/home/rsamuel/Documents/yield-analysis-2021/Figures/FieldBoundaries/Alvin_French C41 2020 boundary.shp"
boundaryPath2 <- "/home/rsamuel/Documents/yield-analysis-2021/Figures/FieldBoundarySegments/Alvin_French C41 2020 boundary.shp"

#Also try: Gibbons Bill Visser 2017 (bad residual structure)

dat <- read.csv(csvPath,stringsAsFactors=TRUE,fileEncoding='latin1') 

if(nrow(dat)>nSubSamp){
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
  mutate(r=1:n()) %>% #row number
  
  #Remove extreme values
  filter(pArea>quantile(pArea,0.05),pArea<quantile(pArea,0.95), #Filter large/small pArea
         DryYield>quantile(DryYield,0.05),DryYield<quantile(DryYield,0.95), #Filter extreme yields
         Speed>quantile(Speed,0.05),Speed<quantile(Speed,0.95), #Filter high and low speeds
         ) %>%
  
  ## Keep extreme values
  # filter(!(pArea>quantile(pArea,0.025)&pArea<quantile(pArea,0.975)& #Filter large/small pArea
  #        DryYield>quantile(DryYield,0.025)&DryYield<quantile(DryYield,0.975)& #Filter extreme yields
  #        Speed>quantile(Speed,0.025)&Speed<quantile(Speed,0.975)), #Filter high and low speeds
  # ) %>% 
  st_centroid() %>% #Convert back to point
  mutate(Pass=factor(seqGroup(ID,FALSE))) %>% 
  group_by(Pass) %>% mutate(rGroup=1:n()) %>% ungroup() %>% 
  mutate(E=st_coordinates(.)[,1],N=st_coordinates(.)[,2]) %>% 
  mutate(E=(E-mean(E)),N=(N-mean(N))) #Center coordinates
# ggplot(dat)+geom_sf()

fieldEdge <- read_sf(boundaryPath) #Read in boundary polygon
fieldEdgeType <- read_sf(boundaryPath2) #Read in boundary type linestrings

#Get distance and type of closest boundary
dat <- dat %>% 
  bind_cols(.,data.frame(t(apply(st_distance(dat,fieldEdgeType),1,function(x) c(dist=min(x,na.rm=TRUE),boundaryType=fieldEdgeType$type[which.min(x)]))))) %>% 
  mutate(dist=as.numeric(dist),boundaryType=factor(boundaryType))

ggplot(dat)+geom_sf(alpha=0.3)+geom_sf(data=fieldEdge,col='red',fill=NA)

#Quantiles
qs <- c(0.01,0.025,0.05,0.10,0.25,0.5,0.75,0.9,0.95,0.975,0.99)
quantile(dat$DryYield,qs); hist(dat$DryYield)
quantile(dat$pArea,qs); hist(dat$pArea)
quantile(dat$Speed,qs); hist(dat$Speed)

#Calculate limits
lyield <- log(dat$DryYield) #Log-transform
lyieldMean <- mean(lyield)
lyieldSD <- sd(lyield)
lyieldLims <- lyieldMean+lyieldSD*c(-1,1)*3 #3 x SD limits
lyieldLims2 <- unname(quantile(lyield,c(0.025,0.975))) #1st and 99th percentiles
windowQuantile <- function(i,x,quant=c(0.025,0.975),win=101){
  if(win%%2==0) stop('Window must be an odd number')
  iLwr <- i-(win-1)/2; iUpr <- i+(win-1)/2
  if(iLwr<0) iLwr <- 0
  if(iUpr>length(x)) iUpr <- length(x)
  return(unname(quantile(x[iLwr:iUpr],probs=quant)))
}
lyieldLims3 <- t(sapply(1:length(lyield),windowQuantile,x=lyield,quant=c(0.01,0.99),win=201))


ylimits <- range(lyield)
par(mfrow=c(2,3))
plot(lyield,type='l',main='3 x SD limits',ylim=ylimits)
abline(h=lyieldLims[1],col='red')
abline(h=lyieldLims[2],col='red')

plot(lyield,type='l',main='98% quantile limits',ylim=ylimits)
abline(h=lyieldLims2[1],col='red')
abline(h=lyieldLims2[2],col='red')

plot(lyield,type='l',main='Windowed filter',ylim=ylimits)
lines(lyieldLims3[,1],col='red')
lines(lyieldLims3[,2],col='red')

#Show data removed
plot(lyield[lyield>lyieldLims[1] & lyield<lyieldLims[2]],type='l',main='3 x SD limits',ylim=ylimits)
plot(lyield[lyield>lyieldLims2[1] & lyield<lyieldLims2[2]],type='l',main='95% quantile limits',ylim=ylimits)
plot(lyield[lyield>lyieldLims3[,1] & lyield<lyieldLims3[,2]],type='l',main='Windowed filter',ylim=ylimits)
sum(!(lyield>lyieldLims3[,1] & lyield<lyieldLims3[,2]))/length(lyield) #Proportion of data removed

#QQnorm plots
par(mfrow=c(2,2))
qqnorm(lyield,main='Raw'); qqline(lyield)
qqnorm(lyield[lyield>lyieldLims[1] & lyield<lyieldLims[2]],main='3xSD'); qqline(lyield[lyield>lyieldLims[1] & lyield<lyieldLims[2]])
qqnorm(lyield[lyield>lyieldLims2[1] & lyield<lyieldLims2[2]],main='95% quantiles'); qqline(lyield[lyield>lyieldLims2[1] & lyield<lyieldLims2[2]])
qqnorm(lyield[lyield>lyieldLims[1] & lyield<lyieldLims[2]],main='Windowed filter'); qqline(lyield[lyield>lyieldLims[1] & lyield<lyieldLims[2]])


#Fit non-isotropic GAM:
kPar <- c(12,60,60,12,60,60)
f <- sqrt(DryYield) ~ s(dist,k=kPar[1]) + s(E,N,k=kPar[2]) + s(r,k=kPar[3]) + log(pArea) #Mean model
f2 <- ~ s(dist,k=kPar[4]) + s(E,N,k=kPar[5]) + s(r,k=kPar[6]) + log(pArea) #Variance model
flist <- list(f,f2) #List of model formulae

a <- Sys.time() 
mod2 <- gam(flist,data=dat,family=gaulss())
fitTime <- paste(as.character(round(Sys.time()-a,2)),units(Sys.time()-a)) #Time taken to fit model
# save(mod,file =  "/home/rsamuel/Documents/yield-analysis-2021/Data/unfilteredMod.Rdata")
# save(mod2,file =  "/home/rsamuel/Documents/yield-analysis-2021/Data/filteredMod.Rdata")

load("/home/rsamuel/Documents/yield-analysis-2021/Data/unfilteredMod.Rdata")
load("/home/rsamuel/Documents/yield-analysis-2021/Data/filteredMod.Rdata")

plot(mod,pages=1,scheme=2)
plot(mod2,pages=1,scheme=2)

par(mfrow=c(2,4))
gam.check(mod); abline(0,1,col='red')
gam.check(mod2); abline(0,1,col='red')

pred <- dat %>% mutate(pArea=median(pArea)) %>% #Set pArea to its median value
  predict(mod,newdata = .,type='response',exclude='s(r)') #Exclude s(r)
pred2 <- dat %>% mutate(pArea=median(pArea)) %>% #Set pArea to its median value
  predict(mod2,newdata = .,type='response',exclude='s(r)') #Exclude s(r)

p1 <- dat %>% mutate(yieldMean=pred[,1],yieldSD=1/pred[,2]) %>% 
  ggplot()+geom_sf(aes(col=yieldMean^2),alpha=0.3)+
  labs(title='Mod1: Predicted yield | Combine Speed, Point Order',fill='Yield (T per Ha)',col='Yield (T per Ha)')+
  scale_colour_distiller(type='div',palette = "Spectral",direction=-1)+theme(legend.position='bottom')

p2 <- dat %>% mutate(yieldMean=pred[,1],yieldSD=1/pred[,2]) %>% 
  ggplot()+geom_sf(aes(col=yieldSD^2),alpha=0.3)+
  labs(title='Mod1: Yield SD | Combine Speed, Point Order',fill='SD Yield ',col='SD Yield')+
  scale_colour_distiller(type='div',palette = "Spectral",direction=-1)+theme(legend.position='bottom')

p3 <- dat %>% mutate(yieldMean=pred2[,1]) %>% 
  ggplot()+geom_sf(aes(col=yieldMean^2),alpha=0.3)+
  labs(title='Mod2: Predicted yield | Combine Speed, Point Order',fill='Yield (T per Ha)',col='Yield (T per Ha)')+
  scale_colour_distiller(type='div',palette = "Spectral",direction=-1)+theme(legend.position='bottom')

p4 <- dat %>% mutate(yieldSD=1/pred2[,2]) %>% 
  ggplot()+geom_sf(aes(col=yieldSD^2),alpha=0.3)+
  labs(title='Mod2: Yield SD | Combine Speed, Point Order',fill='SD Yield ',col='SD Yield')+
  scale_colour_distiller(type='div',palette = "Spectral",direction=-1)+theme(legend.position='bottom')

ggarrange(p1,p2,p3,p4,ncol=2,nrow=2)

#Still heavy space-time dependence in residuals

pred2 <- dat %>% mutate(pArea=median(pArea)) %>% #Set pArea to its median value
  predict(mod2,newdata = .,type='response') 

dat %>% mutate(p=pred2[,1],actual=p+resid(mod2)) %>% #Residuals across time
  mutate(across(c(p,actual),~.^2)) %>% 
  filter(r<2500) %>% 
  ggplot(aes(x=r))+
  geom_point(aes(y=actual),col='red')+
  geom_point(aes(y=p))

library(gstat)
library(spacetime) 
dat_sp <- dat %>% mutate(resid=resid(mod2)) %>% select(ID,resid) %>% 
  as_Spatial()
v <- variogram(resid~1, data=dat_sp,width=5,cutoff=300) #Variogram
p5 <- v %>% 
  ggplot()+geom_point(aes(x=dist,y=gamma))+
  labs(title='Residual variogram',x='Distance',y='Semivariance')
res_acf <- acf(dat_sp$resid,lag.max=300,type='correlation',plot=FALSE)
p6 <- with(res_acf,data.frame(acf=acf,lag=lag)) %>% 
  ggplot(aes(x=lag,y=acf))+geom_col()+
  labs(x='Time Lag',y='Autocorrelation',title='Residual autocorrelation')

  
