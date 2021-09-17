#Seems to be a consensus that filter/cleaning is necessary for yield data. Trying to see how much of a difference this makes for our models

# Load everything -----------------------------------------------

library(tidyverse)
theme_set(theme_classic())
library(sf)
library(mgcv)
library(beepr)
library(ggpubr)
setwd("~/Documents/yield-analysis-2021")
source('./helperFunctions.R')

datSource <- read.csv('./Data/datSource.csv') %>% #Read previous datasource file
  mutate(boundaryComplete=file.exists(boundaryPath)) %>% #Has boundary been made already?
  mutate(modelComplete1=file.exists(modelPath1)) %>% #Has model1 already been run?
  mutate(modelComplete2=file.exists(modelPath2)) %>% #Has model2 already been run?
  mutate(modelComplete0=file.exists(modelPath0)) #Has model0 already been run?

nSubSamp <- 50000 #Number of points to start with

fieldName <- 'Alvin_French C41 2020' #Canola - extremely messy - 1st filter seems better, but hard to tell
# fieldName <- 'Gibbons Bill Visser 2017' #Peas - 2nd filter seems better
# fieldName <- 'Trent_Clark E 06 2018' #Canola - chunks missing - 2nd filter seems better
# fieldName <- 'Dean_Hubbard E_21_11_25 2018' #Wheat
# fieldName <- 'Dean_Hubbard E_21_11_25 2020' #Wheat

rownum <- which(datSource$filename==fieldName)
csvPath <- datSource$dataPath[rownum] #Get file paths
boundaryPath <- datSource$boundaryPath[rownum]
boundaryPath2 <- datSource$boundaryPath2[rownum]
(cropType <- datSource$crop[rownum])

dat <- read.csv(csvPath,stringsAsFactors=TRUE,fileEncoding='latin1') 

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
  st_centroid() %>% #Convert back to point
  mutate(Pass=factor(seqGroup(ID,FALSE))) %>% 
  group_by(Pass) %>% mutate(rGroup=1:n()) %>% ungroup() %>% 
  mutate(E=st_coordinates(.)[,1],N=st_coordinates(.)[,2]) %>% 
  mutate(E=(E-mean(E)),N=(N-mean(N))) #Center coordinates

if(nrow(dat)>nSubSamp){ #If too many samples
  dat <- dat %>% slice(round(seq(1,nrow(dat),length.out=nSubSamp)))
}

fieldEdge <- read_sf(boundaryPath) #Read in boundary polygon
fieldEdgeType <- read_sf(boundaryPath2) #Read in boundary type linestrings

#Get distance and type of closest boundary
dat <- dat %>% 
  bind_cols(.,data.frame(t(apply(st_distance(dat,fieldEdgeType),1,function(x) c(dist=min(x,na.rm=TRUE),boundaryType=fieldEdgeType$type[which.min(x)]))))) %>% 
  mutate(dist=as.numeric(dist),boundaryType=factor(boundaryType))

# ggplot(dat)+geom_sf(aes(col=DryYield),alpha=0.3)+
#   geom_sf(data=fieldEdge,col='red',fill=NA)+
#   scale_colour_distiller(type='div',palette = "Spectral",direction=1)

dat2 <- dat %>% #Trying out different filters
  #Filter 1: Retains only 90th percentils of pArea, dryYield, and speed
  mutate(filt=pArea<quantile(pArea,0.05)|pArea>quantile(pArea,0.95)| #Large/small pArea
           DryYield<quantile(DryYield,0.05)|DryYield>quantile(DryYield,0.95)| #Filter extreme yields
           Speed<quantile(Speed,0.05)|Speed>quantile(Speed,0.95) #Filter high and low speeds
         )
dat3 <- dat %>% 
  #Filter 2: Uses angle and yield differences
  mutate(AngleDiff=abs(bearingDiff(lag(TrackAngle),TrackAngle))) %>% #Absolute value of bearing angle difference (turning)
  mutate(YieldDiff1=abs((lag(DryYield)-DryYield))) %>% #Lagged difference in dry yield
  filter(!is.na(YieldDiff1)) %>%
  mutate(filt=YieldDiff1>quantile(YieldDiff1,0.97)|AngleDiff>quantile(AngleDiff,0.97)|DryYield>quantile(DryYield,0.99)|DryYield<quantile(DryYield,0.01))

p1 <- dat2 %>% #Lineplot
  filter(r>23000&r<24000) %>%
  ggplot(aes(x=r,y=DryYield))+
  geom_line()+
  geom_point(aes(col=filt))
p2 <- dat3 %>% #Lineplot
  filter(r>23000&r<24000) %>%
  ggplot(aes(x=r,y=DryYield))+
  geom_line()+
  geom_point(aes(col=filt))
ggarrange(p1,p2,ncol=1,nrow=2)


p1 <- dat2 %>% filter(!filt) %>% #Yield data after filtering
  ggplot()+
  geom_sf(aes(col=DryYield),alpha=0.3)+
  geom_sf(data=fieldEdge,col='red',fill=NA)+
  scale_colour_distiller(type='div',palette = "Spectral",direction=1)
p2 <- dat3 %>% filter(!filt) %>% #Yield data after filtering
  ggplot()+
  geom_sf(aes(col=DryYield),alpha=0.3)+
  geom_sf(data=fieldEdge,col='red',fill=NA)+
  scale_colour_distiller(type='div',palette = "Spectral",direction=1)
ggarrange(p1,p2,ncol=2,nrow=1)

p1 <- dat2 %>% 
  ggplot()+ #Where are filtered data located in the field?
  geom_sf(aes(col=filt),alpha=0.3)+
  geom_sf(data=fieldEdge,col='red',fill=NA)+
  scale_color_manual(values=c('yellow','black'))
p2 <- dat3 %>% 
  ggplot()+ #Where are filtered data located in the field?
  geom_sf(aes(col=filt),alpha=0.3)+
  geom_sf(data=fieldEdge,col='red',fill=NA)+
  scale_color_manual(values=c('yellow','black'))
ggarrange(p1,p2,ncol=2,nrow=1)

# Look at data --------------------------------------

# #Quantiles
# qs <- c(0.01,0.025,0.05,0.10,0.25,0.5,0.75,0.9,0.95,0.975,0.99)
# quantile(dat$DryYield,qs); hist(dat$DryYield)
# quantile(dat$pArea,qs); hist(dat$pArea)
# quantile(dat$Speed,qs); hist(dat$Speed)
# 
# #Calculate limits
# lyield <- log(dat$DryYield) #Log-transform
# lyieldMean <- mean(lyield)
# lyieldSD <- sd(lyield)
# lyieldLims <- lyieldMean+lyieldSD*c(-1,1)*3 #3 x SD limits
# lyieldLims2 <- unname(quantile(lyield,c(0.025,0.975))) #1st and 99th percentiles
# windowQuantile <- function(i,x,quant=c(0.025,0.975),win=101){
#   if(win%%2==0) stop('Window must be an odd number')
#   iLwr <- i-(win-1)/2; iUpr <- i+(win-1)/2
#   if(iLwr<0) iLwr <- 0
#   if(iUpr>length(x)) iUpr <- length(x)
#   return(unname(quantile(x[iLwr:iUpr],probs=quant)))
# }
# lyieldLims3 <- t(sapply(1:length(lyield),windowQuantile,x=lyield,quant=c(0.01,0.99),win=201))
# 
# 
# ylimits <- range(lyield)
# par(mfrow=c(2,3))
# plot(lyield,type='l',main='3 x SD limits',ylim=ylimits)
# abline(h=lyieldLims[1],col='red')
# abline(h=lyieldLims[2],col='red')
# 
# plot(lyield,type='l',main='98% quantile limits',ylim=ylimits)
# abline(h=lyieldLims2[1],col='red')
# abline(h=lyieldLims2[2],col='red')
# 
# plot(lyield,type='l',main='Windowed filter',ylim=ylimits)
# lines(lyieldLims3[,1],col='red')
# lines(lyieldLims3[,2],col='red')
# 
# #Show data removed
# plot(lyield[lyield>lyieldLims[1] & lyield<lyieldLims[2]],type='l',main='3 x SD limits',ylim=ylimits)
# plot(lyield[lyield>lyieldLims2[1] & lyield<lyieldLims2[2]],type='l',main='95% quantile limits',ylim=ylimits)
# plot(lyield[lyield>lyieldLims3[,1] & lyield<lyieldLims3[,2]],type='l',main='Windowed filter',ylim=ylimits)
# sum(!(lyield>lyieldLims3[,1] & lyield<lyieldLims3[,2]))/length(lyield) #Proportion of data removed
# 
# #QQnorm plots
# par(mfrow=c(2,2))
# qqnorm(lyield,main='Raw'); qqline(lyield)
# qqnorm(lyield[lyield>lyieldLims[1] & lyield<lyieldLims[2]],main='3xSD'); qqline(lyield[lyield>lyieldLims[1] & lyield<lyieldLims[2]])
# qqnorm(lyield[lyield>lyieldLims2[1] & lyield<lyieldLims2[2]],main='95% quantiles'); qqline(lyield[lyield>lyieldLims2[1] & lyield<lyieldLims2[2]])
# qqnorm(lyield[lyield>lyieldLims[1] & lyield<lyieldLims[2]],main='Windowed filter'); qqline(lyield[lyield>lyieldLims[1] & lyield<lyieldLims[2]])


#Fit model --------------------------------

# #Fit non-isotropic GAM:
# kPar <- c(12,60,60,12,60,60)
# f <- sqrt(DryYield) ~ s(dist,k=kPar[1]) + s(E,N,k=kPar[2]) + s(r,k=kPar[3]) + log(pArea) #Mean model
# f2 <- ~ s(dist,k=kPar[4]) + s(E,N,k=kPar[5]) + s(r,k=kPar[6]) + log(pArea) #Variance model
# flist <- list(f,f2) #List of model formulae
# 
# a <- Sys.time() 
# mod2 <- gam(flist,data=dat,family=gaulss())
# fitTime <- paste(as.character(round(Sys.time()-a,2)),units(Sys.time()-a)) #Time taken to fit model
# # save(mod,file =  "/home/rsamuel/Documents/yield-analysis-2021/Data/unfilteredMod_3.Rdata")
# save(mod2,file =  "/home/rsamuel/Documents/yield-analysis-2021/Data/filteredMod_3.Rdata")

# load("/home/rsamuel/Documents/yield-analysis-2021/Data/unfilteredMod.Rdata") #mod
# load("/home/rsamuel/Documents/yield-analysis-2021/Data/filteredMod.Rdata") #mod2
# 
# load("/home/rsamuel/Documents/yield-analysis-2021/Data/unfilteredMod_2.Rdata") #mod
# load("/home/rsamuel/Documents/yield-analysis-2021/Data/filteredMod_2.Rdata") #mod2
# 
# load("/home/rsamuel/Documents/yield-analysis-2021/Data/unfilteredMod_3.Rdata") #mod
# load("/home/rsamuel/Documents/yield-analysis-2021/Data/filteredMod_3.Rdata") #mod2

plot(mod,pages=1,scheme=2) #Some differences between filtered/unfiltered models
plot(mod2,pages=1,scheme=2)

#Very similar shape of distance curves - intercepts shifted up and down
with(dat,data.frame(pArea=median(pArea),dist=seq(min(dist),max(dist),length.out=100),r=median(r),E=median(E),N=median(N))) %>% 
  mutate(predU_mean_mu=predict(mod,newdata=.,se.fit = TRUE)$fit[,1],
         predU_mean_se=predict(mod,newdata=.,se.fit = TRUE)$se.fit[,1],
         predU_logSD_mu=predict(mod,newdata=.,se.fit = TRUE)$fit[,2],
         predU_logSD_se=predict(mod,newdata=.,se.fit = TRUE)$se.fit[,2],
         predF_mean_mu=predict(mod2,newdata=.,se.fit = TRUE)$fit[,1],
         predF_mean_se=predict(mod2,newdata=.,se.fit = TRUE)$se.fit[,1],
         predF_logSD_mu=predict(mod2,newdata=.,se.fit = TRUE)$fit[,2],
         predF_logSD_se=predict(mod2,newdata=.,se.fit = TRUE)$se.fit[,2]) %>% 
  pivot_longer(contains('pred')) %>% separate(name,c('modType','param','moment'),sep='_') %>% 
  pivot_wider(names_from=moment,values_from = value) %>% mutate(upr=mu+se*1.96,lwr=mu-se*1.96) %>% 
  ggplot(aes(x=dist,y=mu))+geom_ribbon(aes(ymax=upr,ymin=lwr,fill=modType),alpha=0.3)+
  geom_line(aes(col=modType))+
  facet_wrap(~param,scales='free')

#Plot of mean field + distance effect
with(dat,data.frame(pArea=median(pArea),dist=dist,r=median(r),E=E,N=N)) %>%
  sample_n(10000) %>% 
  mutate(predU=predict(mod,newdata=.)[,1],predF=predict(mod2,newdata=.)[,1]) %>% 
  pivot_longer(contains('pred')) %>% 
  ggplot(aes(x=E,y=N))+geom_point(aes(col=value))+
  facet_wrap(~name,scales='free')+labs(title='Mean spatial + dist smoothers')+
  scale_colour_distiller(type='div',palette = "Spectral",direction=1)

#Plot of SD field + distance effect
with(dat,data.frame(pArea=median(pArea),dist=dist,r=median(r),E=E,N=N)) %>%
  sample_n(10000) %>% 
  mutate(predU=predict(mod,newdata=.)[,2],predF=predict(mod2,newdata=.)[,2]) %>% 
  pivot_longer(contains('pred')) %>% 
  mutate(value=exp(value)) %>% 
  ggplot(aes(x=E,y=N))+geom_point(aes(col=value),size=2)+
  facet_wrap(~name,scales='free')+labs(title='SD spatial + dist smoothers')+
  scale_colour_distiller(type='div',palette = "Spectral",direction=-1)

par(mfrow=c(2,4))
gam.check(mod); abline(0,1,col='red')
gam.check(mod2); abline(0,1,col='red')


#Still heavy space-time dependence in residuals

# pred2 <- dat %>% mutate(pArea=median(pArea)) %>% #Set pArea to its median value
#   predict(mod2,newdata = .,type='response') 
# 
# dat %>% mutate(p=pred2[,1],actual=p+resid(mod2)) %>% #Residuals across time
#   mutate(across(c(p,actual),~.^2)) %>% 
#   filter(r<2500) %>% 
#   ggplot(aes(x=r))+
#   geom_point(aes(y=actual),col='red')+
#   geom_point(aes(y=p))

library(gstat)
library(spacetime) 
dat_sp <- dat %>% mutate(resid=resid(mod)) %>% select(ID,resid) %>% 
  as_Spatial()
v <- variogram(resid~1, data=dat_sp,width=5,cutoff=300) #Variogram
p5 <- v %>% 
  ggplot()+geom_point(aes(x=dist,y=gamma))+
  labs(title='Residual variogram',x='Distance',y='Semivariance')
res_acf <- acf(dat_sp$resid,lag.max=300,type='correlation',plot=FALSE)
p6 <- with(res_acf,data.frame(acf=acf,lag=lag)) %>% 
  ggplot(aes(x=lag,y=acf))+geom_col()+
  labs(x='Time Lag',y='Autocorrelation',title='Residual autocorrelation')

  
