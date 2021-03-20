#Try loading yield data into R, and run a simple GAM


# Load everything ---------------------------------------------------------

library(tidyverse)
theme_set(theme_bw())
library(mgcv)
library(sf)
library(beepr)

# datLoc <- "C:\\Users\\Samuel\\Documents\\Ag Leader Technology\\SMS\\Export\\Trent Clark\\2019\\SW02.csv" #Multivac
datLoc <- "/media/rsamuel/Storage/geoData/Rasters/yieldData/csv files/Trent Clark/2019/SW 02.csv" #Galpern machine

#Name of field
fieldName <- unlist(strsplit(datLoc,split='/'))
fieldName <- paste(fieldName[c((length(fieldName)-2):length(fieldName))],collapse='_')
fieldName <- gsub('\\.csv','',fieldName)

source('helperFunctions.R')

dat <- read.csv(datLoc,stringsAsFactors=TRUE,fileEncoding='latin1') %>% 
  rename_with(.fn = ~gsub('..L.ha.$','_lHa',.x)) %>%
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
  
fieldEdge <- dat %>% st_union() %>% st_buffer(dist=10) %>% #10m buffer around polygons
  st_cast('LINESTRING')

dat <- dat %>%  #Distance from edge of field
  mutate(dist=as.numeric(st_distance(.,fieldEdge))[1:nrow(.)]) %>% 
  mutate(dist=dist-min(dist)) #Shrink to 0

yieldMap <- ggplot(dat)+
  geom_sf(aes(fill=log(DryYield),col=log(DryYield)))+
  geom_sf(data=fieldEdge,col='red')+
  labs(title=fieldName)
ggsave(paste0('./Figures/YieldMaps/',fieldName,'.png'),yieldMap,height=8,width=8)  

# Look at data ------------------------------------------------------------

dat %>% #Mostly done on Nov 4 2019
  ggplot()+geom_sf(aes(col=DryYield),alpha=0.7)+ #Yield over entire field
  facet_wrap(~Date,ncol=1)+scale_colour_continuous(type='viridis')

#Dates are split between two places in data 
dat %>% ggplot(aes(x=r,y=Date))+geom_point()

#Was this done on 2 combines?
dat %>% ggplot()+
  geom_sf(aes(col=r),alpha=0.7) + #geom_path()+
  facet_wrap(~Date,ncol=1)

dat %>% ggplot(aes(x=r,y=ID))+geom_point()+ #ID number seems useful, not Pass.Number
  facet_wrap(~Date)

dat %>% filter(Date!='2019-10-31') %>% ggplot()+
  geom_sf(aes(col=Pass,fill=Pass),alpha=0.7) #+ geom_path()

dat %>% ggplot(aes(x=r,y=DryYield,col=Pass))+geom_point()+
  facet_wrap(~Date,scales='free_x')

#Distance from edge of field
dat %>% 
  ggplot()+geom_sf(aes(col=sqrt(dist),fill=sqrt(dist)))+
  geom_sf(data=fieldEdge,col='red')



# Additive model ----------------------------------------------------------

library(parallel)
detectCores()
cl <- makeCluster(12)

#Isotropic model:

#DryYield ~ distance from edge, polygon area (turning + speed), geographic smoother
f <- sqrt(DryYield) ~ s(dist,k=10) + s(E,N,k=60) + s(r,k=30) + log(pArea) #Mean model
m1 <- bam(f,cluster=cl,data=dat)

summary(m1)
par(mfrow=c(2,2)); gam.check(m1); abline(0,1,col='red'); par(mfrow=c(1,1))
plot(m1,scheme=2,all.terms=TRUE,too.far=0.01,pages=1)

stopCluster(cl) 

#Non-isotropic model - much better than first one
f2 <- ~ s(dist,k=6) + s(E,N,k=60) + s(r,k=60) + log(pArea) #Variance model
flist <- list(f,f2)

a <- Sys.time()
m2 <- gam(flist,data=dat,  #Gaussian location-scale
          family=gaulss())
Sys.time()-a #Takes about 2 mins
beep(1)

summary(m2)
plot(m2,scheme=2,too.far=0.01,pages=1,all.terms=TRUE)
par(mfrow=c(2,2)); gam.check(m2); abline(0,1,col='red'); par(mfrow=c(1,1))

dat$res2 <- resid(m2) #Add residuals to df
dat$pred <- m2$fit[,1] #Predicted value
dat$sePred <- m2$fit[,2] #SE of prediction


#OK plots, but still messed up by small yield polygons
p1 <- ggplot(dat)+geom_sf(aes(fill=pred,col=pred))+labs(title='Predicted')
p2 <- ggplot(dat)+geom_sf(aes(fill=sePred,col=sePred))+labs(title='SE of Prediction')
p3 <- ggplot(dat)+geom_sf(aes(fill=res2,col=res2))+labs(title='Residuals')
ggarrange(p1,p2,p3)

#Create prediction grid
fieldEdge2 <- st_cast(fieldEdge,'POLYGON') %>% st_buffer(dist=-10) #Shrink to edge of polygons 
grid <-  st_make_grid(fieldEdge2,what='centers',cellsize=10,square=TRUE) 
grid <- grid[st_covers(fieldEdge2,grid)[[1]],] %>% ggplot()+geom_sf() #Chop out grid locations that aren't in field boundary

st_distance(st_cast(fieldEdge2,'LINESTRING'),grid)



#Temporal plot of residuals
ggplot(dat)+geom_point(aes(x=ID,y=res2))+labs(title='Residuals')
acf(dat$res2) #Bunch of autocorrelation in residuals







#Other distributions

a <- Sys.time() 
m3 <- gam(flist,data=dat, #Gamma location-scale 
          family=gammals())
Sys.time()-a #Takes about 1.4 mins
beep(2)

summary(m3)
plot(m3,scheme=2,too.far=0.01,pages=1,all.terms=TRUE)
par(mfrow=c(2,2)); gam.check(m3); abline(0,1,col='red'); par(mfrow=c(1,1))
 
# samp <- sample(1:nrow(dat),round(nrow(dat)/3)) #Use only random 1/3rd of data
# a <- Sys.time() #Gumbel
# m4 <- gam(flist,data=dat[samp,], 
#           family=gumbls())
# Sys.time()-a  #Went over 20mins without converging, even with only 1/3rd of data (4103 rows). Not sure this is going to work
# beep(2)

# samp <- sort(sample(1:nrow(dat),round(nrow(dat)/4))) #Use only random 1/4th of data (3077 samples)
# sampDat <- dat[samp,]
# flist2 <- list(f,f2,~1) #Location, scale, and shape (~1) parameter
# a <- Sys.time() 
# m4 <- gam(flist2,data=dat, #Generalized extreme value location-scale
#           family=gevlss())
# Sys.time()-a  #20 mins on entire dataset
# beep(2)
# summary(m4)
# plot(m4,scheme=2,too.far=0.01,pages=1,all.terms=TRUE)
# par(mfrow=c(2,2)); gam.check(m4); abline(0,1,col='red'); par(mfrow=c(1,1))


#Compare different models

c(m1$aic,m2$aic) #Variance very clearly non-constant. (sqrt)-Gaussian still appears best for the moment
#AIC scores for sqrt data:  55.32474 (gaussian) -21245.06753 (gaulss) -18560.69299 (gammals) -17427.36314 (gevlss)

#Bottom line: 
# - non-isotropic model required
# - Gaussian and Gamma both have long tails for residuals. gumbls didn't converge, even after 20 mins.

# save(m1,m2,file=paste('./Models/',fieldName,'.Rdata'))







