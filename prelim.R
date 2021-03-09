#Try loading yield data into R, and run a simple GAM


# Load everything ---------------------------------------------------------

library(tidyverse)
theme_set(theme_bw())
library(mgcv)
library(sf)

# datLoc <- "C:\\Users\\Samuel\\Documents\\Ag Leader Technology\\SMS\\Export\\Trent Clark\\2019\\SW02.csv" #Multivac
datLoc <- "/media/rsamuel/Storage/geoData/Rasters/yieldData/csv files/Trent Clark/SW02.csv" #Galpern machine

source('helperFunctions.R')

dat <- read.csv(datLoc,stringsAsFactors=TRUE) %>% 
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
dat %>% ggplot()+geom_sf(aes(col=dist,fill=dist))+geom_sf(data=fieldEdge,col='red')

# Deal with overlapping point data ----------------------------------------

#Problem: some points are registered as 0 yield because combine drove over areas that were already harvested.
#Possible solutions: 

dat %>% filter(Date=='2019-10-31',Pass==1) %>% 
  ggplot() + geom_sf()
  # ggplot(aes(x=Lon,y=Lat))+
  # geom_point(alpha=0.7) + geom_path()

dat %>% filter(Date=='2019-10-31',Pass==1) %>% 
  slice(1:3) %>% ggplot()+geom_sf(alpha=0.3)

dat %>% filter(Date=='2019-10-31') %>%
  ggplot()+geom_sf(alpha=0.5)

dat %>% filter(Date=='2019-10-31') %>% 
  select(ID,Pass,DryYield) %>% 
  ggplot()+geom_sf(aes(fill=DryYield))


# Additive model ----------------------------------------------------------

library(parallel)
detectCores()
cl <- makeCluster(12)

#Isotropic model:

#DryYield ~ distance from edge, polygon area (turning + speed), geographic smoother
f <- sqrt(DryYield) ~ s(dist,k=10) + s(E,N,k=60) + s(r,k=30) + log(pArea) 
m1 <- bam(f,cluster=cl,data=dat)

summary(m1)
par(mfrow=c(2,2)); gam.check(m1); abline(0,1,col='red'); par(mfrow=c(1,1))
plot(m1,scheme=2,all.terms=TRUE,too.far=0.01,pages=1)

#Non-isotropic model - much better than first one
f2 <- ~ s(pArea) + s(dist,k=6) + s(E,N,k=60) + s(r,k=60) #Variance model
flist <- list(f,f2)

m2 <- gam(flist,data=dat, #Takes about 60 seconds to fit
          family=gaulss())
summary(m2)
plot(m2,scheme=2,too.far=0.01,pages=1,all.terms=TRUE)
par(mfrow=c(2,2)); gam.check(m2); abline(0,1,col='red'); par(mfrow=c(1,1))

c(m1$aic,m2$aic) #Variance very clearly non-constant

#Bottom line: 
# - non-isotropic model required
# - need to examine model family more. gaussian has long tails for residuals, so perhaps gammals or gevlss?

stopCluster(cl)


