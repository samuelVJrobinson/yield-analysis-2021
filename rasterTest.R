#Test analysis using discrete (rasterized) data
# Reduces processing time from ~20 mins to ~3 mins

library(tidyverse)
theme_set(theme_bw())
library(ggpubr)
library(mgcv)
library(sf)
library(beepr)
source('helperFunctions.R')

datSource <- read.csv('./Data/datSource.csv') %>% #Read previous datasource file
  mutate(boundaryComplete=file.exists(boundaryPath)) %>% #Has boundary been made already?
  mutate(modelComplete1=file.exists(modelPath1)) %>% #Has model1 already been run?
  mutate(modelComplete2=file.exists(modelPath2)) %>% #Has model2 already been run?
  mutate(modelComplete0=file.exists(modelPath0)) %>% #Has model0 already been run?
  mutate(location=ifelse(grepl('(Alvin|Dean)',filename),'Southern','Central')) #Field locations (southern or central AB)

use <- which(datSource$filename == 'Trent_Clark W 34 2014')
load(datSource$modelPath2[use]) #Load model

#Spatial smoothers

library(sf)
fieldBoundary <- st_read(datSource$boundaryPath[use]) #Get boundary
fieldBoundaryType <- st_read(datSource$boundaryPath2[use])# %>% mutate(type=factor(type)) #Get boundary
crs <- st_crs(fieldBoundary) #Save CRS
fieldBoundary <- st_sfc(st_polygon(lapply(fieldBoundary$geometry,function(x) st_coordinates(x)[,c('X','Y')])),crs=crs) #Fix hole geometry
hexGrid <- st_make_grid(fieldBoundary,square=FALSE,n=75)  #Make hexagonal grid
hexGrid <- hexGrid[sapply(st_within(st_centroid(hexGrid),fieldBoundary),function(x) length(x)==1)] %>%  #Strip out points outside polygon
  st_sf() #Set as sf object
hexGrid %>% ggplot()+geom_sf()+geom_sf(data=fieldBoundary,fill=NA,col='red') #Looks OK

dat <- read.csv(datSource$dataPath[use],stringsAsFactors=TRUE,fileEncoding='latin1') %>% 
  st_as_sf(coords=c('Longitude','Latitude')) %>% #Add spatial feature info
  st_set_crs(4326) %>% st_transform(3401) %>% #Lat-lon -> UTM
  transmute(Speed=Speed.km.h.,DryYield=Yld.Mass.Dry..tonne.ha.,Dist=Distance.m.,Swath=Swth.Wdth.m.,
            pArea=Swath*Speed*1/3.6,PassNum=Pass.Num,ID=Obj..Id,
            bearingDiff=Track.deg.-lag(Track.deg.),
            E=st_coordinates(.)[,1],N=st_coordinates(.)[,2],
            r=1:n()) #row number
meanE <- mean(dat$E); meanN <- mean(dat$N)
dat <- dat %>% mutate(E=E-meanE,N=N-meanN) %>% #Center coordinates
  filter(pArea>quantile(pArea,0.05),pArea<quantile(pArea,0.95), #Filter large/small pArea
         DryYield>quantile(DryYield,0.05),DryYield<quantile(DryYield,0.95), #Filter extreme yields
         Speed>quantile(Speed,0.05),Speed<quantile(Speed,0.95) #Filter high and low speeds
  )

temp <- st_within(dat,hexGrid,sparse=TRUE)
temp[sapply(temp,length)==0] <- 0
temp <- unlist(temp)

#Yield
temp2 <- tapply(dat$DryYield,temp,function(x) if(length(x)==0) 0 else median(x))
hexGrid$DryYield <- 0.1
hexGrid$DryYield[as.numeric(names(temp2))] <- unname(temp2)
#Order
temp2 <- tapply(dat$r,temp,function(x) if(length(x)==0) NA else median(x,na.rm=TRUE))
hexGrid$r <- NA
hexGrid$r[as.numeric(names(temp2))] <- unname(temp2)

# ggplot(hexGrid)+geom_sf(aes(fill=DryYield),col=NA)
ggplot(hexGrid)+geom_sf(aes(fill=r),col=NA)

# fieldEdge <- read_sf(boundaryPath) %>% st_cast('MULTILINESTRING') #Read in boundary file
boundaryTypes <- unique(fieldBoundaryType$type) #Get boundary types

#Get distance and type of closest boundary
hexGrid <- hexGrid %>% #st_centroid() %>% 
  bind_cols(.,data.frame(t(apply(st_distance(st_centroid(.),fieldBoundaryType),1,function(x) c(dist=min(x,na.rm=TRUE),boundaryType=fieldBoundaryType$type[which.min(x)]))))) %>%
  mutate(dist=as.numeric(dist),boundaryType=factor(boundaryType)) %>% 
  mutate(E=st_coordinates(st_centroid(.))[,1],N=st_coordinates(st_centroid(.))[,2]) %>% 
  mutate(E=E-meanE,N=N-meanN) #Center coordinates

#Make model formula for non-isotropic gam
kPar <- c(12,60,60,12,60,60)

f2 <- paste0('~ s(dist,k=',kPar[1],',bs="ts",by=boundaryType) + s(E,N,k=',kPar[2],') + s(r,k=',kPar[3],')')
f <- paste0('sqrt(DryYield)~ s(dist,k=',kPar[4],',bs="ts",by=boundaryType) + s(E,N,k=',kPar[5],') + s(r,k=',kPar[6],')')
flist <- list(as.formula(f),as.formula(f2)) #List of model formulae

a <- Sys.time() #Takes about 3 mins
mod <- gam(flist,data=hexGrid,family=gaulss())
fitTime <- paste(as.character(round(Sys.time()-a,2)),units(Sys.time()-a)) #Time taken to fit model

summary(mod)

