#Tests of R-INLA models for yield, using SPDE approx of Matern covariance

#Models to try:
# 0A. Yield ~ s(dist) - done
# 0B. Yield ~ SPDE(E,N)
# 0C. Yield ~ SPDE(Time)
# 1A Yield ~ s(dist) + SPDE(E,N) + SPDE(Time)
# 1B Yield ~ s(dist) + SPDE(E,N) + SPDE(Time), Var ~ s(dist) + SPDE(E,N) + SPDE(Time)
# 1C. Yield ~ s(dist,by=Type) + SPDE(E,N) + SPDE(Time), Var ~ s(dist,by=Type) + SPDE(E,N) + SPDE(Time)

#Load everything --------------------

library(tidyverse)
theme_set(theme_classic())
library(sf)
library(mgcv)
library(INLA)
library(TMB)
setwd("~/Documents/yield-analysis-2021/gam_TMB")
source('../helperFunctions.R')

nSubSamp <- 20000 #Number of points to limit to
csvPath <- "/media/rsamuel/Storage/geoData/Rasters/yieldData/csv files/Alvin_French C41 2020.csv"
boundaryPath <- "/home/rsamuel/Documents/yield-analysis-2021/Figures/FieldBoundaries/Alvin_French C41 2020 boundary.shp"

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
  st_centroid() %>% #Convert back to point
  mutate(r=1:n()) %>% #row number
  mutate(Pass=factor(seqGroup(ID,FALSE))) %>% 
  group_by(Pass) %>% mutate(rGroup=1:n()) %>% ungroup() %>% 
  mutate(E=st_coordinates(.)[,1],N=st_coordinates(.)[,2]) %>% 
  mutate(E=(E-mean(E)),N=(N-mean(N))) #Center coordinates
# ggplot(dat)+geom_sf()

fieldEdge <- read_sf(boundaryPath) %>% st_cast('MULTILINESTRING') #Read in boundary file

#Distance from boundary
dat <- dat %>% mutate(dist=as.numeric(st_distance(.,fieldEdge))[1:nrow(.)]) 

# Model 0A. Yield ~ s(Dist) -----------------------------------------------------

#Get smoother specifications
s0 <- gam(DryYield ~ s(dist,bs='cs',k=10), data = dat, fit = FALSE)

#Input
modMat <- s0$X[,-1] #Model matrix (no intercept)
penMat <- s0$smooth[[1]]$S[[1]] #Penalty matrix (sparse)
penMat <- as(penMat,'sparseMatrix')
penDim <- nrow(penMat) #Dimensions of penalty matrix

predDist <- with(dat,seq(min(dist),max(dist),length.out=50)) #Distances to predict for
predModMat <- PredictMat(s0$smooth[[1]],data = data.frame(dist=predDist)) #New smoothing model matrix

#Input for model
datList <- list(yield=sqrt(dat$DryYield), smoothMat=modMat, penaltyMat=penMat, penaltyDim=penDim, newSmoothMat=predModMat) #Data
parList <- list(b0=0, smoothCoefs=rep(0,ncol(modMat)), log_lambda=0, log_sigma=0) #Parameters

dyn.unload(dynlib("mod0A"))
compile("mod0A.cpp")
dyn.load(dynlib("mod0A"))
obj <- MakeADFun(data = datList, parameters = parList, random=c('smoothCoefs'), DLL = "mod0A")

obj$par
obj$fn() #Works
opt <- nlminb(obj$par,obj$fn,obj$gr) #Runs
report <- sdreport(obj) #Estimates and SEs

report #Looks OK
report$value
report$par.fixed
report$par.random

#Check residuals - very non-normal, even with sqrt-transformation
qqnorm(obj$report()$residuals); abline(0,1)

#Predictions
data.frame(dist=predDist,fit=report$value,se.fit=report$sd) %>% 
  mutate(upr=fit+1.96*se.fit,lwr=fit-1.96*se.fit) %>% 
  mutate(across(c(fit,upr,lwr),~.x^2)) %>% 
  ggplot(aes(x=dist))+geom_ribbon(aes(ymax=upr,ymin=lwr),alpha=0.3)+
  geom_line(aes(y=fit))+labs(x='Distance from edge (m)', y='Dry Yield (t/ha)')

# Model 0B. Yield ~ SPDE(E,N) --------------------------------------------------

fieldEdgePoly <- read_sf(boundaryPath) #Polygon version of boundary
# ggplot()+geom_sf(dat=fieldEdgePoly[1,])

outerBoundaries <- fieldEdgePoly %>% st_union() #Get outermost boundaries
outerBoundaries2 <- outerBoundaries %>% st_buffer(dist = 200) #100m buffer around outermost boundaries
# ggplot()+geom_sf(dat=dat)+
#   geom_sf(dat=outerBoundaries,fill='red',alpha=0.3)+geom_sf(dat=outerBoundaries2,fill='pink',alpha=0.3)

#Get interior polygons
innerPoly <- sapply(st_covered_by(fieldEdgePoly,outerBoundaries),function(x) length(x)==1)&!sapply(st_equals(fieldEdgePoly,outerBoundaries),function(x) length(x)==1)
if(any(innerPoly)){
  innerBoundaries <- fieldEdgePoly[innerPoly,] %>% st_union()
  # ggplot(innerBoundaries)+geom_sf()
  outerBoundaries <- st_difference(outerBoundaries,innerBoundaries)
} 

#Convert to SpatialPolygons, then inla mesh segments
outerBoundaries <- as(outerBoundaries,'Spatial') %>% as.inla.mesh.segment()
outerBoundaries2 <- as(outerBoundaries2,'Spatial') %>% as.inla.mesh.segment()

mesh  <- inla.mesh.2d( #Construct triangulated mesh
  # loc=as(dat,'Spatial'),
  boundary = list(outerBoundaries,outerBoundaries2),
  max.edge=c(30,60),
  min.angle = c(21)
)
plot(mesh) #Looks OK

A <- inla.spde.make.A(mesh,loc=as(dat,'Spatial')) #Spatial interpolation matrix (match data to mesh)
spde <- inla.spde2.matern(mesh, alpha=2) #Create Matern SPDE model from mesh
spdeMatrices <- spde$param.inla[c("M0","M1","M2")] #Matrices needed within TMB

# image(A)
# image(spdeMatrices[[1]])
# image(spdeMatrices[[2]])
# image(spdeMatrices[[3]])

datList <- list(yield = sqrt(dat$DryYield), A = A, spdeMatrices = spdeMatrices)
parList <- list(b0 = 1.9, spdeCoefs = rep(0.0, nrow(datList$spdeMatrices$M0)), log_tau   = 3.13, log_kappa = -2.84, log_sigma = -0.13)

dyn.unload(dynlib("mod0B"))
compile("mod0B.cpp")
dyn.load(dynlib("mod0B"))
obj <- MakeADFun(data = datList, parameters = parList, random=c('spdeCoefs'), DLL = "mod0B") #Takes about 5 mins to run if nSubSamp == 10000, crashes if 20000

obj$par
obj$fn() #Works
a <- Sys.time()
opt <- nlminb(obj$par,obj$fn,obj$gr) 
Sys.time()-a #Takes 13 mins with 20000 samples
report <- sdreport(obj) #Estimates and SEs - takes 3 mins

report #Looks OK
report$value #Correlation drops below 0.1 at ~50 meters
report$par.fixed
report$par.random

#Residuals are still horrible
qqnorm(obj$report()$residuals); abline(0,1)

report$par.random


projMesh  <- inla.mesh.projector(mesh)
latentFieldMAP  <- report$par.random/exp(report$par.fixed[which(names(report$par.fixed)=="log_tau")])
fields::image.plot(projMesh$x,projMesh$y, inla.mesh.project(projMesh, latentFieldMAP),col =  colorRampPalette(c("red","purple", "blue"))(12),
           xlab = 'Easting', ylab = 'Northing',
           main = "MAP estimate of the spatial latent field",
           cex.lab = 1.1,cex.axis = 1.1, cex.main=1, cex.sub= 1.1)
# contour(projMesh$x, projMesh$y,inla.mesh.project(projMesh, latentFieldMAP) ,add = T,labcex  = 1,cex = 1)

with(outerBoundaries,{
  for(i in 1:nrow(loc)) lines(loc[idx[i,],],col='black')   
})
