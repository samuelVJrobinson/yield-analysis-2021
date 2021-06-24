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
boundaryPath2 <- "/home/rsamuel/Documents/yield-analysis-2021/Figures/FieldBoundarySegments/Alvin_French C41 2020 boundary.shp"

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

fieldEdge <- read_sf(boundaryPath) #Read in boundary polygon
fieldEdgeType <- read_sf(boundaryPath2) #Read in boundary type linestrings

#Get distance and type of closest boundary
dat <- dat %>% 
  bind_cols(.,data.frame(t(apply(st_distance(dat,fieldEdgeType),1,function(x) c(dist=min(x,na.rm=TRUE),boundaryType=fieldEdgeType$type[which.min(x)]))))) %>% 
  mutate(dist=as.numeric(dist),boundaryType=factor(boundaryType))

# Model 0A: Yield ~ s(Dist) -----------------------------------------------------

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
datList <- list(yield=sqrt(dat$DryYield), logArea=log(dat$pArea), meanLogArea = mean(log(dat$pArea)),
                smoothMat=modMat, penaltyMat=penMat, penaltyDim=penDim, newSmoothMat=predModMat) #Data
parList <- list(b0=0, b_area = 0, smoothCoefs=rep(0,ncol(modMat)), log_lambda=0, log_sigma=0) #Parameters

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
acf(obj$report()$residuals) #Autocorrelation in residuals 
dat %>% mutate(resid=obj$report()$residuals) %>% #Still some spatial autocorrelation, but not as bad as before
  select(.,ID,resid) %>% 
  # slice(round(seq(1,nrow(.),length.out=3000))) %>% 
  as_Spatial(.) %>% 
  gstat::variogram(resid~1, data=.,width=2,cutoff=100) %>% #Variogram
  ggplot(aes(x=dist,y=gamma))+geom_point()+geom_line()+
  labs(title='Residual variogram',x='Distance',y='Semivariance')

#Predictions
data.frame(dist=predDist,fit=report$value,se.fit=report$sd) %>% 
  mutate(upr=fit+1.96*se.fit,lwr=fit-1.96*se.fit) %>% 
  mutate(across(c(fit,upr,lwr),~.x^2)) %>% 
  ggplot(aes(x=dist))+geom_ribbon(aes(ymax=upr,ymin=lwr),alpha=0.3)+
  geom_line(aes(y=fit))+labs(x='Distance from edge (m)', y='Dry Yield (t/ha)')

#Compare to GAM model
m1 <- gam(sqrt(DryYield)~log(pArea) + s(dist,bs='cr',k=12) + s(E,N,k=60) + s(r,k=60) ,data=dat)
summary(m1)
res <- residuals(m1,type='response')
qqnorm(res); qqline(res)
acf(res) #Autocorrelation in residuals 
dat %>% mutate(resid=res) %>% #Weird semivariance plot
  select(.,ID,resid) %>% as_Spatial(.) %>% 
  gstat::variogram(resid~1, data=.,width=5,cutoff=100) %>% #Variogram
  ggplot(aes(x=dist,y=gamma))+geom_point()+geom_line()+
  labs(title='Residual variogram',x='Distance',y='Semivariance')

#Fit non-isotropic GAM as test
kPar <- c(12,60,60,12,60,60)
f <- sqrt(DryYield) ~ s(dist,k=kPar[1]) + s(E,N,k=kPar[2]) + s(r,k=kPar[3]) + log(pArea) #Mean model
f2 <- ~ s(dist,k=kPar[4]) + s(E,N,k=kPar[5]) + s(r,k=kPar[6]) + log(pArea) #Variance model
flist <- list(f,f2) #List of model formulae

a <- Sys.time() #Takes about 8 mins for a 50000 point model
mod <- gam(flist,data=dat,family=gaulss())
fitTime <- paste(as.character(round(Sys.time()-a,2)),units(Sys.time()-a)) #Time taken to fit model

debugonce(residuals)
res <- residuals(mod,type='response')
# res <- (mod$fitted.values[,1]-mod$y)*mod$fitted.values[,2]
qqnorm(res); qqline(res)
acf(res) #Autocorrelation in residuals 
p1 <- dat %>% mutate(resid=res) %>% 
  select(.,ID,resid) %>% as_Spatial(.) %>% 
  gstat::variogram(resid~1, data=.,width=10,cutoff=200) %>% #Variogram
  ggplot(aes(x=dist,y=gamma))+geom_point()+geom_line()+
  labs(title='Response residuals',x='Distance',y='Semivariance')

res <- residuals(mod)
p2 <- dat %>% mutate(resid=res) %>% 
  select(.,ID,resid) %>% as_Spatial(.) %>% 
  gstat::variogram(resid~1, data=.,width=10,cutoff=200) %>% #Variogram
  ggplot(aes(x=dist,y=gamma))+geom_point()+geom_line()+
  labs(title='Deviance residuals',x='Distance',y='Semivariance')
library(ggpubr)
ggarrange(p1,p2)

# Model 0AA: Yield ~ s(Dist,by='type') -----------------------------------------

#Get smoother specifications
s0 <- gam(DryYield ~ s(dist,bs='cs',by=boundaryType,k=10), data = dat, fit = FALSE)

#Input
modMat <- s0$X[,-1] #Model matrix (no intercept)
penMat <- s0$smooth[[1]]$S[[1]] #Penalty matrix (sparse) - all smoother use the same penalty matrix, so only need first one
penMat <- as(penMat,'sparseMatrix')
penDim <- nrow(penMat) #Dimensions of penalty matrix
numSmooths <- length(s0$smooth) #Number of smoothers - 4 boundary types, so 4 smoothers

predDist <- lapply(round(with(dat,tapply(dist,boundaryType,max))),function(x) seq(0,x,by=5))
predDist <- data.frame(DryYield=-1,dist=unname(unlist(predDist)),boundaryType=factor(rep(names(predDist),sapply(predDist,length)),level=levels(dat$boundaryType)))
predModMat <- gam(formula(s0),data = predDist,fit=FALSE)$X[,-1] #New smoothing model matrix

#Input for model
datList <- list(yield=sqrt(dat$DryYield), logArea=log(dat$pArea), meanLogArea = mean(log(dat$pArea)),
                smoothMat=modMat, penaltyMat=penMat, penaltyDim=penDim, numSmooths=numSmooths, newSmoothMat=predModMat) #Data
parList <- list(b0=0, b_area = 0, smoothCoefs=rep(0,ncol(modMat)), log_lambdas=rep(0,numSmooths), log_sigma=0) #Parameters

dyn.unload(dynlib("mod0AA"))
compile("mod0AA.cpp")
dyn.load(dynlib("mod0AA"))
obj <- MakeADFun(data = datList, parameters = parList, random=c('smoothCoefs'), DLL = "mod0AA")

obj$par
obj$fn() #Works
opt <- nlminb(obj$par,obj$fn,obj$gr) #Runs
report <- sdreport(obj) #Estimates and SEs

report #Looks OK
report$value
report$par.fixed
report$par.random

#Predictions
predDist %>% select(-DryYield) %>% 
  mutate(fit=report$value,se.fit=report$sd) %>% 
  mutate(upr=fit+1.96*se.fit,lwr=fit-1.96*se.fit) %>% 
  mutate(across(c(fit,upr,lwr),~.x^2)) %>% 
  ggplot(aes(x=dist))+geom_ribbon(aes(ymax=upr,ymin=lwr),alpha=0.3)+
  facet_wrap(~boundaryType)+
  geom_line(aes(y=fit))+labs(x='Distance from edge (m)', y='Dry Yield (t/ha)')

# Model 0B: Yield ~ SPDE(E,N) --------------------------------------------------

fieldEdgePoly <- read_sf(boundaryPath) #Polygon version of boundary
# ggplot()+geom_sf(dat=fieldEdgePoly[1,])

outerBoundaries <- fieldEdgePoly %>% st_union() #Get outermost boundaries
outerBoundaries2 <- outerBoundaries %>% st_buffer(dist = 150) #150m buffer around outermost boundaries
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

# meshbuilder()

mesh  <- inla.mesh.2d( #Construct triangulated mesh
  # loc=as(dat,'Spatial'),
  boundary = list(outerBoundaries,outerBoundaries2),
  max.edge=c(30,100),
  cutoff=c(20,50),
  min.angle = c(21)
)
plot(mesh) #Looks OK - this makes a matrix of about N x 3000 
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

a <- Sys.time()
obj <- MakeADFun(data = datList, parameters = parList, random=c('spdeCoefs'), DLL = "mod0B") #Takes about 15 secs to run if A = 20000 x 3000
Sys.time()-a
obj$par
obj$fn() #Works
a <- Sys.time()
opt <- nlminb(obj$par,obj$fn,obj$gr)
b <- Sys.time()
report <- sdreport(obj) #Estimates and SEs 
b-a #2 mins to fit
Sys.time()-b #34 secs to report
save(obj,opt,report,file = 'mod0B.Rdata')
# load('mod0B.Rdata')
# dyn.load(dynlib("mod0B"))

opt$convergence == 0 #Convergence?
report #Looks OK
report$value #Correlation drops below 0.1 at ~60 meters
report$par.fixed
report$par.random

#Residuals are still horrible
qqnorm(obj$report()$residuals); abline(0,1)
acf(obj$report()$residuals) #Temporal autocorrelation in residuals still evident

dat %>% mutate(resid=obj$report()$residuals) %>% #Still some spatial autocorrelation, but not as bad as before
  select(.,ID,resid) %>% 
  slice(round(seq(1,nrow(.),length.out=3000))) %>% 
  as_Spatial(.) %>% 
  gstat::variogram(resid~1, data=.,width=2,cutoff=300) %>% #Variogram
  ggplot(aes(x=dist,y=gamma))+geom_point()+geom_line()+
  labs(title='Residual variogram',x='Distance',y='Semivariance')

#Look at spatial field
projMesh  <- inla.mesh.projector(mesh)

par(mfrow=c(1,2))
latentFieldMAP  <- (report$par.random/exp(report$par.fixed[which(names(report$par.fixed)=="log_tau")])+report$par.fixed[which(names(report$par.fixed)=="b0")])^2
fields::image.plot(projMesh$x,projMesh$y, inla.mesh.project(projMesh, latentFieldMAP),col =  colorRampPalette(c("red","blue","yellow"))(12),
           xlab = 'Easting', ylab = 'Northing',
           main = "MAP estimate of the spatial latent field",
           cex.lab = 1.1,cex.axis = 1.1, cex.main=1, cex.sub= 1.1)
# contour(projMesh$x, projMesh$y,inla.mesh.project(projMesh, latentFieldMAP) ,add = T,labcex  = 1,cex = 1)
with(outerBoundaries, for(i in 1:nrow(loc)) lines(loc[idx[i,],],col='white'))

latentFieldSD <- sqrt(report$diag.cov.random)/exp(report$par.fixed[which(names(report$par.fixed)=="log_tau")])
fields::image.plot(projMesh$x,projMesh$y, inla.mesh.project(projMesh, latentFieldSD),col =  colorRampPalette(c("red","blue", "yellow"))(12),
           xlab = 'Easting', ylab = 'Northing',
           main = "Standard deviation of the estimated spatial latent field",
           cex.lab = 1.1,cex.axis = 1.1, cex.main=1, cex.sub= 1.1)
with(outerBoundaries, for(i in 1:nrow(loc)) lines(loc[idx[i,],],col='white'))

# Model 0C: Yield ~ SPDE(Time) --------------------------------------------

# #First approach - this fits, but sdreport never finishes (too many random effects)
# #Input for model
# datList <- list(yield=sqrt(dat$DryYield)) #Data
# parList <- list(b0=0, u=rep(0,nrow(dat)), atanh_phi=0, log_sigmaU = 0, log_sigma = 0) #Parameters
# 
# dyn.unload(dynlib("mod0C"))
# compile("mod0C.cpp")
# dyn.load(dynlib("mod0C"))
# obj <- MakeADFun(data = datList, parameters = parList, random=c('u'), DLL = "mod0C")
# 
# obj$par
# obj$fn() #Works
# opt <- nlminb(obj$par,obj$fn,obj$gr) #Takes a few seconds
# report <- sdreport(obj) #Estimates and SEs - this takes wayy to long

bLocs <- seq(1,nrow(dat),length.out=10000) #Locations of basis splines
mesh  <- inla.mesh.1d(bLocs,degree=2)

A <- inla.spde.make.A(mesh,loc=1:nrow(dat)) #Interpolation matrix (match data to mesh)
spde <- inla.spde2.matern(mesh, alpha=2) #Create Matern SPDE model from mesh
spdeMatrices <- spde$param.inla[c("M0","M1","M2")] #Matrices needed within TMB

datList <- list(yield = sqrt(dat$DryYield), A = A, spdeMatrices = spdeMatrices)
parList <- list(b0 = 0, spdeCoefs = rep(0.0, nrow(datList$spdeMatrices$M0)), log_tau   = 0, log_kappa = 0, log_sigma = 0)

dyn.unload(dynlib("mod0B")) #I think I can just use the same model for both, since they both just use GMRF with a Q matrix
compile("mod0B.cpp")
dyn.load(dynlib("mod0B"))

obj <- MakeADFun(data = datList, parameters = parList, random=c('spdeCoefs'), DLL = "mod0B") 
obj$par
obj$fn() #Works
opt <- nlminb(obj$par,obj$fn,obj$gr)  
report <- sdreport(obj) #Estimates and SEs

opt$convergence==0 #Did model converge?
report #Looks OK
report$value #Correlation drops below 0.1 at ~30 points
report$par.fixed
report$par.random

#Residuals are still horrible
qqnorm(obj$report()$residuals); abline(0,1)
acf(obj$report()$residuals) #Temporal autocorrelation in residuals still evident
dat %>% mutate(resid=obj$report()$residuals) %>% #Still some spatial autocorrelation, but not as bad as before
  select(.,ID,resid) %>% 
  slice(round(seq(1,nrow(.),length.out=3000))) %>% 
  as_Spatial(.) %>% 
  gstat::variogram(resid~1, data=.,width=2,cutoff=300) %>% #Variogram
  ggplot(aes(x=dist,y=gamma))+geom_point()+geom_line()+
  labs(title='Residual variogram',x='Distance',y='Semivariance')

#Look at temporal field
projMesh  <- inla.mesh.projector(mesh)
#Predictions from latent field map
latentFieldMAP  <- (report$par.random/exp(report$par.fixed[which(names(report$par.fixed)=="log_tau")])+report$par.fixed[which(names(report$par.fixed)=="b0")])^2
latentFieldSD <- sqrt(report$diag.cov.random)/exp(report$par.fixed[which(names(report$par.fixed)=="log_tau")])
par(mfrow=c(2,1))
plot(projMesh$x,inla.mesh.project(projMesh, latentFieldMAP),type='l',xlab='Index',ylab='Temporal Effect')
plot(projMesh$x,inla.mesh.project(projMesh, latentFieldSD),type='l',xlab='Index',ylab='SD Temporal Effect')
plot(latentFieldMAP,latentFieldSD) #Not sure why this occurs

# Model 1A: Yield ~ s(dist) + SPDE(E,N) + SPDE(Time) ------------------------------

#Get smoother parameters:
s0 <- gam(DryYield ~ s(dist,bs='cs',k=10), data = dat, fit = FALSE)
modMat <- s0$X[,-1] #Model matrix (no intercept)
penMat <- s0$smooth[[1]]$S[[1]] #Penalty matrix (sparse)
penMat <- as(penMat,'sparseMatrix')
penDim <- nrow(penMat) #Dimensions of penalty matrix
predDist <- with(dat,seq(min(dist),max(dist),length.out=50)) #Distances to predict for
predModMat <- PredictMat(s0$smooth[[1]],data = data.frame(dist=predDist)) #New smoothing model matrix

#Get spatial SPDE parameters:
fieldEdgePoly <- read_sf(boundaryPath) #Polygon version of boundary
outerBoundaries <- fieldEdgePoly %>% st_union() #Get outermost boundaries
outerBoundaries2 <- outerBoundaries %>% st_buffer(dist = 300) #300m buffer around outermost boundaries
#Interior polygons
innerPoly <- sapply(st_covered_by(fieldEdgePoly,outerBoundaries),function(x) length(x)==1)&!sapply(st_equals(fieldEdgePoly,outerBoundaries),function(x) length(x)==1)
if(any(innerPoly)){
  innerBoundaries <- fieldEdgePoly[innerPoly,] %>% st_union()
  outerBoundaries <- st_difference(outerBoundaries,innerBoundaries)
} 
#Convert to SpatialPolygons, then inla mesh segments
outerBoundaries <- as(outerBoundaries,'Spatial') %>% as.inla.mesh.segment()
outerBoundaries2 <- as(outerBoundaries2,'Spatial') %>% as.inla.mesh.segment()
#Construct triangulated spatial Mesh
spatialMesh  <- inla.mesh.2d(boundary = list(outerBoundaries,outerBoundaries2),
  max.edge=c(30,100), cutoff=c(20,50), min.angle = c(21))
plot(spatialMesh) #Looks OK
spatialA <- inla.spde.make.A(spatialMesh,loc=as(dat,'Spatial')) #Spatial interpolation matrix (match data to spatialMesh)
spatialSPDE <- inla.spde2.matern(spatialMesh, alpha=2) #Create Matern SPDE model from spatialMesh
spatialSPDEMatrices <- spatialSPDE$param.inla[c("M0","M1","M2")] #Matrices needed within TMB
dim(spatialA) #Roughly 3300 vertices

#Get temporal SPDE parameters:
bLocs <- seq(1,nrow(dat),length.out=5000) #Locations of basis splines
temporalMesh  <- inla.mesh.1d(bLocs,degree=2)
temporalA <- inla.spde.make.A(temporalMesh,loc=1:nrow(dat)) #Interpolation matrix (match data to mesh)
temporalSPDE <- inla.spde2.matern(temporalMesh, alpha=2) #Create Matern SPDE model from mesh
temporalSPDEMatrices <- temporalSPDE$param.inla[c("M0","M1","M2")] #Matrices needed within TMB

#Input for model:
#Data
datList <- list(yield=sqrt(dat$DryYield), #Transformed Data
                logArea=log(dat$pArea), meanLogArea = mean(log(dat$pArea)), #log(area) and mean(log(area))
                smoothMat=modMat, penaltyMat=penMat, penaltyDim=penDim, newSmoothMat=predModMat, #Smoother input
                spatialSPDEMatrices=spatialSPDEMatrices, spatialA = spatialA, #Spatial input
                temporalSPDEMatrices = temporalSPDEMatrices, temporalA = temporalA) #Temporal input
#Parameters
parList <- list(b0=0, #Intercept
                b_area=0,
                smoothCoefs=rep(0,ncol(modMat)), log_lambda=0, #Smooth coefficients and penalization parameter
                spdeCoefs_spatial = rep(0.0, nrow(datList$spatialSPDEMatrices$M0)), log_tau_spatial = 0, log_kappa_spatial = 0,
                spdeCoefs_temporal = rep(0.0, nrow(datList$temporalSPDEMatrices$M0)), log_tau_temporal = 0, log_kappa_temporal = 0,
                log_sigma=0 #Error variance
                ) #Parameters

dyn.unload(dynlib("mod1A")) 
compile("mod1A.cpp")
dyn.load(dynlib("mod1A"))

obj <- MakeADFun(data = datList, parameters = parList, random=c('smoothCoefs','spdeCoefs_spatial','spdeCoefs_temporal'), DLL = "mod1A") 
obj$par
# obj$fn() #Works
a <- Sys.time()
opt <- nlminb(obj$par,obj$fn,obj$gr)  #Takes about 9 mins
b <- Sys.time()
report <- sdreport(obj) #Estimates and SEs - takes about 2 mins
b-a
Sys.time()-b 

opt$convergence==0 #Did model converge?
report #Looks OK
report$value #Spatial autocorr drops at about 125m, temporal at 11 points
report$par.fixed
report$par.random

#Residuals are still horrible
res <- obj$report()$residuals#/exp(report$par.fixed[names(report$par.fixed)=='log_sigma'])
qqnorm(res); qqline(res)
acf(res) #Negative autocorrelation
dat %>% mutate(resid=res) %>% #Still some spatial autocorrelation, but not as bad as before
  select(.,ID,resid) %>% 
  slice(round(seq(1,nrow(.),length.out=3000))) %>% 
  as_Spatial(.) %>% 
  gstat::variogram(resid~1, data=.,width=2,cutoff=300) %>% #Variogram
  ggplot(aes(x=dist,y=gamma))+geom_point()+geom_line()+
  labs(title='Residual variogram',x='Distance',y='Semivariance')

#Look at spatial field
projMesh  <- inla.mesh.projector(spatialMesh)
parFixEff <- report$par.fixed
fixedNames <- names(parFixEff)
parRanEff <- report$par.random
randomNames <- names(parRanEff)

latentFieldMAP  <- ((parRanEff[randomNames=='spdeCoefs_spatial']/exp(parFixEff[fixedNames=='log_tau_spatial']))+parFixEff[fixedNames=='b0'])^2
fields::image.plot(projMesh$x,projMesh$y, inla.mesh.project(projMesh, latentFieldMAP),col =  colorRampPalette(c("red","blue","yellow"))(12),
                   xlab = 'Easting', ylab = 'Northing',
                   main = "MAP estimate of the spatial latent field",
                   cex.lab = 1.1,cex.axis = 1.1, cex.main=1, cex.sub= 1.1)
# contour(projMesh$x, projMesh$y,inla.mesh.project(projMesh, latentFieldMAP) ,add = T,labcex  = 1,cex = 1)
with(outerBoundaries, for(i in 1:nrow(loc)) lines(loc[idx[i,],],col='white'))

