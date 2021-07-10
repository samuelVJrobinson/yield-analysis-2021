library(mgcv)
library(TMB)
library(ggplot2)

set.seed(1)

#Make 20000 sampling locations
N <- 20000
E <- runif(N,-10,10)
N <- runif(N,-10,10)
dat <- data.frame(y=0,E,N) #Fill in y with zeros for now

s0 <- gam(y~s(E,N,bs='ts',k=60),data=dat,fit=FALSE) #Create model matrix
modMat <- s0$X[,-1] #Model matrix (no intercept)
penMat <- s0$smooth[[1]]$S[[1]] #Penalty matrix
penMat <- as(penMat,'sparseMatrix') #Change to sparse matrix
penDim <- nrow(penMat) #Dimensions of penalty matrix

b0 <- 1 #Intercept
(smoothCoefs <- rnorm(ncol(modMat))) #Smoothing coefficients
sigma <- 2 #Residual SD
yhat <- b0 + modMat %*% smoothCoefs #Mean values at each location
dat$y <- rnorm(N,yhat,sigma) #Generate data
dat$field <- modMat %*% smoothCoefs #Spatial field

ggplot(dat)+geom_point(aes(x=E,y=N,col=field)) #Plot spatial field (no noise)
ggplot(dat)+geom_point(aes(x=E,y=N,col=y)) #Plot data

#Input for model
datList <- list(y=dat$y, smoothMat=modMat, penaltyMat=penMat, penaltyDim=penDim) #Data
parList <- list(b0=0, smoothCoefs=rep(0,ncol(modMat)), log_lambda=0, log_sigma=0) #Parameters

compile("gam_tmb_tp.cpp")
dyn.load(dynlib("gam_tmb_tp"))
obj <- MakeADFun(data = datList, parameters = parList, random=c('smoothCoefs'), DLL = "gam_tmb_tp")

#Memory usage spikes from 4.5GB to 16GB during Optimizing Tape stage
opt <- nlminb(obj$par,obj$fn,obj$gr) #Run model
report <- sdreport(obj)

opt$convergence==0 #Did model converge?
report #Looks OK. b0 and sigma are very close to actual values (1 and 2)

#Estimates smoothing coefficients fairly well
ggplot(data.frame(predicted=report$par.random,actual=smoothCoefs))+
  geom_point(aes(x=predicted,y=actual))+geom_abline(intercept = 0,slope=1)+
  labs(title='Smoothing coefficients')

