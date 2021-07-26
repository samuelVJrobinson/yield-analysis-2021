library(mgcv)
library(TMB)
library(ggplot2)
library(ggpubr)
set.seed(1)
setwd("~/Documents/yield-analysis-2021/gam_TMB")

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
smoothCoefs <- rnorm(ncol(modMat)) #Smoothing coefficients

sigma <- 2 #Residual SD
yhat <- b0 + modMat %*% smoothCoefs #Mean values at each location
dat$y <- rnorm(N,yhat,sigma) #Generate data
dat$field <- modMat %*% smoothCoefs #Spatial field

ggplot(dat)+geom_point(aes(x=E,y=N,col=field)) #Plot spatial field (no noise)
ggplot(dat)+geom_point(aes(x=E,y=N,col=y)) #Plot data

#Input for model

#Fit model using 5 "batches" of data
nbatches <- 5
batches <- rep(1:nbatches,each=nrow(dat)/nbatches)

#Data, only using first batch
datList <- list(y=dat$y[batches==1], smoothMat=modMat[batches==1,], penaltyMat=penMat, penaltyDim=penDim,usePenalty=1.0) 
parList <- list(b0=0, smoothCoefs=rep(0,ncol(modMat)), log_lambda=0, log_sigma=1) #Parameters

compile("gam_tmb_batches2.cpp")
dyn.load(dynlib("gam_tmb_batches2"))
obj <- MakeADFun(data = datList, parameters = parList, random=c('smoothCoefs'), DLL = "gam_tmb_batches2")

obj$fn() #Objective function works. Spike in memory during Optimizing Tape much reduced from before

#Make new objective function that accumulates values from all batches
oFun <- function(params){ 
  nll <- obj$fn(params)  #NLL from from first batch
  
  for(i in 2:nbatches){ #For each subsequent batch
    obj$env$data$y <- dat$y[batches==i] #Update data
    obj$env$data$smoothMat <- modMat[batches==i,] 
    obj$env$data$usePenalty <- 0
    nll <- nll + obj$fn(params) #Increment NLL
  }
  obj$env$data$y <- dat$y[batches==1] #Reset data to batch 1
  obj$env$data$smoothMat <- modMat[batches==1,] 
  obj$env$data$usePenalty <- 1.0
  return(nll)
}

#Make new gradient function that accumulates values from all batches
oGrad <- function(params){
  grad <- obj$gr(params)  #NLL from from first batch
  for(i in 2:nbatches){ #For each subsequent batch
    obj$env$data$y <- dat$y[batches==i] #Update data
    obj$env$data$smoothMat <- modMat[batches==i,] 
    obj$env$data$usePenalty <- 0 #Turn off penalization
    grad <- grad + obj$gr(params) #Increment gradient
  }
  obj$env$data$y <- dat$y[batches==1] #Reset data to batch 1
  obj$env$data$smoothMat <- modMat[batches==1,] 
  obj$env$data$usePenalty <- 1.0 #Turn penalization back on
  return(grad)
}

opt <- nlminb(obj$par,oFun,oGrad) #Run model using obj$par and new objective/gradient functions
report <- sdreport(obj)

opt$convergence==0 #Did model converge?
report #Looks OK. b0 and logsigma are very close to actual values (1 and 0.693)

#Estimates smoothing coefficients fairly well
p1 <- ggplot(data.frame(var=c('b0','logLambda','logSigma'),estimate=report$par.fixed,stdev=sqrt(diag(report$cov.fixed)),actual=c(b0,NA,log(sigma))),aes(x=var,y=estimate))+
  geom_pointrange(aes(ymax=estimate+stdev*1.96,ymin=estimate-stdev*1.96))+
  geom_point(aes(y=actual),col='red',size=2)+
  facet_wrap(~var,scales='free')+labs(x='fixed effect')+theme(axis.text.x = element_blank())

p2 <- ggplot(data.frame(predicted=report$par.random,actual=smoothCoefs))+
  geom_point(aes(x=predicted,y=actual))+geom_abline(intercept = 0,slope=1)+
  labs(title='Smoothing coefficients')

(p <- ggarrange(p1,p2,ncol=2))

#Compare to non-batch model

#Input for model
datList2 <- list(y=dat$y, smoothMat=modMat, penaltyMat=penMat, penaltyDim=penDim,flag=1) #Data
parList2 <- list(b0=0, smoothCoefs=rep(0,ncol(modMat)), log_lambda=0, log_sigma=1) #Parameters

compile("gam_tmb_tp.cpp")
dyn.load(dynlib("gam_tmb_tp"))
obj2 <- MakeADFun(data = datList2, parameters = parList2, random=c('smoothCoefs'), DLL = "gam_tmb_tp")

#Memory usage spikes from 4.5GB to 16GB during Optimizing Tape stage
opt2 <- nlminb(obj2$par,obj2$fn,obj2$gr) #Run model
report2 <- sdreport(obj2)

opt2$convergence==0 #Did model converge?
report2 #Looks OK. b0 and sigma are very close to actual values (1 and 2)

orig <- data.frame(par=names(opt$par),est=c(1,NA,log(2)),sd=NA,upr=NA,lwr=NA,type='actual') 
pars <- data.frame(par=names(opt$par),est=opt$par,sd=sqrt(diag(report$cov.fixed))) %>% mutate(upr=est+sd*1.96,lwr=est-sd*1.96,type='batches')
pars2 <- data.frame(par=names(opt$par),est=opt2$par,sd=sqrt(diag(report2$cov.fixed))) %>% mutate(upr=est+sd*1.96,lwr=est-sd*1.96,type='no_batches')

rbind(orig,pars,pars2) %>% mutate(par=factor(par,levels=c('b0','log_lambda','log_sigma'))) %>% 
  ggplot(aes(x=par,y=est,col=type))+geom_pointrange(aes(ymax=upr,ymin=lwr),position=position_dodge(width=0.4)) +
  facet_wrap(~par,scales='free')+labs(y='Estimate')+
  theme(axis.text.x = element_blank(), axis.title.x = element_blank(), axis.ticks.x = element_blank())+
  scale_colour_manual(values=c('black','red','blue'))



