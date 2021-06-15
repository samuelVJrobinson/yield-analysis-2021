#GAM testing using TMB

setwd("~/Documents/yield-analysis-2021/gam_TMB")

library(Matrix) #Use sparse matrices
library(TMB)
library(mgcv)
library(tidyverse)
library(beepr)
theme_set(theme_classic())

library(gamair)
data(engine)

#Thin plate spline  - single smoother  ---------------

s0 <- gam(wear~s(size,k=8,bs='ts'),data=engine,fit=FALSE) #Problems start to occur if k<=6; Hessian non-positive definite because lambda gradient = 0
# s00 <- gam(wear~s(size,k=8,bs='ts'),data=engine,fit=FALSE)

#Input
modMat <- s0$X[,-1] #Model matrix (no intercept)
penMat <- s0$smooth[[1]]$S[[1]] #Penalty matrix (sparse)
penMat <- as(penMat,'sparseMatrix')
penDim <- nrow(penMat) #Dimensions of penalty matrix
# 
# s0$X[,-1]-s00$X[,-1] #No difference in model matrices b/w tp and ts
# s0$smooth[[1]]$S[[1]]-s00$smooth[[1]]$S[[1]]
# 
# s0$smooth[[1]]
# s00$smooth[[1]]


#For making predictions:
predSize <- with(engine,seq(min(size),max(size),length.out=50)) #Size of engines to predict for
predModMat <- PredictMat(s0$smooth[[1]],data = data.frame(size=predSize)) #New model matrix to predict with (no intercept)

#Input for model
datList <- list(Y=engine$wear, X=modMat, S1=penMat, S1dim=penDim, newX=predModMat, original=1) #Data
parList <- list(b0=0, smoothCoefs=rep(0,ncol(modMat)), log_lambda=0, log_sigma=0) #Parameters

compile("gam1.cpp")
beep(1)
dyn.load(dynlib("gam1"))
obj <- MakeADFun(data = datList, parameters = parList, random=c('b0','smoothCoefs'), DLL = "gam1")

obj$par
obj$fn() #Works
opt <- nlminb(obj$par,obj$fn,obj$gr) #Runs into NaNs
rep <- sdreport(obj) #Estimates and SEs
rep #Looks OK
rep$value
rep$par.fixed
rep$par.random

data.frame(name=names(obj$par), Est=opt$par, final_gradient=as.vector(obj$gr(opt$par)))

#Get coefs from model
m1Tmb <- rbind(with(rep,data.frame(coef=names(par.fixed),est=par.fixed,sd=sqrt(diag(cov.fixed)),row.names=NULL)),
      with(rep,data.frame(coef=names(par.random),est=par.random,sd=sqrt(diag.cov.random))))

#Compare to GAM fit
s1 <- gam(G=s0,method='REML') #Get s0 and actually fit it this time
summary(s1)

#Intercepts are identical, smoothing coefs are not
data.frame(gamCoefs=coef(s1),tmbCoefs=m1Tmb$est[grepl('(b0|smoothCoefs)',m1Tmb$coef)]) %>% 
  rownames_to_column(var = 'coef') %>% mutate(coef=ifelse(grepl('Intercept',coef),'Intercept','Smoother\nCoefs')) %>% 
  ggplot()+geom_abline(intercept = 0, slope=1)+
  geom_point(aes(x=gamCoefs,y=tmbCoefs,col=coef))+
  scale_colour_manual(values=c('red','black'))+
  labs(x='GAM Coefs',y='TMB Coefs')

gamPred <- predict(s1,newdata=data.frame(size=predSize),se.fit=TRUE) %>% bind_cols() %>% 
  mutate(upr=fit+1.96*se.fit,lwr=fit-1.96*se.fit)
tmbPred <- data.frame(fit=rep$value,se.fit=rep$sd) %>% 
  mutate(upr=fit+1.96*se.fit,lwr=fit-1.96*se.fit)

(p <- bind_rows(gamPred,tmbPred) %>% mutate(mod=rep(c('GAM','TMB'),c(nrow(gamPred),nrow(tmbPred)))) %>% 
    mutate(size=rep(predSize,2)) %>% 
    ggplot(aes(x=size,y=fit))+
    geom_point(data=engine,aes(x=size,y=wear))+
    geom_ribbon(aes(ymax=upr,ymin=lwr,fill=mod),alpha=0.3)+
    geom_line(aes(col=mod),size=1)+
    scale_fill_manual(values=c('red','blue'))+
    scale_colour_manual(values=c('red','blue'))+
    labs(x='Size',y='Wear',col='Model',fill='Model'))
ggsave('./gamTMB_tp.png',p)

# Cubic regression spline - as above ---------------

s0 <- gam(wear~s(size,k=8,bs='cr'),data=engine,fit=FALSE) 

#Input
modMat <- s0$X[,-1] #Model matrix (no intercept)
penMat <- s0$smooth[[1]]$S[[1]] #Penalty matrix (sparse)
penMat <- as(penMat,'sparseMatrix')
penDim <- nrow(penMat) #Dimensions of penalty matrix

#For making predictions:
predSize <- with(engine,seq(min(size),max(size),length.out=50)) #Size of engines to predict for
predModMat <- PredictMat(s0$smooth[[1]],data = data.frame(size=predSize)) #New model matrix to predict with (no intercept)

#Input for model
datList <- list(Y=engine$wear, X=modMat, S1=penMat, S1dim=penDim, newX=predModMat) #Data
parList <- list(b0=0, smoothCoefs=rep(0,ncol(modMat)), log_lambda=0, log_sigma=0) #Parameters

compile("gam1.cpp")
beep(1)
dyn.load(dynlib("gam1"))
obj <- MakeADFun(data = datList, parameters = parList, random=c('b0','smoothCoefs'), DLL = "gam1")

obj$par
obj$fn() #Works
opt <- nlminb(obj$par,obj$fn,obj$gr) #Runs into NaNs
rep <- sdreport(obj) #Estimates and SEs
rep #Looks OK
rep$value
rep$par.fixed
rep$par.random

#Compare to GAM fit
s1 <- gam(G=s0,method='REML') #Get s0 and actually fit it this time
summary(s1)

#Intercepts are identical, smoothing coefs are not
m1Tmb <- rbind(with(rep,data.frame(coef=names(par.fixed),est=par.fixed,sd=sqrt(diag(cov.fixed)),row.names=NULL)),
               with(rep,data.frame(coef=names(par.random),est=par.random,sd=sqrt(diag.cov.random))))

data.frame(gamCoefs=coef(s1),tmbCoefs=m1Tmb$est[grepl('(b0|smoothCoefs)',m1Tmb$coef)]) %>% 
  rownames_to_column(var = 'coef') %>% mutate(coef=ifelse(grepl('Intercept',coef),'Intercept','Smoother\nCoefs')) %>% 
  ggplot()+geom_abline(intercept = 0, slope=1)+
  geom_point(aes(x=gamCoefs,y=tmbCoefs,col=coef))+
  scale_colour_manual(values=c('red','black'))+
  labs(x='GAM Coefs',y='TMB Coefs')

gamPred <- predict(s1,newdata=data.frame(size=predSize),se.fit=TRUE) %>% bind_cols() %>% 
  mutate(upr=fit+1.96*se.fit,lwr=fit-1.96*se.fit)
tmbPred <- data.frame(fit=rep$value,se.fit=rep$sd) %>% 
  mutate(upr=fit+1.96*se.fit,lwr=fit-1.96*se.fit)

(p <- bind_rows(gamPred,tmbPred) %>% mutate(mod=rep(c('GAM','TMB'),c(nrow(gamPred),nrow(tmbPred)))) %>% 
  mutate(size=rep(predSize,2)) %>% 
  ggplot(aes(x=size,y=fit))+
  geom_point(data=engine,aes(x=size,y=wear))+
  geom_ribbon(aes(ymax=upr,ymin=lwr,fill=mod),alpha=0.3)+
  geom_line(aes(col=mod),size=1)+
  scale_fill_manual(values=c('red','blue'))+
  scale_colour_manual(values=c('red','blue'))+
  labs(x='Size',y='Wear',col='Model',fill='Model'))
ggsave('./gamTMB_cr.png',p)
