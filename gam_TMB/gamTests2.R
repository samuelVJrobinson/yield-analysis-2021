#GAM testing using TMB
# Thin plate spline - 2 smoothers 

setwd("~/Documents/yield-analysis-2021/gam_TMB")

library(Matrix) #Use sparse matrices
library(TMB)
library(mgcv)
library(tidyverse)
theme_set(theme_classic())

s0 <- gam(Volume~s(Girth,k=6,bs='tp')+s(Height,k=6,bs='tp'),data=trees) 

#Input
modMat <- model.matrix(s0)[,-1] #Model matrix (no intercept)
penMat1 <- s0$smooth[[1]]$S[[1]] #Penalty matrices (sparse)
penMat1 <- as(penMat1,'sparseMatrix')
penMat2 <- s0$smooth[[2]]$S[[1]] 
penMat2 <- as(penMat2,'sparseMatrix')
penDim1 <- nrow(penMat1) #Dimensions of penalty matrices
penDim2 <- nrow(penMat2) 

#For making predictions:
predSize <- with(trees,rbind(expand.grid(Girth=mean(Girth),Height=seq(min(Height),max(Height),1)),
                 expand.grid(Girth=seq(min(Girth),max(Girth),0.5),Height=mean(Height))))
predMat1 <- PredictMat(s0$smooth[[1]],data = predSize) 
predMat2 <- PredictMat(s0$smooth[[2]],data = predSize) 
predModMat <- cbind(predMat1,predMat2) #New model matrix to predict with (no intercept)

#Input for model
datList <- list(Y=trees$Volume, X=modMat, 
                S1=penMat1, S2=penMat2, S1dim=penDim1, S2dim=penDim2, newX=predModMat) #Data
parList <- list(b0=coef(s0)[1], smoothCoefs=coef(s0)[-1], log_lambda1=1, log_lambda2=1, log_sigma=1) #Parameters

compile("gam2.cpp")
dyn.load(dynlib("gam2"))
obj <- MakeADFun(data = datList, parameters = parList, random=c('b0','smoothCoefs'), DLL = "gam2")

obj$par
obj$fn() #Works
opt <- nlminb(obj$par,obj$fn,obj$gr)
rep <- sdreport(obj) #Estimates and SEs
rep #Looks OK
rep$value