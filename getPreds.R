#Script to get parameter estimates from saved GAM models, and fit meta-model that incorporates uncertainty at the lower levels
# Written by Sam Robinson for Lan Nguyen, Summer 2021

# Load libraries and functions --------------------------

library(mgcv)
library(tidyverse)
theme_set(theme_classic())

# Convenience function to get parameters from models
# d = data frame of new parameters to use. Useful if you want predictions of 1 thing at a time (e.g. Distance at specific Lat/Lon values)
#   Can also be: 'original' uses original data, 'distance' uses unique distances and median row/col values.
getPars <- function(path,d='distance'){
  load(path) #Load file into environment
  beta <- coefficients(mod) #Coefficients
  sigma <- sqrt(diag(vcov(mod))) #Standard errors
  
  #Range of predictors in original model
  predRange <- cbind(summary(with(mod$smooth[[1]],data.frame(Dist=as.vector(Xu)+as.vector(shift)))), 
                     summary(with(mod$smooth[[2]],data.frame(Row=Xu[,1]+shift[1],Col=Xu[,2]+shift[2]))))
  predNames <- gsub('\\s*(:|:\\s).*','',predRange[,1])
  predRange <- apply(predRange,2,function(x) as.numeric(gsub('.+(:|:\\s)','',x)))
  
  nUnique <- c(with(mod$smooth[[1]],length(unique(Xu))), #Number of unique values
               with(mod$smooth[[2]],apply(Xu,2,function(x) length(unique(x)))))
  predRange <- rbind(predRange,nUnique)
  rownames(predRange) <- c(predNames,'Nunique')
  
  if(is.null(d)){
    mm <- model.matrix(mod)
    d <- NULL
    } else if(is.data.frame(d)) {
      mm <- predict(mod,type='lpmatrix',newdata=d) #Model matrix from new dataframe
    } else if(d=='distance'){
      d <- data.frame(Dist=with(mod$smooth[[1]],sort(unique(as.vector(Xu)+as.vector(shift)))), #Unique distances
                      Row=with(mod$smooth[[2]],median(Xu[,1]+shift[1])), #Median row
                      Col=with(mod$smooth[[2]],median(Xu[,2]+shift[2]))) #Median column
      mm <- predict(mod,type='lpmatrix',newdata=d) #Model matrix from new dataframe
      mm[,!grepl('Dist',colnames(mm))] <- 0 #Set any non-distance column to 0
    }
  
  rm(mod); gc(); #Run garbage collector
  ret <- list(beta=beta,sigma=sigma,mm=mm,predRange=predRange,dat=d) #List of objects to return
  return(ret)
}


#Make distance smoother from model parameters
#if sim = TRUE, draws a single value from posterior distribution
distSmoother <- function(modPar,sim=FALSE){ 
  parInd <- attr(modPar$mm,'lpi') #Parameter index (are parameters for mean or variance?)
  
  if(sim){ #Simulate new coefficients using mean/SE
    betaMu <- rnorm(length(parInd[[1]]),modPar$beta[parInd[[1]]],modPar$sigma[parInd[[1]]])
    betaSigma <- rnorm(length(parInd[[2]]),modPar$beta[parInd[[2]]],modPar$sigma[parInd[[2]]])
  } else {
    betaMu <- modPar$beta[parInd[[1]]] #Use estimate
    betaSigma <- modPar$beta[parInd[[2]]]
  }
  #Model matrices
  XMu <- modPar$mm[,parInd[[1]]]
  XSigma <- modPar$mm[,parInd[[2]]]
  
  #Predictions
  data.frame(Dist=modPar$dat$Dist,
             predMean=XMu %*% betaMu,
             predSD= XSigma %*% betaSigma)
}

# Run models -------------------------

#Get paths for models from directory storing all models
paths <- dir('/home/rsamuel/Downloads/GAMs',pattern='*.RData',full.names = TRUE)

#Get parameters from all models (takes a few seconds for 10 models on my computer)
modParList <- lapply(paths,getPars) 

#Plot of mean smoothers for each field with overall smoothers
#Problem: this treats smoothers like data, when in reality they have uncertainty in them
predSmooths <- lapply(modParList,distSmoother) %>% 
  bind_rows(.id='Field') %>% rename_with(.fn=function(x) gsub('pred','',x),contains('pred')) %>% 
  mutate(SD=exp(SD)) %>% #Back-transform
  pivot_longer(Mean:SD,names_to = 'var')

ggplot(predSmooths,aes(x=Dist,y=value))+geom_line(aes(group=Field))+
  geom_smooth(method='gam',formula=y ~ s(x),col='red',fill='red',alpha=0.3)+
  facet_wrap(~var,ncol=1,scales='free_y')+
  labs(x='Distance from Boundary(?)',y='Parameter')

#Solution: 1) simulate smoothers for each field, 2) fit a GAM to these simulated smoothers (I call this a meta-model), 
# 3) Repeat many times, 4) Compare predictions of meta-models


makeSimSmooth <- function(bs){ #Argument does nothing
  require(mgcv)
  simSmooths <- do.call('rbind',lapply(modParList,distSmoother,sim=TRUE)) #Draws a set of random smooths
  m1 <- gam(predMean~s(Dist,k=20),data=simSmooths) #Fits a model of the mean
  m2 <- gam(predSD~s(Dist,k=20),data=simSmooths) #Fits a model of the log(SD)
  nd <- data.frame(Dist=sort(unique(simSmooths$Dist))) #Get unique distances
  nd$Mean <- predict(m1,newdata=nd) #mean model results
  nd$SD <- predict(m2,newdata=nd) #log(SD) model results
  return(nd)
}
  
#This takes a long time, but can be done in parallel (see below)
Nrep <- 1000
a <- lapply(1:Nrep,makeSimSmooth)

# library(parallel) #Parallel version
# Nrep <- 1000
# detectCores() #How many cores are on your machine? Keep 1 core unused for background processes
# cluster <- makeCluster(5) #Make connection to CPU clusters
# clusterExport(cl = cluster,varlist=c('distSmoother','modParList'), envir = .GlobalEnv) #Export function to clusters
# a <- parLapply(cl=cluster,1:Nrep,fun=makeSimSmooth) #Run makeSimSmooth in parallel
# stopCluster(cluster) #Stop cluster connection

#Get results and plot them
#Red line is the median smoother, pink areas are 90% credible intervals (i.e. 90% of the simulations were between these lines)
#Grey lines are original field level smoothers
a %>% bind_rows(.id = 'sim') %>% 
  mutate(across(contains('SD'),exp)) %>% #Back-transform SD 
  group_by(Dist) %>% #Group by distances
  summarize(across(c(Mean,SD),list(med=median,upr=~quantile(.x,0.95),lwr=~quantile(.x,0.05)))) %>% #Median, upper, and lower quantiles
  ungroup() %>% pivot_longer(-Dist) %>% separate(name,c('var','range'),sep='_') %>%  #Reshape data
  pivot_wider(names_from=range,values_from=value) %>% 
  ggplot(aes(x=Dist))+ #Make plot
  geom_ribbon(aes(ymax=upr,ymin=lwr),alpha=0.3,fill='red')+
  geom_line(data=predSmooths,aes(x=Dist,y=value,group=Field),alpha=0.3)+ #Original field-level smoothers
  geom_line(aes(y=med),col='red')+
  facet_wrap(~var,ncol=1,scales = 'free_y')+
  labs(x='Distance from Boundary(?)',y='Parameter')
