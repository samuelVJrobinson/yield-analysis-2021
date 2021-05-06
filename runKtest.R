#SCRIPT TO RUN TEST OF K-PARAMETER VALUES BETWEEN MODELS, TO GAUGE SENSITIVITY OF RESULTS


# Run models --------------------------------------------------------------

library(tidyverse)
theme_set(theme_bw())
library(ggpubr)
library(mgcv)
library(sf)
library(beepr)
source('helperFunctions.R')

set.seed(1)

#Get paths to data
datSource <- read.csv('./Data/datSource.csv') %>% #Read previous datasource file
  # select(dataPath:boundaryPath,use) %>% 
  filter(use) %>% #Sample 3 fields from each grower
  group_by(grower) %>% 
  slice_sample(n=3) %>% 
  mutate(boundaryComplete=file.exists(boundaryPath)) %>% #Has boundary been made already?
  mutate(modelComplete=FALSE) #Has model1 already been run? Set to false
  

#List of k-parameters to use
kPars<-list(
  lwrSpace=c(12,30,60,12,30,60),
  uprSpace=c(12,120,60,12,120,60),
  lwrTime=c(12,60,30,12,60,30),
  uprTime=c(12,60,120,12,60,120)
  # original=c(12,60,60,12,60,60),
    )

# #Make directories
# for(kp in 2:length(kPars)) dir.create(paste0('./Data/kTestResults/',names(kPars)[kp]))

library(parallel)

for(kp in 2:length(kPars)){
  outputDir <- paste0('./Data/kTestResults/',names(kPars)[kp])
  # runModI(i=1,dS=datSource,kPar=kPars[[kp]],modelCheckDir=outputDir,resultsDir=outputDir) #Test on 1 field
  cluster <- makeCluster(15) #10 procs max - uses about 90% of memory
  parLapply(cl=cluster,1:nrow(datSource),runModI,
            dS=datSource,kPar=kPars[[kp]],modelCheckDir=outputDir,resultsDir=outputDir) 
  stopCluster(cluster)
}







