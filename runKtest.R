#SCRIPT TO RUN TEST OF K-PARAMETER VALUES BETWEEN MODELS, TO GAUGE SENSITIVITY OF RESULTS


# Setup --------------------------------------------------------------

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

# datSource$modelComplete <- 
  file.exists(datSource$modelPath)
  

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

# Run models -------------------------------

library(parallel)

for(kp in 1:length(kPars)){
  outputDir <- paste0('./Data/kTestResults/',names(kPars)[kp])
  # runModI(i=1,dS=datSource,kPar=kPars[[kp]],modelCheckDir=outputDir,resultsDir=outputDir) #Test on 1 field
  cluster <- makeCluster(15) #10 procs max - uses about 90% of memory
  parLapply(cl=cluster,1:nrow(datSource),runModI,
            dS=datSource,kPar=kPars[[kp]],modelCheckDir=outputDir,resultsDir=outputDir) 
  stopCluster(cluster)
}

# Get results -------------------------------

#K-test results
resultPaths <- dir(paste0('./Data/kTestResults'),pattern='*.Rdata',recursive=TRUE,full.names=TRUE)
resultNames <- sapply(strsplit(resultPaths,'/'),function(x){l <- length(x); paste(x[(l-1):l],collapse='__')}) %>% gsub(' modList.Rdata','',.)

#"Original" results
resultPaths2 <- paste0('./Figures/ModelCheck/',unique(sapply(strsplit(resultNames,'__'),function(x) x[2])),' modList.Rdata')
resultNames2 <- sapply(strsplit(resultPaths2,'/'),function(x){l <- length(x); x[l]}) %>% gsub(' modList.Rdata','',.) %>% paste0('original__',.)

resultPaths <- c(resultPaths,resultPaths2); resultNames <- c(resultNames,resultNames2) #Join
rm(resultPaths2,resultNames2) #Cleanup
file.exists(resultPaths)


#MISMATCH BETWEEN ORIGINAL AND KTEST DIST RANGES, ONLY AT TRENT CLARK'S FIELDS

#Trent Clark
load("./Figures/ModelCheck/Trent Clark 2019 W 34 modList.Rdata") #Original 
modList$distRange
load("./Data/kTestResults/lwrSpace/Trent Clark 2019 W 34 modList.Rdata") #Ktest 
modList$distRange

# #Alvin French
# load("./Figures/ModelCheck/Alvin French 2020 Al_Jr modList.Rdata") #Original
# load("./Data/kTestResults/lwrSpace/Alvin French 2020 Al_Jr modList.Rdata") #Ktest

results <- lapply(resultPaths,function(p){
  load(p)
  getSmooths(smoothLabel='s(dist)',modList,xvals=seq(modList$distRange[1],modList$distRange[2],by=10))}) %>% 
  set_names(resultNames) %>% bind_rows(.id='name') %>% 
  separate(name,c('parSet','field'),sep='__') %>% mutate(dist=round(dist))

# temp <- 
  left_join(filter(results,parSet!='original'),filter(results,parSet=='original'),by=c('field','dist')) %>% 
    filter(is.na(pred.y))
    # arrange(dist,field) 
  
filter(results,field=='Trent Clark 2019 W 34')  %>% group_by(parSet) %>% slice(1:5) %>% data.frame()
