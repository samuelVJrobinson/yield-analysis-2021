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
modListPaths <- dir(paste0('./Data/kTestResults'),pattern='*.Rdata',recursive=TRUE,full.names=TRUE)
resultNames <- sapply(strsplit(modListPaths,'/'),function(x){l <- length(x); paste(x[(l-1):l],collapse='__')}) %>% gsub(' modList.Rdata','',.)

#"Original" results
modListPaths2 <- paste0('./Figures/ModelCheck/',unique(sapply(strsplit(resultNames,'__'),function(x) x[2])),' modList.Rdata')
resultNames2 <- sapply(strsplit(modListPaths2,'/'),function(x){l <- length(x); x[l]}) %>% 
  gsub(' modList.Rdata','',.) %>% paste0('original__',.)

modListPaths <- c(modListPaths,modListPaths2); resultNames <- c(resultNames,resultNames2) #Join
rm(modListPaths2,resultNames2) #Cleanup
file.exists(modListPaths)

results <- lapply(modListPaths,function(p){
  load(p)
  getSmooths(smoothLabel='s(dist)',modList,xvals=seq(modList$distRange[1],modList$distRange[2],by=10))}) %>% 
  set_names(resultNames) %>% bind_rows(.id='name') %>% separate(name,c('parSet','field'),sep='__') %>% 
  mutate(dist=round(dist)) 

resultsOrig <- filter(results,parSet=='original')
resultsNew <- filter(results,parSet!='original')
# head(resultsOrig)
# head(resultsNew)
results <- left_join(resultsNew,resultsOrig,by=c('field','dist'),suffix=c('new','orig')) %>% 
  select(-parSetorig) %>% mutate(diff=predorig-prednew) %>% 
  mutate(parSetnew=factor(parSetnew,labels=c('Less k: space','Less k: time','More k: space','More k: time')))
  # filter(field=='Alvin French 2020 Al_Jr',dist==5)

p1 <- ggplot(results)+geom_point(aes(x=predorig,y=prednew))+
  facet_wrap(~parSetnew)+
  labs(x='Original effect',y='Updated effect')+
  geom_abline(intercept=0,slope=1,col='red',linetype='dashed')

p2 <- ggplot(results)+geom_line(aes(x=dist,y=diff,group=field))+
  facet_wrap(~parSetnew)+
  labs(x='Distance from boundary',y='Difference in effects')

p3 <- ggplot(results)+
  geom_line(aes(x=dist,y=predorig,group=field),col='black',alpha=0.5)+
  geom_line(aes(x=dist,y=prednew,group=field),col='red',alpha=0.5)+
  facet_wrap(~parSetnew)+
  labs(x='Distance from boundary',y='Boundary effect')

(p <- ggarrange(p3,p1,p2,ncol=3))
ggsave('./Data/kTestResults/kTestResults.png',p,height=4,width=12,dpi=250)  

#How long did models take to run ----------------

resultPaths <- gsub(' modList.Rdata',' results.txt',modListPaths)
modelInfo <- sapply(resultPaths,getModelInfo)

getModelInfo(resultPaths[50])


for(i in 1:length(resultPaths)){
  getModelInfo(resultPaths[i])
}


'Dean Hubbard 2016 E_21_11_25 results.txt'