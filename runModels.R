#Run non-isotropic GAMs on multiple field years

# Load everything ---------------------------------------------------------

library(tidyverse)
theme_set(theme_bw())
library(ggpubr)
library(mgcv)
library(sf)
library(beepr)
source('helperFunctions.R')

# Get paths to data 

# #Generate new dataSource dataframe
# rootPath <- "/media/rsamuel/Storage/geoData/Rasters/yieldData/csv files" #Path to csv files
# datSource <- data.frame(dataPath=dir(rootPath,pattern=".csv",recursive=TRUE,full.names = TRUE)) %>%
#   mutate(path=gsub('/media/rsamuel/Storage/geoData/Rasters/yieldData/csv files/','',dataPath)) %>%
#   separate(path,c('grower','field','year'),sep=" ",remove=FALSE) %>% select(-path) %>%
#   mutate(year=gsub('\\.csv','',year)) %>%
#   unite(filename,c(grower:year),sep=' ',remove=FALSE) %>%
#   mutate(boundaryPath=paste0('./Figures/FieldBoundaries/',filename,' boundary.shp')) %>% #Paths to boundary shapefiles
#   mutate(boundaryPath2=paste0('./Figures/FieldBoundarySegments/',filename,' boundary.shp')) %>%
#   mutate(modelPath=paste0('./Figures/ModelCheck/',filename,' modList.Rdata')) %>% #Paths to saved model files
#   mutate(modelPath2=paste0('./Figures/ModelCheck2/',filename,' modList.Rdata'))
# 
# #Read 10000 lines from each csv, and get counts of products - takes about 20 seconds
# cropType <- sapply(datSource$dataPath,function(path){
#   table(read.csv(path,nrows=10000,fileEncoding='latin1')$Product)}) %>% 
#   set_names(datSource$filename) %>% 
#   sapply(.,function(x) names(x)[1])
# 
# datSource$crop <- case_when(grepl('(CANOLA|InVigor|Dekalb 74-54|Dekalb 75-42|Brett Young 6056|Dupont D3152|VT 9562)',cropType) ~ 'Canola',
#                             grepl('(WHEAT|Wheat|CPS AC Foremost|CWRS|CWSW)',cropType) ~ 'Wheat',
#                             grepl('Peas',cropType) ~ 'Peas',
#                             grepl('(BARLEY|Barley)',cropType) ~ 'Barley',
#                             grepl('Beans',cropType) ~ 'Beans',
#                             grepl('Oats',cropType) ~ 'Oats',
#                             grepl('Rye',cropType) ~ 'Rye',
#                             TRUE ~ cropType) 
# write.csv(datSource,'./Data/datSource_new.csv',row.names = FALSE)

datSource <- read.csv('./Data/datSource.csv') %>% #Read previous datasource file
  mutate(boundaryComplete=file.exists(boundaryPath)) %>% #Has boundary been made already?
  mutate(modelComplete=file.exists(modelPath)) %>% #Has model1 already been run?
  mutate(modelComplete2=file.exists(modelPath2)) #Has model2 already been run?

#Get summary results from models
datSource %>% pull(use) %>% sum()
datSource %>% count(grower,year) #How many fields per grower per year
datSource %>% group_by(grower) %>% summarize(nLoc=length(unique(field)),nYears=length(unique(year)),n=n()) #Number of distinct fields 
datSource %>% count(crop) #Enough to do separate analyses on canola and wheat, maybe peas

# #First set of models - no boundary types
# resultPaths <- gsub(' modList.Rdata',' results.txt',datSource$modelPath)[datSource$use]
# modInfo <- lapply(resultPaths,getModelInfo)
# sapply(modInfo,function(x) as.numeric(x$timeTaken,units='hours')) %>% summary() #Hours taken
# sapply(modInfo,function(x) x$percDevianceExplained) %>% summary() #Explained deviance

# #Second set of models - boundary types included
# resultPaths <- gsub(' modList.Rdata',' results.txt',datSource$modelPath2)[datSource$use]
# modInfo <- lapply(resultPaths,getModelInfo)
# sum(!file.exists(datSource$modelPath2[datSource$use])) #24 models not completed
# sapply(modInfo,function(x) as.numeric(x$timeTaken,units='hours')) %>% summary() #Hours taken
# sapply(modInfo,function(x) x$percDevianceExplained) %>% summary() #Explained deviance

# Run first set of models - no boundary type --------------------------------------------------------------

# a <- Sys.time() #Test
# runModI(51,dS=datSource)
# Sys.time()-a
# beep(1)

library(parallel)
cluster <- makeCluster(15) #10 procs max - uses about 90% of memory
parLapply(cl=cluster,1:nrow(datSource),runModI,dS=datSource) #Takes about 42 mins for running 10 procs. Some seem to take longer than others (weird shaped fields?)
beep(1)
stopCluster(cluster)
Sys.time() #Takes ~ 7 hrs

#Need to identify boundary types at each field

# Get smoother info from first set of models -------------------------------------------

# #Estimate at field 1
# est <- getPreds(paste0('./Figures/ModelCheck/',datSource$filename[1],' modList.Rdata'),margInt=c(FALSE,TRUE,TRUE),samp=FALSE)$distDat
# #Samples around estimate at field 1
# Nsamp <- 100 #Number of samples
# samp <- lapply(1:Nsamp,
#           function(x) getPreds(paste0('./Figures/ModelCheck/',datSource$filename[1],' modList.Rdata'),
#                                  margInt=c(FALSE,TRUE,TRUE),samp=TRUE)$distDat) %>% 
#   do.call('rbind',.) %>% mutate(N=rep(1:length(unique(dist)),each=Nsamp))
# ggplot()+geom_line(data=test,aes(x=dist,y=mean,group=N),alpha=0.3)+geom_line(data=est,aes(x=dist,y=mean),col='red',size=3)

#Get smoother from each field
#Takes about 10 seconds
getFiles <- datSource$filename[datSource$use]
allSmooths <- lapply(paste0('./Figures/ModelCheck/',getFiles,' modList.Rdata'),
                     getPreds,margInt=c(FALSE,TRUE,TRUE))
names(allSmooths) <- gsub(' ','-',getFiles)
names(allSmooths[[1]]) #Variables to get from allSmooths

#Get df of predictions for each variable
allEff <- lapply(names(allSmooths[[1]]),function(y){
  lapply(allSmooths,function(x) x[[y]]) %>% 
    do.call('rbind',.) %>% 
    rownames_to_column('field') %>% 
    mutate(field=gsub('\\.\\d{1,3}','',field))
})

#Plot with individual field-level smoothers

tHa2buAc <- 17.0340 #tonnes per hectare to bushels per acre (https://www.agrimoney.com/calculators/calculators)
ylimVals <- list(mean=c(-1.2,1.2)*tHa2buAc,sd=c(0,5)) #Y limits
ylabMean <- 'Mean Yield (bushels/acre)'
ylabSD <- 'log(SD Yield)'

p1 <- allEff[[1]] %>% ggplot(aes(x=pArea,y=mean*tHa2buAc))+geom_line(aes(group=field),alpha=0.3)+
  geom_smooth(method='lm',formula=y~log2(x),col='blue',se=FALSE,n=500)+
  labs(x='Polygon Area',y=ylabMean)+
  coord_cartesian(xlim = c(0,200))
p2 <- allEff[[1]] %>% ggplot(aes(x=pArea,y=logSD*tHa2buAc))+geom_line(aes(group=field),alpha=0.3)+
  geom_smooth(method='lm',formula=y~log(x),col='blue',se=FALSE,n=500)+
  labs(x='Polygon Area',y=ylabSD)+coord_cartesian(xlim = c(0,200))

p3 <- allEff[[2]] %>% ggplot(aes(x=dist,y=mean*tHa2buAc))+geom_line(aes(group=field),alpha=0.3)+
  labs(x='Boundary Distance',y=ylabMean)+
  geom_hline(yintercept = 0,linetype='dashed',col='red')+
  geom_smooth(method='gam',formula=y~s(x,k=10),col='blue',se=TRUE)+
  coord_cartesian(xlim = c(0,400),ylim=ylimVals$mean)
p4 <- allEff[[2]] %>% ggplot(aes(x=dist,y=log(exp(logSD)*tHa2buAc)))+geom_line(aes(group=field),alpha=0.3)+
  labs(x='Boundary Distance',y=ylabSD)+
  # geom_hline(yintercept = 0,linetype='dashed',col='red')+
  geom_smooth(method='gam',formula=y~s(x),col='blue',se=FALSE)+
  coord_cartesian(xlim = c(0,400),ylim=ylimVals$sd)

p5 <- allEff[[3]] %>% ggplot(aes(x=r,y=mean*tHa2buAc))+geom_line(aes(group=field),alpha=0.3)+
  labs(x='Point Order',y=ylabMean)+
  # geom_hline(yintercept = 0,linetype='dashed',col='red')+
  geom_smooth(method='gam',formula=y~s(x),col='blue',se=FALSE)+
  coord_cartesian(ylim=ylimVals$mean)
p6 <- allEff[[3]] %>% ggplot(aes(x=r,y=log(exp(logSD)*tHa2buAc)))+geom_line(aes(group=field),alpha=0.3)+
  labs(x='Point Order',y=ylabSD)+
  # geom_hline(yintercept = 0,linetype='dashed',col='red')+
  geom_smooth(method='gam',formula=y~s(x),col='blue',se=FALSE)+
  coord_cartesian(ylim=ylimVals$sd)

(p <- ggarrange(p1,p3,p5,p2,p4,p6,ncol=3,nrow=2))
ggsave(paste0('./Figures/ModelSummary.png'),p,height=6,width=12,dpi=300)  

#Get smoother from each field

getFiles <- datSource$filename[datSource$use] #Field names
files <- paste0('./Figures/ModelCheck/',getFiles,' modList.Rdata') #File paths

#Function to get predicted smoother values from list of functions, fit models to these, and return results
# Used for getting 
# f = files to draw from
# s = sample from posterior
# m = marginalize across intercept
sampleSmooth <- function(f,s,m,...){
  require(mgcv)
  require(tidyverse)
  
  allSmooths <- lapply(f,getPreds,margInt=m,samp=s)
  names(allSmooths) <- gsub(' ','-',f)
  names(allSmooths[[1]]) #Variables to get from allSmooths
  
  #Get df of predictions for each variable
  allEff <- lapply(names(allSmooths[[1]]),function(y){
    lapply(allSmooths,function(x) x[[y]]) %>% 
      do.call('rbind',.) %>% 
      rownames_to_column('field') %>% 
      mutate(field=gsub('\\.\\d{1,3}','',field))
  })
  
  #Fit models of sampled estimates
  m1mean <- lm(mean~log(pArea),data=allEff[[1]])
  m1sd <- lm(logSD~log(pArea),data=allEff[[1]])
  m2mean <- gam(mean~s(dist,k=30),data=allEff[[2]])
  m2sd <- gam(logSD~s(dist,k=30),data=allEff[[2]])
  m3mean <- gam(mean~s(r,k=30),data=allEff[[3]])
  m3sd <- gam(logSD~s(r,k=30),data=allEff[[3]])
  
  #Predictions from models
  datList <- list(
    #200 locations along log-pArea line
    pArea=data.frame(pArea=exp(seq(min(log(allEff[[1]]$pArea)),max(log(allEff[[1]]$pArea)),length.out=200))) %>%
      mutate(mean=predict(m1mean,newdata=.),logSD=predict(m1sd,newdata=.)),
    #Unique rounded distance measurements
    dist=data.frame(dist=sort(unique(round(allEff[[2]]$dist)))) %>% 
      mutate(mean=predict(m2mean,newdata=.),logSD=predict(m2sd,newdata=.)), 
    #Unique rounded sequence measurements
    r=data.frame(r=round(seq(min(allEff[[3]]$r),max(allEff[[3]]$r),length.out=1000))) %>% 
      mutate(mean=predict(m3mean,newdata=.),logSD=predict(m3sd,newdata=.))
  )
  rm(list=ls()[ls()!='datList']) #Cleanup
  gc()
  return(datList)
}

temp <- sampleSmooth(f=files,s=FALSE,m=c(TRUE,TRUE,TRUE)) #Mean smoother

# Nrep <- 3
# # temp2 <- replicate(Nrep,sampleSmooth(f=files,s=TRUE),simplify=FALSE)
# debugonce(sampleSmooth)
# debugonce(getPreds)
# temp2 <- sampleSmooth(f=files,s=TRUE)
# 
# lapply(temp2,function(x) x$r) %>% set_names(paste0('s',1:Nrep)) %>% 
#   do.call('rbind',.) %>%  rownames_to_column('rep') %>%   mutate(rep=gsub('\\.\\d{1,3}','',rep)) %>% 
#   group_by(rep) %>% slice(1:3)

library(parallel) #Parallel version
Nrep <- 200
detectCores()
cluster <- makeCluster(15) 
clusterExport(cl = cluster,varlist='getPreds', envir = .GlobalEnv) #Export function to clusters
tempSamp <- parLapply(cl=cluster,1:Nrep,fun=sampleSmooth,f=files,s=TRUE,m=c(TRUE,TRUE,TRUE))
beep(1)
stopCluster(cluster)

#Get range of variability for models 

ylabMean <- 'Mean Yield (bushels/acre)'
ylabSD <- 'log(SD Yield)'

#Low variability in pArea models
p1 <- lapply(tempSamp,function(x) x$pArea) %>% set_names(paste0('s',1:Nrep)) %>% 
  do.call('rbind',.) %>%  rownames_to_column('rep') %>%   mutate(rep=gsub('\\.\\d{1,3}','',rep)) %>% 
  ggplot(aes(x=pArea,y=mean))+geom_line(aes(group=rep),alpha=0.3)+
  geom_line(data=temp$pArea,col='blue')+labs(x='Polygon Area',y=ylabMean)
p2 <- lapply(tempSamp,function(x) x$pArea) %>% set_names(paste0('s',1:Nrep)) %>% 
  do.call('rbind',.) %>%  rownames_to_column('rep') %>%   mutate(rep=gsub('\\.\\d{1,3}','',rep)) %>% 
  ggplot(aes(x=pArea,y=logSD))+geom_line(aes(group=rep),alpha=0.3)+
  geom_line(data=temp$pArea,col='blue')+labs(x='Polygon Area',y=ylabSD)


#Higher variability in dist model
p3 <- lapply(tempSamp,function(x) x$dist) %>% set_names(paste0('s',1:Nrep)) %>% 
  do.call('rbind',.) %>%  rownames_to_column('rep') %>%   mutate(rep=gsub('\\.\\d{1,3}','',rep)) %>% 
  ggplot(aes(x=dist,y=mean))+geom_line(aes(group=rep),alpha=0.3)+
  geom_line(data=temp$dist,col='blue',size=2)+labs(x='Boundary Distance',y=ylabMean)
p4 <- lapply(tempSamp,function(x) x$dist) %>% set_names(paste0('s',1:Nrep)) %>% 
  do.call('rbind',.) %>%  rownames_to_column('rep') %>%   mutate(rep=gsub('\\.\\d{1,3}','',rep)) %>% 
  ggplot(aes(x=dist,y=logSD))+
  geom_line(aes(group=rep),alpha=0.3)+
  geom_line(data=temp$dist,col='blue',size=2)+labs(x='Boundary Distance',y=ylabSD)

#Low variability in r model
p5 <- lapply(tempSamp,function(x) x$r) %>% set_names(paste0('s',1:Nrep)) %>% 
  do.call('rbind',.) %>%  rownames_to_column('rep') %>%   mutate(rep=gsub('\\.\\d{1,3}','',rep)) %>% 
  ggplot(aes(x=r,y=mean))+
  geom_line(aes(group=rep),alpha=0.3)+
  geom_line(data=temp$r,col='blue',size=2)+labs(x='Point Order',y=ylabMean)
p6 <- lapply(tempSamp,function(x) x$r) %>% set_names(paste0('s',1:Nrep)) %>% 
  do.call('rbind',.) %>%  rownames_to_column('rep') %>%   mutate(rep=gsub('\\.\\d{1,3}','',rep)) %>% 
  ggplot(aes(x=r,y=logSD))+
  geom_line(aes(group=rep),alpha=0.3)+
  geom_line(data=temp$r,col='blue',size=2)+labs(x='Point Order',y=ylabSD)
(p <- ggarrange(p1,p3,p5,p2,p4,p6,ncol=3,nrow=2))
ggsave(paste0('./Figures/ModelSummary1a.png'),p,height=6,width=12,dpi=300)  

# Run second set of models - boundary type included ---------------------------

# a <- Sys.time()
# runModII(295,dS=datSource) 
# Sys.time()-a
# # beep(1)
# debugonce(runModII)

library(parallel)
detectCores()
cluster <- makeCluster(15) #10 procs uses about 30% of memory - could probably max it out
a <- Sys.time()
parLapply(cl=cluster,1:nrow(datSource),runModII,dS=datSource,useClosest=TRUE) #All models
beep(1)
stopCluster(cluster)
Sys.time()-a #Takes ~ 10 hrs for 50 models at 15 procs
#Takes 2.364893 days to do all fields at 15 procs

# Get smoother info from second set of models ---------------------------

# #Estimate at single field
# est <- getPreds(paste0('./Figures/ModelCheck2/',datSource$filename[2],' modList.Rdata'),
#                 margInt=c(FALSE,TRUE,TRUE),
#                 samp=FALSE)$distDat %>%
#   bind_rows(.id='dist_type')
# 
# #Draw posterior samples at field 1
# Nsamp <- 100 #Number of samples
# # samp <- lapply(1:Nsamp,
# #           function(x){
# #             getPreds(paste0('./Figures/ModelCheck2/',datSource$filename[2],' modList.Rdata'),margInt=c(FALSE,TRUE,TRUE),samp=TRUE)$distDat %>%
# #               bind_rows(.id='dist_type')
# #           }) %>% bind_rows(.id='N')
# library(parallel)
# cluster <- makeCluster(15) #10 procs uses about 30% of memory - could probably max it out
# clusterExport(cluster,c('datSource'))
# samp <- parLapply(cl=cluster,1:Nsamp,fun=function(x){
#   require(tidyverse)
#   source('helperFunctions.R')
#   getPreds(paste0('./Figures/ModelCheck2/',datSource$filename[2],' modList.Rdata'),margInt=c(FALSE,TRUE,TRUE),samp=TRUE)$distDat %>%
#     bind_rows(.id='dist_type')}) %>% bind_rows(.id='N') 
# stopCluster(cluster)
# 
# ggplot()+ #Visualize draws
#   geom_line(data=samp,aes(x=dist,y=mean,group=N),alpha=0.3)+
#   geom_line(data=est,aes(x=dist,y=mean),col='blue',size=1)+
#   geom_hline(yintercept=0,col='red',linetype='dashed')+
#   facet_wrap(~dist_type)
# 
# est %>% group_by(dist_type) %>% slice(1)

#Get smoother from each field
getFiles <- datSource %>% filter(use,modelComplete2) #Completed models
allSmooths <- lapply(getFiles$modelPath2,getPreds,margInt=c(FALSE,TRUE,TRUE)) %>% #Takes about 20 seconds
  set_names(getFiles$filename)

#Get df of predictions for each variable
allEff <- lapply(names(allSmooths[[1]]),function(y){
  if(class(allSmooths[[1]][[y]])=='data.frame'){
    lapply(allSmooths,function(x) x[[y]]) %>% 
      bind_rows(.id='field') 
  } else {
    lapply(allSmooths,function(x) x[[y]] %>% bind_rows(.id = 'dist_type')) %>% 
      bind_rows(.id='field')
  }
}) %>% set_names(names(allSmooths[[1]]))

allEff[[2]] <- mutate(allEff[[2]],dist_type=gsub('dist\\:boundaryType','',dist_type)) #Change names of boundary types

tHa2buAc <- 17.0340 #tonnes per hectare to bushels per acre (https://www.agrimoney.com/calculators/calculators)
# ylimVals <- list(mean=c(-1.2,1.2)*tHa2buAc,sd=c(0,5)) #Y limits
ylabMean <- 'Mean Yield (bushels/acre)'
ylabSD <- 'log(SD Yield)'
alphaVal <- 0.1

p1 <- allEff[[1]] %>% ggplot(aes(x=pArea,y=mean*tHa2buAc))+geom_line(aes(group=field),alpha=alphaVal)+
  geom_smooth(method='lm',formula=y~log(x),col='blue',se=FALSE,n=500)+
  labs(x='Polygon Area',y=ylabMean)+
  coord_cartesian(xlim = c(0,200))
p2 <- allEff[[1]] %>% ggplot(aes(x=pArea,y=logSD*tHa2buAc))+geom_line(aes(group=field),alpha=alphaVal)+
  geom_smooth(method='lm',formula=y~log(x),col='blue',se=FALSE,n=500)+
  labs(x='Polygon Area',y=ylabSD)+coord_cartesian(xlim = c(0,200))

p3 <- allEff[[2]] %>% ggplot(aes(x=dist,y=mean*tHa2buAc))+geom_line(aes(group=field),alpha=alphaVal)+
  facet_wrap(~dist_type)+
  labs(x='Boundary Distance',y=ylabMean)+
  geom_hline(yintercept = 0,linetype='dashed',col='red')+
  geom_smooth(method='gam',formula=y~s(x,k=10),col='blue',se=TRUE)+
  coord_cartesian(xlim = c(0,400),ylim=c(-15,15))
p4 <- allEff[[2]] %>% ggplot(aes(x=dist,y=log(exp(logSD)*tHa2buAc)))+geom_line(aes(group=field),alpha=alphaVal)+
  facet_wrap(~dist_type)+
  labs(x='Boundary Distance',y=ylabSD)+
  # geom_hline(yintercept = 0,linetype='dashed',col='red')+
  geom_smooth(method='gam',formula=y~s(x),col='blue',se=FALSE)+
  coord_cartesian(xlim = c(0,400),ylim=c(0,5))

p5 <- allEff[[3]] %>% ggplot(aes(x=r,y=mean*tHa2buAc))+geom_line(aes(group=field),alpha=alphaVal)+
  labs(x='Point Order',y=ylabMean)+
  # geom_hline(yintercept = 0,linetype='dashed',col='red')+
  geom_smooth(method='gam',formula=y~s(x),col='blue',se=FALSE)
p6 <- allEff[[3]] %>% ggplot(aes(x=r,y=log(exp(logSD)*tHa2buAc)))+geom_line(aes(group=field),alpha=alphaVal)+
  labs(x='Point Order',y=ylabSD)+
  # geom_hline(yintercept = 0,linetype='dashed',col='red')+
  geom_smooth(method='gam',formula=y~s(x),col='blue',se=FALSE)

(p <- ggarrange(p1,p3,p5,p2,p4,p6,ncol=3,nrow=2))
ggsave(paste0('./Figures/ModelSummary2.png'),p,height=6,width=12,dpi=300)
(p <- ggarrange(p3,p4,ncol=1,nrow=2))
ggsave(paste0('./Figures/ModelSummary2a.png'),p,height=6,width=12,dpi=300)  

#Draw posterior samples from each field and fit curve at each field

# Function to sample from posterior of each prediction, amalgamate, fit meta-model, return results of meta-model
# My way of getting around the problem of random effects
samplePreds <- function(a,ds=datSource,nX=c(100,100,100)){ #a ignored, ds = datSource csv, nX = number of samples along range of X
  require(tidyverse); require(mgcv)
  source('helperFunctions.R')
  
  getFiles <- ds %>% filter(use,modelComplete2) #Completed models
  
  #Sample from posterior, get new lines for each field
  allSmooths <- lapply(getFiles$modelPath2,getPreds,margInt=c(FALSE,TRUE,TRUE),samp=TRUE) %>% set_names(getFiles$filename)
  
  #log(Area) meta-models
  pAreaDf <- lapply(allSmooths,function(x) x$pAreaDat) %>% bind_rows(.id = 'field')
  
  m1 <- lm(mean~log(pArea),data=pAreaDf) 
  m2 <- lm(logSD~log(pArea),data=pAreaDf)
  
  # pAreaDf %>% mutate(pred=predict(m1,se.fit = TRUE)$fit,se=predict(m1,se.fit = TRUE)$se) %>% 
  #   arrange(pArea) %>% ggplot(aes(x=pArea,y=mean))+geom_line(aes(group=field),alpha=0.3)+
  #   geom_ribbon(aes(ymax=pred+2*se,ymin=pred-2*se),fill='red',alpha=0.5)+
  #   geom_line(aes(y=pred),col='red')
  #   
  # pAreaDf %>% mutate(pred=predict(m2,se.fit = TRUE)$fit,se=predict(m2,se.fit = TRUE)$se) %>% 
  #   arrange(pArea) %>% ggplot(aes(x=pArea,y=logSD))+geom_line(aes(group=field),alpha=0.3)+
  #   geom_ribbon(aes(ymax=pred+2*se,ymin=pred-2*se),fill='red',alpha=0.5)+
  #   geom_line(aes(y=pred),col='red')
  
  #cover class meta-models
  
  #Get cover class names
  coverNames <- sort(unique(unlist(lapply(allSmooths, function(x) names(x$distDat)))))
  
  #Which models have which cover classes
  containsCover <- sapply(coverNames,function(z) sapply(lapply(allSmooths, function(x) names(x$distDat)), function(y) z %in% y))
  
  #Lists for cover distance models
  meanModList <- logSDModList <- vector(mode='list',length=length(coverNames)) %>% set_names(coverNames)
  
  for(i in 1:ncol(containsCover)){
    # i <- 5 #Debugging
    useThese <- unname(which(containsCover[,i]))
    cname <- colnames(containsCover)[i]
    coverDF <- lapply(allSmooths[useThese],function(x) x$distDat[[cname]]) %>% bind_rows(.id = 'field')
    meanModList[[i]] <- gam(mean~s(dist,k=20),data=coverDF)
    logSDModList[[i]] <- gam(logSD~s(dist,k=20),data=coverDF)
    # gam.check(meanModList[[i]])
    # gam.check(logSDModList[[i]])
  }
  
  #r (point order) meta-models
  rDatDf <- lapply(allSmooths,function(x) x$rDat) %>% bind_rows(.id = 'field')
  
  r1 <- gam(mean~s(r,k=50),data=rDatDf)
  r2 <- gam(logSD~s(r,k=50),data=rDatDf)
  # gam.check(r1)
  # gam.check(r2)
  # plot(r1)
  # plot(r2)
  
  #Assemble predictions
  pAreaPred <- data.frame(pArea=exp(seq(log(min(pAreaDf$pArea)),log(max(pAreaDf$pArea)),length.out=nX[1]))) %>% 
    mutate(predMean=predict(m1,newdata=.),predLogSD=predict(m2,newdata=.))
  
  coverPred <- vector(mode='list',length=length(coverNames)) %>% set_names(gsub('dist:boundaryType','',coverNames))
  for(i in 1:length(coverNames)){
    
    useThese <- unname(which(containsCover[,i]))
    cname <- colnames(containsCover)[i]
    coverDF <- lapply(allSmooths[useThese],function(x) x$distDat[[cname]]) %>% bind_rows(.id = 'field')
    
    coverPred[[i]] <- data.frame(dist=seq(min(coverDF$dist),max(coverDF$dist),length.out=nX[2])) %>% 
      mutate(predMean=predict(meanModList[[i]],newdata=.),predLogSD=predict(logSDModList[[i]],newdata=.))
  }
  
  rPred <- data.frame(r=seq(min(rDatDf$r),max(rDatDf$r),length.out=nX[3])) %>% 
    mutate(predMean=predict(r1,newdata=.),predLogSD=predict(r2,newdata=.))
  
  retList <- list(pArea=pAreaPred,coverDist=coverPred,r=rPred) #Return list of results
  
  rm(list=ls()[!ls() %in% 'retList']) #Remove everything except retList
  gc()
  return(retList) 
}

# debugonce(samplePreds)
# samplePreds() #Test

## Create samples from scratch
# Nsamp <- 100 #Number of samples
# library(parallel)
# cluster <- makeCluster(15) 
# clusterExport(cluster,c('datSource'))
# samp <- parLapply(cl=cluster,1:Nsamp,fun=samplePreds)
# #Takes about 10 mins for 100 samples
# Sys.time()
# stopCluster(cluster)
# save(samp,file='./Data/postSamples2.Rdata')

# ## Add to current samples
# Nsamp <- 700 #Number of samples
# library(parallel)
# cluster <- makeCluster(15)
# clusterExport(cluster,c('datSource'))
# a <- Sys.time()
# samp2 <- parLapply(cl=cluster,1:Nsamp,fun=samplePreds)
# #Takes about 10 mins for 100 samples
# Sys.time()-a
# stopCluster(cluster)
# load('./Data/postSamples.Rdata')
# samp <- c(samp,samp2)
# save(samp,file='./Data/postSamples.Rdata')
# rm(samp2)

# Load saved posterior samples
load('./Data/postSamples.Rdata')

tHa2buAc <- 17.0340 #tonnes per hectare to bushels per acre (https://www.agrimoney.com/calculators/calculators)
# ylimVals <- list(mean=c(-1.2,1.2)*tHa2buAc,sd=c(0,5)) #Y limits
ylabMean <- 'Mean Yield (bushels/acre)'
ylabSD <- 'log(SD Yield)'
alphaVal <- 0.1
qs <- c(0.1,0.5,0.9) #Quantiles

p1 <- lapply(samp,function(x) x$pArea) %>% bind_rows(.id = 'sample') %>% 
  group_by(pArea) %>% 
  summarise(predMean = quantile(predMean, qs), q = qs) %>% 
  ungroup() %>% mutate(q=factor(q,labels=c('lwr','med','upr'))) %>% 
  pivot_wider(names_from=q,values_from=predMean) %>% 
  ggplot(aes(x=pArea))+
  geom_line(data=allEff[[1]],aes(x=pArea,y=mean,group=field),alpha=0.3)+
  geom_ribbon(aes(ymax=upr,ymin=lwr),alpha=0.3,fill='cyan')+
  geom_line(aes(y=med),col='cyan',size=1)+
  labs(x='Polygon Area',y=ylabMean)

p2 <- lapply(samp,function(x) x$pArea) %>% bind_rows(.id = 'sample') %>% 
  group_by(pArea) %>% 
  summarise(predLogSD = quantile(predLogSD, qs), q = qs) %>% 
  ungroup() %>% mutate(q=factor(q,labels=c('lwr','med','upr'))) %>% 
  pivot_wider(names_from=q,values_from=predLogSD) %>% 
  ggplot(aes(x=pArea))+
  geom_line(data=allEff[[1]],aes(x=pArea,y=logSD,group=field),alpha=0.3)+
  geom_ribbon(aes(ymax=upr,ymin=lwr),alpha=0.3,fill='cyan')+
  geom_line(aes(y=med),col='cyan',size=1)+
  labs(x='Polygon Area',y=ylabSD)

p3 <- lapply(samp,function(x) x$coverDist %>% bind_rows(.id = 'dist_type')) %>% 
  bind_rows(.id = 'sample') %>% group_by(dist_type,dist) %>% 
  summarise(predMean = quantile(predMean, qs), q = qs) %>% 
  ungroup() %>% mutate(q=factor(q,labels=c('lwr','med','upr'))) %>% 
  pivot_wider(names_from=q,values_from=predMean) %>% 
  ggplot(aes(x=dist)) + 
  geom_line(data=allEff[[2]],aes(x=dist,y=mean,group=field),alpha=0.1) + #"Raw" smoothers
  geom_ribbon(aes(ymax=upr,ymin=lwr),alpha=0.3,fill='cyan') + 
  geom_line(aes(y=med),col='cyan') + #Meta-model
  geom_hline(yintercept = 0,col='red',linetype='dashed')+
  facet_wrap(~dist_type) + labs(x='Distance',y=ylabMean) +
  coord_cartesian(xlim=c(0,400),ylim=c(-5,5))

p4 <- lapply(samp,function(x) x$coverDist %>% bind_rows(.id = 'dist_type')) %>% 
  bind_rows(.id = 'sample') %>% group_by(dist_type,dist) %>% 
  summarise(predLogSD = quantile(predLogSD, qs), q = qs) %>% 
  ungroup() %>% mutate(q=factor(q,labels=c('lwr','med','upr'))) %>% 
  pivot_wider(names_from=q,values_from=predLogSD) %>% 
  ggplot(aes(x=dist)) + 
  geom_line(data=allEff[[2]],aes(x=dist,y=logSD,group=field),alpha=0.1) + #"Raw" smoothers
  geom_ribbon(aes(ymax=upr,ymin=lwr),alpha=0.3,fill='cyan') + geom_line(aes(y=med),col='cyan') + #Meta-model
  geom_hline(yintercept = 0,col='red',linetype='dashed')+
  facet_wrap(~dist_type) + labs(x='Distance',y=ylabSD) +
  coord_cartesian(xlim=c(0,400),ylim=c(-5,5))

p5 <- lapply(samp,function(x) x$r) %>% bind_rows(.id = 'sample') %>% 
  group_by(r) %>% 
  summarise(predMean = quantile(predMean, qs), q = qs) %>% 
  ungroup() %>% mutate(q=factor(q,labels=c('lwr','med','upr'))) %>% 
  pivot_wider(names_from=q,values_from=predMean) %>% 
  ggplot(aes(x=r))+
  geom_line(data=allEff[[3]],aes(x=r,y=mean,group=field),alpha=0.1) + #"Raw" smoothers
  geom_ribbon(aes(ymax=upr,ymin=lwr),alpha=0.3,fill='cyan') + geom_line(aes(y=med),col='cyan')+
  geom_hline(yintercept = 0,col='red',linetype='dashed')+
  labs(x='Point order',y=ylabMean)

p6 <- lapply(samp,function(x) x$r) %>% bind_rows(.id = 'sample') %>% 
  group_by(r) %>% 
  summarise(predLogSD = quantile(predLogSD, qs), q = qs) %>% 
  ungroup() %>% mutate(q=factor(q,labels=c('lwr','med','upr'))) %>% 
  pivot_wider(names_from=q,values_from=predLogSD) %>% 
  ggplot(aes(x=r))+
  geom_line(data=allEff[[3]],aes(x=r,y=logSD,group=field),alpha=0.1) + #"Raw" smoothers
  geom_ribbon(aes(ymax=upr,ymin=lwr),alpha=0.3,fill='cyan') + geom_line(aes(y=med),col='cyan')+
  geom_hline(yintercept = 0,col='red',linetype='dashed')+
  labs(x='Point order',y=ylabMean)

(p <- ggarrange(p1,p3,p5,p2,p4,p6,ncol=3,nrow=2))
ggsave(paste0('./Figures/ModelSummary3.png'),p,height=6,width=12,dpi=300)
(p <- ggarrange(p3,p4,ncol=1,nrow=2))
ggsave(paste0('./Figures/ModelSummary3a.png'),p,height=6,width=12,dpi=300)  
