#Run non-stationary GAMs on multiple field years

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
#   mutate(modelPath1=paste0('./Figures/ModelCheck1/',filename,' modList.Rdata')) %>% #Paths to saved model files
#   mutate(modelPath2=paste0('./Figures/ModelCheck2/',filename,' modList.Rdata')) %>% 
#   mutate(modelPath0=paste0('./Figures/ModelCheck0/',filename,' modList.Rdata')) 
# 
# #Read 10000 lines from each csv, and get counts of products - takes about 20 seconds
# cropType <- sapply(datSource$dataPath,function(path){
#   table(read.csv(path,nrows=10000,fileEncoding='latin1')$Product)}) %>% 
#   set_names(datSource$filename) %>% 
#   sapply(.,function(x) names(x)[1])
# 
# 
# #Get crop types
# datSource$crop <- case_when(grepl('(CANOLA|InVigor|Dekalb 74-54|Dekalb 75-42|Brett Young 6056|Dupont D3152|VT 9562)',cropType) ~ 'Canola',
#                             grepl('(WHEAT|Wheat|CPS AC Foremost|CWRS|CWSW)',cropType) ~ 'Wheat',
#                             grepl('Peas',cropType) ~ 'Peas',
#                             grepl('(BARLEY|Barley)',cropType) ~ 'Barley',
#                             grepl('Beans',cropType) ~ 'Beans',
#                             grepl('Oats',cropType) ~ 'Oats',
#                             grepl('Rye',cropType) ~ 'Rye',
#                             TRUE ~ cropType) 
# datSource$npoints <- sapply(datSource$dataPath,function(x) length(readLines(x)))
# write.csv(datSource,'./Data/datSource_new.csv',row.names = FALSE)

# #Get length of boundary features at each field
# bLengths <- lapply(datSource$boundaryPath2,function(x){
#   bTypes <- c('BARE','GRASSLAND','OTHERCROP','SHELTERBELT','STANDARD','WETLAND')
#   if(!file.exists(x)){
#     return(data.frame(type=bTypes,len=NA))
#   }
#   st_read(x,quiet=TRUE) %>%
#     mutate(len=as.numeric(st_length(.))) %>%
#     filter(len>0) %>%
#     mutate(type=factor(type,levels=bTypes)) %>%
#     st_drop_geometry() %>% group_by(type,.drop=FALSE) %>%
#     summarize(len=sum(len))}) %>%
#   set_names(nm = datSource$filename) %>%
#   bind_rows(.id='filename') %>% group_by(filename) %>% mutate(prop=len/sum(len)) %>%
#   pivot_wider(names_from=type,values_from=c(len,prop)) %>%
#   data.frame()
# write.csv(bLengths,file = './Data/boundaryLengths.csv')

datSource <- read.csv('./Data/datSource.csv') %>% #Read previous datasource file
  mutate(boundaryComplete=file.exists(boundaryPath)) %>% #Has boundary been made already?
  mutate(modelComplete1=file.exists(modelPath1)) %>% #Has model1 already been run?
  mutate(modelComplete2=file.exists(modelPath2)) %>% #Has model2 already been run?
  mutate(modelComplete0=file.exists(modelPath0)) %>% #Has model0 already been run?
  mutate(location=ifelse(grepl('(Alvin|Dean)',filename),'Southern','Central')) #Field locations (southern or central AB)

#Get summary results from models
datSource %>% pull(use) %>% sum()
datSource %>% filter(use) %>% count(grower,year) #How many fields per grower per year
datSource %>% group_by(grower) %>% summarize(nLoc=length(unique(field)),nYears=length(unique(year)),n=n()) #Number of distinct fields 
datSource %>% filter(use) %>% count(crop) #Enough to do separate analyses on canola and wheat, maybe peas

# #First set of models - no boundary types
# resultPaths <- gsub(' modList.Rdata',' results.txt',datSource$modelPath1)[datSource$use]
# modInfo <- lapply(resultPaths,getModelInfo)
# sapply(modInfo,function(x) as.numeric(x$timeTaken,units='hours')) %>% summary() #Hours taken
# sapply(modInfo,function(x) x$percDevianceExplained) %>% summary() #Explained deviance

# #Second set of models - boundary types included
# resultPaths <- gsub(' modList.Rdata',' results.txt',datSource$modelPath2)[datSource$use]
# modInfo <- lapply(resultPaths,getModelInfo)
# sum(!file.exists(datSource$modelPath2[datSource$use])) #24 models not completed
# sapply(modInfo,function(x) as.numeric(x$timeTaken,units='hours')) %>% summary() #Hours taken
# sapply(modInfo,function(x) x$percDevianceExplained) %>% summary() #Explained deviance


# Run "zero-th" set of models - no boundary -------------------------------

# a <- Sys.time() #Test
# runMod0(1,dS=datSource)
# Sys.time()-a
# beep(1)

library(parallel)
cluster <- makeCluster(15) #10 procs max - uses about 90% of memory
a <- Sys.time()
parLapply(cl=cluster,1:nrow(datSource),runMod0,dS=datSource) 
stopCluster(cluster)
Sys.time() - a

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

# # Test with single model
# which(datSource$filename=='Gibbons Lafonda 2018') #18 mins
# which(datSource$filename=='Alvin_French Al_Jr 2020') #12 mins
# a <- Sys.time()
# runModII(113,dS=datSource,kPar=c(5,60,5,60))
# Sys.time()-a
# # beep(1)
# debugonce(runModII)

library(parallel)
detectCores()
cluster <- makeCluster(15) #10 procs uses about 30% of memory - could probably max it out
a <- Sys.time()
runWhich <- c(1:nrow(datSource))[with(datSource,use&!modelComplete2&grepl('(Canola|Wheat|Peas)',crop))] #Models to use, 
parLapply(cl=cluster,runWhich,runModII,dS=datSource,useClosest=TRUE) #All models
beep(1)
stopCluster(cluster)
Sys.time()-a #Takes ~ 10 hrs for 50 models at 15 procs
#Takes 2.364893 days to do all fields at 15 procs

# Get smoother info from second set of models ---------------------------

#Estimate at single field
# debugonce(getPreds)
est <- getPreds(paste0('./Figures/ModelCheck2/',datSource$filename[13],' modList.Rdata'),
                margInt=c(FALSE,TRUE,TRUE),
                samp=FALSE)$distDat %>% bind_rows(.id='dist_type')

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

#Get dataframe of predictions for each variable
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

# Smoother info from second set of models, with sub-model error included ------------------

#Draw posterior samples from each field and fit curve at each field

# Get crop-specific smoother info from second set of models 

# # Add to current samples - canola
# isCanola <- which(with(datSource,use & crop=='Canola' & modelComplete2)) #Canola crops only
# Nsamp <- 500 #Number of samples
# library(parallel)
# cluster <- makeCluster(15) #Memory usage is OK, so could probably max it out
# clusterExport(cluster,c('datSource'))
# a <- Sys.time()
# samp <- parLapply(cl=cluster,1:Nsamp,fun=samplePreds,useRows=isCanola,margInt=c(FALSE,FALSE,FALSE),kPar=10)
# Sys.time()-a
# stopCluster(cluster)
# save(samp,file='./Data/postSamples_canola.Rdata')
# 
# # Add to current samples - wheat
# isWheat <- which(datSource$crop=='Wheat') #Wheat
# Nsamp <- 500 #Number of samples
# library(parallel)
# cluster <- makeCluster(15)
# clusterExport(cluster,c('datSource'))
# a <- Sys.time()
# samp <- parLapply(cl=cluster,1:Nsamp,fun=samplePreds,useRows=isWheat,margInt=c(FALSE,FALSE,FALSE),kPar=10)
# Sys.time()-a
# stopCluster(cluster)
# save(samp,file='./Data/postSamples_wheat.Rdata')
# 
# # Add to current samples - peas
# isPeas <- which(datSource$crop=='Peas') #Peas
# Nsamp <- 500 #Number of samples
# library(parallel)
# cluster <- makeCluster(15)
# clusterExport(cluster,c('datSource'))
# a <- Sys.time()
# samp <- parLapply(cl=cluster,1:Nsamp,fun=samplePreds,useRows=isPeas,margInt=c(FALSE,FALSE,FALSE),kPar=10)
# Sys.time()-a
# stopCluster(cluster)
# save(samp,file='./Data/postSamples_peas.Rdata')

#Get samples from storage
croptype <- c('canola','wheat','peas')
samp2 <- vector(mode='list',length=3) %>% set_names(croptype)
for(i in 1:length(croptype)){
  load(paste0('./Data/postSamples_',croptype[i],'.Rdata'))
  samp <- lapply(samp,backTrans) #Back-transform
  samp2[[i]] <- samp
} #Store in nested lists
names(samp) <- croptype
samp <- samp2; rm(samp2)

#Get baseline smoother fit to sampled data ("mean" effect)
samp_mean <- vector(mode='list',length=3) %>% set_names(croptype)
for(i in 1:length(croptype)){
  isCrop <- which(tolower(datSource$crop)==croptype[i] & datSource$use)
  samp_mean[[i]] <- backTrans(samplePreds(1,useRows=isCrop,samp=FALSE))
} #Takes about 10 secs

# tHa2buAc <- 17.0340 #tonnes per hectare to bushels per acre (https://www.agrimoney.com/calculators/calculators)
ylabMean <- 'Average Yield (t/ha)'
ylabSD <- 'Yield variablity (log(t/ha)) '

distLims <- c(0,200)
alphaVal <- 0.1
qs <- c(0.05,0.5,0.95) #Quantiles
colshade <- 'black'
fillshade <- 'black'

#Number of fields that have at least _some_ boundary type
nContains <- lapply(c('Canola','Wheat','Peas'),function(cr){
  isCrop <- which(with(datSource,cr==crop & use & modelComplete2))
  allSmooths <- lapply(datSource$modelPath2[isCrop],getPreds,samp=FALSE,reps=1)
  
  #Get cover class names
  coverNames <- sort(unique(unlist(lapply(allSmooths, function(x) names(x$distDat)))))
  
  #Which models have which cover classes
  containsCover <- sapply(coverNames,function(z) sapply(lapply(allSmooths, function(x) names(x$distDat)), function(y) z %in% y))
  ret <- data.frame(type = set_names(gsub('dist:boundaryType','',coverNames)),nfields=apply(containsCover,2,sum))
  return(ret)
}) %>% set_names(nm = c('canola','wheat','peas'))

#Within those fields, what makes up proportion of boundary?
propBoundary <- lapply(c('Canola','Wheat','Peas'),function(cr){
  isCrop <- which(with(datSource,cr==crop & use & modelComplete2))
  
  ret <- lapply(datSource$boundaryPath2[isCrop],function(x){
    st_read(x,quiet=TRUE) %>% mutate(len=as.numeric(st_length(.))) %>% 
    filter(len>0) %>% st_drop_geometry() %>% #Get rid of zero-length segments
    group_by(type) %>% summarize(sumLength=sum(len)) %>% data.frame()
  }) %>% bind_rows(.id='field') %>% mutate(type=factor(type)) %>% group_by(field) %>% mutate(total=sum(sumLength)) %>% 
    ungroup %>% mutate(propLen=sumLength/total) %>% 
    group_by(type) %>% summarize(avgProp=mean(propLen)) %>% 
    data.frame()
  
  return(ret)
}) %>% set_names(nm = c('canola','wheat','peas'))

#Assemble figures
figList <- vector(mode = 'list',length = 3) %>% set_names(croptype)

for(i in 1:length(croptype)){
  
  if(croptype[i]=='canola'){
    ylimsMean <- c(1,15)
    ylimsSD <- c(-7,0)
  } else if(croptype[i]=='wheat'){
    ylimsMean <- c(2,8)
    ylimsSD <- c(-6,1)
  } else if(croptype[i]=='peas'){
    ylimsMean <- c(2,10)
    ylimsSD <- c(-7,2)
  }
  
  meanSmooth <- samp_mean[[i]]$coverDist %>% bind_rows(.id = 'type') %>% 
    left_join(nContains[[i]],by='type') %>%
    left_join(propBoundary[[i]],by='type') %>% 
    mutate(type=paste0(type,' (N=',nfields,', %=',round(avgProp*100,1),')'))
  
  (p1 <-lapply(samp[[i]],function(x) x$coverDist %>% bind_rows(.id = 'type')) %>% 
    bind_rows(.id = 'sample') %>% group_by(type,dist) %>% 
    summarise(predMean = quantile(predMean, qs), q = qs) %>% 
    ungroup() %>% mutate(q=factor(q,labels=c('lwr','med','upr'))) %>% 
    pivot_wider(names_from=q,values_from=predMean) %>% 
    filter(dist>=distLims[1],dist<=distLims[2]) %>% 
    left_join(nContains[[i]],by='type') %>%
    left_join(propBoundary[[i]],by='type') %>% 
    mutate(type=paste0(type,' (N=',nfields,', %=',round(avgProp*100,1),')')) %>% 
    ggplot(aes(x=dist)) + 
    geom_ribbon(aes(ymax=upr,ymin=lwr),alpha=0.3,fill=fillshade) + 
    geom_line(aes(y=med),col=colshade) + #Meta-model median
    # geom_line(data=meanSmooth,aes(y=predMean))+
    facet_wrap(~type,nrow=1) + labs(x='Distance (m)',y=ylabMean) 
    # +
    # coord_cartesian(xlim=distLims,ylim=ylimsMean)
    )
  
  (p2 <- lapply(samp[[i]],function(x) x$coverDist %>% bind_rows(.id = 'type')) %>% 
      bind_rows(.id = 'sample') %>% group_by(type,dist) %>% 
      summarise(predLogSD = quantile(predLogSD, qs), q = qs) %>% 
      ungroup() %>% mutate(q=factor(q,labels=c('lwr','med','upr'))) %>% 
      pivot_wider(names_from=q,values_from=predLogSD) %>% 
      filter(dist>=distLims[1],dist<=distLims[2]) %>% 
      left_join(nContains[[i]],by='type') %>%
      left_join(propBoundary[[i]],by='type') %>% 
      mutate(type=paste0(type,' (N=',nfields,', %=',round(avgProp*100,1),')')) %>% 
      ggplot(aes(x=dist)) + 
      geom_ribbon(aes(ymax=upr,ymin=lwr),alpha=0.3,fill=fillshade) + 
      geom_line(aes(y=med),col=colshade) + #Meta-model
      # geom_line(data=meanSmooth,aes(y=predLogSD))+
      facet_wrap(~type,nrow=1) + labs(x='Distance (m)',y=ylabSD)
    # +
    #   coord_cartesian(xlim=distLims,ylim=ylimsSD)
    )
  
  figList[[i]] <- list(p1,p2)
}


# (p <- ggarrange(p1,p3,p5,p2,p4,p6,ncol=3,nrow=2)) #Plot of everything
# ggsave(paste0('./Figures/ModelSummary3_',croptype,'.png'),p,height=6,width=16,dpi=350)

#Write to files for each crop types
for(i in 1:length(croptype)){
  p1 <- figList[[i]][[1]]
  p2 <- figList[[i]][[2]]
    
  (p <- ggarrange(p1,p2,ncol=1,nrow=2)) #Edge distance smoothers only
  ggsave(paste0('./Figures/ModelSummary3a_',croptype[i],'.png'),p,height=6,width=12,dpi=350)    
}


# Example figures from single field (Trent Clark Johnson 2014) ---------------------------------------------------------

# use <- which(datSource$filename == 'Trent_Clark W 34 2014')
use <- which(datSource$filename == 'Trent_Clark JOHNSON 2014') #Could also use this

load(datSource$modelPath2[use]) #Load model
# load("C:\\Users\\Samuel\\Documents\\Projects\\UofC postdoc\\yield-analysis-2021\\Figures\\ExamplePlots\\Trent_Clark JOHNSON 2014 modList.Rdata") #Path on Multivac

#Spatial smoothers

library(sf)
fieldBoundary <- st_read(datSource$boundaryPath[use]) #Get boundary
fieldBoundaryType <- st_read(datSource$boundaryPath2[use]) %>% mutate(type=factor(type)) #Get boundary
crs <- st_crs(fieldBoundary) #Save CRS
fieldBoundary <- st_sfc(st_polygon(lapply(fieldBoundary$geometry,function(x) st_coordinates(x)[,c('X','Y')])),crs=crs) #Fix hole geometry
hexGrid <- st_make_grid(fieldBoundary,square=FALSE,n=75)  #Make hexagonal grid
hexGrid <- hexGrid[sapply(st_within(st_centroid(hexGrid),fieldBoundary),function(x) length(x)==1)] %>%  #Strip out points outside polygon
  st_sf() #Set as sf object
# hexGrid %>% ggplot()+geom_sf()+geom_sf(data=fieldBoundary,fill=NA,col='red') #Looks OK

dat <- read.csv(datSource$dataPath[use],stringsAsFactors=TRUE,fileEncoding='latin1') %>% 
  st_as_sf(coords=c('Longitude','Latitude')) %>% #Add spatial feature info
  st_set_crs(4326) %>% st_transform(3401) %>% #Lat-lon -> UTM
  transmute(Speed=Speed.km.h.,DryYield=Yld.Mass.Dry..tonne.ha.,Dist=Distance.m.,Swath=Swth.Wdth.m.,
            pArea=Swath*Speed*1/3.6,PassNum=Pass.Num,ID=Obj..Id,
            bearingDiff=Track.deg.-lag(Track.deg.),
            E=st_coordinates(.)[,1],N=st_coordinates(.)[,2])
meanE <- mean(dat$E); meanN <- mean(dat$N)
dat <- dat %>% mutate(E=E-meanE,N=N-meanN) %>% #Center coordinates
  filter(pArea>quantile(pArea,0.05),pArea<quantile(pArea,0.95), #Filter large/small pArea
         DryYield>quantile(DryYield,0.05),DryYield<quantile(DryYield,0.95), #Filter extreme yields
         Speed>quantile(Speed,0.05),Speed<quantile(Speed,0.95) #Filter high and low speeds
  )

hexGrid <- hexGrid %>%  #Center coordinates
  transmute(E=st_coordinates(st_centroid(.))[,1],N=st_coordinates(st_centroid(.))[,2]) %>% 
  mutate(E=E-meanE,N=N-meanN)

theme_set(theme_bw())

# palette <- 'RdYlBu'
# palette <- 'Spectral'
palette <- 'YlGnBu'

#Raw data + 
(p1 <-ggplot()+
    geom_sf(data=dat,aes(geometry=geometry,col=DryYield),size=0.5,show.legend = 'point')+
    geom_sf(data=fieldBoundaryType,aes(geometry=geometry),show.legend= 'line')+
    scale_colour_distiller(type='div',palette = palette, direction = -1)+
    labs(col='Yield (t/ha)',title='Raw Yield Data')+
    theme(legend.position='bottom',axis.text = element_blank(),axis.ticks = element_blank())
  )

# temp <- lapply(levels(fieldBoundaryType$type),function(x) filter(fieldBoundaryType,type==x)) %>% 
#   set_names(nm=levels(fieldBoundaryType$type))
# 
# p+
#   geom_sf(data=temp$GRASSLAND,aes(geometry=geometry),col='green')+
#   geom_sf(data=temp$OTHERCROP,aes(geometry=geometry),col='black')+
#   geom_sf(data=temp$WETLAND,aes(geometry=geometry),col='darkgreen')+
#   geom_sf(data=temp$SHELTERBELT,aes(geometry=geometry),col='brown')+
# ggsave(paste0('./Figures/ExamplePlots/rawData.png'),p,height=8,width=8,dpi=350)  
  
sapply(modList$smooths,function(x) x$label)

useSmooths <- which(grepl('(E,N)',sapply(modList$smooths,function(x) x$label),fixed=TRUE)) #Spatial smoothers only

#Spatial smoothers
(p2 <- hexGrid %>% 
    bind_cols(getSmooths(smoothLabel=modList$smooths[[useSmooths[1]]]$label,modList=modList,xvals=st_drop_geometry(hexGrid[,c('E','N')]),
                         noIntercept=FALSE,returnSE = TRUE)[,c('pred','se')]) %>% 
    mutate(pred=pred^2) %>% #Back-transform
    ggplot()+geom_sf(col=NA,aes(fill=pred))+
    geom_sf(data=fieldBoundaryType,aes(geometry=geometry),show.legend= 'line')+
    labs(fill='Yield (t/ha)',title='Yield Average Smoother')+
    scale_fill_distiller(type='div',palette = palette, direction = -1) +
    theme(legend.position='bottom',axis.text = element_blank(),axis.ticks = element_blank())
  )

(p3 <- hexGrid %>% 
  bind_cols(getSmooths(smoothLabel=modList$smooths[[useSmooths[2]]]$label,modList=modList,xvals=st_drop_geometry(hexGrid[,c('E','N')]),noIntercept=FALSE,returnSE = TRUE)[,c('pred','se')]) %>% 
  mutate(pred=log(exp(pred)^2)) %>% #Back-transform
  ggplot()+geom_sf(col=NA,aes(fill=pred))+
    geom_sf(data=fieldBoundaryType,aes(geometry=geometry),show.legend= 'line')+
    labs(fill='log Yield SD',title='Yield Variability Smoother')+
    scale_fill_distiller(type='div',palette = palette, direction = -1) +
    theme(legend.position='bottom',axis.text = element_blank(),axis.ticks = element_blank())
  )

(p <- ggarrange(p1,p2,p3,ncol=3,nrow=1))
ggsave(paste0('./Figures/ExamplePlots/spatialSmooths.png'),p,height=6,width=12,dpi=350)  

#Distance smoothers
useSmooths <- grepl('dist',sapply(modList$smooths,function(x) x$label)) #Distance smoothers only

smoothDat <- lapply(modList$smooths[useSmooths],function(x){
  xvals <- x$Xu[,1]+as.vector(x$shift)
  xvals <- seq(min(xvals),max(xvals),length.out=150)
  getSmooths(x$label,modList,xvals,noIntercept = FALSE, returnSE = TRUE) #1D smooth  
}) %>% bind_rows() %>% 
  mutate(type=rep(ifelse(grepl('s.1',sapply(modList$smooths,function(x) x$label)[useSmooths],fixed = TRUE),'logSD','mean'),each=150)) %>% 
  mutate(upr=pred+se*1.96,lwr=pred-se*1.96) %>% 
  filter(dist<200)

(p1 <- smoothDat %>% filter(type=='mean') %>% mutate(across(c(pred,upr,lwr),~.x^2)) %>% 
  ggplot(aes(x=dist))+geom_ribbon(aes(ymax=upr,ymin=lwr),alpha=0.3)+geom_line(aes(y=pred))+
  facet_wrap(~boundaryType,nrow=1)+
  labs(x='Distance (m)',y='Yield (t/ha)'))

(p2 <- smoothDat %>% filter(type=='logSD') %>% mutate(across(c(pred,upr,lwr),~log(exp(.x)^2))) %>% 
  ggplot(aes(x=dist))+geom_ribbon(aes(ymax=upr,ymin=lwr),alpha=0.3)+geom_line(aes(y=pred))+
  facet_wrap(~boundaryType,nrow=1)+
  labs(x='Distance (m)',y='log Yield SD'))

(p <- ggarrange(p1,p2,ncol=1))
ggsave(paste0('./Figures/ExamplePlots/distSmooths.png'),p,height=6,width=12,dpi=350)  

#Plots of basis functions

# sapply(modList$smooths,function(x) x$term)

colnum <- 15 #How many bases to show?
dat <- expand.grid(E=seq(-400,400,length.out=21),N=seq(-400,400,length.out=21))
Z <- PredictMat(modList$smooths[[5]],data=dat)[,1:colnum] #Columns 1-colnum
colnames(Z) <- paste0('Z',formatC(1:colnum,width=2,flag='0'))

data.frame(dat,Z) %>% 
  pivot_longer(cols=contains('Z')) %>% 
  ggplot(aes(x=E,y=N,fill=value))+geom_raster(show.legend = FALSE)+
  facet_wrap(~name,nrow=3)+
  scale_fill_distiller(type='div',palette = palette, direction = -1) 

# Compare model types -----------------------------------------------------

# Get model info
mod0info <- lapply(gsub('modList.Rdata','results.txt',datSource$modelPath0),getModelInfo)
mod1info <- lapply(gsub('modList.Rdata','results.txt',datSource$modelPath1),getModelInfo)
mod2info <- lapply(gsub('modList.Rdata','results.txt',datSource$modelPath2),getModelInfo)

datSource %>% mutate(reml0=sapply(mod0info,function(x) x$REML),
                     reml1=sapply(mod1info,function(x) x$REML),
                     reml2=sapply(mod2info,function(x) x$REML),
                     diff1=reml0-reml1,diff2=reml0-reml2)

load(datSource$modelPath0[1])
# Simulate different types of boundaries at the edge of a quarter-section ---------------

#Get samples from storage
croptype <- c('canola','wheat','peas')
samp2 <- vector(mode='list',length=3) %>% set_names(croptype)
for(i in 1:length(croptype)){
  load(paste0('./Data/postSamples_',croptype[i],'.Rdata'))
  # samp <- lapply(samp,backTrans) #Back-transform
  samp2[[i]] <- samp
} #Store in nested lists
names(samp) <- croptype
samp <- samp2; rm(samp2)

#Get baseline smoother fit to sampled data ("mean" effect)
samp_mean <- vector(mode='list',length=3) %>% set_names(croptype)
for(i in 1:length(croptype)){
  isCrop <- which(tolower(datSource$crop)==croptype[i] & datSource$use)
  # samp_mean[[i]] <- backTrans(samplePreds(1,useRows=isCrop,samp=FALSE))
  samp_mean[[i]] <- samplePreds(1,useRows=isCrop,samp=FALSE)
} #Takes about 10 secs



lapply(samp$wheat,function(x) x$coverDist[c('STANDARD','SHELTERBELT','WETLAND')]) %>% set_names(paste0('samp',1:length(.))) %>% 
  lapply(.,function(x) bind_rows(x,.id='boundary') %>% filter(dist<400)) %>% 
  bind_rows(.id='samp') %>% 
  pivot_longer(c(predMean,predLogSD)) %>% 
  ggplot(aes(x=dist))+
  geom_line(aes(y=value,group=samp),alpha=0.1)+
  facet_grid(name~boundary,scales='free')


#Create field boundary polygon
qSecDist <- 805/4 #Quarter section = 805 x 805 m
fieldPoly <- st_polygon(x=list(matrix(c(0,0,0,qSecDist,qSecDist,qSecDist,qSecDist,0,0,0),nrow=5,ncol=2,byrow = TRUE))) #Square field
fieldBoundary <- lapply(1:(length(st_coordinates(fieldPoly)[, 1]) - 1), #4 boundary segments
                        function(i){
                          rbind(
                            as.numeric(st_coordinates(fieldPoly)[i, 1:2]),
                            as.numeric(st_coordinates(fieldPoly)[i + 1, 1:2])
                          )
                        }) %>%
    lapply(.,function(i) st_linestring(i)) %>% 
    st_sfc() %>% 
    st_sf(a=c('a','b','c','d'),geom=.)
sampGrid <- st_make_grid(fieldPoly,cellsize = qSecDist/100,square=TRUE) #100 x 100 cells

#Problem: cells are equidistant from some edges, and no systematic way to make sure edges are equally represented.
#Solution: add very small rotation to the grid - "nearest edge" changes but not distance
rotGrd <- function(g,a){ #Function to rotate grid g by a degrees
  rot <- function(aa) matrix(c(cos(aa), sin(aa), -sin(aa), cos(aa)), 2, 2)  
  (g - st_centroid(st_union(g))) * rot(a * pi / 180) + st_centroid(st_union(g))
}
rotSampGrid <- rotGrd(sampGrid,0.1) #Rotate by 0.1 degrees
cSamp <- st_centroid(sampGrid) #Centroid of each cell 
cSampRot <- st_centroid(rotSampGrid) #Centroid of rotated cells

#Sampling grid, with nearest edges coded as a-d
sampGrid <- st_sf(nearest=letters[1:4][apply(st_distance(cSampRot,fieldBoundary),1,function(x) which.min(x)[1])],
      dist=apply(st_distance(cSamp,fieldBoundary),1,min),
      geom=st_geometry(sampGrid))
sampGrid %>% ggplot()+geom_sf(aes(fill=nearest)) #Looks OK
rm(qSecDist,fieldPoly,fieldBoundary)

# #Single sample
# samp1 <- samp$canola[[1]]$coverDist
# grid1 <- sampGrid %>% mutate(m=NA,sd=NA,y=NA) %>% st_drop_geometry()
# 
# #Looks OK
# samp1 %>% bind_rows(.id='boundary') %>% filter(dist<400) %>% 
#   pivot_longer(predMean:predLogSD) %>% 
#   ggplot(aes(x=dist,y=value))+geom_point()+
#   facet_grid(name~boundary,scales='free_y')

#Function to simulate yield given:
# boundaries = dataframe of boundaries: nearest = letters[1:4], boundary = corresponding name of boundary type (eg "STANDARD")
# crop = 'canola','wheat','peas'
# samp = posterior samples from storage
# grid = sampling 
simYield <- function(boundaries,crop='canola',samp,grid,Nreps=1){
  cropSamp <- lapply(1:length(samp[[crop]]), function(i) samp[[crop]][[i]]$coverDist) #Get samples from specific crop
  f1 <- function(s1,b){ #Replace boundary of grid, fill with means
    g1 <- grid
    for(i in 1:nrow(b)){
      isb <- grid$nearest==b$nearest[i]
      g1$nearest[isb] <- b$boundary[i]
      g1$m[isb] <- approx(x = s1[[i]]$dist, y = s1[[i]]$predMean, xout = g1$dist[isb])$y
      g1$sd[isb] <- exp(approx(x = s1[[i]]$dist, y = s1[[i]]$predLogSD, xout = g1$dist[isb])$y)
    }
    g1$y <- (rnorm(nrow(g1),g1$m,g1$sd))^2 #Simulate yield, then transform to normal scale
    return(mean(g1$y)) #All cells are the same size, so mean is OK
  }
  simY <- replicate(Nreps, #Replicate Nreps replicates
                    sapply(cropSamp, f1, b=boundaries), #Get simulation for each cropSamp
                    simplify = FALSE)
  return(unlist(simY))
}

#Test
test <- data.frame(b=rep(c('STANDARD','SHELTERBELT','WETLAND'),each=1000),
           y=c(simYield(data.frame(nearest = letters[1:4], boundary= rep('STANDARD',4)),'canola',samp,sampGrid,Nreps=2),
               simYield(data.frame(nearest = letters[1:4], boundary= rep('SHELTERBELT',4)),'canola',samp,sampGrid,Nreps=2),
               simYield(data.frame(nearest = letters[1:4], boundary= rep('WETLAND',4)),'canola',samp,sampGrid,Nreps=2))) 
test %>% ggplot(aes(x=y*17.84,fill=b))+geom_density(alpha=0.2)+xlim(NA,60)+labs(x='Bu/acre canola',fill='Boundary')

test <- data.frame(b=rep(c('STANDARD','SHELTERBELT','WETLAND'),each=1000),
                   y=c(simYield(data.frame(nearest = letters[1:4], boundary= rep('STANDARD',4)),'wheat',samp,sampGrid,Nreps=2),
                       simYield(data.frame(nearest = letters[1:4], boundary= rep('SHELTERBELT',4)),'wheat',samp,sampGrid,Nreps=2),
                       simYield(data.frame(nearest = letters[1:4], boundary= rep('WETLAND',4)),'wheat',samp,sampGrid,Nreps=2))) 
test %>% ggplot(aes(x=y*14.85,fill=b))+geom_density(alpha=0.2)+labs(x='Bu/acre wheat',fill='Boundary')
 
# #Sets of boundaries to use
# boundaries <- lapply(0:4,function(i) data.frame(nearest = letters[1:4], boundary= c(rep('STANDARD',4-i),rep('SHELTERBELT',i))))
# 
# #Shelterbelt
# sb_canola <- lapply(boundaries,simYield,crop='canola',samp=samp,grid=sampGrid,Nreps=10) %>% 
#   set_names(paste('shelterbelt',0:4,sep='_')) %>% 
#   bind_rows(.id = 'scenario') %>% 
#   pivot_longer(everything()) %>% mutate(crop='canola')
# sb_wheat <- lapply(boundaries,simYield,crop='wheat',samp=samp,grid=sampGrid,Nreps=10) %>% 
#   set_names(paste('shelterbelt',0:4,sep='_')) %>% 
#   bind_rows(.id = 'scenario') %>% 
#   pivot_longer(everything()) %>% mutate(crop='wheat')
# 
# list(sb_canola,sb_wheat) %>% set_names(c('canola','wheat')) %>% 
#   bind_rows(.id='croptype') %>% 
#   ggplot(aes(x=value))+geom_density(aes(fill=name))+
#   facet_wrap(~croptype,ncol=1,scales='free')
# 
# 
