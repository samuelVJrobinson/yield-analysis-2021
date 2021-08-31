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
  mutate(modelComplete0=file.exists(modelPath0)) #Has model0 already been run?

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

# #Test with single model
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

# Smoother info from second set of models, with sub-model error included ------------------

#Draw posterior samples from each field and fit curve at each field

# Function to sample from posterior of each prediction, amalgamate, fit meta-model, return results of meta-model
# My way of getting around the problem of random effects
# a = does nothing (mandatory argument for parLapply), useRows = rows of ds to be used, ds = datSource csv, nX = number of samples along range of X
# rCutoff = restrict point number (some fields have hundreds of thousands of points, but all data uses approximately the same size yield rectangles, so cuts off meta-smoother at "about" the maximum field size - i.e. the most common number of points per field)
samplePreds <- function(a=NULL,useRows=NULL,ds=datSource,nX=c(30,100,500),rCutoff=50000,margInt=c(FALSE,TRUE,TRUE)){ 
  if(is.null(useRows)){
    warning('Rows not specified. Using entire dataset.')
    useRows <- 1:nrow(datSource)
  } 
  require(tidyverse); require(mgcv); require(lme4)
  source('helperFunctions.R')
  getFiles <- ds %>% slice(useRows) %>% filter(use,modelComplete2) #Completed models
  
  #Sample from posterior, get new lines for each field
  allSmooths <- lapply(getFiles$modelPath2,getPreds,l=nX,margInt=margInt,samp=TRUE) %>% set_names(getFiles$filename)
  
  #log(Area) meta-models
  pAreaDf <- lapply(allSmooths,function(x) x$pAreaDat) %>% bind_rows(.id = 'field')
  
  #Adds a small amount of noise so model fit isn't singular
  m1 <- pAreaDf %>% 
    mutate(mean=mean+rnorm(n(),0,(max(mean)-min(mean))*0.005)) %>% 
    lmer(mean~log(pArea)+(log(pArea)|field),data=.) 
  m2 <- pAreaDf %>% mutate(logSD=logSD+rnorm(n(),0,(max(logSD)-min(logSD))*0.005)) %>% 
    lmer(logSD~log(pArea)+(log(pArea)|field),data=.)
  
  #cover class meta-models
  
  #Get cover class names
  coverNames <- sort(unique(unlist(lapply(allSmooths, function(x) names(x$distDat)))))
  
  #Which models have which cover classes
  containsCover <- sapply(coverNames,function(z) sapply(lapply(allSmooths, function(x) names(x$distDat)), function(y) z %in% y))
  
  #Lists for cover distance models
  meanModList <- logSDModList <- vector(mode='list',length=length(coverNames)) %>% set_names(coverNames)
  
  for(i in 1:ncol(containsCover)){
    useThese <- unname(which(containsCover[,i]))
    cname <- colnames(containsCover)[i]
    coverDF <- lapply(allSmooths[useThese],function(x) x$distDat[[cname]]) %>% 
      bind_rows(.id = 'field') %>% mutate(field=factor(field))
    meanModList[[i]] <- gam(mean~s(dist,k=20)+s(field,bs='re'),data=coverDF)
    logSDModList[[i]] <- gam(logSD~s(dist,k=20)+s(field,bs='re'),data=coverDF)
  }
  
  #r (point order) meta-models
  rDatDf <- lapply(allSmooths,function(x) x$rDat) %>% 
    bind_rows(.id = 'field') %>% mutate(field=factor(field)) %>% 
    filter(r<=rCutoff) #Filters out point values above cutoff
  
  r1 <- gam(mean~s(r,k=50)+s(field,bs='re'),data=rDatDf)
  r2 <- gam(logSD~s(r,k=50)+s(field,bs='re'),data=rDatDf)
  
  #Assemble predictions
  pAreaPred <- data.frame(pArea=exp(seq(log(min(pAreaDf$pArea)),log(max(pAreaDf$pArea)),length.out=nX[1]))) %>% 
    mutate(predMean=predict(m1,newdata=.,re.form=~0),predLogSD=predict(m2,newdata=.,re.form=~0))
  
  # coverPred <- lapply(1:length(coverNames),function(i){
  #   useThese <- unname(which(containsCover[,i]))
  #   cname <- colnames(containsCover)[i]
  #   coverDF <- lapply(allSmooths[useThese],function(x) x$distDat[[cname]]) %>% 
  #     bind_rows(.id = 'field') %>% mutate(field=factor(field))
  #   
  #   data.frame(dist=seq(min(coverDF$dist),max(coverDF$dist),length.out=nX[2]),field=coverDF$field[1]) %>% 
  #     mutate(predMean=predict(meanModList[[i]],newdata=.,exclude="s(field)"),predLogSD=predict(logSDModList[[i]],newdata=.,exclude="s(field)"))
  # }) %>% set_names(gsub('dist:boundaryType','',coverNames))

  #Empty list for results  
  coverPred <- vector(mode='list',length=length(coverNames)) %>% set_names(gsub('dist:boundaryType','',coverNames))

  for(i in 1:length(coverNames)){

    useThese <- unname(which(containsCover[,i]))
    cname <- colnames(containsCover)[i]
    coverDF <- lapply(allSmooths[useThese],function(x) x$distDat[[cname]]) %>%
      bind_rows(.id = 'field') %>%  mutate(field=factor(field))

    coverPred[[i]] <- data.frame(dist=seq(min(coverDF$dist),max(coverDF$dist),length.out=nX[2]),field=coverDF$field[1]) %>%
      mutate(predMean=predict(meanModList[[i]],newdata=.,exclude="s(field)"),predLogSD=predict(logSDModList[[i]],newdata=.,exclude="s(field)"))
  }
  
  rPred <- data.frame(r=seq(min(rDatDf$r),max(rDatDf$r),length.out=nX[3]),field=rDatDf$field[1]) %>% 
    mutate(predMean=predict(r1,newdata=.,exclude="s(field)"),predLogSD=predict(r2,newdata=.,exclude="s(field)"))
  
  retList <- list(pArea=pAreaPred,coverDist=coverPred,r=rPred) #Return list of results
  
  rm(list=ls()[!ls() %in% 'retList']) #Remove everything except retList
  gc()
  return(retList) 
}

debugonce(samplePreds)
samplePreds() #Test

# Get crop-specific smoother info from second set of models 

## Add to current samples - canola
# isCanola <- which(datSource$crop=='Canola') #Canola crops only
# Nsamp <- 1000 #Number of samples
# library(parallel)
# cluster <- makeCluster(15) #Memory usage is OK, so could probably max it out
# clusterExport(cluster,c('datSource'))
# a <- Sys.time() #53 seconds for 15 reps, using 15 cores - 45 mins for 1000 samples
# samp <- parLapply(cl=cluster,1:Nsamp,fun=samplePreds,useRows=isCanola,margInt=c(FALSE,TRUE,TRUE))
# Sys.time()-a
# stopCluster(cluster)
# # load('./Data/postSamples_canola.Rdata')
# # samp <- c(samp,samp2)
# save(samp,file='./Data/postSamples_canola.Rdata')

# # Add to current samples - wheat
# isWheat <- which(datSource$crop=='Wheat') #Wheat
# Nsamp <- 1000 #Number of samples
# library(parallel)
# cluster <- makeCluster(15)
# clusterExport(cluster,c('datSource'))
# a <- Sys.time() #38 mins for 1000 reps
# samp <- parLapply(cl=cluster,1:Nsamp,fun=samplePreds,useRows=isWheat,margInt=c(TRUE,TRUE,TRUE))
# Sys.time()-a
# stopCluster(cluster)
# # load('./Data/postSamples_wheat.Rdata')
# # samp <- c(samp,samp2)
# save(samp,file='./Data/postSamples_wheat.Rdata')
# rm(samp2)

# croptype <- 'canola'
croptype <- 'wheat'

load(paste0('./Data/postSamples_',croptype,'.Rdata'))

# tHa2buAc <- 17.0340 #tonnes per hectare to bushels per acre (https://www.agrimoney.com/calculators/calculators)
ylabMean <- 'Mean Yield (T/ha)'
ylabSD <- 'log SD Yield'
ylims <- c(-1.5,1.5)
distLims <- c(0,200)
alphaVal <- 0.1
qs <- c(0.1,0.5,0.9) #Quantiles
colshade <- 'black'
fillshade <- 'black'

p1 <- lapply(samp,function(x) x$pArea) %>% bind_rows(.id = 'sample') %>% 
  group_by(pArea) %>% 
  summarise(predMean = quantile(predMean, qs), q = qs) %>% 
  ungroup() %>% mutate(q=factor(q,labels=c('lwr','med','upr'))) %>% 
  pivot_wider(names_from=q,values_from=predMean) %>% 
  ggplot(aes(x=pArea))+
  # geom_line(data=allEff[[1]],aes(x=pArea,y=mean,group=field),alpha=0.3)+
  geom_ribbon(aes(ymax=upr,ymin=lwr),alpha=0.3,fill=fillshade)+
  geom_line(aes(y=med),col=colshade,size=1)+
  labs(x='Polygon Area (m^2)',y=ylabMean)

p2 <- lapply(samp,function(x) x$pArea) %>% bind_rows(.id = 'sample') %>% 
  group_by(pArea) %>% 
  summarise(predLogSD = quantile(predLogSD, qs), q = qs) %>% 
  ungroup() %>% mutate(q=factor(q,labels=c('lwr','med','upr'))) %>% 
  pivot_wider(names_from=q,values_from=predLogSD) %>% 
  ggplot(aes(x=pArea))+
  # geom_line(data=allEff[[1]],aes(x=pArea,y=logSD,group=field),alpha=0.3)+
  geom_ribbon(aes(ymax=upr,ymin=lwr),alpha=0.3,fill=fillshade)+
  geom_line(aes(y=med),col=colshade,size=1)+
  labs(x='Polygon Area (m^2)',y=ylabSD)

p3 <- lapply(samp,function(x) x$coverDist %>% bind_rows(.id = 'dist_type')) %>% 
  bind_rows(.id = 'sample') %>% group_by(dist_type,dist) %>% 
  summarise(predMean = quantile(predMean, qs), q = qs) %>% 
  ungroup() %>% mutate(q=factor(q,labels=c('lwr','med','upr'))) %>% 
  pivot_wider(names_from=q,values_from=predMean) %>% 
  filter(dist>=distLims[1],dist<=distLims[2]) %>% 
  ggplot(aes(x=dist)) + 
  # geom_line(data=allEff[[2]],aes(x=dist,y=mean,group=field),alpha=0.1) + #"Raw" smoothers
  geom_ribbon(aes(ymax=upr,ymin=lwr),alpha=0.3,fill=fillshade) + 
  geom_line(aes(y=med),col=colshade) + #Meta-model
  geom_hline(yintercept = 0,col='red',linetype='dashed')+
  facet_wrap(~dist_type,scales='free') + labs(x='Distance (m)',y=ylabMean)+ 
  coord_cartesian(xlim=distLims,ylim=ylims)

p4 <- lapply(samp,function(x) x$coverDist %>% bind_rows(.id = 'dist_type')) %>% 
  bind_rows(.id = 'sample') %>% group_by(dist_type,dist) %>% 
  summarise(predLogSD = quantile(predLogSD, qs), q = qs) %>% 
  ungroup() %>% mutate(q=factor(q,labels=c('lwr','med','upr'))) %>% 
  pivot_wider(names_from=q,values_from=predLogSD) %>% 
  filter(dist>=distLims[1],dist<=distLims[2]) %>% 
  ggplot(aes(x=dist)) + 
  # geom_line(data=allEff[[2]],aes(x=dist,y=logSD,group=field),alpha=0.1) + #"Raw" smoothers
  geom_ribbon(aes(ymax=upr,ymin=lwr),alpha=0.3,fill=fillshade) + geom_line(aes(y=med),col=colshade) + #Meta-model
  geom_hline(yintercept = 0,col='red',linetype='dashed')+
  facet_wrap(~dist_type) + labs(x='Distance (m)',y=ylabSD) +
  coord_cartesian(xlim=distLims,ylim=ylims)

p5 <- lapply(samp,function(x) x$r) %>% bind_rows(.id = 'sample') %>% 
  group_by(r) %>% 
  summarise(predMean = quantile(predMean, qs), q = qs) %>% 
  ungroup() %>% mutate(q=factor(q,labels=c('lwr','med','upr'))) %>% 
  pivot_wider(names_from=q,values_from=predMean) %>% 
  ggplot(aes(x=r))+
  # geom_line(data=allEff[[3]],aes(x=r,y=mean,group=field),alpha=0.1) + #"Raw" smoothers
  geom_ribbon(aes(ymax=upr,ymin=lwr),alpha=0.3,fill=fillshade) + geom_line(aes(y=med),col=colshade)+
  geom_hline(yintercept = 0,col='red',linetype='dashed')+
  labs(x='Harvest order',y=ylabMean)

p6 <- lapply(samp,function(x) x$r) %>% bind_rows(.id = 'sample') %>% 
  group_by(r) %>% 
  summarise(predLogSD = quantile(predLogSD, qs), q = qs) %>% 
  ungroup() %>% mutate(q=factor(q,labels=c('lwr','med','upr'))) %>% 
  pivot_wider(names_from=q,values_from=predLogSD) %>% 
  ggplot(aes(x=r))+
  # geom_line(data=allEff[[3]],aes(x=r,y=logSD,group=field),alpha=0.1) + #"Raw" smoothers
  geom_ribbon(aes(ymax=upr,ymin=lwr),alpha=0.3,fill=fillshade) + geom_line(aes(y=med),col=colshade)+
  geom_hline(yintercept = 0,col='red',linetype='dashed')+
  labs(x='Harvest order',y=ylabSD)

(p <- ggarrange(p1,p3,p5,p2,p4,p6,ncol=3,nrow=2)) #Plot of everything
ggsave(paste0('./Figures/ModelSummary3_',croptype,'.png'),p,height=6,width=16,dpi=350)

(p <- ggarrange(p3+facet_wrap(~dist_type,nrow=1),p4+facet_wrap(~dist_type,nrow=1),ncol=1,nrow=2)) #Edge distance smoothers only
ggsave(paste0('./Figures/ModelSummary3a_',croptype,'.png'),p,height=6,width=12,dpi=350)  

(p <- ggarrange(p1,p5,p2,p6,ncol=2,nrow=2)) #Polygon area and harvest order only
ggsave(paste0('./Figures/ModelSummary3b_',croptype,'.png'),p,height=6,width=8,dpi=350)  

# Smoother info from 2nd set of models, using combine ground speed rather than polygon area --------------------------------------------------------------

#Get coefs for converting pArea (width x distance) to ground speed - takes about 1 minute
convertArea <- lapply(datSource$dataPath, function(path){
  d <- read.csv(path,nrows=10000,fileEncoding='latin1')
  #Get statistics on swath width and speed
  minSwathWidth <- min(d$Swth.Wdth.m.,na.rm=TRUE)
  maxSwathWidth <- max(d$Swth.Wdth.m.,na.rm=TRUE)
  # minSpeed <- min(d$Speed.km.h.[d$Speed.km.h.>0])
  minSpeed <- unname(quantile(d$Speed.km.h.[d$Speed.km.h.>0],0.1)) #10th percentile instead of min
  maxSpeed <- unname(quantile(d$Speed.km.h.,0.9)) #90th percentile of speed instead of max
  d <- subset(d,d$Swth.Wdth.m.==maxSwathWidth) #Retain measurements at max swath width
  if(nrow(d)==0) return(list(maxSwathWidth=NA,dist2Speed=NA,r2=NA))
  m <- lm(Speed.km.h.~Distance.m.-1,data=d) #Predict speed (km/hr) using distance (m)
  dist2Speed <- unname(coef(m)) #Slope coefficient to convert distance to speed
  retList <- list(minSwathWidth=minSwathWidth,maxSwathWidth=maxSwathWidth,
                  minSpeed=minSpeed,maxSpeed=maxSpeed,
                  dist2Speed=dist2Speed, #Slope of speed:distance line
                  r2=summary(m)$r.squared)
  rm(list=ls()[!ls() %in% 'retList']) #Remove everything except retList
  gc()
  return(retList)
})

#Shorter version that only gets estimates of log(pArea) relationship and turns them into speed measurements
# Combine ground speed is talked about a lot, so this is more appropriate for growers/agronomists

getPredsSpeed <- function(path,l=30,margInt=FALSE,samp=FALSE, minSpeed = NA, maxSpeed = NA, width = NA, speed2Distance = NA){
  
  # #Debugging
  # l <- 30
  # # path <- paste0('./Figures/ModelCheck1/',datSource$filename[2],' modList.Rdata') #Model type 1
  # path <- paste0('./Figures/ModelCheck2/',datSource$filename[2],' modList.Rdata') #Model type 2
  # minSpeed <- convertArea[[2]]$minSpeed
  # maxSpeed <- convertArea[[2]]$maxSpeed
  # margInt <- FALSE
  # samp <- FALSE
  # width <- convertArea[[2]]$maxSwathWidth
  # speed2Distance <- 1/convertArea[[2]]$dist2Speed
  # rm(l,path,minSpeed,maxSpeed,margInt,samp,width,speed2Distance,logSpeed,speed,WAS,logWAS)
  
  if(any(is.na(c(minSpeed,maxSpeed,width,speed2Distance)))) stop('Missing term')
  load(path) #Load data
  
  #Converts speed to area term for use with original coefficients (holding width at max value)
  logSpeed <- seq(log(minSpeed),log(maxSpeed),length.out=l) 
  speed <- exp(logSpeed) #Speed measurements to use
  WAS <- width*speed*speed2Distance #Width * speed * alpha = creates proxy for pArea
  logWAS <- log(WAS) #proxy for log(pArea)
  
  meanVars <- which(grepl('(^\\(Intercept\\)$|^log\\(pArea\\)$)',names(modList$coefs)))
  sdVars <- which(grepl('(^\\(Intercept\\)\\.1$|^log\\(pArea\\)\\.1$)',names(modList$coefs)))
  
  #Coefficients (intercept and slope) for mean and logSD relationship
  if(samp){ #Sample from posterior
    meanCoefs <- rnorm(rep(1,length(meanVars)),modList$coefs[meanVars],sqrt(diag(modList$vcv)[meanVars]))
    sdCoefs <- rnorm(rep(1,length(meanVars)),modList$coefs[sdVars],sqrt(diag(modList$vcv)[sdVars]))
  } else {
    meanCoefs <- modList$coefs[meanVars] 
    sdCoefs <- modList$coefs[sdVars]
  }
  
  if(margInt){ #Marginalize across intercept (set intercept coef to 0)
    meanCoefs[1] <- 0 
    sdCoefs[1] <- 0  
  }
  
  speedDat <- data.frame(speed, 
                         mean=cbind(rep(1,length(logSpeed)),logWAS) %*% meanCoefs,
                         logSD=cbind(rep(1,length(logSpeed)),logWAS) %*% sdCoefs)
  
  rm(list=ls()[!ls() %in% 'speedDat']) #Remove everything except speedDat
  gc()
  return(speedDat)
} 

# Same function as samplePreds, but only gets ground speed measurements (after converting to pArea)
# a = does nothing (argument for parLapply), ds = datSource csv, nX = number of samples along range of X
samplePredsSpeed <- function(a=NULL,ds=datSource,nX=100,minSpeeds=NA,maxSpeeds=NA,maxWidth=NA,s2d=NA,margInt=FALSE,useLmer=TRUE){ 
  
  # #Debugging
  # ds <- datSource; nX <- 100
  # useThese <- with(ds,use&modelComplete2&crop=='Wheat')
  # minSpeeds <- sapply(convertArea,function(x) x$minSpeed)[useThese]
  # maxSpeeds <- sapply(convertArea,function(x) x$maxSpeed)[useThese]
  # maxWidth <- sapply(convertArea,function(x) x$maxSwathWidth)[useThese]
  # s2d <- 1/sapply(convertArea,function(x) x$dist2Speed)[useThese]
  # ds <- ds %>% filter(useThese) #Completed models
  # # rm(ds,nX,useThese,minSpeeds,maxSpeeds,maxWidth,s2d,ds,allSmooths,speedDf)

  require(tidyverse)
  source('helperFunctions.R')
  
  #Sample from posterior, get new lines for each field
  allSmooths <- mapply(FUN = getPredsSpeed,path=ds$modelPath2,minSpeed=minSpeeds,maxSpeed=maxSpeeds,
                       width=maxWidth,speed2Distance=s2d,MoreArgs=list(l=30,margInt=margInt,samp=TRUE),
                       SIMPLIFY=FALSE) %>% set_names(nm=ds$filename)
  
  #log(Speed) meta-models
  speedDf <- allSmooths %>% bind_rows(.id = 'field')
  
  # speedDf %>% filter(speed<15) %>% #Plot of field-level estimates
  #   pivot_longer(cols=mean:logSD) %>% 
  #   ggplot(aes(x=speed,y=value,group=field))+geom_line(size=0.1)+
  #   facet_wrap(~name,scales='free',ncol=1)
  if(useLmer){
    library(lme4)
    m1 <- speedDf %>% mutate(mean=mean+rnorm(length(logSD),0,(max(mean)-min(mean))*0.01)) %>% 
      lmer(mean~log(speed)+(log(speed)|field),data=.) 
    m2 <- speedDf %>% mutate(logSD=logSD+rnorm(length(logSD),0,(max(logSD)-min(logSD))*0.01)) %>% 
      lmer(logSD~log(speed)+(log(speed)|field),data=.)
    
    #Assemble predictions
    speedPred <- data.frame(speed=exp(seq(log(min(speedDf$speed)),log(max(speedDf$speed)),length.out=nX))) %>% 
      mutate(predMean=predict(m1,newdata=.,re.form=~0),predLogSD=predict(m2,newdata=.,re.form=~0))
    
  } else {
    m1 <- lm(mean~log(speed),data=speedDf) 
    m2 <- lm(logSD~log(speed),data=speedDf)
    
    #Assemble predictions
    speedPred <- data.frame(speed=exp(seq(log(min(speedDf$speed)),log(max(speedDf$speed)),length.out=nX))) %>% 
      mutate(predMean=predict(m1,newdata=.),predLogSD=predict(m2,newdata=.))
  }
  rm(list=ls()[!ls() %in% 'speedPred']) #Remove everything except speedPred
  gc()
  return(speedPred) 
}

# library(parallel)
# cluster <- makeCluster(15)
# clusterExport(cluster,c('getPredsSpeed'))
# samp <- parLapply(cl=cluster,1:Nsamp,fun=samplePredsSpeed,ds=tempDat,minSpeeds=minSpeeds,maxSpeeds=maxSpeeds,maxWidth=maxWidth,s2d=s2d)
# #Takes about 10 mins for 100 samples
# stopCluster(cluster)
# save(samp,file='./Data/postSamplesSpeed.Rdata')

# Add to current samples - canola
Nsamp <- 30 #Number of samples
useThese <- with(datSource,use&modelComplete2&crop=='Canola')
tempDat <- datSource %>% filter(useThese)
minSpeeds <- sapply(convertArea,function(x) x$minSpeed)[useThese]
maxSpeeds <- sapply(convertArea,function(x) x$maxSpeed)[useThese]
maxWidth <- sapply(convertArea,function(x) x$maxSwathWidth)[useThese]
s2d <- 1/sapply(convertArea,function(x) x$dist2Speed)[useThese]
library(parallel)
cluster <- makeCluster(15)
clusterExport(cluster,c('getPredsSpeed'))
samp <- parLapply(cl=cluster,1:Nsamp,fun=samplePredsSpeed,ds=tempDat,minSpeeds=minSpeeds,maxSpeeds=maxSpeeds,maxWidth=maxWidth,s2d=s2d,margInt=FALSE)
stopCluster(cluster)
# load('./Data/postSamplesSpeed_canola.Rdata')
# samp <- c(samp,samp2)
save(samp,file='./Data/postSamplesSpeed_canola.Rdata')
 
# #TEST - negative relationship b/w speed and mean/SD, but only shows up if margInt = FALSE. Think a bit about how to deal with this/
# #Solution: use lmer to deal with this
# useThese <- with(datSource,use&modelComplete2&crop=='Canola')
# tempDat <- datSource %>% filter(useThese)
# minSpeeds <- sapply(convertArea,function(x) x$minSpeed)[useThese]
# maxSpeeds <- sapply(convertArea,function(x) x$maxSpeed)[useThese]
# maxWidth <- sapply(convertArea,function(x) x$maxSwathWidth)[useThese]
# s2d <- 1/sapply(convertArea,function(x) x$dist2Speed)[useThese]
# debugonce(samplePredsSpeed)
# test <- samplePredsSpeed(1,ds=tempDat,minSpeeds=minSpeeds,maxSpeeds=maxSpeeds,maxWidth=maxWidth,s2d=s2d,margInt=FALSE)

# Make samples - wheat
Nsamp <- 30 #Number of samples
useThese <- with(datSource,use&modelComplete2&crop=='Wheat')
tempDat <- datSource %>% filter(useThese)
minSpeeds <- sapply(convertArea,function(x) x$minSpeed)[useThese]
maxSpeeds <- sapply(convertArea,function(x) x$maxSpeed)[useThese]
maxWidth <- sapply(convertArea,function(x) x$maxSwathWidth)[useThese]
s2d <- 1/sapply(convertArea,function(x) x$dist2Speed)[useThese]
library(parallel)
cluster <- makeCluster(15)
clusterExport(cluster,c('getPredsSpeed'))
samp <- parLapply(cl=cluster,1:Nsamp,fun=samplePredsSpeed,ds=tempDat,minSpeeds=minSpeeds,maxSpeeds=maxSpeeds,maxWidth=maxWidth,s2d=s2d,margInt=FALSE)
stopCluster(cluster)
# load('./Data/postSamplesSpeed_wheat.Rdata')
# samp <- c(samp,samp2)
save(samp,file='./Data/postSamplesSpeed_wheat.Rdata')
rm(samp,samp2)


#Make ground speed:yield plots
croptype <- 'canola'
croptype <- 'wheat'

load(paste0('./Data/postSamplesSpeed_',croptype,'.Rdata'))

qs <- c(0.1,0.5,0.9) #Quantiles
p1 <- samp %>% bind_rows(.id = 'sample') %>% group_by(speed) %>% 
  summarise(predMean = quantile(predMean, qs), q = qs) %>% 
  ungroup() %>% mutate(q=factor(q,labels=c('lwr','med','upr'))) %>% 
  pivot_wider(names_from=q,values_from=predMean) %>% 
  ggplot(aes(x=speed)) + geom_ribbon(aes(ymax=upr,ymin=lwr),alpha=0.3) +
  geom_line(aes(y=med),size=1) + labs(x='Ground Speed (km/hr)',y='Mean yield (T/ha)')

p2 <- samp %>% bind_rows(.id = 'sample') %>% group_by(speed) %>% 
  summarise(predLogSD = quantile(predLogSD, qs), q = qs) %>% 
  ungroup() %>% mutate(q=factor(q,labels=c('lwr','med','upr'))) %>% 
  pivot_wider(names_from=q,values_from=predLogSD) %>% 
  ggplot(aes(x=speed))+
  geom_ribbon(aes(ymax=upr,ymin=lwr),alpha=0.3)+
  geom_line(aes(y=med),size=1)+
  labs(x='Ground Speed (km/hr)',y='log SD yield')

(p <- ggarrange(p1,p2,ncol=1))
ggsave(paste0('./Figures/groundSpeed_',croptype,'.png'),p,height=8,width=4,dpi=350)  


# 
# #Function to map f(pArea) to f(speed|width,alpha); b0 = intercept, b1=slope of y~log(pArea), alpha=convert from distance to speed
# function(b0,b1,width,alpha,speed) b0 + b1*log(width*alpha*speed)

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