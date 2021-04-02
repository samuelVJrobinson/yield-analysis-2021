# Get boundary shapefiles from fields

library(tidyverse)
library(sf)

rootPath <- "/media/rsamuel/Storage/geoData/Rasters/yieldData/csv files"
datSource <- data.frame(path=dir(rootPath,pattern=".csv",recursive=TRUE)) %>% 
  separate(path,c('grower','year','field'),sep="/",remove=FALSE) %>% 
  mutate(field=gsub('\\.csv','',field)) %>% unite(filename,c(grower:field),sep=' ',remove=FALSE) %>% 
  mutate(completed=filename %in% gsub(' boundary.shp','',dir('./Figures/FieldBoundaries',pattern=".shp",recursive=TRUE))) #Have files already been processed?
  # mutate(completed=FALSE)

# datFiles <- dir('./Figures/ModelCheck/',pattern='.Rdata',full.names=TRUE)

#Function to get boundary from ith field
makeBoundary <- function(i,dS,rP,nSubSamp=50000){ 
  # i <- 3 #Debugging
  # dS <- datSource
  # rP <- rootPath
  # nSubSamp <- 50000
  require(tidyverse)
  require(sf)
  
  source('helperFunctions.R')
  
  csvPath <- paste(rP,dS$path[i],sep='/')
  fieldName <- with(dS[i,],paste(grower,year,field))
  if(dS$completed[i]){ #If shapefile already exists
    return('File already exists')
  }
  
  print('Reading in data')
  dat <- read.csv(csvPath,stringsAsFactors=TRUE,fileEncoding='latin1') 
  
  #Takes 40 seconds with subsamp of 50000 using full 800000 samples from Alvin French's Al Jr Field
  if(nrow(dat)>nSubSamp){
    #Limit to nSubSamp sequential samples
    dat <- dat %>% slice(round(seq(1,nrow(dat),length.out=nSubSamp)))
  }
  
  dat <- dat %>% rename_with(.fn = ~gsub('..L.ha.$','_lHa',.x)) %>%
    rename_with(.fn = ~gsub('..tonne.ha.$','_tHa',.x)) %>%
    rename_with(.fn = ~gsub('..m.s.$','_ms',.x)) %>%
    rename_with(.fn = ~gsub('.km.h.$','_kmh',.x)) %>%
    rename_with(.fn = ~gsub('..tonne.h.$','_th',.x)) %>%
    rename_with(.fn = ~gsub('.ha.h.','_hah',.x)) %>%
    rename_with(.fn = ~gsub('.m.','_m',.x)) %>%
    rename_with(.fn = ~gsub('.deg.','Angle',.x)) %>%
    rename_with(.fn = ~gsub('\\.','',.x)) %>%
    rename('ID'='ObjId','DryYield'='YldMassDry_tHa','Lon'='Longitude','Lat'='Latitude','Pass'='PassNum','Speed'='Speed__m') %>%
    st_as_sf(coords=c('Lon','Lat')) %>% #Add spatial feature info
    st_set_crs(4326) %>% st_transform(3401) %>% #Lat-lon -> UTM
    makePolys(width='SwthWdth_m',dist='Distance_m',angle='TrackAngle') %>%    
    mutate(pArea=as.numeric(st_area(.))) %>% #Area of polygon
    st_centroid() %>% #Convert back to point
    mutate(r=1:n()) %>% #row number
    mutate(Pass=factor(seqGroup(ID,FALSE))) %>% 
    group_by(Pass) %>% mutate(rGroup=1:n()) %>% ungroup() %>% 
    mutate(E=st_coordinates(.)[,1],N=st_coordinates(.)[,2]) %>% 
    mutate(E=E-mean(E),N=N-mean(N)) #Center coordinates
  
  print('Calculating field edge')
  fieldEdge <- dat %>% st_union() %>% st_buffer(dist=10) %>% st_cast('MULTILINESTRING') 
  
  st_write(fieldEdge,paste0('./Figures/FieldBoundaries/',fieldName,' boundary.shp'),append=FALSE)
  
  gc()
}
# makeBoundary(6,datSource,rootPath)

library(parallel)
library(beepr)
cluster <- makeCluster(8) #8 procs max
parLapply(cl=cluster,1:nrow(datSource),makeBoundary,dS=datSource,rP=rootPath) 
beep(1)
stopCluster(cluster)
