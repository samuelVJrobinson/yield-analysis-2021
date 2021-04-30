# Get boundary shapefiles from fields

library(tidyverse)
library(sf)

setwd("~/Documents/yield-analysis-2021")

rootPath <- "/media/rsamuel/Storage/geoData/Rasters/yieldData/csv files"
datSource <- data.frame(path=dir(rootPath,pattern=".csv",recursive=TRUE)) %>% 
  separate(path,c('grower','year','field'),sep="/",remove=FALSE) %>% 
  mutate(field=gsub('\\.csv','',field)) %>% unite(filename,c(grower:field),sep=' ',remove=FALSE) %>% 
  mutate(completed=filename %in% gsub(' boundary.shp','',dir('./Figures/FieldBoundaries',pattern=".shp",recursive=TRUE))) #Have files already been processed?
  # mutate(completed=FALSE)

# datFiles <- dir('./Figures/ModelCheck/',pattern='.Rdata',full.names=TRUE)

#Function to get boundary from ith field
makeBoundary <- function(i,dS,rP,outerOnly=TRUE,overwrite=FALSE,nSubSamp=50000){ 
  # #Test
  # # i <- 1 #Debugging
  # i <- 2 #Debugging #2 = 2-part field
  # dS <- datSource
  # rP <- rootPath
  # nSubSamp <- 50000
  # outerOnly <- TRUE
  
  require(tidyverse)
  require(sf)
  
  source('helperFunctions.R')
  
  csvPath <- paste(rP,dS$path[i],sep='/')
  fieldName <- with(dS[i,],paste(grower,year,field))
  if(dS$completed[i] & !overwrite){ #If shapefile already exists and overwrite=FALSE
    return('File already exists')
  }
  
  cat('Reading in data.')
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
  
  cat(' Calculating field edge.')
  #Get all boundaries (including internal holes)
  fieldEdge <- dat %>% st_union() %>% 
    st_buffer(dist=10) %>% st_buffer(dist=-5) #Dilate by 10 then erode by 5
  
  if(outerOnly){ #Gets only outer boundary
    if(st_geometry_type(fieldEdge)=='POLYGON'){
      fieldEdge2 <- lapply(fieldEdge, function(x) st_polygon(x[1]))
    } else if(st_geometry_type(fieldEdge)=='MULTIPOLYGON'){
      fieldEdge2 <- lapply(fieldEdge, function(x) lapply(x,function(y) st_polygon(y[1])))
      fieldEdge2 <- fieldEdge2[[1]] %>% do.call('st_union',.) #Union of polygons
    } else {
      stop('Polygon structure error')
    }
    fieldEdge2 <- st_sfc(fieldEdge2,crs=st_crs(fieldEdge)) #Turn into sfc
  } else { #Gets all boundaries (including internal)
    fieldEdge2 <- fieldEdge
  }
  
  #Write to shapefile
  cat('. Writing to shapefile.')
  st_write(fieldEdge2,paste0('./Figures/FieldBoundaries/',fieldName,' boundary.shp'),append=FALSE,driver='ESRI Shapefile')
  rm(dat,fieldEdge,fieldEdge2); gc() #Garbage collection
}

# makeBoundary(58,datSource,rootPath,overwrite=TRUE) #Test
# debugonce(makeBoundary)


library(parallel)
library(beepr)
cluster <- makeCluster(8) #8 procs max
parLapply(cl=cluster,1:nrow(datSource),makeBoundary,dS=datSource,rP=rootPath,outerOnly=TRUE) 
beep(1)
stopCluster(cluster)

# Turn polygons into linestrings ------------------------------------------
library(tidyverse)
library(sf)

setwd("~/Documents/yield-analysis-2021")

rootPath <- "/media/rsamuel/Storage/geoData/Rasters/yieldData/csv files"
# datSource <- data.frame(path=dir(rootPath,pattern=".csv",recursive=TRUE)) %>% 
#   separate(path,c('grower','year','field'),sep="/",remove=FALSE) %>% 
#   mutate(field=gsub('\\.csv','',field)) %>% unite(filename,c(grower:field),sep=' ',remove=FALSE) %>% 
#   mutate(completed=filename %in% gsub(' boundary.shp','',dir('./Figures/FieldBoundaries',pattern=".shp",recursive=TRUE)))

datSource <- read.csv('./Data/datSource.csv') #Read previous datsource file
datSource <- datSource %>% mutate(boundaryComplete=file.exists(boundaryPath)) %>% #Has boundary been made already?
  mutate(modelComplete=file.exists(modelPath)) #Has model already been run?
# shpFiles <- paste0('./Figures/FieldBoundaries/',dir('./Figures/FieldBoundaries',pattern=".shp",recursive=TRUE)) #Paths to shapefiles
shpFiles <- unique(datSource$boundaryPath[datSource$use]) #Only unique shapefiles used (some fields use a single year's boundary)

makeLinestring <- function(shp){
  # shp <-shpFiles[1] 
  fieldName <- strsplit(shp,'/')[[1]][length(strsplit(shp,'/')[[1]])]
  #Polygon
  p <- read_sf(shp) %>% st_cast('MULTILINESTRING') %>% mutate(type='STANDARD') %>% st_simplify(dTolerance=0.5)
  st_write(p,paste0('./Figures/FieldBoundaryLines/',fieldName),append=FALSE,driver='ESRI Shapefile')
}
makeLinestring(shpFiles[1])
sapply(shpFiles,makeLinestring)
