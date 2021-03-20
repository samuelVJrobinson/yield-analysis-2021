#Run non-isotropic GAMs on multiple field years

# Load required libraries ---------------------------------------------------------

library(tidyverse)
theme_set(theme_bw())
library(mgcv)
library(sf)

source('helperFunctions.R')

rootPath <- "/media/rsamuel/Storage/geoData/Rasters/yieldData/csv files"
datSource <- data.frame(path=dir(rootPath,pattern=".csv",recursive=TRUE)) %>% 
  separate(path,c('grower','year','field'),sep="/",remove=FALSE) %>% 
  mutate(field=gsub('\\.csv','',field))

datSource %>% count(grower,year) 

# for(i in 1:nrow(datSource)){ #Looped version
#   
#   
# }





  
