#Try loading yield data into R, and run a simple GAM


# Load everything ---------------------------------------------------------

library(tidyverse)
theme_set(theme_bw())
library(mgcv)
library(sf)

# datLoc <- "C:\\Users\\Samuel\\Documents\\Ag Leader Technology\\SMS\\Export\\Trent Clark\\2019\\SW02.csv" #Multivac
datLoc <- "/media/rsamuel/Storage/geoData/Rasters/yieldData/csv files/Trent Clark/SW02.csv" #Galpern machine

source('helperFunctions.R')

dat <- read.csv(datLoc,stringsAsFactors=TRUE) %>% 
  rename_with(.fn = ~gsub('..L.ha.$','_lHa',.x)) %>%
  rename_with(.fn = ~gsub('..tonne.ha.$','_tHa',.x)) %>%
  rename_with(.fn = ~gsub('..m.s.$','_ms',.x)) %>%
  rename_with(.fn = ~gsub('.km.h.$','_kmh',.x)) %>%
  rename_with(.fn = ~gsub('..tonne.h.$','_th',.x)) %>%
  rename_with(.fn = ~gsub('.ha.h.','_hah',.x)) %>%
  rename_with(.fn = ~gsub('.m.','_m',.x)) %>%
  rename_with(.fn = ~gsub('.deg.','Angle',.x)) %>%
  rename_with(.fn = ~gsub('\\.','',.x)) %>%
  rename('ID'='ObjId','DryYield'='YldMassDry_tHa','Lon'='Longitude','Lat'='Latitude','Pass'='PassNum') %>%
  # select(Lon:Distance,TrackAngle,Pass.Num,DryYield,Date) %>% 
  st_as_sf(coords=c('Lon','Lat')) %>% #Add spatial feature info
  st_set_crs(4326) %>% st_transform(3401) %>% #Lat-lon -> UTM
  mutate(logYield=log(DryYield)) %>% 
  mutate(r=1:n()) %>% #row number
  mutate(Pass=factor(seqGroup(ID))) %>% 
  group_by(Pass) %>% mutate(rGroup=1:n()) %>% ungroup() 

# Look at data ------------------------------------------------------------

dat %>% #Mostly done on Nov 4 2019
  ggplot()+geom_sf(aes(col=logYield),alpha=0.7)+ #Yield over entire field
  facet_wrap(~Date,ncol=1)+scale_colour_continuous(type='viridis')

#Dates are split between two places in data 
dat %>% ggplot(aes(x=r,y=Date))+geom_point()

#Was this done on 2 combines?
dat %>% ggplot()+
  geom_sf(aes(col=r),alpha=0.7) + #geom_path()+
  facet_wrap(~Date,ncol=1)

dat %>% ggplot(aes(x=r,y=ID))+geom_point()+ #ID number seems useful, not Pass.Number
  facet_wrap(~Date)

dat %>% filter(Date!='2019-10-31') %>% 
  ggplot()+
  geom_sf(aes(col=Pass),alpha=0.7) #+ geom_path()

dat %>% ggplot(aes(x=r,y=logYield,col=Pass))+geom_point()+
  facet_wrap(~Date,scales='free_x')

# Deal with overlapping point data ----------------------------------------

#Problem: some points are registered as 0 yield because combine drove over areas that were already harvested.
#Possible solutions: 

dat %>% filter(Date=='2019-10-31',Pass==1) %>% 
  ggplot() + geom_sf()

  # ggplot(aes(x=Lon,y=Lat))+
  # geom_point(alpha=0.7) + geom_path()

dat %>% filter(Date=='2019-10-31',Pass==1) %>% 
  slice(1:3) 

dat %>% filter(Date=='2019-10-31') %>%
  makePolys(width='SwthWdth_m',dist='Distance_m',angle='TrackAngle') %>%
  ggplot()+geom_sf(alpha=0.5)


dat %>% filter(Date=='2019-10-31') %>% 
  makePolys(width='SwthWdth_m',dist='Distance_m',angle='TrackAngle') %>% 
  ggplot()+geom_sf(alpha=0.5)
  # st_intersection() %>%
  # filter(r==1|r==2|r==3) %>% 
  # group_by(r) %>% 
  # st_union() %>% ggplot()+geom_sf(alpha=0.5)#+
  # st_combine() %>% ggplot()+geom_sf(alpha=0.5)#+


  ggplot()+geom_sf(aes(fill=DryYield,col=DryYield),alpha=0.5)#+
  # geom_sf(dat=filter(dat,Date=='2019-10-31'),col='red')

dat %>% filter(Date=='2019-10-31') %>% 
  makePolys(width='SwthWdth_m',dist='Distance_m',angle='TrackAngle') %>% 
  mergePoly()
  # st_intersection() %>% 
  ggplot()+geom_sf(aes(fill=DryYield,col=DryYield),alpha=0.5)


# Additive model ----------------------------------------------------------

library(parallel)
detectCores()
cl <- makeCluster(8)
m1 <- gam(list(sqrt(DryYield)~s(Lon,Lat,k=40)+Pass+s(r,k=40),
               ~ s(Lon,Lat,k=40)+Pass+s(r,k=40)),
          data=dat, #cluster=cl,
          family=gaulss())
summary(m1)
plot(m1,scheme=2,too.far=0.01,pages=1)
par(mfrow=c(2,2)); gam.check(m1); abline(0,1,col='red'); par(mfrow=c(1,1))

#Look at basis dimensions
select(dat,Lon,Lat) %>% bind_cols(as.data.frame(model.matrix(m1))) %>% 
  select(-contains('Intercept')) %>% 
  # select(Lon,Lat,matches('\\.[1-9]$')) %>% #First 9 dimensions
  select(Lon,Lat,matches('\\.([1-3]$|2[4-9]$)')) %>% #First 3, last 6
  pivot_longer(-Lon:-Lat) %>% 
  ggplot(aes(x=Lon,y=Lat,col=value))+geom_point()+
  facet_wrap(~name,ncol=3)
  
stopCluster(cl)
