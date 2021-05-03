#TAKING A LOOK AT YIELD DATA FROM TRENT CLARK - PRELIMINARY EXPLORATION

library(tidyverse)
library(sf)
library(mgcv)
theme_update(panel.border=element_rect(size=1,fill=NA),axis.line=element_blank()) #Better for maps

#L/bu
35.239
#Ac/ha
2.471

#L/ha -> bu/ac
L2bu <- 2.471/35.239

#Canola yield data from JOHNSON field only
yieldDat <- read.csv("/media/samuel/Storage/Eyeshigh/Rasters/yieldData/TrentClark/Converted/2019/all2019.csv",header=T,fileEncoding="latin1") %>% 
  st_as_sf(coords=c('Longitude','Latitude'),crs=4326) %>% 
  filter(Product=='CANOLA') %>% filter(Field=='JOHNSON') %>%
  # select(-Area.Count:-Crop.Flw.M..tonne.h.,-contains('Wet')) %>% 
  st_transform(crs=3403) %>% 
  mutate(xloc=st_coordinates(.)[,1],yloc=st_coordinates(.)[,2]) %>% 
  rename('Speed'='Speed.km.h.','YieldDryTHa'='Yld.Mass.Dry..tonne.ha.','YieldDryLHa'='Yld.Vol.Dry..L.ha.') %>% 
  rename('ID'='Obj..Id') %>% 
  mutate(Yield=YieldDryTHa*1000) %>% 
  mutate(Yield=ifelse(Yield>quantile(Yield,probs=c(0.975)),NA,Yield)) %>% 
  mutate(sYield=(Yield-mean(Yield,na.rm=T))/sd(Yield,na.rm=T)) %>% #Scaled yield
  mutate(chunk=factor(cumsum(as.numeric((ID-lag(ID))<0&!is.na(lag(ID))))+1))  #"Chunks" caused by different passes (new machine?)
  
#Average yield
yieldDat %>% st_drop_geometry() %>% 
  mutate(tonnesYield=YieldDryTHa*Swth.Wdth.m.*Distance.m./10000) %>% 
  group_by(Product,Field) %>% 
  summarize(tonnes=sum(tonnesYield)) %>% 
  data.frame()

#Yield
yieldDat %>% ggplot() + geom_histogram(aes(Yield))+geom_vline(xintercept=2326.2,col='red')+labs(x='Yield (kg/ha)')

yieldDat %>% mutate(Yield=ifelse(Yield>quantile(Yield,probs=c(0.975)),NA,Yield)) %>%
  ggplot() + geom_histogram(aes(Yield))+geom_vline(xintercept=2326.2,col='red')+labs(x='Yield (kg/ha)')

yieldDat %>% mutate(Yield=ifelse(Yield>quantile(Yield,probs=c(0.975)),NA,Yield)) %>%
  ggplot() + geom_sf(aes(col=Yield),size=1)+
  scale_colour_gradient(low='blue',high='red')+
  facet_wrap(~Dataset)

#Speed
yieldDat %>% ggplot() + geom_histogram(aes(Speed))+labs(x='Speed km/hr')

yieldDat %>% ggplot() + geom_sf(aes(col=Speed),size=1) +
  scale_colour_gradient(low='blue',high='red')

#Speed/Yield  
yieldDat %>% ggplot(aes(Speed,Yield))+geom_point()
  # scale_x_log10()+scale_y_log10()
  # scale_x_sqrt()+scale_y_sqrt()

#ObjectID - order of points?
yieldDat %>% 
  ggplot()+geom_sf(aes(col=ID))+
  facet_grid(chunk~Dataset)

yieldDat %>% mutate(rownum=1:nrow(.)) %>% 
  ggplot()+geom_point(aes(rownum,ID))+
  facet_wrap(~Dataset,ncol=1)

#Test models 
mod1 <- gam(sYield~Speed+te(xloc,yloc,k=12)+s(ID,by=chunk,k=10),data=yieldDat,subset=!is.na(Yield))
summary(mod1)
plot(mod1,se=F,scheme=2,pages=1)

yieldDat %>% select(Speed,xloc,yloc,sYield,ID,chunk) %>% filter(!is.na(sYield)) %>% 
  mutate(Speed=mean(Speed)) %>% 
  mutate(pred=predict(mod1,newdata=.),resid=resid(mod1)) %>% st_drop_geometry() %>% 
  pivot_longer(cols = c('sYield','pred','resid'),names_to='type') %>% 
  st_as_sf(coords=c('xloc','yloc'),crs=3403) %>% 
  ggplot() + geom_sf(aes(col=value))+
  facet_wrap(~type,ncol=1)

yieldDat %>% select(Speed,xloc,yloc,sYield) %>% filter(!is.na(sYield)) %>% 
  mutate(Speed=mean(Speed)) %>% 
  mutate(pred=predict(mod1,newdata=.)) %>% 
  ggplot()+geom_sf(aes(col=pred))


