#Try loading yield data into R, and run a simple GAM


# Load everything ---------------------------------------------------------

library(tidyverse)
theme_set(theme_bw())
library(mgcv)

# datLoc <- "C:\\Users\\Samuel\\Documents\\Ag Leader Technology\\SMS\\Export\\Trent Clark\\2019\\SW02.csv" #Multivac
datLoc <- "/media/rsamuel/Storage/geoData/Rasters/yieldData/csv files/Trent Clark/SW02.csv" #Multivac


seqGroup <- function(x){ #Turns sequential IDs into numbered groups
  xl <- (x-lag(x))!=1
  cumsum(ifelse(is.na(xl),FALSE,xl))+1
}

dat <- read.csv(datLoc,stringsAsFactors=TRUE) %>% 
  rename('DryYield'='Yld.Mass.Dry..tonne.ha.','Lon'='Longitude','Lat'='Latitude') %>%
  rename('ID'='Obj..Id') %>% 
  mutate(across(Lon:Lat,~as.vector(scale(.x)))) %>% 
  mutate(logYield=log(DryYield)) %>% 
  mutate(r=1:n()) %>% #row number
  mutate(Pass=factor(seqGroup(ID))) %>% 
  group_by(Pass) %>% mutate(rGroup=1:n()) %>% ungroup()

# Look at data ------------------------------------------------------------

dat %>% #Mostly done on Nov 4 2019
  ggplot(aes(x=Lon,y=Lat,col=logYield))+geom_point(alpha=0.7)+ #Yield over entire field
  facet_wrap(~Date,ncol=1)+scale_colour_continuous(type='viridis')

#Dates are split between two places in data 
dat %>% ggplot(aes(x=r,y=Date))+geom_point()

#Was this done on 2 combines?
dat %>% ggplot(aes(x=Lon,y=Lat,col=r))+
  geom_point(alpha=0.7) + geom_path()+
  facet_wrap(~Date,ncol=1)

dat %>% ggplot(aes(x=r,y=ID))+geom_point()+ #ID number seems useful, not Pass.Number
  facet_wrap(~Date)

dat %>% filter(Date!='2019-10-31') %>% 
  ggplot(aes(x=Lon,y=Lat,col=Pass))+
  geom_point(alpha=0.7) + geom_path()

dat %>% ggplot(aes(x=r,y=logYield,col=Pass))+geom_point()+
  facet_wrap(~Date,scales='free_x')


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
