#Try loading yield data into R, and run a simple GAM

library(tidyverse)
library(mgcv)

dat <- read.csv("C:\\Users\\Samuel\\Documents\\Ag Leader Technology\\SMS\\Export\\Alvin French\\D42.csv") %>% 
  rename('DryYield'='Yld.Mass.Dry..tonne.ha.','Lon'='Longitude','Lat'='Latitude') %>%
  mutate(across(Lon:Lat,~as.vector(scale(.x)))) %>% 
  mutate(logYield=log(DryYield))

ggplot(dat,aes(x=Lon,y=Lat,col=logYield))+geom_point()


library(parallel)
detectCores()
cl <- makeCluster(6)
m1 <- bam(DryYield~s(Lon,Lat),data=dat,cluster=cl)
plot(m1,scheme=2)
