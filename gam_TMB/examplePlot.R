#Example plot of basis functions

library(mgcv)
library(tidyverse)

N <- 100
dat <- data.frame(x=runif(N,-10,10)) %>% 
  mutate(y=1-0.5*x+0.2*x^2+rnorm(N)) %>% 
  arrange(x)

plot(y~x,data=dat)

m1 <- gam(y~s(x,k=6,bs = 'cr'),data=dat)

par(mfrow=c(2,1),mar=c(4,4.1, 1, 1))
plot(dat$x,model.matrix(m1)[,2],col='red',type='l',ylim=c(-0.6,1),xlab='Distance',ylab='Z')
lines(dat$x,model.matrix(m1)[,3],col='blue')
lines(dat$x,model.matrix(m1)[,4],col='green')
lines(dat$x,model.matrix(m1)[,5],col='magenta')
lines(dat$x,model.matrix(m1)[,6],col='grey40')

plot(dat$x,model.matrix(m1)[,c(2:6)]%*%coef(m1)[2:6],type='l',lwd=2,ylim=c(-10,20),
     xlab='Distance',ylab='Zu')
lines(dat$x,model.matrix(m1)[,2]*coef(m1)[2],col='red')
lines(dat$x,model.matrix(m1)[,3]*coef(m1)[3],col='blue')
lines(dat$x,model.matrix(m1)[,4]*coef(m1)[4],col='green')
lines(dat$x,model.matrix(m1)[,5]*coef(m1)[5],col='magenta')
lines(dat$x,model.matrix(m1)[,6]*coef(m1)[6],col='grey40')


m1$smooth[[1]]

Z <- m1$smooth[[1]]$UZ
