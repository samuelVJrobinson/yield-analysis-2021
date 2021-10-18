
# First attempt: make functions from scratch ------------------------------

dist <- 1:100

expFun1 <- function(x,a,b){ #Yield loss around field edges
  a*exp(-(1/x)^b)
}

curve(expFun1(x,1,3),0,100,ylim=c(0,1))

par(mfrow=c(1,1))

#Beneficial insect presence - works for -ive relationships, but not positive

expFun2 <- function(x,a,b){ 
  a*exp(-(x/b))
}
curve(expFun2(x,a=0.5,b=100),0,100,ylim=c(0,1),ylab='Bugs',xlab='Dist')
curve(expFun2(x,a=0.3,b=100),0,100,add=TRUE)
curve(expFun2(x,a=0.2,b=-100),0,100,add=TRUE)

par(mfrow=c(2,1))
#Scenario 1 - beneficial bugs decrease with distance into field
curve(expFun1(x,a=1,b=1)*(expFun2(x,a=0.5,b=100)+0.5),0,100,ylim=c(0,1),xlim=range(dist),
      lwd=4,ylab='Yield',main='Scenario 1')
curve(expFun1(x,a=1,b=1),0,100,ylim=c(0,1),xlim=range(dist),col='red',add=TRUE,lwd=2)
curve(expFun2(x,a=0.5,b=100)+0.5,0,100,ylim=c(0,1),xlim=range(dist),col='blue',add=TRUE,lwd=2)
text(5,0.9,'Bugs',col='blue'); text(80,0.9,'Environment',col='red'); text(80,0.6,'Yield')

#Scenario 2 - beneficial bugs increase with distance into field
curve(expFun1(x,a=1,b=3)*(expFun2(x,a=0.5,b=-100)+0.5),0,100,xlim=range(dist),
      lwd=4,ylab='Yield',main='Scenario 2')
curve(expFun1(x,a=1,b=3),0,100,ylim=c(0,1),xlim=range(dist),col='red',add=TRUE,lwd=2)
curve(expFun2(x,a=0.5,b=-100)+0.5,0,100,ylim=c(0,1),xlim=range(dist),col='blue',add=TRUE,lwd=2)
text(5,0.9,'Bugs',col='blue'); text(80,0.9,'Environment',col='red'); text(80,0.6,'Yield')
par(mfrow=c(1,1))


# Second attempt: use logistic function -----------------------------------

logit <- function(x) -log((1/x) -1)
invLogit <- function(x) 1/(1+exp(-x))
leq <- function(x,a,b) 1/(1+exp(-(a+b*x)))
leq2 <- function(x,a1,b1,a2,b2){
  # exp(x*b1+x*b2+a1+a2)/(exp(a1+b1*x)+1)*(exp(a2+b2*x)+1)
  leq(x,a1,b1)*leq(x,a2,b2)
}

curve(invLogit,0,100)
curve(leq(x,-1,exp(-2)),0,100,ylim=c(0,1)) #Environment
curve(leq(x,2,-0.1),0,100,ylim=c(0,1)) #Insects
curve(leq(x,1,exp(-2))*leq(x,1,-0.1),0,100,ylim=c(0,1)) #Both
curve(leq2(x,1,exp(-2),1,-0.1),0,100,ylim=c(0,1)) #Both


png(paste0('./Figures/ExamplePlots/hypotheses.png'),width=10,height=10,units='in',res=200)

# png("C:\\Users\\Samuel\\Desktop\\hypothesisPlot.png",width=12,height=6,units='in',res=200) #Path on Multivac

savepar <- par(no.readonly = TRUE)
par(mfrow=c(2,2),mai=c(1,0.5,0.8,0.4))


#Scenario: no effect environmental effect
a1 <- 0; b1 <- 0
plot(0,0,type='n',xlim=c(0,100),ylim=c(0,1),axes=FALSE,
     ylab=NA,xlab='Distance from Edge',main='Scenario 1',xaxt='n')
axis(2,at = c(0,1),labels = c('Low','High')); axis(1,at=seq(0,100,by=20))
curve(leq(x,a1,b1)+0.05,0,100,xlim=range(dist),col='red',add=TRUE,lwd=2)
curve(leq(x,a1,b1)-0.05,0,100,xlim=range(dist),col='black',add=TRUE,lwd=2)
text(70,0.6,'Edge effect',col='red'); text(70,0.4,'Crop Yield')

#Scenario: negative environmental effect
a1 <- 0; b1 <- 0.1
plot(0,0,type='n',xlim=c(0,100),ylim=c(0,1),axes=FALSE,
     ylab=NA,xlab='Distance from Edge',main='Scenario 2',xaxt='n')
axis(2,at = c(0,1),labels = c('Low','High')); axis(1,at=seq(0,100,by=20))
curve(leq(x,a1,b1)-0.1,0,100,xlim=range(dist),col='red',add=TRUE,lwd=2)
curve(leq(x,a1,b1)-0.2,0,100,xlim=range(dist),col='black',add=TRUE,lwd=2)
text(70,0.95,'Edge effect',col='red'); text(70,0.7,'Crop Yield')

#Scenario: positive environmental effect
a1 <- 0; b1 <- -0.01
plot(0,0,type='n',xlim=c(0,100),ylim=c(0,1),axes=FALSE,
     ylab=NA,xlab='Distance from Edge',main='Scenario 3',xaxt='n')
axis(2,at = c(0,1),labels = c('Low','High')); axis(1,at=seq(0,100,by=20))

curve(leq(x,a1,b1)+0.2,0,100,xlim=range(dist),col='blue',add=TRUE,lwd=2)
curve(leq(x,a1,b1)+0.1,0,100,xlim=range(dist),col='black',add=TRUE,lwd=2)
text(70,0.65,'Ecosystem\nservice effect',col='blue'); text(70,0.35,'Crop Yield')


#Scenario 3 - beneficial bugs decrease with distance into field
# a1 <- -1; b1 <- 0.2; a2 <- 3; b2 <- -0.02
# 
# plot(0,0,type='n',xlim=c(0,100),ylim=c(0,1),axes=FALSE,
#      ylab=NA,xlab='Distance from Edge',main='Scenario 4',xaxt='n')
# axis(2,at = c(0,1),labels = c('Low','High')); axis(1,at=seq(0,100,by=20))
# 
# curve(leq2(x,a1,b1,a2,b2)*0.95,0,100,lwd=4,add=TRUE)#,main='Scenario 3: bugs - with dist')
# curve(leq(x,a1,b1),0,100,xlim=range(dist),col='red',add=TRUE,lwd=2)
# curve(leq(x,a2,b2),0,100,xlim=range(dist),col='blue',add=TRUE,lwd=2)
# text(10,0.83,'Ecosystem\nservice\neffect',col='blue'); text(70,0.95,'Edge effect',col='red'); text(80,0.6,'Crop Yield')

a1 <- 0; b1 <- 0.1; a2 <- 0; b2 <- -0.01

plot(0,0,type='n',xlim=c(0,100),ylim=c(0,1),axes=FALSE,
     ylab=NA,xlab='Distance from Edge',main='Scenario 4',xaxt='n')
axis(2,at = c(0,1),labels = c('Low','High')); axis(1,at=seq(0,100,by=20))

curve(leq2(x,a1,b1,a2,b2)+0.1,0,100,lwd=2,add=TRUE)#,main='Scenario 3: bugs - with dist')
curve(leq(x,a1,b1)-0.1,0,100,xlim=range(dist),col='red',add=TRUE,lwd=2)
curve(leq(x,a2,b2)+0.2,0,100,xlim=range(dist),col='blue',add=TRUE,lwd=2)
text(70,0.65,'Ecosystem\nservice effect',col='blue'); text(70,0.95,'Edge effect',col='red'); text(70,0.3,'Crop Yield')
arrows(25,0.3,25,0.43); text(25,0.2,'Intermediate yield\nincrease')

par(savepar)

dev.off()

#Approximately how many basis functions needed to represent this?
library(tidyverse)
library(mgcv)

a1 <- -2; b1 <- 0.2; a2 <- 3; b2 <- -0.03
scen1 <- data.frame(dist=1:100,yield=leq2(1:100,a1,b1,a2,b2)) %>% mutate(y=yield+rnorm(n(),0,0.1))
a2 <- 0; b2 <- 0.03
scen2 <- data.frame(dist=1:100,yield=leq2(1:100,a1,b1,a2,b2)) %>% mutate(y=yield+rnorm(n(),0,0.1))
a2 <- 1; b2 <- 0
scen3 <- data.frame(dist=1:100,yield=leq2(1:100,a1,b1,a2,b2)) %>% mutate(y=yield+rnorm(n(),0,0.1))

m1 <- gam(y~s(dist,k=5,bs='cs'),data=scen1) #Looks like about 5 is necessary
m2 <- gam(y~s(dist,k=5),data=scen2)
m3 <- gam(y~s(dist,k=5),data=scen3)
summary(m1)

scen1 %>% mutate(pred=predict(m1)) %>% 
  ggplot(aes(x=dist,y=y))+geom_point()+
  geom_line(aes(y=yield))+geom_line(aes(y=pred),col='red')

scen2 %>% mutate(pred=predict(m2)) %>% 
  ggplot(aes(x=dist,y=y))+geom_point()+
  geom_line(aes(y=yield))+geom_line(aes(y=pred),col='red')

scen3 %>% mutate(pred=predict(m3)) %>% 
  ggplot(aes(x=dist,y=y))+geom_point()+
  geom_line(aes(y=yield))+geom_line(aes(y=pred),col='red')

#Simulate data using logit-normal distribution --------------

a1 <- -2; b1 <- 0.2; a2 <- 3; b2 <- -0.03; sig <- 0.5
dist <- rep(1:100,each=3)

par(mfrow=c(1,1))
yhat <- leq2(dist,a1,b1,a2,b2)
y <- invLogit(logit(yhat)+rnorm(length(dist),0,sig))
plot(dist,yhat,type='l')
points(dist,y)

f1 <- function(p,d){ #log-likelihood function
  y <- d[,1] #Data
  dist <- d[,2] #Covariate
  yhat <- leq2(dist,p[1],p[2],p[3],p[4]) #Predicted
  -sum(dnorm(y,yhat,exp(p[5]),log=TRUE))
}

mod1 <- optim(par=unlist(pars),fn=f1,d=cbind(y,dist))

#Bootstrapping
boot <- lapply(1:1000,function(x){
  o <- optim(par=unlist(pars),fn=f1,d=cbind(y,dist)[sample(1:length(dist),length(dist),TRUE),])
})

#Parameter ranges
parRange <- apply(sapply(boot,function(x) x$par),1,function(l) quantile(l,c(0.1,0.5,0.9)))

#Correlation b/w parameters - really weird sets of correlations. Looks like Neil's Funnel.
pairs(t(sapply(boot,function(x) x$par)))

#Prediction ranges
predRange <- apply(sapply(boot,function(x) leq2(dist,x$par[1],x$par[2],x$par[3],x$par[4])),1,
           function(l) quantile(l,c(0.1,0.5,0.9)))

plot(dist,y,pch=19,cex=0.8,ylab='Yield')
# lines(dist,yhat,lwd=2) #Actual process
lines(dist,predRange[2,],col='red',lwd=3) #Recovered process
lines(dist,predRange[1,],col='red',lty=2);lines(dist,predRange[3,],col='red',lty=2)

par(mfrow=c(3,2))
pars <- list(a1=a1,b1=b1,a2=a2,b2=b2,logsig=log(sig))
for(i in 1:length(pars)){
  plot(c(1,2),c(pars[[i]],parRange[2,1]),ylim=range(c(pars[[i]],parRange[,1])),xlim=c(0,3),
       col=c('black','red'),pch=19,main=names(pars)[i],ylab='Estimate',xlab=NA)
  lines(c(2,2),parRange[c(1,3),1],col='red')
}; par(mfrow=c(1,1))


# Simulate data using beta distribution ----------------------------------------

logit <- function(x) -log((1/x) -1)
invLogit <- function(x) 1/(1+exp(-x))
# leq <- function(x,a,b) 1/(1+exp(-(a+b*x)))
leq2 <- function(x,a1,b1,a2,b2){
  # l <- (a1*a2)+(a1*b2+a2*b1)*x+b1*b2*x^2
  # invLogit(l)
  1/((1+exp(-(a1+b1*x)))*(1+exp(-(a2+b2*x))))
  # y <- (a1+b1*x)*(a2+b2*x)
  # 1/(1+exp(-y))
}

curve(invLogit(x),-5,5)
curve(leq2(x,-2,0.2,3,-0.02),0,100)


a1 <- -2; b1 <- 0.2; a2 <- 3; b2 <- -0.02; sig <- 0.1

#Reparameterize beta in terms of mean and sd
rbeta2 <- function(N,mean,sd){
  a <- (((1-mean)/sd^2)-(1/mean))*mean^2
  b <- a*((1/mean)-1)
  rbeta(N,a,b)
}
dbeta2 <- function(x,mean,sd,l=FALSE){
  a <- (((1-mean)/sd^2)-(1/mean))*mean^2
  b <- a*((1/mean)-1)
  dbeta(x,a,b,log=l)
}

#Parameters
a1 <- -2; b1 <- 0.2; a2 <- 3; b2 <- -0.02; sig <- 0.1
pars <- list(a1=a1,b1=b1,a2=a2,b2=b2,logsig=log(sig))
dist <- rep(1:100,each=3)

par(mfrow=c(1,1))
yhat <- leq2(dist,a1,b1,a2,b2)
y <- rbeta2(1:length(dist),yhat,sig)
plot(dist,yhat,type='l',ylim=c(0,1))
points(dist,y)

f1 <- function(p,d){ #log-likelihood function
  y <- d[,1] #Data
  dist <- d[,2] #Covariate
  yhat <- leq2(dist,p[1],p[2],p[3],p[4]) #Predicted
  -sum(dbeta2(y,yhat,exp(p[5]),l=TRUE))
}

f1(unlist(pars),cbind(y,dist)) #Works

mod1 <- optim(par=unlist(pars),fn=f1,d=cbind(y,dist))
data.frame(unlist(pars),mod1$par) #Works OK

#Bootstrapping
boot <- lapply(1:500,function(x){
  o <- optim(par=unlist(pars),fn=f1,d=cbind(y,dist)[sample(1:length(dist),length(dist),TRUE),])
})

#Parameter ranges
parRange <- apply(sapply(boot,function(x) x$par),1,function(l) quantile(l,c(0.1,0.5,0.9)))

#Correlation b/w parameters - intercepts & slopes correlated for each relationship set
pairs(t(sapply(boot,function(x) x$par)))

#Prediction ranges
predRange <- apply(sapply(boot,function(x) leq2(dist,x$par[1],x$par[2],x$par[3],x$par[4])),1,
                   function(l) quantile(l,c(0.1,0.5,0.9)))

plot(dist,y,pch=19,cex=0.8,ylab='Yield')
# lines(dist,yhat,lwd=2) #Actual process
lines(dist,predRange[2,],col='red',lwd=3) #Recovered process
lines(dist,predRange[1,],col='red',lty=2);lines(dist,predRange[3,],col='red',lty=2)

par(mfrow=c(3,2))
pars <- list(a1=a1,b1=b1,a2=a2,b2=b2,logsig=log(sig))
for(i in 1:length(pars)){
  plot(c(1,2),c(pars[[i]],parRange[2,1]),ylim=range(c(pars[[i]],parRange[,1])),xlim=c(0,3),
       col=c('black','red'),pch=19,main=names(pars)[i],ylab='Estimate',xlab=NA)
  lines(c(2,2),parRange[c(1,3),1],col='red')
}; par(mfrow=c(1,1))
