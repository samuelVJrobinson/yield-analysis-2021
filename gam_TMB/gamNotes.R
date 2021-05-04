#Code from Ch.2

llm <- function(theta,X,Z,y){
  ##untransform parameters...
  sigma.b <- exp(theta[1])
  sigma <- exp(theta[2])
  
  ##extract dimensions
  n <- length(y) #Number of data points
  pr <- ncol(Z) #Number of random effects
  pf <- ncol(X) #Number of fixed effects
  
  ##obtain beta-hat, b-hat
  X1 <- cbind(X,Z) #Put fixef and ranef matrices together
  #Vector of 0s (fixed effects) + Vector of 1/sigma.b^2 (random effects)
  ipsi <- c(rep(0,pf),rep(1/sigma.b^2,pr)) 
  b1 <- solve(crossprod(X1)/sigma^2+diag(ipsi),t(X1) %*% y/sigma^2) #random effects coefficients
  
  ##compute log|Z'Z/sigma^2 + I/sigma.b^2
  ldet <- sum(log(diag(chol(crossprod(Z)/sigma^2 + diag(ipsi[-(1:pf)])))))
  
  ##compute log profile likelihood
  l <- (-sum((y-X1 %*%b1)^2)/sigma^2 - sum(b1^2*ipsi) - n*log(sigma^2) - pr*log(sigma.b^2) - 2*ldet -n*log(2*pi))/2
  attr(l,"b") <- as.numeric(b1) #random effects coefficients
  
  return(-l) ##return beta-hat (fixed effect coefs) and b-hat (random effect variance)
}

library(nlme)
options(contrasts=c('contr.treatment','contr.treatment'))
Z <- model.matrix(~ Rail$Rail -1) #Random effects model matrix
X <- matrix(1,18,1) #Fixed effects model matrix (Intercepts)

# llm(c(0,0),X=X,Z=Z,y=Rail$travel) #Test log-lik function
rail.mod <- optim(c(0,0),llm,hessian=TRUE,X=X,Z=Z,y=Rail$travel)
exp(rail.mod$par) ## variance components
solve(rail.mod$hessian) ##approximate covariance matrix for theta
attr(llm(rail.mod$par,X,Z,Rail$travel),"b") #Random effect

llm(rail.mod$par,X,Z,Rail$travel) #Log-likelihood at convergence
debugonce(llm)



#Code from Ch.4 of Wood 2016 GAMs

library(gamair)
data(engine); attach(engine)

plot(size,wear,xlab='Engine capacity',ylab='Wear index')

tf <- function(x,xj,j){
  #generate jth tent function from set defined by knots xj
  dj <- xj*0
  dj[j] <- 1
  approx(xj,dj,x)$y #THIS is the actual basis function. In this case, it's just a series of triangles
}

tf.X <- function(x,xj){
  #tent function basis matrix given data x
  # and knot sequence xj
  nk <- length(xj)
  n <- length(x)
  X <- matrix(NA,n,nk)
  for(j in 1:nk) X[,j] <- tf(x,xj,j)
  X
}

(sj <- seq(min(size),max(size),length=6)) #Generate knot positions along range of "size" data
(X <- tf.X(size,sj)) #Get model matrix -> each row is a set of weights for each knot
(b <- lm(wear~X-1)) #Fit model - get coefficients to multiply onto model matrix
s <- seq(min(size),max(size),length=200) #New data to fit - across range of "size" data
(Xp <- tf.X(s,sj)) #make prediction matrix - get "weights" for new data
image(x=1:ncol(Xp),y=1:nrow(Xp),t(Xp),ylim=c(nrow(Xp),1),xlab='Xp columns',ylab='Xp rows')
plot(size,wear,ylim=c(0,5))
lines(s,Xp %*% coef(b)) #Multiply coefficients by model matrix to get predicted values
for(i in 1:ncol(Xp)){ 
  lines(s,Xp[,i],col=palette()[i],lty=2) #Basis functions (dashed)
  lines(s,Xp[,i]*coef(b)[i],col=palette()[i],lty=1) #Basis function * coef (solid)
}


#This works, but how to choose knot number?
# Could do backwards model selection, but this doesn't work because models aren't nested, perform poorly, and depend heavily on knot position

#Better method: penalized smoothing

prs.fit <- function(y,x,xj,sp){
  #y = y-value, thing we want to predict
  #x = x-value, predictor
  #xj = knot positions
  #sp = smoothing penalty - fixed for now
  
  (X <- tf.X(x,xj)) #model matrix - from tent function basis 
  (D <- diff(diag(length(xj)),differences=2)) #penalty matrix 
  # This penalizes differences in coefficients that are directly next to each other
  # eg. coef1 - 2*coef2 + coef3: if coef1 and 2 are very different from each other, penalty is large
  (D <- sqrt(sp)*D) #Multiplies sqrt of sp (penalty term) by penalty matrix
  
  #Uses augmented model matrix - this is a way to shrink X %*% coefs close to Y, and penalty close to 0
  
  (X <- rbind(X,D)) #augmented model matrix - binds penalty matrix _below_ model matrix
  (y <- c(y,rep(0,nrow(D)))) # augmented data - binds column of zeros _below_ data
  
  (lm(y ~ X-1)) #penalized least square fit - fit using LM
}

(sj <- seq(min(size),max(size),length=20)) #Knot positions
# debugonce(prs.fit)
(b <- prs.fit(y=wear,x=size,xj=sj,sp=.1)) #Coefficients for smoothing matrix, using penalty sp = 0.1
plot(size,wear)
Xp <- tf.X(s,sj) #Get model matrix for points s, given knots sj
lines(s,Xp %*% coef(b)) #Looks pretty good, actually
points(sj,tf.X(sj,sj) %*% coef(b),pch=19) #knot points

#This looks pretty good, but how to choose sp (penalty term)?
#Better method: generalized cross-validation - similar to leave-one-out 

#Direct search for GCV optimal smoothing parameter
rho <- seq(-9,11,length=90) #sp parameters to try out (log-space)
n <- length(wear)
V <- rep(NA,90) #Storage for GCV score
for(i in 1:90){
  b <- prs.fit(wear,size,sj,exp(rho[i])) #fit model
  trF <- sum(influence(b)$hat[1:n]) #extract EDF
  rss <- sum((wear-fitted(b)[1:n])^2) #residual SS
  V[i] <- n*rss/(n-trF)^2 #GCV score
}
plot(rho,V,type='l',xlab=expression(log(lambda)),main='GCV score')
(sp <- exp(rho[V==min(V)])) #optimal sp term, 18.3 = exp(2.9)
b <- prs.fit(wear,size,sj,sp) #refit model using new sp term
plot(size,wear,main='GCV optimal fit')
lines(s,Xp %*% coef(b))
(b <- prs.fit(y=wear,x=size,xj=sj,sp=0.1)) #Older coefficients (where sp=0.1)
lines(s,Xp %*% coef(b),col='red') #Very similar to optimal model


#Smoothing penalties are introduced because we believe that smoothness is more likely than wigglyness
#This suggests a Bayesian form for the smoother

X0 <- tf.X(size,sj) #X in original parameterization
D <- rbind(0,0,diff(diag(20),difference=2)) #Penalty matrix
diag(D) <- 1 #Augmented D - change diagonal elements to 1
X <- t(backsolve(t(D),t(X0))) #re-parameterized X
Z <- X[,-c(1,2)] #Mixed model matrices
X <- X[,1:2]

#estimate smoothing and variance parameters
m <- optim(c(0,0),llm,method="BFGS",X=X,Z=Z,y=wear)
exp(m$par) #random effect sigma, and overall sigma
b <- attr(llm(m$par,X,Z,wear),"b") #Get coefficients

#plot results
plot(size,wear)
Xp1 <- t(backsolve(t(D),t(Xp))) ##re-parameterized prediction matrix
lines(s,Xp1 %*% as.numeric(b),col='grey') #ML estimate

#Same thing using lme
library(nlme)
g <- factor(rep(1,nrow(X)))  #dummy factor (all ones)
m <- lme(wear~X-1, random = list(g =pdIdent(~Z-1)))
lines(s,Xp1 %*% as.numeric(coef(m))) #REML estimate - more variable than ML
