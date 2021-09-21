#Script from Vega et al 2019 to remove outliers and "inliers" - https://link.springer.com/article/10.1007/s11119-018-09632-8

#Remove "inliers" (points that have bad spatial autocorrelation)
require(sf)
require(sp)
require(spdep)
require(tidyverse)
theme_set(theme_classic())

# mydat <- dat
basedatos <- read.table("C:\\Users\\Samuel\\Desktop\\dataset1.txt",header=TRUE)

mydat <- basedatos #Simple df

neighbor <- dis.neighbors <- 30
yielcol <- 'Yield'

cord <- coordinates(mydat[,1:2])

buscar<-function(x){ #Makes neighbourhood weights
  gri <- dnearneigh(cord,0,(dis.neighbors+x))
  strsplit(as.character(try(nb2listw(gri, style = "W"))),
           "[ ]")[[1]][1]=="Error"
}

ErrorNo4<-function (i) {
  if(buscar(i)==TRUE) {
    stop("no hay vecinos")
  } else {
    list(i,dnearneigh(cord,0,(dis.neighbors+i))) #Returns list with neighbour indices for each point
  }
}

gri1<-list()
for (i in c(0,10,20,30,40,50,60)) { #Looks like this tries various neighbourhoods until one of them works
  B <- tryCatch(ErrorNo4(i), error=function (e) {})
  if(is.null(B)){
    next
  } else {
    gri1<-B
    break
  }
}

gri <- gri1[[2]] #Select only the neighbourhood list
lw <- nb2listw(gri, style = "W") #Adds spatial weights to neighbourhood list; W = row-standardized: each neighbour given equal weight (i.e. sum neighbour weights for each point = 1, sum(lw$weights[[1]])=1)

# Local Morans Index for each neighbourhood - uses weights from lw above, and spits out 5 numbers for each neighbourhood:
# Ii = local (empirical) Moran statistic
# E.Ii = expected value
# Var.Ii = variance
# Z.Ii = Z-score
# Pr(z) = p-value
LM <- localmoran(mydat[,yielcol],lw,p.adjust.method="bonferroni",alternative ="less")

# Moran Plot of Yield data against spatially lagged values
MP <- moran.plot(mydat[,yielcol], lw,quiet=F,labels=TRUE,zero.policy=FALSE,
                 xlab=yielcol, ylab="Spatially Lagged")

# Identification of inliers
mydat <- data.frame(mydat,LM,MP) #Joins to dataframe

#Gets rid of points with a negative Local Moran's I (low yield point surrounded by high yield)

mydat$keepThese <- mydat$Ii > 0 | #Local Moran's I 
  mydat$Pr.z...E.Ii.. > 0.05 #Z-score, but I think this should be p-value instead
mydat$r <- 1:nrow(mydat)

mydat %>% 
  ggplot(aes(x=r,y=Yield))+
  geom_line()+
  geom_point(aes(col=keepThese),alpha=0.3)+
  coord_cartesian(
    # xlim=c(5000,5100),
                  ylim=c(0,30))
  
mydata1 <- filter(mydat,keepThese) 

nrow(mydat)-nrow(mydata1) #Chops out about 78 data points

# I have no idea what this part is supposed to do... all of these variables are in MP not mydata1, and are numeric not logical
# # Elimination of inliers identified with Moran Plot
# eoyeoi <- mydata1[mydata1$dfb.1_ == FALSE & mydata1$dfb.x == FALSE & mydata1$dffit == FALSE
#                   & mydata1$cov.r == FALSE & mydata1$cook.d  == FALSE & mydata1$hat == FALSE, ]
# eoyeoi[,1:5]

