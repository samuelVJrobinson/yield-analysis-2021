#Script from Vega et al 2019 to remove outliers and "inliers" - https://link.springer.com/article/10.1007/s11119-018-09632-8

#install.packages("BiocParallel")## For install go to https://www.bioconductor.org/install/
#install.packages("sp")
#install.packages("spdep")

# only for download dataset for test and updates of function "paradepu"
url <- "https://drive.google.com/drive/folders/142nA5wW3CwUwjOOLSakgis8DqZDPM4I_?usp=sharing"
browseURL(url, browser = getOption("browser"),encodeIfNeeded = FALSE)

#questions please write to: andresvega@agro.unc.edu.ar

## to create object of the function paradepu
paradepu<-function(file.input,file.out,col.yiel=3,col.x=1,col.y=2,dis.buffer=30,times.sd=2.5,dis.neighbors=30,N_CPU=4){
  library("BiocParallel")
  carp<-list.files(file.input)
  
  depu<-function(x){ 
    file.input="C:/Users/veget/Desktop/for test/"
    files<-file.input
    library("sp")
    library("spdep")
    # basedatos<-read.table(paste0(file.input,"/",carp[x]),header=TRUE)
    basedatos <- read.table("C:\\Users\\Samuel\\Desktop\\dataset1.txt",header=TRUE)
    
    ##################### borde #############################
    redu<-function(mapa,discor=dis.buffer){
      require("ggplot2")
      require("raster")
      require("SDMTools")
      require("rgeos")
      borde<-mapa[chull(mapa[,1:2]),1:2]
      pol=Polygons(list(Polygon(cbind(borde[,1],
                                      borde[,2]))),as.character(dim(borde)[1]))
      sr=SpatialPolygons(list(pol))
      buff<-buffer(sr, width = -discor, dissolve = T)
      df.pl <- fortify(buff)[, 1:2]
      clapoli<-pnt.in.poly(mapa[,1:2],df.pl)
      sali<-cbind(mapa,"clapoli"=clapoli[,3])
      return(sali)
    }
    
    #Remove outliers (2.5 SD away from mean)
    sinoutli<-function(mydat,times_sd=times.sd,yielcol=col.yiel){
      nom<-c(colnames(mydat))
      vari<-nom[yielcol]
      Mean <- mean(mydat[,yielcol])
      stde <- sd(mydat[,yielcol])
      Lower <- Mean-times_sd*stde
      Upper <- Mean+times_sd*stde
      mydat[,yielcol][Upper<mydat[,yielcol] | mydat[,yielcol]<Lower] <-NA
      sinout <- mydat[which(!is.na(mydat[,yielcol])),]
      sinout
    }
    

    #Remove "inliers" (points that have bad spatial autocorrelation)
    sininli<-function(mydat,neighbor=dis.neighbors,yielcol=col.yiel){
      
      nom<-c(colnames(mydat))
      vari<-nom[yielcol]
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
          list(i,dnearneigh(cord,0,(dis.neighbors+i)))
        }
      }
      
      gri1<-list()
      for (i in c(0,10,20,30,40,50,60)) {
        B<-  tryCatch(ErrorNo4(i), error=function (e) {})
        if(is.null(B)){next
        } else {
          gri1<-B
          break}
      }
      gri <- gri1[[2]]
      lw <- nb2listw(gri, style = "W")
      # Local Morans Index
      LM <- localmoran(mydat[,yielcol],lw,p.adjust.method="bonferroni",alternative ="less")
      # Moran Plot
      MP<-moran.plot(mydat[,yielcol], lw,quiet=T,labels=F,zero.policy=F,xlab=vari, ylab="Spatially Lagged")
      # Identification of inliers
      Influ <- MP$is.inf
      mydat <- data.frame(mydat,LM,Influ)
      mydata1 <- subset(mydat,mydat[,dim(basedatos)[2]+1] > 0 | mydat[,dim(basedatos)[2]+4] > 0.05 )
      
      # Elimination of inliers identified with Moran Plot
      eoyeoi <- mydata1[mydata1$dfb.1_ == FALSE & mydata1$dfb.x == FALSE & mydata1$dffit == FALSE
                        & mydata1$cov.r == FALSE & mydata1$cook.d  == FALSE & mydata1$hat == FALSE, ]
      eoyeoi[,1:5]
    }
    
    base1<-redu(basedatos,dis.buffer)
    base<-base1[which(base1$clapoli==1),]
    out<- sinoutli(base,times_sd=times.sd,yielcol=5)
    final<-sininli(out,neighbor=dis.neighbors)
    final
    write.table(final,paste0(file.out,"/",carp[x]),row.names = FALSE)
  }
  
  #Run depu function (automatic cleaning) in parallel for each file in list (carp = file list)
  bplapply(1:length(carp),depu, 
           BPPARAM=SnowParam(workers=N_CPU, progressbar=TRUE, type="SOCK"))
} ## so far is the function paradepu

###
#file.input:folder where only the databases are downloaded, see example of how to write
#file.put: folder where you will save the cleaned databases, see example of how to write
#col.yiel: Yield column.
#col.x and col.y: coordinates columns only points UTM.
#dis.buffer: distance for remove border effects for defect is 30 .
#times.sd: To identiﬁed automatically and remove yield values that are outside the mean ± times.sd for defect is 2.5
#dis.neighbors: To calculate W_ij, neighbour points were deﬁned as contiguous points.
#N_CPU: It depends on the number of cores that your computer can have. 
#to know the number of cores of your computer see in resource monitor the exact number of CPU, in the link you can see an example
#for this it is recommended to use ncores-1
# "http://www.thewindowsclub.com/use-resource-monitor-windows-10"
### Run the function.
paradepu(file.input="C:/Users/veget/Desktop/for test/",file.out="C:/Users/veget/Desktop/new/",
         col.yiel=5,col.x=1,col.y=2,dis.buffer=30,times.sd=3,dis.neighbors=30,N_CPU=5)