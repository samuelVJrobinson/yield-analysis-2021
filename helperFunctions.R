# Helper functions

seqGroup <- function(x,singles=TRUE){ #Turns sequential IDs into numbered groups.
  #Singles = TRUE, sequential values must be exactly 1 greater (gaps matter)
  #Singles = FALSE, sequential values must be greater (gaps don't matter)
  if(singles){
    xl <- (x-lag(x))!=1
  } else {
    xl <- (x-lag(x))<0
  }
  cumsum(ifelse(is.na(xl),FALSE,xl))+1
}

makePolys <- function(dat,width='w',dist='d',angle='a',backwards=FALSE){
  # Make polygons from width, dist, and angle measurements, centered on location from dat
  
  gType <- dat %>% st_geometry_type(FALSE)
  
  if(gType!='POINT') warning(paste('Input data type is',gType,'not POINT',sep=' '))
  
  rectFun <- function(x,y,w,d,a,b){
    #Function to create corners of rotated rectangle from:
    #   starting location (x,y),
    #   width (w), distance (d), and angle (a)
    #   rotate rectangles 180 degrees? (b)
    rotate <- ifelse(b,90,270)
    a <- (rotate-a)*pi/180 #90-angle in radians. 
    v <- c(x,y) #Starting point
    v1 <- c(d*cos(a),d*sin(a)) #Vectors to add together to get corners
    v2 <- c(-(w/2)*sin(a),(w/2)*cos(a))
    return(rbind(v+v2,v+v1+v2,v+v1-v2,v-v2,v+v2)) #Corners of rotated rectangle
  }
  datCRS <- st_crs(dat) #Coordinate system from dat
  #Apply rectFun to all rows in dat. Probably could be done through map or some other purrr related thing.
  xcoord <- st_coordinates(dat)[,1]
  ycoord  <- st_coordinates(dat)[,2]
  polys <- lapply(1:nrow(dat),function(i){
    r <- rectFun(x = xcoord[i], y = ycoord[i],  w = dat[[width]][i],  d = dat[[dist]][i],
                 a = dat[[angle]][i], b = backwards)
    st_polygon(list(r))
  })
  #Combine dat with new polygon geometry, and add original CRS
  dat2 <- st_sf(st_drop_geometry(dat),st_sfc(polys)) %>% st_set_crs(datCRS)
  return(dat2)
}

mergePoly <- function(dat,fList=NULL){
  #Function to merge together completely overlapping yield polygons (cover or perfect overlap)
  # takes a list of functions (fList)
  # columns not mentioned are dropped
  
  
  vList <- as.list(names(fList))
  # doesContain <- unique(st_contains(dat)) #Deals with identical overlap
  doesContain <- st_covers(dat) 
  #Find polygons that are overlapped by others
  remove <- unique(unlist(lapply(doesContain,function(x) if(length(x)>1) x[2:length(x)] else logical(0)))) 
  doesContain <- doesContain[-remove] #Remove polygons from list
  
  # cMat <- table(rep(1:length(doesContain),sapply(doesContain,length)),unlist(doesContain))
  # diag(cMat) <- 0
  # doesContain <- doesContain[colSums(cMat)[1:nrow(cMat)]==0] #Strips out inverse overlap cases
  # any(sapply(doesContain,length)>1) #Do any polygons overlap completely?
  
  newPolys <- lapply(doesContain,function(i){
    if(length(i)>1){ 
      p1 <- dat[i,] #Dataframe with overlapping polygons
      newData <- map2(vList,fList, ~ st_drop_geometry(p1) %>% summarize(across(all_of(.x),.y))) %>%   
        reduce(data.frame)
      rownames(newData) <- paste(i,collapse='.')
      newGeom <- st_union(p1) %>% st_geometry()
      return(list(data=newData,geom=newGeom))
    } else if(i>0) {
      newData <- dat[i,unlist(vList)] %>% st_drop_geometry() 
      newGeom <- dat[i,] %>% st_geometry()
      return(list(data=newData,geom=newGeom))
    }
  }) 
  
  dat2 <- do.call('rbind',lapply(newPolys,function(x) x$data)) %>%
    st_set_geometry(st_sfc(sapply(newPolys,function(x) x$geom),crs=st_crs(dat)))
  
  return(dat2)
}

#Yield conversions 
# Requires starting and ending units
convertYield <- function(x,inUnits=NULL,outUnits=NULL){
  
  if(identical(inUnits,outUnits)) return(x)
  
  #Convert to g/m2
  x <- switch(inUnits,
              tpha = x/100, #Tonnes/ha
              gpm2 = x #g per m2
              )
  
  #Convert to output units
  xOut <- switch(outUnits,
                 tpha = x*100,
                 gpm2 = x)
  
  return(xOut)
} 
