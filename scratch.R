#Scratch space

library(tidyverse)
theme_set(theme_bw())
library(sf)

# #Test polygon set (only deals with overlap of 2, plus exact overlap)
# m <- rbind(c(1,0), c(1,1), c(0,1), c(0,0),c(1,0))
# p <- st_polygon(list(m))
# dat <- vector("list", 2)
# dat[[1]] <- p
# dat[[2]] <- p + c(0.65,0.75)
# dat[[3]] <- p + 3
# dat[[4]] <- p*0.25+0.25
# dat[[5]] <- dat[[3]] #5 and 6 exactly overlap 3
# dat[[6]] <- dat[[3]]
# dat[[7]] <- p + c(3,2)
# dat <- st_sf(speed=rep(1,length(dat)),yield=runif(length(dat),0,3),st_sfc(dat)) %>%
#   st_set_crs(3401)
# dat %>% ggplot()+geom_sf(aes(fill=yield),alpha=0.6)

#Test polygon set (deals with overlaps of 3)
m <- rbind(c(1,0), c(1,1), c(0,1), c(0,0),c(1,0))
p <- st_polygon(list(m))
dat <- vector("list", 2)
dat[[1]] <- p
dat[[2]] <- p + 0.75
dat[[3]] <- p + 3
dat[[4]] <- p*0.25+0.25
dat[[5]] <- p*0.5+0.2
dat[[6]] <- p+c(0.75,0)
dat[[7]] <- p*c(2,0.3) + c(2.5,3.5)
dat[[8]] <- p + c(2.1,3.1)
dat[[9]] <- p + c(3,1)
dat <- st_sf(speed=rep(1,length(dat)),yield=runif(length(dat),0,3),st_sfc(dat)) %>%
  st_set_crs(3401)
dat %>% ggplot()+geom_sf(aes(fill=yield),alpha=0.6)

# #Triangulate turns a polygon into a set of triangles. 
# dat2 <- rbind(c(1,0), c(1,1), c(0.5,1.5),c(0,1), c(0,0),c(0.5,-0.5),c(1,0))
# dat2 <- st_polygon(list(dat2))
# ggplot(dat2)+geom_sf(fill=NA)+
#   geom_sf(data=st_triangulate(dat2),col='red',fill=NA)

# #Intersecting + non-intersecting bits
# st_intersection(dat) %>% 
#   ggplot()+geom_sf(aes(fill=n.overlaps),alpha=0.5)+
#   facet_wrap(~n.overlaps)
# 
# #Non-intersecting bits
# st_difference(dat) %>% mutate(n=1:n()) %>% 
#   ggplot()+geom_sf(aes(fill=n),alpha=0.5)
# 
# #Turns all polygons into single multipolygon (keeps them distinct)
# st_combine(dat) %>% ggplot()+geom_sf(alpha=0.5)
# 
# #Turns all polygons into single multipolygon (merges ones that overlap)
# st_union(dat) %>% ggplot()+geom_sf(alpha=0.5) 
# 
# #Add buffer to polygons
# st_buffer(dat,dist=0.1) %>% ggplot()+geom_sf(alpha=0.5) +geom_sf(data=dat,col='red',alpha=0.3)
# 
# #Boundary of a geometry
# st_boundary(dat) %>% ggplot+geom_sf() #Returns the same thing. I think this more useful for multipoints
# 
# st_boundary(st_union(dat)) %>% ggplot+geom_sf() #Again, returns the same thing.
# 
# #Centroids for each polygon
# st_centroid(dat) %>% ggplot() + geom_sf(alpha=0.5)
# 
# #Triangulate (make triangles) between points
# st_triangulate(dat) %>% ggplot() + geom_sf(alpha=0.5)
# 
# st_intersection(dat) %>% filter(n.overlaps==2) %>% #Only between intersecting points
#   st_cast('MULTIPOINT') %>% 
#   st_triangulate() %>% ggplot() + geom_sf(alpha=0.5) +
#   geom_sf(data=st_difference(dat),fill=NA)
# 
# st_intersection(dat) %>% filter(n.overlaps==2) %>% #Only between intersecting points
#   st_triangulate() %>% st_collection_extract(type='POLYGON') %>% 
#   st_centroid() %>% 
#   ggplot() + geom_sf(alpha=0.5)
# 
# 
# #Not sure this will work for dividing up intersection areas. Perhaps voronoi polygons:
# int <- st_intersection(dat) %>% filter(n.overlaps==2) #Intersecting polygons
#   # st_triangulate() %>% st_collection_extract(type='POLYGON') %>% 
#   # st_centroid()


#Idea for merging polygons:
# 1. Merge polygons that completely overlap others

source('helperFunctions.R')

#Deals with complete overlapping cases
dat2 <- dat %>% 
  mutate(mass=as.numeric(yield*st_area(.))) %>% #Turn yield into raw mass (difficult to deal with areas inside of mergePoly)
  mergePoly(fList=lst(speed= mean, mass= sum)) %>% 
  mutate(yield=as.numeric(mass/st_area(.))) #Turn raw mass back into yield
# debugonce(mergePoly)
  
#Looks OK
dat2 %>% mutate(r=1:n()) %>% ggplot()+geom_sf(aes(fill=yield),alpha=0.5)

# 2. For intersecting areas, merge together overlapping ares
# Voronai polygons work in simple cases (2 polygons), but not in multiple, which is the majority of yield polygons
weldPoly <- function(dat,weldVars){
  #Function to "weld" overlapping polygons
  # weldVars = vars to spatially average
  # columns not mentioned are kept the same
  
  #Convenience function to calculate weighted sum for each "welded" polygon
  # x = variable of interest (vector) 1xN
  # A = Area of each polygon (vector) 1xN
  # Ad = Area of overlap, split between each polygon (vector) 1xN
  weldSum <- function(x,A,Ad){
    A <- as.numeric(A); Ad <- as.numeric(Ad)
    xd <- sum(x*(sum(Ad)/A)) #Variable in overlapping area
    x - (x* sum(Ad)/A) + xd*(Ad/sum(Ad)) #Mass in new polygon
  }

  # intPoly <- st_intersects(dat) #Polygons that intersect
  # intPoly <- unique(c(st_intersects(dat),st_overlaps(dat)))
  # intPoly <- intPoly[sapply(intPoly,length)!=0] #Get rid of zero overlap cases
  
  dat <- dat2 #Prototyping
  
  #All cases
  doesIntersect <- unique(st_intersects(dat))
  # 
  # st_intersects(dat)
  # st_disjoint(dat)
  # st_touches(dat)
  # st_crosses(dat)
  # st_within(dat)
  # st_contains(dat)
  # st_overlaps(dat)
  # st_covers(dat)
  
  
  # #Pairwise cases
  # intPoly <- unique(c(st_overlaps(dat))) #Get overlapping polygons
  # intPoly[sapply(intPoly,length)==0] <- c(1:length(intPoly))[sapply(intPoly,length)==0] #Get rid of zero overlap cases
  # intPoly <- intPoly[order(sapply(intPoly,length),decreasing=TRUE)] #Sort
  
  # for(i in 1:length(intPoly)){
  
  i <- 1 #Eventually inside a for/lapply loop
  
  polySet <- doesIntersect[[i]] #Set of polygons to use
  # if(length(polySet)>1){
  p1 <- dat[polySet,] #Dataframe with overlapping polygons
  
  p1 %>% rownames_to_column('poly') %>% 
    ggplot()+geom_sf(aes(col=poly),fill=NA)+
    geom_sf_text(aes(label=poly,col=poly))
      
  # #Calculate new geometry - pairwise
  # #This works only for pairwise overlap
  # pDiff <- st_sym_difference(st_geometry(p1[1,]), #Polygons with intersection clipped
  #                            st_geometry(p1[2,])) %>% st_cast('POLYGON')
  # pInt <- st_intersection(st_geometry(p1[1,]), #Only intersection
  #                         st_geometry(p1[2,]))
  
  # st_intersection(p1) %>% ggplot()+geom_sf(alpha=0.3)
      
  intPoly <- st_intersection(p1)# %>% st_geometry() #Polygon intersections
      
  p1 %>% ggplot()+geom_sf(alpha=0.3)
  
  #Decompose into list of polygons. Each level has polygons of 1,2,3 sets of overlap
  pList <- lapply(1:max(intPoly$n.overlaps),function(x){
    st_cast(st_cast(intPoly[intPoly$n.overlaps==x,],'MULTIPOLYGON'),'POLYGON') #%>% 
      # select(-n.overlaps,-contains('origins')) #%>% st_geometry()
  })
      
  for(j in 1:(length(pList)-1)){ #Remove difference from upper-level intersections (ie 3-overlap separate from 2, separate from 1)
    pList[[j]] <- st_difference(pList[[j]],st_union(do.call('rbind',pList[(j+1):length(pList)])))
      # select(-matches('\\.\\d'))
  }
      
  # Bits of the intersecting polygons
  intPoly <- do.call('rbind',pList) %>% mutate(area=as.numeric(st_area(.)))
  
  intPoly %>% #Subdivisions
    mutate(r=1:n()) %>%
    ggplot()+
    geom_sf(alpha=0.3)+
    geom_sf_label(aes(label=r))
  
  

  # #Original 3 polygons
  # p1 %>% mutate(r=1:n()) %>% 
  #   ggplot()+
  #   geom_sf(alpha=0.3)+
  #   geom_sf_label(aes(label=r))
  
  #Features to merge together
  st_intersection(intPoly,st_geometry(p1[1,])) %>%
    filter(grepl("POLYGON", st_geometry_type(.))) %>% 
    select(-matches('\\.\\d')) %>%
    mutate(area=as.numeric(st_area(.))) %>% 
    st_drop_geometry()
    # ggplot()+geom_sf(alpha=0.3)
  
  st_overlaps(intPoly,p1[1,])
  
  st_intersection(intPoly,p1[1,]) 
      
  intPoly %>% mutate(r=1:n()) %>% 
    ggplot()+geom_sf(alpha=0.5,fill='red')+
    geom_sf_text(aes(label=n.overlaps))+
    geom_sf(data=st_intersection(p1),fill=NA)+
    facet_wrap(~r)
      
      
      
      # st_difference(p1) %>% ggplot()+geom_sf(alpha=0.3)
      # 
      # st_difference(p1) %>% mutate(r=1:n()) %>% 
      #   ggplot()+geom_sf(alpha=0.3)+
      #   facet_wrap(~r)
      # 
      # 
      # st_intersection(p1) %>% #filter(n.overlaps>1) %>% 
      #   # bind_rows(st_difference(p1)) %>% 
      #   mutate(r=1:n()) %>% 
      #   ggplot()+geom_sf(alpha=0.3)+facet_wrap(~r)
      # # st_sym_difference %>% ggplot()+geom_sf(aes(fill=n.overlaps),alpha=0.3)
      # 
      # st_difference(p1) %>% st_cast('MULTIPOLYGON') %>% 
      #   st_cast('POLYGON') %>% 
      #   mutate(r=1:n()) %>% 
      #   ggplot()+geom_sf(alpha=0.3)+facet_wrap(~r)
      
      
      
      
      
      st_union(st_cast(st_boundary(p1),'POLYGON'))
      
      intBorder <- pInt %>% st_cast('LINESTRING') #Border of intersection
      intCenter <- pInt %>% st_centroid() #Center of intersection
      diffCenters <- st_geometry(pDiff) %>% st_centroid() #Centers of difference polygons
      intPoints <- st_union(diffCenters,intCenter) %>%
        st_cast('LINESTRING') %>% #Lines b/w center of intersection to center of differences
        st_intersection(st_cast(pDiff,'LINESTRING')) #Intersection with border of difference
      
      vPoly <- st_voronoi(st_union(intPoints),envelope=st_union(p1)) %>% st_cast()  #Voronoi polygons
      vPoly <- st_crop(vPoly,pInt) #Clip to intersecting box
      
      #Geometry of new polygons
      newGeom <- sapply(1:length(pDiff),function(j) st_geometry(st_union(vPoly[j],pDiff[j]))) %>%
        st_sfc() #%>% st_set_crs(st_crs(dat))
      
      pArea <- st_area(dat[polySet,]) #Area of original polygons
      pOverlap <- st_area(vPoly) #Area of overlap polygons
      
      #Associated data
      newData <- st_drop_geometry(p1) %>% #select({{weldVars}}) %>%
        mutate(across({{weldVars}},~weldSum(.x,A=pArea,Ad=pOverlap)))
      
      #Replace data from dat
      st_sfc(newData,newGeom)
      
    } 
  
  
  ## This doesn't work if doing pairwise joins (only works iteratively). Either figure out iterative process or do joins non-pairwise
    
  
  # #Join together list into spatial df 
  # datOut <- do.call('rbind',lapply(newPolys,function(x) x$data)) %>% 
  #   st_set_geometry(do.call('c',sapply(newPolys,function(x) x$geom))) %>% 
  #   st_set_crs(st_crs(dat))
  
  return(datOut)
}

#Deals with incomplete overlap
debugonce(weldPoly)
dat3 <- dat2 %>% weldPoly(mass)

#Works for case of only 2 overlapping polygons
dat3 %>% mutate(r=1:n()) %>% ggplot()+geom_sf(aes(fill=yield),alpha=0.5)


