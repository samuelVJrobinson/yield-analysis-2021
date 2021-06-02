#GAM testing using TMB

library(TMB)
library(mgcv)
library(tidyverse)
library(gamair)
data(engine)

#Cubic splines, as in example ------------------------

s1 <- smoothCon(object=s(wear,bs='cr'),data=engine)

s1[[1]]$X #Model matrix
s1[[1]]$S #Penalty matrix

# Thin-plate splines ------------------------------------------------------

s1 <- smoothCon(object=s(wear,bs='cs'),data=engine)

s1[[1]]$X #Model matrix
s1[[1]]$S #Penalty matrix #Two zero columns. Not sure this is going to work

# Thin-plate splines ------------------------------------------------------

s1 <- smoothCon(object=s(wear,bs='tp'),data=engine)

s1[[1]]$X #Model matrix
s1[[1]]$S #Penalty matrix #Two zero columns. Not sure this is going to work
