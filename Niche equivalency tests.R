
library(ecospat)   
library(raster)     
library(SDMTools)   
library(dismo)      
library(phyloclim)
library(ade4)


###########################################################################################
###### PAIRWISE COMPARISONS WITH A SINGLE COMBINED "M" (SYMPATRIC DISTRIBUTION) ###########
###########################################################################################

############################################################
########################## Colinus ##########################
############################################################

setwd("WORKING DIRECTORY")

# Se leen los puntos de presencia de las  especies (previamente depurados, filtrados espacial y/o ambientalmente, etc.)
sp <- read.csv("~/Colinus_nigrogularis.csv", header=TRUE)
ev <- read.csv("~/YP_VS.csv", header=TRUE)

# Variables climaticas que se usaran para generar los componentes principales
# Las capas climaticas debieron ser previamente recortadas solo para el o las areas geograficas relevante para la comparacion de los nichos ecologicos, en este caso Norteamerica.
setwd("~/colinus_yp/asc")
varclim <- stack(list.files(pattern = "*.asc$",full.names = T))
resol <- res(varclim)[1]

climpunto <- rasterToPoints(varclim[[1]], fun=NULL, spatial=TRUE)
clim <- raster::extract(varclim, climpunto)

clim <- data.frame(coordinates(climpunto),clim)
clim <- subset(clim, !is.na(bio3) & !is.na(bio8) & !is.na(bio9) & !is.na(bio11) & !is.na(bio18))

occ.sp1 <- sp[2:3]
occ.sp2 <- ev[2:3]

#Anadir variables climaticas a datos
occ_sp1<-na.exclude(ecospat.sample.envar(dfsp=occ.sp1,colspxy=1:2, 
                                         colspkept=1:2,dfvar=clim, colvarxy=1:2,colvar="all",resolution=0.0416666))

occ_sp2<-na.exclude(ecospat.sample.envar(dfsp=occ.sp2,colspxy=1:2, 
                                         colspkept=1:2,dfvar=clim, colvarxy=1:2,colvar="all",resolution=0.0416666))

# N?mero de iteraciones
iterations<-1000
R=1000

########################################################################
############################ PCA-AMBIENTE ##############################
########################################################################
data <-rbind(clim[,3:7],occ_sp1[,3:7],occ_sp2[,3:7]) 
data <- subset(data, !is.na(bio3) & !is.na(bio8) & !is.na(bio9) & !is.na(bio11) & !is.na(bio18))

# vector de peso 0 para las ocurrencias y 1 para todos los sitios del area de estudio
w<-c(rep(1,nrow(clim)),rep(0,nrow(occ_sp1)),rep(0,nrow(occ_sp2))) 

pca.cal <-dudi.pca(data, row.w = w, center = T, scale = T, scannf = F, nf = 2)

# Filas en que estan los datos de clim y de cada una de las especies:
row.clim <- 1:nrow(clim)
row.sp1 <-  (1+nrow(clim)):(nrow(clim) + nrow(occ_sp1))
row.sp2 <-  (1 + nrow(clim) + nrow(occ_sp1)) : (nrow(clim) + nrow(occ_sp1) + nrow(occ_sp2))

scores.clim <- pca.cal$li[row.clim,] 
scores.sp1  <- pca.cal$li[row.sp1,]   
scores.sp2  <- pca.cal$li[row.sp2,]   

#Contribucion de cada variable a cada componente del PCA
ecospat.plot.contrib(contrib=pca.cal$co, eigen=pca.cal$eig)


########################################################################
########### SUPERFICIE DE DENSIDAD DE OCURRENCIAS ######################
########################################################################
z1<-ecospat.grid.clim.dyn (scores.clim, scores.clim, scores.sp1, R=200) 
z2<-ecospat.grid.clim.dyn (scores.clim, scores.clim, scores.sp2, R=200) 

#Metricos de sobrelape de nicho observado (D - el metrico de Schoener e I - el metrico de Warren)
ecospat.niche.overlap (z1=z1, z2=z2, cor=TRUE)

########################################################################
################### TEST DE EQUIVALENCIA DE NICHO ######################
########################################################################
a.dyn<-ecospat.niche.equivalency.test(z1=z1 , z2=z2, rep=1000) 

#GRAFICO PARA EL TEST DE EQUIVALENCIA DE NICHO
ecospat.plot.overlap.test(a.dyn,"D","Equivalency Colinus")





############################################################
########################## Campylorhynchus ##########################
############################################################

setwd("WORKING DIRECTORY")

# Se leen los puntos de presencia de las  especies (previamente depurados, filtrados espacial y/o ambientalmente, etc.)
sp <- read.csv("~/Campylorhynchus_yucatanicus.csv", header=TRUE)
ev <- read.csv("~/YP_VS.csv", header=TRUE)

# Variables climaticas que se usaran para generar los componentes principales
# Las capas climaticas debieron ser previamente recortadas solo para el o las areas geograficas relevante para la comparacion de los nichos ecologicos, en este caso Norteamerica.
setwd("~/campylorhynchus_yp/asc")
varclim <- stack(list.files(pattern = "*.asc$",full.names = T))
resol <- res(varclim)[1]

climpunto <- rasterToPoints(varclim[[1]], fun=NULL, spatial=TRUE)
clim <- raster::extract(varclim, climpunto)

clim <- data.frame(coordinates(climpunto),clim)
clim <- subset(clim, !is.na(bio8) & !is.na(bio14) & !is.na(bio18) & !is.na(bio19))

occ.sp1 <- sp[2:3]
occ.sp2 <- ev[2:3]

#Anadir variables climaticas a datos
occ_sp1<-na.exclude(ecospat.sample.envar(dfsp=occ.sp1,colspxy=1:2, 
                                         colspkept=1:2,dfvar=clim, colvarxy=1:2,colvar="all",resolution=0.0416666))

occ_sp2<-na.exclude(ecospat.sample.envar(dfsp=occ.sp2,colspxy=1:2, 
                                         colspkept=1:2,dfvar=clim, colvarxy=1:2,colvar="all",resolution=0.0416666))

# N?mero de iteraciones
iterations<-1000
R=1000

########################################################################
############################ PCA-AMBIENTE ##############################
########################################################################
data <-rbind(clim[,3:6],occ_sp1[,3:6],occ_sp2[,3:6]) 
data <- subset(data, !is.na(bio8) & !is.na(bio14) & !is.na(bio18) & !is.na(bio19))

# vector de peso 0 para las ocurrencias y 1 para todos los sitios del area de estudio
w<-c(rep(1,nrow(clim)),rep(0,nrow(occ_sp1)),rep(0,nrow(occ_sp2))) 

pca.cal <-dudi.pca(data, row.w = w, center = T, scale = T, scannf = F, nf = 2)

# Filas en que estan los datos de clim y de cada una de las especies:
row.clim <- 1:nrow(clim)
row.sp1 <-  (1+nrow(clim)):(nrow(clim) + nrow(occ_sp1))
row.sp2 <-  (1 + nrow(clim) + nrow(occ_sp1)) : (nrow(clim) + nrow(occ_sp1) + nrow(occ_sp2))

scores.clim <- pca.cal$li[row.clim,] 
scores.sp1  <- pca.cal$li[row.sp1,]   
scores.sp2  <- pca.cal$li[row.sp2,]   

#Contribucion de cada variable a cada componente del PCA
ecospat.plot.contrib(contrib=pca.cal$co, eigen=pca.cal$eig)


########################################################################
########### SUPERFICIE DE DENSIDAD DE OCURRENCIAS ######################
########################################################################
z1<-ecospat.grid.clim.dyn (scores.clim, scores.clim, scores.sp1, R=200) 
z2<-ecospat.grid.clim.dyn (scores.clim, scores.clim, scores.sp2, R=200) 

#Metricos de sobrelape de nicho observado (D - el metrico de Schoener e I - el metrico de Warren)
ecospat.niche.overlap (z1=z1, z2=z2, cor=TRUE)

########################################################################
################### TEST DE EQUIVALENCIA DE NICHO ######################
########################################################################
a.dyn<-ecospat.niche.equivalency.test(z1=z1 , z2=z2, rep=1000) 

#GRAFICO PARA EL TEST DE EQUIVALENCIA DE NICHO
ecospat.plot.overlap.test(a.dyn,"D","Equivalency Campylorhynchus")




############################################################
########################## Geococcyx ##########################
############################################################

setwd("WORKING DIRECTORY")

# Se leen los puntos de presencia de las  especies (previamente depurados, filtrados espacial y/o ambientalmente, etc.)
sp <- read.csv("~/Geococcyx_velox.csv", header=TRUE)
ev <- read.csv("~/YP_MPS_VS.csv", header=TRUE)

# Variables climaticas que se usaran para generar los componentes principales
# Las capas climaticas debieron ser previamente recortadas solo para el o las areas geograficas relevante para la comparacion de los nichos ecologicos, en este caso Norteamerica.
setwd("~/geococcyx_ypmps/asc")
varclim <- stack(list.files(pattern = "*.asc$",full.names = T))
resol <- res(varclim)[1]

climpunto <- rasterToPoints(varclim[[1]], fun=NULL, spatial=TRUE)
clim <- raster::extract(varclim, climpunto)

clim <- data.frame(coordinates(climpunto),clim)
clim <- subset(clim, !is.na(bio2) & !is.na(bio3) & !is.na(bio9) & !is.na(bio14) & !is.na(bio15) & !is.na(bio18) & !is.na(bio19))

occ.sp1 <- sp[2:3]
occ.sp2 <- ev[2:3]

#Anadir variables climaticas a datos
occ_sp1<-na.exclude(ecospat.sample.envar(dfsp=occ.sp1,colspxy=1:2, 
                                         colspkept=1:2,dfvar=clim, colvarxy=1:2,colvar="all",resolution=0.0416666))

occ_sp2<-na.exclude(ecospat.sample.envar(dfsp=occ.sp2,colspxy=1:2, 
                                         colspkept=1:2,dfvar=clim, colvarxy=1:2,colvar="all",resolution=0.0416666))

# N?mero de iteraciones
iterations<-1000
R=1000

########################################################################
############################ PCA-AMBIENTE ##############################
########################################################################
data <-rbind(clim[,3:9],occ_sp1[,3:9],occ_sp2[,3:9]) 
data <- subset(data, !is.na(bio2) & !is.na(bio3) & !is.na(bio9) & !is.na(bio14) & !is.na(bio15) & !is.na(bio18) & !is.na(bio19))

# vector de peso 0 para las ocurrencias y 1 para todos los sitios del area de estudio
w<-c(rep(1,nrow(clim)),rep(0,nrow(occ_sp1)),rep(0,nrow(occ_sp2))) 

pca.cal <-dudi.pca(data, row.w = w, center = T, scale = T, scannf = F, nf = 2)

# Filas en que estan los datos de clim y de cada una de las especies:
row.clim <- 1:nrow(clim)
row.sp1 <-  (1+nrow(clim)):(nrow(clim) + nrow(occ_sp1))
row.sp2 <-  (1 + nrow(clim) + nrow(occ_sp1)) : (nrow(clim) + nrow(occ_sp1) + nrow(occ_sp2))

scores.clim <- pca.cal$li[row.clim,] 
scores.sp1  <- pca.cal$li[row.sp1,]   
scores.sp2  <- pca.cal$li[row.sp2,]   

#Contribucion de cada variable a cada componente del PCA
ecospat.plot.contrib(contrib=pca.cal$co, eigen=pca.cal$eig)


########################################################################
########### SUPERFICIE DE DENSIDAD DE OCURRENCIAS ######################
########################################################################
z1<-ecospat.grid.clim.dyn (scores.clim, scores.clim, scores.sp1, R=200) 
z2<-ecospat.grid.clim.dyn (scores.clim, scores.clim, scores.sp2, R=200) 

#Metricos de sobrelape de nicho observado (D - el metrico de Schoener e I - el metrico de Warren)
ecospat.niche.overlap (z1=z1, z2=z2, cor=TRUE)

########################################################################
################### TEST DE EQUIVALENCIA DE NICHO ######################
########################################################################
a.dyn<-ecospat.niche.equivalency.test(z1=z1 , z2=z2, rep=1000) 

#GRAFICO PARA EL TEST DE EQUIVALENCIA DE NICHO
ecospat.plot.overlap.test(a.dyn,"D","Equivalency Geococcyx")




############################################################
########################## Amazilia ##########################
############################################################

setwd("WORKING DIRECTORY")

# Se leen los puntos de presencia de las  especies (previamente depurados, filtrados espacial y/o ambientalmente, etc.)
sp <- read.csv("~/Amazilia_rutila.csv", header=TRUE)
ev <- read.csv("~/YP_MPS_VS.csv", header=TRUE)

# Variables climaticas que se usaran para generar los componentes principales
# Las capas climaticas debieron ser previamente recortadas solo para el o las areas geograficas relevante para la comparacion de los nichos ecologicos, en este caso Norteamerica.
setwd("~/amazilia_ypmps/asc")
varclim <- stack(list.files(pattern = "*.asc$",full.names = T))
resol <- res(varclim)[1]

climpunto <- rasterToPoints(varclim[[1]], fun=NULL, spatial=TRUE)
clim <- raster::extract(varclim, climpunto)

clim <- data.frame(coordinates(climpunto),clim)
clim <- subset(clim, !is.na(bio2) & !is.na(bio3) & !is.na(bio9) & !is.na(bio13) & !is.na(bio14) & !is.na(bio15) & !is.na(bio18) & !is.na(bio19))

occ.sp1 <- sp[2:3]
occ.sp2 <- ev[2:3]

#Anadir variables climaticas a datos
occ_sp1<-na.exclude(ecospat.sample.envar(dfsp=occ.sp1,colspxy=1:2, 
                                         colspkept=1:2,dfvar=clim, colvarxy=1:2,colvar="all",resolution=0.0416666))

occ_sp2<-na.exclude(ecospat.sample.envar(dfsp=occ.sp2,colspxy=1:2, 
                                         colspkept=1:2,dfvar=clim, colvarxy=1:2,colvar="all",resolution=0.0416666))

# N?mero de iteraciones
iterations<-1000
R=1000

########################################################################
############################ PCA-AMBIENTE ##############################
########################################################################
data <-rbind(clim[,3:10],occ_sp1[,3:10],occ_sp2[,3:10]) 
data <- subset(data, !is.na(bio2) & !is.na(bio3) & !is.na(bio9) & !is.na(bio13) & !is.na(bio14) & !is.na(bio15) & !is.na(bio18) & !is.na(bio19))

# vector de peso 0 para las ocurrencias y 1 para todos los sitios del area de estudio
w<-c(rep(1,nrow(clim)),rep(0,nrow(occ_sp1)),rep(0,nrow(occ_sp2))) 

pca.cal <-dudi.pca(data, row.w = w, center = T, scale = T, scannf = F, nf = 2)

# Filas en que estan los datos de clim y de cada una de las especies:
row.clim <- 1:nrow(clim)
row.sp1 <-  (1+nrow(clim)):(nrow(clim) + nrow(occ_sp1))
row.sp2 <-  (1 + nrow(clim) + nrow(occ_sp1)) : (nrow(clim) + nrow(occ_sp1) + nrow(occ_sp2))

scores.clim <- pca.cal$li[row.clim,] 
scores.sp1  <- pca.cal$li[row.sp1,]   
scores.sp2  <- pca.cal$li[row.sp2,]   

#Contribucion de cada variable a cada componente del PCA
ecospat.plot.contrib(contrib=pca.cal$co, eigen=pca.cal$eig)


########################################################################
########### SUPERFICIE DE DENSIDAD DE OCURRENCIAS ######################
########################################################################
z1<-ecospat.grid.clim.dyn (scores.clim, scores.clim, scores.sp1, R=200) 
z2<-ecospat.grid.clim.dyn (scores.clim, scores.clim, scores.sp2, R=200) 

#Metricos de sobrelape de nicho observado (D - el metrico de Schoener e I - el metrico de Warren)
ecospat.niche.overlap (z1=z1, z2=z2, cor=TRUE)

########################################################################
################### TEST DE EQUIVALENCIA DE NICHO ######################
########################################################################
a.dyn<-ecospat.niche.equivalency.test(z1=z1 , z2=z2, rep=1000) 

#GRAFICO PARA EL TEST DE EQUIVALENCIA DE NICHO
ecospat.plot.overlap.test(a.dyn,"D","Equivalency Amazilia")









############################################################
########################## Forpus ##########################
############################################################

setwd("WORKING DIRECTORY")

# Se leen los puntos de presencia de las  especies (previamente depurados, filtrados espacial y/o ambientalmente, etc.)
sp <- read.csv("~/Forpus_mexicanus.csv", header=TRUE)
ev <- read.csv("~/MPS_VS.csv", header=TRUE)

# Variables climaticas que se usaran para generar los componentes principales
# Las capas climaticas debieron ser previamente recortadas solo para el o las areas geograficas relevante para la comparacion de los nichos ecologicos, en este caso Norteamerica.
setwd("~/forpus_mps/asc")
varclim <- stack(list.files(pattern = "*.asc$",full.names = T))
resol <- res(varclim)[1]

climpunto <- rasterToPoints(varclim[[1]], fun=NULL, spatial=TRUE)
clim <- raster::extract(varclim, climpunto)

clim <- data.frame(coordinates(climpunto),clim)
clim <- subset(clim, !is.na(bio2) & !is.na(bio3) & !is.na(bio14) & !is.na(bio15) & !is.na(bio18) & !is.na(bio19))

occ.sp1 <- sp[2:3]
occ.sp2 <- ev[2:3]

#Anadir variables climaticas a datos
occ_sp1<-na.exclude(ecospat.sample.envar(dfsp=occ.sp1,colspxy=1:2, 
                                         colspkept=1:2,dfvar=clim, colvarxy=1:2,colvar="all",resolution=0.0416666))

occ_sp2<-na.exclude(ecospat.sample.envar(dfsp=occ.sp2,colspxy=1:2, 
                                         colspkept=1:2,dfvar=clim, colvarxy=1:2,colvar="all",resolution=0.0416666))

# N?mero de iteraciones
iterations<-1000
R=1000

########################################################################
############################ PCA-AMBIENTE ##############################
########################################################################
data <-rbind(clim[,3:8],occ_sp1[,3:8],occ_sp2[,3:8])
data <- subset(data, !is.na(bio2) & !is.na(bio3) & !is.na(bio14) & !is.na(bio15) & !is.na(bio18) & !is.na(bio19))

# vector de peso 0 para las ocurrencias y 1 para todos los sitios del area de estudio
w<-c(rep(1,nrow(clim)),rep(0,nrow(occ_sp1)),rep(0,nrow(occ_sp2))) 

pca.cal <-dudi.pca(data, row.w = w, center = T, scale = T, scannf = F, nf = 2)

# Filas en que estan los datos de clim y de cada una de las especies:
row.clim <- 1:nrow(clim)
row.sp1 <-  (1+nrow(clim)):(nrow(clim) + nrow(occ_sp1))
row.sp2 <-  (1 + nrow(clim) + nrow(occ_sp1)) : (nrow(clim) + nrow(occ_sp1) + nrow(occ_sp2))

scores.clim <- pca.cal$li[row.clim,] 
scores.sp1  <- pca.cal$li[row.sp1,]   
scores.sp2  <- pca.cal$li[row.sp2,]   

#Contribucion de cada variable a cada componente del PCA
ecospat.plot.contrib(contrib=pca.cal$co, eigen=pca.cal$eig)


########################################################################
########### SUPERFICIE DE DENSIDAD DE OCURRENCIAS ######################
########################################################################
z1<-ecospat.grid.clim.dyn (scores.clim, scores.clim, scores.sp1, R=200) 
z2<-ecospat.grid.clim.dyn (scores.clim, scores.clim, scores.sp2, R=200) 

#Metricos de sobrelape de nicho observado (D - el metrico de Schoener e I - el metrico de Warren)
ecospat.niche.overlap (z1=z1, z2=z2, cor=TRUE)

########################################################################
################### TEST DE EQUIVALENCIA DE NICHO ######################
########################################################################
a.dyn<-ecospat.niche.equivalency.test(z1=z1 , z2=z2, rep=1000) 

#GRAFICO PARA EL TEST DE EQUIVALENCIA DE NICHO
ecospat.plot.overlap.test(a.dyn,"D","Equivalency Forpus")





############################################################
########################## Momotus ##########################
############################################################

setwd("WORKING DIRECTORY")

# Se leen los puntos de presencia de las  especies (previamente depurados, filtrados espacial y/o ambientalmente, etc.)
sp <- read.csv("~/Momotus_mexicanus.csv", header=TRUE)
ev <- read.csv("~/MPS_VS.csv", header=TRUE)

# Variables climaticas que se usaran para generar los componentes principales
# Las capas climaticas debieron ser previamente recortadas solo para el o las areas geograficas relevante para la comparacion de los nichos ecologicos, en este caso Norteamerica.
setwd("~/momotus_mps/asc")
varclim <- stack(list.files(pattern = "*.asc$",full.names = T))
resol <- res(varclim)[1]

climpunto <- rasterToPoints(varclim[[1]], fun=NULL, spatial=TRUE)
clim <- raster::extract(varclim, climpunto)

clim <- data.frame(coordinates(climpunto),clim)
clim <- subset(clim, !is.na(bio2) & !is.na(bio3) & !is.na(bio4) & !is.na(bio14) & !is.na(bio15) & !is.na(bio18) & !is.na(bio19))

occ.sp1 <- sp[2:3]
occ.sp2 <- ev[2:3]

#Anadir variables climaticas a datos
occ_sp1<-na.exclude(ecospat.sample.envar(dfsp=occ.sp1,colspxy=1:2, 
                                         colspkept=1:2,dfvar=clim, colvarxy=1:2,colvar="all",resolution=0.0416666))

occ_sp2<-na.exclude(ecospat.sample.envar(dfsp=occ.sp2,colspxy=1:2, 
                                         colspkept=1:2,dfvar=clim, colvarxy=1:2,colvar="all",resolution=0.0416666))

# N?mero de iteraciones
iterations<-1000
R=1000

########################################################################
############################ PCA-AMBIENTE ##############################
########################################################################
data <-rbind(clim[,3:9],occ_sp1[,3:9],occ_sp2[,3:9])
data <- subset(data, !is.na(bio2) & !is.na(bio3) & !is.na(bio4) & !is.na(bio14) & !is.na(bio15) & !is.na(bio18) & !is.na(bio19))

# vector de peso 0 para las ocurrencias y 1 para todos los sitios del area de estudio
w<-c(rep(1,nrow(clim)),rep(0,nrow(occ_sp1)),rep(0,nrow(occ_sp2))) 

pca.cal <-dudi.pca(data, row.w = w, center = T, scale = T, scannf = F, nf = 2)

# Filas en que estan los datos de clim y de cada una de las especies:
row.clim <- 1:nrow(clim)
row.sp1 <-  (1+nrow(clim)):(nrow(clim) + nrow(occ_sp1))
row.sp2 <-  (1 + nrow(clim) + nrow(occ_sp1)) : (nrow(clim) + nrow(occ_sp1) + nrow(occ_sp2))

scores.clim <- pca.cal$li[row.clim,] 
scores.sp1  <- pca.cal$li[row.sp1,]   
scores.sp2  <- pca.cal$li[row.sp2,]   

#Contribucion de cada variable a cada componente del PCA
ecospat.plot.contrib(contrib=pca.cal$co, eigen=pca.cal$eig)


########################################################################
########### SUPERFICIE DE DENSIDAD DE OCURRENCIAS ######################
########################################################################
z1<-ecospat.grid.clim.dyn (scores.clim, scores.clim, scores.sp1, R=200) 
z2<-ecospat.grid.clim.dyn (scores.clim, scores.clim, scores.sp2, R=200) 

#Metricos de sobrelape de nicho observado (D - el metrico de Schoener e I - el metrico de Warren)
ecospat.niche.overlap (z1=z1, z2=z2, cor=TRUE)

########################################################################
################### TEST DE EQUIVALENCIA DE NICHO ######################
########################################################################
a.dyn<-ecospat.niche.equivalency.test(z1=z1 , z2=z2, rep=1000) 

#GRAFICO PARA EL TEST DE EQUIVALENCIA DE NICHO
ecospat.plot.overlap.test(a.dyn,"D","Equivalency Momotus")




###############################################################################################
####### PAIRWISE COMPARISONS BETWEEN VS WITH DISTINCT "M" (ALLOPATRIC DISTRIBUTION) ###########
###############################################################################################

#######################################################################
################### MPS vs YP ################################
#######################################################################
##1. Leer los archivos .CSV de los puntos de ocurrencia de las especies.
setwd("WORKING DIRECTORY")
VS_1 <- read.csv("MPS_VS.csv", header = T, sep = ",")
VS_1$lat <- as.numeric(as.character(VS_1$lat))
VS_1$lon <- as.numeric(as.character(VS_1$lon))

VS_2 <- read.csv("YP_VS.csv", header = T, sep = ",")
VS_2$lat <- as.numeric(as.character(VS$lat))
VS_2$lon <- as.numeric(as.character(VS$lon))

##2. Leer los archivos de las condiciones climaticas de la M de cada especie.
## especie virtual 1
setwd("~/MPS/YP") 
variables_species1 <- list.files(".",pattern = "*.asc$",full.names = T) 
varclim_VS_1 <- stack(variables_species1)

## especie virtual 2
setwd("~/YP/MPS") 
variables_species2 <- list.files(".",pattern = "*.asc$",full.names = T)
varclim_VS_2 <- stack(variables_species2)

##3. Crear un data frame de las condiciones climaticas del area M de cada especie.
clim_punto_sp1 <- rasterToPoints(varclim_VS_1[[1]], fun=NULL, spatial=TRUE)
clim_punto_sp2 <- rasterToPoints(varclim_VS_2[[1]], fun=NULL, spatial=TRUE)

##4. Extraer los valores de todas las variables ambientales para cada punto de la M

## virtualspecies 1
clima_VS1 <- extract(varclim_VS_1, clim_punto_sp1)
clima_VS1 <- data.frame(coordinates(clim_punto_sp1),clima_VS1)
clima_VS1 <- subset(clima_VS1,!is.na(bio2) & !is.na(bio3) & !is.na(bio14) & !is.na(bio18) & !is.na(bio19))

## virtualspecies 2
clima_VS2 <- extract(varclim_VS_2, clim_punto_sp2)
clima_VS2 <- data.frame(coordinates(clim_punto_sp2),clima_VS2)
clima_VS2 <- subset(clima_VS2, !is.na(bio2) & !is.na(bio3) & !is.na(bio14) & !is.na(bio18) & !is.na(bio19))

#5. Seleccionar los datos de presencia de las especies pero considerando solo las coordenadas geograficas.
occ.sp1 <- VS_1[2:3]
occ.sp2 <- VS_2[2:3]

#6.Integrar la informacion de datos de ocurrencia y los datos del background para cada especie. Adem?s limpiar los datos a la resoluci?n espacial deseada (en este caso 1km)
occ_vs1<-na.exclude(ecospat.sample.envar(dfsp=occ.sp1,colspxy=1:2, 
                                         colspkept=1:2,dfvar=clima_VS1,
                                         colvarxy=1:2,colvar="all",resolution= 0.041665))

occ_vs2<-na.exclude(ecospat.sample.envar(dfsp=occ.sp2,colspxy=1:2,
                                         colspkept=1:2,dfvar=clima_VS2, 
                                         colvarxy=1:2,colvar="all",resolution= 0.041665))


#7. Agregar una nueva columna para diferenciar los datos de presencia (1) vs. el background (0) para luego unir los dos set de datos para cada especie.
# VS 1
occ_vs1 <- cbind(occ_vs1,species_occ = 1)
clima_VS1 <- cbind(clima_VS1,species_occ = 0)

names(clima_VS1)[1] = "lon" 
names(clima_VS1)[2] = "lat" 

data_vs_1 <- rbind(occ_vs1, clima_VS1)

# VS 2
occ_vs2 <- cbind(occ_vs2,species_occ = 1)
clima_VS2 <- cbind(clima_VS2,species_occ = 0)

names(clima_VS2)[1] = "lon" 
names(clima_VS2)[2] = "lat" 

data_vs_2 <- rbind(occ_vs2, clima_VS2)


##8. Calcula los valores de PCA para las variables que se estan comparando al analizar los nichos de las dos especies
pca.env <- dudi.pca(rbind(data_vs_1,data_vs_2)[,3:7],center = T, scale = T, scannf=F,nf=2)
ecospat.plot.contrib(contrib=pca.env$co, eigen=pca.env$eig)

##Rescata los valores "scores" de los ejes para los componentes obtenidos
scores.globclim <- pca.env$li

scores.vs1 <- suprow(pca.env,data_vs_1[which(data_vs_1[,8]==1),3:7])$li
scores.vs2 <- suprow(pca.env,data_vs_2[which(data_vs_2[,8]==1),3:7])$li
scores.clim.vs1 <- suprow(pca.env,data_vs_1[which(data_vs_1[,8]==0),3:7])$li
scores.clim.vs2 <- suprow(pca.env,data_vs_2[which(data_vs_2[,8]==0),3:7])$li

#9. Calcular las gr?ficas de densidad de ocurrencia 
#Species1
grid.clim.vs1 <- ecospat.grid.clim.dyn(glob= scores.globclim,
                                       glob1= scores.clim.vs1,
                                       sp=scores.sp, R=200,
                                       th.sp=0)
#Species2
grid.clim.vs2 <- ecospat.grid.clim.dyn(glob=scores.globclim,
                                       glob1=scores.clim.vs2,
                                       sp=scores.vs, R=200,
                                       th.sp=0)

###Calculate Niche Overlap with ecospat.niche.overlap()Compute Schoener's D, index of niche overlap
D.overlap <- ecospat.niche.overlap (grid.clim.vs1, grid.clim.vs2, cor= T)
D.overlap

#10. Delimitat los nichos de las especies y cuantificar el grado de sobrelape.
niche.dyn <- ecospat.niche.dyn.index (grid.clim.vs1, grid.clim.vs2, intersection = NA)

########################################################################
################### TEST DE EQUIVALENCIA DE NICHO ######################
########################################################################
a.dyn<-ecospat.niche.equivalency.test(z1=grid.clim.vs1 , z2=grid.clim.vs2, rep=1000) #m?nimo son 100 en art?culos, est?ndar son 1000

#GRAFICO PARA EL TEST DE EQUIVALENCIA DE NICHO
ecospat.plot.overlap.test(a.dyn,"D","Equivalency YP-VS vs MPS-VS")





#######################################################################
################### MPS vs YP/MPS ################################
#######################################################################
##1. Leer los archivos .CSV de los puntos de ocurrencia de las especies.
setwd("WORKING DIRECTORY")
VS_1 <- read.csv("YP_MPS_VS.csv", header = T, sep = ",")
VS_1$lat <- as.numeric(as.character(VS_1$lat))
VS_1$lon <- as.numeric(as.character(VS_1$lon))

VS_2 <- read.csv("MPS_VS.csv", header = T, sep = ",")
VS_2$lat <- as.numeric(as.character(VS$lat))
VS_2$lon <- as.numeric(as.character(VS$lon))

##2. Leer los archivos de las condiciones climaticas de la M de cada especie.
## especie virtual 1
setwd("~/YP_MPS/MPS") 
variables_species1 <- list.files(".",pattern = "*.asc$",full.names = T) 
varclim_VS_1 <- stack(variables_species1)

## especie virtual 2
setwd("~/MPS/YP_MPS") 
variables_species2 <- list.files(".",pattern = "*.asc$",full.names = T)
varclim_VS_2 <- stack(variables_species2)

##3. Crear un data frame de las condiciones climaticas del area M de cada especie.
clim_punto_sp1 <- rasterToPoints(varclim_VS_1[[1]], fun=NULL, spatial=TRUE)
clim_punto_sp2 <- rasterToPoints(varclim_VS_2[[1]], fun=NULL, spatial=TRUE)

##4. Extraer los valores de todas las variables ambientales para cada punto de la M
## virtualspecies 1
clima_VS1 <- extract(varclim_VS_1, clim_punto_sp1)
clima_VS1 <- data.frame(coordinates(clim_punto_sp1),clima_VS1)
clima_VS1 <- subset(clima_species,!is.na(bio2) & !is.na(bio3) & !is.na(bio14) & !is.na(bio15) & !is.na(bio18) & !is.na(bio19))

## virtualspecies 2
clima_VS2 <- extract(varclim_VS_2, clim_punto_sp2)
clima_VS2 <- data.frame(coordinates(clim_punto_sp2),clima_VS2)
clima_VS2 <- subset(clima_VS, !is.na(bio2) & !is.na(bio3) & !is.na(bio14) & !is.na(bio15) & !is.na(bio18) & !is.na(bio19))

#5. Seleccionar los datos de presencia de las especies pero considerando solo las coordenadas geograficas.
occ.sp1 <- VS_1[2:3]
occ.sp2 <- VS_2[2:3]

#6.Integrar la informacion de datos de ocurrencia y los datos del background para cada especie. Adem?s limpiar los datos a la resoluci?n espacial deseada (en este caso 1km)
occ_vs1<-na.exclude(ecospat.sample.envar(dfsp=occ.sp1,colspxy=1:2, 
                                         colspkept=1:2,dfvar=clima_VS1,
                                         colvarxy=1:2,colvar="all",resolution= 0.041665))

occ_vs2<-na.exclude(ecospat.sample.envar(dfsp=occ.sp2,colspxy=1:2,
                                         colspkept=1:2,dfvar=clima_VS2, 
                                         colvarxy=1:2,colvar="all",resolution= 0.041665))


#7. Agregar una nueva columna para diferenciar los datos de presencia (1) vs. el background (0) para luego unir los dos set de datos para cada especie.
# VS 1
occ_vs1 <- cbind(occ_vs1,species_occ = 1)
clima_VS1 <- cbind(clima_VS1,species_occ = 0)

names(clima_VS1)[1] = "lon" 
names(clima_VS1)[2] = "lat" 

data_vs_1 <- rbind(occ_vs1, clima_VS1)

# VS 2
occ_vs2 <- cbind(occ_vs2,species_occ = 1)
clima_VS2 <- cbind(clima_VS2,species_occ = 0)

names(clima_VS2)[1] = "lon" 
names(clima_VS2)[2] = "lat" 

data_vs_2 <- rbind(occ_vs2, clima_VS2)


##8. Calcula los valores de PCA para las variables que se estan comparando al analizar los nichos de las dos especies
pca.env <- dudi.pca(rbind(data_vs_1,data_vs_2)[,3:8],center = T, scale = T, scannf=F,nf=2)
ecospat.plot.contrib(contrib=pca.env$co, eigen=pca.env$eig)

##Rescata los valores "scores" de los ejes para los componentes obtenidos
scores.globclim <- pca.env$li

scores.vs1 <- suprow(pca.env,data_vs_1[which(data_vs_1[,9]==1),3:8])$li
scores.vs2 <- suprow(pca.env,data_vs_2[which(data_vs_2[,9]==1),3:8])$li
scores.clim.vs1 <- suprow(pca.env,data_vs_1[which(data_vs_1[,9]==0),3:8])$li
scores.clim.vs2 <- suprow(pca.env,data_vs_2[which(data_vs_2[,9]==0),3:8])$li

#9. Calcular las gr?ficas de densidad de ocurrencia 
#Species1
grid.clim.vs1 <- ecospat.grid.clim.dyn(glob= scores.globclim,
                                       glob1= scores.clim.vs1,
                                       sp=scores.sp, R=200,
                                       th.sp=0)
#Species2
grid.clim.vs2 <- ecospat.grid.clim.dyn(glob=scores.globclim,
                                       glob1=scores.clim.vs2,
                                       sp=scores.vs, R=200,
                                       th.sp=0)

###Calculate Niche Overlap with ecospat.niche.overlap()Compute Schoener's D, index of niche overlap
D.overlap <- ecospat.niche.overlap (grid.clim.vs1, grid.clim.vs2, cor= T)
D.overlap

#10. Delimitat los nichos de las especies y cuantificar el grado de sobrelape.
niche.dyn <- ecospat.niche.dyn.index (grid.clim.vs1, grid.clim.vs2, intersection = NA)

########################################################################
################### TEST DE EQUIVALENCIA DE NICHO ######################
########################################################################
a.dyn<-ecospat.niche.equivalency.test(z1=grid.clim.vs1 , z2=grid.clim.vs2, rep=1000) #m?nimo son 100 en art?culos, est?ndar son 1000

#GRAFICO PARA EL TEST DE EQUIVALENCIA DE NICHO
ecospat.plot.overlap.test(a.dyn,"D","Equivalency YP/MPS-VS vs MPS-VS")




#######################################################################
################### YP vs YP/MPS ################################
#######################################################################
##1. Leer los archivos .CSV de los puntos de ocurrencia de las especies.
setwd("WORKING DIRECTORY")
VS_1 <- read.csv("YP_VS.csv", header = T, sep = ",")
VS_1$lat <- as.numeric(as.character(VS_1$lat))
VS_1$lon <- as.numeric(as.character(VS_1$lon))

VS_2 <- read.csv("YP_MPS_VS.csv", header = T, sep = ",")
VS_2$lat <- as.numeric(as.character(VS$lat))
VS_2$lon <- as.numeric(as.character(VS$lon))

##2. Leer los archivos de las condiciones climaticas de la M de cada especie.
## especie virtual 1
setwd("~/YP/YP_MPS") 
variables_species1 <- list.files(".",pattern = "*.asc$",full.names = T) 
varclim_VS_1 <- stack(variables_species1)

## especie virtual 2
setwd("~/YP_MPS/YS") 
variables_species2 <- list.files(".",pattern = "*.asc$",full.names = T)
varclim_VS_2 <- stack(variables_species2)

##3. Crear un data frame de las condiciones climaticas del area M de cada especie.
clim_punto_sp1 <- rasterToPoints(varclim_VS_1[[1]], fun=NULL, spatial=TRUE)
clim_punto_sp2 <- rasterToPoints(varclim_VS_2[[1]], fun=NULL, spatial=TRUE)

##4. Extraer los valores de todas las variables ambientales para cada punto de la M
## virtualspecies 1
clima_VS1 <- extract(varclim_VS_1, clim_punto_sp1)
clima_VS1 <- data.frame(coordinates(clim_punto_sp1),clima_VS1)
clima_VS1 <- subset(clima_species,!is.na(bio2) & !is.na(bio3) & !is.na(bio8) & !is.na(bio9) & !is.na(bio14) & !is.na(bio18) & !is.na(bio19))

## virtualspecies 2
clima_VS2 <- extract(varclim_VS_2, clim_punto_sp2)
clima_VS2 <- data.frame(coordinates(clim_punto_sp2),clima_VS2)
clima_VS2 <- subset(clima_VS, !is.na(bio2) & !is.na(bio3) & !is.na(bio8) & !is.na(bio9) & !is.na(bio14) & !is.na(bio18) & !is.na(bio19))

#5. Seleccionar los datos de presencia de las especies pero considerando solo las coordenadas geograficas.
occ.sp1 <- VS_1[2:3]
occ.sp2 <- VS_2[2:3]

#6.Integrar la informacion de datos de ocurrencia y los datos del background para cada especie. Adem?s limpiar los datos a la resoluci?n espacial deseada (en este caso 1km)
occ_vs1<-na.exclude(ecospat.sample.envar(dfsp=occ.sp1,colspxy=1:2, 
                                         colspkept=1:2,dfvar=clima_VS1,
                                         colvarxy=1:2,colvar="all",resolution= 0.041665))

occ_vs2<-na.exclude(ecospat.sample.envar(dfsp=occ.sp2,colspxy=1:2,
                                         colspkept=1:2,dfvar=clima_VS2, 
                                         colvarxy=1:2,colvar="all",resolution= 0.041665))


#7. Agregar una nueva columna para diferenciar los datos de presencia (1) vs. el background (0) para luego unir los dos set de datos para cada especie.
# VS 1
occ_vs1 <- cbind(occ_vs1,species_occ = 1)
clima_VS1 <- cbind(clima_VS1,species_occ = 0)

names(clima_VS1)[1] = "lon" 
names(clima_VS1)[2] = "lat" 

data_vs_1 <- rbind(occ_vs1, clima_VS1)

# VS 2
occ_vs2 <- cbind(occ_vs2,species_occ = 1)
clima_VS2 <- cbind(clima_VS2,species_occ = 0)

names(clima_VS2)[1] = "lon" 
names(clima_VS2)[2] = "lat" 

data_vs_2 <- rbind(occ_vs2, clima_VS2)


##8. Calcula los valores de PCA para las variables que se estan comparando al analizar los nichos de las dos especies
pca.env <- dudi.pca(rbind(data_vs_1,data_vs_2)[,3:9],center = T, scale = T, scannf=F,nf=2)
ecospat.plot.contrib(contrib=pca.env$co, eigen=pca.env$eig)

##Rescata los valores "scores" de los ejes para los componentes obtenidos
scores.globclim <- pca.env$li

scores.vs1 <- suprow(pca.env,data_vs_1[which(data_vs_1[,10]==1),3:9])$li
scores.vs2 <- suprow(pca.env,data_vs_2[which(data_vs_2[,10]==1),3:9])$li
scores.clim.vs1 <- suprow(pca.env,data_vs_1[which(data_vs_1[,10]==0),3:9])$li
scores.clim.vs2 <- suprow(pca.env,data_vs_2[which(data_vs_2[,10]==0),3:9])$li

#9. Calcular las gr?ficas de densidad de ocurrencia 
#Species1
grid.clim.vs1 <- ecospat.grid.clim.dyn(glob= scores.globclim,
                                       glob1= scores.clim.vs1,
                                       sp=scores.sp, R=200,
                                       th.sp=0)
#Species2
grid.clim.vs2 <- ecospat.grid.clim.dyn(glob=scores.globclim,
                                       glob1=scores.clim.vs2,
                                       sp=scores.vs, R=200,
                                       th.sp=0)

###Calculate Niche Overlap with ecospat.niche.overlap()Compute Schoener's D, index of niche overlap
D.overlap <- ecospat.niche.overlap (grid.clim.vs1, grid.clim.vs2, cor= T)
D.overlap

#10. Delimitat los nichos de las especies y cuantificar el grado de sobrelape.
niche.dyn <- ecospat.niche.dyn.index (grid.clim.vs1, grid.clim.vs2, intersection = NA)

########################################################################
################### TEST DE EQUIVALENCIA DE NICHO ######################
########################################################################
a.dyn<-ecospat.niche.equivalency.test(z1=grid.clim.vs1 , z2=grid.clim.vs2, rep=1000) #m?nimo son 100 en art?culos, est?ndar son 1000

#GRAFICO PARA EL TEST DE EQUIVALENCIA DE NICHO
ecospat.plot.overlap.test(a.dyn,"D","Equivalency YP-VS vs YP/MPS-VS")

#FIN

