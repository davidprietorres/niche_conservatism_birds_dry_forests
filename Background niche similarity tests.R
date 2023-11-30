
library(ecospat)  
library(raster)     
library(SDMTools) 
library(dismo)      
library(phyloclim)
library(ade4)
library(knitr)
library(spThin)
library(humboldt)


setwd("WORKING DIRECTORY")

##########################################################################################
###### PAIRWISE COMPARISONS WITH A SINGLE COMBINED "M" (SYMPATRIC DISTRIBUTION) ###########
##########################################################################################


################################
######### COLINUS ############
################################
sp <- read.csv("~/Colinus_nigrogularis_model.csv", header=TRUE)
ev <- read.csv("~/YP_VS.csv", header=TRUE)

# Variables climaticas que se usaran para generar los componentes principales
# Las capas climaticas debieron ser previamente recortadas solo para el o las areas geograficas relevante para la comparacion de los nichos ecologicos, en este caso Norteamerica.
setwd("~/M/colinus_yp/asc")
varclim <- stack(list.files(pattern = "*.asc$",full.names = T)) 
resol <- res(varclim)[1]

# Generar una tabla que para cada pixel (coordenada x, y) me diga el valor de todas las variables ambientales
climpunto <- rasterToPoints(varclim[[1]], fun=NULL, spatial=TRUE)

# Extraemos datos ambientales para cada capa en varclim de cada coordenada en climpunto
clim <- raster::extract(varclim, climpunto)

# Formateamos clim para que sea una tabla normal con dos primeras columnas x y y
clim <- data.frame(coordinates(climpunto),clim)
clim <- subset(clim, !is.na(bio3) & !is.na(bio8) & !is.na(bio9) & !is.na(bio11) & !is.na(bio18))

occ.sp1 <- sp[2:3]
occ.sp2 <- ev[2:3]

#Anadir variables climaticas a datos
occ_sp1<-na.exclude(ecospat.sample.envar(dfsp=occ.sp1,colspxy=1:2, 
                                         colspkept=1:2,dfvar=clim, colvarxy=1:2,colvar="all",resolution=0.0416666))

occ_sp2<-na.exclude(ecospat.sample.envar(dfsp=occ.sp2,colspxy=1:2, 
                                         colspkept=1:2,dfvar=clim, colvarxy=1:2,colvar="all",resolution=0.0416666))

iterations<-1000
R=1000

########################################################################
############################ PCA-AMBIENTE ##############################
########################################################################
data <-rbind(clim[,3:7],occ_sp1[,3:7],occ_sp2[,3:7]) # para dos especies
data <- subset(data, !is.na(bio3) & !is.na(bio8) & !is.na(bio9) & !is.na(bio11) & !is.na(bio18))

# vector de peso 0 para las ocurrencias y 1 para todos los sitios del area de estudio
w<-c(rep(1,nrow(clim)),rep(0,nrow(occ_sp1)),rep(0,nrow(occ_sp2))) 

pca.cal <-dudi.pca(data, row.w = w, center = T, scale = T, scannf = F, nf = 2)

# Filas en que estan los datos de clim y de cada una de las especies:
row.clim <- 1:nrow(clim)
row.sp1 <-  (1+nrow(clim)):(nrow(clim) + nrow(occ_sp1))
row.sp2 <-  (1 + nrow(clim) + nrow(occ_sp1)) : (nrow(clim) + nrow(occ_sp1) + nrow(occ_sp2))

# las coordenadas en cada uno de los ejes del PCA para los datos de toda el area de estudio y de cada una de las especies
scores.clim <- pca.cal$li[row.clim,] 
scores.sp1  <- pca.cal$li[row.sp1,]   
scores.sp2  <- pca.cal$li[row.sp2,]   

#Contribucion de cada variable a cada componente del PCA
ecospat.plot.contrib(contrib=pca.cal$co, eigen=pca.cal$eig)

####################
######HUMBOLDT######
####################

z_1<- humboldt.grid.espace(scores.clim,
                           scores.clim,
                           scores.sp1,
                           kern.smooth=1,R=200)

z_2<- humboldt.grid.espace(scores.clim,
                           scores.clim,
                           scores.sp2,
                           kern.smooth=1,R=200)

par(mfrow=c(1,2))
plot.sp<- humboldt.plot.niche(z_1, title = "", name.axis1 = "PC1", name.axis2 = "PC2", 
                              correct.env = F, color.ramp = 1)

plot.vs <- humboldt.plot.niche(z_2, title = "", name.axis1 = "PC1", name.axis2 = "PC2", 
                               correct.env = F, color.ramp = 1)

zz <- humboldt.g2e(env1= clim, 
                   env2= clim, 
                   sp1= sp, 
                   sp2= ev, 
                   reduce.env = 2,
                   reductype = "PCA",
                   env.trim =F,
                   rarefy.dist = 5,
                   rarefy.units = "km",
                   env.reso = 0.041665)

humboldt.plot.overlap(in.g2e = zz,
                      pdf.out = F,
                      pcx = 1,
                      pcy = 2,
                      swap=F)

########################################################################
########### SUPERFICIE DE DENSIDAD DE OCURRENCIAS ######################
########################################################################
z1<-ecospat.grid.clim.dyn (scores.clim, scores.clim, scores.sp1, R=200) 
z2<-ecospat.grid.clim.dyn (scores.clim, scores.clim, scores.sp2, R=200) 

ecospat.niche.overlap (z1=z1, z2=z2, cor=TRUE)

########################################################################
#################### Background niche similarity tests ######################
########################################################################
b.dyn_phor_phai<-ecospat.niche.similarity.test(z1=z1 , z2=z2, rep=1000)
b.dyn_phor_phai2<-ecospat.niche.similarity.test(z1=z2 , z2=z1, rep=1000)

par(mfrow=c(1,2))
ecospat.plot.overlap.test(b.dyn_phor_phai,"D","C. nigrogularis vs. YP-VS")
ecospat.plot.overlap.test(b.dyn_phor_phai2,"D", "YP-VS vs. C. nigrogularis")




######################################
######### CAMPYLORHYNHUS #############
#####################################
sp <- read.csv("~/Campylorhynchus_yucatanicus.csv", header=TRUE)
ev <- read.csv("~/YP_VS.csv", header=TRUE)

# Variables climaticas que se usaran para generar los componentes principales
# Las capas climaticas debieron ser previamente recortadas solo para el o las areas geograficas relevante para la comparacion de los nichos ecologicos, en este caso Norteamerica.
setwd("~/campylorhynchus_yp/asc")
varclim <- stack(list.files(pattern = "*.asc$",full.names = T)) 
resol <- res(varclim)[1]

# Generar una tabla que para cada pixel (coordenada x, y) me diga el valor de todas las variables ambientales
climpunto <- rasterToPoints(varclim[[1]], fun=NULL, spatial=TRUE)

# Extraemos datos ambientales para cada capa en varclim de cada coordenada en climpunto
clim <- raster::extract(varclim, climpunto)

# Formateamos clim para que sea una tabla normal con dos primeras columnas x y y
clim <- data.frame(coordinates(climpunto),clim)
clim <- subset(clim, !is.na(bio8) & !is.na(bio14) & !is.na(bio18) & !is.na(bio19))

occ.sp1 <- sp[2:3]
occ.sp2 <- ev[2:3]

#Anadir variables climaticas a datos
occ_sp1<-na.exclude(ecospat.sample.envar(dfsp=occ.sp1,colspxy=1:2, 
                                         colspkept=1:2,dfvar=clim, colvarxy=1:2,colvar="all",resolution=0.041666))

occ_sp2<-na.exclude(ecospat.sample.envar(dfsp=occ.sp2,colspxy=1:2, 
                                         colspkept=1:2,dfvar=clim, colvarxy=1:2,colvar="all",resolution=0.041666))

iterations<-1000
R=1000

########################################################################
############################ PCA-AMBIENTE ##############################
########################################################################
data<-rbind(clim[,3:6],occ_sp1[,3:6],occ_sp2[,3:6]) # para dos especies

data <- subset(data, !is.na(bio8) & !is.na(bio14) & !is.na(bio18) & !is.na(bio19))

# vector de peso 0 para las ocurrencias y 1 para todos los sitios del area de estudio
w<-c(rep(1,nrow(clim)),rep(0,nrow(occ_sp1)),rep(0,nrow(occ_sp2))) # para tres especies

pca.cal <-dudi.pca(data, row.w = w, center = T, scale = T, scannf = F, nf = 2)

# Filas en que estan los datos de clim y de cada una de las especies:
row.clim <- 1:nrow(clim)
row.sp1 <-  (1+nrow(clim)):(nrow(clim) + nrow(occ_sp1))
row.sp2 <-  (1 + nrow(clim) + nrow(occ_sp1)) : (nrow(clim) + nrow(occ_sp1) + nrow(occ_sp2))

# las coordenadas en cada uno de los ejes del PCA para los datos de toda el area de estudio y de cada una de las especies
scores.clim <- pca.cal$li[row.clim,] 
scores.sp1  <- pca.cal$li[row.sp1,]   
scores.sp2  <- pca.cal$li[row.sp2,]  

#Contribucion de cada variable a cada componente del PCA
par(mfrow=c(1,1))
ecospat.plot.contrib(contrib=pca.cal$co, eigen=pca.cal$eig)

####################
######HUMBOLDT######
####################

z_1<- humboldt.grid.espace(scores.clim,
                           scores.clim,
                           scores.sp1,
                           kern.smooth=1,R=200)

z_2<- humboldt.grid.espace(scores.clim,
                           scores.clim,
                           scores.sp2,
                           kern.smooth=1,R=200)

par(mfrow=c(1,2))
plot.sp<- humboldt.plot.niche(z_1, title = "", name.axis1 = "PC1", name.axis2 = "PC2", 
                              correct.env = F, color.ramp = 1)

plot.vs <- humboldt.plot.niche(z_2, title = "", name.axis1 = "PC1", name.axis2 = "PC2", 
                               correct.env = F, color.ramp = 1)

zz <- humboldt.g2e(env1= clim, 
                   env2= clim, 
                   sp1= sp, 
                   sp2= ev, 
                   reduce.env = 2,
                   reductype = "PCA",
                   env.trim =F,
                   rarefy.dist = 5,
                   rarefy.units = "km",
                   env.reso = 0.041665)

humboldt.plot.overlap(in.g2e = zz,
                      pcx = 1,
                      pcy = 2,
                      swap=F)

########################################################################
########### SUPERFICIE DE DENSIDAD DE OCURRENCIAS ######################
########################################################################
z1<-ecospat.grid.clim.dyn (scores.clim, scores.clim, scores.sp1, R=200) 
z2<-ecospat.grid.clim.dyn (scores.clim, scores.clim, scores.sp2, R=200) 

ecospat.niche.overlap (z1=z1, z2=z2, cor=TRUE)

########################################################################
#################### Background niche similarity tests ######################
########################################################################
b.dyn_phor_phai<-ecospat.niche.similarity.test(z1=z1 , z2=z2, rep=1000)
b.dyn_phor_phai2<-ecospat.niche.similarity.test(z1=z2 , z2=z1, rep=1000)

par(mfrow=c(1,2))
ecospat.plot.overlap.test(b.dyn_phor_phai,"D","C. yucatanicus vs. YP-VS")
ecospat.plot.overlap.test(b.dyn_phor_phai2,"D", "YP-VS vs. C. yucatanicus")




################################
######### GEOCOCCYX ############
################################
sp <- read.csv("~/Geococcyx_velox_model.csv", header=TRUE)
ev <- read.csv("~/YP_MPS_VS.csv", header=TRUE)

# Variables climaticas que se usaran para generar los componentes principales
# Las capas climaticas debieron ser previamente recortadas solo para el o las areas geograficas relevante para la comparacion de los nichos ecologicos, en este caso Norteamerica.
setwd("~/geococcyx_ypmps/asc")
varclim <- stack(list.files(pattern = "*.asc$",full.names = T))s
resol <- res(varclim)[1]

# Generar una tabla que para cada pixel (coordenada x, y) me diga el valor de todas las variables ambientales
climpunto <- rasterToPoints(varclim[[1]], fun=NULL, spatial=TRUE)

# Extraemos datos ambientales para cada capa en varclim de cada coordenada en climpunto
clim <- raster::extract(varclim, climpunto)

# Formateamos clim para que sea una tabla normal con dos primeras columnas x y y
clim <- data.frame(coordinates(climpunto),clim)
clim <- subset(clim, !is.na(bio2) & !is.na(bio3) & !is.na(bio9) & !is.na(bio14) & !is.na(bio15) & !is.na(bio18) & !is.na(bio19))

occ.sp1 <- sp[2:3]
occ.sp2 <- ev[2:3]

#Anadir variables climaticas a datos
occ_sp1<-na.exclude(ecospat.sample.envar(dfsp=occ.sp1,colspxy=1:2, 
                                         colspkept=1:2,dfvar=clim, colvarxy=1:2,colvar="all",resolution=0.0416666))

occ_sp2<-na.exclude(ecospat.sample.envar(dfsp=occ.sp2,colspxy=1:2, 
                                         colspkept=1:2,dfvar=clim, colvarxy=1:2,colvar="all",resolution=0.0416666))

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

# las coordenadas en cada uno de los ejes del PCA para los datos de toda el area de estudio y de cada una de las especies
scores.clim <- pca.cal$li[row.clim,] 
scores.sp1  <- pca.cal$li[row.sp1,]   
scores.sp2  <- pca.cal$li[row.sp2,]   

#Contribucion de cada variable a cada componente del PCA
ecospat.plot.contrib(contrib=pca.cal$co, eigen=pca.cal$eig)

####################
######HUMBOLDT######
####################

z_1<- humboldt.grid.espace(scores.clim,
                           scores.clim,
                           scores.sp1,
                           kern.smooth=1,R=200)

z_2<- humboldt.grid.espace(scores.clim,
                           scores.clim,
                           scores.sp2,
                           kern.smooth=1,R=200)

par(mfrow=c(1,2))
plot.sp<- humboldt.plot.niche(z_1, title = "", name.axis1 = "PC1", name.axis2 = "PC2", 
                              correct.env = F, color.ramp = 1)

plot.vs <- humboldt.plot.niche(z_2, title = "", name.axis1 = "PC1", name.axis2 = "PC2", 
                               correct.env = F, color.ramp = 1)

zz <- humboldt.g2e(env1= clim, 
                   env2= clim, 
                   sp1= sp, 
                   sp2= ev, 
                   reduce.env = 2,
                   reductype = "PCA",
                   env.trim =F,
                   rarefy.dist = 5,
                   rarefy.units = "km",
                   env.reso = 0.041665)

humboldt.plot.overlap(in.g2e = zz,
                      pdf.out = T,
                      pcx = 1,
                      pcy = 2,
                      swap=F)


########################################################################
########### SUPERFICIE DE DENSIDAD DE OCURRENCIAS ######################
########################################################################
z1<-ecospat.grid.clim.dyn (scores.clim, scores.clim, scores.sp1, R=200) 
z2<-ecospat.grid.clim.dyn (scores.clim, scores.clim, scores.sp2, R=200) 

ecospat.niche.overlap (z1=z1, z2=z2, cor=TRUE)

########################################################################
#################### Background niche similarity tests ######################
########################################################################
b.dyn_phor_phai<-ecospat.niche.similarity.test(z1=z1 , z2=z2, rep=1000)
b.dyn_phor_phai2<-ecospat.niche.similarity.test(z1=z2 , z2=z1, rep=1000)

par(mfrow=c(1,2))
ecospat.plot.overlap.test(b.dyn_phor_phai,"D","G. velox vs. YP/MPS-VS")
ecospat.plot.overlap.test(b.dyn_phor_phai2,"D", "YP/MPS-VS vs. G. velox")




#############################
######### AMAZILIA ############
#############################
sp <- read.csv("~/Amazilia_rutila_model.csv", header=TRUE)
ev <- read.csv("~/YP_MPS_VS.csv", header=TRUE)

# Variables climaticas que se usaran para generar los componentes principales
# Las capas climaticas debieron ser previamente recortadas solo para el o las areas geograficas relevante para la comparacion de los nichos ecologicos, en este caso Norteamerica.
setwd("~/amazilia_ypmps/asc")
varclim <- stack(list.files(pattern = "*.asc$",full.names = T)) 
resol <- res(varclim)[1]

# Generar una tabla que para cada pixel (coordenada x, y) me diga el valor de todas las variables ambientales
climpunto <- rasterToPoints(varclim[[1]], fun=NULL, spatial=TRUE)

# Extraemos datos ambientales para cada capa en varclim de cada coordenada en climpunto
clim <- raster::extract(varclim, climpunto)

# Formateamos clim para que sea una tabla normal con dos primeras columnas x y y
clim <- data.frame(coordinates(climpunto),clim)
clim <- subset(clim, !is.na(bio2) & !is.na(bio3) & !is.na(bio9) & !is.na(bio13) & !is.na(bio14) & !is.na(bio15) & !is.na(bio18) & !is.na(bio19))

occ.sp1 <- sp[2:3]
occ.sp2 <- ev[2:3]

#Anadir variables climaticas a datos
occ_sp1<-na.exclude(ecospat.sample.envar(dfsp=occ.sp1,colspxy=1:2, 
                                         colspkept=1:2,dfvar=clim, colvarxy=1:2,colvar="all",resolution=0.0416666))

occ_sp2<-na.exclude(ecospat.sample.envar(dfsp=occ.sp2,colspxy=1:2, 
                                         colspkept=1:2,dfvar=clim, colvarxy=1:2,colvar="all",resolution=0.0416666))

iterations<-1000
R=1000

########################################################################
############################ PCA-AMBIENTE ##############################
########################################################################
data <-rbind(clim[,3:10],occ_sp1[,3:10],occ_sp2[,3:10]) 
data <- subset(data, !is.na(bio2) & !is.na(bio3) & !is.na(bio9) & !is.na(bio13) & !is.na(bio14) & !is.na(bio15) & !is.na(bio18) & !is.na(bio19))

# vector de peso 0 para las ocurrencias y 1 para todos los sitios del area de estudio
w<-c(rep(1,nrow(clim)),rep(0,nrow(occ_sp1)),rep(0,nrow(occ_sp2))) # para tres especies

pca.cal <-dudi.pca(data, row.w = w, center = T, scale = T, scannf = F, nf = 2)

# Filas en que estan los datos de clim y de cada una de las especies:
row.clim <- 1:nrow(clim)
row.sp1 <-  (1+nrow(clim)):(nrow(clim) + nrow(occ_sp1))
row.sp2 <-  (1 + nrow(clim) + nrow(occ_sp1)) : (nrow(clim) + nrow(occ_sp1) + nrow(occ_sp2))

# las coordenadas en cada uno de los ejes del PCA para los datos de toda el area de estudio y de cada una de las especies
scores.clim <- pca.cal$li[row.clim,] 
scores.sp1  <- pca.cal$li[row.sp1,]   
scores.sp2  <- pca.cal$li[row.sp2,]   

#Contribucion de cada variable a cada componente del PCA
ecospat.plot.contrib(contrib=pca.cal$co, eigen=pca.cal$eig)

####################
######HUMBOLDT######
####################

z_1<- humboldt.grid.espace(scores.clim,
                           scores.clim,
                           scores.sp1,
                           kern.smooth=1,R=200)

z_2<- humboldt.grid.espace(scores.clim,
                           scores.clim,
                           scores.sp2,
                           kern.smooth=1,R=200)

par(mfrow=c(1,2))
plot.sp<- humboldt.plot.niche(z_1, title = "", name.axis1 = "PC1", name.axis2 = "PC2", 
                              correct.env = F, color.ramp = 1)

plot.vs <- humboldt.plot.niche(z_2, title = "", name.axis1 = "PC1", name.axis2 = "PC2", 
                               correct.env = F, color.ramp = 1)

zz <- humboldt.g2e(env1= clim, 
                   env2= clim, 
                   sp1= sp, 
                   sp2= ev, 
                   reduce.env = 2,
                   reductype = "PCA",
                   env.trim =F,
                   rarefy.dist = 5,
                   rarefy.units = "km",
                   env.reso = 0.041665)

humboldt.plot.overlap(in.g2e = zz,
                      pcx = 1,
                      pcy = 2,
                      swap=F)

########################################################################
########### SUPERFICIE DE DENSIDAD DE OCURRENCIAS ######################
########################################################################
z1<-ecospat.grid.clim.dyn (scores.clim, scores.clim, scores.sp1, R=200) #spp sp
z2<-ecospat.grid.clim.dyn (scores.clim, scores.clim, scores.sp2, R=200) #spp ev

ecospat.niche.overlap (z1=z1, z2=z2, cor=TRUE)

########################################################################
#################### Background niche similarity tests ######################
########################################################################
b.dyn_phor_phai<-ecospat.niche.similarity.test(z1=z1 , z2=z2, rep=1000)
b.dyn_phor_phai2<-ecospat.niche.similarity.test(z1=z2 , z2=z1, rep=1000)

par(mfrow=c(1,2))
ecospat.plot.overlap.test(b.dyn_phor_phai,"D","A. rutila vs. YP/MPS-VS")
ecospat.plot.overlap.test(b.dyn_phor_phai2,"D","YP/MPS-VS vs. A. rutila")




#############################
######### FORPUS ############
#############################
sp <- read.csv("~/Forpus_cyanopygius_model.csv", header=TRUE)
ev <- read.csv("~/MPS_VS.csv", header=TRUE)

# Variables climaticas que se usaran para generar los componentes principales
# Las capas climaticas debieron ser previamente recortadas solo para el o las areas geograficas relevante para la comparacion de los nichos ecologicos, en este caso Norteamerica.
setwd("~/forpus_mps/asc")
varclim <- stack(list.files(pattern = "*.asc$",full.names = T)) 
resol <- res(varclim)[1]

# Generar una tabla que para cada pixel (coordenada x, y) me diga el valor de todas las variables ambientales
climpunto <- rasterToPoints(varclim[[1]], fun=NULL, spatial=TRUE)

# Extraemos datos ambientales para cada capa en varclim de cada coordenada en climpunto
clim <- raster::extract(varclim, climpunto)

# Formateamos clim para que sea una tabla normal con dos primeras columnas x y y
clim <- data.frame(coordinates(climpunto),clim)
clim <- subset(clim, !is.na(bio2) & !is.na(bio3) & !is.na(bio14) & !is.na(bio15) & !is.na(bio18) & !is.na(bio19))

occ.sp1 <- sp[2:3]
occ.sp2 <- ev[2:3]

#Anadir variables climaticas a datos
occ_sp1<-na.exclude(ecospat.sample.envar(dfsp=occ.sp1,colspxy=1:2, 
                                         colspkept=1:2,dfvar=clim, colvarxy=1:2,colvar="all",resolution=0.0416666))

occ_sp2<-na.exclude(ecospat.sample.envar(dfsp=occ.sp2,colspxy=1:2, 
                                         colspkept=1:2,dfvar=clim, colvarxy=1:2,colvar="all",resolution=0.0416666))

iterations<-1000
R=1000

########################################################################
############################ PCA-AMBIENTE ##############################
########################################################################
data <-rbind(clim[,3:8],occ_sp1[,3:8],occ_sp2[,3:8]) 

data <- subset(data, !is.na(bio2) & !is.na(bio3) & !is.na(bio14) & !is.na(bio15) & !is.na(bio18) & !is.na(bio19))

# vector de peso 0 para las ocurrencias y 1 para todos los sitios del area de estudio
w<-c(rep(1,nrow(clim)),rep(0,nrow(occ_sp1)),rep(0,nrow(occ_sp2))) # para tres especies

pca.cal <-dudi.pca(data, row.w = w, center = T, scale = T, scannf = F, nf = 2)

# Filas en que estan los datos de clim y de cada una de las especies:
row.clim <- 1:nrow(clim)
row.sp1 <-  (1+nrow(clim)):(nrow(clim) + nrow(occ_sp1))
row.sp2 <-  (1 + nrow(clim) + nrow(occ_sp1)) : (nrow(clim) + nrow(occ_sp1) + nrow(occ_sp2))

# las coordenadas en cada uno de los ejes del PCA para los datos de toda el area de estudio y de cada una de las especies
scores.clim <- pca.cal$li[row.clim,] 
scores.sp1  <- pca.cal$li[row.sp1,]   
scores.sp2  <- pca.cal$li[row.sp2,]   

#Contribucion de cada variable a cada componente del PCA
par(mfrow=c(1,1))
ecospat.plot.contrib(contrib=pca.cal$co, eigen=pca.cal$eig)


####################
######HUMBOLDT######
####################

z_1<- humboldt.grid.espace(scores.clim,
                           scores.clim,
                           scores.sp1,
                           kern.smooth=1,R=200)

z_2<- humboldt.grid.espace(scores.clim,
                           scores.clim,
                           scores.sp2,
                           kern.smooth=1,R=200)

par(mfrow=c(1,2))
plot.sp<- humboldt.plot.niche(z_1, title = "", name.axis1 = "PC1", name.axis2 = "PC2", 
                              correct.env = F, color.ramp = 1)

plot.vs <- humboldt.plot.niche(z_2, title = "", name.axis1 = "PC1", name.axis2 = "PC2", 
                               correct.env = F, color.ramp = 1)

zz <- humboldt.g2e(env1= clim, 
                   env2= clim, 
                   sp1= sp, 
                   sp2= ev, 
                   reduce.env = 2,
                   reductype = "PCA",
                   env.trim =F,
                   rarefy.dist = 5,
                   rarefy.units = "km",
                   env.reso = 0.041665)

humboldt.plot.overlap(in.g2e = zz,
                      pcx = 1,
                      pcy = 2,
                      swap=F)

########################################################################
########### SUPERFICIE DE DENSIDAD DE OCURRENCIAS ######################
########################################################################
z1<-ecospat.grid.clim.dyn (scores.clim, scores.clim, scores.sp1, R=200) 
z2<-ecospat.grid.clim.dyn (scores.clim, scores.clim, scores.sp2, R=200) 

ecospat.niche.overlap (z1=z1, z2=z2, cor=TRUE)

########################################################################
#################### Background niche similarity tests ######################
########################################################################
b.dyn_phor_phai<-ecospat.niche.similarity.test(z1=z1 , z2=z2, rep=1000)
b.dyn_phor_phai2<-ecospat.niche.similarity.test(z1=z2 , z2=z1, rep=1000)

par(mfrow=c(1,2))
ecospat.plot.overlap.test(b.dyn_phor_phai,"D","F. cyanopgygius vs. MPS-VS")
ecospat.plot.overlap.test(b.dyn_phor_phai2,"D","MPS-VS vs. F.cyanopygius")




#############################
######### MOMOTUS ############
#############################
sp <- read.csv("~/Momotus_mexicanus_model.csv", header=TRUE)
ev <- read.csv("~/MPS_VS.csv", header=TRUE)

# Variables climaticas que se usaran para generar los componentes principales
# Las capas climaticas debieron ser previamente recortadas solo para el o las areas geograficas relevante para la comparacion de los nichos ecologicos, en este caso Norteamerica.
setwd("~/momotus_mps/asc")
varclim <- stack(list.files(pattern = "*.asc$",full.names = T)) 
resol <- res(varclim)[1]

# Generar una tabla que para cada pixel (coordenada x, y) me diga el valor de todas las variables ambientales
climpunto <- rasterToPoints(varclim[[1]], fun=NULL, spatial=TRUE)

# Extraemos datos ambientales para cada capa en varclim de cada coordenada en climpunto
clim <- raster::extract(varclim, climpunto)

# Formateamos clim para que sea una tabla normal con dos primeras columnas x y y
clim <- data.frame(coordinates(climpunto),clim)
clim <- subset(clim, !is.na(bio2) & !is.na(bio3) & !is.na(bio4) & !is.na(bio14) & !is.na(bio15) & !is.na(bio18) & !is.na(bio19))

occ.sp1 <- sp[2:3]
occ.sp2 <- ev[2:3]

#Anadir variables climaticas a datos
occ_sp1<-na.exclude(ecospat.sample.envar(dfsp=occ.sp1,colspxy=1:2, 
                                         colspkept=1:2,dfvar=clim, colvarxy=1:2,colvar="all",resolution=0.0416666))

occ_sp2<-na.exclude(ecospat.sample.envar(dfsp=occ.sp2,colspxy=1:2, 
                                         colspkept=1:2,dfvar=clim, colvarxy=1:2,colvar="all",resolution=0.0416666))

iterations<-1000
R=1000

########################################################################
############################ PCA-AMBIENTE ##############################
########################################################################
data <-rbind(clim[,3:9],occ_sp1[,3:9],occ_sp2[,3:9]) # para dos especies

data <- subset(data, !is.na(bio2) & !is.na(bio3) & !is.na(bio4) & !is.na(bio14) & !is.na(bio15) & !is.na(bio18) & !is.na(bio19))

# vector de peso 0 para las ocurrencias y 1 para todos los sitios del area de estudio
w<-c(rep(1,nrow(clim)),rep(0,nrow(occ_sp1)),rep(0,nrow(occ_sp2)))

pca.cal <-dudi.pca(data, row.w = w, center = T, scale = T, scannf = F, nf = 2)

# Filas en que estan los datos de clim y de cada una de las especies:
row.clim <- 1:nrow(clim)
row.sp1 <-  (1+nrow(clim)):(nrow(clim) + nrow(occ_sp1))
row.sp2 <-  (1 + nrow(clim) + nrow(occ_sp1)) : (nrow(clim) + nrow(occ_sp1) + nrow(occ_sp2))

# las coordenadas en cada uno de los ejes del PCA para los datos de toda el area de estudio y de cada una de las especies
scores.clim <- pca.cal$li[row.clim,] 
scores.sp1  <- pca.cal$li[row.sp1,]   
scores.sp2  <- pca.cal$li[row.sp2,]   

#Contribucion de cada variable a cada componente del PCA
ecospat.plot.contrib(contrib=pca.cal$co, eigen=pca.cal$eig)

####################
######HUMBOLDT######
####################

z_1<- humboldt.grid.espace(scores.clim,
                           scores.clim,
                           scores.sp1,
                           kern.smooth=1,R=200)

z_2<- humboldt.grid.espace(scores.clim,
                           scores.clim,
                           scores.sp2,
                           kern.smooth=1,R=200)

par(mfrow=c(1,2))
plot.sp<- humboldt.plot.niche(z_1, title = "", name.axis1 = "PC1", name.axis2 = "PC2", 
                              correct.env = F, color.ramp = 1)

plot.vs <- humboldt.plot.niche(z_2, title = "", name.axis1 = "PC1", name.axis2 = "PC2", 
                               correct.env = F, color.ramp = 1)

zz <- humboldt.g2e(env1= clim, 
                   env2= clim, 
                   sp1= sp, 
                   sp2= ev, 
                   reduce.env = 2,
                   reductype = "PCA",
                   env.trim =F,
                   rarefy.dist = 5,
                   rarefy.units = "km",
                   env.reso = 0.041665)

humboldt.plot.overlap(in.g2e = zz,
                      pcx = 1,
                      pcy = 2,
                      swap=F)

########################################################################
########### SUPERFICIE DE DENSIDAD DE OCURRENCIAS ######################
########################################################################
z1<-ecospat.grid.clim.dyn (scores.clim, scores.clim, scores.sp1, R=200) 
z2<-ecospat.grid.clim.dyn (scores.clim, scores.clim, scores.sp2, R=200) 

ecospat.niche.overlap (z1=z1, z2=z2, cor=TRUE)

########################################################################
#################### Background niche similarity tests ######################
########################################################################
b.dyn_phor_phai<-ecospat.niche.similarity.test(z1=z1 , z2=z2, rep=1000)
b.dyn_phor_phai2<-ecospat.niche.similarity.test(z1=z2 , z2=z1, rep=1000)

par(mfrow=c(1,2))
ecospat.plot.overlap.test(b.dyn_phor_phai,"D","M. mexicanus vs. MPS-VS")
ecospat.plot.overlap.test(b.dyn_phor_phai2,"D","MPS-VS vs. M. mexicanus")




###############################################################################################
####### PAIRWISE COMPARISONS BETWEEN VS WITH DISTINCT "M" (ALLOPATRIC DISTRIBUTION) ###########
###############################################################################################

#######################################################################
################### YP vs MPS ################################
#######################################################################
##1. Leer los archivos .CSV de los puntos de ocurrencia de las especies.

setwd("WORKING DIRECTORY")
VS_1<- read.csv("MPS_VS.csv", header = T, sep = ",")
VS_1$lat <- as.numeric(as.character(VS_1$lat))
VS_1$lon <- as.numeric(as.character(VS_1$lon))

VS_2 <- read.csv("YP_VS.csv", header = T, sep = ",")
VS_2$lat <- as.numeric(as.character(VS_2$lat))
VS_2$lon <- as.numeric(as.character(VS_2$lon))

##2. Leer los archivos de las condiciones climaticas de la M de cada especie.
##especie virtual 1
setwd("~/MPS/YP") 
variables_VS1 <- list.files(".",pattern = "*.asc$",full.names = T)
varclim_vs1 <- stack(variables_VS1)

##especie virtual 2
setwd("~/YP/MPS") 
variables_VS2 <- list.files(".",pattern = "*.asc$",full.names = T) 
varclim_vs2 <- stack(variables_VS2)

##3. Crear un data frame de las condiciones climaticas del area M de cada especie.
clim_punto_sp1 <- rasterToPoints(varclim_vs1[[1]], fun=NULL, spatial=TRUE)
clim_punto_sp2 <- rasterToPoints(varclim_vs2[[1]], fun=NULL, spatial=TRUE)

##4. Extraer los valores de todas las variables ambientales para cada punto de la M
##species 1
clima_VS_1 <- extract(varclim_vs1, clim_punto_sp1)
clima_VS_1 <- data.frame(coordinates(clim_punto_sp1),clima_VS_1)
clima_VS_1 <- subset(clima_VS_1,!is.na(bio2) & !is.na(bio3) & !is.na(bio14) & !is.na(bio18) & !is.na(bio19))

##species 2
clima_VS_2 <- extract(varclim_vs2, clim_punto_sp2) 
clima_VS_2 <- data.frame(coordinates(clim_punto_sp2),clima_VS)
clima_VS_2 <- subset(clima_VS_2, !is.na(bio2) & !is.na(bio3) & !is.na(bio14) & !is.na(bio18) & !is.na(bio19))

#5. Seleccionar los datos de presencia de las especies pero considerando solo las coordenadas geogr?ficas.
occ.sp1 <- VS_1[2:3]
occ.sp2 <- VS_2[2:3]

#6.Integrar la informaci?n de datos de ocurrencia y los datos del background para cada especie. Adem?s limpiar los datos a la resoluci?n espacial deseada (en este caso 1km)
occ_vs1<-na.exclude(ecospat.sample.envar(dfsp=occ.sp1,colspxy=1:2, 
                                        colspkept=1:2,dfvar=clima_VS_1,
                                        colvarxy=1:2,colvar="all",resolution= 0.041665))

occ_vs2<-na.exclude(ecospat.sample.envar(dfsp=occ.sp2,colspxy=1:2,
                                        colspkept=1:2,dfvar=clima_VS_2, 
                                        colvarxy=1:2,colvar="all",resolution= 0.041665))


#7. Agregar una nueva columna para diferenciar los datos de presencia (1) vs. el background (0) para luego unir los dos set de datos para cada especie.

occ_vs1 <- cbind(occ_vs1,species_occ = 1)
clima_VS_1 <- cbind(clima_VS_1,species_occ = 0)

names(clima_VS_1)[1] = "lon" 
names(clima_VS_1)[2] = "lat" 
data_VS_1 <- rbind(occ_vs1, clima_VS_1)

#VS
occ_vs2 <- cbind(occ_vs2,species_occ = 1)
clima_VS_2 <- cbind(clima_VS_2,species_occ = 0)

names(clima_VS_2)[1] = "lon" 
names(clima_VS_2)[2] = "lat" 
data_VS_2 <- rbind(occ_vs2, clima_VS_2)


##8. Calcula los valores de PCA para las variables que se est?n comparando al analizar los nichos de las dos especies
pca.env <- dudi.pca(rbind(data_VS_1,data_VS_2)[,3:7],center = T, scale = T, scannf=F,nf=2)
ecospat.plot.contrib(contrib=pca.env$co, eigen=pca.env$eig)

##Rescata los valores "scores" de los ejes para los componentes obtenidos
scores.globclim <- pca.env$li# PCA scores para toda el ?rea de estudio.

scores.vs1 <- suprow(pca.env,data_VS_1[which(data_VS_1[,8]==1),3:7])$li
scores.vs2 <- suprow(pca.env,data_VS_1[which(data_VS_2[,8]==1),3:7])$li
scores.clim.vs1 <- suprow(pca.env,data_VS_1[which(data_VS_1[,8]==0),3:7])$li
scores.clim.vs2 <- suprow(pca.env,data_VS_2[which(data_VS_2[,8]==0),3:7])$li

####################
######HUMBOLDT######
####################

z_1<- humboldt.grid.espace(scores.globclim,
                           scores.clim.vs1,
                           scores.vs1,
                           kern.smooth=1,R=200)

z_2<- humboldt.grid.espace(scores.globclim,
                           scores.clim.vs2,
                           scores.vs2,
                           kern.smooth=1,R=200)

par(mfrow=c(1,2))
plot.vs1 <- humboldt.plot.niche(z_1, title = "MPS-VS", name.axis1 = "PC1", name.axis2 = "PC2", 
                              correct.env = F, color.ramp = 1)

plot.vs2 <- humboldt.plot.niche(z_2, title = "YP-VS", name.axis1 = "PC1", name.axis2 = "PC2", 
                               correct.env = F, color.ramp = 1)

#9. Calcular las gr?ficas de densidad de ocurrencia 
#Species1
grid.clim.vs1 <- ecospat.grid.clim.dyn(glob= scores.globclim,
                                      glob1= scores.clim.vs1,
                                      sp=scores.vs1, R=200,
                                      th.sp=0)
#Species2
grid.clim.vs2 <- ecospat.grid.clim.dyn(glob=scores.globclim,
                                      glob1=scores.clim.vs2,
                                      sp=scores.vs2, R=200,
                                      th.sp=0)

###Calculate Niche Overlap with ecospat.niche.overlap()Compute Schoener's D, index of niche overlap
D.overlap <- ecospat.niche.overlap (grid.clim.vs1, grid.clim.vs2, cor= T)
D.overlap

#10. Delimitat los nichos de las especies y cuantificar el grado de sobrelape.
niche.dyn <- ecospat.niche.dyn.index (grid.clim.vs1, grid.clim.vs2, intersection = NA)

##ver las graficas de los nichos por separado
ecospat.plot.niche (grid.clim.vs1, title="MPS-VS", name.axis1="PC1", name.axis2="PC2", cor=F)
ecospat.plot.niche (grid.clim.vs2, title="YP-VS", name.axis1="PC1", name.axis2="PC2", cor=T)

########################################################################
################### Background niche similarity tests ######################
########################################################################
sim.test_sp1 <- ecospat.niche.similarity.test(grid.clim.vs1, grid.clim.vs2,
                                              rep = 1000,
                                              rand.type=1)

sim.test_sp2 <- ecospat.niche.similarity.test(grid.clim.vs2, grid.clim.vs1,
                                              rep = 1000,
                                              rand.type=1)
par(mfrow=c(1,2))
ecospat.plot.overlap.test(sim.test_sp1,"D","MPS-VS vs. YP-VS")
ecospat.plot.overlap.test(sim.test_sp2,"D", "YP-VS vs. MPS-VS")



#######################################################################
################### YP_MPS VS YP ################################
#######################################################################
##1. Leer los archivos .CSV de los puntos de ocurrencia de las especies.
setwd("~")
species<- read.csv("YP_MPS_VS.csv", header = T, sep = ",")
species$lat <- as.numeric(as.character(species$lat))
species$lon <- as.numeric(as.character(species$lon))

VS <- read.csv("YP_VS.csv", header = T, sep = ",")
VS$lat <- as.numeric(as.character(VS$lat))
VS$lon <- as.numeric(as.character(VS$lon))

##2. Leer los archivos de las condiciones clim?ticas de la M de cada especie.
##especie 
setwd("~/MPS_YP/YP") 
variables_VS1 <- list.files(".",pattern = "*.asc$",full.names = T)
varclim_vs1 <- stack(variables_VS1)

##especie virtual
setwd("~/YP/MPS_YP") 
variables_VS2 <- list.files(".",pattern = "*.asc$",full.names = T) 
varclim_vs2 <- stack(variables_VS2)

##3. Crear un data frame de las condiciones clim?ticas del ?rea M de cada especie.
clim_punto_sp1 <- rasterToPoints(varclim_vs1[[1]], fun=NULL, spatial=TRUE)
clim_punto_sp2 <- rasterToPoints(varclim_vs2[[1]], fun=NULL, spatial=TRUE)

##4. Extraer los valores de todas las variables ambientales para cada punto de la M
##species 1
clima_VS_1 <- extract(varclim_vs1, clim_punto_sp1)
clima_VS_1 <- data.frame(coordinates(clim_punto_sp1),clima_VS_1)
clima_VS_1 <- subset(clima_VS_1,!is.na(bio2) & !is.na(bio3) & !is.na(bio8) & !is.na(bio9) & !is.na(bio14) & !is.na(bio18) & !is.na(bio19))

##species 2
clima_VS_2 <- extract(varclim_vs2, clim_punto_sp2) 
clima_VS_2 <- data.frame(coordinates(clim_punto_sp2),clima_VS_2)
clima_VS_2 <- subset(clima_VS_2, !is.na(bio2) & !is.na(bio3) & !is.na(bio8) & !is.na(bio9) & !is.na(bio14) & !is.na(bio18) & !is.na(bio19))

#5. Seleccionar los datos de presencia de las especies pero considerando solo las coordenadas geogr?ficas.
occ.sp1 <- VS_1[2:3]
occ.sp2 <- VS_2[2:3]

#6.Integrar la informaci?n de datos de ocurrencia y los datos del background para cada especie. Adem?s limpiar los datos a la resoluci?n espacial deseada (en este caso 1km)
occ_vs1<-na.exclude(ecospat.sample.envar(dfsp=occ.sp1,colspxy=1:2, 
                                         colspkept=1:2,dfvar=clima_VS_1,
                                         colvarxy=1:2,colvar="all",resolution= 0.041665))

occ_vs2<-na.exclude(ecospat.sample.envar(dfsp=occ.sp2,colspxy=1:2,
                                         colspkept=1:2,dfvar=clima_VS_2, 
                                         colvarxy=1:2,colvar="all",resolution= 0.041665))


#7. Agregar una nueva columna para diferenciar los datos de presencia (1) vs. el background (0) para luego unir los dos set de datos para cada especie.

occ_vs1 <- cbind(occ_vs1,species_occ = 1)
clima_VS_1 <- cbind(clima_VS_1,species_occ = 0)

names(clima_VS_1)[1] = "lon" 
names(clima_VS_1)[2] = "lat" 
data_VS_1 <- rbind(occ_vs1, clima_VS_1)

#VS
occ_vs2 <- cbind(occ_vs2,species_occ = 1)
clima_VS_2 <- cbind(clima_VS_2,species_occ = 0)

names(clima_VS_2)[1] = "lon" 
names(clima_VS_2)[2] = "lat" 
data_VS_2 <- rbind(occ_vs2, clima_VS_2)


##8. Calcula los valores de PCA para las variables que se est?n comparando al analizar los nichos de las dos especies
pca.env <- dudi.pca(rbind(data_VS_1,data_VS_2)[,3:9],center = T, scale = T, scannf=F,nf=2)
ecospat.plot.contrib(contrib=pca.env$co, eigen=pca.env$eig)

##Rescata los valores "scores" de los ejes para los componentes obtenidos
scores.globclim <- pca.env$li

scores.vs1 <- suprow(pca.env,data_VS_1[which(data_VS_1[,10]==1),3:9])$li
scores.vs2 <- suprow(pca.env,data_VS_2[which(data_VS_2[,10]==1),3:9])$li
scores.clim.vs1 <- suprow(pca.env,data_VS_1[which(data_VS_1[,10]==0),3:9])$li
scores.clim.vs2 <- suprow(pca.env,data_VS_2[which(data_VS_2[,10]==0),3:9])$li

####################
######HUMBOLDT######
####################

z_1<- humboldt.grid.espace(scores.globclim,
                           scores.clim.vs1,
                           scores.vs1,
                           kern.smooth=1,R=200)

z_2<- humboldt.grid.espace(scores.globclim,
                           scores.clim.vs2,
                           scores.vs2,
                           kern.smooth=1,R=200)

par(mfrow=c(1,2))
plot.vs1 <- humboldt.plot.niche(z_1, title = "YP/MPS-VS", name.axis1 = "PC1", name.axis2 = "PC2", 
                                correct.env = F, color.ramp = 1)

plot.vs2 <- humboldt.plot.niche(z_2, title = "YP-VS", name.axis1 = "PC1", name.axis2 = "PC2", 
                                correct.env = F, color.ramp = 1)

#9. Calcular las gr?ficas de densidad de ocurrencia 
#Species1
grid.clim.vs1 <- ecospat.grid.clim.dyn(glob= scores.globclim,
                                       glob1= scores.clim.vs1,
                                       sp=scores.vs1, R=200,
                                       th.sp=0)
#Species2
grid.clim.vs2 <- ecospat.grid.clim.dyn(glob=scores.globclim,
                                       glob1=scores.clim.vs2,
                                       sp=scores.vs2, R=200,
                                       th.sp=0)

###Calculate Niche Overlap with ecospat.niche.overlap()Compute Schoener's D, index of niche overlap
D.overlap <- ecospat.niche.overlap (grid.clim.vs1, grid.clim.vs2, cor= T)
D.overlap

#10. Delimitat los nichos de las especies y cuantificar el grado de sobrelape.
niche.dyn <- ecospat.niche.dyn.index (grid.clim.vs1, grid.clim.vs2, intersection = NA)

##ver las graficas de los nichos por separado
ecospat.plot.niche (grid.clim.vs1, title="YP/MPS-VS", name.axis1="PC1", name.axis2="PC2", cor=F)
ecospat.plot.niche (grid.clim.vs2, title="YP-VS", name.axis1="PC1", name.axis2="PC2", cor=T)

########################################################################
################### Background niche similarity tests ######################
########################################################################
sim.test_sp1 <- ecospat.niche.similarity.test(grid.clim.vs1, grid.clim.vs2,
                                              rep = 1000,
                                              rand.type=1)

sim.test_sp2 <- ecospat.niche.similarity.test(grid.clim.vs2, grid.clim.vs1,
                                              rep = 1000,
                                              rand.type=1)

par(mfrow=c(1,2))
ecospat.plot.overlap.test(sim.test_sp1,"D","YP/MPS-VS vs. YP-VS")
ecospat.plot.overlap.test(sim.test_sp2,"D", "YP-VS vs. YP/MPS-VS")




#######################################################################
################### YP_MPS VS MPS ################################
#######################################################################
##1. Leer los archivos .CSV de los puntos de ocurrencia de las especies.
setwd("~")
VS_1 <- read.csv("YP_MPS_VS.csv", header = T, sep = ",")
VS_1$lat <- as.numeric(as.character(VS_1$lat))
VS_1$lon <- as.numeric(as.character(VS_1$lon))

VS_2 <- read.csv("MPS_VS.csv", header = T, sep = ",")
VS_2$lat <- as.numeric(as.character(VS_2$lat))
VS_2$lon <- as.numeric(as.character(VS_2$lon))

##2. Leer los archivos de las condiciones clim?ticas de la M de cada especie.
##especie 
setwd("~/MPS_YP/MPS") 
variables_VS1 <- list.files(".",pattern = "*.asc$",full.names = T)
varclim_vs1 <- stack(variables_VS1)

##especie virtual
setwd("~/MPS/MPS_YP") 
variables_VS2 <- list.files(".",pattern = "*.asc$",full.names = T) 
varclim_vs2 <- stack(variables_VS2)

##3. Crear un data frame de las condiciones clim?ticas del ?rea M de cada especie.
clim_punto_sp1 <- rasterToPoints(varclim_vs1[[1]], fun=NULL, spatial=TRUE)
clim_punto_sp2 <- rasterToPoints(varclim_vs2[[1]], fun=NULL, spatial=TRUE)

##4. Extraer los valores de todas las variables ambientales para cada punto de la M
##species 1
clima_VS_1 <- extract(varclim_vs1, clim_punto_sp1)
clima_VS_1 <- data.frame(coordinates(clim_punto_sp1),clima_VS_1)
clima_VS_1 <- subset(clima_species,!is.na(bio2) & !is.na(bio3) & !is.na(bio14) & !is.na(bio15) & !is.na(bio18) & !is.na(bio19))

##species 2
clima_VS_2 <- extract(varclim_vs2, clim_punto_sp2)
clima_VS_2 <- data.frame(coordinates(clim_punto_sp2),clima_VS2)
clima_VS_2 <- subset(clima_VS_2, !is.na(bio2) & !is.na(bio3) & !is.na(bio14) & !is.na(bio15) & !is.na(bio18) & !is.na(bio19))

#5. Seleccionar los datos de presencia de las especies pero considerando solo las coordenadas geogr?ficas.
occ.sp1 <- VS_1[2:3]
occ.sp2 <- VS_2[2:3]

#6.Integrar la informaci?n de datos de ocurrencia y los datos del background para cada especie. Adem?s limpiar los datos a la resoluci?n espacial deseada (en este caso 1km)
occ_vs1<-na.exclude(ecospat.sample.envar(dfsp=occ.sp1,colspxy=1:2, 
                                         colspkept=1:2,dfvar=clima_VS_1,
                                         colvarxy=1:2,colvar="all",resolution= 0.041665))

occ_vs2<-na.exclude(ecospat.sample.envar(dfsp=occ.sp2,colspxy=1:2,
                                         colspkept=1:2,dfvar=clima_VS_2, 
                                         colvarxy=1:2,colvar="all",resolution= 0.041665))


#7. Agregar una nueva columna para diferenciar los datos de presencia (1) vs. el background (0) para luego unir los dos set de datos para cada especie.

occ_vs1 <- cbind(occ_vs1,species_occ = 1)
clima_VS_1 <- cbind(clima_VS_1,species_occ = 0)

names(clima_VS_1)[1] = "lon" 
names(clima_VS_1)[2] = "lat" 
data_VS_1 <- rbind(occ_vs1, clima_VS_1)

#VS
occ_vs2 <- cbind(occ_vs2,species_occ = 1)
clima_VS_2 <- cbind(clima_VS_2,species_occ = 0)

names(clima_VS_2)[1] = "lon" 
names(clima_VS_2)[2] = "lat" 
data_VS_2 <- rbind(occ_vs2, clima_VS_2)

##8. Calcula los valores de PCA para las variables que se est?n comparando al analizar los nichos de las dos especies
pca.env <- dudi.pca(rbind(data_VS_1,data_VS_2)[,3:8],center = T, scale = T, scannf=F,nf=2)
ecospat.plot.contrib(contrib=pca.env$co, eigen=pca.env$eig)

##Rescata los valores "scores" de los ejes para los componentes obtenidos
scores.globclim <- pca.env$li
scores.vs1 <- suprow(pca.env,data_VS_1[which(data_VS_1[,9]==1),3:8])$li
scores.vs2 <- suprow(pca.env,data_VS_2[which(data_VS_2[,9]==1),3:8])$li
scores.clim.vs1 <- suprow(pca.env,data_VS_1[which(data_VS_1[,9]==0),3:8])$li
scores.clim.vs2 <- suprow(pca.env,data_VS_2[which(data_VS_2[,9]==0),3:8])$li

####################
######HUMBOLDT######
####################

z_1<- humboldt.grid.espace(scores.globclim,
                           scores.clim.vs1,
                           scores.vs1,
                           kern.smooth=1,R=200)

z_2<- humboldt.grid.espace(scores.globclim,
                           scores.clim.vs2,
                           scores.vs2,
                           kern.smooth=1,R=200)

par(mfrow=c(1,2))
plot.sp<- humboldt.plot.niche(z_1, title = "YP/MPS-VS", name.axis1 = "PC1", name.axis2 = "PC2", 
                              correct.env = F, color.ramp = 1)

plot.vs <- humboldt.plot.niche(z_2, title = "MPS-VS", name.axis1 = "PC1", name.axis2 = "PC2", 
                               correct.env = F, color.ramp = 1)

#9. Calcular las gr?ficas de densidad de ocurrencia 
#Species1
grid.clim.vs1 <- ecospat.grid.clim.dyn(glob= scores.globclim,
                                       glob1= scores.clim.vs1,
                                       sp=scores.vs1, R=200,
                                       th.sp=0)
#Species2
grid.clim.vs2 <- ecospat.grid.clim.dyn(glob=scores.globclim,
                                       glob1=scores.clim.vs2,
                                       sp=scores.vs2, R=200,
                                       th.sp=0)

###Calculate Niche Overlap with ecospat.niche.overlap()Compute Schoener's D, index of niche overlap
D.overlap <- ecospat.niche.overlap (grid.clim.vs1, grid.clim.vs2, cor= T)
D.overlap

#10. Delimitat los nichos de las especies y cuantificar el grado de sobrelape.
niche.dyn <- ecospat.niche.dyn.index (grid.clim.vs1, grid.clim.vs2, intersection = NA)

########################################################################
################### Background niche similarity tests ######################
########################################################################
sim.test_sp1 <- ecospat.niche.similarity.test(grid.clim.vs1, grid.clim.vs2,
                                              rep = 1000, overlap.alternative= "higher",
                                              rand.type=1)

sim.test_sp2 <- ecospat.niche.similarity.test(grid.clim.vs2, grid.clim.vs1,
                                              rep = 1000, overlap.alternative ="higher",
                                              rand.type=1)
par(mfrow=c(1,2))
ecospat.plot.overlap.test(sim.test_sp1,"D","YP/MPS-VS vs. MPS-VS")
ecospat.plot.overlap.test(sim.test_sp2,"D", "MPS-VS vs. YP/MPS-VS")

###FIN

