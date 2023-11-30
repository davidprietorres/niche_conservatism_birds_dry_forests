
# Installing and loading packages
if(!require(devtools)){
  install.packages("devtools")
}

if(!require(kuenm)){
  devtools::install_github("marlonecobos/kuenm")
}

library(kuenm)
library(dismo)
library(biomod2)
library(sp)
library(raster)
library(rgeos)
library(maptools)
library(rgdal)
library(usdm)
library(ENMeval)
library(foreign)
library(rJava)
library(spocc)
library(corrplot)
library(usdm)
library(XML)
library(ecospat)
library(dplyr)
library(reshape)

setwd("WORKING DIRECTORY") 

## Same script for Campylorhynchus yucatanicus, Geocooccyx velox, Amazilia rutila, Forpus cyanopygius y Momotus mexicanus
## change species name to run the model
data<- read.csv2("species.csv", sep = ",", header = TRUE)

species<-data$species 
lat<-as.numeric(data$lat)
lon<-as.numeric(data$lon)
datos<-data.frame(species,lon,lat) 


## Create layer stack
setwd("~")
pca_path <- list.files(".",pattern = "*.asc$",full.names = T)###crea el stack de las 19 variables climaticas del presente
capas_presente<- stack(pca_path)

## Extract environmental values for each occurence point
presencias_clima <- data.frame(raster::extract(capas_presente,datos[,2:3]))
presencias_clima2<-data.frame(datos,presencias_clima)
presencias_clima3 <- na.omit(presencias_clima2)

## Download world shp
data(wrld_simpl)
plot(wrld_simpl, xlim = c(-135,-75), ylim = c(25,70), axes = TRUE, col = "light blue")
points(presencias_clima3$lon, presencias_clima3$lat, col = "red", pch=20, cex= 0.9)


###################################################################
########### ENVIRONMENTAL VARIABLES SELECTION #####################
###################################################################

## Collinearity matrix and VIF estimation 
cormatriz <- cor(presencias_clima3[,4:22])
corrplot(cormatriz, outline = T, tl.col = "black", mar = c(2,0,1,1.5), title = "species", method = "shade")###gr?fica las correlaciones entre variables. OJO no guarda el archivo, hay que darle en la ventana y mandar a guardarlo

vif(presencias_clima3[,4:22])

## Omit correlated variables r < 0.8
nocor <-vifcor(presencias_clima3[,4:22], th= 0.8)

## Omit variables with VIF < 10
species<- vifstep(presencias_clima3[,4:22], th=10)

################################################################################
################################# DATA PARTITION ###############################
################################################################################

setwd("WORKING DIRECTORY") 

## Same script for Campylorhynchus yucatanicus, Geocooccyx velox, Amazilia rutila, Forpus cyanopygius y Momotus mexicanus
## change species name to run the model

data1<- read.csv2("species.csv", sep=",", header=TRUE)
names(data1)
data1$lat <- as.numeric(as.character(data1$lat))
data1$lon <- as.numeric(as.character(data1$lon))

species<-data1$species 
lat<-data1$lat
lon<-data1$lon
data2<-data.frame(species,lon,lat) 

todos <- unique(data2)
todos$sp <- paste(todos[,2], todos[,3], sep = "_")

###Data partition in 3/4 parts(=75%)
train <- todos[sample(nrow(todos), round((length(todos[,1])/4 *3))), ] 
test_ind <- todos[!todos[,4] %in% train[,4], ]

setwd("~")
write.csv(test_ind[,1:3], "species_ind.csv", row.names = FALSE)
write.csv(train[,1:3], "species_joint.csv", row.names = FALSE)
###nombre con el que guardará el archivo .csv que contiene el 80% de datos originales

#### First set (joint) division in test and training
rm(list=ls()) 
data1<- read.csv2("~/species_joint.csv", sep=",", header=TRUE) ###seleccionar el archivo que ocntine los datos de presencia de la especie a modelar
names(data1)

species<-data1$species 
lat<-data1$lat
lon<-data1$lon
data2<-data.frame(species,lon,lat)
data2$lat <- as.numeric(as.character(data2$lat))
data2$lon <- as.numeric(as.character(data2$lon))

## Random division of occurences data
todos <- unique(data2)
todos$sp <- paste(todos[,2], todos[,3], sep = "_")
train <- todos[sample(nrow(todos), round((length(todos[,1])/4 *3))), ]
test_test <- todos[!todos[,4] %in% train[,4], ]

setwd("~")

write.csv(train[,1:3], "~/species_train.csv", row.names = FALSE)###nombre con el que guardará el archivo .csv que contiene el 20% de datos independientes de evaluación
write.csv(test_test[,1:3], "~/species_test.csv", row.names = FALSE)###nombre con el que guardará el archivo .csv que contiene el 80% de datos originales

rm(list=ls()) 

############################################
############ kuenm #########################
############################################
setwd("~")### directorio donde estan todos los archivos .csv por separado de las especies con TODOS los datos originales
# If you have your own data and they are organized as in the first part of Figure 1, change 
# your directory and follow the instructions below.


############################################
######### The model calibration ############
############################################
# Variables with information to be used as arguments. Change "YOUR/DIRECTORY" by your actual directory.
occ_joint <- "species_joint.csv"
occ_tra <- "species_train.csv"
M_var_dir <- "M_variables"
batch_cal <- "Candidate_models"
out_dir <- "Candidate_models"
reg_mult <- c(seq(0.4, 2, 0.2), seq(2, 6, 0.5),8)
f_clas <- "all" 
args <- NULL # e.g., "maximumbackground=20000" for increasing the number of pixels in the bacground or


# note that some arguments are fixed in the function and should not be changed
maxent_path <- "~/R/win-library/4.3/dismo/java"
wait <- TRUE
run <- TRUE
kuenm_cal(occ.joint = occ_joint, occ.tra = occ_tra, M.var.dir = M_var_dir, batch = batch_cal,
          out.dir = out_dir, reg.mult = reg_mult, f.clas = f_clas, args = args,
          maxent.path = maxent_path, wait = wait, run = run)
warnings()


############################################
######### The model evaluation #############
############################################
occ_test <- "species_test.csv"
out_eval <- "Calibration_results"
threshold <- 5
rand_percent <- 50
iterations <- 500 
kept <- TRUE
selection <- "OR_AICc"
paral_proc <- FALSE # make this true to perform pROC calculations in parallel, recommended
# only if a powerfull computer is used (see function's help)
# Note, some of the variables used here as arguments were already created for previous function

cal_eval <- kuenm_ceval(path = out_dir, occ.joint = occ_joint, occ.tra = occ_tra, occ.test = occ_test, batch = batch_cal, out.eval = out_eval, threshold = threshold, rand.percent = rand_percent, iterations = iterations, kept = kept, selection = selection, parallel.proc = paral_proc)


############################################
######### The model creation ###############
############################################
batch_fin <- "Final_models2"
mod_dir <- "Final_Models"
rep_n <- 10
rep_type <- "Bootstrap"
jackknife <- FALSE
out_format <- "cloglog" 
project <- TRUE
G_var_dir <- "G_variables"
ext_type <- "all" 
write_mess <- TRUE
write_clamp <- TRUE
wait1 <- FALSE
run1 <- TRUE
args <- NULL # e.g., "maximumbackground=20000" for increasing the number of pixels in the bacground or
# "outputgrids=false" which avoids writing grids of replicated models and only writes the 
# summary of them (e.g., average, median, etc.) when rep.n > 1
# note that some arguments are fixed in the function and should not be changed
# Again, some of the variables used here as arguments were already created for previous functions


kuenm_mod(occ.joint = occ_joint, M.var.dir = M_var_dir, out.eval = out_eval, 
          batch = batch_fin, rep.n = rep_n,rep.type = rep_type, jackknife = jackknife,
          out.dir = mod_dir, out.format = out_format, project = project,
          G.var.dir = G_var_dir, ext.type = ext_type, write.mess = write_mess, 
          write.clamp = write_clamp, maxent.path = maxent_path, args = args,
          wait = wait1, run = run1)

############################################
####### The model evaluation final #########
############################################
occ_ind <- "species_ind.csv"
replicates <- TRUE
threshold <- 10
out_feval <- "Final_Models_evaluation"
# Most of the variables used here as arguments were already created for previous functions
fin_eval <- kuenm_feval(path = mod_dir, occ.joint = occ_joint, occ.ind = occ_ind, replicates = replicates,
                        out.eval = out_feval, threshold = threshold, rand.percent = rand_percent,
                        iterations = iterations, parallel.proc = paral_proc)

best <- read.csv("Calibration_results/best_candidate_models_OR_AICc.csv")
knitr::kable(best, caption = "Models selected based on significance, omission rates, and AICc, in that order.")
