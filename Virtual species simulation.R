
library(dismo)
library(raster)
library(rworldmap)
library(ade4)
library(maptools)
library(virtualspecies)
library(sp)
library(sf)


setwd("WORKING DIRECTORY")    

##### 1. CREATE STACK AND ENVIRONMENTAL SUBSETS #####
pca_path <- list.files(".",pattern = "*.asc$",full.names = T)
capas<- stack(pca_path)

wc <- getData("worldclim", var = "bio", res = 2.5)
wc <- crop(capas,extent(-116.3689,-78.71014, 8.53292, 32.71804)) 

wc1 <- wc[[c("BIO2","BIO3","BIO4","BIO9","BIO14","BIO15","BIO18","BIO19")]] 



##### REMOVE COLLINEARITY #####
removeCollinearity(wc1,
                   multicollinearity.cutoff = 0.8,
                   select.variables = T,
                   sample.points = F,
                   plot = T,
                   method = "pearson")


##### 2. ENVIRONMENTAL VALUES EXTRACTED FROM 1000 POINTS IN MESOAMERICA SDTF #####
param <- formatFunctions(BIO2 = c(fun = "dnorm", mean = 132, sd = 24.21),
                         BIO3 = c(fun = "dnorm", mean = 67, sd = 7.07),
                         BIO9 = c(fun = "dnorm", mean = 232, sd = 29.01),
                         BIO14 = c(fun = "dnorm", mean = 9, sd = 11.16),
                         BIO15 = c(fun = "dnorm", mean = 93, sd = 19.59),
                         BIO18 = c(fun = "dnorm", mean = 353, sd = 146.54))


##### 3. GENERATE SPECIES FROM FUNCTION #####
VS1 <- generateSpFromFun(raster.stack = wc1[[c("BIO2","BIO3","BIO9","BIO14","BIO15","BIO18")]],  
                         parameters = param, rescale= T, rescale.each.response = T)

plot(VS1)
plotResponse(VS1)


##### 4. PROBABILITY CONVERSION #####
sp1 <- convertToPA(VS1,
                   species.prevalence = 0.3,
                   plot = TRUE)
sp1
plotSuitabilityToProba(sp1)


##### 5. LIMIT DISTRIBUTION #####
#### LIMIT TO YUCATAN PENINSULA DRY FORESTS
Yucatan <- st_read("~/virtual species/YP.shp") 
YP <- as_Spatial(Yucatan)
YP_SDTF <- limitDistribution(sp1,
                             geographical.limit = "polygon",
                             area = Yucatan)

#### LIMIT TO YUCATAN PENINSULA AND MESOAMERICAN PACIFIC SLOPE DRY FORESTS
YPMPS <- st_read("~/virtual species/MPS_YP.shp") 
YP_MPS <- as_Spatial(YPMPS)
YP_MPS_SDTF<- limitDistribution(sp1,
                          geographical.limit = "polygon",
                          area = Bosques)

#### LIMIT TO MESOAMERICAN PACIFIC SLOPE DRY FORESTS
Pacific <- st_read("~/virtual species/MPS.shp") 
MPS <- as_Spatial(Pacific)
MPS_SDTF <- limitDistribution(sp1,
                             geographical.limit = "polygon",
                             area = Pacific)

color <- colorRampPalette(c("dark grey","green"))
plot("DESIRED PATCH", add= T)
par(mfrow = c(1, 2))
plot(sp1$pa.raster, main = "Theoretical distribution")
plot("DESIRED PATCH"$occupied.area, main = "Realised distribution", col = color(20))



##### 6. SAMPLING #####
PA.points <- sampleOccurrences(sp1,
                               n = 1000,
                               type = "presence-absence",
                               sampling.area = "DESIRED PATCH") ## Sampling area to YP, YP_MPS or MPS shapefile

PA.points$sample.points

occ <- PA.points$sample.points
points(occ[occ$Observed == 1, c("x", "y")], pch = 16, cex = .8)
points(occ[occ$Observed == 0, c("x", "y")], pch = 1, cex = .8)


##### 7. WRITE CSV FILE #####
setwd("MAIN DIRECTORY")
write.csv(PA.points$sample.points, file= "SDTF_patch.csv") 

