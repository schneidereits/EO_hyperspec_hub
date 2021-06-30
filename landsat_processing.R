# 23.6.2021
# Data processing 

# Library ----

library(tidyverse)
library(rgdal)
library(raster)
library(randomForest)


## hyper spec ----

### Import and manipulation ----

# hyperspectral data 
hyperspec <- raster("O:/SS21_EO/Hyperspectral/data/enmap/2013SU_BA_ENMAP_L2SIM.bsq")
plot(hyperspec)

# extract training data point from hyperspec stack
hyperspec_df <- raster::extract(hyperspec, training_data, sp=T) 

# save as df and remove coords
df <- as.data.frame(hyperspec_df) %>% 
  dplyr::select(-coords.x1, -coords.x2)# %>% mutate(classID = as.factor(classID))

# creating narrow band indices ----


## multi-temporal data ----

# import files
landsat_files <- list.files("O:/SS21_EO/Hyperspectral/data/landsat/")
landsat_BOA_files <- landsat_files[grepl("BOA", landsat_files)]
landsat_QAI_files <- landsat_files[grepl("QAI", landsat_files)]

# QAI raster 
QAI <- raster(paste0("O:/SS21_EO/Hyperspectral/data/landsat/",c("20130102"), "_LEVEL2_LND07_QAI",".tif"))
plot(QAI)
(QAI_count <- freq(QAI) %>% as.data.frame() %>% arrange(-count))


# stack all the BOA
stack <- stack(paste0("O:/SS21_EO/Hyperspectral/data/landsat/", landsat_BOA_files[1:2])) # change to 60 or n 
plot(stack)



# Extent cropping 
(ext_scene_one <- extent(stack))
(ext_scene_two <- extent(hyperspec))
(common.extent <- intersect(ext_scene_one, ext_scene_two))
roi <- c(496665, 568665 , 4261665 , 4276665 )
multispec <- crop(stack, roi)
plot(multispec)


# Creating a cloud mask ----

# Define function to find fill values from Landsat BQA
fill_pixels <- function(x) {intToBits(x)[1] == T}

# high confidence clouds or high confidence cloud shadows or fill values
a_pixels <- function(x) {
  all(intToBits(x)[c(1)] == T) | 
    all(intToBits(x)[c(6,7)] == T) | 
    all(intToBits(x)[c(8,9)] == T)}

# high and medium confidence clouds or high and medium confidence cloud shadows or fill values 
b_pixels <- function(x) {
  all(intToBits(x)[c(1)] == T) | 
    all(intToBits(x)[c(6,7)] == T) | 
    intToBits(x)[c(7)] == T |
    all(intToBits(x)[c(8,9)] == T) | 
    intToBits(x)[c(9)] == T}

mask_a <- calc(QAI, fun = a_pixels)
mask_b <- calc(QAI, fun = b_pixels)

plot(mask_a)
plot(mask_b)


# synthetic endmember mixing

# Generating synthetic training data using the spectral lib




































