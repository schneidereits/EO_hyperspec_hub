# 23.6.2021
# Landsat data processing 

# Library ----

library(tidyverse)
library(rgdal)
library(raster)
library(randomForest)

# data manipulation ----

landsat_files <- list.files("O:/SS21_EO/Hyperspectral/data/landsat/")
landsat_BOA_files <- landsat_files[grepl("BOA", landsat_files)]

# stack all the temporal metrics
stack <- stack(paste0("O:/SS21_EO/Hyperspectral/data/landsat/",c("20130102"), "_LEVEL2_LND07_BOA",".tif"))
plot(stack)



