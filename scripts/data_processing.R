
# Library ----

library(tidyverse)
library(rgdal)
library(raster)
library(randomForest)
library(e1071)

# functions ----

# unclear if needed for cloudmasking
mask_highconf <- function(x){
  bs <- intToBits(x)
  return ( ((bs[1]) | (bs[6] & bs[7]) | (bs[8] & bs[9])) == T)
}

# function used to convert QA band into bits 
fill_pixels <- function(x) {((intToBits(x)[1] == T))}

# 
# for (i in length(QAI)) {
#   mask <- calc(QAI[[i]], fun = pixels) # function should probably be mask_highconf???
# }

## hyper spec ----

### Import and manipulation ----

# hyperspectral data 
#hyperspec <- stack("data/2013SU_BA_ENMAP_L2sim.bsq")
# create vector of names
names <- names(hyperspec)
# use sub to extract standard ".000000.Nanometers." string at end of col name
bands <- sub(".000000.Nanometers.", "", names) 
# remove string before standard '.*bsq...' string at begining of col name
bands <- paste0("wavelength_", c(sub('.*bsq...', '', bands)))

# assign the correct band names to the hyperspec stack
names(hyperspec) <- bands

writeRaster(hyperspec, "data/2013SU_BA_ENMAP_L2sim_renamed.bsq",
            # overwrite=TRUE,
            format='GTiff')

hyperspec <- stack("data/2013SU_BA_ENMAP_L2sim_renamed.tif")

plot(hyperspec)

# import our training data that we collected in QGIS
training_data <- readOGR(dsn='spectral_library/spectral_library.shp')

# extract training data point from hyperspec stack
hyperspec_df <- raster::extract(hyperspec, training_data, sp=T) 

# save as df and remove coords
df <- as.data.frame(hyperspec_df) %>% 
  dplyr::select(-coords.x1, -coords.x2)# %>% mutate(classID = as.factor(classID))

# creating narrow band indices ----


## multi-temporal data ----

# import files
# create list of files in the data directory 
landsat_BOA_files <- list.files("data/BOA/")[-54]
landsat_QAI_files <- list.files("data/QAI/")[-54]


# define extent
# extent of landsat stack
#old path on HU desktop
#stack <- stack(paste0("O:/SS21_EO/Hyperspectral/data/landsat/", landsat_QAI_files[1]))

# new path on local
QAI_stack <- stack(paste0("data/QAI/", QAI_files[i]))

(ext_scene_one <- extent(QAI_stack))
# extent of hyperspec stack
(ext_scene_two <- extent(hyperspec))
(common.extent <- intersect(ext_scene_one, ext_scene_two))
roi <- c(496665, 568665 , 4261665 , 4276665)

#loop for QAI raster
for (i in c(1:59)) {
  QAI_stack <- stack(paste0("data/QAI/", landsat_QAI_files[i]))
  QAI_stack <- crop(QAI_stack, roi)
  mask <- calc(QAI_stack, fun = fill_pixels)

  if (i==1) {
    final_mask_stack <- mask

  } else {

    final_mask_stack <- stack(final_mask_stack, mask)
    print(paste("moin layer", i, "is done"))

  }
}

plot(final_mask_stack)

writeRaster(final_mask_stack, "data/landsat_cloud_mask_20072021",
          # overwrite=TRUE,
            format='GTiff')

# import created final mask stack
# original computed on HU desktop
final_mask_stack <- stack("data/landsat_cloud_mask.tif")
# secound run computed on local machine 
final_mask_stack_2 <- stack("data/landsat_cloud_mask_20072021.tif")
# visual assesment of masks
plot(final_mask_stack[[45:55]])
plot(final_mask_stack_2[[45:55]])


# loop used to mask landsat images
for (i in c(1:59)) {
  # stack each image in directory 
  BOA_stack <- stack(paste0("data/BOA/", landsat_BOA_files[i]))

  # crop size
  BOA_stack <- crop(BOA_stack, roi)
  
  if(i==1)  {
    # apply mask
    final_masked_BOA <- mask(BOA_stack,
                             mask = final_mask_stack[[i]], 
                             maskvalue = 1) #takes everything that is a 1 in our mask and makes an NA
    
  } else { 
    
    masked_BOA <-  mask(BOA_stack,
                        mask = final_mask_stack[[i]], 
                        maskvalue =1) #takes everything that is a 1 in our mask and makes an NA
    # add mask of individual image to stack of already masked BOAs
    final_masked_BOA <- stack(final_masked_BOA, masked_BOA)
    print(paste("layer", i, "is done"))
    
  }
}

writeRaster(final_masked_BOA, "data/masked_BOA_21072021",
            #overwrite=TRUE,
            format='GTiff')

landsat_BOA <- stack("data/masked_BOA_21072021.tif")

# filter reluctance values from BOA images 
 landsat_BOA[(landsat_BOA>10000) | (landsat_BOA<0)] <- NA

writeRaster(landsat_BOA, "data/landsat_BOA_filtered_21072021.tif",
            #overwrite=TRUE,
            format='GTiff')

landsat_BOA <- stack("data/landsat_BOA_filtered_21072021.tif")

plot(landsat_BOA[[264:276]])

# relic code used to count the frequency of of QA pixel types
# QAI <- stack(paste0("O:/SS21_EO/Hyperspectral/data/landsat/", landsat_QAI_files[59:60]))
# plot(QAI)
# (QAI_count <- freq(QAI) %>% as.data.frame() %>% arrange(-count))
# 

# tassel cap transformation ----

# TC Coefficients from Crist (1985) (COEF NEED TO BE CHECKED )
tcc <- matrix(c( 0.2043,  0.4158,  0.5524, 0.5741,  0.3124,  0.2303, 
                 -0.1603, -0.2819, -0.4934, 0.7940, -0.0002, -0.1446,
                 0.0315,  0.2021,  0.3102, 0.1594, -0.6806, -0.6109), 
              # assign names to bands 
              dimnames = list(
                c('blue', 'green', 'red', 'nIR', 'swIR1', 'swIR2'),
                c('brightness', 'greenness', 'wetness')), ncol = 3)
# visual check
print(tcc)

# reduced dataset for testing
# landsat_reduced <- landsat_BOA[[1:12]]
# 
# landsat_reduced <- stack(I_imgs[[1]])
# hist(landsat_reduced) # for outlier check
#


# loop for tcb

for (i in c(1:59)) {
  
  # separate handeling to create final final_landsat_BOA_tcb stack
  # note the loop is performed in iterations of 6 for each landsat image (once for each band)
  if (i==1) {
    
    loop_images <- c(1:6)
    
    tcb <- sum(landsat_BOA[[loop_images]] * tcc[,1])
    
    
    final_landsat_BOA_tcb <- stack(tcb)
    
    print("tasseled cap brightness 1 is done")
    
  } else {
    # logic for loading correct bands (per image) to be converted to tc
    loop_images <-c((((i-1)*6)+1):((i)*6))
    
    tcb <- sum(landsat_BOA[[loop_images]] * tcc[,1])
    
    
    landsat_BOA_tcb <- stack(tcb)
    
    final_landsat_BOA_tcb <- stack(final_landsat_BOA_tcb, landsat_BOA_tcb)
    
    print(paste0("tasseled cap brightness ", i,  "is done"))
    
  }
}

writeRaster(final_landsat_BOA_tcb, "data/tcb_21072021.tif",
            #overwrite=TRUE,
            format='GTiff')

# loop for tcg

for (i in c(1:59)) {
  
  if (i==1) {
    
    loop_images <- c(1:6)
    
    tcg <- sum(landsat_BOA[[loop_images]] * tcc[,2])
    
    
    final_landsat_BOA_tcg <- stack(tcg)
    
    print("tassel cap greeness 1 is done")
    
  } else {
    
    loop_images <-c((((i-1)*6)+1):((i)*6))
    
    tcg <- sum(landsat_BOA[[loop_images]] * tcc[,2])
    
    
    landsat_BOA_tcg <- stack(tcg)
    
    final_landsat_BOA_tcg <- stack(final_landsat_BOA_tcg, landsat_BOA_tcg)
    
    print(paste0("tasseled cap greeness ", i,  "is done"))
    
  }
}

writeRaster(final_landsat_BOA_tcg, "data/tcg_21072021.tif",
            #overwrite=TRUE,
            format='GTiff')
# loop for tcw

for (i in c(1:59)) {
  
  if (i==1) {
    
    
    
    loop_images <- c(1:6)
    
    tcw <- sum(landsat_BOA[[loop_images]] * tcc[,3])
    
    
    final_landsat_BOA_tcw <- stack(tcw)
    
    print("tassel cap wettness 1 is done")
    
  } else {
    
    loop_images <-c((((i-1)*6)+1):((i)*6))
    
    tcw <- sum(landsat_BOA[[loop_images]] * tcc[,3])
    
    
    landsat_BOA_tcw <- stack(tcw)
    
    final_landsat_BOA_tcw <- stack(final_landsat_BOA_tcw, landsat_BOA_tcw)
    
    print(paste0("tassel cap wettness ", i,  "is done"))
    
  }
}

writeRaster(final_landsat_BOA_tcw, "data/tcw_21072021.tif",
            #overwrite=TRUE,
            format='GTiff')



# Creating spectral temporal metrics ----

final_landsat_BOA_tcb <- stack("data/tcb_21072021.tif")
final_landsat_BOA_tcg <- stack("data/tcg_21072021.tif")
final_landsat_BOA_tcw <- stack("data/tcw_21072021.tif")

# convert to matrix, inorder to have a final output of one layer per STM
tcb_matrix <- as.matrix(final_landsat_BOA_tcb)
tcg_matrix <- as.matrix(final_landsat_BOA_tcg)
tcw_matrix <- as.matrix(final_landsat_BOA_tcw) 


# tcb stm ----

# Calculate IQR across rows in matrix
# na.rm = T in order to prevent all values to become NA
tcb_matrix_IQR <- apply(tcb_matrix,1, FUN=IQR, na.rm=T)

# Write results to empty raster
tcb_matrix_IQR_raster <- raster(nrows=final_landsat_BOA_tcb@nrows, 
                                ncols=final_landsat_BOA_tcb@ncols, 
                                crs=final_landsat_BOA_tcb@crs, 
                                vals=tcb_matrix_IQR,
                                ext=extent(final_landsat_BOA_tcb))

# Calculate mean across rows in matrix
tcb_matrix_mean <- apply(tcb_matrix,1, FUN=mean, na.rm=T)

# Write results to empty raster
tcb_matrix_mean_raster <- raster(nrows=final_landsat_BOA_tcb@nrows, 
                                 ncols=final_landsat_BOA_tcb@ncols, 
                                 crs=final_landsat_BOA_tcb@crs, 
                                 vals=tcb_matrix_mean,
                                 ext=extent(final_landsat_BOA_tcb))

# Calculate median across rows in matrix
tcb_matrix_median <- apply(tcb_matrix,1, FUN=median, na.rm=T)

# Write results to empty raster
tcb_matrix_median_raster <- raster(nrows=final_landsat_BOA_tcb@nrows, 
                                   ncols=final_landsat_BOA_tcb@ncols, 
                                   crs=final_landsat_BOA_tcb@crs, 
                                   vals=tcb_matrix_median,
                                   ext=extent(final_landsat_BOA_tcb))


# Calculate SD across rows in matrix
tcb_matrix_sd <- apply(tcb_matrix,1, FUN=sd, na.rm=T)

# Write results to empty raster
tcb_matrix_sd_raster <- raster(nrows=final_landsat_BOA_tcb@nrows, 
                                   ncols=final_landsat_BOA_tcb@ncols, 
                                   crs=final_landsat_BOA_tcb@crs, 
                                   vals=tcb_matrix_sd,
                                   ext=extent(final_landsat_BOA_tcb))

# Calculate max across rows in matrix
tcb_matrix_max <- apply(tcb_matrix,1, FUN=max, na.rm=T)

# Write results to empty raster
tcb_matrix_max_raster <- raster(nrows=final_landsat_BOA_tcb@nrows, 
                                ncols=final_landsat_BOA_tcb@ncols, 
                                crs=final_landsat_BOA_tcb@crs, 
                                vals=tcb_matrix_max,
                                ext=extent(final_landsat_BOA_tcb))

# Calculate min across rows in matrix
tcb_matrix_min <- apply(tcb_matrix,1, FUN=min, na.rm=T)

# Write results to empty raster
tcb_matrix_min_raster <- raster(nrows=final_landsat_BOA_tcb@nrows, 
                                ncols=final_landsat_BOA_tcb@ncols, 
                                crs=final_landsat_BOA_tcb@crs, 
                                vals=tcb_matrix_min,
                                ext=extent(final_landsat_BOA_tcb))


tcb_stm <- stack(tcb_matrix_IQR_raster, tcb_matrix_mean_raster,
                 tcb_matrix_median_raster, tcb_matrix_sd_raster,
                 tcb_matrix_max_raster, tcb_matrix_min_raster)

writeRaster(tcb_stm, "data/tcb_stm_21072021.tif",
                       #overwrite=TRUE,
                       format='GTiff')

tcb_stm_stack <- stack("data/tcb_stm_21072021.tif")
# assign STM names to layers
names(tcb_stm_stack) <- c("IQR", "mean", "median", 
                          "Std.dev", "max", "min")
plot(tcb_stm_stack)

# tcg stm ----
# Calculate p25 across rows in matrix
tcg_matrix_p25 <- apply(tcg_matrix,1, FUN=p25, na.rm=T)

# Write results to empty raster
tcb_matrix_p25_raster <- raster(nrows=final_landsat_BOA_tcb@nrows, 
                                ncols=final_landsat_BOA_tcb@ncols, 
                                crs=final_landsat_BOA_tcb@crs, 
                                vals=tcb_matrix_p25,
                                ext=extent(final_landsat_BOA_tcb))

# Calculate p75 across rows in matrix
tcg_matrix_p75 <- apply(tcg_matrix,1, FUN=p75, na.rm=T)

# Write results to empty raster
tcg_matrix_p75_raster <- raster(nrows=final_landsat_BOA_tcb@nrows, 
                                ncols=final_landsat_BOA_tcb@ncols, 
                                crs=final_landsat_BOA_tcb@crs, 
                                vals=tcg_matrix_p75,
                                ext=extent(final_landsat_BOA_tcb))

# Calculate IQR across rows in matrix
tcg_matrix_IQR <- apply(tcg_matrix,1, FUN=IQR, na.rm=T)

# Write results to empty raster
tcg_matrix_IQR_raster <- raster(nrows=final_landsat_BOA_tcb@nrows, 
                                ncols=final_landsat_BOA_tcb@ncols, 
                                crs=final_landsat_BOA_tcb@crs, 
                                vals=tcg_matrix_IQR,
                                ext=extent(final_landsat_BOA_tcb))

# Calculate mean across rows in matrix
tcg_matrix_mean <- apply(tcg_matrix,1, FUN=mean, na.rm=T)

# Write results to empty raster
tcg_matrix_mean_raster <- raster(nrows=final_landsat_BOA_tcb@nrows, 
                                 ncols=final_landsat_BOA_tcb@ncols, 
                                 crs=final_landsat_BOA_tcb@crs, 
                                 vals=tcg_matrix_mean,
                                 ext=extent(final_landsat_BOA_tcb))

# Calculate median across rows in matrix
tcg_matrix_median <- apply(tcg_matrix,1, FUN=median, na.rm=T)

# Write results to empty raster
tcg_matrix_median_raster <- raster(nrows=final_landsat_BOA_tcb@nrows, 
                                   ncols=final_landsat_BOA_tcb@ncols, 
                                   crs=final_landsat_BOA_tcb@crs, 
                                   vals=tcg_matrix_median,
                                   ext=extent(final_landsat_BOA_tcb))

# Calculate SD across rows in matrix
tcg_matrix_sd <- apply(tcg_matrix,1, FUN=sd, na.rm=T)

# Write results to empty raster
tcg_matrix_sd_raster <- raster(nrows=final_landsat_BOA_tcb@nrows, 
                               ncols=final_landsat_BOA_tcb@ncols, 
                               crs=final_landsat_BOA_tcb@crs, 
                               vals=tcg_matrix_sd,
                               ext=extent(final_landsat_BOA_tcb))

# Calculate max across rows in matrix
tcg_matrix_max <- apply(tcg_matrix,1, FUN=max, na.rm=T)

# Write results to empty raster
tcg_matrix_max_raster <- raster(nrows=final_landsat_BOA_tcb@nrows, 
                                ncols=final_landsat_BOA_tcb@ncols, 
                                crs=final_landsat_BOA_tcb@crs, 
                                vals=tcg_matrix_max,
                                ext=extent(final_landsat_BOA_tcb))


# Calculate min across rows in matrix
tcg_matrix_min <- apply(tcg_matrix,1, FUN=min, na.rm=T)

# Write results to empty raster
tcg_matrix_min_raster <- raster(nrows=final_landsat_BOA_tcb@nrows, 
                                ncols=final_landsat_BOA_tcb@ncols, 
                                crs=final_landsat_BOA_tcb@crs, 
                                vals=tcg_matrix_min,
                                ext=extent(final_landsat_BOA_tcb))


tcg_stm <- stack(tcg_matrix_IQR_raster, tcg_matrix_mean_raster,
                 tcg_matrix_median_raster, tcg_matrix_sd_raster,
                 tcg_matrix_max_raster, tcg_matrix_min_raster)

writeRaster(tcg_stm, "data/tcg_stm_21072021.tif",
                       #overwrite=TRUE,
                       format='GTiff')

tcg_stm_stack <- stack("data/tcg_stm_21072021.tif")
names(tcg_stm_stack) <- c("IQR", "mean", "median", 
                          "Std.dev", "max", "min")
plot(tcg_stm_stack)

# tcw stm ----
# Calculate p25 across rows in matrix
tcw_matrix_p25 <- apply(tcw_matrix,1, FUN=p25, na.rm=T)

# Write results to empty raster
tcw_matrix_p25_raster <- raster(nrows=final_landsat_BOA_tcb@nrows, 
                                ncols=final_landsat_BOA_tcb@ncols, 
                                crs=final_landsat_BOA_tcb@crs, 
                                vals=tcw_matrix_p25,
                                ext=extent(final_landsat_BOA_tcb))

# Calculate p75 across rows in matrix
tcw_matrix_p75 <- apply(tcw_matrix,1, FUN=p75, na.rm=T)

# Write results to empty raster
tcw_matrix_p75_raster <- raster(nrows=final_landsat_BOA_tcb@nrows, 
                                ncols=final_landsat_BOA_tcb@ncols, 
                                crs=final_landsat_BOA_tcb@crs, 
                                vals=tcw_matrix_p75,
                                ext=extent(final_landsat_BOA_tcb))

# Calculate IQR across rows in matrix
tcw_matrix_IQR <- apply(tcw_matrix,1, FUN=IQR, na.rm=T)

# Write results to empty raster
tcw_matrix_IQR_raster <- raster(nrows=final_landsat_BOA_tcb@nrows, 
                                ncols=final_landsat_BOA_tcb@ncols, 
                                crs=final_landsat_BOA_tcb@crs, 
                                vals=tcw_matrix_IQR,
                                ext=extent(final_landsat_BOA_tcb))

# Calculate mean across rows in matrix
tcw_matrix_mean <- apply(tcw_matrix,1, FUN=mean, na.rm=T)

# Write results to empty raster
tcw_matrix_mean_raster <- raster(nrows=final_landsat_BOA_tcb@nrows, 
                                 ncols=final_landsat_BOA_tcb@ncols, 
                                 crs=final_landsat_BOA_tcb@crs, 
                                 vals=tcw_matrix_mean,
                                 ext=extent(final_landsat_BOA_tcb))

# Calculate median across rows in matrix
tcw_matrix_median <- apply(tcw_matrix,1, FUN=median, na.rm=T)

# Write results to empty raster
tcw_matrix_median_raster <- raster(nrows=final_landsat_BOA_tcb@nrows, 
                                   ncols=final_landsat_BOA_tcb@ncols, 
                                   crs=final_landsat_BOA_tcb@crs, 
                                   vals=tcw_matrix_median,
                                   ext=extent(final_landsat_BOA_tcb))


# Calculate SD across rows in matrix
tcw_matrix_sd <- apply(tcw_matrix,1, FUN=sd, na.rm=T)

# Write results to empty raster
tcw_matrix_sd_raster <- raster(nrows=final_landsat_BOA_tcb@nrows, 
                                   ncols=final_landsat_BOA_tcb@ncols, 
                                   crs=final_landsat_BOA_tcb@crs, 
                                   vals=tcw_matrix_sd,
                                   ext=extent(final_landsat_BOA_tcb))

# Calculate max across rows in matrix
tcw_matrix_max <- apply(tcw_matrix,1, FUN=max, na.rm=T)

# Write results to empty raster
tcw_matrix_max_raster <- raster(nrows=final_landsat_BOA_tcb@nrows, 
                                ncols=final_landsat_BOA_tcb@ncols, 
                                crs=final_landsat_BOA_tcb@crs, 
                                vals=tcw_matrix_max,
                                ext=extent(final_landsat_BOA_tcb))

# Calculate min across rows in matrix
tcw_matrix_min <- apply(tcw_matrix,1, FUN=min, na.rm=T)

# Write results to empty raster
tcw_matrix_min_raster <- raster(nrows=final_landsat_BOA_tcb@nrows, 
                                ncols=final_landsat_BOA_tcb@ncols, 
                                crs=final_landsat_BOA_tcb@crs, 
                                vals=tcw_matrix_min,
                                ext=extent(final_landsat_BOA_tcb))


tcw_stm <- stack(tcw_matrix_IQR_raster, tcw_matrix_mean_raster,
                 tcw_matrix_median_raster, tcw_matrix_sd_raster,
                 tcw_matrix_max_raster, tcw_matrix_min_raster)

writeRaster(tcw_stm, "data/tcw_stm_21072021.tif",
                      # overwrite=TRUE,
                       format='GTiff')

tcw_stm_stack <- stack("data/tcw_stm_21072021.tif")
names(tcw_stm_stack) <- c("IQR", "mean", "median", 
                          "Std.dev", "max", "min")
plot(tcw_stm_stack)

full_stm_stack <- stack(tcb_stm_stack, tcg_stm_stack, tcw_stm_stack)



writeRaster(full_stm_stack, "data/full_stm_stack.tif",
            #overwrite=TRUE,
            format='GTiff')

full_stm_stack <- stack("data/full_stm_stack.tif")
names <-  c("tcb_IQR", "tcb_mean", "tcb_median", 
            "tcb_Std.dev", "tcb_max", "tcb_min",
            "tcg_IQR", "tcg_mean", "tcg_median", 
            "tcg_Std.dev", "tcg_max", "tcg_min",
            "tcw_IQR", "tcw_mean", "tcw_median", 
            "tcw_Std.dev", "tcw_max", "tcw_min")

names(full_stm_stack) <- names

plot(full_stm_stack)

# synthetic endmember mixing ----

# read in training data for spectral library (collected in QGIS)
training_data <- readOGR(dsn='spectral_library/spectral_library.shp')


# Use extract() to create a data.frame with training points as rows, and class labels 
# (classID) as well as the spectral bands of your composites as columns. 
# Remove the day of year and year flags (band 7 and 8) for the next steps.

# extracting hyperspec points
hyperspec_df <- raster::extract(hyperspec, training_data, sp=T) 

df <- as.data.frame(hyperspec_df) %>% 
  na.omit() %>% 
  # remove coordinates 
  dplyr::select(-coords.x1, -coords.x2)
# colnames(df[,2:198]) <- str_sub(colnames(df[,2:198]), -23, -18)   

#sub(".Nanometers.", "", colnames(df))
#sub(".*...", "", colnames(df))

write.csv(df, file = "spectral_library/spectral_library_hyperspec_extended")

# synthmix_hyperspec.csv created with python function "snythmix"
sli_hyperspec <- read.csv("spectral_library/spectral_library_hyperspec_extended.csv")

sli_hyperspec_df <- sli_hyperspec %>% 
  dplyr::select(-class_ID) %>% 
  rename(class_ID = Unnamed..0)
 
# for landsat



# extracting hyperspec points
landsat_df <- raster::extract(full_stm_stack, training_data, sp=T) 

df <- as.data.frame(landsat_df) %>% 
  na.omit() %>% 
  # remove coordinates 
  dplyr::select(-coords.x1, -coords.x2)
# colnames(df[,2:198]) <- str_sub(colnames(df[,2:198]), -23, -18)   

#sub(".Nanometers.", "", colnames(df))
#sub(".*...", "", colnames(df))

write.csv(df, file = "spectral_library/spectral_library_landsat")

# synthmix_hyperspec.csv created with python function "snythmix"
sli_landsat <- read.csv("spectral_library/synthmix_landsat.csv")

sli_landsat_df <- sli_landsat %>% 
  dplyr::select(-class_ID) %>% 
  rename(class_ID = Unnamed..0)



#  Modeling fractional cover ----

# hyperspec SVM ----

# Define accuracy from 10-fold cross-validation as optimization measure
cv <- tune.control(cross = 10) # change for final 

# Use tune.svm() for a grid search of the gamma and cost parameters


sli_hyperspec_grouped <- sli_hyperspec_df %>% 
  group_by(class_ID) %>% 
  group_split() %>% 
 map(., dplyr::select, -"class_ID")


# old code for only class_ID
# sli_hyperspec_NV <- sli_hyperspec_df %>% 
#   filter(class_ID == 5) %>% 
#   dplyr::select(-class_ID)

# svm_NV_test <-  svm(fraction~., 
#                     data = sli_hyperspec_NV, 
#                     kernel = 'radial',
#                     gamma = 1, 
#                     cost = 1, 
#                     epsilon = 0.001,
#                     tunecontrol = cv)

svm_list <- sli_hyperspec_grouped %>% map(~ svm(fraction~., 
                                              data        = .,
                                              kernel      = 'radial',
                                              gamma       = 1, 
                                              cost        = 1, 
                                              epsilon     = 0.001,
                                              tunecontrol = cv))


svm.tune_list <-  sli_hyperspec_grouped %>% map(~ tune.svm(fraction~., 
                                                           data        = as.data.frame(.), 
                                                           kernel      = 'radial', 
                                                           gamma       = (0.01:100), 
                                                           cost        = 10^(-2:2), 
                                                           epsilon     = 0.001,
                                                           tunecontrol = cv))

# Store the best model in a new object
svm.best <- svm.tune$best.model

# Which parameters performed best?
print(svm.best$gamma)
print(svm.best$cost)

classes <- list("conifer", "decid", "shrub", "NW", "NV")
for (i in 1:length(classes)) {
  
  svm <- svm_list[[i]]
  prediction_hyperspec <- predict(hyperspec, svm)
  classes[[i]] <- prediction_hyperspec
  print(paste("SVM Layer", i, "is done"))
}

prediction_conifer <- classes[[1]]
prediction_decid <- classes[[2]]
prediction_shrub <- classes[[3]]
prediction_NW <- classes[[4]]
prediction_NV <- classes[[5]]

prediction_stack_hyperspec <- stack(prediction_conifer, prediction_decid, prediction_shrub, 
                                    prediction_NW, prediction_NV)

writeRaster(prediction_stack_hyperspec, "data/prediction_hyperspec_extended.tif", 
           # datatype="INT1S", 
            #overwrite=T
           )

prediction_stack_hyperspec <- stack("data/prediction_hyperspec_extended.tif")
plot(prediction_stack_hyperspec)

# hyperspec random forest ----

# Train a randomForest() classification model with the data.frame created in the 
# prior step. Make sure to include only useful predictors.

sli_hyperspec_grouped <- sli_hyperspec_df %>% 
  group_by(class_ID) %>% 
  group_split() %>% 
  map(., dplyr::select, -"class_ID")

rf_test <- randomForest(fraction~., 
                          data = sli_hyperspec_grouped[[1]],
                          ntree= 1000)

rf_list <- sli_hyperspec_grouped %>% map(~ randomForest(fraction~., 
                                                        data = (.),
                                                        ntree= 1000))

saveRDS(rf_list, file = "data/hyperspec_rf_extended.rds")
readRDS(file = "data/hyperspec_rf_extended.rds")


# Define accuracy from 5-fold cross-validation as optimization measure
cv <- tune.control(cross = 10) 

# Use tune.randomForest to assess the optimal combination of ntree and mtry
rf.tune_list <- sli_hyperspec_grouped %>% map(~ randomForest(fraction~., 
                                                             data        = (.),
                                                             ntree       = 1000, 
                                                             mtry        = c(64:65), 
                                                             tunecontrol = cv))
saveRDS(rf.tune_list, file = "data/hyperspec_rf_tuned_extended.rds")
readRDS(file = "data/hyperspec_rf_tuned_extended.rds")

# Store the best model in a new object for further use
rf.best <- rf.tune$library(e1071)


# Is the parametrization and/or different from your previous model?
print(rf.tune)

# RF predictions

classes_rf <- list("conifer", "decid", "shrub", "NW", "NV")
for (i in 1:length(classes_rf)) {
  
  rf <- rf.tune_list[[i]]
  prediction_hyperspec <- predict(hyperspec, rf)
  classes_rf[[i]] <- prediction_hyperspec
  print(paste("rf", i, "is done"))
  
}

prediction_rf_conifer <- classes_rf[[1]]
prediction_rf_decid <- classes_rf[[2]]
prediction_rf_shrub <- classes_rf[[3]]
prediction_rf_NW <- classes_rf[[4]]
prediction_rf_NV <- classes_rf[[5]]

prediction_rf_tune_stack_hyperspec <- stack(prediction_rf_conifer, prediction_rf_decid, prediction_rf_shrub, 
                                       prediction_rf_NW, prediction_rf_NV)

writeRaster(prediction_rf_tune_stack_hyperspec, "data/prediction_rf_tune_hyperspec_extended.tif", 
            # datatype="INT1S", 
            overwrite=T
            )

prediction_rf_tune_stack_hyperspec_extended <- stack("data/prediction_rf_tune_hyperspec_extended.tif")
plot(prediction_rf_tune_stack_hyperspec_extended)

# landsat SVM ----

sli_landsat_grouped <- sli_landsat_df %>% 
  group_by(class_ID) %>% 
  group_split() %>% 
  map(., dplyr::select, -"class_ID")


# Define accuracy from 10-fold cross-validation as optimization measure
cv <- tune.control(cross = 10) 
svm_list <- sli_landsat_grouped %>% map(~ svm(fraction~., 
                                                data        = .,
                                                kernel      = 'radial',
                                                gamma       = 1, 
                                                cost        = 1, 
                                                epsilon     = 0.001,
                                                tunecontrol = cv))


svm.tune_list <-  sli_landsat_grouped %>% map(~ tune.svm(fraction~., 
                                                           data        = as.data.frame(.), 
                                                           kernel      = 'radial', 
                                                           gamma       = (0.01:100), 
                                                           cost        = 10^(-2:2), 
                                                           epsilon     = 0.001,
                                                           tunecontrol = cv))

# Store the best model in a new object
svm.best <- svm.tune$best.model

# Which parameters performed best?
print(svm.best$gamma)
print(svm.best$cost)


# SVM predictions
classes <- list("conifer", "decid", "shrub", "NW", "NV")
for (i in 1:length(classes)) {
  
  svm <- svm_list[[i]]
  prediction_landsat <- predict(full_stm_stack, svm)
  classes[[i]] <- prediction_landsat
  print(paste("SVM for class", i, "is done"))
}

prediction_conifer <- classes[[1]]
prediction_decid <- classes[[2]]
prediction_shrub <- classes[[3]]
prediction_NW <- classes[[4]]
prediction_NV <- classes[[5]]

prediction_stack_landsat <- stack(prediction_conifer, prediction_decid, prediction_shrub, 
                                    prediction_NW, prediction_NV)

writeRaster(prediction_stack_landsat, "data/prediction_landsat_extended.tif", 
            # datatype="INT1S", 
           # overwrite=T
            )

prediction_stack_landsat <- stack("data/prediction_landsat_extended.tif")
plot(prediction_stack_landsat)

# landsat random forest ----

# Train a randomForest() classification model with the data.frame created in the 
# prior step. Make sure to include only useful predictors.

sli_landsat_grouped <- sli_landsat_df %>% 
  group_by(class_ID) %>% 
  group_split() %>% 
  map(., dplyr::select, -"class_ID")

rf_test <- randomForest(fraction~., 
                        data = sli_landsat_grouped[[1]],
                        ntree= 1000)

rf_list <- sli_landsat_grouped %>% map(~ randomForest(fraction~., 
                                                        data = (.),
                                                        ntree= 1000))

saveRDS(rf_list, file = "data/landsat_rf_extended.rds")
readRDS(file = "data/landsat_rf.rds")


# Define accuracy from 5-fold cross-validation as optimization measure
cv <- tune.control(cross = 10) 

# Use tune.randomForest to assess the optimal combination of ntree and mtry
rf.tune_list <- sli_landsat_grouped %>% map(~ randomForest(fraction~., 
                                                             data        = (.),
                                                             ntree       = 1000, 
                                                             mtry        = c(17:18), 
                                                             tunecontrol = cv))

saveRDS(rf.tune_list, file = "data/landsat_rf_tuned_extended.rds")
readRDS(file = "data/landsat_rf_tuned.rds")

# Store the best model in a new object for further use
rf.best <- rf.tune$library(e1071)


# Is the parametrization and/or different from your previous model?
print(rf.tune)

# RF predictions

classes_rf <- list("conifer", "decid", "shrub", "NW", "NV")
for (i in 1:length(classes_rf)) {
  
  rf <- rf.tune_list[[i]]
  prediction_landsat <- predict(full_stm_stack, rf)
  classes_rf[[i]] <- prediction_landsat
  print(paste("rf", i, "is done"))
  
}

prediction_rf_conifer <- classes_rf[[1]]
prediction_rf_decid <- classes_rf[[2]]
prediction_rf_shrub <- classes_rf[[3]]
prediction_rf_NW <- classes_rf[[4]]
prediction_rf_NV <- classes_rf[[5]]

prediction_rf_tune_stack_landsat <- stack(prediction_rf_conifer, prediction_rf_decid, prediction_rf_shrub, 
                                            prediction_rf_NW, prediction_rf_NV)

writeRaster(prediction_rf_tune_stack_landsat, "data/prediction_rf_tuned_landsat_extended.tif", 
            # datatype="INT1S", 
            overwrite=T
            )

prediction_rf_tuned_landsat_extended <- stack("data/prediction_rf_tuned_landsat_extended.tif")
plot(prediction_rf_tuned_landsat_extended)

# difference betweem hyperspec and landsat tuned rf prediction
diff <- overlay(prediction_rf_tune_stack_hyperspec,
                prediction_rf_tune_stack_landsat,
        fun=function(r1, r2){return(abs(r1-r2))})
diff_mean <- mean(diff)
plot(diff)
plot(diff_mean)

writeRaster(diff, "data/prediction_rf_tuned_diff_extended.tif", 
            # datatype="INT1S", 
            #overwrite=T
            )

writeRaster(diff_mean, "data/prediction_rf_tuned_diff_mean_extended.tif", 
            # datatype="INT1S", 
            #overwrite=T
            )




# load in validation points from cooper et al 2020
validation <- readOGR(dsn='validation/validation_poly.shp')

prediction_landsat <- svm_list %>% map(~ full_stm_stack, predict((.)[1]))
plot(prediction_landsat)


writeRaster(prediction_NPV, "data/gcg_eo_s09/prediction_NPV.tif", 
            datatype="INT1S", 
            #overwrite=T
            )

# Based on this, can you identify fields in different stages of the crop 
# phenology?

# Based on this we can only easy identify crops in the senescence phenological
# stage, or that have already been harvested. The other colored areas could be a 
# mix of soil or photosynthetic vegetation. 


#Evaluation of fractional cover est----


reference_data <- readOGR("data/centroids.shp") # gis calculated centroids from sams validation polygons

transformed_data <- as.data.frame(reference_data) %>%  na.omit() %>% 
  mutate(across(c(Conifer:Other), as.numeric)) %>% 
#  filter(IgnoreSite != "y") %>% 
  mutate(MedShr = MedShr+ Ag,
         NW = UplGra+MngGra,
         NV = Soil+Imperv+Other) %>% 
  rename(conifer  = Conifer,
         decid = Broadleaf,
         shrub= MedShr) %>% 
  select(conifer, decid, shrub, NW, NV)

reference_data@data <- transformed_data


# extracting our reference data points from prediction_NPV
predictions_stack <- raster::extract(prediction_rf_tune_stack_hyperspec, reference_data, sp=TRUE) 

valdation_df <- as.data.frame(predictions_stack) %>% 
  rename(conifer_predict = prediction_rf_tune_hyperspec.1,
         decid_predict = prediction_rf_tune_hyperspec.2,
         shrub_predict = prediction_rf_tune_hyperspec.3,
         NW_predict = prediction_rf_tune_hyperspec.4,
         NV_predict = prediction_rf_tune_hyperspec.5) %>% 
  mutate(across(c(conifer_predict:NV_predict), ~ .*100),
         across(c(conifer_predict:NV_predict), round)) %>% 
  select(-coords.x1, -coords.x2)

valdation_df$dens_conifer <- get_density(x=valdation_df$conifer, y=valdation_df$conifer_predict, n=100)
valdation_df$dens_decid <- get_density(x=valdation_df$decid, y=valdation_df$decid_predict, n=100)
valdation_df$dens_shrub <- get_density(x=valdation_df$shrub, y=valdation_df$shrub_predict, n=100)
valdation_df$dens_NW <- get_density(x=valdation_df$NW, y=valdation_df$NW_predict, n=100)
valdation_df$dens_NV <- get_density(x=valdation_df$NV, y=valdation_df$NV_predict, n=100)

plot_validation <- function(data, class, class_predict){

  ggplot(data, aes(x=substitute(class), y=substitute(class_predict))) +
    geom_point()+
    xlim(0,100) +
    ylim(0,100) +
    scale_color_viridis()+
    geom_abline(aes(slope=1, intercept=0))+
    geom_smooth(method="lm")+
    annotate("text", x=1, y=100,
             label= paste0("R2 =",round((cor(substitute(data$class),substitute(data$class_predict), use="complete.obs")^2),2)),
             size=4.5, hjust=0)+
    annotate("text", x= 1, y=95, hjust=0, size=4.5,
             label= paste0("RMSE ==",round(sqrt(mean((substitute(data$class)-substitute(data$class_predict))^2, na.rm=TRUE)),2)),
             parse=TRUE)+
    annotate("text", x=1, y= 90, size=4.5, hjust=0,
             label=paste0("Bias = ", round((mean(substitute(data$class), na.rm=TRUE) - mean(substitute(data$class_predict),na.rm=TRUE)),2)))+
    annotate("text", x=1, y=85,size=4.5,hjust=0, label=paste0("MAE = ", round(mean(abs(substitute(data$class)- substitute(data$class_predict))),2))) 
  
  
  
  
}

plot_validation(valdation_df, conifer, conifer_predict)

 ggplot(valdation_df, aes(x=conifer, y=conifer_predict)) +
  geom_point()+
  xlim(0,100) +
  ylim(0,100) +
  scale_color_viridis()+
  geom_abline(aes(slope=1, intercept=0))+
  geom_smooth(method="lm")+
  annotate("text", x=1, y=100,
           label= paste0("R2 =",round((cor(valdation_df$conifer,valdation_df$conifer_predict, use="complete.obs")^2),2)),
           size=4.5, hjust=0)+
  annotate("text", x= 1, y=95, hjust=0, size=4.5,
           label= paste0("RMSE ==",round(sqrt(mean((valdation_df$conifer-valdation_df$conifer_predict)^2, na.rm=TRUE)),2)),
           parse=TRUE)+
  annotate("text", x=1, y= 90, size=4.5, hjust=0,
           label=paste0("Bias = ", round((mean(valdation_df$conifer, na.rm=TRUE) - mean(valdation_df$conifer_predict,na.rm=TRUE)),2)))+
  annotate("text", x=1, y=85,size=4.5,hjust=0, label=paste0("MAE = ", round(mean(abs(valdation_df$conifer - valdation_df$conifer_predict)),2))) 

ggplot(valdation_df, aes(x=decid, y=decid_predict)) +
  geom_point(#aes(color=dens_PV)
  )+
  xlim(0,100) +
  ylim(0,100) +
  scale_color_viridis()+
  geom_abline(aes(slope=1, intercept=0))+
  geom_smooth(method="lm")+
  annotate("text", x=1, y=100,
           label= paste0("R2 =",round((cor(valdation_df$decid,valdation_df$decid_predict, use="complete.obs")^2),2)),
           size=4.5, hjust=0.5)+
  annotate("text", x= -0.05, y=1, hjust=0, size=4.5,
           label= paste0("RMSE ==",round(sqrt(mean((valdation_df$decid-valdation_df$decid_predict)^2, na.rm=TRUE)),2)),
           parse=TRUE)+
  annotate("text", x=-0.05, y= 0.9, size=4.5, hjust=0,
           label=paste0("Bias = ", round((mean(valdation_df$decid, na.rm=TRUE) - mean(valdation_df$decid_predict,na.rm=TRUE)),2)))+
  annotate("text", x=-0.05, y=0.8,size=4.5,hjust=0, label=paste0("MAE = ", round(mean(abs(valdation_df$decid - valdation_df$decid_predict)),2))) 

ggplot(valdation_df, aes(x=shrub, y=shrub_predict)) +
  geom_point(#aes(color=dens_PV)
  )+
  scale_color_viridis()+
  geom_abline(aes(slope=1, intercept=0))+
  geom_smooth(method="lm")+
  annotate("text", x=-0.05, y=1.1,
           label= paste0("R2 =",round((cor(valdation_df$shrub,valdation_df$shrub_predict, use="complete.obs")^2),2)),
           size=4.5, hjust=0.5)+
  annotate("text", x= -0.05, y=1, hjust=0, size=4.5,
           label= paste0("RMSE ==",round(sqrt(mean((valdation_df$shrub-valdation_df$shrub_predict)^2, na.rm=TRUE)),2)),
           parse=TRUE)+
  annotate("text", x=-0.05, y= 0.9, size=4.5, hjust=0,
           label=paste0("Bias = ", round((mean(valdation_df$shrub, na.rm=TRUE) - mean(valdation_df$shrub_predict,na.rm=TRUE)),2)))+
  annotate("text", x=-0.05, y=0.8,size=4.5,hjust=0, label=paste0("MAE = ", round(mean(abs(valdation_df$shrub - valdation_df$shrub_predict)),2))) 

ggplot(valdation_df, aes(x=NW, y=NW_predict)) +
  geom_point(#aes(color=dens_PV)
  )+
  xlim(0,100) +
  ylim(0,100) +
  scale_color_viridis()+
  geom_abline(aes(slope=1, intercept=0))+
  geom_smooth(method="lm")+
  annotate("text", x=1, y=100,
           label= paste0("R2 =",round((cor(valdation_df$NW,valdation_df$NW_predict, use="complete.obs")^2),2)),
           size=4.5, hjust=0.5)+
  annotate("text", x= -0.05, y=1, hjust=0, size=4.5,
           label= paste0("RMSE ==",round(sqrt(mean((valdation_df$NW-valdation_df$NW_predict)^2, na.rm=TRUE)),2)),
           parse=TRUE)+
  annotate("text", x=-0.05, y= 0.9, size=4.5, hjust=0,
           label=paste0("Bias = ", round((mean(valdation_df$NW, na.rm=TRUE) - mean(valdation_df$NW_predict,na.rm=TRUE)),2)))+
  annotate("text", x=-0.05, y=0.8,size=4.5,hjust=0, label=paste0("MAE = ", round(mean(abs(valdation_df$NW - valdation_df$NW_predict)),2))) 

ggplot(valdation_df, aes(x=NV, y=NV_predict)) +
  geom_point(#aes(color=dens_PV)
  )+
  xlim(0,100) +
  ylim(0,100) +
  scale_color_viridis()+
  geom_abline(aes(slope=1, intercept=0))+
  geom_smooth(method="lm")+
  annotate("text", x=1, y=100,
           label= paste0("R2 =",round((cor(valdation_df$NV,valdation_df$NV_predict, use="complete.obs")^2),2)),
           size=4.5, hjust=0.5)+
  annotate("text", x= 1, y=80, hjust=0, size=4.5,
           label= paste0("RMSE ==",round(sqrt(mean((valdation_df$NV-valdation_df$NV_predict)^2, na.rm=TRUE)),2)),
           parse=TRUE)+
  annotate("text", x=-0.05, y= 0.9, size=4.5, hjust=0,
           label=paste0("Bias = ", round((mean(valdation_df$NV, na.rm=TRUE) - mean(valdation_df$NV_predict,na.rm=TRUE)),2)))+
  annotate("text", x=-0.05, y=0.8,size=4.5,hjust=0, label=paste0("MAE = ", round(mean(abs(valdation_df$NV - valdation_df$NV_predict)),2))) 



ggplot(valdation_df, aes(conifer, conifer_predict)) +
  geom_point() +
  xlim(0,100) +
  ylim(0,100) +
  geom_abline(intercept = 0, slope = 1) + 
  theme_classic()

ggplot(valdation_df, aes(decid, decid_predict)) +
  geom_point() +
  xlim(0,100) +
  ylim(0,100) +
  geom_abline(intercept = 0, slope = 1) + 
  theme_classic()

ggplot(valdation_df, aes(NV, NV_predict)) +
  geom_point() +
  xlim(0,100) +
  ylim(0,100) +
  geom_abline(intercept = 0, slope = 1) + 
  theme_classic()

library("gridExtra") 
grid.arrange(p1,p2,p3)

# Function that returns Root Mean Squared Error
rmse <- function(error)
{
  sqrt(mean(error^2))
}

# Function that returns Mean Absolute Error
mae <- function(error)
{
  mean(abs(error))
}

m1 <- lm(data = predictions_df, NPV_val ~ NPV_predict)
summary(m1)
# R2 0.6431 
error <- predictions_df$NPV_val - predictions_df$NPV_predict
error <-  na.omit(error)
unique(is.na(error))

(mae <- abs(mean(error))) #0.01277667

(rmse <- sqrt(mean(error^2))) # 0.1276346

m2 <- lm(data = predictions_df, PV_val ~ PV_predict)
summary(m2)
# R2 0.8933 
error <- predictions_df$PV_val - predictions_df$PV_predict
error <-  na.omit(error)
unique(is.na(error))

(mae <- abs(mean(error))) #0.03547932

(rmse <- sqrt(mean(error^2))) # 0.09200218

m3 <- lm(data = predictions_df, soil_val ~ soil_predict)
summary(m3)
# R2 0.6522 
error <- predictions_df$soil_val - predictions_df$soil_predict
error <-  na.omit(error)
unique(is.na(error))

(mae <- abs(mean(error))) #0.1728898

(rmse <- sqrt(mean(error^2))) # 0.2090371

# You can use this function to calculate the simple linear regression coefficients 
# of observed and predicted data and annotate it in your plots
lm_eqn <- function(pred_x,obs_y,df){
  m <- lm(obs_y ~ pred_x, df);
  eq <- substitute(y ==  a + b * x,
                   list(a = format(as.numeric(round(coef(m)[1]/100,2))),
                        b = format(as.numeric(round(coef(m)[2],2)))))
  as.character(as.expression(eq));
}

# For a better visualization, you can use this function to calculate the point
# density and display it as the color in your scatterplots:
get_density <- function(x, y,...) { # set x and y to predicted and observed values and n = 100
  dens <- MASS::kde2d(x, y, ...)
  ix <- findInterval(x, dens$x)
  iy <- findInterval(y, dens$y)
  ii <- cbind(ix, iy)
  return(dens$z[ii])
}



































































