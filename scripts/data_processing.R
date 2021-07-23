
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

write.csv(df, file = "spectral_library/spectral_library_hyperspec")

# synthmix_hyperspec.csv created with python function "snythmix"
sli_hyperspec <- read.csv("spectral_library/synthmix_hyperspec.csv")

sli_hyperspec_df <- sli_hyperspec %>% 
  dplyr::select(-class_ID) %>% 
  rename(class_ID = Unnamed..0)
 


# df <- as.data.frame(hyperspec_df) %>% 
#   select(-1) %>% 
#   na.omit() %>% 
#   select(-coords.x1, -coords.x2) 

names(df)


sli_snythmix <- sli
head(sli)
sli$wavelength <- c(493, 560, 665, 704, 740, 783, 883, 865, 1610, 2190)
sli_long <- sli %>% 
  pivot_longer(names_to = "type", values_to = "reflectance", cols = c(2:5))


ggplot(sli_long, aes(x=wavelength, y= reflectance, color=type)) +
  geom_line() +
  theme_classic()


#  Modeling fractional cover ----


# Define accuracy from 10-fold cross-validation as optimization measure
cv <- tune.control(cross = 3) # change for final 

# Use tune.svm() for a grid search of the gamma and cost parameters

# test code

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

# load in validation points from cooper et al 2020
validation <- readOGR(dsn='validation/validation_poly.shp')

prediction_hyperspec_NV <- predict(hyperspec, svm_NV_test)
plot(prediction_hyperspec_NV)

prediction_PV <- predict(s2, svm_test_PV)
plot(prediction_PV)

prediction_soil <- predict(s2, svm_test_soil)
plot(prediction_soil)

writeRaster(prediction_NPV, "data/gcg_eo_s09/prediction_NPV.tif", 
            datatype="INT1S", 
            overwrite=T)

# Based on this, can you identify fields in different stages of the crop 
# phenology?

# Based on this we can only easy identify crops in the senescence phenological
# stage, or that have already been harvested. The other colored areas could be a 
# mix of soil or photosynthetic vegetation. 


#############################################################################
# 4) Evaluation of fractional cover
#############################################################################

reference_data <- readOGR("data/gcg_eo_s09/s09_validation/Validation_scaled_20190726.shp")
prediction_stack <-  (stack(c(prediction_NPV, prediction_PV, prediction_soil)))
plot(prediction_stack)

# extracting our reference data points from prediction_NPV
predictions_stack <- raster::extract(prediction_stack, reference_data, sp=TRUE)
# nothing is negative or above one, so no adjustemts are needed here
predictions_df <- as.data.frame(predictions_stack) %>% 
  rename(NPV_predict = layer.1,
         PV_predict = layer.2,
         soil_predict = layer.3) %>% 
  mutate(total_cover = NPV_predict + PV_predict + soil_predict)

p1 <- ggplot(predictions_df, aes(NPV_predict, NPV_val)) +
  geom_point() +
  geom_abline(intercept = 0, slope = 1)

p2 <- ggplot(predictions_df, aes(PV_predict, PV_val)) +
  geom_point() +
  geom_abline(intercept = 0, slope = 1)

p3 <- ggplot(predictions_df, aes(soil_predict, soil_val)) +
  geom_point() +
  geom_abline(intercept = 0, slope = 1)

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



































































