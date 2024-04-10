##### Producing HK-wide rasters of LiDAR metrics (2020 data) #####
# Guidance for producing own metrics: https://rdrr.io/cran/lidR/src/R/metrics_stdmetrics.R

library(rgdal)
library(raster)
library(tmap)
library(tmaptools)
library(lidR)
library(terra)
library(stringr)

#for mosaicking rasters
library(raster)
library(terra)
library(tictoc)
library(dplyr)
library(tictoc)

##### *2020 LiDAR metrics - unthinned* #####

# use the hashed below to chunk tiles into overlapping segments and produce normalised las files

las_cat <- readLAScatalog("Z:/Data/2020_LiDAR/")
#las_cat <- readLAScatalog("Z:/Data/2020_LiDAR/")
HK1980 <- st_crs(2326) # set projection ('2326' is the HK1980 EPSG)
st_crs(las_cat) <- HK1980 # apply projection to las file
summary(las_cat)
las_check(las_cat)
plot(las_cat)

# multiples of 150 work for the HK LiDAR data. Enough overlap between tiles so that there aren't edge effects
#opt_chunk_size(las_cat) <- 150 # need to pick a value that is small enough for overlapping between tiles, but doesn't bust the compute power.
# opt_chunk_size(las_cat) <- 0 # need to pick a value that is small enough for overlapping between tiles, but doesn't bust the compute power.
# opt_chunk_buffer(las_cat) <- 15
# plot(las_cat, chunk_pattern = TRUE)
# 
# # parallel processing
# library(future)
# plan(multisession)
# 
# # normalise las files
# # first, create DTM
opt_output_files(las_cat) <- "D:/HK_LiDAR/2020/DTM/{ORIGINALFILENAME}" # Martha's Maxtor hard drive
#opt_output_files(las_cat) <- "/Users/marthaledger/OneDrive - The University Of Hong Kong/HKU/Butterflies/Data/HK_RS_data/2020/L2Processed_Martha/H4/Normal/norm_{XLEFT}_{YBOTTOM}"

opt_stop_early(las_cat) <- FALSE # bypasses errors
  #las_cat_normal <- normalize_height(las_cat, res = 1, knnidw(k = 10, p = 2), keep_lowest = FALSE)
  #las_cat_normal <- las_cat - rasterize_terrain(las_cat, res = 1, knnidw())
  las_cat_DTM <- rasterize_terrain(las_cat, res = 1, knnidw())
  # can opt to restart at the chunk which failed: requires being present at script though.
# opt_restart(ctg) <- 9

  path = "D:/HK_LiDAR/2020/DTM/"

  tiffcat <- function(f){
    do.call(rbind,
            lapply(f, function(r){
              ras = rast(r)
              poly = as.polygons(ext(ras), crs=crs(ras))
              poly$path = r
              poly
            })
    )
  }

  files = list.files(path, "*.tif", full=TRUE)
  tc = tiffcat(files)
#   
# library(stringr)
# DTM_names <- list.files("D:/HK_LiDAR/2020/DTM/")
# DTM_names <- str_replace(DTM_names, ".tif", "")
# remove .tif suffix

##### create normalised las files
# multiples of 150 work for the HK LiDAR data. Enough overlap between tiles so that there aren't edge effects
# #opt_chunk_size(las_cat) <- 150 # need to pick a value that is small enough for overlapping between tiles, but doesn't bust the compute power.
# opt_chunk_size(las_cat) <- 0 # need to pick a value that is small enough for overlapping between tiles, but doesn't bust the compute power.
# opt_chunk_buffer(las_cat) <- 15
# plot(las_cat, chunk_pattern = TRUE)
# 
# # parallel processing
# library(future)
# plan(multisession)
# 
# # normalise las files first for CHM
# # first, create DTM
# opt_output_files(las_cat) <- "D:/HK_LiDAR/2020/Normalised/{ORIGINALFILENAME}" # Martha's Maxtor hard drive
# opt_stop_early(las_cat) <- FALSE # bypasses errors
# las_cat_normal <- normalize_height(las_cat, res = 1, knnidw(k = 10, p = 2), keep_lowest = FALSE)

##### Normalised point cloud using DTM raster #####

#insert a loop which skips las files that throw an error
DTM_names <- list.files("D:/HK_LiDAR/2020/DTM/")
DTM_names <- str_replace(DTM_names, ".tif", "")
  
for(i in DTM_names){
  print(i)
  
  # if(file.exists(paste("D:/HK_LiDAR/2020/Normalised/", i, ".las", sep = ""))){
  # 
  #   }
  # 
  # else{
  
    if(file.exists(paste("Z:/Data/2020_LiDAR/", i, ".las", sep = ""))){
      
      if(file.exists(paste("D:/HK_LiDAR/2020/DTM/", i, ".tif", sep = ""))){
      
        las <- readLAS(paste("Z:/Data/2020_LiDAR/", i, ".las", sep = ""), select = "xyzc", filter = "-keep_class 2L 3L 4L 5L")
        
        dtm <- raster(paste("D:/HK_LiDAR/2020/DTM/", i, ".tif", sep = ""))
      
        las_normal <- las - dtm
      
        writeLAS(las_normal, paste("D:/HK_LiDAR/2020/Normalised/", i, ".las", sep = ""))
      }
   }
  }
# }
  
##### vertical metrics (below 0.2, 0.2-1, 1-5, 5-20, >20) #####

DTM_names <- list.files("D:/HK_LiDAR/2020/Normalised/")
DTM_names <- str_replace(DTM_names, ".las", "")

for(i in DTM_names){
  print(i)
  
  if(file.exists(paste("D:/HK_LiDAR/2020/Normalised/", i, ".las", sep = ""))){
      
    las <- readLAS(paste("D:/HK_LiDAR/2020/Normalised/", i, ".las", sep = ""), select = "xyzc") 

    las_veg <- filter_poi(las, Classification == 2L | Classification == 3L | Classification == 4L | Classification == 5L)
    
    myMetrics = function(z){
      
      npoints <- length(z)
      v_b02 <- length(z[z < 0.2])
      v_02_1 <- length(z[z > 0.2 & z < 1])
      v_1_5 <- length(z[z > 1 & z < 5])
      v_5_20 <- length(z[z > 5 & z < 20])
      v_a20 <- length(z[z > 20])
      
      below_02 <- v_b02/npoints*100
      betw_02_1 <- v_02_1/npoints*100
      betw_1_5 <- v_1_5/npoints*100
      betw_5_20 <- v_5_20/npoints*100
      above_20 <- v_a20/npoints*100
      
      # User's metrics
      metrics <- list(
        below_02 =  below_02,
        betw_02_1 = betw_02_1,
        betw_1_5 = betw_1_5,
        betw_5_20 = betw_5_20,
        above_20 = above_20
      )
      
      return(c(metrics))
    }
    
    m_vert <- pixel_metrics(las_veg, ~myMetrics(Z), res = 1)
    writeRaster(m_vert, paste("D:/HK_LiDAR/2020/Vertical_metrics/", i, ".tif", sep = ""))
#    writeRaster(m_vert, paste("/Users/marthaledger/OneDrive - The University Of Hong Kong/HKU/Butterflies/Data/HK_RS_data/2020/L2Processed_Martha/H4/", "v_metrics_", i, ".tif", sep = ""))
    
  }
}

# time to mosaic all vertical rasters...

path <- "D:/HK_LiDAR/2020/Vertical_metrics/"
list <- list.files(path, "*.tif", full=TRUE)

ListRasters <- function(list) {
  raster_list <- list() # initialise the list of rasters
  for (i in 1:(length(list))){ 
    grd_name <- list[i] # list_names contains all the names of the images in .grd format
    raster_file <- raster(grd_name, band = 3)
  }
  raster_list <- append(raster_list, raster_file) # update raster_list at each iteration
}

tic()
raster.list <- sapply(list, FUN = ListRasters)
toc()

names(raster.list) <- NULL

raster.list$fun <- mean
tic()
mos <- do.call(mosaic, raster.list)
toc()

plot(mos)

crs(mos) <- '+init=EPSG:2326'
terra::writeRaster(mos, "C:/Users/TCB-Martha/OneDrive - The University Of Hong Kong/HKU/Butterflies/Data/HK_RS_data/2020/L2Processed_Martha/H4/vertical_metrics_betw_15_HK_test.tif", overwrite = T)


##### all other rasters... #####
# height, total roughness, low roughness, elevation, slope, open area, open patches, edge extent, urban land cover %

library(raster)
library(terra)
library(tictoc)
library(dplyr)

##### canopy height model (vegetation and ground points only) #####

# multiples of 150 work for the HK LiDAR data. Enough overlap between tiles so that there aren't edge effects
#opt_chunk_size(las_cat) <- 150 # need to pick a value that is small enough for overlapping between tiles, but doesn't bust the compute power.
las_cat_normal <- readLAScatalog("D:/HK_LiDAR/2020/Normalised")
opt_chunk_size(las_cat_normal) <- 0 # need to pick a value that is small enough for overlapping between tiles, but doesn't bust the compute power.
opt_chunk_buffer(las_cat_normal) <- 15
plot(las_cat_normal, chunk_pattern = TRUE)

# parallel processing
library(future)
plan(multisession)

# normalise las files
# first, create DTM
opt_output_files(las_cat_normal) <- "D:/HK_LiDAR/2020/CHM/{ORIGINALFILENAME}" # Martha's Maxtor hard drive
opt_stop_early(las_cat_normal) <- FALSE # bypasses errors

#las_cat_CHM <- rasterize_canopy(las_cat_normal, res = 1, pitfree(thresholds = c(0, 10, 20), max_edge = c(0, 1.5))) - pitfree can't cope with vegetation points only
las_cat_CHM <- grid_metrics(las_cat_normal, ~quantile(Z, 0.90), res=1)

# mosaic CHM raster
path <- "D:/HK_LiDAR/2020/CHM/"
list <- list.files(path, "*.tif", full=TRUE)

ListRasters <- function(list) {
  raster_list <- list() # initialise the list of rasters
  for (i in 1:(length(list))){ 
    grd_name <- list[i] # list_names contains all the names of the images in .grd format
    raster_file <- raster(grd_name, band = 1)
  }
  raster_list <- append(raster_list, raster_file) # update raster_list at each iteration
}

tic()
raster.list <- sapply(list, FUN = ListRasters)
toc()

names(raster.list) <- NULL

raster.list$fun <- mean
tic()
mos <- do.call(mosaic, raster.list)
toc()

plot(mos)

crs(mos) <- '+init=EPSG:2326'
terra::writeRaster(mos, "C:/Users/TCB-Martha/OneDrive - The University Of Hong Kong/HKU/Butterflies/Data/HK_RS_data/2020/L2Processed_Martha/H4/2020_CHM.tif", overwrite = T)


##### total veg roughness #####

CHM_names <- list.files("D:/HK_LiDAR/2020/CHM/")

for(i in CHM_names){
  print(i)
  
  if(file.exists(paste("D:/HK_LiDAR/2020/CHM/", i, sep = ""))){
    
    chm <- raster(paste("D:/HK_LiDAR/2020/CHM/", i, sep = ""))
    
    tri_dsm <- terrain(chm, opt = "roughness", neighbors=8)
    
    terra::writeRaster(tri_dsm, paste("D:/HK_LiDAR/2020/Total_roughness/", i, sep = ""), overwrite = T)
  }
}


# mosaic CHM raster
path <- "D:/HK_LiDAR/2020/Total_roughness/"
list <- list.files(path, "*.tif", full=TRUE)

ListRasters <- function(list) {
  raster_list <- list() # initialise the list of rasters
  for (i in 1:(length(list))){ 
    grd_name <- list[i] # list_names contains all the names of the images in .grd format
    raster_file <- raster(grd_name, band = 1)
  }
  raster_list <- append(raster_list, raster_file) # update raster_list at each iteration
}

tic()
raster.list <- sapply(list, FUN = ListRasters)
toc()

names(raster.list) <- NULL

raster.list$fun <- mean
tic()
mos <- do.call(mosaic, raster.list)
toc()

plot(mos)

crs(mos) <- '+init=EPSG:2326'
terra::writeRaster(mos, "C:/Users/TCB-Martha/OneDrive - The University Of Hong Kong/HKU/Butterflies/Data/HK_RS_data/2020/L2Processed_Martha/H4/2020_total_roughness.tif", overwrite = T)


##### low veg roughness #####

# multiples of 150 work for the HK LiDAR data. Enough overlap between tiles so that there aren't edge effects
las_cat_normal <- readLAScatalog("D:/HK_LiDAR/2020/Normalised")
opt_chunk_size(las_cat_normal) <- 0 # need to pick a value that is small enough for overlapping between tiles, but doesn't bust the compute power.
opt_chunk_buffer(las_cat_normal) <- 15
plot(las_cat_normal, chunk_pattern = TRUE)

# parallel processing
library(future)
plan(multisession)

opt_output_files(las_cat_normal) <- "D:/HK_LiDAR/2020/DSM_B1/{ORIGINALFILENAME}" # Martha's Maxtor hard drive
opt_stop_early(las_cat_normal) <- FALSE # bypasses errors

dsm_b1 <- grid_metrics(las_cat_normal, ~mean(Z[Z<1]), res=1)

#

dsm_b1_names <- list.files("D:/HK_LiDAR/2020/DSM_B1/")

for(i in dsm_b1_names[2363:3268]){
  print(i)
  
  if(file.exists(paste("D:/HK_LiDAR/2020/DSM_B1/", i, sep = ""))){
    
    dsm_b1 <- raster(paste("D:/HK_LiDAR/2020/DSM_B1/", i, sep = ""))
    
    tri_dsm_b1 <- terrain(dsm_b1, opt = "roughness", neighbors=8)
    
    terra::writeRaster(tri_dsm_b1, paste("D:/HK_LiDAR/2020/Low_roughness/", i, sep = ""), overwrite = T)
  }
}


# mosaic CHM raster
path <- "D:/HK_LiDAR/2020/Low_roughness/"
list <- list.files(path, "*.tif", full=TRUE)

ListRasters <- function(list) {
  raster_list <- list() # initialise the list of rasters
  for (i in 1:(length(list))){ 
    grd_name <- list[i] # list_names contains all the names of the images in .grd format
    raster_file <- raster(grd_name, band = 1)
  }
  raster_list <- append(raster_list, raster_file) # update raster_list at each iteration
}

tic()
raster.list <- sapply(list, FUN = ListRasters)
toc()

names(raster.list) <- NULL

raster.list$fun <- mean
tic()
mos <- do.call(mosaic, raster.list)
toc()

plot(mos)

crs(mos) <- '+init=EPSG:2326'
terra::writeRaster(mos, "C:/Users/TCB-Martha/OneDrive - The University Of Hong Kong/HKU/Butterflies/Data/HK_RS_data/2020/L2Processed_Martha/H4/2020_low_roughness.tif", overwrite = T)


##### elevation #####

path <- "D:/HK_LiDAR/2020/DTM/"
list <- list.files(path, "*.tif", full=TRUE)

ListRasters <- function(list) {
  raster_list <- list() # initialise the list of rasters
  for (i in 1:(length(list))){ 
    grd_name <- list[i] # list_names contains all the names of the images in .grd format
    raster_file <- raster(grd_name)
  }
  raster_list <- append(raster_list, raster_file) # update raster_list at each iteration
}

tic()
raster.list <- sapply(list, FUN = ListRasters)
toc()

names(raster.list) <- NULL

raster.list$fun <- mean
tic()
mos <- do.call(mosaic, raster.list)
toc()

plot(mos)

crs(mos) <- '+init=EPSG:2326'
terra::writeRaster(mos, "C:/Users/TCB-Martha/OneDrive - The University Of Hong Kong/HKU/Butterflies/Data/HK_RS_data/2020/L2Processed_Martha/H4/2020_DTM.tif", overwrite = T)

 
##### slope #####

#7NE15A(e843n833,e844n833)

DTM_names <- list.files("D:/HK_LiDAR/2020/DTM/")

for(i in DTM_names[2349:3253]){
  print(i)
      
      if(file.exists(paste("D:/HK_LiDAR/2020/DTM/", i, sep = ""))){
        
        dtm <- raster(paste("D:/HK_LiDAR/2020/DTM/", i, sep = ""))
        
        dtm_slope <- terrain(dtm, v = "slope", unit = "degrees", neighbors = 4)
        
        terra::writeRaster(dtm_slope, paste("D:/HK_LiDAR/2020/Slope/", i, sep = ""), overwrite = T)
    }
  }
# }  

# mosaic slope
path <- "D:/HK_LiDAR/2020/Slope/"
list <- list.files(path, "*.tif", full=TRUE)

ListRasters <- function(list) {
  raster_list <- list() # initialise the list of rasters
  for (i in 1:(length(list))){ 
    grd_name <- list[i] # list_names contains all the names of the images in .grd format
    raster_file <- raster(grd_name)
  }
  raster_list <- append(raster_list, raster_file) # update raster_list at each iteration
}

tic()
raster.list <- sapply(list, FUN = ListRasters)
toc()

names(raster.list) <- NULL

raster.list$fun <- mean
tic()
mos <- do.call(mosaic, raster.list)
toc()

plot(mos)

crs(mos) <- '+init=EPSG:2326'
terra::writeRaster(mos, "C:/Users/TCB-Martha/OneDrive - The University Of Hong Kong/HKU/Butterflies/Data/HK_RS_data/2020/L2Processed_Martha/H4/2020_slope.tif", overwrite = T)


##### urbanicity #####
# https://stackoverflow.com/questions/65796267/extract-mean-value-of-raster-with-buffer-condition-on-second-layer-attribute 

LUMHK <- raster("C:/Users/TCB-Martha/OneDrive - The University Of Hong Kong/D24H/Rat_project/Data/LUMHK_RasterGrid_2021/LUMHK_RasterGrid_2021/LUM_end2021.tif")
dtm <- rast("C:/Users/TCB-Martha/OneDrive - The University Of Hong Kong/HKU/Butterflies/Data/HK_RS_data/2020/L2Processed_Martha/H4/2020_DTM.tif")

# reclassify the values into two groups (urban, not urban)
m <- c(1, 1,
       2, 1,
       3, 1,
       11, 1,
       21, 1,
       22, 1,
       23, 1,
       31, 1,
       32, 1,
       41, 1,
       42, 1,
       43, 1,
       44, 1,
       51, 1,
       52, 1,
       53, 1,
       54, 1,
       61, 0,
       62, 0,
       71, 0,
       72, 0,
       73, 0,
       74, 0,
       81, 0,
       83, 0,
       91, 0,
       92, 0)
rclmat <- matrix(m, ncol=2, byrow=TRUE)

new_LUMHK <- reclassify(LUMHK, rclmat) # the binary urban raster (urban = 1, non-urban = 0)
new_LUMHK <- rast(new_LUMHK)

r <- dtm > -Inf
p <- as(terra::as.polygons(r, dissolve = T, values = F), "Spatial") #To create a mask
m <- raster::crop(raster::mask(raster(new_LUMHK), p), p) #To mask the DTM

# try focal statistics (terra::focal) - Adam to the rescue!

getCircleMatrix <- function(d){  # d stands for dimension
  
  if(!as.logical(d%%2)){
    # calculate centre point
    centre_pt_coordinates <- d/2 +0.5
    
    # create coordinate matrix
    coor_matrix <- matrix(c(rep(c(1:d,0), d), rep(0,d+1)),nrow = d+1, byrow = TRUE)
    
    transposed_coor_matrix <- t(coor_matrix)
    
    valid_matrix <- (transposed_coor_matrix > 0) *1
    
    vert_distance <- (transposed_coor_matrix- centre_pt_coordinates)**2 * valid_matrix
    hort_distance <- (coor_matrix- centre_pt_coordinates)**2*valid_matrix
    
    distance_matrix <- vert_distance+hort_distance
    
    circle_radius <- (centre_pt_coordinates)**2 # minus 1 is optional - if you remove it, the circle will fit right to the edge.
    
    circle_matrix <- (distance_matrix <= circle_radius) *1 *valid_matrix
    
    return(circle_matrix)
    
  }
  
  if(as.logical(d%%2)){
    # calculate centre point
    centre_pt_coordinates <- (d+1)/2
    
    # create matrix
    coor_matrix <- matrix(rep(1:d, d),nrow = d, byrow = TRUE)
    transposed_coor_matrix <- t(coor_matrix)
    
    vert_distance <- (transposed_coor_matrix- centre_pt_coordinates)**2 
    hort_distance <- (coor_matrix- centre_pt_coordinates)**2
    
    distance_matrix <- vert_distance+hort_distance
    
    circle_radius <- (centre_pt_coordinates-1)**2
    
    circle_matrix <- (distance_matrix <= circle_radius) *1
    
    return(circle_matrix)
  }
  
}
circle_matrix <- getCircleMatrix(10)

plan(multisession)
final_raster <- terra::focal(m, w=circle_matrix, fun="mean", na.policy="omit")

crs(final_raster) <- '+init=EPSG:2326'
writeRaster(final_raster, "C:/Users/TCB-Martha/OneDrive - The University Of Hong Kong/HKU/Butterflies/Data/HK_RS_data/2020/L2Processed_Martha/H4/2020_urbanicity.tif", overwrite = T)


##### open area #####
# how to use focal statistics with this?
# function within terra::focal - If 1, sum. If 0, 0.
library(landscapemetrics)

chm <- raster("C:/Users/TCB-Martha/OneDrive - The University Of Hong Kong/HKU/Butterflies/Data/HK_RS_data/2020/L2Processed_Martha/H4/2020_CHM.tif")

height_class <- reclassify(chm, c(-Inf,1,1, # this raster is enough for creating plot averages, but not for predictive modelling using rasters...
                                  1,5,0,
                                  5,Inf,0))

# write the raster for model training
crs(height_class) <- '+init=EPSG:2326'
terra::writeRaster(height_class, "C:/Users/TCB-Martha/OneDrive - The University Of Hong Kong/HKU/Butterflies/Data/HK_RS_data/2020/L2Processed_Martha/H4/2020_openarea.tif", overwrite = T)

# create circular moving window

getCircleMatrix <- function(d){  # d stands for dimension
  
  if(!as.logical(d%%2)){
    # calculate centre point
    centre_pt_coordinates <- d/2 +0.5
    
    # create coordinate matrix
    coor_matrix <- matrix(c(rep(c(1:d,0), d), rep(0,d+1)),nrow = d+1, byrow = TRUE)
    
    transposed_coor_matrix <- t(coor_matrix)
    
    valid_matrix <- (transposed_coor_matrix > 0) *1
    
    vert_distance <- (transposed_coor_matrix- centre_pt_coordinates)**2 * valid_matrix
    hort_distance <- (coor_matrix- centre_pt_coordinates)**2*valid_matrix
    
    distance_matrix <- vert_distance+hort_distance
    
    circle_radius <- (centre_pt_coordinates)**2 # minus 1 is optional - if you remove it, the circle will fit right to the edge.
    
    circle_matrix <- (distance_matrix <= circle_radius) *1 *valid_matrix
    
    return(circle_matrix)
    
  }
  
  if(as.logical(d%%2)){
    # calculate centre point
    centre_pt_coordinates <- (d+1)/2
    
    # create matrix
    coor_matrix <- matrix(rep(1:d, d),nrow = d, byrow = TRUE)
    transposed_coor_matrix <- t(coor_matrix)
    
    vert_distance <- (transposed_coor_matrix- centre_pt_coordinates)**2 
    hort_distance <- (coor_matrix- centre_pt_coordinates)**2
    
    distance_matrix <- vert_distance+hort_distance
    
    circle_radius <- (centre_pt_coordinates-1)**2
    
    circle_matrix <- (distance_matrix <= circle_radius) *1
    
    return(circle_matrix)
  }
  
}
circle_matrix <- getCircleMatrix(200)


library(future)
plan(multisession)
final_raster <- terra::focal(rast(height_class), w=circle_matrix, fun="sum", na.rm = T) 
final_raster_ha <- final_raster*0.0001 # convert area to hectares.

crs(final_raster) <- '+init=EPSG:2326'
crs(final_raster_ha) <- '+init=EPSG:2326'
terra::writeRaster(final_raster, "C:/Users/TCB-Martha/OneDrive - The University Of Hong Kong/HKU/Butterflies/Data/HK_RS_data/2020/L2Processed_Martha/H4/2020_openarea_predict.tif", overwrite = T)
terra::writeRaster(final_raster_ha, "C:/Users/TCB-Martha/OneDrive - The University Of Hong Kong/HKU/Butterflies/Data/HK_RS_data/2020/L2Processed_Martha/H4/2020_openarea_predict_ha.tif", overwrite = T)


##### open patch number #####

library(landscapemetrics)

chm <- raster("C:/Users/TCB-Martha/OneDrive - The University Of Hong Kong/HKU/Butterflies/Data/HK_RS_data/2020/L2Processed_Martha/H4/2020_CHM.tif")
chm <- raster("/Users/marthaledger/OneDrive - The University Of Hong Kong/HKU/Butterflies/Data/HK_RS_data/2020/L2Processed_Martha/H4/2020_CHM.tif")

height_class <- reclassify(chm, c(-Inf,1,1, # this raster is enough for creating plot averages, but not for predictive modelling using rasters...
                                  1,5,0,
                                  5,Inf,0))


# resample the open area raster to 2m (R can't handle 1m resolution for this function)
height_class_5 <- aggregate(height_class, fact=5, fun=modal)

#height_class <- rast(height_class)

# write the raster for model training
# crs(height_class) <- '+init=EPSG:2326'
# terra::writeRaster(height_class, "C:/Users/TCB-Martha/OneDrive - The University Of Hong Kong/HKU/Butterflies/Data/HK_RS_data/2020/L2Processed_Martha/H4/2020_openarea.tif", overwrite = T)

# create circular moving window

getCircleMatrix <- function(d){  # d stands for dimension
  
  if(!as.logical(d%%2)){
    # calculate centre point
    centre_pt_coordinates <- d/2 +0.5
    
    # create coordinate matrix
    coor_matrix <- matrix(c(rep(c(1:d,0), d), rep(0,d+1)),nrow = d+1, byrow = TRUE)
    
    transposed_coor_matrix <- t(coor_matrix)
    
    valid_matrix <- (transposed_coor_matrix > 0) *1
    
    vert_distance <- (transposed_coor_matrix- centre_pt_coordinates)**2 * valid_matrix
    hort_distance <- (coor_matrix- centre_pt_coordinates)**2*valid_matrix
    
    distance_matrix <- vert_distance+hort_distance
    
    circle_radius <- (centre_pt_coordinates)**2 # minus 1 is optional - if you remove it, the circle will fit right to the edge.
    
    circle_matrix <- (distance_matrix <= circle_radius) *1 *valid_matrix
    
    return(circle_matrix)
    
  }
  
  if(as.logical(d%%2)){
    # calculate centre point
    centre_pt_coordinates <- (d+1)/2
    
    # create matrix
    coor_matrix <- matrix(rep(1:d, d),nrow = d, byrow = TRUE)
    transposed_coor_matrix <- t(coor_matrix)
    
    vert_distance <- (transposed_coor_matrix- centre_pt_coordinates)**2 
    hort_distance <- (coor_matrix- centre_pt_coordinates)**2
    
    distance_matrix <- vert_distance+hort_distance
    
    circle_radius <- (centre_pt_coordinates-1)**2
    
    circle_matrix <- (distance_matrix <= circle_radius) *1
    
    return(circle_matrix)
  }
  
}
circle_matrix <- getCircleMatrix(40)

# function for open patches
open_patches <- function(x){
  y <- lsm_c_np(x, directions = 8)
  return(y$value[2])
}

library(future)
plan(multisession)
final_raster <- terra::focal(rast(height_class_5), w=circle_matrix, fun=open_patches) 

crs(final_raster) <- '+init=EPSG:2326'
terra::writeRaster(final_raster, "C:/Users/TCB-Martha/OneDrive - The University Of Hong Kong/HKU/Butterflies/Data/HK_RS_data/2020/L2Processed_Martha/H4/2020_openpatches_predict.tif", overwrite = T)
terra::writeRaster(final_raster, "/Users/marthaledger/OneDrive - The University Of Hong Kong/HKU/Butterflies/Data/HK_RS_data/2020/L2Processed_Martha/H4/2020_openpatches_predict.tif", overwrite = T)


##### edge extent (need to make a custom function somehow. Landscape metrics does not work within terra::focal) #####

library(landscapemetrics)

chm <- raster("C:/Users/TCB-Martha/OneDrive - The University Of Hong Kong/HKU/Butterflies/Data/HK_RS_data/2020/L2Processed_Martha/H4/2020_CHM.tif")
chm <- raster("/Users/marthaledger/OneDrive - The University Of Hong Kong/HKU/Butterflies/Data/HK_RS_data/2020/L2Processed_Martha/H4/2020_CHM.tif")

height_class <- reclassify(chm, c(-Inf,1,1, # this raster is enough for creating plot averages, but not for predictive modelling using rasters...
                                  1,5,0,
                                  5,Inf,0))


# resample the open area raster to 2m (R can't handle 1m resolution for this function)
height_class_5 <- aggregate(height_class, fact=5, fun=modal)

#height_class <- rast(height_class)

# write the raster for model training
# crs(height_class) <- '+init=EPSG:2326'
# terra::writeRaster(height_class, "C:/Users/TCB-Martha/OneDrive - The University Of Hong Kong/HKU/Butterflies/Data/HK_RS_data/2020/L2Processed_Martha/H4/2020_openarea.tif", overwrite = T)

# create circular moving window

getCircleMatrix <- function(d){  # d stands for dimension
  
  if(!as.logical(d%%2)){
    # calculate centre point
    centre_pt_coordinates <- d/2 +0.5
    
    # create coordinate matrix
    coor_matrix <- matrix(c(rep(c(1:d,0), d), rep(0,d+1)),nrow = d+1, byrow = TRUE)
    
    transposed_coor_matrix <- t(coor_matrix)
    
    valid_matrix <- (transposed_coor_matrix > 0) *1
    
    vert_distance <- (transposed_coor_matrix- centre_pt_coordinates)**2 * valid_matrix
    hort_distance <- (coor_matrix- centre_pt_coordinates)**2*valid_matrix
    
    distance_matrix <- vert_distance+hort_distance
    
    circle_radius <- (centre_pt_coordinates)**2 # minus 1 is optional - if you remove it, the circle will fit right to the edge.
    
    circle_matrix <- (distance_matrix <= circle_radius) *1 *valid_matrix
    
    return(circle_matrix)
    
  }
  
  if(as.logical(d%%2)){
    # calculate centre point
    centre_pt_coordinates <- (d+1)/2
    
    # create matrix
    coor_matrix <- matrix(rep(1:d, d),nrow = d, byrow = TRUE)
    transposed_coor_matrix <- t(coor_matrix)
    
    vert_distance <- (transposed_coor_matrix- centre_pt_coordinates)**2 
    hort_distance <- (coor_matrix- centre_pt_coordinates)**2
    
    distance_matrix <- vert_distance+hort_distance
    
    circle_radius <- (centre_pt_coordinates-1)**2
    
    circle_matrix <- (distance_matrix <= circle_radius) *1
    
    return(circle_matrix)
  }
  
}
circle_matrix <- getCircleMatrix(40)

# function for edge extent
edge_extent <- function(x){
  y <- lsm_c_te(x, directions = 8)
  return(y$value[2])
}

library(future)
plan(multisession)
final_raster <- terra::focal(rast(height_class_5), w=circle_matrix, fun=edge_extent) 

crs(final_raster) <- '+init=EPSG:2326'
terra::writeRaster(final_raster, "C:/Users/TCB-Martha/OneDrive - The University Of Hong Kong/HKU/Butterflies/Data/HK_RS_data/2020/L2Processed_Martha/H4/2020_edgeextent_predict.tif", overwrite = T)
terra::writeRaster(final_raster, "/Users/marthaledger/OneDrive - The University Of Hong Kong/HKU/Butterflies/Data/HK_RS_data/2020/L2Processed_Martha/H4/2020_edgeextent_predict.tif", overwrite = T)

##### 2010 LiDAR metrics #####

# use the hashed below to chunk tiles into overlapping segments and produce normalised las files

las_cat_2010 <- readLAScatalog("Z:/Data/2010_LiDAR/Input")
HK1980 <- st_crs(2326) # set projection ('2326' is the HK1980 EPSG)
st_crs(las_cat_2010) <- HK1980 # apply projection to las file
summary(las_cat_2010)
las_check(las_cat_2010)
plot(las_cat_2010)

# multiples of 150 work for the HK LiDAR data. Enough overlap between tiles so that there aren't edge effects
#opt_chunk_size(las_cat) <- 150 # need to pick a value that is small enough for overlapping between tiles, but doesn't bust the compute power.
opt_chunk_size(las_cat_2010) <- 0 # need to pick a value that is small enough for overlapping between tiles, but doesn't bust the compute power.
opt_chunk_buffer(las_cat_2010) <- 15
plot(las_cat_2010, chunk_pattern = TRUE)

# parallel processing
library(future)
plan(multisession)

# normalise las files
# first, create DTM
opt_output_files(las_cat_2010) <- "D:/HK_LiDAR/2010/DTM/{ORIGINALFILENAME}" # Martha's Maxtor hard drive
#opt_output_files(las_cat) <- "/Users/marthaledger/OneDrive - The University Of Hong Kong/HKU/Butterflies/Data/HK_RS_data/2020/L2Processed_Martha/H4/Normal/norm_{XLEFT}_{YBOTTOM}"

opt_stop_early(las_cat_2010) <- FALSE # bypasses errors
#las_cat_normal <- normalize_height(las_cat, res = 1, knnidw(k = 10, p = 2), keep_lowest = FALSE)
#las_cat_normal <- las_cat - rasterize_terrain(las_cat, res = 1, knnidw())
las_cat_DTM <- rasterize_terrain(las_cat_2010, res = 1, knnidw())
# can opt to restart at the chunk which failed: requires being present at script though.
# opt_restart(ctg) <- 9

path = "D:/HK_LiDAR/2010/DTM/"

tiffcat <- function(f){
  do.call(rbind,
          lapply(f, function(r){
            ras = rast(r)
            poly = as.polygons(ext(ras), crs=crs(ras))
            poly$path = r
            poly
          })
  )
}  

files = list.files(path, "*.tif", full=TRUE)
tc = tiffcat(files)  

library(stringr)
DTM_names <- list.files("D:/HK_LiDAR/2010/DTM/")
DTM_names <- str_replace(DTM_names, ".tif", "")
# remove .tif suffix

# insert a loop which skips las files that throw an error
for(i in DTM_names){
  print(i)
  
  if(file.exists(paste("Z:/Data/2010_LiDAR/Input/", i, ".las", sep = ""))){
    
    if(file.exists(paste("D:/HK_LiDAR/2010/DTM/", i, ".tif", sep = ""))){
      
      las <- readLAS(paste("Z:/Data/2010_LiDAR/Input/", i, ".las", sep = ""), select = "xyzc")
      
      dtm <- raster(paste("D:/HK_LiDAR/2010/DTM/", i, ".tif", sep = ""))
      
      las_normal <- las - dtm
      
      writeLAS(las_normal, paste("D:/HK_LiDAR/2010/Normalised/", i, ".las", sep = ""))
    }
  }
}

# now produce metrics!

DTM_names <- list.files("D:/HK_LiDAR/2010/Normalised/")
DTM_names <- str_replace(DTM_names, ".las", "")

for(i in DTM_names){
  print(i)
  
  if(file.exists(paste("D:/HK_LiDAR/2010/Normalised/", i, ".las", sep = ""))){
    
    las <- readLAS(paste("D:/HK_LiDAR/2010/Normalised/", i, ".las", sep = ""), select = "xyzc") 
    
    las_veg <- filter_poi(las, Classification == 2L | Classification == 3L | Classification == 4L | Classification == 5L)
    
    myMetrics = function(z){
      
      npoints <- length(z)
      v_b02 <- length(z[z < 0.2])
      v_02_1 <- length(z[z > 0.2 & z < 1])
      v_1_5 <- length(z[z > 1 & z < 5])
      v_5_20 <- length(z[z > 5 & z < 20])
      v_a20 <- length(z[z > 20])
      
      below_02 <- v_b02/npoints*100
      betw_02_1 <- v_02_1/npoints*100
      betw_1_5 <- v_1_5/npoints*100
      betw_5_20 <- v_5_20/npoints*100
      above_20 <- v_a20/npoints*100
      
      # User's metrics
      metrics <- list(
        below_02 =  below_02,
        betw_02_1 = betw_02_1,
        betw_1_5 = betw_1_5,
        betw_5_20 = betw_5_20,
        above_20 = above_20
      )
      
      return(c(metrics))
    }
    
    m_vert <- pixel_metrics(las_veg, ~myMetrics(Z), res = 1)
    writeRaster(m_vert, paste("D:/HK_LiDAR/2010/Vertical_metrics/", i, ".tif", sep = ""))
    #    writeRaster(m_vert, paste("/Users/marthaledger/OneDrive - The University Of Hong Kong/HKU/Butterflies/Data/HK_RS_data/2020/L2Processed_Martha/H4/", "v_metrics_", i, ".tif", sep = ""))
    
  }
}

# time to mosaic all rasters...!

library(raster)
library(terra)
library(tictoc)
library(dplyr)
library(tictoc)

path <- "D:/HK_LiDAR/2010/Vertical_metrics/"
list <- list.files(path, "*.tif", full=TRUE)

ListRasters <- function(list) {
  raster_list <- list() # initialise the list of rasters
  for (i in 1:(length(list))){ 
    grd_name <- list[i] # list_names contains all the names of the images in .grd format
    raster_file <- raster(grd_name, band = 5)
  }
  raster_list <- append(raster_list, raster_file) # update raster_list at each iteration
}

tic()
raster.list <- sapply(list, FUN = ListRasters)
toc()

names(raster.list) <- NULL

raster.list$fun <- mean
tic()
mos <- do.call(mosaic, raster.list)
toc()

plot(mos)

crs(mos) <- '+init=EPSG:2326'
terra::writeRaster(mos, "C:/Users/TCB-Martha/OneDrive - The University Of Hong Kong/HKU/Butterflies/Data/HK_RS_data/2010/vertical_metrics_2010_above_20_HK.tif", overwrite = T)
