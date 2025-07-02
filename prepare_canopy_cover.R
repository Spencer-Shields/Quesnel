library(tidyverse)
library(terra)
library(lidR)
library(pbapply)
library(future)
library(future.apply)
# library(lasR)
library(Rcpp)
source('helper_functions.R')

parent_dir = 'data/quesnel_thinning_las'

pre_dir = paste0(parent_dir,'/pre')
pre_blockdirs = list.dirs(pre_dir, recursive=F)

post_dir = paste0(parent_dir,'/post')
post_blockdirs = list.dirs(post_dir, recursive=F)

#----inspect point clouds----

pre_harvest_las = pblapply(pre_blockdirs, function(x)readLAScatalog(file.path(x, "input/las/norm")))
names(pre_harvest_las) = basename(pre_blockdirs)

post_harvest_las = pblapply(post_blockdirs, function(x)readLAScatalog(file.path(x, "input/las/norm")))
names(post_harvest_las) = basename(post_blockdirs)

#----view point density of pre-harvest and post harvest point clouds----
# Set resolution
res <- 0.25

# Function to compute multiband density rasters for one LAScatalog
compute_density_pixel_metrics <- function(lascat, resolution = res) {
  # This calculates: total points, and counts per return number (1st, 2nd, 3rd, 4th)
  metric_fun <- ~ list(
    total = length(Z),
    return_1 = sum(ReturnNumber == 1),
    return_2 = sum(ReturnNumber == 2),
    return_3 = sum(ReturnNumber == 3),
    return_4 = sum(ReturnNumber == 4)
  )
  
  # Compute metrics
  r <- pixel_metrics(lascat, func = metric_fun, res = resolution)
  return(r)  # convert to terra raster
}

plan('multisession', workers=8)
pre_density_rasters <- pblapply(pre_harvest_las, compute_density_pixel_metrics)
names(pre_density_rasters) <- names(pre_harvest_las)

post_density_rasters <- pblapply(post_harvest_las, compute_density_pixel_metrics)
names(post_density_rasters) <- names(post_harvest_las)

density_rast_dir = file.path(parent_dir, 'point_density_rasters')
dir.check(density_rast_dir)
pre_density_dir = file.path(density_rast_dir, 'pre')
dir.check(pre_density_dir)
post_density_dir = file.path(density_rast_dir, 'post')
dir.check(post_density_dir)

pblapply(1:length(pre_density_rasters), function(i){
  block_id = names(pre_density_rasters)[i]
  pre_filename = paste0(pre_density_dir,'/',block_id,'.tif')
  if(!file.exists(pre_filename)){writeRaster(pre_density_rasters[[block_id]], pre_filename)}
  post_filename = paste0(post_density_dir,'/',block_id,'.tif')
  if(!file.exists(post_filename)){writeRaster(post_density_rasters[[block_id]], post_filename)}
})
plan('sequential')

#----decimate post-harvest lidar----
# target_density = 5
# res_value = 0.25
# 
# plan('multisession', workers=4)
# post_harvest_decimated = pblapply(post_blockdirs, function(x){ #x is the name of the thinning block
#   dir <- file.path(x, "input/las/norm")
#   ctg <- readLAScatalog(dir)
#   
#   # Set processing options
#   opt_chunk_size(ctg) <- 250
#   opt_chunk_buffer(ctg) <- res_value  # Set explicitly instead of relying on algorithm env
#   opt_progress(ctg) <- TRUE
#   opt_output_files(ctg) <- tempfile(pattern = "decimated_", fileext = "_{ID}")
#   
#   #decimate pointcloud
#   # d_lidr = decimate_points(ctg, random_per_voxel(res=res_value, n=target_density))
#   
#   # Create algorithm and manually assign res into its environment
#   algo <- random_per_voxel(res = res_value, n = target_density)
#   environment(algo)$res <- res_value  # <<-- This is the critical line!
#   
#   # Now call decimate_points
#   decimated <- decimate_points(ctg, algo)
#   
#   puls = retrieve_pulses(ctg)
#   decimated2 = decimate_points(ctg, homogenize(density = 30, res = 1))
#   
#   # Canopy cover: proportion of first returns above 2m
#   # Return a raster with values 0–1 (fraction)
#   # canopy_cover <- grid_metrics(ctg, ~sum(ReturnNumber == 1 & Z > 2) / sum(ReturnNumber == 1), res = 0.25)
#   
#   return(decimated)
# }
# # ,cl = 'future'
# # ,future.seed=TRUE
# )
# names(post_harvest_decimated) = basename(pre_blockdirs)

pre_harvest_pulses_df = tribble(
  ~block_id, ~pulses_per_m2, ~points_per_m2,
  '12L_B8C3', 11.2, 25,
  '12L_C4', 21.9, 47.7,
  '12L_C5', 11.5, 24.9,
  '12L_C7', 11.5, 24.9,
  '12L_D345', 11.2, 24,
  '12N_1X1W', 11.4, 24.6,
  '12N_T3', 11.3, 25.5,
  '12N_V1', 11.5, 24.9
)

####random decimation
post_decimated_dir = paste0(parent_dir,'/post_decimated_random')
dir.check(post_decimated_dir)

# plan('multisession', workers=4)
pblapply(post_blockdirs, function(x){
  
  #make output directory
  block_id = basename(x)
  dec_block_dir = paste0(post_decimated_dir,'/',block_id)
  dir.check(dec_block_dir)
  
  #get target pulse density for block
  block_pulses = pre_harvest_pulses_df$pulses_per_m2[pre_harvest_pulses_df$block_id==block_id]
  
  #get list of las files in input directory
  las_list = list.files(file.path(x, "input/las/norm")
                        , pattern = '\\.(laz|las)$', ignore.case = TRUE, recursive = T, full.names = T)
  #decimate every laz in the list, save to output directory
  pblapply(las_list, function(y){
    filename = paste0(dec_block_dir,'/',basename(y))
    if(!file.exists(filename)){
      las = readLAS(y)
      las_puls = retrieve_pulses(las)
      decimated = decimate_points(las_puls, random(density = block_pulses, use_pulse = T))
      writeLAS(decimated, filename, index=T)
      gc()
    }
  }
  # ,cl='future'
  )
})
# plan(future::sequential())

###homogenized decimation
post_decimated_dir = paste0(parent_dir,'/post_decimated_homogenize')
dir.check(post_decimated_dir)

# plan('multisession', workers=4)
pblapply(post_blockdirs, function(x){
  
  #make output directory
  block_id = basename(x)
  dec_block_dir = paste0(post_decimated_dir,'/',block_id)
  dir.check(dec_block_dir)
  
  #get target pulse density for block
  block_pulses = pre_harvest_pulses_df$pulses_per_m2[pre_harvest_pulses_df$block_id==block_id]
  
  #get list of las files in input directory
  las_list = list.files(file.path(x, "input/las/norm")
                        , pattern = '\\.(laz|las)$', ignore.case = TRUE, recursive = T, full.names = T)
  #decimate every laz in the list, save to output directory
  pblapply(las_list, function(y){
    filename = paste0(dec_block_dir,'/',basename(y))
    if(!file.exists(filename)){
      las = readLAS(y)
      las_puls = retrieve_pulses(las)
      decimated = decimate_points(las_puls, homogenize(density = block_pulses, res=5, use_pulse = T))
      writeLAS(decimated, filename, index=T)
      # gc()
    }
  }
  # ,cl='future'
  )
})
# plan(future::sequential())

###voxel decimation


#----produce canopy cover rasters----

pre_CC_l = pblapply(pre_blockdirs, function(x){
  dir <- file.path(x, "input/las/norm")
  ctg <- readLAScatalog(dir)
  
  # Set processing options
  opt_chunk_size(ctg) <- 250
  opt_chunk_buffer(ctg) <- 25
  # opt_cores(ctg) <- 4
  opt_progress(ctg) <- TRUE
  
  # Canopy cover: proportion of first returns above 2m
  # Return a raster with values 0–1 (fraction)
  canopy_cover <- grid_metrics(ctg, ~sum(ReturnNumber == 1 & Z > 2) / sum(ReturnNumber == 1), res = 0.25)
  
  return(canopy_cover)
}
# ,cl = 'future'
# ,future.seed=TRUE
)
names(pre_CC_l) = basename(pre_blockdirs)

# post_CC_l = pblapply(post_blockdirs, function(x){
#   dir <- file.path(x, "input/las/norm")
#   ctg <- readLAScatalog(dir)
#   
#   # Set processing options
#   opt_chunk_size(ctg) <- 200
#   opt_chunk_buffer(ctg) <- 25
#   # opt_cores(ctg) <- 4
#   opt_progress(ctg) <- TRUE
#   
#   # Canopy cover: proportion of first returns above 2m
#   # Return a raster with values 0–1 (fraction)
#   canopy_cover <- grid_metrics(ctg, ~sum(ReturnNumber == 1 & Z > 2) / sum(ReturnNumber == 1), res = 0.25)
#   
#   return(canopy_cover)
# }
# # ,cl='future'
# # ,future.seed=TRUE
# )
# names(post_CC_l) = basename(post_blockdirs) 


#get canopy cover of decimated point clouds
post_dec_blockdirs = list.dirs(post_decimated_dir)[list.dirs(post_decimated_dir)!=post_decimated_dir]

post_dec_CC_l = pblapply(post_dec_blockdirs, function(x){
  ctg = readLAScatalog(x)
  # Set processing options
  opt_chunk_size(ctg) <- 250
  opt_chunk_buffer(ctg) <- 25
  # opt_cores(ctg) <- 4
  opt_progress(ctg) <- TRUE
  
  # Canopy cover: proportion of first returns above 2m
  # Return a raster with values 0–1 (fraction)
  canopy_cover <- grid_metrics(ctg, ~sum(ReturnNumber == 1 & Z > 2) / sum(ReturnNumber == 1), res = 0.25)
  
  return(canopy_cover)
})
names(post_dec_CC_l) = basename(post_blockdirs)



# # plan('sequential')
# change_cc_l = pblapply(names(pre_CC_l), function(x){
#   dcc = post_CC_l[[x]] - pre_CC_l[[x]]
# })
# names(change_cc_l) = names(pre_CC_l)
# 
# plot(change_cc_l$`12N_T3`)
# 
# cc_dir = paste0(parent_dir,'/canopy_cover_change')
# dir.check(cc_dir)
# pblapply(1:length(change_cc_l), function(i){
#   id = names(change_cc_l)[i]
#   filename = paste0(cc_dir,'/',id,'.tif')
#   if(!file.exists(filename)){writeRaster(change_cc_l[[i]], filename)}
# })


#calculate and save canopy cover change rasters using decimated point clouds
change_dec_cc_l = pblapply(names(pre_CC_l), function(x){
  dcc = post_dec_CC_l[[x]] - pre_CC_l[[x]]
})
names(change_dec_cc_l) = names(pre_CC_l)

plot(change_dec_cc_l$`12N_T3`)

cc_dec_dir = paste0(parent_dir,'/canopy_cover_change_Decimated')
dir.check(cc_dec_dir)
pblapply(1:length(change_dec_cc_l), function(i){
  id = names(change_dec_cc_l)[i]
  filename = paste0(cc_dec_dir,'/',id,'.tif')
  if(!file.exists(filename)){writeRaster(change_dec_cc_l[[i]], filename)}
})

