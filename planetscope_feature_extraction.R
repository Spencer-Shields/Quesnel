library(tidyverse)
library(terra)
source('helper_functions.R')
source('https://raw.githubusercontent.com/Spencer-Shields/planetscope_radiometric_correction/refs/heads/main/softmax.R')
library(future)
library(future.apply)
library(glcm)
library(pbapply)
library(tools)
library(jsonify)
# library(microbenchmark)
library(RStoolbox)

#----load data----

#set directories
list.dirs('data')
ps_dir = 'data/planet_scenes'
raw_dir = paste0(ps_dir,'/raw')

#get dataframe of scene metadata
meta_df = ps_meta(ps_dir)

#filter dataframe to only include scenes which are more than 95% clear

clear_threshold = 95

clear_df = meta_df %>% 
  filter(clear_percent >= clear_threshold) %>% #get data where the percent of clear pixels is greater than or equal to a threshold
  filter(quality_category == 'standard') #get standard quality data

#see temporal distribution of PlanetScope scenes
clear_df[['acquisition_date']] = sapply(clear_df$acquired, function(x) (str_split_1(x, 'T')[1]))

clear_df = clear_df %>%
  mutate(acquisition_year = year(acquisition_date),
         acquisition_month = month(acquisition_date),
         acquisition_day = day(acquisition_date)) %>%
  mutate(acquisition_month_day = paste0(acquisition_month,'-',acquisition_day))

clear_df_dates_bymonth = clear_df %>%
  group_by(acquisition_year, acquisition_month) %>%
  summarise(Number_of_scenes = n())

ggplot()+
  geom_col(data = clear_df_dates_bymonth, aes(x = acquisition_month, y = Number_of_scenes))+
  facet_wrap(vars(acquisition_year), ncol = 1)


#----apply UDM2 mask to scenes that need it----

dir_string = 'UDMmasked'
udm_masked_dir = paste0(ps_dir, '/',dir_string)
dir.check(udm_masked_dir)

indir = ps_dir
outdir = udm_masked_dir

rasters = list.files(indir, full.names = T, recursive = T, pattern = '\\.tif$')

plan(multisession, workers = 8)


pblapply(1:length(clear_df$id), function(i){
  
  id = clear_df$id[i] #get planet id number
  filename = paste0(outdir,'/',id,'_',dir_string,'.tif')
  if(!file.exists(filename)){
    
    #if(clear_df$clear_percent[i] != 100){
      id_rasters = rasters[str_detect(rasters, id)]
      udm = rast(id_rasters[str_detect(id_rasters, 'udm')])
      scene = rast(id_rasters[!str_detect(id_rasters, 'udm')])
      
      #apply mask
      scene_m = mask(scene, is.na(udm[['clear']]))
      
      #save file
      writeRaster(scene_m, filename = filename)
    #}
  }
}
, cl = 'future'
)

#plan(sequential)

#----calculate glcm textures for Disturbance index for a range of window sizes----

# windows = c(3,5,7,9,11,13,15,17,21,23,25,27,29,31)
# textures = c(
#   "mean"
#   # , "variance"
#   # , "homogeneity"
#   , "contrast"
#   , "dissimilarity"
#   , "entropy"
#   # , "second_moment"
#   , "correlation"
# )
# 
# dir_string = 'GLCM_DIfullscene'
# GLCMDI_dir = paste0(ps_dir, '/',dir_string)
# dir.check(GLCMDI_dir)
# 
# #plan('multisession', workers = 12)
# 
# indir = DI_dir
# rasters = list.files(indir, full.names = T, recursive = T, pattern = '\\.tif')
# 
# library(raster) #apparently needed for glcm function
# 
# pblapply(1:length(windows), function(j){
#   
#   window_size = c(windows[j], windows[j])
#   window_size_string = paste0(window_size[1],'x',window_size[2])
#   
#   window_dir = paste0(GLCMDI_dir,'/',window_size[1],'x',window_size[2])
#   dir.check(window_dir)
#   
#   indir = DI_dir
#   outdir = window_dir
#   
#   future_lapply(1:nrow(clear_df), function(i){
#     id = clear_df$id[i]
#     filename = paste0(outdir,'/',id,'_',dir_string,'_',window_size_string,'.tif')
#     if(!file.exists(filename)){
#       scene = raster(rasters[str_detect(rasters, id)])
#       
#       #calculate glcm textures based on the disturbance index
#       
#       scene_g = glcm(scene, n_grey = 32, window = window_size, statistics = textures, na_opt = 'center')
#       scene_g = terra::rast(scene_g)
#       names(scene_g) = paste0(names(scene_g),'_',basename(indir),'_',window_size_string)
#       
#       #save file
#       terra::writeRaster(scene_g, filename = filename)
#       
#     }
#     #mark progress
#     print(paste0(basename(filename), ' done, ',i,'/',nrow(clear_df),' for this window size.'))
#     
#   })
# })
# 
# # plan('sequential')
# 
# detach('package:raster', unload = T)
# library(terra)
# 
# # #Fix names
# # glcm_files = list.files(GLCMDI_dir, full.names = T, recursive = T, pattern = '\\.tif$')
# # 
# # pblapply(glcm_files, function(x) {
# #   r = rast(x)
# #   window_size = basename(dirname(x))
# #   names(r) = paste0(textures, '_', DI_dir, '_', window_size)
# #   
# #   # Write to a temporary file first
# #   temp_x <- tempfile(fileext = ".tif")
# #   terra::writeRaster(r, filename = temp_x, overwrite = TRUE)
# #   
# #   # Replace the original file
# #   file.rename(temp_x, x)
# # }
# # ,cl = 'future'
# # )
# 
# plan('sequential')

#----Calculate a bunch of spectral indices----
dir_string = 'VegIndices'
VI_dir = paste0(ps_dir, '/', dir_string)
dir.check(VI_dir)

# plan('multisession', workers = 12)

indir = udm_masked_dir
outdir = VI_dir

rasters = list.files(indir, full.names = T, recursive = T, pattern = '\\.tif')

pblapply(1:nrow(clear_df), function(i){
  id = clear_df$id[i]
  filename = paste0(outdir,'/',id,'_',dir_string,'.tif')
  if(!file.exists(filename)){
    scene = rast(rasters[str_detect(rasters, id)])
    
    #get maximum pixel value iacross all bands to use as scale factor
    
    maxes = global(scene,'max',na.rm=T)
    sf = max(maxes[['max']]) #maximum value across all bands in a scene
    
    sf = 65535#maximum theoretical value for a 16-bit pixel
    
    #calculate a bunch of vegetation indices
    scene_vi = spectralIndices(img = scene, blue = 'blue', green = 'green', red = 'red', nir = 'nir', redEdge1 = 'rededge'
                               ,scaleFactor = sf, skipRefCheck = T)
    
    #save file
    terra::writeRaster(scene_vi, filename = filename)
  }
}
, cl = 'future'
)

plan('sequential')


#---- normalize vegetation indices using z-score

#----Calculate change in DI metrics ----


# dir_string = 'DeltaDIfullscene'
# DeltaDI_dir = paste0(ps_dir, '/', dir_string)
# dir.check(DeltaDI_dir)
# 
# indirs = c(DI_dir, GLCMDI_dir, VI_dir)
# outdir = DeltaDI_dir
# 
# rasters = unlist(lapply(indirs, list.files, full.names = TRUE, recursive = T, pattern = '\\.tif$'))
# 
# 
# plan('multisession', workers = 4)
# 
# future_lapply(1:nrow(clear_df), function(i){
#   id = clear_df$id[i]
#   filename = paste0(outdir,'/',id,'_',dir_string,'.tif')
#   
#   reference_id = clear_df$id[1]
#   # refscene1 = rast(rasters1[str_detect(rasters1, reference_id)])
#   # refscene2 = rast(rasters2[str_detect(rasters2, reference_id)])
#   
#   # reference_scene = rast(rasters[str_detect(rasters, reference_id)])
#   ref_l = lapply(rasters[str_detect(rasters, reference_id)], rast)
#   reference_scene = rast(ref_l)
#   
#   if(!file.exists(filename)){
#     
#     r_l = lapply(rasters[str_detect(rasters, id)], rast)
#     scene = rast(r_l)
#     
#     
#     #calculate delta metrics
#     delta_scene = reference_scene - scene
#     
#     #save file
#     terra::writeRaster(delta_scene, filename = filename)
#   }
# })
# 
# plan('sequential')
# 
# l = list.files(DeltaDI_dir, full.names = T)
# r = rast(l[length(l)-1])



