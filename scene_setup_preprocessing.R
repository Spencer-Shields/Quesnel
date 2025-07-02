#----packages----
library(tidyverse)
library(terra)
library(RStoolbox)
library(tools)
library(sf)
library(pbapply)
library(parallel)
library(future)
library(future.apply)
library(tidyterra)
library(Rcpp)
library(tictoc)
library(arrow)
library(data.table)
library(car)
library(nortest)
library(pROC)
library(varSel)
library(tseries)
library(lme4)
library(glmmTMB)
#extra functions
source('helper_functions.R')
#functions for checking or visualizing radiometric consistency
source('https://raw.githubusercontent.com/Spencer-Shields/planetscope_radiometric_correction/refs/heads/main/check_radiometric_consistency.R')


#----load data----

#set directories
list.dirs('data')
ps_dir = 'data/planet_scenes'
raw_dir = paste0(ps_dir,'/raw')

#load thinning block vector data
blocks = st_read('data/Quesnel_thinning/12l_12n_bdy.geojson')

thinning_block_ids = blocks$BLOCKNUM #get ids of thinning blocks

block_ids = blocks$BLOCKNUM

#load harvest dates
harvest_dates_file = 'data/Quesnel_thinning/harvest_dates.csv'
harvest_dates_df = read_csv(harvest_dates_file)

#load Otsu thresholds
otsu_df = read_csv("data/Quesnel_thinning/Otsu_change_thresholds.csv")

#load LiDAR data
lidar_dir = 'data/Quesnel_thinning/chm_change_scenes'
lidar_files = list.files(lidar_dir, full.names = T, recursive = T)
lidar_ids = basename(file_path_sans_ext(lidar_files))
names(lidar_files) = lidar_ids

#load non-vegetation mask data
nonveg_mask_dir = 'data/Quesnel_thinning/nonveg_mask_pre=0.5_post=0.3'
nonveg_mask_files = list.files(nonveg_mask_dir, full.names = T)

#get dataframe of scene metadata
meta_df = ps_meta(raw_dir) |> distinct()
meta_df[['acquisition_date']] = as.Date(sapply(meta_df$acquired, function(x){
  s = str_split_1(x, 'T')[1]
  # d = as.Date(s, format = "%Y-%m-%d")
  return(s)
}))
head(meta_df |> select(acquired, acquisition_date))
meta_df = meta_df |>
  mutate(acquisition_year = year(acquisition_date),
         acquisition_month = month(acquisition_date),
         acquisition_day = day(acquisition_date)) %>%
  mutate(acquisition_month_day = paste0(acquisition_month,'-',acquisition_day))


#filter dataframe to remove undesirable scenes

meta_df_filtered = meta_df %>%
  # filter(clear_percent >= clear_threshold) %>% #get data where the percent of clear pixels is greater than or equal to a threshold
  filter(cloud_cover <= 60) |>
  filter(quality_category == 'standard') |> #get standard quality data
  filter(acquisition_year > 2020)


# #see temporal distribution of filtered images
{
  # filtered_dates_bymonth = meta_df_filtered %>%
  #   group_by(acquisition_year, acquisition_month) %>%
  #   summarise(Number_of_scenes = n())
  # 
  # ggplot()+
  #   geom_col(data = filtered_dates_bymonth, aes(x = acquisition_month, y = Number_of_scenes))+
  #   facet_wrap(vars(acquisition_year), ncol = 1) +
  #   ggtitle(paste0('Scenes with >=',clear_threshold,'% clear pixels'))
  # 
  # #see temporal distribution of ALL PlanetScope scenes
  # meta_df_dates_bymonth = meta_df %>%
  #   group_by(acquisition_year, acquisition_month) %>%
  #   summarise(Number_of_scenes = n())
  # 
  # ggplot()+
  #   geom_col(data = meta_df_dates_bymonth, aes(x = acquisition_month, y = Number_of_scenes))+
  #   facet_wrap(vars(acquisition_year), ncol = 1) +
  #   ggtitle('All scenes')
  }

#----preprocess scenes (avoid saving intermediate data)----

ids = meta_df_filtered$id
raw_rasters = list.files(raw_dir, recursive = T, full.names = T, pattern = '\\.tif$')

if(crs(blocks) != crs(rast(raw_rasters[1]))){ #reproject blocks if necessary to match raster crs
  blocks = st_transform(blocks, crs = crs(rast(raw_rasters[1])))
}

scale_factor = 65535 #scale factor (max value of 16-bit pixel)

scale_up_factor = 10000 #value to multiply indices by in order to avoid floating point errors

raw_bands = names(rast(raw_rasters[!str_detect(raw_rasters,'udm')][1])) #get names of bands in scene rasters (not indices)
indices_to_use = c('BI', 'MSAVI','SR','Hue','NDVI') #indices to test, determined by corr matrix analysis and PCA (below)
feats_to_use = c(raw_bands, indices_to_use)

nonnorm_string = 'Non-normalized'
nonnorm_dir = paste0(ps_dir,'/', nonnorm_string,'_REDO')
dir.check(nonnorm_dir)

z_string = 'Z'
z_dir = paste0(ps_dir,'/', z_string,'_REDO')
dir.check(z_dir)

zr_string = 'Zrobust'
zr_dir = paste0(ps_dir, '/', zr_string,'_REDO')
dir.check(zr_dir)

sm_string = 'SM'
sm_dir = paste0(ps_dir,'/',sm_string,'_REDO')
dir.check(sm_dir)

dirs = c(z_dir, zr_dir, sm_dir, nonnorm_dir)

# dir_strings = c('Z', 'Zrobust', 'SM', 'Non-normalized')
# dirs = paste0(ps_dir,'/',dir_strings,'_REDO')
# lapply(dirs,dir.check)

pblapply(dirs, function(x){ #make subdirectories to store scenes by block
  pblapply(thinning_block_ids, function(y){
    d = paste0(x,'/',y)
    dir.check(d)
  })
})

if(length(list.files(dirs,recursive = T,pattern='\\.tif$')) < length(dirs)*length(thinning_block_ids)*length(ids)){
  #process rasters
  
  cl = makeCluster(8)
  plan('cluster', workers = cl)
  
  future_lapply(1:length(ids), function(i){
    
    id = ids[i]
    
    pblapply(1:length(thinning_block_ids), function(j){
      
      block_dir = thinning_block_ids[j]
      
      #Non-normalized raster
      nonnorm_file = paste0(nonnorm_dir,'/',block_dir,'/',id,'.tif')
      print(paste('Processing', nonnorm_file))
      if(!file.exists(nonnorm_file)){
        
        #load raster
        f_r = raw_rasters[str_detect(raw_rasters, id) & !str_detect(raw_rasters, 'udm')][1]
        r = rast(f_r)
        r = r[[feats_to_use]] #subset raster to use features described above
        
        #mask using UDM
        f_udm = raw_rasters[str_detect(raw_rasters, id) & str_detect(raw_rasters, 'udm')][1]
        udm = rast(f_udm)
        
        r = mask(r, udm[['cloud']]) #remove cloud
        r = mask(r, udm[['shadow']]) #remove shadow
        r = mask(r, udm[['haze_heavy']]) #remove heavy haze
        
        #crop to block
        block = blocks |> filter(BLOCKNUM == block_dir)
        r = crop(r, block, mask = T)
        
        #use nonveg mask
        nvf = nonveg_mask_files[which(str_detect(nonveg_mask_files,block_dir))]
        nvm = rast(nvf)
        nvm_rs = resample(nvm,r)
        r = mask(r, nvm_rs)
        
        
        #calculate vegetation indices
        
        si = spectralIndices(img = r, blue = 'blue', green = 'green', red = 'red', nir = 'nir', redEdge1 = 'rededge'
                             ,scaleFactor = scale_factor, skipRefCheck = T)
        
        si = si*scale_up_factor
        
        vi = rast.batch.functions(r, include_input = F, fl = c(
          hue.vi,
          gcc.vi,
          ndgr.vi,
          bi.vi,
          ci.vi,
          cri550.vi,
          gli.vi,
          # tvi.vi,
          varig.vi
        ))
        
        vi = vi*scale_up_factor
        
        r = c(r, si, vi)
        
        #save file
        terra::writeRaster(r, nonnorm_file, datatype='FLT8S')
        
      } else {
        r = rast(nonnorm_file)
      }
      
      #Z-score raster
      z_file = paste0(z_dir,'/',block_dir,'/',id,'.tif')
      # print(paste('Processing',z_file))
      if(!file.exists(z_file)){
        z_r = rast(lapply(r, z_rast))
        writeRaster(z_r, z_file)
      }
      
      #robust Z-score raster
      zr_file = paste0(zr_dir,'/',block_dir,'/',id,'.tif')
      # print(paste('Processing',zr_file))
      if(!file.exists(zr_file)){
        zr_r = rast(lapply(r, function(x)z_rast(x,robust = T)))
        writeRaster(zr_r, zr_file)
      }
      
      #softmax raster
      sm_file = paste0(sm_dir,'/',block_dir,'/',id,'.tif')
      # print(paste('Processing',sm_file))
      if(!file.exists(sm_file)){
        
        #scale bands for use with softmax
        rs = rast(lapply(1:nlyr(r), function(k){
          b = r[[k]]
          nom = names(r)[k]
          if(nom %in% raw_bands){
            b = b/scale_factor
          } else {
            b 
          }
          return(b)
        }))
        
        sm_r = softmax(rs, append_name = F)
        sm_r = sm_r * scale_up_factor
        writeRaster(sm_r, sm_file, datatype='FLT8S')
      }
    })
    
  }
  # ,cl = 'future'
  # ,future.seed=T
  )
  
  # plan('sequential')
}

