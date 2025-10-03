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
library(Rcpp)
# library(tictoc)
# library(arrow)
# library(data.table)
# library(car)
# library(nortest)
# library(pROC)
# library(varSel)
# library(tseries)
# library(lme4)
# library(glmmTMB)
#extra functions
source('helper_functions.R')
#functions for checking or visualizing radiometric consistency
source('https://raw.githubusercontent.com/Spencer-Shields/planetscope_radiometric_correction/refs/heads/main/check_radiometric_consistency.R')


#----load data----

#set directories
ps_dir = 'data/planet_scenes'
raw_dir = paste0(ps_dir,'/raw')

#load thinning block vector data
blocks = st_read('data/Quesnel_thinning/12l_12n_bdy.geojson')

thinning_block_ids = blocks$BLOCKNUM #get ids of thinning blocks

block_ids = blocks$BLOCKNUM

#load harvest dates
harvest_dates_file = 'data/Quesnel_thinning/harvest_dates.csv'
harvest_dates_df = read_csv(harvest_dates_file)

if(!file.exists('data/Quesnel_thinning/harvest_dates_plus.csv')){
  harvest_dates_plus = blocks |>
    rename(block_id = BLOCKNUM) |>
    left_join(harvest_dates_df) |>
    mutate(area_m2 = st_area(geometry)) |>
    mutate(area_ha = units::set_units(area_m2, 'ha')) |>
    mutate(area_km2 = units::set_units(area_m2, 'km^2'))
  
  write_csv(harvest_dates_plus |> st_drop_geometry(), 'data/Quesnel_thinning/harvest_dates_plus.csv')
}

#load Otsu thresholds
otsu_df = read_csv("data/Quesnel_thinning/Otsu_change_thresholds.csv")

#load LiDAR data
lidar_dir = 'data/Quesnel_thinning/chm_change_scenes'
lidar_files = list.files(lidar_dir, full.names = T, recursive = T)
lidar_ids = basename(file_path_sans_ext(lidar_files))
names(lidar_files) = lidar_ids

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
  # filter(cloud_cover <= 60) |>
  filter(quality_category == 'standard') |> #get standard quality data
  filter(acquisition_year >= 2020)


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

#get ids of raw rasters and list of raw raster files
ids = meta_df_filtered$id
raw_rasters = list.files(raw_dir, recursive = T, full.names = T, pattern = '\\.tif$')

#reproject blocks
if(terra::crs(vect(blocks)) != terra::crs(rast(raw_rasters[1]))){ #reproject blocks if necessary to match raster crs
  blocks = st_transform(blocks, crs = terra::crs(rast(raw_rasters[1])))
}

#make blocks spatvector and wrap them for parallelization
blocks_wrapped = blocks |> vect() |> wrap()

scale_factor = 10000 #what the reflectance values get divided by to facilitate calculating indices and softmax values

#create directories for each type of dataset to store the processed rasters
nonnorm_string = 'Non-normalized'
nonnorm_dir = paste0(ps_dir,'/', nonnorm_string,'_20250910')
dir.check(nonnorm_dir)

z_string = 'Z'
z_dir = paste0(ps_dir,'/', z_string,'_20250910')
dir.check(z_dir)

zr_string = 'Zrobust'
zr_dir = paste0(ps_dir, '/', zr_string,'_20250910')
dir.check(zr_dir)

sm_string = 'SM'
sm_dir = paste0(ps_dir,'/',sm_string,'_20250910')
dir.check(sm_dir)

dirs = c(z_dir, zr_dir, sm_dir, nonnorm_dir)

#make subdirectories within dataset directories to store scenes for each thinning block
pblapply(dirs, function(x){
  pblapply(thinning_block_ids, function(y){
    d = paste0(x,'/',y)
    dir.check(d)
  })
})

#check if there are fewer files in the output folders than expected, engage in preprocessing if yes
if(length(list.files(dirs,recursive = T,pattern='\\.tif$')) < length(dirs)*length(thinning_block_ids)*length(ids)){
  
  #define dummy raster to save if the actual raster ends up being blank or invalid (saves time reloading data later)
  r_dummy = rast(nrow=1,ncol=1)
  values(r_dummy) = NA #r_dummy has a single cell, the value is NA
  r_dummy_wrapped = wrap(r_dummy)
  
  #process rasters
  
  clust = makeCluster(16)
  plan('cluster', workers = clust)
  
  future_lapply(1:length(ids), function(i){
  # pblapply(1:length(ids), function(i){
    
    
    #get file id
    id = ids[i]
    
    #unwrap r_dummy for workers
    r_dummy = unwrap(r_dummy_wrapped)
    
    pblapply(1:length(thinning_block_ids), function(j){
      
      block_dir = thinning_block_ids[j]
      
      #Non-normalized raster
      nonnorm_file = paste0(nonnorm_dir,'/',block_dir,'/',id,'.tif')
      nonnorm_dummy_file = paste0(file_path_sans_ext(nonnorm_file),'_NULL.tif')
      cat('Processing',i,'---',j,'; ', nonnorm_file, '\n')
      if(!file.exists(nonnorm_file) & !file.exists(nonnorm_dummy_file)){
        
        #load raster
        f_r = raw_rasters[str_detect(raw_rasters, id) & !str_detect(raw_rasters, 'udm')][1]
        r = rast(f_r)
        # r = r[[feats_to_use]] #subset raster to use features described above
        
        
        #mask using UDM
        f_udm = raw_rasters[str_detect(raw_rasters, id) & str_detect(raw_rasters, 'udm')][1]
        udm = rast(f_udm)
        
        r = mask(r, udm[['cloud']]) #remove cloud
        r = mask(r, udm[['shadow']]) #remove shadow
        r = mask(r, udm[['haze_heavy']]) #remove heavy haze
        
        #load thinning block boundary
        blocks = unwrap(blocks_wrapped)
        block = blocks[blocks$BLOCKNUM==block_dir,]
        
        
        #check if thinning block and raster overlap, save dummy file if they do not, proceed otherwise
        if(!any(relate(r,block,'intersects'))){
          
          terra::writeRaster(r_dummy, nonnorm_dummy_file)
          
        } else {
          
          #crop to block
          r = crop(r, block, mask = T)
          
          #check if raster only contains NA values, save dummy file if yes, proceed with preprocessing otherwise
          g = global(r, 'isNA', na.rm=T)
          
          if(g[1,1]==ncell(r)){
            
            terra::writeRaster(r_dummy, nonnorm_dummy_file)
            
          } else {
            
            #calculate vegetation indices, combine all features into single spatraster
            {
              #rescale raster
              r = r/scale_factor
              
              #calculate vegetation indices
              si = spectralIndices(img = r, blue = 'blue', green = 'green', red = 'red', nir = 'nir', redEdge1 = 'rededge'
                                   # , skipRefCheck = T
              )
              #remove EVI2 (since spectralIndices doesn't calculate it well for some reason)
              si = subset(si, 'EVI2', negate=T)
              
              #calculate additional indices with custom functions
              vi = rast.batch.functions(r, include_input = F,
                                        #band mapping
                                        coastal_blue=1, blue=2, green_i=3, green=4, yellow=5, red=6, rededge=7, nir=8, 
                                        #list of indices to calculate
                                        fl = c(
                                          hue.vi,
                                          gcc.vi,
                                          ndgr.vi,
                                          bi.vi,
                                          ci.vi,
                                          cri550.vi,
                                          gli.vi,
                                          # tvi.vi, #no TVI since spectralIndices handles this
                                          varig.vi,
                                          evi2.vi,
                                          yndvi.vi
                                        ))
              
              r = c(r, si, vi)
            }
            
            #save file
            terra::writeRaster(r, nonnorm_file, datatype='FLT8S')
          }
        }
        
      } else { #if the file doesn't exist make it, if it does exist load it for further processing checks
        
        #load the dummy file if nonnorm file does not exist
        if(file.exists(nonnorm_file)){
                   r = rast(nonnorm_file)
        } else {
                  r = rast(nonnorm_dummy_file)
          }
      }
      
      #Z-score raster
      z_file = paste0(z_dir,'/',block_dir,'/',id,'.tif')
      z_dummy_file = paste0(file_path_sans_ext(z_file),'_NULL.tif')
      
      # print(paste('Processing',z_file))
      if(!file.exists(z_file)&!file.exists(z_dummy_file)){
        
        if(ncell(r) == ncell(r_dummy)){
          writeRaster(r, z_dummy_file)
        } else {
          z_r = z_rast(r)
          writeRaster(z_r, z_file)
        }
      }
      
      #robust Z-score raster
      zr_file = paste0(zr_dir,'/',block_dir,'/',id,'.tif')
      zr_dummy_file = paste0(file_path_sans_ext(zr_file),'_NULL.tif')
      
      # print(paste('Processing',zr_file))
      if(!file.exists(zr_file)&!file.exists(zr_dummy_file)){
        
        if(ncell(r) == ncell(r_dummy)){
          writeRaster(r, zr_dummy_file)
        } else {
          zr_r = z_rast(r,robust = T)
          writeRaster(zr_r, zr_file)
        }
      }
      
      #softmax raster
      sm_file = paste0(sm_dir,'/',block_dir,'/',id,'.tif')
      sm_dummy_file = paste0(file_path_sans_ext(sm_file),'_NULL.tif')
      
      # print(paste('Processing',sm_file))
      if(!file.exists(sm_file)&!file.exists(sm_dummy_file)){
        
        if(ncell(r) == ncell(r_dummy)){
          writeRaster(r, sm_dummy_file)
        } else{
          
          sm_r = softmax(r, append_name = F)
          # sm_r = sm_r * scale_up_factor
          writeRaster(sm_r, sm_file, datatype='FLT8S')
        }
      }
      
    })
    
  }
  # ,cl = 'future'
  # ,future.seed=T
  )
  
  stopCluster(clust)
  plan('sequential')
}

#resample data to ensure consistent geometry (because PlanetScope analytic scenes don't have consistent resolution for some damn reason)
{
  
  #figure out what spatial resolution is most common
  {
  # #get list of spatial resolutions
  # fl = list.files(dirs, recursive = T, full.names = T)
  # fl = fl[!str_detect(fl, 'NULL')]
  # 
  # # clust = makeCluster(4)
  # # plan('cluster', workers = clust)
  # res_l = pblapply(1:length(fl), function(i){
  #   cat('\r',i)
  #   r = rast(fl[i])
  #   res = terra::res(r)
  #   return(res)
  # }
  # # ,cl='future'
  # )
  # # stopCluster(clust)
  # # plan('sequential')
  # 
  # res_df = tibble(x = pbsapply(res_l, function(r)r[1]), 
  #                 y = pbsapply(res_l, function(r)r[2]))
  # 
  # res_summ = res_df |>
  #   group_by(x, y) |>
  #   dplyr::summarise(n = n())
  # 
  # ggplot(res_summ |> pivot_longer(cols = c('x', 'y'), names_to = 'dim', values_to = 'val'))+
  #   geom_col(aes(x = val, y = n, fill = dim), alpha = 0.5) +
  #   facet_grid(rows=vars(dim))
  # 
  # 
  # # reproject rasters and check again
  # fl2 = fl[str_detect(fl, 'Non-normalized')]
  # 
  # clust = makeCluster(4)
  # plan('cluster', workers = clust)
  # res_l = pblapply(1:length(fl2), function(i){
  #   cat('\r',i)
  #   r = rast(fl2[i])
  #   r = project(r, 'EPSG:3005')
  #   res = terra::res(r)
  #   return(res)
  # }
  # ,cl='future'
  # )
  # stopCluster(clust)
  # plan('sequential')
  # 
  # res_df = tibble(x = pbsapply(res_l, function(r)r[1]), 
  #                 y = pbsapply(res_l, function(r)r[2]))
  # 
  # res_summ = res_df |>
  #   group_by(x, y) |>
  #   dplyr::summarise(n = n())
  # 
  # ggplot(res_summ |> pivot_longer(cols = c('x', 'y'), names_to = 'dim', values_to = 'val'))+
  #   geom_col(aes(x = val, y = n, fill = dim), alpha = 0.5) +
  #   facet_grid(rows=vars(dim))
  # 
  # 
  }
  
  #define target spatial resolution (m, BC Albers)
  # target_res = 3.764306e-05
  target_res=3
  target_crs = 'EPSG:3005'
  
  #reproject blocks
  blocks_p = blocks |>
    vect()|>
    project(target_crs)
  
  res_dirs = paste0(dirs, '_',target_res,'m')
  
  #check if 
  if(length(list.files(res_dirs,recursive = T)) < length(list.files(dirs, recursive=T))){
    
    #make template raster for resampling scenes
    aoi_ext = ext(vect(blocks) |> project(target_crs))
    
    snap_extent <- function(ext, res) {
      xmin <- floor(ext[1] / res) * res
      xmax <- ceiling(ext[2] / res) * res
      ymin <- floor(ext[3] / res) * res
      ymax <- ceiling(ext[4] / res) * res
      ext(xmin, xmax, ymin, ymax)
    }
    
    aoi_ext_snapped <- snap_extent(aoi_ext, target_res)
    
    aoi_rast <- rast(aoi_ext_snapped, resolution = target_res, crs = target_crs)
    aoi_rast_wrapped = terra::wrap(aoi_rast)
    # res(aoi_rast)  # should be exactly 3, 3
    
    #check if dataset directories exist, make them if they don't
    lapply(res_dirs, dir.check)
    
    #check if subdirectories for blocks exist for each dataset, make them if they don't
    lapply(res_dirs, function(d){
      lapply(block_ids, function(b){
        dir = paste0(d,'/',b)
        dir.check(dir)
      })
    })
    
    #define dummy raster to save if the actual raster ends up being blank or invalid (saves time reloading data later)
    r_dummy = rast(nrow=1,ncol=1)
    values(r_dummy) = NA #r_dummy has a single cell, the value is NA
    r_dummy_wrapped = wrap(r_dummy)
    
    #wrap blocks for parallel
    blocks_p_wrapped = blocks_p |>
      wrap()
      
    
    
    #reproject rasters to match the template
    
    fl = list.files(dirs, recursive = T, full.names = T , pattern = '\\.tif$')
    
    clust = makeCluster(8)
    plan('cluster', workers = clust)
    
    future_lapply(1:length(fl), function(i){
      f = fl[i]
      
      reprojected_f = str_replace(f, dirs[which(startsWith(f, dirs))]
                                  ,res_dirs[which(startsWith(f, dirs))])
      if(!file.exists(reprojected_f)){
        if(str_detect(basename(reprojected_f), 'NULL')){
          writeRaster(unwrap(r_dummy_wrapped), reprojected_f)
        } else {
          
          r = rast(f)
          r_p = project(r, terra::unwrap(aoi_rast_wrapped))
          
          b = basename(dirname(f))
          blocks_ = unwrap(blocks_p_wrapped)
          block = blocks_[blocks_$BLOCKNUM==b,]
          r_p = crop(r_p, block, mask=T)
          
          writeRaster(r_p, reprojected_f)
        }
      }
    }
    # ,cl = 'future'
    )
    
    stopCluster(clust)
    plan('sequential')
    
  } 
  
  #redefine dirs as res_dirs in order to propagate changes through the processing 
  dirs = res_dirs
  
  
  
}

