library(tidyverse)
library(terra)
library(sf)
source('helper_functions.R')
library(future)
library(future.apply)
library(tools)
library(data.table)

#----load directories, block data----

ps_dir = 'data/planet_scenes'

blocks = st_read('data/Quesnel_thinning/12l_12n_bdy.gpkg')
no_change = st_read('Arc project/NoChangeStands.shp')

no_change = no_change %>% rename(geom = geometry)

blocks = rbind(blocks, no_change)

#
#---- process z rasters: crop to block, get global data, get timeseries stats ----
{
  dir_string = 'Z-blocks' #name of new directory
  indir = 'data/planet_scenes/Z-score' #name of directory with input data
  csv_string = 'z_score_timeseries_data_byblock.csv' #directory for dataframe of output statistics
  
  outdir = paste0(ps_dir,'/', dir_string)
  dir.check(outdir)
  
  data_filename = paste0(outdir,'/',csv_string)
  if(!file.exists(data_filename)){
    
    files = list.files(indir, full.names = T)
    
    blocks_p = st_transform(blocks, crs = crs(rast(files[1]))) #reproject blocks to match rasters
    
    plan('multisession', workers = 8)
    
    results_df_l = pblapply(1:nrow(blocks_p), function(i){
      #iterate over blocks, make directory for output
      block = blocks_p[i,]
      block_id = block$BLOCKNUM
      block_dir = paste0(outdir, '/',block_id)
      dir.check(block_dir)
      
      block_v = vect(block)
      
      #process raster data
      b_df_l = future_lapply(1:length(files), function(j){
        id = basename(file_path_sans_ext(files[j]))
        filename = paste0(block_dir,'/',id,'_',block_id,'.tif')
        
        r = rast(files[j]) #load raster
        r_c = crop(r, block) #crop to block
        
        if(!file.exists(filename)){writeRaster(r_c, filename)} #save cropped file if it doesn't exist
        
        d = global(r_c, fun = c("max", "min", "mean", "sum", "range", "rms", "sd", "std", "isNA", "notNA"), na.rm=T) #calculate global stats for each band
        d[['rast_id']] = paste0(id,'_',block_id) #assign raster/site id to each row of stats
        return(d)
      })
      b_df = bind_rows(b_df_l)
      return(b_df)
    })
    
    plan('sequential')
    
    results_df = bind_rows(results_df_l)
    
    write.csv(results_df, data_filename)
  }
}

z_data = fread(data_filename)
z_data[['dataset']] = 'Z-score'
#----process z_score data----

# z_data = fread(data_filename)
# z_data[['block_id']] = sapply(z_data$rast_id, function(x)str_split_1(x, '_Z-score_'))[2]
# z_data[['r_meta']] = sapply(z_data$rast_id, function(x)str_split_1(x, '_Z-score_'))[1]
# 
# z_data = z_data %>%
#   mutate(
#     year = str_sub(r_meta, 1, 4),
#     month = str_sub(r_meta, 5, 6),
#     day = str_sub(r_meta, 7, 8)
#   ) %>%
#   mutate(acquisition_date = as.Date(paste0(year,'-',month,'-',day))) %>%
#   mutate(julian_day = yday(acquisition_date)) %>%
#   mutate(var_name = str_remove(V1, "_Z-score.*"))
# 
# ggplot(z_data %>% filter())


#---- process softmax rasters: crop to block, get global data, get timeseries stats ----
{
  dir_string = 'softmax-blocks' #name of new directory
  indir = 'data/planet_scenes/softmax' #name of directory with input data
  csv_string = 'softmax_timeseries_data_byblock.csv' #directory for dataframe of output statistics
  
  outdir = paste0(outdir,'/', dir_string)
  dir.check(outdir)
  
  data_filename = paste0(outdir,'/',csv_string)
  if(!file.exists(data_filename)){
    
    files = list.files(indir, full.names = T)
    
    blocks_p = st_transform(blocks, crs = crs(rast(files[1]))) #reproject blocks to match rasters
    
    plan('multisession', workers = 8)
    
    results_df_l = pblapply(1:nrow(blocks_p), function(i){
      #iterate over blocks, make directory for output
      block = blocks_p[i,]
      block_id = block$BLOCKNUM
      block_dir = paste0(outdir, '/',block_id)
      dir.check(block_dir)
      
      block_v = vect(block)
      
      #process raster data
      b_df_l = future_lapply(1:length(files), function(j){
        id = basename(file_path_sans_ext(files[j]))
        filename = paste0(block_dir,'/',id,'_',block_id,'.tif')
        
        r = rast(files[j]) #load raster
        r_c = crop(r, block) #crop to block
        
        if(!file.exists(filename)){writeRaster(r_c, filename)} #save cropped file if it doesn't exist
        
        d = global(r_c, fun = c("max", "min", "mean", "sum", "range", "rms", "sd", "std", "isNA", "notNA"), na.rm=T) #calculate global stats for each band
        d[['rast_id']] = paste0(id,'_',block_id) #assign raster/site id to each row of stats
        return(d)
      })
      b_df = bind_rows(b_df_l)
      return(b_df)
    })
    
    plan('sequential')
    
    results_df = bind_rows(results_df_l)
    
    write.csv(results_df, data_filename)
  }
}

softmax_data = fread(data_filename)
softmax_data[['dataset']] = 'softmax'

#----process non-normalized time series----
{
  dir_string = 'nonorm-blocks' #name of new directory
  indir = 'data/planet_scenes/VegIndices' #name of directory with input data
  csv_string = 'nonnorm_timeseries_data_byblock.csv' #directory for dataframe of output statistics
  
  outdir = paste0(outdir,'/', dir_string)
  dir.check(outdir)
  
  data_filename = paste0(outdir,'/',csv_string)
  if(!file.exists(data_filename)){
    
    files = list.files(indir, full.names = T)
    
    blocks_p = st_transform(blocks, crs = crs(rast(files[1]))) #reproject blocks to match rasters
    
    plan('multisession', workers = 8)
    
    results_df_l = pblapply(1:nrow(blocks_p), function(i){
      #iterate over blocks, make directory for output
      block = blocks_p[i,]
      block_id = block$BLOCKNUM
      block_dir = paste0(outdir, '/',block_id)
      dir.check(block_dir)
      
      block_v = vect(block)
      
      #process raster data
      b_df_l = future_lapply(1:length(files), function(j){
        id = basename(file_path_sans_ext(files[j]))
        filename = paste0(block_dir,'/',id,'_',block_id,'.tif')
        
        r = rast(files[j]) #load raster
        r_c = crop(r, block) #crop to block
        
        if(!file.exists(filename)){writeRaster(r_c, filename)} #save cropped file if it doesn't exist
        
        d = global(r_c, fun = c("max", "min", "mean", "sum", "range", "rms", "sd", "std", "isNA", "notNA"), na.rm=T) #calculate global stats for each band
        d[['rast_id']] = paste0(id,'_',block_id) #assign raster/site id to each row of stats
        return(d)
      })
      b_df = bind_rows(b_df_l)
      return(b_df)
    })
    
    plan('sequential')
    
    results_df = bind_rows(results_df_l)
    
    write.csv(results_df, data_filename)
  }


nonnorm_data = fread(data_filename)
nonnorm_data[['dataset']] = 'VegIndices'
}

#----process combined data----

combined_df = bind_rows(list(z_data, softmax_data, nonnorm_data))

combined_df = combined_df %>%
  mutate(
    block_id = str_remove(rast_id, paste0(".*_", dataset, "_")),
    r_meta = str_remove(rast_id,  paste0("_", dataset, ".*"))
  ) %>%
  mutate(
    year = str_sub(r_meta, 1, 4),
    month = str_sub(r_meta, 5, 6),
    day = str_sub(r_meta, 7, 8)
  ) %>%
  mutate(acquisition_date = as.Date(paste0(year,'-',month,'-',day))
    ) %>%
  mutate(julian_day = yday(acquisition_date)) %>%
  mutate(var_name = str_remove(V1, paste0("_", dataset, ".*"))) %>%
  mutate(var_name = str_remove(var_name, "\\.\\.\\..*"))

var_summ = combined_df %>%
  group_by(var_name, dataset) %>%
  summarise(count = n())
  
#----plotting vegetation indices over time----

#correlation matrix of vegetation indices


inds = c(unique(combined_df$var_name[!combined_df$var_name %in% c('EVI2','SR')])) #vector of indices to look at
blocks_ = c("12L_C5",    "12L_D345",  "12L_C4",    "12L_C7",    "12L_B8C3",  "12N_T3",    "12N_1X1W",  "12N_V1",   
             "NoChange1")

subset_df = combined_df %>% filter(var_name %in% inds, block_id %in% blocks_)
ggplot(subset_df, aes(x = acquisition_date, y = mean, color = block_id)) +
  geom_point()+
  geom_line()+
  # geom_ribbon(aes(x = acquisition_date, ymin = mean-sds, ymax = mean+sds))+
  facet_grid(rows = vars(dataset), cols = vars(var_name), scales = 'free')+
  scale_x_date(
    date_breaks = "1 year",       # Major tick marks every year
    date_labels = "%Y",           # Year labels on major ticks
    date_minor_breaks = "1 month" # Minor tick marks every month
  )




