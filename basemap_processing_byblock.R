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
#extra functions
source('helper_functions.R')
#functions for checking or visualizing radiometric consistency
source('https://raw.githubusercontent.com/Spencer-Shields/planetscope_radiometric_correction/refs/heads/main/check_radiometric_consistency.R')

#set up number of cores to be used for parallel processing
cores = detectCores()

#----data preprocessing
{
  #----load and examine basemap data----
  bm_dir = 'data/planet_basemaps/global_monthly'
  
  raw_dir = paste0(bm_dir, '/raw')
  # indir = 'data/planet_basemaps/ps_monthly_sen2_normalized_analytic_8b_sr/raw'
  # indir = 'data/planet_basemaps/global_quarterly/raw'
  
  rast_files = list.files(raw_dir, full.names = T, recursive = T
                          , pattern = '\\.tif$'
  )
  r_l = lapply(rast_files,rast) #get rasters
  
  
  #----get thinning blocks and aoi----
  ps_dir = 'data/planet_scenes'
  blocks = st_read('data/Quesnel_thinning/12l_12n_bdy.gpkg')
  no_change = st_read('Arc project/NoChangeStands_conifer.shp')
  no_change = no_change %>% rename(geom = geometry)
  blocks = rbind(blocks, no_change)
  blocks_p = st_transform(blocks, crs = crs(r_l[[1]]))
  
  aoi = st_read('data/AOI_fullsite_wgs84.kml') %>%
    st_transform(crs = crs(r_l[[1]]))
  
  ggRGB(r_l[[1]])+
    geom_sf(data = blocks_p, aes(color = BLOCKNUM))
  # 
  # ggRGB(r_l[[2]]) +
  #   geom_sf(data = blocks_p, aes(color = BLOCKNUM))
  
  #----mosaic basemap quads from each month, crop to aoi----
  
  #get subdirectories that contain monthly basemap quads
  month_dirs = list.dirs(raw_dir, full.names = T)
  month_dirs = month_dirs[month_dirs != raw_dir] #exclude parent directory
  
  desired_months = c('01','02','03','04', '05','06', '07', '08', '09', '10', '11','12')
  month_str = paste0('_',desired_months,'_')
  month_pattern = paste0("_(", paste(desired_months, collapse = "|"), ")_")
  
  desired_month_dirs = month_dirs[str_detect(month_dirs, month_pattern)]
  
  mosaic_string = 'CroppedMosaics'
  outdir = paste0(bm_dir, '/', mosaic_string)
  dir.check(outdir)
  
  
  # cl = makeCluster(ceiling(cores/3))
  # plan('cluster', workers = cl)
  
  future_lapply(desired_month_dirs, function(d){
    output_file = paste0(outdir, '/', basename(d),'.tif')
    if(!file.exists(output_file)){
      f = list.files(d, full.names = T, recursive = T, pattern = '\\.tif$')
      rl = lapply(f,rast)
      sc = sprc(rl)
      m = mosaic(sc)
      mc = crop(m,aoi)
      writeRaster(mc, output_file)
    }
  })
  print('Mosaic crop done')
  # stopCluster(cl)
  # plan('sequential')
  
  #----calculate vegetation indices----
  
  indir = outdir
  files_list = list.files(indir, full.names = T, recursive = T, pattern = '\\.tif$')
  
  vi_string = 'Indices'
  outdir = paste0(indir,'_', vi_string)
  dir.check(outdir)
  
  # cl = makeCluster(ceiling(cores/4))
  # plan('cluster', workers = cl)
  
  future_lapply(files_list, function(x){
    output_file = paste0(outdir, '/', basename(x),'.tif')
    if(!file.exists(output_file)){
      
      #load raster
      r = rast(x)
      
      #calculate vegetation indices, include original raster bands in the output
      vi = rast.batch.functions(r, include_input = T, fl = c(
        hue.vi,
        gcc.vi,
        ndgr.vi,
        bi.vi,
        ci.vi,
        cri550.vi,
        gli.vi,
        tvi.vi,
        varig.vi
      ))
      names(vi[[1:4]]) = c('blue', 'green', 'red', 'max_DN') #rename original bands
      
      #save raster
      writeRaster(vi, output_file)
      gc()
    }
  })
  
  # print('Calculate VegIndices done')
  # stopCluster(cl)
  # plan('sequential')
  
  
  #----mask data using NTEMS forest cover mask----
  
  # indir = outdir
  # rasters = list.files(indir, full.names = T, pattern = '\\.tif$')
  # 
  # fc_cropped_string = 'ntemsFC'
  # outdir = paste0(indir, '_',fc_cropped_string)
  # dir.check(outdir)
  # 
  # #load NTEMS data
  # ntems_prepped_file = 'data/spencer-ntems/lc_2021_prepped.tif'
  # if(!file.exists(ntems_prepped_file)){
  #   ntems_fc = rast("data/spencer-ntems/spencer-lc") %>%
  #     project(crs(aoi)) %>% #reproject
  #     crop(aoi) #crop to match aoi
  #   ntems_fc_cat = as.factor(ntems_fc) #make raster categorical
  #   ntems_fc_r = resample(ntems_fc_cat, rast(rasters[1]), method = 'near') #resample to match PS basemaps
  #   terra::writeRaster(ntems_fc_r, ntems_prepped_file)
  # }else{
  #   ntems_fc_r = rast(ntems_prepped_file)
  # }
  # 
  # #make mask
  # forest_classes = c('81','210','220','230') #landcover classes which represent forest cover
  # fc_mask = ifel(ntems_fc_r %in% forest_classes, 1, NA) #assign NA values to non-forest landcover types
  # 
  # pblapply(1:length(rasters), function(i){ #not parallelized since spatraster is being passed to function
  #   id = basename(rasters[i])
  #   filename = paste0(outdir, '/', id)
  #   if(!file.exists(filename)){
  #     r = rast(rasters[i])
  #     r_m = mask(r, fc_mask)
  #     terra::writeRaster(r_m, filename)
  #   }
  # })
  # 
  # print('NTEMS mask done')
  
  #----crop data to blocks----
  
  indir = outdir
  rasters = list.files(indir, pattern = '\\.tif$', full.names = T)
  
  clipped_string = 'BlockClipped'
  outdir = paste0(indir, '_', clipped_string)
  nonorm_dir = outdir #save this directory name for later
  dir.check(outdir)
  
  lapply(1:nrow(blocks_p), function(i){ #create subdirectories for every block
    d = paste0(outdir,'/',blocks_p$BLOCKNUM[i])
    dir.check(d)
  })
  
  pblapply(1:length(rasters), function(i){
    r_id = basename(file_path_sans_ext(rasters[i]))
    future_lapply(1:nrow(blocks_p), function(j){
      block_id = blocks_p$BLOCKNUM[j]
      filename = paste0(outdir,'/',block_id,'/',r_id)
      print(paste0('Processing ',filename))
      if(!file.exists(filename)){
        r = rast(rasters[i])
        b = blocks_p[j,]
        r_c = crop(r,b,mask=T)
        terra::writeRaster(r_c, filename)
      }
    })
  })
  
  
  #----Calculate z-scores, softmax, and robust z-scores using cropped data----
  
  indir = outdir
  rasters = list.files(indir, pattern = '\\.tif$', full.names = T, recursive = T)
  
  z_dir_string = 'Z'
  z_dir = paste0(indir, '_', z_dir_string)
  dir.check(z_dir)
  
  lapply(1:nrow(blocks_p), function(i){ #create subdirectories for every block
    d = paste0(z_dir,'/',blocks_p$BLOCKNUM[i])
    dir.check(d)
  })
  
  sm_dir_string = 'SM'
  sm_dir = paste0(indir, '_', sm_dir_string)
  dir.check(sm_dir)
  
  lapply(1:nrow(blocks_p), function(i){ #create subdirectories for every block
    d = paste0(sm_dir,'/',blocks_p$BLOCKNUM[i])
    dir.check(d)
  })
  
  zr_dir_string = 'Zrobust'
  zr_dir = paste0(indir, '_', zr_dir_string)
  dir.check(zr_dir)
  
  lapply(1:nrow(blocks_p), function(i){ #create subdirectories for every block
    d = paste0(zr_dir,'/',blocks_p$BLOCKNUM[i])
    dir.check(d)
  })
  
  sf = 255 #scale factor to use for calculating softmax (i.e. maximum possible band pixel value)
  
  bands = c('blue', 'green', 'red')
  
  cl = makeCluster(ceiling(cores/4))
  plan('cluster', workers = cl)
  
  #calculate z-scores on cropped data
  future_lapply(1:length(rasters), function(i){
    filename = str_replace(rasters[i], #get filename by modifying the name of the parent directory
                           basename(indir), 
                           paste0(basename(indir),'_',z_dir_string)) 
    if(!file.exists(filename)){
      #load scene
      scene = rast(rasters[i])
      
      #calculate and save z-score
      scene_z_l = lapply(scene, z_rast)
      scene_z = rast(scene_z_l)
      names(scene_z) = paste0(names(scene),'_',z_dir_string)
      terra::writeRaster(scene_z, filename, overwrite = T)
    }
  })
  
  #calculate robust z-scores on cropped data
  future_lapply(1:length(rasters), function(i){
    filename = str_replace(rasters[i], #get filename by modifying the name of the parent directory
                           basename(indir), 
                           paste0(basename(indir),'_',zr_dir_string))  
    if(!file.exists(filename)){
      #load scene
      scene = rast(rasters[i])
      
      #calculate and save robust z-score
      scene_z_l = lapply(scene, function(x)z_rast(x, robust = T))
      scene_z = rast(scene_z_l)
      names(scene_z) = paste0(names(scene),'_',zr_dir_string)
      terra::writeRaster(scene_z, filename, overwrite = T)
    }
  })
  
  #calculate softmax on cropped data
  future_lapply(1:length(rasters), function(i){
    filename = str_replace(rasters[i], #get filename by modifying the name of the parent directory
                           basename(indir), 
                           paste0(basename(indir),'_',sm_dir_string))  
    if(!file.exists(filename)){
      #load scene
      scene = rast(rasters[i])
      
      #calculate and save softmax
      
      s_l = lapply(1:nlyr(scene), function(i){ #scale the raw bands for use with softmax
        b = scene[[i]]
        nom = names(scene)[i]
        if(nom %in% bands){
          b = b/sf
        } else {
          b 
        }
        return(b)
      })
      scene_s = rast(s_l)
      
      scene_sm = softmax(scene_s, append_name = T)
      terra::writeRaster(scene_sm, filename, overwrite = T)
    }
  })
  
  stopCluster(cl)
  plan('sequential')
  
  #----Calculate difference between time t and time t-1 for whole time series----
  
  #----z_score delta
  
  #create output dir
  dz_dir = paste0(z_dir,'_delta')
  dir.check(dz_dir)
  
  lapply(1:nrow(blocks_p), function(i){ #create subdirectories for every block
    d = paste0(dz_dir,'/',blocks_p$BLOCKNUM[i])
    dir.check(d)
  })
  
  indir = z_dir #input files
  
  in_subdirs = list.dirs(indir)
  in_subdirs = in_subdirs[in_subdirs != indir]
  
  #calculate difference rasters
  lapply(in_subdirs, function(d){
    #load input block subdirectory
    d_id = basename(d)
    files = list.files(d, full.names = T, pattern = '\\.tif$')
    df = tibble(
      files
    ) %>%
      #get dates from filepaths
      mutate(year_month = str_extract(files, "\\d{4}_\\d{2}")) %>%
      mutate(year_month_day = paste0(year_month,'_01')) %>%
      mutate(date = as.Date(str_replace_all(year_month_day,'_','-'))) %>%
      #sort filepaths by increasing order of date
      arrange(date)
    
    files_sorted = df$files
    
    #calculate rasters from files in subdirectory
    pblapply(2:length(files_sorted), function(i){ #start at 
      
      f_id = file_path_sans_ext(basename(files[i]))
      filename = paste0(dz_dir,'/',d_id,'/',f_id,'_delta.tif')
      if(!file.exists(filename)){
        rd = rast(files_sorted[i]) - rast(files_sorted[i-1])
        terra::writeRaster(rd, filename = filename)
      }
    })
  })
  
  #----not normalized delta
  
  #create output dir
  dn_dir = paste0(nonorm_dir,'_delta')
  dir.check(dn_dir)
  
  lapply(1:nrow(blocks_p), function(i){ #create subdirectories for every block
    d = paste0(dn_dir,'/',blocks_p$BLOCKNUM[i])
    dir.check(d)
  })
  
  indir = nonorm_dir #input files
  
  in_subdirs = list.dirs(indir)
  in_subdirs = in_subdirs[in_subdirs != indir]
  
  #calculate difference rasters
  lapply(in_subdirs, function(d){
    #load input block subdirectory
    d_id = basename(d)
    files = list.files(d, full.names = T, pattern = '\\.tif$')
    df = tibble(
      files
    ) %>%
      #get dates from filepaths
      mutate(year_month = str_extract(files, "\\d{4}_\\d{2}")) %>%
      mutate(year_month_day = paste0(year_month,'_01')) %>%
      mutate(date = as.Date(str_replace_all(year_month_day,'_','-'))) %>%
      #sort filepaths by increasing order of date
      arrange(date)
    
    files_sorted = df$files
    
    #calculate rasters from files in subdirectory
    pblapply(2:length(files_sorted), function(i){ #start at 
      
      f_id = file_path_sans_ext(basename(files[i]))
      filename = paste0(dn_dir,'/',d_id,'/',f_id,'_delta.tif')
      if(!file.exists(filename)){
        rd = rast(files_sorted[i]) - rast(files_sorted[i-1])
        terra::writeRaster(rd, filename = filename)
      }
    })
  })
  
  #----softmax delta
  
  #create output dir
  dsm_dir = paste0(sm_dir,'_delta')
  dir.check(dsm_dir)
  
  lapply(1:nrow(blocks_p), function(i){ #create subdirectories for every block
    d = paste0(dsm_dir,'/',blocks_p$BLOCKNUM[i])
    dir.check(d)
  })
  
  indir = sm_dir #input files
  
  in_subdirs = list.dirs(indir)
  in_subdirs = in_subdirs[in_subdirs != indir]
  
  #calculate difference rasters
  lapply(in_subdirs, function(d){
    #load input block subdirectory
    d_id = basename(d)
    files = list.files(d, full.names = T, pattern = '\\.tif$')
    df = tibble(
      files
    ) %>%
      #get dates from filepaths
      mutate(year_month = str_extract(files, "\\d{4}_\\d{2}")) %>%
      mutate(year_month_day = paste0(year_month,'_01')) %>%
      mutate(date = as.Date(str_replace_all(year_month_day,'_','-'))) %>%
      #sort filepaths by increasing order of date
      arrange(date)
    
    files_sorted = df$files
    
    #calculate rasters from files in subdirectory
    pblapply(2:length(files_sorted), function(i){ #start at 
      
      f_id = file_path_sans_ext(basename(files[i]))
      filename = paste0(dsm_dir,'/',d_id,'/',f_id,'_delta.tif')
      if(!file.exists(filename)){
        rd = rast(files_sorted[i]) - rast(files_sorted[i-1])
        terra::writeRaster(rd, filename = filename)
      }
    })
  })
  
  #----robust z_score delta
  
  #create output dir
  dzr_dir = paste0(zr_dir,'_delta')
  dir.check(dzr_dir)
  
  lapply(1:nrow(blocks_p), function(i){ #create subdirectories for every block
    d = paste0(dzr_dir,'/',blocks_p$BLOCKNUM[i])
    dir.check(d)
  })
  
  indir = zr_dir #input files
  
  in_subdirs = list.dirs(indir)
  in_subdirs = in_subdirs[in_subdirs != indir]
  
  #calculate difference rasters
  lapply(in_subdirs, function(d){
    #load input block subdirectory
    d_id = basename(d)
    files = list.files(d, full.names = T, pattern = '\\.tif$')
    df = tibble(
      files
    ) %>%
      #get dates from filepaths
      mutate(year_month = str_extract(files, "\\d{4}_\\d{2}")) %>%
      mutate(year_month_day = paste0(year_month,'_01')) %>%
      mutate(date = as.Date(str_replace_all(year_month_day,'_','-'))) %>%
      #sort filepaths by increasing order of date
      arrange(date)
    
    files_sorted = df$files
    
    #calculate rasters from files in subdirectory
    pblapply(2:length(files_sorted), function(i){ #start at 
      
      f_id = file_path_sans_ext(basename(files[i]))
      filename = paste0(dzr_dir,'/',d_id,'/',f_id,'_delta.tif')
      if(!file.exists(filename)){
        rd = rast(files_sorted[i]) - rast(files_sorted[i-1])
        terra::writeRaster(rd, filename = filename)
      }
    })
  })
  
  
  
}

#----timeseries analysis by block
{
  #----set up for processing----
  
  #get list of all files
  dirs = c(
    #change data
    dz_dir,
    dzr_dir,
    dsm_dir,
    dn_dir
    #no-change data
    ,z_dir
    ,zr_dir
    ,sm_dir
    ,nonorm_dir
  )
  
  all_files <- unlist(sapply(dirs, list.files, pattern = "\\.tif$", full.names = TRUE, recursive = TRUE))
  
  #get vector of all possible dates in timeseries
  start_date <- as.Date("2021-01-01")
  end_date <- as.Date("2024-12-31")
  date_vector <- seq.Date(from = start_date, to = end_date, by = "month")
  date_vector <- format(date_vector, "%Y-%m")
  date_vector = as.character(date_vector)
  date_vector = str_replace_all(date_vector,'-','_')
  
  #all possible datasets (strings to look for in filepaths to tell what dataset a raster belongs to)
  dataset_strings <- str_extract(dirs, "BlockClipped.*") %>% paste0('/') #add '/' to avoid confusion where some datasets are substrings of others
  
  #block ids
  block_ids = blocks_p$BLOCKNUM
  
  #----Turn each raster into a dataframe with global statistics----
  
  data_string = 'Global_Stats_nomask'
  data_filename = paste0(bm_dir,'/',data_string,'.csv')
  
  if(!file.exists(data_filename)){
    
    # plan('multisession', workers = 8)
    
    df_l = pblapply(all_files, function(x){
      
      # print(paste('Processing', x)) #uncomment to identify files where processing fails
      
      r = rast(x)
      d = global(r, fun = c("max", "min", "mean", "sum", "range", "rms", "sd", "std", "isNA", "notNA"), na.rm=T)
      
      #get dataset
      dataset = find_substring(x,dataset_strings)
      d[['dataset']] = dataset
      
      #get block
      block = find_substring(x,block_ids)
      d[['block_id']] = block
      
      #get date
      date_ = find_substring(x, date_vector)
      d[['acquisition_date']] = date_
      
      #filepath
      d[['path']] = x
      
      #rownames
      d[['v1']] = rownames(d)
      
      #delta
      d[['delta']] = str_detect(x, '_delta')
      
      return(d)
    })
    
    results_df = bind_rows(df_l)
    
    write.csv(results_df, data_filename)
    
    # plan('sequential')
  }
  
  #----process combined data----
  results_df = read_csv(data_filename) %>%
    #reformat dates
    mutate(
      year = str_remove(acquisition_date, '_.*'), 
      month = str_remove(acquisition_date, '.*_')
    ) %>%
    mutate(acquisition_date2 = as.Date(paste0(year,'-',month,'-01'))
    ) %>%
    mutate(julian_day = yday(acquisition_date2)
    ) %>%
    #modify variable names
    mutate(var_name = str_remove(v1, '_.*')) %>%
    # mutate(var_name = str_remove(var_name, "\\.\\.\\..*")) %>%
    filter(var_name != 'max_DN') %>%
    # fix dataset names
    mutate(block_type = ifelse(str_detect(block_id, 'NoChange'), 'Control', 'Thinning')) %>%
    mutate(dataset = str_replace(dataset, 'BlockClipped','')) %>%
    mutate(dataset = ifelse(str_detect(dataset, 'Z|SM'), dataset, paste0('Non-normalized', dataset))) %>%
    mutate(dataset = str_remove(dataset, "^_")) %>%
    mutate(dataset = str_remove(dataset, '/$'))
  
  
  #----plotting vegetation indices over time----
  
  {
    indices = c(
      # "blue", "green", "red"
      # ,
      "Hue"
      # ,
      # "GCC"
      ,
      "NDGR"
      # ,
      # "BI",        
      # "CI", "CRI550",   
      # "GLI", "TVI"
      # , "VARIgreen"
    ) #vector of indices to look at
    
    blocks_ = c(
      # "12L_C5",
      "12L_D345"
      # ,  "12L_C4",    "12L_C7",    "12L_B8C3",  "12N_T3",    "12N_1X1W",  "12N_V1"
      ,"NoChange1.1_conifer", "NoChange1.2_conifer",
      "NoChange4_conifer",   "NoChange5_conifer",  
      "NoChange6_conifer",   "NoChange7_conifer")
    
    datasets_ = c(
      'Z','Z_delta','Zrobust', 'Zrobust_delta'
    )
    
    months_ = c(
      '01','02','03',
      '04',
      '05','06','07','08','09','10','11'
      ,'12'
    )
    
    harvest_dates_df = tibble(
      block = blocks_[!str_detect(blocks_, 'NoChange')],
      harvest_month = c('2023-02'
                        ,'2022-07'
                        ,'2023-02'
                        ,'2023-02'
                        ,'2023-02'
                        ,'2023-02'
                        ,'2024-03'
                        ,'2024-03')
      ,note = c('-','-','-','-','-','partial, complete 2023-03','-','partial, complete 2024-04')
    )
    
    
    
    subset_df = results_df %>% 
      filter(var_name %in% indices
             , str_detect(block_id, paste(blocks_, collapse = "|"))
             , dataset %in% datasets_
             , month %in% months_)
    
    ggplot(subset_df, aes(x = acquisition_date2, y = mean)) +
      geom_point(aes(color = block_type))+
      geom_line(aes(color = block_id))+
      # geom_ribbon(aes(x = acquisition_date, ymin = mean-sds, ymax = mean+sds))+
      scale_x_date(
        date_breaks = "1 year",       # Major tick marks every year
        date_labels = "%Y %m",           # Year labels on major ticks
        date_minor_breaks = "1 month" # Minor tick marks every month
      )+
      facet_grid(rows = vars(dataset), cols = vars(var_name), scales = 'free') +
      geom_vline(xintercept = as.Date('2022-07-01')) +
      ggtitle('Change from previous image')
  }
  
}

#----timeseries analysis by thinned vs not thinned pixels ----
{

  harvest_threshold = -1 #the drop in height (in m) for a pixel to be considered "harvested" 
  
  data_string = paste0('GlobalStats_byblock_ThinnedVsNotVsTotal_HarvestThreshold=',harvest_threshold,'m')
  data_filename = paste0(bm_dir,'/',data_string,'.csv')
  
  if(!file.exists(data_filename)){
    
    #----load LiDAR data----
    lidar_dir = 'data/Quesnel_thinning/chm_change'
    lidar_files = list.files(lidar_dir, full.names = T, recursive = T)
    
    lidar_ids = basename(file_path_sans_ext(lidar_files))
    lidar_l = pblapply(lidar_files, rast)
    lidar_l = pblapply(lidar_l, function(x)if(crs(x)!=crs(blocks_p)){project(x, crs(blocks_p))})
    
    
    #----create harvest mask----
    
    harvest_masks = pblapply(lidar_l, function(x)ifel(x <= harvest_threshold, 1, 0)) #1==thinned, 0==not thinned
    names(harvest_masks) = lidar_ids
    # list_plot(harvest_masks)
    
    #----crop masks to harvest blocks----
    
    harvest_masks_clipped = pblapply(1:length(harvest_masks), function(i){
      block = blocks_p[blocks_p$BLOCKNUM == lidar_ids[i],]
      hm = crop(harvest_masks[[i]], block,mask = T)
    })
    names(harvest_masks_clipped) = lidar_ids
    # list_plot(harvest_masks_clipped)
    
    #----mask PS rasters, save data in dataframe----
    
    #remove NoChange rasters from list of all files
    all_files_thinning = all_files[!str_detect(all_files, 'NoChange')]
    
    # plan('multisession', workers = 8)
    
    
    df_l = pblapply(all_files_thinning, function(x){
      
      # print(paste('Processing', x)) #uncomment to identify files where processing fails
      
      #get raster
      r = rast(x)
      
      #select appropriate mask layer
      block = find_substring(x, lidar_ids)
      m = harvest_masks_clipped[[block]]
      
      #resample mask to match raster
      m_rs = resample(m, r, method = 'mode')
      
      #calculate global stats for thinning pixels
      m_rs1 = ifel(m_rs == 1, 1,NA) #make mask for thinning
      r_m = mask(r, m_rs1)
      d1 = global(r_m, fun = c("max", "min", "mean", "sum", "range", "rms", "sd", "std", "isNA", "notNA"), na.rm=T)
      d1[['block_pixel_stratum']] = 'Thinned'
      
      #calculate global stats for thinning pixels
      m_rs2 = ifel(m_rs == 0, 1,NA) #make mask for thinning
      r_m = mask(r, m_rs2)
      d2 = global(r_m, fun = c("max", "min", "mean", "sum", "range", "rms", "sd", "std", "isNA", "notNA"), na.rm=T)
      d2[['block_pixel_stratum']] = 'Not_thinned'
      
      #total block stats
      d3 = global(r, fun = c("max", "min", "mean", "sum", "range", "rms", "sd", "std", "isNA", "notNA"), na.rm=T)
      d3[['block_pixel_stratum']] = 'Total'
      
      #combine thinning, non-thinning, and total stats
      d = rbind(d1,d2,d3)
      
      #get dataset
      dataset = find_substring(x,dataset_strings)
      d[['dataset']] = dataset
      
      #get block
      d[['block_id']] = block
      
      #get date
      date_ = find_substring(x, date_vector)
      d[['acquisition_date']] = date_
      
      #filepath
      d[['file_path']] = x
      
      #rownames
      d[['v1']] = rownames(d)
      
      #delta
      d[['delta']] = str_detect(x, '_delta')
      
      return(d)
    })
    
    results_df = bind_rows(df_l)
    
    write.csv(results_df, data_filename)
    
    # plan('sequential')
  }
  
  data_string = paste0('GlobalStats_byblock_ThinnedVsNotVsTotal_HarvestThreshold=',harvest_threshold,'m')
  data_filename = paste0(bm_dir,'/',data_string,'.csv')
  
  if(!file.exists(data_filename)){
    
    #----load LiDAR data----
    lidar_dir = 'data/Quesnel_thinning/chm_change'
    lidar_files = list.files(lidar_dir, full.names = T, recursive = T)
    
    lidar_ids = basename(file_path_sans_ext(lidar_files))
    lidar_l = pblapply(lidar_files, rast)
    lidar_l = pblapply(lidar_l, function(x)if(crs(x)!=crs(blocks_p)){project(x, crs(blocks_p))})
    
    
    #----create harvest mask----
    
    harvest_masks = pblapply(lidar_l, function(x)ifel(x <= harvest_threshold, 1, 0)) #1==thinned, 0==not thinned
    names(harvest_masks) = lidar_ids
    # list_plot(harvest_masks)
    
    #----crop masks to harvest blocks----
    
    harvest_masks_clipped = pblapply(1:length(harvest_masks), function(i){
      block = blocks_p[blocks_p$BLOCKNUM == lidar_ids[i],]
      hm = crop(harvest_masks[[i]], block,mask = T)
    })
    names(harvest_masks_clipped) = lidar_ids
    # list_plot(harvest_masks_clipped)
    
    #----mask PS rasters, save data in dataframe----
    
    #remove NoChange rasters from list of all files
    all_files_thinning = all_files[!str_detect(all_files, 'NoChange')]
    
    # plan('multisession', workers = 8)
    
    
    df_l = pblapply(all_files_thinning, function(x){
      
      print(paste('Processing', x)) #uncomment to identify files where processing fails
      
      #get raster
      r = rast(x)
      
      #select appropriate mask layer
      block = find_substring(x, lidar_ids)
      m = harvest_masks_clipped[[block]]
      
      #resample mask to match raster
      m_rs = resample(m, r, method = 'mode')
      
      #calculate global stats for thinning pixels
      m_rs1 = ifel(m_rs == 1, 1,NA) #make mask for thinning
      r_m = mask(r, m_rs1)
      d1 = global(r_m, fun = c("max", "min", "mean", "sum", "range", "rms", "sd", "std", "isNA", "notNA"), na.rm=T)
      d1[['block_pixel_stratum']] = 'Thinned'
      
      #calculate global stats for thinning pixels
      m_rs2 = ifel(m_rs == 0, 1,NA) #make mask for thinning
      r_m = mask(r, m_rs2)
      d2 = global(r_m, fun = c("max", "min", "mean", "sum", "range", "rms", "sd", "std", "isNA", "notNA"), na.rm=T)
      d2[['block_pixel_stratum']] = 'Not_thinned'
      
      #total block stats
      d3 = global(r, fun = c("max", "min", "mean", "sum", "range", "rms", "sd", "std", "isNA", "notNA"), na.rm=T)
      d3[['block_pixel_stratum']] = 'Total'
      
      #combine thinning, non-thinning, and total stats
      d = rbind(d1,d2,d3)
      
      #get dataset
      dataset = find_substring(x,dataset_strings)
      d[['dataset']] = dataset
      
      #get block
      d[['block_id']] = block
      
      #get date
      date_ = find_substring(x, date_vector)
      d[['acquisition_date']] = date_
      
      #filepath
      d[['file_path']] = x
      
      #rownames
      d[['v1']] = rownames(d)
      
      #delta
      d[['delta']] = str_detect(x, '_delta')
      
      return(d)
    })
    
    results_df = bind_rows(df_l)
    
    write.csv(results_df, data_filename)
    
    # plan('sequential')
  }
  
  #----process combined data----
  results_df = read_csv(data_filename) %>%
    #reformat dates
    mutate(
      year = str_remove(acquisition_date, '_.*'), 
      month = str_remove(acquisition_date, '.*_')
    ) %>%
    mutate(acquisition_date2 = as.Date(paste0(year,'-',month,'-01'))
    ) %>%
    mutate(julian_day = yday(acquisition_date2)
    ) %>%
    #modify variable names
    mutate(var_name = str_remove(v1, '_.*')) %>%
    # mutate(var_name = str_remove(var_name, "\\.\\.\\..*")) %>%
    filter(var_name != 'max_DN') %>%
    # fix dataset names
    mutate(block_type = ifelse(str_detect(block_id, 'NoChange'), 'Control', 'Thinning')) %>%
    mutate(dataset = str_replace(dataset, 'BlockClipped','')) %>%
    mutate(dataset = ifelse(str_detect(dataset, 'Z|SM'), dataset, paste0('Non-normalized', dataset))) %>%
    mutate(dataset = str_remove(dataset, "^_")) %>%
    mutate(dataset = str_remove(dataset, '/$'))
  
  
  
  #----plotting vegetation indices over time----
  
  {
    indices = c(
      # "blue", "green", "red"
      # ,
      # "Hue"
      # ,
      "GCC"
      # ,
      # "NDGR"
      # ,
      # "BI"
      # ,        
      # "CI"
      # , 
      # "CRI550"
      # ,   
      # "GLI"
      # , 
      # "TVI"
      # , 
      # "VARIgreen"
    ) #vector of indices to look at
    
    blocks_ = c(
      "12L_C5"
      ,
      "12L_D345"
      ,  "12L_C4",    "12L_C7",    "12L_B8C3",  "12N_T3",    "12N_1X1W",  "12N_V1"
      # ,"NoChange1.1_conifer", "NoChange1.2_conifer",
      # "NoChange4_conifer",   "NoChange5_conifer",  
      # "NoChange6_conifer",   "NoChange7_conifer"
      )
    
    datasets_ = c(
      'Z'
      # ,
      # 'Z_delta',
      # 'Zrobust'
      # , 'Zrobust_delta'
    )
    
    months_ = c(
      '01','02','03',
      '04',
      '05','06','07','08','09','10','11'
      ,'12'
    )
    
    harvest_dates_df = tibble(
      block_id = blocks_p$BLOCKNUM[!str_detect(blocks_p$BLOCKNUM, 'NoChange')],
      harvest_month = c('2023-02'
                        ,'2022-07'
                        ,'2023-02'
                        ,'2023-02'
                        ,'2023-02'
                        ,'2023-02'
                        ,'2024-03'
                        ,'2024-03')
      ,harvest_note = c('-','-','-','-','-','partial, complete 2023-03','-','partial, complete 2024-04')
    )
    
    
    
    subset_df = results_df %>% 
      filter(var_name %in% indices
             , str_detect(block_id, paste(blocks_, collapse = "|"))
             , dataset %in% datasets_
             , month %in% months_) %>%
      left_join(harvest_dates_df %>% mutate(harvest_date = as.Date(paste0(harvest_month,'-01'))))
    
    ggplot(subset_df, aes(x = acquisition_date2, y = mean)) +
      geom_point(aes(color = block_pixel_stratum))+
      geom_line(aes(color = block_pixel_stratum))+
      geom_ribbon(aes(ymin = mean - sd, ymax = mean + sd, fill = block_pixel_stratum), alpha = 0.2, color = NA) +
      # geom_ribbon(aes(x = acquisition_date, ymin = mean-sds, ymax = mean+sds))+
      scale_x_date(
        date_breaks = "1 year",       # Major tick marks every year
        date_labels = "%Y %m",           # Year labels on major ticks
        date_minor_breaks = "1 month" # Minor tick marks every month
      )+
      facet_grid(rows = vars(block_id), cols = vars(dataset), scales = 'free') +
      # geom_vline(xintercept = as.Date('2022-07-01')) + #12L_D345
      geom_vline(aes(xintercept = harvest_date, linetype = 'Harvest date')) + #12L_C4, 12L_C5
      ggtitle(unique(subset_df$var_name))
  }
}








#save raw pixel values