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
#extra functions
source('helper_functions.R')
#functions for checking or visualizing radiometric consistency
source('https://raw.githubusercontent.com/Spencer-Shields/planetscope_radiometric_correction/refs/heads/main/check_radiometric_consistency.R')

#set up number of cores to be used for parallel processing
cores = detectCores()

# cl = makeCluster(ceiling(cores/4))
# plan('multisession', workers = 8)

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
  
  scale_factor = 255 #scale factor to use for calculating softmax (i.e. maximum possible band pixel value)
  
  bands = c('blue', 'green', 'red')
  
  # cl = makeCluster(ceiling(cores/4))
  plan('multisession', workers = 12)
  
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
      # names(scene_z) = paste0(names(scene),'_',z_dir_string)
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
      # names(scene_z) = paste0(names(scene),'_',zr_dir_string)
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
          b = b/scale_factor
        } else {
          b 
        }
        return(b)
      })
      scene_s = rast(s_l)
      
      scene_sm = softmax(scene_s, append_name = F)
      terra::writeRaster(scene_sm, filename, overwrite = T)
    }
  })
  
  # stopCluster(cl)
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
    
    results1_df = bind_rows(df_l)
    
    write.csv(results1_df, data_filename)
    
    # plan('sequential')
  }
  
  #----process combined data----
  results1_df = read_csv(data_filename) %>%
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
    
    
    
    subset_df = results1_df %>% 
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

#----timeseries analysis by block and by thinned vs not thinned pixels
{
  harvest_threshold = -5.2 #the drop in height (in m) for a pixel to be considered "harvested". -5.2 = average Otsu threshold based on 0.25m raster
  
  data_string = paste0('GlobalStats_byblock_ThinnedVsNotVsTotal_HarvestThreshold=',harvest_threshold,'m')
  data_filename = paste0(bm_dir,'/',data_string,'.csv')
  
  #run if summary data file does not exist
  if(!file.exists(data_filename)){
    
    #----load LiDAR data----
    lidar_dir = 'data/Quesnel_thinning/chm_change_basemap'
    lidar_files = list.files(lidar_dir, full.names = T, recursive = T)
    lidar_ids = basename(file_path_sans_ext(lidar_files))
    names(lidar_files) = lidar_ids
    
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
    
    #----mask PS rasters, save global values in dataframe----
    
    #remove NoChange rasters from list of all files
    all_files_thinning = all_files[!str_detect(all_files, 'NoChange')]
    
    plan('multisession', workers = 10)
    
    df_l = pblapply(all_files_thinning, function(x){
      
      # print(paste('Processing', x)) #uncomment to identify files where processing fails
      
      #get raster
      r = rast(x)
      
      #select appropriate lidar layer, create binary thinning mask
      block = find_substring(x, lidar_ids)
      lid_file = lidar_files[block]
      lid_r = rast(lid_file)
      m = ifel(lid_r <= harvest_threshold, 1, 0)
      
      #calculate global stats for thinning pixels
      m1 = ifel(m == 1, 1,NA) #make mask for thinning
      r_m = mask(r, m1)
      d1 = global(r_m, fun = c("max", "min", "mean", "sum", "range", "rms", "sd", "std", "isNA", "notNA"), na.rm=T)
      d1[['block_pixel_stratum']] = 'Thinned'
      d1[['v1']] = rownames(d1)
      
      #calculate global stats for thinning pixels
      m2 = ifel(m == 0, 1,NA) #make mask for thinning
      r_m = mask(r, m2)
      d2 = global(r_m, fun = c("max", "min", "mean", "sum", "range", "rms", "sd", "std", "isNA", "notNA"), na.rm=T)
      d2[['block_pixel_stratum']] = 'Not_thinned'
      d2[['v1']] = rownames(d2)
      
      #total block stats
      d3 = global(r, fun = c("max", "min", "mean", "sum", "range", "rms", "sd", "std", "isNA", "notNA"), na.rm=T)
      d3[['block_pixel_stratum']] = 'Total'
      d3[['v1']] = rownames(d3)
      
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
      
      # #rownames
      # d[['v1']] = rownames(d)
      
      #delta
      d[['delta']] = str_detect(x, '_delta')
      
      return(d)
    }
    ,cl = 'future'
    )
    
    results_df = bind_rows(df_l)
    
    write.csv(results_df, data_filename)
    
    plan('sequential')
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
  
  
  
  #----check results----
  results_summ = results_df %>%
    group_by(dataset, block_pixel_stratum, block_id, var_name) %>%
    summarise(count = n()) %>%
    ungroup()%>%
    complete(dataset, block_pixel_stratum, block_id, var_name, fill = list(count = 0))
  # view(results_summ)
  
  #----rescale variable values so that they can be plotted together----
  
  #get min and max means and standard deviations
  vars_summ = results_df %>%
    group_by(var_name) %>%
    summarise(max_mean = max(mean),
              min_mean = min(mean),
              mean_mean = mean(mean),
              max_sd = max(std),
              min_sd = min(std),
              mean_sd = mean(std))
  
  #----plotting vegetation indices over time----
  {
    indices = c(
      # "blue", "green", "red"
      # ,
      "Hue"
      ,
      "GCC"
      ,
      "NDGR"
      ,
      "BI"
      # ,
      # "CI", "CRI550",
      # "GLI", "TVI"
      # , "VARIgreen"
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
      # 'Z_delta'
      # ,
      # 'Zrobust'
      # , 'Zrobust_delta'
      # ,
      # 'SM'
      # ,
      # 'SM_delta'
      # 'Non-normalized'
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
    
    pixel_stratum = c('Not_thinned', 'Thinned')
    
    subset_df = results_df %>% 
      filter(var_name %in% indices
             , str_detect(block_id, paste(blocks_, collapse = "|"))
             , dataset %in% datasets_
             , month %in% months_
             , block_pixel_stratum %in% pixel_stratum) %>%
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
      facet_grid(rows = vars(block_id), cols = vars(var_name), scales = 'free') +
      # geom_vline(xintercept = as.Date('2022-07-01')) + #12L_D345
      geom_vline(aes(xintercept = harvest_date, linetype = 'Harvest date')) + #get harvest date for each plot
      ggtitle(unique(subset_df$dataset))
  }
}

#----timeseries analysis using raw pixel values
{
  data_string = 'RawBasemapValues_withCHM_cleaned.parquet'
  data_filename = paste0(bm_dir,'/',data_string)
  if(!file.exists(data_filename)){
    #----set up for processing----
    #load LiDAR data
    lidar_dir = 'data/Quesnel_thinning/chm_change_basemap'
    lidar_files = list.files(lidar_dir, full.names = T, recursive = T)
    lidar_ids = basename(file_path_sans_ext(lidar_files))
    names(lidar_files) = lidar_ids
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
    
    #----combine planetscope raster stacks with canopy height, save raw pixel values in giant data table----
    
    #remove NoChange rasters from list of all files
    all_files_thinning = all_files[!str_detect(all_files, 'NoChange')]
    
    plan('multisession', workers = 10)
    
    df_l = pblapply(all_files_thinning, function(x){
      
      # print(paste('Processing', x)) #uncomment to identify files where processing fails
      
      #get raster
      r = rast(x)
      
      #select appropriate lidar layer, create binary thinning mask
      block = find_substring(x, lidar_ids)
      lid_file = lidar_files[block]
      lid_r = rast(lid_file)
      names(lid_r) = 'CHM_change'
      
      #append chm to planetscope stack
      r_c = c(r, lid_r)
      
      #extract values, convert to data.table
      v = values(r_c, dataframe=T)
      d = setDT(v)
      
      # Remove 'max_DN' column (in-place)
      d[, max_DN := NULL]
      
      # Assign 'id' to each pixel (efficient with data.table)
      d[, id := .I]
      
      # Get 'dataset' with string manipulations in one step
      d[, dataset := find_substring(x, dataset_strings)]
      d[, dataset := ifelse(str_detect(dataset, 'SM|Z'), dataset, paste0('Non-normalized', dataset))]
      d[, dataset := str_replace(dataset, 'BlockClipped', '')]
      d[, dataset := str_remove(dataset, "^_")]
      d[, dataset := str_remove(dataset, '/$')]
      
      # Assign 'block_id' directly (assuming 'block' is a pre-defined vector)
      d[, block_id := block]
      
      # Get 'acquisition_date' with string manipulations and conversion to Date
      d[, acquisition_date := find_substring(x, date_vector)]
      d[, acquisition_date := str_replace(acquisition_date, '_', '-')]
      d[, acquisition_date := paste0(acquisition_date, '-01')]
      d[, acquisition_date := as.Date(acquisition_date)]
      
      # Calculate Julian day and assign it
      d[, julian_day := yday(acquisition_date)]
      
      # d = d |>
      #   pivot_longer(cols = ps_feats, names_to = 'var_name', values_to = 'value')
      
      return(d)
    }
    ,cl = 'future'
    )
    
    results2_df = rbindlist(df_l)
    
    write_parquet(results2_df, data_filename)
    
    plan('sequential')
  }
  
  #----separate into multiple parquet files based on dataset to speed up processing----
  rawvals_dir = paste0(bm_dir,'/','RawValues_parquets_cleaned')
  if(!dir.exists(rawvals_dir)){
    dir.check(rawvals_dir)
    results2_df = open_dataset(data_filename) |>
      group_by(dataset) |>
      write_dataset(path = rawvals_dir, format = 'parquet')
  }
  results2_df = open_dataset(rawvals_dir)
  
  
  
  #----test if there's a significant difference between thinning and not thinning for each block and each date----
  
  
  all_cols = names(results2_df)
  non_ps_cols = c('block_id', 'dataset', 'id', 'acquisition_date', 'julian_day', 'CHM_change')
  ps_feats = all_cols[!all_cols %in% non_ps_cols]
  ps_feats = ps_feats[ps_feats != 'VARIgreen'] #remove this index since it is undefined for many dates
  
  acquisition_date <- results2_df %>%
    select(acquisition_date) %>%
    collect() %>%
    unique() %>%
    pull(acquisition_date)  # Extract as a vector
  
  dataset <- results2_df %>%
    select(dataset) %>%
    collect() %>%
    unique() %>%
    pull(dataset)
  
  statistical_significance_summ = CJ(
    var_name = ps_feats,
    acquisition_date = acquisition_date,
    dataset = dataset)
  
  harvest_threshold = -5.2 #mean Otsu's threshold applied to the change in CHM layers
  
  # plan('multisession', workers = 10)
  
  stat_results = pblapply(1:nrow(statistical_significance_summ), function(i){
  # stat_results = pblapply(1:10, function(i){
    
    # ad_thinned_p = c()
    # ad_notthinned_p = c()
    # levene_p = c()
    # test_type = c()
    # difference_p = c()
    
    #initialize list where results will go
    results_list = list(NA, NA, NA, NA, NA)
    
    print(paste0('Processing ',i,'/',nrow(statistical_significance_summ)))
    
    sub_df = results2_df |>
      filter(dataset == statistical_significance_summ$dataset[i],
             acquisition_date == statistical_significance_summ$acquisition_date[i]) |>
      mutate(thinned = case_when(
        CHM_change <= harvest_threshold ~ 'Thinned',
        TRUE ~ 'Not_thinned'
      )) |>
      select(statistical_significance_summ$var_name[i], 'thinned') |>
      filter(!is.na(!!sym(statistical_significance_summ$var_name[i]))) |>
      collect()
    
    #proceed if there are sufficient data (since not ever combinaiton of block, acquisition date, and varaible might exist)
    if(nrow(sub_df)>0){
      
      thinned = sub_df[[statistical_significance_summ$var_name[i]]][sub_df$thinned == 'Thinned']
      not_thinned = sub_df[[statistical_significance_summ$var_name[i]]][sub_df$thinned == 'Not_thinned']
      
      ad_t = ad.test(thinned)
      results_list[[1]] = ad_t[['p.value']]
      ad_nt = ad.test(not_thinned)
      results_list[[2]] = ad_nt[['p.value']]
      
      formula = as.formula(paste(statistical_significance_summ$var_name[i], '~ thinned'))
      
      sub_df$thinned = as.factor(sub_df$thinned)
      lev = leveneTest(formula, data = sub_df)
      results_list[[3]] = lev$`Pr(>F)`[1]
      
      # Test selection logic
      # Significance level
      alpha <- 0.05
      
      # Check normality and variance conditions
      is_normal_thinned <- results_list[[1]] > alpha
      is_normal_not_thinned <- results_list[[2]] > alpha
      is_variance_equal <- results_list[[3]] > alpha
      
      # Select appropriate test
      if (is_normal_thinned && is_normal_not_thinned) {
        if (is_variance_equal) {
          # Student's t-test (equal variances)
          t_test_result <- t.test(thinned, not_thinned, var.equal = TRUE)
          results_list[[4]] <- "Student's t-test"
        } else {
          # Welch's t-test (unequal variances)
          t_test_result <- t.test(thinned, not_thinned, var.equal = FALSE)
          results_list[[4]] <- "Welch's t-test"
        }
        results_list[[5]] <- t_test_result$p.value
      } else {
        # Non-parametric test if normality assumption is violated
        wilcox_result <- wilcox.test(thinned, not_thinned)
        results_list[[4]] <- "Mann-Whitney U"
        results_list[[5]] <- wilcox_result$p.value
      }
      
    }
    
    # results_list = list(
    #   ad_thinned_p,
    #   ad_notthinned_p,
    #   levene_p,
    #   test_type,
    #   difference_p
    # )
    
    names(results_list) = c('Anderson-Darling thinned p-value',
                            'Anderson-Darling notthinned p-value',
                            'Levene p-value',
                            'Test type',
                            'Difference p-value')
    
    return(results_list)
    
  }
  )
  
  #add stats
  statistical_significance_summ[['anderson_darling_thinned_p']] = sapply(stat_results, function(x)x[[1]])
  statistical_significance_summ[['anderson_darling_notthinned_p']] = sapply(stat_results, function(x)x[[2]])
  statistical_significance_summ[['levene_p']] = sapply(stat_results, function(x)x[[3]])
  statistical_significance_summ[['test_type']] = sapply(stat_results, function(x)x[[4]])
  statistical_significance_summ[['difference_p']] = sapply(stat_results, function(x)x[[5]])
  
  #----save results----
  
  stat_sig_file = paste0(bm_dir, '/StatDiff_ThinnedvsNotthinned.csv')
  
  if(!file.exists(stat_sig_file)){
    write_csv(statistical_significance_summ, stat_sig_file)
  }
  
  statistical_significance_summ = fread(stat_sig_file)
  
  #filter out NAs
  statistical_significance_summ = statistical_significance_summ[!is.na(difference_p)]
}
