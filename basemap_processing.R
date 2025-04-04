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
library(moments)
library(mgcv)
library(fitdistrplus)
library(betareg)
#extra functions
source('helper_functions.R')
#functions for checking or visualizing radiometric consistency
source('https://raw.githubusercontent.com/Spencer-Shields/planetscope_radiometric_correction/refs/heads/main/check_radiometric_consistency.R')

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
  # plan('multisession', workers = 12)
  
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
  # plan('sequential')
  
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
  
  data_string = paste0('GlobalStats_byblock_ThinnedVsNotVsTotal_HarvestThreshold=',harvest_threshold,'m','_withTests_EqualClasses')
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
    
    #remove NoChange rasters from list of all files
    all_files_thinning = all_files[!str_detect(all_files, 'NoChange')]
    
    #----get spatial sample for each study area to make thinning/non-thinning classes same size----
    
    #load lidar CHM change layers
    lid_masks_equalclasses = pblapply(lidar_files, function(x){
      #load lidar raster
      # print(paste('Processing',x))
      r = rast(x)
      
      #create raster stack with thinning/nonthinning mask layers
      m_t = ifel(r <= harvest_threshold, 1, NA) #create thinning mask layer
      m_nt = ifel(r > harvest_threshold, 0, NA) #create non-thinning mask layer
      m = c(m_t, m_nt) #stack layers
      names(m) = c('Thinning', 'Non_thinning') #name layers
      
      #get number of non-NA pixels in each layer of the raster stack
      gvals = global(m, 'notNA')
      
      n_t = gvals$notNA[rownames(gvals)=='Thinning'] #number of non-NA thinning pixels
      n_nt = gvals$notNA[rownames(gvals)=='Non_thinning'] #number of non-NA non-thinning pixels
      
      #if one layer has more pixels than the other, do a subsample so that they have the same number and replace the raster layer with the subsampled layer
      if(n_t > n_nt){
        target_size = n_nt
        #sample thinning mask
        m_s_xy = terra::spatSample(m_t, size = target_size, xy=T, method='random', na.rm=T, replace=F, exhaustive = T)
        m_s = rast(m_s_xy, type = 'xyz')
        m[[1]] = m_s |> resample(m_t)
      }
      if(n_t < n_nt){
        target_size = n_t
        #sample thinning mask
        m_s_xy = terra::spatSample(m_nt, size = target_size, xy=T, method='random', na.rm=T, replace=F, exhaustive = T)
        m_s = rast(m_s_xy, type = 'xyz')
        # m_s = ifel(m_s == 1, 0,NA) #reset values to zero since they become one for some reason
        m[[2]] = m_s |> resample(m_nt)
      }
      
      names(m) = c('Thinning', 'Non-thinning')
      m = wrap(m) #wrap for parallel processing
      return(m)
    })
    
    #----mask PS rasters, save global values in dataframe----
    
    plan('multisession', workers = 10)
    
    # df_l = pblapply(1:10, function(i){
    df_l = pblapply(1:length(all_files_thinning), function(i){
      x = all_files_thinning[i]
      # print(paste0('Processing ', x,', ',i,'/',length(all_files_thinning))) #uncomment to identify files where processing fails
      
      #get raster
      r = rast(x)
      
      # #old code used when lidar is not sampled
      # #select appropriate lidar layer, create binary thinning mask
      # block = find_substring(x, lidar_ids)
      # lid_file = lidar_files[block]
      # lid_r = rast(lid_file)
      # m = ifel(lid_r <= harvest_threshold, 1, 0)
      # 
      # #make mask for thinning and non-thinng pixels
      # mt = ifel(m == 1, 1,NA) #make mask for thinning
      # r_t = mask(r, mt) #use mask to isolate thinning pixels in raster
      # 
      # mnt = ifel(m == 0, 1,NA) #make mask for non-thinning
      # r_nt = mask(r, mnt) #use mask to isolate non-thinning pixels in raster
      
      #get lidar mask
      block = find_substring(x, lidar_ids)
      lid_mask_wrapped = lid_masks_equalclasses[[block]]
      m = unwrap(lid_mask_wrapped)
      
      #mask planetscope raster using thinning and non-thinning masks
      r_t = mask(r, m[[1]])
      r_nt = mask(r, m[[2]])
      
      #calculate global stats for thinning pixels
      d1 = global(r_t, fun = c("max", "min", "mean", "sum", "range", "rms", "sd", "std", "isNA", "notNA"), na.rm=T)
      d1 = cbind(d1, 
                 global(r_t, median, na.rm=T) |> rename(median = global)
      )
      d1[['block_pixel_stratum']] = 'Thinned'
      d1[['v1']] = rownames(d1)
      
      #calculate global stats for non-thinning pixels
      d2 = global(r_nt, fun = c("max", "min", "mean", "sum", "range", "rms", "sd", "std", "isNA", "notNA"), na.rm=T)
      d2 = cbind(d2, 
                 global(r_nt, median, na.rm=T) |> rename(median = global)
      )
      d2[['block_pixel_stratum']] = 'Not_thinned'
      d2[['v1']] = rownames(d2)
      
      #total block stats
      d3 = global(r, fun = c("max", "min", "mean", "sum", "range", "rms", "sd", "std", "isNA", "notNA"), na.rm=T)
      d3 = cbind(d3, 
                 global(r, median, na.rm=T) |> rename(median = global)
      )
      d3[['block_pixel_stratum']] = 'Total'
      d3[['v1']] = rownames(d3)
      
      #combine thinning, non-thinning, and global stats
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
      
      #delta
      d[['delta']] = str_detect(x, '_delta')
      
      ###statistical measures between groups
      
      v_t = values(r_t, dataframe=T) |> mutate(thinning = 'Thinned')
      v_nt = values(r_nt, dataframe=T) |> mutate(thinning = 'Not_thinned') 
      v_df = rbind(v_t, v_nt)
      
      feats = unique(d$v1)
      feats = feats[!feats %in% c('max_DN', 'VARIgreen')] #don't run tests on max_DN
      
      tests_l = pblapply(feats, function(y){
        
        print(paste('Processing',y))
        #subset values to test from dataframe
        v = v_df |> select(all_of(c(y, 'thinning')))
        v = v[!is.na(v[[y]]),] #remove NA values
        v$thinning = as.factor(v$thinning)
        
        #initialize tibble to store results
        test_df = tibble(
          v1 = y,
          levene_p = NA,
          anderson_darling_p_thinning = NA,
          anderson_darling_p_notthinning = NA,
          anderson_darling_p_pooled = NA,
          student_t_p = NA,
          welch_t_p = NA,
          mann_whitney_p = NA,
          anova_p = NA,
          kruskal_wallis_p = NA,
          logistic_p = NA,
          logistic_aic = NA,
          logistic_AUC = NA,
          fisher_discriminant_ratio = NA,
          cohen_d = NA,
          jeffries_matusita_dist = NA,
          bhattacharya_dist = NA
        )
        
        if(nrow(v)>0){ #proceed if there are rows in the dataframe (i.e. values that aren't na)
          ##perform tests, store results in test_df
          
          #create formula and linear model
          formula = as.formula(paste(y,'~ thinning'))
          v_lm = lm(formula, v)
          
          #levene test
          levene = leveneTest(formula, v)
          test_df$levene_p = levene$`Pr(>F)`[1]
          
          #anderson-darling tests (if else statements handle cases when there is one unique value)
          t = v[[y]][v$thinning == 'Thinned']
          test_df$anderson_darling_p_thinning = if(length(unique(t)|length(t)<8) == 1) 0 else ad.test(t)$p.value
          nt = v[[y]][v$thinning == 'Not_thinned']
          test_df$anderson_darling_p_notthinning = if(length(unique(nt)|length(t)<8) == 1) 0 else ad.test(nt)$p.value
          total = v[[y]]
          test_df$anderson_darling_p_pooled = if(length(unique(total)|length(total)<8) == 1) 0 else ad.test(total)$p.value
          
          #Students t-test
          s_t = t.test(formula, v)
          test_df$student_t_p = s_t$p.value
          
          #Welch's t-test
          w_t = t.test(formula, v, var.equal = F)
          test_df$welch_t_p = w_t$p.value
          
          #Mann-Whitney U test
          mw = wilcox.test(formula, v)
          test_df$mann_whitney_p = mw$p.value
          
          #ANOVA
          anv = anova(v_lm)
          test_df$anova_p = anv$`Pr(>F)`[1]
          
          #Kruskal-Wallis
          kw = kruskal.test(formula, v)
          test_df$kruskal_wallis_p = kw$p.value
          
          #logistic regression model
          log_form = paste('thinning ~', y)
          log_m = glm(log_form, v, family = 'binomial')
          log_m_summ = summary(log_m)
          
          test_df$logistic_p = log_m_summ$coefficients[,4][[y]]
          test_df$logistic_aic = log_m_summ$aic
          
          log_p = v |> mutate(predictions = predict(log_m, v, type = 'response'))
          log_roc = roc(thinning ~ predictions, log_p, quiet=T)
          test_df$logistic_AUC = as.numeric(auc(log_roc))
          
          #Fisher discriminant ratio
          fdr = ((d1$mean[d1$v1==y]-d2$mean[d2$v1==y])^2)/(d1$sd[d1$v1==y]^2 + d2$sd[d2$v1==y]^2)
          test_df$fisher_discriminant_ratio = fdr
          
          #Cohen's d
          pooled_sd = sqrt(((d1$notNA[d1$v1 == y] - 1) * d1$sd[d1$v1 == y]^2 +
                              (d2$notNA[d2$v1 == y] - 1) * d2$sd[d2$v1 == y]^2) /
                             (d1$notNA[d1$v1 == y] + d2$notNA[d2$v1 == y] - 2))
          cohen_d = (d1$mean[d1$v1==y] - d2$mean[d2$v1==y])/pooled_sd
          test_df$cohen_d = cohen_d
          
          #Jeffries-Matusita distance
          jmd = JMdist(g = v$thinning, X = tibble(v[[y]]))
          test_df$jeffries_matusita_dist = jmd$jmdist
          
          #Bhattacharya distance
          bhat = BHATdist(g = v$thinning, X = tibble(v[[y]]))
          test_df$bhattacharya_dist = bhat$bhatdist
          
          #I 
        }
        
        #return results dataframe
        return(test_df)
      })
      test_results = bind_rows(tests_l)
      
      ###attach statistical test results to the summary stats dataframe
      
      d = d %>% left_join(test_results, by = 'v1')
      
      ####return final result
      return(d)
    }
    ,cl = 'future'
    )
    
    results_df = bind_rows(df_l)
    
    write_csv(results_df, data_filename)
    
    plan('sequential')
  }
  #----process combined data----
  {
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
    
    #identify significant values
    alpha = 0.05
    
    results_df = results_df |>
      mutate(mann_whitney_sig = ifelse(mann_whitney_p < alpha, 1, 0),
             logistic_sig = ifelse(logistic_p < alpha, 1, 0))
    #add harvest dates
    harvest_dates_df = tribble(
      ~block_id, ~harvest_month, ~note,
      '12L_C5', '2022-08', '-',
      '12L_D345','2022-07', '-',
      '12L_C4', '2022-11', '-',
      '12L_C7', '2022-09', '-',
      '12L_B8C3', '2022-09', '-',
      '12N_T3', '2023-02', 'partial, complete 2023-03',
      '12N_1X1W', '2024-03', '-',
      '12N_V1', '2024-03', 'partial, complete 2024-04'
    )
    results_df = results_df |>
      left_join(harvest_dates_df, by = 'block_id') |>
      mutate(harvest_date = as.Date(paste0(harvest_month, '-01')))
    #calculate time since harvest
    results_df = results_df |>
      mutate(months_since_harvest = as.numeric(acquisition_date2 - harvest_date)/30.44)
    
  }
  #----plotting vegetation indices over time----
  {
    indices = c(
      # "blue", "green", "red"
      # ,
      # "Hue"
      # ,
      # "GCC"
      # ,
      # "NDGR"
      # ,
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
    
    tests_ = c(
      'levene_p',
      'anderson_darling_p_thinning',
      'anderson_darling_p_notthinning',
      'anderson_darling_p_pooled',
      'student_t_p',
      'welch_t_p',
      'mann_whitney_p',
      'anova_p',
      'kruskal_wallis_p',
      'logistic_p',
      'logistic_aic',
      'logistic_AUC',
      'fisher_discriminant_ratio',
      'cohen_d',
      'jeffries_matusita_dist'
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
    
    # harvest_dates_df = tibble(
    #   block_id = blocks_p$BLOCKNUM[!str_detect(blocks_p$BLOCKNUM, 'NoChange')],
    #   harvest_month = c('2023-02'
    #                     ,'2022-07'
    #                     ,'2023-02'
    #                     ,'2023-02'
    #                     ,'2023-02'
    #                     ,'2023-02'
    #                     ,'2024-03'
    #                     ,'2024-03')
    #   ,harvest_note = c('-','-','-','-','-','partial, complete 2023-03','-','partial, complete 2024-04')
    # )
    
    pixel_stratum = c('Not_thinned', 'Thinned', 'Total')
    
    subset_df = results_df %>% 
      filter(var_name %in% indices
             , str_detect(block_id, paste(blocks_, collapse = "|"))
             , dataset %in% datasets_
             , month %in% months_
             , block_pixel_stratum %in% pixel_stratum)
    
    #harvest pixel values
    pixels_p = ggplot(subset_df, aes(x = acquisition_date2, y = mean)) +
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
      ggtitle(unique(subset_df$dataset))+
      theme_classic()
    
    #seperability test results
    log_p = ggplot(subset_df |> filter(block_pixel_stratum == 'Thinned') #remove redundant info (since thinning, not_thinning, and total have same stats info)
                   , aes(x = acquisition_date2, y = logistic_AUC)) +
      geom_point()+
      geom_line()+
      geom_point(data = subset_df|>
                   filter(logistic_sig == 1)|>
                   mutate(label = paste0('p<',alpha)), 
                 aes(x = acquisition_date2, y =1, color=label), shape = 8)+ #significance of model
      # scale_shape_manual(name = "Model Significance", values = c("Significant Model" = 8))+
      scale_x_date(
        date_breaks = "1 year",       # Major tick marks every year
        date_labels = "%Y %m",           # Year labels on major ticks
        date_minor_breaks = "1 month" # Minor tick marks every month
      )+
      facet_grid(rows = vars(block_id), cols = vars(var_name)
                 , scales = 'free'
      ) +
      geom_vline(aes(xintercept = harvest_date, linetype = 'Harvest date')) + #get harvest date for each plot
      ggtitle(unique(subset_df$dataset))+
      theme_classic()
    
    library(patchwork)
    pixels_p+log_p+plot_layout(guides = 'collect')
  }
  
  #----use GAM to assess seperability before and after harvesting----
  
  # #check what distribution best fits the data
  # samp = results_df |> filter(dataset == 'Z',
  #                             var_name == 'blue') 
  # fit_beta = fitdistrplus::fitdist(samp$logistic_AUC, 'beta') #beta is the choice since the range of possible values are [0-1]
  # fit_gamma = fitdistrplus::fitdist(samp$logistic_AUC, 'gamma')
  # 
  # ps_feats = unique(results_df$var_name)
  # ps_feats = ps_feats[!ps_feats %in% c('VARIgreen', 'max_DN')]
  # datasets = unique(results_df$dataset)
  # 
  # # plan('multisession', workers=8)
  # dataset_gam_l = pblapply(datasets, function(x){
  #   d = results_df |> filter(dataset == x)
  #   feat_gam_l = pblapply(ps_feats, function(y){
  #     d_ = d |> filter(var_name == y)
  #     d_$block_id = as.factor(d_$block_id)
  #     mod = gam(logistic_AUC ~ s(months_since_harvest) + s(block_id, bs = "re")
  #               ,data = d_, family = betar())
  #     return(mod)
  #   })
  #   names(feat_gam_l) = ps_feats
  #   return(feat_gam_l)
  # }
  # # ,cl = 'future'
  # )
  # # plan('sequential')
  # 
  # names(dataset_gam_l) = datasets
  # 
  # #compare with glms
  # plan('multisession', workers=8)
  # dataset_glm_l = pblapply(datasets, function(x){
  #   d = results_df |> filter(dataset == x)
  #   feat_glm_l = pblapply(ps_feats, function(y){
  #     d_ = d |> filter(var_name == y)
  #     d_$block_id = as.factor(d_$block_id)
  #     mod = betareg(logistic_AUC ~ months_since_harvest + (1|block_id)
  #               ,data = d_)
  #     return(mod)
  #   })
  #   names(feat_glm_l) = ps_feats
  #   return(feat_glm_l)
  # }
  # ,cl = 'future'
  # )
  # plan('sequential')
  # names(dataset_glm_l) = datasets
}

#----timeseries analysis using raw pixel values
{
  data_string = 'RawBasemapValues_withCHM_withXY_cleaned.parquet'
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
    
    plan('multisession', workers = 4)
    
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
      
      # # Assign 'id' to each pixel (efficient with data.table)
      # d[, id := .I]
      
      # Get x and y coordinates for each cell
      d[, X := xFromCell(r,1:ncell(r))]
      d[, Y := yFromCell(r,1:ncell(r))]
      
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
  
  stat_sig_file = paste0(bm_dir, '/StatDiff_ThinnedvsNotthinned.csv')
  
  if(!file.exists(stat_sig_file)){
    
    all_cols = names(results2_df)
    non_ps_cols = c('block_id', 'dataset', 'id', 'acquisition_date', 'julian_day', 'CHM_change')
    ps_feats = all_cols[!all_cols %in% non_ps_cols]
    ps_feats = ps_feats[ps_feats != 'VARIgreen'] #remove this index since it is undefined for many dates
    
    
    #make dataframe with every combination of block_id, dataset, acquisition_date, and var_name
    statistical_significance_summ = results2_df |>
      select(acquisition_date, dataset, block_id) |>
      distinct() |>
      collect() |>
      crossing(var_name = ps_feats) |>
      filter(var_name != 'VARIgreen')
    
    harvest_threshold = -5.2 #mean Otsu's threshold applied to the change in CHM layers
    
    # plan('multisession', workers = 10)
    
    stat_results = pblapply(1:nrow(statistical_significance_summ), function(i){
      # stat_results = pblapply(1:10, function(i){
      
      #initialize list where results will go
      results_list = list(NA, NA, NA, NA, NA)
      
      print(paste0('Processing ',i,'/',nrow(statistical_significance_summ)))
      
      sub_df = results2_df |>
        filter(dataset == statistical_significance_summ$dataset[i],
               acquisition_date == statistical_significance_summ$acquisition_date[i],
               block_id == statistical_significance_summ$block_id[i]) |>
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
        
        if(length(unique(thinned))==1){
          results_list[[1]] = 0
        } else{
          ad_t = ad.test(thinned)
          results_list[[1]] = ad_t[['p.value']]
        }
        if(length(unique(not_thinned))==1){
          results_list[[2]] = 0
        } else{
          ad_nt = ad.test(not_thinned)
          results_list[[2]] = ad_nt[['p.value']]
        }
        
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
    
    
    write_csv(statistical_significance_summ, stat_sig_file)
  }
  
  
  statistical_significance_summ = fread(stat_sig_file)
  
  #filter out NAs
  statistical_significance_summ = statistical_significance_summ[!is.na(difference_p)]
  
  #----join results with results_df (ie global stats for thinning vs not-thinning by block) ----
  
  alpha = 0.001
  
  
  results3_df <- inner_join( #join statistical significance results to global stats thinnedvsnotthinned
    results_df,
    statistical_significance_summ |> 
      mutate(acquisition_date2 = as.Date(acquisition_date)),
    by = c('dataset', 'acquisition_date2', 'block_id', 'var_name')
  ) |>
    mutate(significant_dif = ifelse(difference_p < alpha, 1, 0))
  
  #----plotting vegetation indices over time----
  {
    indices = c(
      "blue", "green", "red"
      ,
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
    
    
    pixel_stratum = c('Not_thinned', 'Thinned')
    
    subset_df = results3_df %>% 
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
      #add significant difference signs
      geom_point(data = subset_df %>% filter(significant_dif==1), aes(x = acquisition_date2, y = -2
                                                                      , shape = paste0('Difference (p<',alpha,')')))+
      guides(
        color = guide_legend(title = NULL),
        fill = guide_legend(title = NULL),
        linetype = guide_legend(title = NULL),
        shape = guide_legend(title = NULL)
      ) +
      ggtitle(unique(subset_df$dataset))
  }
  
  
  #----PCA----
  
  # datasets = results2_df |> 
  #   select(dataset) |>
  #   unique()|>
  #   collect()
  # datasets = datasets$dataset
  # 
  # ps_feats = c('blue', 'green', 'red', 'Hue', 'GCC', 'NDGR', 'BI', 'CI', 'CRI550', 'GLI', 'TVI')
  # 
  # pca_l = pblapply(datasets, function(x){
  #   pc = results2_df |>
  #     filter(dataset == x)|>
  #     select(all_of(ps_feats)) |>
  #     collect() |>
  #     drop_na()|>
  #     princomp()
  #   return(pc)
  # })
  # 
  # names(pca_l)= datasets
  # 
  # z = pca_l$Z
  
  
  #----fit mixed effects models to the timeseries for each feature and dataset----
  
  #create directory to store models
  lm_dir_string = 'linear_mixed_models_bydataset_byfeat_binaryharvest_binarythinning'
  lm_dir = paste0(bm_dir, '/',lm_dir_string)
  dir.check(lm_dir)
  
  ##using arrow dataset (i.e. library of parquet files, quite slow)
  {
  # #get vector of dataset names
  # datasets = results2_df |>
  #   dplyr::select(dataset) |>
  #   distinct()|>
  #   collect()
  # datasets = datasets$dataset
  # 
  # #vector of features to test
  # ps_feats = c('blue', 'green', 'red', 'Hue', 'GCC', 'NDGR', 'BI', 'CI', 'CRI550', 'GLI', 'TVI')
  # 
  # #reformat harvest dates
  # harvest_dates_df = harvest_dates_df |>
  #   mutate(harvest_date = as.Date(paste0(harvest_month,'-01')))
  # 
  # #create directory to store models
  # lm_dir_string = 'linear_mixed_models_bydataset_byfeat_binaryharvest_binarythinning'
  # lm_dir = paste0(bm_dir, '/',lm_dir_string)
  # dir.check(lm_dir)
  # 
  # #create and save models
  # 
  # 
  # # feat_lm_l = 
  # pblapply(datasets, function(x){
  #   
  #   d = results2_df |>
  #     filter(dataset == x)|>
  #     #get time since harvest
  #     left_join(harvest_dates_df, by = 'block_id') |>
  #     # collect()|>
  #     # mutate(months_since_harvest = as.numeric(acquisition_date - harvest_date, units='days')/30.44) |>
  #     mutate(harvested = ifelse(acquisition_date >= harvest_date, 1,0))|> #binary harvesting variable
  #     #classify CHM data using harvest threshold
  #     mutate(thinning_pixel = ifelse(is.na(CHM_change), NA,
  #                                    ifelse(CHM_change < harvest_threshold, 1, 0))) #binary thinning variable
  #   subdir = paste0(lm_dir,'/',x)
  #   dir.check(subdir)
  #   
  #   # f_l = 
  #   pblapply(ps_feats, function(y){
  #     
  #     print(paste0('Processing ',x,'-',y))
  #     filename = paste0(subdir,'/',y,'.rds')
  #     if(!file.exists(filename)){
  #       d_ = d |>
  #         filter(!is.na(!!sym(y))) |>
  #         dplyr::select(all_of(c(y, 'thinning_pixel', 'block_id', 'harvested'))) |>
  #         collect()
  #       
  #       # formula = as.formula(paste(y,'~ thinning_pixel * harvested * block_id + (1|block_id)'))
  #       # mod = lmer(formula, data = d)
  #       
  #       formula = as.formula(paste( y,'~ thinning_pixel* harvested + thinning_pixel + harvested + (1|block_id)'))
  #       mod = lmer(formula, data = d_)
  #       saveRDS(mod, file = filename)
  #     }
  #     
  #     # return(mod)
  #   }
  #   )
  #   # names(f_l) = ps_feats
  #   # return(f_l)
  # })
  }
  
  
  #using single parquet file
  results_pixels_df = setDT(read_parquet(data_filename)) #load and make data.table
  datasets = unique(results_pixels_df$dataset)
  ps_feats = c('blue', 'green', 'red', 'Hue', 'GCC', 'NDGR', 'BI', 'CI', 'CRI550', 'GLI', 'TVI')
  
  
  
  # cl = makeCluster(2)
  pblapply(datasets, function(x){
    
    d = results_pixels_df[dataset == x][
      setDT(harvest_dates_df), on = "block_id"
    ][
      , harvested := as.integer(acquisition_date >= harvest_date)
    ][
      , thinning_pixel := fifelse(
        is.na(CHM_change), 
        NA_integer_, 
        fifelse(CHM_change < harvest_threshold, 1L, 0L)
      )
    ]
    
    subdir = paste0(lm_dir,'/',x)
    dir.check(subdir)
    
    pblapply(ps_feats, function(y){
      
      print(paste0('Processing ',x,'-',y))
      filename = paste0(subdir,'/',y,'.rds')
      if(!file.exists(filename)){
        d_ = d |>
          filter(!is.na(!!sym(y))) |>
          dplyr::select(all_of(c(y, 'thinning_pixel', 'block_id', 'harvested')))
        # formula = as.formula(paste(y,'~ thinning_pixel * harvested * block_id + (1|block_id)'))
        # mod = lmer(formula, data = d)
        
        formula = as.formula(paste( y,'~ thinning_pixel* harvested + thinning_pixel + harvested + (1|block_id)'))
        mod = lmer(formula, data = d_)
        saveRDS(mod, file = filename)
      }
    }
    # ,cl = cl
    )
  })
  # stopCluster(cl)
  
  #load linear models for Z-score dataset
  lm_subdirs = list.dirs(lm_dir)
  lm_subdirs = lm_subdirs[lm_subdirs != lm_dir]
  datset_names = basename(lm_subdirs)
  
  z_lm_l = pblapply(lm_subdirs[basename(lm_subdirs) %in% c('Z', 'Non-normalized')], function(x){
    files_list = list.files(x, full.names = T, recursive = T)
    file_names = basename(file_path_sans_ext(files_list))
    
    m_l = pblapply(files_list, readRDS)
    names(m_l) = file_names
    return(m_l)
  })
  
  
}
