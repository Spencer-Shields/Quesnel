#----packages
{
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
  library(patchwork)
  library(caret)
  #extra functions
  source('helper_functions.R')
  #functions for checking or visualizing radiometric consistency
  source('https://raw.githubusercontent.com/Spencer-Shields/planetscope_radiometric_correction/refs/heads/main/check_radiometric_consistency.R')
}

#----data preprocessing, directory setup
{
  #----set global options so that terra saves rasters as FLT8S----
  terraOptions(datatype = 'FLT8S')
  
  #----harvest dates----
  # harvest_dates_df = tribble(
  #   ~block_id, ~harvest_month, ~note,
  #   '12L_C5', '2022-08', '-',
  #   '12L_D345','2022-07', '-',
  #   '12L_C4', '2022-11', '-',
  #   '12L_C7', '2022-09', '-',
  #   '12L_B8C3', '2022-09', '-',
  #   '12N_T3', '2023-02', 'partial, complete 2023-03',
  #   '12N_1X1W', '2024-03', '-',
  #   '12N_V1', '2024-03', 'partial, complete 2024-04'
  # )
  
  harvest_dates_file = 'data/Quesnel_thinning/harvest_dates.csv'
  harvest_dates_df = read_csv(harvest_dates_file)
  
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
  blocks = st_read('data/Quesnel_thinning/12l_12n_bdy.gpkg')
  no_change = st_read('Arc project/NoChangeStands_conifer.shp')
  no_change = no_change %>% rename(geom = geometry)
  blocks = rbind(blocks, no_change)
  blocks_p = st_transform(blocks, crs = crs(r_l[[1]]))
  blocks_v = vect(blocks_p)
  blocks_wr = wrap(blocks_v)
  thinning_block_ids = blocks$BLOCKNUM[!str_detect(blocks$BLOCKNUM,'NoChange')]
  
  aoi = st_read('data/AOI_fullsite_wgs84.kml') %>%
    st_transform(crs = crs(r_l[[1]]))
  
  # ggRGB(r_l[[1]])+
  #   geom_sf(data = blocks_p, aes(color = BLOCKNUM))
  # 
  # ggRGB(r_l[[2]]) +
  #   geom_sf(data = blocks_p, aes(color = BLOCKNUM))
  
  #load non-vegetation mask data
  nonveg_mask_dir = 'data/Quesnel_thinning/nonveg_mask_pre=0.5_post=0.3'
  nonveg_mask_files = list.files(nonveg_mask_dir, full.names = T)
  
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
  # plan('multisession',workers=6)
  
  pblapply(desired_month_dirs, function(d){
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
  # plan('multisession', workers = 4)
  
  pblapply(files_list, function(x){
    output_file = paste0(outdir, '/', basename(file_path_sans_ext(x)),'.tif')
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
      writeRaster(vi, output_file, datatype='FLT8S')
      gc()
    }
  })
  
  
  #----crop data to blocks, apply non-vegetation mask----
  
  indir = outdir
  rasters = list.files(indir, pattern = '\\.tif$', full.names = T)
  
  clipped_string = 'BlockClipped'
  outdir = paste0(indir, '_', clipped_string)
  nonorm_dir = outdir #save this directory name for later
  dir.check(outdir)
  
  lapply(1:length(thinning_block_ids), function(i){ #create subdirectories for every block
    d = paste0(outdir,'/',thinning_block_ids[i])
    dir.check(d)
  })
  
  pblapply(1:length(rasters), function(i){
    r_id = basename(file_path_sans_ext(rasters[i]))
    r = rast(rasters[i])
    print(i)
    pblapply(1:length(thinning_block_ids), function(j){
      print(j)
      block_id = thinning_block_ids[j]
      filename = paste0(nonorm_dir,'/',block_id,'/',r_id,'.tif')
      print(paste0('Processing ',filename))
      if(!file.exists(filename)){
        #crop basemap
        bu = unwrap(blocks_wr)
        b = bu[j,]
        r_c = crop(r,b,mask=T)
        #non-vegetation mask
        mf = nonveg_mask_files[which(str_detect(nonveg_mask_files,block_id))]
        m_ = rast(mf)
        m_ = project(m_, r_c)
        r_m = mask(r_c,m_)
        
        #save raster
        terra::writeRaster(r_m, filename)
      }
    })
  })
  
  
  #----Calculate z-scores, softmax, and robust z-scores using cropped data----
  
  indir = outdir
  rasters = list.files(indir, pattern = '\\.tif$', full.names = T, recursive = T)
  
  z_dir_string = 'Z'
  z_dir = paste0(indir, '_', z_dir_string)
  dir.check(z_dir)
  
  lapply(1:length(thinning_block_ids), function(i){ #create subdirectories for every block
    d = paste0(z_dir,'/',blocks_p$BLOCKNUM[i])
    dir.check(d)
  })
  
  sm_dir_string = 'SM'
  sm_dir = paste0(indir, '_', sm_dir_string)
  dir.check(sm_dir)
  
  lapply(1:length(thinning_block_ids), function(i){ #create subdirectories for every block
    d = paste0(sm_dir,'/',blocks_p$BLOCKNUM[i])
    dir.check(d)
  })
  
  zr_dir_string = 'Zrobust'
  zr_dir = paste0(indir, '_', zr_dir_string)
  dir.check(zr_dir)
  
  lapply(1:length(thinning_block_ids), function(i){ #create subdirectories for every block
    d = paste0(zr_dir,'/',blocks_p$BLOCKNUM[i])
    dir.check(d)
  })
  
  scale_factor = 255 #scale factor to use for calculating softmax (i.e. maximum possible band pixel value)
  
  bands = c('blue', 'green', 'red')
  
  # cl = makeCluster(ceiling(cores/4))
  # plan('multisession', workers = 12)
  
  #calculate z-scores on cropped data
  pblapply(1:length(rasters), function(i){
    filename = str_replace(rasters[i], #get filename by modifying the name of the parent directory
                           basename(indir), 
                           paste0(basename(indir),'_',z_dir_string))
    print(paste('Processing',filename))
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
  pblapply(1:length(rasters), function(i){
    filename = str_replace(rasters[i], #get filename by modifying the name of the parent directory
                           basename(indir), 
                           paste0(basename(indir),'_',zr_dir_string))
    print(paste('Processing',filename))
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
  pblapply(1:length(rasters), function(i){
    filename = str_replace(rasters[i], #get filename by modifying the name of the parent directory
                           basename(indir), 
                           paste0(basename(indir),'_',sm_dir_string))
    print(paste('Processing',filename))
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
      scene_sm = scene_sm * 10000 #multiply to avoid errors from having very small floating point numbers
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
  
  lapply(1:length(thinning_block_ids), function(i){ #create subdirectories for every block
    d = paste0(dz_dir,'/',thinning_block_ids[i])
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
      print(paste('Processing',filename))
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
  
  lapply(1:length(thinning_block_ids), function(i){ #create subdirectories for every block
    d = paste0(dn_dir,'/',blocks_p$BLOCKNUM[i])
    dir.check(d)
  })
  
  indir = nonorm_dir #input files
  
  in_subdirs = list.dirs(indir)
  in_subdirs = in_subdirs[in_subdirs != indir]
  
  #calculate difference rasters
  pblapply(in_subdirs, function(d){
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
      print(paste('Processing',filename))
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
  
  lapply(1:length(thinning_block_ids), function(i){ #create subdirectories for every block
    d = paste0(dsm_dir,'/',blocks_p$BLOCKNUM[i])
    dir.check(d)
  })
  
  indir = sm_dir #input files
  
  in_subdirs = list.dirs(indir)
  in_subdirs = in_subdirs[in_subdirs != indir]
  
  #calculate difference rasters
  pblapply(in_subdirs, function(d){
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
      print(paste('Processing',filename))
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
  
  lapply(1:length(thinning_block_ids), function(i){ #create subdirectories for every block
    d = paste0(dzr_dir,'/',blocks_p$BLOCKNUM[i])
    dir.check(d)
  })
  
  indir = zr_dir #input files
  
  in_subdirs = list.dirs(indir)
  in_subdirs = in_subdirs[in_subdirs != indir]
  
  #calculate difference rasters
  pblapply(in_subdirs, function(d){
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
      print(paste('Processing',filename))
      if(!file.exists(filename)){
        rd = rast(files_sorted[i]) - rast(files_sorted[i-1])
        terra::writeRaster(rd, filename = filename)
      }
    })
  })
  
  
  
}