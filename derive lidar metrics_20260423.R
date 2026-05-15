
source('scene_setup_preprocessing_20250909.R')


lid_metrics_dir = paste0('data/lidar_metrics_Res=',target_res)
dir.check(lid_metrics_dir)

lid_metrics_subdirs = list.dirs(lid_metrics_dir)
lid_metrics_subdirs = lid_metrics_subdirs[!str_detect(lid_metrics_subdirs, lid_metrics_dir)]
  
  #----load packages, define directories----


library(tictoc)
  
raw_las_parent_dir = "G:/quesnel_thinning_las"

pre_dir = paste0(raw_las_parent_dir,'/pre')
pre_blockdirs = list.dirs(pre_dir, recursive=F)

post_dir = paste0(raw_las_parent_dir,'/post')
post_blockdirs = list.dirs(post_dir, recursive=F)
  
  #----inspect point clouds----
  # 
  # pre_harvest_las = pblapply(pre_blockdirs, function(x)readLAScatalog(file.path(x, "input/las/norm"), filter='-keep_first'))
  # names(pre_harvest_las) = basename(pre_blockdirs)
  # 
  # post_harvest_las = pblapply(post_blockdirs, function(x)readLAScatalog(file.path(x, "input/las/norm"), filter = '-keep_first'))
  # names(post_harvest_las) = basename(post_blockdirs)
  # 
  #----process point clouds, generate lidar metrics----
  
  #set parameters for building pipelines
  
  metrics = c(
    #height
    'z_max', 'z_min', paste0('z_p',seq(5,95,5)), 'z_mean'
    #canopy cover
    , paste0('z_above',seq(1,5,1))
    #height variability
    ,'z_sd', 'z_cv'
    #number of points in each pixel
    ,'count'
  )
  
  pipe = lasR::reader() + lasR::rasterize(target_res, operators = metrics, filter=lasR::keep_first())
  
  #generate and save pre-harvest metrics
  {
    lid_mets_pre_dir = paste0(lid_metrics_dir,'/pre')
    dir.check(lid_mets_pre_dir)
    
    tic()
    pblapply(pre_blockdirs, function(lc){
      
      #get block from lascatalogue directory
      b = basename(lc)
      cat(b)
      
      #check if file exists
      f = paste0(lid_mets_pre_dir,'/',b,'.tif')
      if(!file.exists(f)){
        
        library(lasR)
        
        lc_ = paste0(lc,'/input/las/norm')
        # set_parallel_strategy(concurrent_files(12))
        r = lasR::exec(pipe
                       , on = lc_, with = list(progress = T))
        if(terra::crs(r)!=terra::crs(blocks_p)){
          r = terra::project(r, terra::crs(blocks_p))
        }
        terra::writeRaster(r, f)
      }
    })
    toc()
  }
  
  #generate and save post-harvest metrics
  {
    lid_mets_post_dir = paste0(lid_metrics_dir,'/post')
    dir.check(lid_mets_post_dir)
    
    tic()
    pblapply(post_blockdirs, function(lc){
      
      #get block from lascatalogue directory
      b = basename(lc)
      cat(b)
      
      #check if file exists
      f = paste0(lid_mets_post_dir,'/',b,'.tif')
      if(!file.exists(f)){
        
        library(lasR)
        
        lc_ = paste0(lc,'/input/las/norm')
        set_parallel_strategy(concurrent_files(12))
        r = lasR::exec(pipe
                       , on = lc_, with = list(progress = T))
        if(terra::crs(r)!=terra::crs(blocks_p)){
          r = terra::project(r, terra::crs(blocks_p))
        }
        terra::writeRaster(r, f)
      }
    })
    toc()
  }
  
  #generate and save post-harvest metrics with decimation
  {
    
    #get pulse densities from pre_harvest and post-harvest data
    pre_harv_pulse_densities = sapply(list.files(lid_mets_pre_dir, full.names = T), function(x){
      r = rast(x, lyrs='count')
      mean_pulse_count = terra::global(r,'mean',na.rm=T)[1,1]
      pixel_area = prod(res(r))
      mean_pulse_density = mean_pulse_count/pixel_area
      return(mean_pulse_density)
    })
    names(pre_harv_pulse_densities) = basename(file_path_sans_ext(list.files(lid_mets_pre_dir, full.names = T)))

    post_harv_pulse_densities = sapply(list.files(lid_mets_post_dir, full.names = T), function(x){
      r = rast(x, lyrs='count')
      mean_pulse_count = terra::global(r,'mean',na.rm=T)[1,1]
      pixel_area = prod(res(r))
      mean_pulse_density = mean_pulse_count/pixel_area
      return(mean_pulse_density)
    })
    names(post_harv_pulse_densities) = basename(file_path_sans_ext(list.files(lid_mets_post_dir, full.names = T)))
    
    pulse_densities_tbl = tibble(
      pre_harv = pre_harv_pulse_densities,
      post_harv = post_harv_pulse_densities,
      block_id = basename(file_path_sans_ext(list.files(lid_mets_pre_dir, full.names = T)))
    ) %>%
      mutate(
        # Target density is pre_harv (points per m²)
        # sample_pixel keeps ~1 point per pixel
        # Density = 1 / (res²)
        # Therefore: res = sqrt(1 / target_density)
        sample_pixel_res = sqrt(1 / pre_harv)
      )
    
    #unload terra to avoid conflicts
    detach("package:terra", unload = TRUE)
    
    
    # target_pulse_density = 25
    
    #set directory
    lid_mets_post_dec_dir = paste0(lid_metrics_dir,'/post_decimated')
    dir.check(lid_mets_post_dec_dir)
    
    tic()
    pblapply(post_blockdirs, function(lc){
      
      #get block from lascatalogue directory
      b = basename(lc)
      cat(b)
      
      #get spatial resolution to match target density
      spat_res = pulse_densities_tbl$sample_pixel_res[pulse_densities_tbl$block_id==b]
      names(spat_res) = NULL
      
      # #set pipe
      pipe2 = lasR::reader(filter = lasR::keep_first()) +  # Apply filter at reader stage
        lasR::sampling_pixel(spat_res, method = 'random') +  # Decimate point cloud (already filtered to first returns)
        lasR::rasterize(target_res, operators = metrics)  # Rasterize (no need to filter again)

      
      
      #check if file exists
      f = paste0(lid_mets_post_dec_dir,'/',b,'.tif')
      if(!file.exists(f)){
        
        #unload terra to avoid conflicts
        detach("package:terra", unload = TRUE)
        library(lasR)
        
        lc_ = paste0(lc,'/input/las/norm')
        
        
        # ctg = readLAScatalog(lc_, filter = "-keep_first")
        # 
        # # Set catalog options
        # opt_chunk_buffer(ctg) = 10  # Buffer to avoid edge effects
        # opt_chunk_size(ctg) = 0     # Process by file
        # opt_progress(ctg) = TRUE
        # opt_output_files(ctg) = paste0(tempdir(), "/{ORIGINALFILENAME}_decimated")
        # 
        # # Decimate using homogenize - keeps only first returns and thins to target density
        # ctg_decimated = lidR::decimate_points(las = ctg, algorithm = lidR::homogenize(density=target_pulse_density, res=10))
        
        
        set_parallel_strategy(concurrent_files(12))
        # r = lasR::exec(pipe, on = ctg_decimated, with = list(progress = T))
        
        r = lasR::exec(pipe2, on = lc_, with = list(progress = T))
        
        if(terra::crs(r)!=terra::crs(blocks_p)){
          r = terra::project(r, terra::crs(blocks_p))
        }
        terra::writeRaster(r, f)
      }
    })
    toc()
  }
  
  #detach lidR and lasR to avoid function namespace conflicts 
  # if('lidR' %in% loadedNamespaces()){detach('package:lidR', unload=T)} #check if package is loaded, remove if it is
  # if('lasR' %in% loadedNamespaces()){detach('package:lasR', unload=T)}
  for (pkg in c("lidR", "lasR")) {
    if (paste0("package:", pkg) %in% search()) {
      detach(paste0("package:", pkg), unload = TRUE, character.only = TRUE)
    } else if (pkg %in% loadedNamespaces()) {
      unloadNamespace(pkg)
    }
  }
  
  #generate and save change metrics
  {
    lid_mets_change_dir = paste0(lid_metrics_dir,'/change_decimated')
    dir.check(lid_mets_change_dir)
    
    lapply(block_ids, function(b){
      #check if file exists
      f = paste0(lid_mets_change_dir,'/',b,'.tif')
      if(!file.exists(f)){
        
        cat(b)
        #load pre and post-thinning metrics
        pre_f = list.files(lid_mets_pre_dir, full.names = T) |> keep(~ str_detect(.x,b))
        pre_r = rast(pre_f)
        
        post_f = list.files(lid_mets_post_dec_dir, full.names = T) |> keep(~ str_detect(.x,b))
        post_r = rast(post_f)
        
        #resample post-harvest data to match extent of pre-harvest, calculate change as post-pre
        post_r = resample(post_r, pre_r)
        
        change = post_r - pre_r
        
        #save
        terra::writeRaster(change,f)
      }
    })
  }
  
  #save change metrics which are cropped to the thinning blocks
  {
  lid_mets_change_cropped_dir = paste0(lid_metrics_dir,'/change_decimated_cropped')
  dir.check(lid_mets_change_cropped_dir)
  
  lapply(block_ids, function(b){
    #check if file exists
    f = paste0(lid_mets_change_cropped_dir,'/',b,'.tif')
    if(!file.exists(f)){
      
      cat(b)
      #load change metrics
      change = list.files(lid_mets_change_dir, full.names = T) |> 
        keep(~ str_detect(.x,b)) |>
        rast()
      block = blocks_p[blocks_p$BLOCKNUM==b,]
      
      change_c = terra::crop(change, block, mask=T)
      
      #save
      terra::writeRaster(change_c,f)
    }
  })
  }
  
  #save pre and post-harvest metrics which are cropped to the thinning blocks
  {
    lid_mets_pre_cropped_dir = paste0(lid_metrics_dir,'/pre_cropped')
    dir.check(lid_mets_pre_cropped_dir)
    
    lapply(block_ids, function(b){
      #check if file exists
      f = paste0(lid_mets_pre_cropped_dir,'/',b,'.tif')
      if(!file.exists(f)){
        
        cat(b)
        #load uncropped metrics
        change = list.files(lid_mets_pre_dir, full.names = T) |> 
          keep(~ str_detect(.x,b)) |>
          rast()
        block = blocks_p[blocks_p$BLOCKNUM==b,]
        
        change_c = terra::crop(change, block, mask=T)
        
        #save
        terra::writeRaster(change_c,f)
      }
    })
    
    lid_mets_post_dec_cropped_dir = paste0(lid_metrics_dir,'/post_decimated_cropped')
    dir.check(lid_mets_post_dec_cropped_dir)
    
    lapply(block_ids, function(b){
      #check if file exists
      f = paste0(lid_mets_post_dec_cropped_dir,'/',b,'.tif')
      if(!file.exists(f)){
        
        cat(b)
        #load uncropped metrics
        change = list.files(lid_mets_post_dec_dir, full.names = T) |> 
          keep(~ str_detect(.x,b)) |>
          rast()
        block = blocks_p[blocks_p$BLOCKNUM==b,]
        
        change_c = terra::crop(change, block, mask=T)
        
        #save
        terra::writeRaster(change_c,f)
      }
    })
  }
  
  library(tidyverse)
  library(terra)
  
  #plotting to check data
  {
  #   fl = list.files(lid_metrics_post, full.names = T)
  #   harvest_metrics = lapply(fl, rast)
  #   names(harvest_metrics) = basename(file_path_sans_ext(fl))
  # 
  #   b = block_ids[i]
  # 
  #   plot(harvest_metrics[[b]][['z_above2']])
  #   plot(blocks_p[blocks_p$BLOCKNUM==b,], add=T, border='red')
  }
  
  # ggplot()+geom_spatraster(data=change_c, aes(fill=z_max))+geom_spatvector(data=blocks_p[blocks_p$BLOCKNUM==b,], color='red', fill=NA)
  

