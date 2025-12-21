source('scene_global_stats_filtering_20250925.R')
# source('scene_global_stats_filtering_20250925.R')

lid_metrics_dir = paste0('data/lidar_metrics_Res=',target_res)
dir.check(lid_metrics_dir)

lid_metrics_subdirs = list.dirs(lid_metrics_dir)
lid_metrics_subdirs = lid_metrics_subdirs[!str_detect(lid_metrics_subdirs, lid_metrics_dir)]
  
  #----load packages, define directories----

#don't load these packages unless you need to (since they mess up terra plotting functions)  
if(length(list.files(lid_metrics_dir, recursive = T)) %% length(block_ids) != 0){
  library(lidR)
  library(lasR)
}
library(tictoc)
  
raw_las_parent_dir = "G:/quesnel_thinning_las"

pre_dir = paste0(raw_las_parent_dir,'/pre')
pre_blockdirs = list.dirs(pre_dir, recursive=F)

post_dir = paste0(raw_las_parent_dir,'/post')
post_blockdirs = list.dirs(post_dir, recursive=F)
  
  #----inspect point clouds----
  
  # pre_harvest_las = pblapply(pre_blockdirs, function(x)readLAScatalog(file.path(x, "input/las/norm")))
  # names(pre_harvest_las) = basename(pre_blockdirs)
  # 
  # post_harvest_las = pblapply(post_blockdirs, function(x)readLAScatalog(file.path(x, "input/las/norm")))
  # names(post_harvest_las) = basename(post_blockdirs)
  
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
    lid_mets_change_dir = paste0(lid_metrics_dir,'/change')
    dir.check(lid_mets_change_dir)
    
    lapply(block_ids, function(b){
      #check if file exists
      f = paste0(lid_mets_change_dir,'/',b,'.tif')
      if(!file.exists(f)){
        
        cat(b)
        #load pre and post-thinning metrics
        pre_f = list.files(lid_mets_pre_dir, full.names = T) |> keep(~ str_detect(.x,b))
        pre_r = rast(pre_f)
        
        post_f = list.files(lid_mets_post_dir, full.names = T) |> keep(~ str_detect(.x,b))
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
  lid_mets_change_cropped_dir = paste0(lid_metrics_dir,'/change_cropped')
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
  

