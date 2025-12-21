source('scene_global_stats_filtering_20250925.R')
library(ranger)

#set up directories
rf_dir = paste0(ps_dir,'/random_forest_models_nopredict_perscene'
                ,'_',harvest_threshold
)
dir.check(rf_dir)

#get full list of files to test
all_files_totest = unique(results_meta_df$file_path)
feats_to_use = unique(results_meta_df$var_name)

pblapply(all_files_totest, function(f){
  
  #get dataset from input filepath, create output directory if necessary
  data_dir = str_replace(dirname(dirname(f)),ps_dir,rf_dir)
  dir.check(data_dir)
  
  #get block id from input filepath, create output directory if necessary
  block_dir = str_replace(dirname(f), ps_dir,rf_dir)
  block_id = basename(block_dir)
  dir.check(block_dir)
  
  #get id of scene from filepath, make output filepath
  id = basename(file_path_sans_ext(f))
  rf_filename = paste0(block_dir,'/',id,'.rds')
  rf_failed_filename = paste0(block_dir,'/',id,'_NULL.rds') #make dummy filepath for if model fails
  
  #run random forest model if output file does not exist
  # cat('\rProcessing ',rf_filename)
  cat('\nProcessing',f)
  if(!file.exists(rf_filename) & !file.exists(rf_failed_filename)){
    
    #load and subset raster
    r = rast(f)
    r = r[[feats_to_use]]
    
    #load chm_change file
    chm = lidar_files[str_detect(lidar_files, block_id)] |> rast() #load height change raster
    ht = otsu_df$chm_change_otsu_threshold[otsu_df$block_id==block_id] #get height change threshold
    thinning_binary = ifel(chm <= ht,1,0) #binarize raster with threshold
    names(thinning_binary) = 'thinned' #name layer
    
    #prepare lidar data
    t_r = resample(thinning_binary,r)
    r_ = c(t_r,r)
    v = values(r_,dataframe=T) |> 
      na.omit() |> 
      mutate(thinned = ifelse(thinned == 1, 'Thinned','Not_thinned'))
    v$thinned = as.factor(v$thinned)
    
    #balance classes
    min_class_size <- v %>%
      count(thinned) %>%
      pull(n) %>%
      min()
    v_b <- v %>%
      group_by(thinned) %>%
      slice_sample(n = min_class_size) %>%
      ungroup()
    
    #fit random forest model, message if error thrown
    rf <- tryCatch({
      ranger(thinned ~ ., data = v_b, importance = 'permutation', write.forest = F)
    }, error = function(e) {
      message('\nError fitting RF for file: ', f)
      message('Error: ', e$message)
      return(NULL)
    })
    
    # save model if successful
    if (!is.null(rf)){ 
      
      saveRDS(rf, rf_filename)
      
    } else { #save dummy file if not successful (so you don't have to try and redo all the processing again)
      dummy_rf <- list(
        model = NULL,
        error = TRUE,
        file_path = f,
        timestamp = Sys.time(),
        class = "failed_ranger")
      
      saveRDS(dummy_rf, rf_failed_filename)
      }
    
    
  }
  cat('\nDone ',match(f,all_files_totest),'/',length(all_files_totest))
})
