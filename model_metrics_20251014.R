source('derive lidar metrics_202501020.R')
library(ranger)
# source('scene_SeparabilityAnalysis_20250925.R')

#----choose pre and post harvest dates----

#get ids which are common to each thinning block
ds_ids = sapply(unique(results_meta_df$dataset), function(d)results_meta_df$id[results_meta_df$dataset==d])
common_ids = Reduce(intersect, ds_ids)

common_ids <- results_meta_df %>%
  group_by(block_id, dataset) %>%
  summarise(ids = list(unique(id)), .groups = "drop") %>%
  summarise(common_ids = list(Reduce(intersect, ids)), .by = block_id)
ids_in_all_blocks <- Reduce(intersect, common_ids$common_ids)


#generate random list of pre-harvest and post-harvest growing season images that are common to each block/dataset
gs_months = c(5,6,7,8,9,10) #define growing season months

n_combos = 19

set.seed(123)
random_pre_harv_ids = sample(size = n_combos, x = unique(results_meta_df$id[results_meta_df$id %in% ids_in_all_blocks &
                                                                              month(as.Date(results_meta_df$acquired)) %in% gs_months &
                                                                              as.Date(results_meta_df$acquired)<as.Date(min(harvest_dates_df$harvest_start_date))]))
random_pre_harv_ids = c("20210713_191755_66_2414", random_pre_harv_ids) #manually add scene which is close to pre-harv id

set.seed(123)
random_post_harv_ids = sample(size = n_combos, x = unique(results_meta_df$id[results_meta_df$id %in% ids_in_all_blocks &
                                                                               month(as.Date(results_meta_df$acquired)) %in% gs_months &
                                                                               as.Date(results_meta_df$acquired)>as.Date(max(harvest_dates_df$harvest_finish_date))]))
random_post_harv_ids = c("20241011_192907_10_251a", random_post_harv_ids) #manually add scene which is close to post-harv id

#load lidar change metrics
print('Load lidar change metrics')
change_lids = lapply(block_ids, function(b){
  fl = list.files(lid_mets_change_cropped_dir, full.names = T)
  f = fl[str_detect(fl,b)]
  r = rast(f)
  return(r)
})
names(change_lids) = block_ids

#generate a lidar change mosaic
lid_change_mosaic = sprc(change_lids) |> merge()

#iterate over date pairs
pbmapply(function(prehid, posthid){
  
  gc()
  
  pre_harv_id = prehid
  post_harv_id = posthid
  
  cat('Processing',pre_harv_id,'---',post_harv_id)
  
  #----load data, make mosaics----
  
  # cat('\nLoad Planetscope data, make mosaics')
  
  #load pre-harvest scenes
  # print('Load pre-harvest scenes')
  pre_harvest_scenes = lapply(unique(results_meta_df$dataset), function(d){
    bl = lapply(block_ids, function(b){
      f = results_meta_df$file_path[results_meta_df$id == pre_harv_id &
                                      results_meta_df$block_id==b &
                                      results_meta_df$dataset==d][1]
      r = rast(f)
      
      r = r[[!names(r) %in% 'WDVI']] #remove WDVI since it's the same as DVI
      
      return(r)
    })
    names(bl) = block_ids
    return(bl)
  })
  names(pre_harvest_scenes) = unique(results_meta_df$dataset)
  
  #generate pre-harvest mosaics (i.e. combine thinning blocks into single raster stack)
  # print('Generate pre-harvest mosaics')
  pre_harv_mosaics = lapply(pre_harvest_scenes, function(x){
    sc = sprc(x)
    r = merge(sc)
  })
  
  #load post-harvest scenes
  # print('Load post harvest scenes')
  post_harvest_scenes = lapply(unique(results_meta_df$dataset), function(d){
    bl = lapply(block_ids, function(b){
      f = results_meta_df$file_path[results_meta_df$id == post_harv_id &
                                      results_meta_df$block_id==b &
                                      results_meta_df$dataset==d][1]
      r = rast(f)
      r = r[[!names(r) %in% 'WDVI']]
      return(r)
    })
    names(bl) = block_ids
    return(bl)
  })
  names(post_harvest_scenes) = unique(results_meta_df$dataset)
  
  # print('Generate mosaics')
  post_harv_mosaics = lapply(post_harvest_scenes, function(x){
    sc = sprc(x)
    r = merge(sc)
  })
  
  #calculate change mosaics
  # print('Generate change mosaics')
  change_mosaics = mapply(function(x, y){
    x-y
  }, pre_harv_mosaics, post_harv_mosaics, SIMPLIFY = F)
  
  
  #----select lidar metrics and planetscope features to use for modelling----
  
  # print('Set PlanetScope features to use')
  #lidar
  mets_to_test = c( 
    "z_above2","z_p95","z_sd"
    # ,
    # "z_max",    
    #"z_min"    "z_p5"     "z_p10"    "z_p15"    "z_p20"    "z_p25"    "z_p30"    "z_p35"    "z_p40"   
    # "z_p45"    "z_p50"    "z_p55"    "z_p60"    "z_p65"    "z_p70"    "z_p75"    "z_p80"    "z_p85"    
    # "z_p90",   
    #    #"z_mean",   
    # "z_above1", "z_above3", "z_above4", "z_above5"
    # , "z_cv" #,"count"
  )
  
  #ps spectral
  spec.feats4 = c("coastal_blue", "blue"        , "green_i"     , "green"      ,  "yellow"       ,"red"         , "rededge"   ,   "nir"         
                  , "CTVI"       ,  "DVI"        ,  "EVI"        ,  "GEMI"      ,   "GNDVI"       , "KNDVI"      ,  "MCARI"    ,    "MSAVI"       
                  , "MSAVI2"      , "NDVI"        , "NRVI"        , "SAVI"       ,  "SR"           ,"TVI"         , "Hue"       ,   "GCC"         
                  , "NDGR"         ,"BI"           ,"CRI550"       ,"GLI"         , "VARIgreen"    ,"EVI2"         ,"YNDVI" )
  
  # #ps texture
  # glcm_mets = c(
  #   "glcm_contrast",
  #   # "glcm_dissimilarity",
  #   #  "glcm_homogeneity", "glcm_ASM",
  #   "glcm_entropy", 
  #   "glcm_mean"
  #   # , "glcm_variance", "glcm_correlation"
  # )
  # glcm_quant_lvls = 32
  # glcm_quant_method = 'range'
  # glcm_windows = c(3,5,7,9,11
  #                  # ,15,21
  #                  )
  
  
  feats_to_test = spec.feats4
  
  #remove unnecessary spectral features from mosaics
  # print('Remove unnecessary features from PlanetScope mosiacs')
  # pre_harv_mosaics = lapply(pre_harv_mosaics, function(x)x[[feats_to_test]])
  # post_harv_mosaics =lapply(post_harv_mosaics, function(x)x[[feats_to_test]])
  change_mosaics = lapply(change_mosaics, function(x)x[[feats_to_test]])
  
  #----generate/load glcm layers----
  # require(GLCMTextures)
  # 
  # glcm_mosaics_dir = paste0(ps_dir,'/GLCM_mosaics')
  # dir.check(glcm_mosaics_dir)
  # 
  # cat('\nGenerating pre-harvest GLCM layers')
  # glcm_pre = pblapply(names(pre_harv_mosaics), function(nom){ #iterate over datasets
  #   
  #   r = pre_harv_mosaics[[nom]][[feats_to_test]]
  #   
  #   wins = pblapply(glcm_windows, function(win){ #iterate over window sizes
  #     feats = pblapply(names(r), function(feat){
  #       
  #       #make directory for textures
  #       g_dir = paste0(glcm_mosaics_dir,'/','feat=',feat,'_w=',win,'_mets=',str_c(str_replace_all(glcm_mets,'glcm_',''),collapse=','),'_quant=',glcm_quant_lvls,'_quantmethod=',glcm_quant_method)
  #       dir.check(g_dir)
  #       
  #       #define file for dataset, check if it exists, make if it doesn't load if it does
  #       ds_dir = paste0(g_dir,'/',nom)
  #       dir.check(ds_dir)
  #       
  #       f = paste0(ds_dir,'/',pre_harv_id,'.tif')
  #       
  #       if(!file.exists(f)){
  #         g = GLCMTextures::glcm_textures(r[[feat]], w=win, n_levels = glcm_quant_lvls, quant_method = glcm_quant_method, metrics = glcm_mets)
  #         names(g) = paste0(names(g),'_',feat,'_',win)
  #         writeRaster(g, f)
  #       } else {
  #         g = rast(f)
  #       }
  #       return(g)
  #     }) |> rast() #make single spatraster from all layers for a single window size
  #     return(feats)
  #   }) |> rast() #make single spatraster for all layers for a single dataset
  #   return(wins)
  # })
  # 
  # cat('\nGenerating post-harvest GLCM layers')  
  # glcm_post = pblapply(names(post_harv_mosaics), function(nom){ #iterate over datasets
  #   
  #   r = post_harv_mosaics[[nom]][[feats_to_test]]
  #   
  #   wins = pblapply(glcm_windows, function(win){ #iterate over window sizes
  #     feats = pblapply(names(r), function(feat){
  #       
  #       #make directory for textures
  #       g_dir = paste0(glcm_mosaics_dir,'/','feat=',feat,'_w=',win,'_mets=',str_c(str_replace_all(glcm_mets,'glcm_',''),collapse=','),'_quant=',glcm_quant_lvls,'_quantmethod=',glcm_quant_method)
  #       dir.check(g_dir)
  #       
  #       #define file for dataset, check if it exists, make if it doesn't load if it does
  #       ds_dir = paste0(g_dir,'/',nom)
  #       dir.check(ds_dir)
  #       
  #       f = paste0(ds_dir,'/',post_harv_id,'.tif')
  #       
  #       if(!file.exists(f)){
  #         g = GLCMTextures::glcm_textures(r[[feat]], w=win, n_levels = glcm_quant_lvls, quant_method = glcm_quant_method, metrics = glcm_mets)
  #         names(g) = paste0(names(g),'_',feat,'_',win)
  #         writeRaster(g, f)
  #       } else {
  #         g = rast(f)
  #       }
  #       return(g)
  #     }) |> rast() #make single spatraster from all layers for a single window size
  #     return(feats)
  #   }) |> rast() #make single spatraster for all layers for a single dataset
  #   return(wins)
  # })
  # 
  # print('Generating GLCM change mosaics')
  # glcm_change = mapply(function(x,y){x-y}, glcm_pre, glcm_post)
  # names(glcm_change) = names(pre_harv_mosaics)
  
  #----write function for modelling lidar for a pre-harvest and post-harvest scene using variable selection and GLCM textures----
  
  
  metric_predictor_fn4 = function(met, ps_scene,
                                  lid=lid_change_mosaic#,
                                  # glcm_mets = c(
                                  #   "glcm_contrast",
                                  #   # "glcm_dissimilarity",
                                  #   #  "glcm_homogeneity", "glcm_ASM",
                                  #   "glcm_entropy", 
                                  #   "glcm_mean"
                                  #   # , "glcm_variance", "glcm_correlation"
                                  # ),
                                  # glcm_rast = NULL, #spatraster with all possible glcm layers, leave as NULL to make this on-the-fly using ps_scene
                                  # glcm_quant_lvls = glcm_quant_lvls,
                                  # glcm_quant_method = glcm_quant_method
                                  # ,
                                  # glcm_window = 3
                                  ){
    
    
    #get lidar change metric
    lid_met = lid[[met]]
    
    #make glcm_textures for the selected variables
    # if(is.null(glcm_rast)){
    #   glcm_r = pblapply(names(ps_scene), function(x){
    #     r = ps_scene[[x]]
    #     g = glcm_textures(r, w = c(glcm_window, glcm_window), n_levels = glcm_quant_lvls, quant_method = glcm_quant_method, metrics = glcm_mets)
    #     names(g) = paste0(x,'_',names(g))
    #     return(g)
    #   }) |> rast()
    # } else {
    #   glcm_lyrs_wanted = as.vector(outer(names(ps_scene), glcm_mets, FUN = function(x, y) paste0(x, "_", y))) #subset the premade raster to only have the desired layers
    #   glcm_r = glcm_rast[[glcm_lyrs_wanted]]
    # }
    # 
    # glcm_lyrs_wanted = as.vector(outer(names(ps_scene), glcm_mets, FUN = function(x, y) paste0(x, "_", y))) #subset the premade raster to only have the desired layers
    # glcm_r = glcm_rast[[glcm_lyrs_wanted]]
    
    # glcm_r = glcm_rast
    
    #stack lidar metric, planetscope scene, and planetscope glcm layers
    if(compareGeom(lid_met, ps_scene, stopOnError=F)==F){ps_scene = resample(ps_scene_scene, lid_met)}
    # if(compareGeom(lid_met, glcm_r, stopOnError=F)==F){glcm_r = resample(glcm_r, lid_met)}
    combined_scene = c(lid_met, ps_scene#, glcm_r
                       )
    
    #get random sample of pixels based on size of smallest block 
    blocks_rast = rasterize(blocks_p, combined_scene, 'BLOCKNUM')
    blocks_ncells_summ = values(blocks_rast, dataframe=T) |>
      group_by(BLOCKNUM) |>
      summarise(n = n())
    
    set.seed(123)
    combined_values_sample = lapply(block_ids, function(b){
      block = blocks_p[blocks_p$BLOCKNUM==b,]
      r_b = crop(combined_scene, block, mask=T)
      s = spatSample(x=r_b, size=min(blocks_ncells_summ$n), method='random', na.rm=T, exhaustive=T)
      return(s)
    }) |> bind_rows()
    
    # #select features using VSURF
    # # tic()
    # set.seed(123)
    # vsurf = VSURF(y=combined_values_sample[,1], x = combined_values_sample[,2:ncol(combined_values_sample)], RFimplem='ranger', parallel=T)
    # # toc()
    # 
    # # Extract the variables selected by VSURF
    # vsurf.selected_indices <- vsurf$varselect.pred
    # vsurf.vars <- names(combined_scene)[vsurf.selected_indices]
    
    #subset sample dataframe to only have columns with selected features
    final_dat = combined_values_sample #|> select(all_of(c(met,vsurf.vars)))
    
    #get final ranger model using VSURF variables
    final_mod = ranger(y = final_dat[,1], x = final_dat[,2:ncol(final_dat), drop=F], importance='permutation')
    
    #save data
    input_vars = c(names(ps_scene)
                   # , names(glcm_r)
                   )
    variables = list(met, input_vars)
    names(variables) = c('Dependent.Var', 'Independent.Vars')
    
    r_l = list(final_mod
               # , vsurf
               , variables, combined_values_sample)
    names(r_l) = c('Final_model'
                   # , 'VSURF'
                   , 'Variables', 
                   'Training.Data')
    return(r_l)
  }
  
  
  #----predict lidar metric changes for each dataset----
  
  # print('Generating models')
  #train models
  ps_lid_change_mod_dir = paste0(ps_dir,'/lidar_change_metrics_models')
  dir.check(ps_lid_change_mod_dir)
  
  rf_vsurf_dir = paste0(ps_lid_change_mod_dir,'/randomforest_vsurf')
  dir.check(rf_vsurf_dir)
  
  feat_set = 'FeatureSet4_NoVarSel'
  
  mod_dir = paste0(rf_vsurf_dir,'/',feat_set,'_Pre=',pre_harv_id,'_Post=',post_harv_id)
  dir.check(mod_dir)
  
  #generate models by feature first
  pblapply(mets_to_test, function(m){
    pblapply(unique(results_meta_df$dataset), function(d){
      
      # print(paste(m,'---',d))
      
      ds_dir = paste0(mod_dir,'/',d)
      dir.check(ds_dir)
      
      f = paste0(ds_dir,'/',m,'.rds')
      cat('\n----Processing',f)
      tic()
      if(!file.exists(f)){
        
        #load appropriate spectral and texture mosaics
        ps_change = change_mosaics[[d]]
        # glcm_c = glcm_change[[d]]
        
        # #subset dataframe to get features appropriate for the dataset and metric
        # feats = lid_cor_summ$var_name[lid_cor_summ$dataset==d & lid_cor_summ$metric==m]
        # ps_change = ps_change[[feats]]
        
        
        #run modelling procedure, save output
        mod = metric_predictor_fn4(met = m, ps_scene = ps_change
                                   # , glcm_rast = glcm_c
                                   )
        write_rds(mod, f)
        
      }
      toc()
    })
  })
  
}, random_pre_harv_ids, random_post_harv_ids)
