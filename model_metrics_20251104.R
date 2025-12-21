source('derive lidar metrics_202501020.R')
library(ranger)
# source('scene_SeparabilityAnalysis_20250925.R')

#----function for assessing model performance----
perf_mets_calculator = function(pred=NULL, obs=NULL){
  
  #vector of performance metrics names
  Metric = c('R2', 'RMSE', 'RRMSE', 'ME', 'RME')
  
  #calculate performance metrics if both predictions and observations are provided
  if(!is.null(pred)&!is.null(obs)){
    
    residuals = obs - pred
    
    rss = sum(residuals^2)
    tss = sum((obs-mean(obs))^2)
    
    R2 = 1-(rss/tss) #coefficient of determination
    RMSE = sqrt(mean((obs-pred)^2)) #root mean squared error
    RRMSE = 100*RMSE/abs(mean(obs)) #relative rmse / percent rmse
    ME = mean(residuals) #mean error / bias
    RME = 100*ME/abs(mean(obs)) #relative mean error
    
    results = tibble(
      Metric = Metric,
      Value = c(R2, RMSE, RRMSE, ME, RME)
    )
    return(results)
  }
  
  #return error if predictions or observations are missing
  if(xor(is.null(pred),is.null(obs))){
    message("Error: need both predictions and observations to calculate metrics")
  }
  
  #return a vector of performance metrics names if neither predictions or observations are provided
  if(is.null(obs)&is.null(pred)){
    return(Metric)
  }
}

#----choose pre and post harvest dates----

#get ids which are common to each thinning block
ds_ids = sapply(unique(results_meta_df$dataset), function(d)results_meta_df$id[results_meta_df$dataset==d])
common_ids = Reduce(intersect, ds_ids)

common_ids <- results_meta_df %>%
  group_by(block_id, dataset) %>%
  summarise(ids = list(unique(id)), .groups = "drop") %>%
  summarise(common_ids = list(Reduce(intersect, ids)), .by = block_id)
ids_in_all_blocks <- Reduce(intersect, common_ids$common_ids)


#select pre-harvest and post-harvest scenes which are common to each block

pre_harv_ids = c("20210726_191611_86_2307", "20210729_191439_38_2274", "20210930_182644_85_2435")
post_harv_ids = c("20241011_192907_10_251a", "20240717_183454_97_24a7", "20240717_192626_64_24c6")

#make combinations of pre and post-harvest metrics to iterate over
grid = expand.grid(
  pre_harv_id = pre_harv_ids,
  post_harv_id = post_harv_ids
  # ,
  # dataset = unique(results_meta_df$dataset),
  # block_id = block_ids,
  # lidar_met = c('z_above2', 'z_p95', 'z_sd')
)

#define spectral features to use and lidar metrics to predict
base_feats = c('coastal_blue', 'blue', 'green_i', 'green', 'yellow', 'red', 'rededge', 'nir'
               ,'GNDVI','Hue', 'GCC','SAVI','NDVI','YNDVI')
mets_to_test = c('z_above2', 'z_p95', 'z_sd', 'z_cv', 'z_mean')

#set up output directories
ps_lid_change_mod_dir = paste0(ps_dir,'/lidar_change_metrics_models')
dir.check(ps_lid_change_mod_dir)

modtype_dir = paste0(ps_lid_change_mod_dir,'/randomforest')
dir.check(modtype_dir)

feat_set = 'BaseFeats_NoGLCM_NoVarsel' #identifier for the feature set used here


#loop over date pairs
pblapply(1:nrow(grid),function(i){
  
  gc()
  
  pre_harv_id = grid$pre_harv_id[i]
  post_harv_id = grid$post_harv_id[i]
  
  
  #----load data, make mosaics----
  
  cat('\nLoad data, make mosaics')
  
  #load lidar change metrics
  print('Load lidar change metrics')
  change_lids = lapply(block_ids, function(b){
    fl = list.files(lid_mets_change_cropped_dir, full.names = T)
    f = fl[str_detect(fl,b)]
    r = rast(f, lyrs=mets_to_test)
    return(r)
  })
  names(change_lids) = block_ids
  
  #generate a lidar change mosaic
  lid_change_mosaic = sprc(change_lids) |> merge()
  
  #load pre-harvest scenes
  print('Load pre-harvest scenes')
  pre_harvest_scenes = lapply(unique(results_meta_df$dataset), function(d){
    bl = lapply(block_ids, function(b){
      f = results_meta_df$file_path[results_meta_df$id == pre_harv_id &
                                      results_meta_df$block_id==b &
                                      results_meta_df$dataset==d][1]
      r = rast(f,lyrs=base_feats)
      
      r = r[[!names(r) %in% 'WDVI']] #remove WDVI since it's the same as DVI
      
      return(r)
    })
    names(bl) = block_ids
    return(bl)
  })
  names(pre_harvest_scenes) = unique(results_meta_df$dataset)
  
  #generate pre-harvest mosaics (i.e. combine thinning blocks into single raster stack)
  print('Generate mosaics')
  pre_harv_mosaics = lapply(pre_harvest_scenes, function(x){
    sc = sprc(x)
    r = merge(sc)
  })
  
  #load post-harvest scenes
  print('Load post harvest scenes')
  post_harvest_scenes = lapply(unique(results_meta_df$dataset), function(d){
    bl = lapply(block_ids, function(b){
      f = results_meta_df$file_path[results_meta_df$id == post_harv_id &
                                      results_meta_df$block_id==b &
                                      results_meta_df$dataset==d][1]
      r = rast(f, lyrs=base_feats)
      r = r[[!names(r) %in% 'WDVI']]
      return(r)
    })
    names(bl) = block_ids
    return(bl)
  })
  names(post_harvest_scenes) = unique(results_meta_df$dataset)
  
  print('Generate mosaics')
  post_harv_mosaics = lapply(post_harvest_scenes, function(x){
    sc = sprc(x)
    r = merge(sc)
  })
  
  #calculate change mosaics
  print('Generate change mosaics')
  change_mosaics = mapply(function(x, y){
    x-y
  }, pre_harv_mosaics, post_harv_mosaics, SIMPLIFY = F)
  
  
  #----write function for modelling lidar for a pre-harvest and post-harvest scene using variable selection and GLCM textures----
  
  
  metric_predictor_fn4 = function(met, ps_scene,
                                  lid=lid_change_mosaic
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
    combined_scene = c(lid_met, ps_scene)
    
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
    
    # # Extract the variables selected by VSURF
    # vsurf.selected_indices <- vsurf$varselect.pred
    # vsurf.vars <- names(combined_scene)[vsurf.selected_indices]
    
    #subset sample dataframe to only have columns with selected features
    final_dat = combined_values_sample #|> select(all_of(c(met,vsurf.vars)))
    
    #get final ranger model using VSURF variables
    final_mod = ranger(y = final_dat[,1], x = final_dat[,2:ncol(final_dat), drop=F], num.trees = 1000, importance='permutation')
    
    #get model performance metrics
    perf_mets = perf_mets_calculator(obs = final_dat[,1], pred = final_mod$predictions)
    
    #get variable importance
    vi = tibble(
      var_name = names(final_mod$variable.importance),
      rf_permutation_importance = final_mod$variable.importance
    )
    
    #save data
    input_vars = c(names(ps_scene))
    variables = list(met, input_vars)
    names(variables) = c('Dependent.Var', 'Independent.Vars')
    
    r_l = list(perf_mets,
               final_mod, 
               # vsurf, 
               variables, combined_values_sample)
    names(r_l) = c('Performance_metrics', 
                   'Final_model', 
                   # 'VSURF', 
                   'Variables', 'Training.Data')
    return(r_l)
  }
  
  
  #----predict lidar metric changes for each dataset----
  
  print('Generating models')
  # #train models
  # ps_lid_change_mod_dir = paste0(ps_dir,'/lidar_change_metrics_models')
  # dir.check(ps_lid_change_mod_dir)
  # 
  # modtype_dir = paste0(ps_lid_change_mod_dir,'/randomforest')
  # dir.check(modtype_dir)
  # 
  # feat_set = 'BaseFeats_NoGLCM_NoVarsel'
  # 
  mod_dir = paste0(modtype_dir,'/',feat_set,'_Pre=',pre_harv_id,'_Post=',post_harv_id)
  dir.check(mod_dir)
  
  #generate models by feature first
  pblapply(mets_to_test, function(m){
    pblapply(unique(results_meta_df$dataset), function(d){
      
      print(paste(m,'---',d))
      
      ds_dir = paste0(mod_dir,'/',d)
      dir.check(ds_dir)
      
      f = paste0(ds_dir,'/',m,'.rds')
      print(f)
      tic()
      if(!file.exists(f)){
        
        #load appropriate spectral and texture mosaics
        ps_change = change_mosaics[[d]]
        # glcm_c = glcm_change[[d]]
        
        # #subset dataframe to get features appropriate for the dataset and metric
        # feats = lid_cor_summ$var_name[lid_cor_summ$dataset==d & lid_cor_summ$metric==m]
        # ps_change = ps_change[[feats]]
        
        
        #run modelling procedure, save output
        mod = metric_predictor_fn4(met = m, ps_scene = ps_change #, glcm_rast = glcm_c
        )
        write_rds(mod, f)
        
      }
      toc()
    })
  })
})

#----extract model results----

mod_fl = list.files(modtype_dir, full.names = T, recursive = T, pattern='\\.tif$')

model_performance_f = file.path(modtype_dir,paste0('model_performance_summary.csv'))
variable_importance_f = file.path(modtype_dir,paste0('variable_importance_summary.csv'))

if(!(file.exists(model_performance_f)&file.exists(variable_importance_f))){
  
  mod_results_l = pblapply(mod_fl, function(f){
    obj = read_rds(f)
    
    dataset = basename(dirname(f))
    lid_met = basename(file_path_sans_ext(f))
    
    #get parameters
    harv_ids_dir = basename(dirname(dirname(f)))
    pre_id = str_extract(harv_ids_dir, "(?<=Pre=).*?(?=_Post)")
    post_id = str_extract(harv_ids_dir, '(?<=Post=).*')
    
    #model performance df
    {
      #load model performance metrics
      mod_perf_df = obj$Performance_metrics
      
      y_bar = mean(obj$Training.Data[,1])
      
      mod_perf_df = mod_perf_df |>
        #fix RRMSE and RME (scale them by absolute mean observations to avoid changing the sign)
        mutate(
          Value = case_match(Metric,
                             c('RRMSE','RME') ~ (Value*y_bar)/abs(y_bar),
                             # 'RME' ~ (Value*y_hat)/abs(y_hat),
                             .default = Value)
        ) |>
      #add out-of-bag r2
        cbind(tibble(Metric = 'OOB_R2', Value = obj$Final_model$r.squared))
      
      #add parameters to model performance dataframe
      mod_perf_df$dataset = dataset
      mod_perf_df$lidar_metric = lid_met
      mod_perf_df$pre_harvest_id = pre_id
      mod_perf_df$post_harvest_id = post_id
    }
    
    #variable importance df
    {
      vi_tbl = tibble(
        var_name = names(obj$Final_model$variable.importance),
        permutation_importance = obj$Final_model$variable.importance,
        dataset = dataset,
        lidar_metric = lid_met,
        pre_harvest_id = pre_id,
        post_harvest_id = post_id
      )
    }
    res_l = list(mod_perf_df, vi_tbl)
    names(res_l) = c('Model_performance', 'Variable_importance')
    return(res_l)
  })
  
  if(!file.exists(model_performance_f)){
    model_perf_df = lapply(mod_results_l, function(x)x[[1]]) |> bind_rows()
    write_csv(model_perf_df, file = model_performance_f)
  }
  
  if(!file.exists(variable_importance_f)){
    var_imp_df = lapply(mod_results_l, function(x)x[[2]]) |> bind_rows()
    write_csv(var_imp_df, file = variable_importance_f)
  }
  
} else {
  
  model_perf_df = read_csv(model_performance_f)
  var_imp_df = read_csv(variable_importance_f)
}

#----plot results----

#code the date combinations
model_perf_df = model_perf_df |> 
  mutate(id_combos = paste0(pre_harvest_id,' - ',post_harvest_id))

#do two way Kruskal Wall predicting R2 

#plot model performance metrics
mod_perf_plot_l = lapply(unique(model_perf_df$Metric), function(met){
  ggplot(model_perf_df |> filter(Metric==met))+
    geom_boxplot(aes(x = dataset, y = Value), outliers=F) +
    geom_jitter(aes(x = dataset, y = Value, color=id_combos), width=0.1, height=0, size=2, alpha = 0.5) +
    facet_grid(rows=vars(lidar_metric), scales='free')+
    theme_bw()
})
names(mod_perf_plot_l) = unique(mod_perf_df$Metric)

# x11()
mod_perf_plot_l$R2

mod_perf_plot_l$RRMSE

#variable importance plot
var_imp_df = var_imp_df |>
  mutate(id_combos = paste0(pre_harvest_id,' - ',post_harvest_id))


var_imp_p = ggplot(var_imp_df)+
  geom_boxplot(aes(x = var_name, y = permutation_importance), outliers=F) +
  geom_jitter(aes(x = var_name, y = permutation_importance, color=id_combos), width=0.1, height=0, size=1, alpha = 0.5) +
  coord_flip()+
  facet_grid(rows=vars(dataset),cols = vars(lidar_metric),scales='free')+
  theme_bw() #+
  # theme(axis.text.x = element_text(angle=90))
var_imp_p

mean_var_imp_df = var_imp_df |>
  group_by(var_name,dataset,lidar_metric) |>
  summarise(mean_importance = mean(permutation_importance))

