source('derive lidar metrics_202501001.R')
# source('scene_SeparabilityAnalysis_20250925.R')

#----choose pre and post harvest dates----

#get ids which are common to each thinning block
ds_ids = sapply(unique(results_meta_df$dataset), function(d)results_meta_df$id[results_meta_df$dataset==d])
common_ids = Reduce(intersect, ds_ids)

results_meta_df |> select(id, file_path, dataset, block_id) |> distinct() |> filter(id == "20241011_192907_10_251a")

pre_harv_date = '20210713'
post_harv_date = '20240818'

pre_harv_id = "20210713_191755_66_2414"
post_harv_id = "20241011_192907_10_251a"


#----load data, make mosaics----

#load lidar change metrics
change_lids = lapply(block_ids, function(b){
  fl = list.files(lid_mets_change_cropped_dir, full.names = T)
  f = fl[str_detect(fl,b)]
  r = rast(f)
  return(r)
})
names(change_lids) = block_ids

#generate a lidar change mosaic
lid_change_mosaic = sprc(change_lids) |> merge()

#load pre-harvest scenes
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

pre_harv_mosaics = pblapply(pre_harvest_scenes, function(x){
  sc = sprc(x)
  r = merge(sc)
})

#load post-harvest scenes
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

post_harv_mosaics = pblapply(post_harvest_scenes, function(x){
  sc = sprc(x)
  r = merge(sc)
})

#calculate change mosaics
change_mosaics = pbmapply(function(x, y){
  x-y
}, pre_harv_mosaics, post_harv_mosaics, SIMPLIFY = F)






#----write function for modelling lidar for a pre-harvest and post-harvest scene using variable selection and GLCM textures----

require(GLCMTextures)
require(VSURF)

# # sample values for debugging
# {met = 'z_above2';  ps_scene = pre_harv_mosaics$`Non-normalized`-post_harv_mosaics$`Non-normalized`; lid=lid_change_mosaic}

glcm_mets = c(
  "glcm_contrast",
  # "glcm_dissimilarity",
  #  "glcm_homogeneity", "glcm_ASM",
  "glcm_entropy", 
  "glcm_mean"
  # , "glcm_variance", "glcm_correlation"
)
glcm_quant_lvls = 32
glcm_quant_method = 'range'
glcm_window = 3


metric_predictor_fn2 = function(met, ps_scene,
                               lid=lid_change_mosaic,
                               glcm_mets = c(
                                 "glcm_contrast",
                                 # "glcm_dissimilarity",
                                 #  "glcm_homogeneity", "glcm_ASM",
                                 "glcm_entropy", 
                                 "glcm_mean"
                                 # , "glcm_variance", "glcm_correlation"
                               ),
                               glcm_rast = NULL, #spatraster with all possible glcm layers, leave as NULL to make this on-the-fly using ps_scene
                               glcm_quant_lvls = 32,
                               glcm_quant_method = 'range',
                               glcm_window = 3){
  
  
  #get lidar change metric
  lid_met = lid[[met]]
  
  #make glcm_textures for the selected variables
  if(is.null(glcm_rast)){
    glcm_r = pblapply(names(ps_scene), function(x){
      r = ps_scene[[x]]
      g = glcm_textures(r, w = c(glcm_window, glcm_window), n_levels = glcm_quant_lvls, quant_method = glcm_quant_method, metrics = glcm_mets)
      names(g) = paste0(x,'_',names(g))
      return(g)
    }) |> rast()
  } else {
    glcm_lyrs_wanted = as.vector(outer(names(ps_scene), glcm_mets, FUN = function(x, y) paste0(x, "_", y))) #subset the premade raster to only have the desired layers
    glcm_r = glcm_rast[[glcm_lyrs_wanted]]
  }
  
  #stack lidar metric, planetscope scene, and planetscope glcm layers
  if(compareGeom(lid_met, ps_scene, stopOnError=F)==F){ps_scene = resample(ps_scene_scene, lid_met)}
  if(compareGeom(lid_met, glcm_r, stopOnError=F)==F){glcm_r = resample(glcm_r, lid_met)}
  combined_scene = c(lid_met, ps_scene, glcm_r)
  
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
  
  #select features using VSURF
  # tic()
  set.seed(123)
  vsurf = VSURF(y=combined_values_sample[,1], x = combined_values_sample[,2:ncol(combined_values_sample)], RFimplem='ranger', parallel=T)
  # toc()
  
  # Extract the variables selected by VSURF
  vsurf.selected_indices <- vsurf$varselect.pred
  vsurf.vars <- names(combined_scene)[vsurf.selected_indices]
  
  #subset sample dataframe to only have columns with selected features
  final_dat = combined_values_sample |> select(all_of(c(met,vsurf.vars)))
  
  #get final ranger model using VSURF variables
  final_mod = ranger(y = final_dat[,1], x = final_dat[,2:ncol(final_dat)], importance='permutation')
  
  #save data
  input_vars = c(names(ps_scene), names(glcm_r))
  variables = list(met, input_vars)
  names(variables) = c('Dependent.Var', 'Independent.Vars')
  
  r_l = list(final_mod, vsurf, variables, combined_values_sample)
  names(r_l) = c('Final_model', 'VSURF', 'Variables', 'Training.Data')
  return(r_l)
  }


#----predict lidar metric changes for each dataset----

#choose metrics to test
mets_to_test = c( 
  "z_above2","z_p95","z_sd","z_max",    #"z_min"    "z_p5"     "z_p10"    "z_p15"    "z_p20"    "z_p25"    "z_p30"    "z_p35"    "z_p40"   
                 # "z_p45"    "z_p50"    "z_p55"    "z_p60"    "z_p65"    "z_p70"    "z_p75"    "z_p80"    "z_p85"    
                 "z_p90",   
                 #    #"z_mean",   
                 "z_above1", "z_above3", "z_above4", "z_above5"
                 # , "z_cv" #,"count"
                 )

#choose spectral features
{
  library(ClustOfVar)
  
  #combine lidar with spectral features
  combined_ps_lid_mosaics = lapply(change_mosaics, function(x)c(lid_change_mosaic[[mets_to_test]], x))
  
  #do heirarchical clustering of variables based on correlation matrics
  ps_cors = pblapply(change_mosaics, function(x)cor(values(x,na.rm=T)))
  var_hclusts = pblapply(ps_cors, hclustvar2)
  
  n_clusts = 8
  
  var_cutrees = pblapply(var_hclusts, function(x){
    c=cutree(x, k=n_clusts)}
    )
  # select_vars_l = pblapply(var_cutrees, function(x){
  #   p = x$
  # })
  
  #make dataframes of correlations with lidar metrics
  combined_ps_lid_dfs = pblapply(combined_ps_lid_mosaics, function(x)values(x,na.rm=T))
  combined_ps_lid_cors= pblapply(combined_ps_lid_dfs, cor)
  lid_cor_df = pblapply(names(combined_ps_lid_cors), function(x){
    corr = combined_ps_lid_cors[[x]]
    corr = abs(corr)
    d = as.data.frame(corr) |> 
      select(all_of(mets_to_test)) |>
      mutate(var_name = rownames(corr), dataset = x) |>
      filter(!var_name %in% mets_to_test)
    
    cut = var_cutrees[[x]]
    cd = data.frame(var_name = names(cut), cluster = cut)
    
    d = d |> left_join(cd)
  }) |> bind_rows()
  
  lid_cor_df_long = lid_cor_df |> pivot_longer(cols = all_of(mets_to_test), values_to = 'abs_cor', names_to = 'metric')
  
  lid_cor_summ = lid_cor_df_long |>
    group_by(dataset, metric, cluster) |>
    arrange(desc(abs_cor), .by_group = T) |>
    mutate(var_rank = row_number(desc(abs_cor))) |>
    ungroup() |>
    filter(var_rank==1
                     # ,abs_cor>=0.05
    )
  
  vars_summ = lid_cor_summ |>
    group_by(var_name,metric)|>
    summarise(mean_abs_cor = mean(abs_cor)) |>
    group_by(metric)|>
    arrange(desc(mean_abs_cor), .by_group = T)
  
  # lid_cor_varstouse = lid_cor_summ |> 
  #   filter(var_rank==1
  #          # ,abs_cor>=0.05
  #          ) |>
  #   pivot_wider(names_from = metric, values_from = abs_cor)
  
  all_spec_vars = unique(lid_cor_summ$var_name) 
    
  print(unique(lid_cor_varstouse$var_name))
  
  # ggplot(lid_cor_varstouse)+
  #   geom_histogram(aes(abs_cor))+
  #   facet_grid(rows=vars(dataset), cols = vars(metric))
}

#choose PlanetScope features to use in the initial model
{
  # spec_feats_1 = c('NDVI', 'GEMI', 
  #                  'GCC', 'GLI', 'Hue',
  #                  'NRVI', 'NDGR', 'SR', 'VARIgreen',
  #                  'coastal_blue', 'blue', 'green_i',
  #                  'green', 'yellow', 'red', 'rededge', 'nir')
}


glcm_mets = c(
  "glcm_contrast",
  # "glcm_dissimilarity",
  #  "glcm_homogeneity", "glcm_ASM",
  "glcm_entropy", 
  "glcm_mean"
  # , "glcm_variance", "glcm_correlation"
)
glcm_quant_lvls = 64
glcm_quant_method = 'range'
glcm_window = 3


#train models
ps_lid_change_mod_dir = paste0(ps_dir,'/lidar_change_metrics_models')
dir.check(ps_lid_change_mod_dir)

rf_vsurf_dir = paste0(ps_lid_change_mod_dir,'/randomforest_vsurf')
dir.check(rf_vsurf_dir)

spec_feats_1 = c('NDVI', 'GEMI', 
                 'GCC', 'GLI', 'Hue',
                 # 'NRVI', 
                 'NDGR', 'SR', 'VARIgreen',
                 'coastal_blue', 'blue', 'green_i',
                 'green', 'yellow', 'red', 'rededge', 'nir')

feat_set = 'FeatureSet2_1stepVSURF'

mod_dir = paste0(rf_vsurf_dir,'/',feat_set,'_Pre=',pre_harv_id,'_Post=',post_harv_id)
dir.check(mod_dir)

# feats_to_test = spec_feats_1
feats_to_test = names(pre_harv_mosaics[[1]])

#subset change mosaics to only have target set of spectral features
change_mosaics = lapply(change_mosaics, function(r)r[[feats_to_test]])

#pre-make glcm layers to save on processing time
glcm_dir = paste0(ps_dir, '/GLCM_','quant=',glcm_quant_lvls,'_quantmethod=',glcm_quant_method,'_w=',glcm_window,'_mets=',str_c(str_replace_all(glcm_mets,'glcm_',''),collapse=','))
dir.check(glcm_dir)

glcm_pre_dir = paste0(glcm_dir,'/',pre_harv_id)
dir.check(glcm_pre_dir)

glcm_pre = pblapply(names(pre_harv_mosaics), function(nom){
  f = paste0(glcm_pre_dir,'/',nom,'.tif')
  if(!file.exists(f)){
    r_ = pre_harv_mosaics[[nom]]
    r_ = r_[[feats_to_test]]
    gl = pblapply(names(r_), function(x){
      s = r_[[x]]
      g = glcm_textures(s, w = c(glcm_window, glcm_window), n_levels = glcm_quant_lvls, quant_method = glcm_quant_method, metrics = glcm_mets)
      names(g) = paste0(x,'_',names(g))
      return(g)
    }) |> rast()
    
    writeRaster(gl, f)
  } else {
    gl = rast(f)
  }
  return(gl)   
})
names(glcm_pre) = names(pre_harv_mosaics)

glcm_post_dir = paste0(glcm_dir,'/',post_harv_id)
dir.check(glcm_post_dir)

glcm_post = pblapply(names(post_harv_mosaics), function(nom){
  f = paste0(glcm_post_dir,'/',nom,'.tif')
  if(!file.exists(f)){
    r_ = unwrap(r_w)
    r_ = r_[[feats_to_test]]
    gl = future_lapply(names(r_), function(x){
      s = unwrap(r_[[x]])
      g = glcm_textures(s, w = c(glcm_window, glcm_window), n_levels = glcm_quant_lvls, quant_method = glcm_quant_method, metrics = glcm_mets)
      names(g) = paste0(x,'_',names(g))
      return(g)
    }) |> rast()
    writeRaster(gl, f)
  } else {
    gl = rast(f)
  }
  return(gl)   
})
names(glcm_post) = names(post_harv_mosaics)

# glcm_change_dir = paste0(glcm_dir,'/change',pre_harv_id,'to',post_harv_id)
# dir.check(glcm_change_dir)

glcm_change = mapply(function(x,y){x-y}, glcm_pre, glcm_post)
names(glcm_change) = names(pre_harv_mosaics)
print(glcm_change)

#generate models by feature first
pblapply(mets_to_test, function(m){
  pblapply(unique(results_meta_df$dataset), function(d){
    
    ds_dir = paste0(mod_dir,'/',d)
    dir.check(ds_dir)
    
    f = paste0(ds_dir,'/',m,'.rds')
    print(f)
    tic()
    if(!file.exists(f)){
      
      tic()
      
      #load appropriate spectral and texture mosaics
      ps_change = change_mosaics[[d]]
      glcm_c = glcm_change[[d]]
      
      #subset dataframe to get features appropriate for the dataset and metric
      feats = lid_cor_summ$var_name[lid_cor_summ$dataset==d & lid_cor_summ$metric==m]
      ps_change = ps_change[[feats]]
      
      
      #run modelling procedure, save output
      mod = metric_predictor_fn2(met = m, ps_scene = ps_change, glcm_rast = glcm_c)
      write_rds(mod, f)
      
      toc()
    }
  })
})

