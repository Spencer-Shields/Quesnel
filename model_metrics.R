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


metric_predictor_fn = function(met, ps_scene,
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
  
  #stack lidar metric with change scene
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
  
  #select spectral features using VSURF
  tic()
  set.seed(123)
  spectral.vsurf = VSURF(y=combined_values_sample[[met]], x = combined_values_sample[names(ps_scene)], RFimplem='ranger', parallel=T)
  toc()
  
  # Get the selected variables from the prediction step (final set)
  spectral.vsurf.selected_indices <- spectral.vsurf$varselect.pred
  spectral.vsurf.vars <- names(ps_scene)[spectral.vsurf.selected_indices]
  
  #remove spectral layers eliminated by vsurf
  combined_scene_vsurf = combined_scene[[c(met,spectral.vsurf.vars)]]
  
  #make glcm_textures for the selected variables
  if(is.null(glcm_rast)){
  glcm_r = pblapply(spectral.vsurf.vars, function(x){
    r = combined_scene[[x]]
    g = glcm_textures(r, w = c(glcm_window, glcm_window), n_levels = glcm_quant_lvls, quant_method = glcm_quant_method, metrics = glcm_mets)
    names(g) = paste0(x,'_',names(g))
    return(g)
  }) |> rast()
  } else {
    glcm_lyrs_wanted = as.vector(outer(spectral.vsurf.vars, glcm_mets, FUN = function(x, y) paste0(x, "_", y)))
    glcm_r = glcm_rast[[glcm_lyrs_wanted]]
  }

  #combine textures with other variables, run vsurf again
  combined_scene_vs_glcm = c(combined_scene_vsurf, glcm_r)
  
  set.seed(123)
  combined_values_glcm_sample = lapply(block_ids, function(b){
    block = blocks_p[blocks_p$BLOCKNUM==b,]
    r_b = crop(combined_scene_vs_glcm, block, mask=T)
    s = spatSample(x=r_b, size=min(blocks_ncells_summ$n), method='random', na.rm=T, exhaustive=T)
    return(s)
  }) |> bind_rows()
  
  spect_text.vsurf = VSURF(y=combined_values_glcm_sample[[met]], x = combined_values_glcm_sample|>select(!all_of(met)), RFimplem='ranger', parallel=T)
  
  # Get the selected variables from the prediction step (final set)
  spect_text.vsurf.selected_indices <- spect_text.vsurf$varselect.pred
  spect_text.vsurf.vars <- names(combined_values_glcm_sample)[spect_text.vsurf.selected_indices]
  
  #get final dataset
  final_dat = combined_values_glcm_sample |> select(all_of(c(met,spect_text.vsurf.vars)))
  
  #get final ranger model
  final_mod = ranger(y = final_dat[,1], x = final_dat[,2:ncol(final_dat)], importance='permutation')
  
  #output results
  r_all = c(combined_scene, glcm_r)
  
  r_l = list(final_mod, spect_text.vsurf, spectral.vsurf)
  names(r_l) = c('Final_model', 'Spectral_VSURF', 'Spect.Text_VSURF')
  return(r_l)
  }


#----predict lidar metric changes for each dataset----

#choose metrics to test
mets_to_test = c( 
  "z_max",    #"z_min"    "z_p5"     "z_p10"    "z_p15"    "z_p20"    "z_p25"    "z_p30"    "z_p35"    "z_p40"   
                 # "z_p45"    "z_p50"    "z_p55"    "z_p60"    "z_p65"    "z_p70"    "z_p75"    "z_p80"    "z_p85"    
                 "z_p90",   
                     #"z_mean",
                 "z_above1", "z_above2", "z_sd","z_p95", "z_above3", "z_above4", "z_above5"
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
  
  n_clusts = 4
  
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
  
  # lid_cor_varstouse = lid_cor_summ |> 
  #   filter(var_rank==1
  #          # ,abs_cor>=0.05
  #          ) |>
  #   pivot_wider(names_from = metric, values_from = abs_cor)
    
  print(unique(lid_cor_varstouse$var_name))
  
  ggplot(lid_cor_varstouse)+
    geom_histogram(aes(abs_cor))+
    facet_grid(rows=vars(dataset), cols = vars(metric))
}

#choose PlanetScope features to use in the initial model
{
  spec_feats_1 = c('NDVI', 'GEMI', 
                   'GCC', 'GLI', 'Hue',
                   'NRVI', 'NDGR', 'SR', 'VARIgreen',
                   'coastal_blue', 'blue', 'green_i',
                   'green', 'yellow', 'red', 'rededge', 'nir')
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

# mod_dir = paste0(rf_vsurf_dir,'/',
#   str_c(feats_to_test,collapse = '_'),'_', 
#                  str_c(glcm_mets, collapse = '_'),'_',
#                  glcm_quant_lvls,'_',
#                  glcm_window,'_',
#                  glcm_quant_method,'_',
#   'PreID=',pre_harv_id,'_',
#   'PostID=', post_harv_id)
# dir.check(mod_dir)

# mod_parameters = paste0(str_c(feats_to_test,collapse = '_'),'_', 
#                     str_c(glcm_mets, collapse = '_'),'_',
#                     glcm_quant_lvls,'_',
#                     glcm_window,'_',
#                     glcm_quant_method)


spec_feats_1 = c('NDVI', 'GEMI', 
                 'GCC', 'GLI', 'Hue',
                 # 'NRVI', 
                 'NDGR', 'SR', 'VARIgreen',
                 'coastal_blue', 'blue', 'green_i',
                 'green', 'yellow', 'red', 'rededge', 'nir')
feat_set = 'FeatureSet1'

mod_dir = paste0(rf_vsurf_dir,'/',feat_set,'_Pre=',pre_harv_id,'_Post=',post_harv_id)
dir.check(mod_dir)

feats_to_test = spec_feats_1

#subset change mosaics to only have target set of spectral features
change_mosaics = lapply(change_mosaics, function(r)r[[feats_to_test]])

#pre-make glcm layers to save on processing time
glcm_pre = pblapply(pre_harv_mosaics, function(r){
  r_ = r[[feats_to_test]]
  gl = pblapply(names(r_), function(x){
    r = r_[[x]]
    g = glcm_textures(r, w = c(glcm_window, glcm_window), n_levels = glcm_quant_lvls, quant_method = glcm_quant_method, metrics = glcm_mets)
    names(g) = paste0(x,'_',names(g))
    return(g)
}) |> rast()
 return(gl)   
})
  
glcm_post = pblapply(post_harv_mosaics, function(r){
  r_ = r[[feats_to_test]]
  gl = pblapply(names(r_), function(x){
    r = r_[[x]]
    g = glcm_textures(r, w = c(glcm_window, glcm_window), n_levels = glcm_quant_lvls, quant_method = glcm_quant_method, metrics = glcm_mets)
    names(g) = paste0(x,'_',names(g))
    return(g)
  }) |> rast()
  return(gl)  
})

glcm_change = mapply(function(x,y){x-y}, glcm_pre, glcm_post)
names(glcm_change) = names(pre_harv_mosaics)
print(glcm_change)

##generate models by dataset first
{
# ps_lid_change_mods = pblapply(unique(results_meta_df$dataset), function(d){
#   
#   ds_dir = paste0(mod_dir,'/',d)
#   dir.check(ds_dir)
#   
#   ps_change = change_mosaics[[d]]
#   # ps_change = ps_change[[feats_to_test]]
#   ml = pblapply(mets_to_test, function(m){
#     
#     f = paste0(ds_dir,'/',m,'.rds')
#     print(f)
#     tic()
#     if(!file.exists(f)){
#     mod = metric_predictor_fn(met = m, ps_scene = ps_change)
#     write_rds(mod, f)
#     } else {
#       mod = read_rds(f)
#     }
#     toc()
#     return(mod)
#   })
#   names(ml) = mets_to_test
#   return(ml)
# })
# names(ps_lid_change_mods) = unique(results_meta_df$dataset)
}

#generate models by feature first
pblapply(mets_to_test, function(m){
  ml = pblapply(unique(results_meta_df$dataset), function(d){
    
    ds_dir = paste0(mod_dir,'/',d)
    dir.check(ds_dir)
    
    f = paste0(ds_dir,'/',m,'.rds')
    print(f)
    tic()
    Sys.time()
    if(!file.exists(f)){
      ps_change = change_mosaics[[d]]
      glcm_c = glcm_change[[d]] 
      mod = metric_predictor_fn(met = m, ps_scene = ps_change, glcm_rast = glcm_c)
      write_rds(mod, f)
    } else {
      mod = read_rds(f)
    }
    toc()
    return(mod)
  })
})




#----choose spectral features to compare----
{
# #get average AUC and JMD for every combination of id, block, and dataset
# vars_summ = results_meta_df |>
#   select(-block_pixel_stratum)|> distinct() |>
#   group_by(id, dataset, var_name) |>
#   summarise(mean_jmd = mean(jmd_scaled),
#             mean_logauc = mean(logistic_AUC))
# 
# #pre-harvest
# {
# pre_harv_vars_summ = vars_summ[vars_summ$id==pre_harv_id,] #|>
#   # pivot_longer(cols = c('mean_jmd', 'mean_logauc'), names_to = 'metric', values_to = 'val')
# 
# ggplot(pre_harv_vars_summ)+
#   geom_col(aes(x = reorder(var_name, -mean_jmd), y = mean_jmd)) +
#   facet_wrap(vars(dataset), nrow = length(unique(pre_harv_vars_summ$dataset)), scales='free')
# 
# #load pre-harvest scenes
# pre_harvest_scenes = lapply(unique(results_meta_df$dataset), function(d){
#   bl = lapply(block_ids, function(b){
#     f = results_meta_df$file_path[results_meta_df$id == pre_harv_id &
#                                     results_meta_df$block_id==b &
#                                     results_meta_df$dataset==d][1]
#     r = rast(f)
#     return(r)
#   })
#   names(bl) = block_ids
#   return(bl)
# })
# names(pre_harvest_scenes) = unique(results_meta_df$dataset)
# 
# pre_harv_mosaics = pblapply(pre_harvest_scenes, function(x){
#   sc = sprc(x)
#   r = merge(sc)
# })
# 
# ggplot()+geom_spatraster(data=pre_harv_mosaics[[1]], aes(fill=SR))
# 
# pre_harv_mosaics_df_l = pblapply(pre_harv_mosaics, function(x){v = values(x, na.rm=T)})
# 
# pre_harv_mosaics_pca = lapply(pre_harv_mosaics_df_l, function(x)princomp(scale(x)))
# 
# corrplot(cor(pre_harv_mosaics_df_l[[1]]), order='hclust')
#   }
# 
# #post-harvest
# {
#   post_harv_vars_summ = vars_summ[vars_summ$id==post_harv_id,] #|>
#   # pivot_longer(cols = c('mean_jmd', 'mean_logauc'), names_to = 'metric', values_to = 'val')
#   
#   ggplot(post_harv_vars_summ)+
#     geom_col(aes(x = reorder(var_name, -mean_jmd), y = mean_jmd)) +
#     facet_wrap(vars(dataset), nrow = length(unique(post_harv_vars_summ$dataset)), scales='free')
#   
#   #load post-harvest scenes
#   post_harvest_scenes = lapply(unique(results_meta_df$dataset), function(d){
#     bl = lapply(block_ids, function(b){
#       f = results_meta_df$file_path[results_meta_df$id == post_harv_id &
#                                       results_meta_df$block_id==b &
#                                       results_meta_df$dataset==d][1]
#       r = rast(f)
#       return(r)
#     })
#     names(bl) = block_ids
#     return(bl)
#   })
#   names(post_harvest_scenes) = unique(results_meta_df$dataset)
#   
#   post_harv_mosaics = pblapply(post_harvest_scenes, function(x){
#     sc = sprc(x)
#     r = merge(sc)
#   })
#   
#   ggplot()+geom_spatraster(data=post_harv_mosaics[[1]], aes(fill=SR))
#   
#   post_harv_mosaics_df_l = pblapply(post_harv_mosaics, function(x){v = values(x, na.rm=T)})
#   
#   post_harv_mosaics_pca = lapply(post_harv_mosaics_df_l, function(x)princomp(scale(x)))
#   
#   corrplot(cor(post_harv_mosaics_df_l[[1]]), order='hclust')
# }

}
# 
# 
# #load pre-harvest scenes
# pre_harvest_scenes = lapply(unique(results_meta_df$dataset), function(d){
#   bl = lapply(block_ids, function(b){
#     f = results_meta_df$file_path[results_meta_df$id == pre_harv_id &
#                                     results_meta_df$block_id==b &
#                                     results_meta_df$dataset==d][1]
#     r = rast(f)
# 
#     r = r[[!names(r) %in% 'WDVI']] #remove WDVI since it's the same as DVI
# 
#     return(r)
#   })
#   names(bl) = block_ids
#   return(bl)
# })
# names(pre_harvest_scenes) = unique(results_meta_df$dataset)
# 
# pre_harv_mosaics = pblapply(pre_harvest_scenes, function(x){
#   sc = sprc(x)
#   r = merge(sc)
# })
# 
# #load post-harvest scenes
# post_harvest_scenes = lapply(unique(results_meta_df$dataset), function(d){
#   bl = lapply(block_ids, function(b){
#     f = results_meta_df$file_path[results_meta_df$id == post_harv_id &
#                                     results_meta_df$block_id==b &
#                                     results_meta_df$dataset==d][1]
#     r = rast(f)
#     r = r[[!names(r) %in% 'WDVI']]
#     return(r)
#   })
#   names(bl) = block_ids
#   return(bl)
# })
# names(post_harvest_scenes) = unique(results_meta_df$dataset)
# 
# post_harv_mosaics = pblapply(post_harvest_scenes, function(x){
#   sc = sprc(x)
#   r = merge(sc)
# })
# 
# #calculate change mosaics
# change_mosaics = pbmapply(function(x, y){
#   x-y
# }, pre_harv_mosaics, post_harv_mosaics, SIMPLIFY = F)
# 
# #change mosaics PCA
# change_mosaics_pca = pblapply(change_mosaics, function(r){
#   v = values(r, na.rm=T)
#   vs = scale(v)
#   p = princomp(vs)
#   return(p)
# })
# 
# #get variance explained for individual PCs
# cm_pc_var_l = pblapply(change_mosaics_pca, function(p){
#   variance_explained = (p$sdev^2)/sum(p$sdev^2) * 100
# }) #first 4 PCs explain most variance
# 
# #extract dataframe of absolute loadings for first 4 PCs
# n = 4
# 
# cm_pc_load_df = pblapply(1:length(change_mosaics_pca), function(i){
#   p = change_mosaics_pca[[i]]
#   ld = p$loadings[,1:n]
#   ld = abs(ld) #use absolute loadings
#   ldf = as.data.frame(ld)
#   ldf$var_name = rownames(ldf)
#   ldf$dataset = names(change_mosaics_pca)[i]
#   return(ldf)
# }) |>
#   bind_rows() |>
#   pivot_longer(cols = paste0('Comp.', 1:n), names_to = 'PC', values_to = 'loading')%>%
#   group_by(dataset, PC) %>%           # group by facets
#   mutate(var_name_f = factor(var_name, levels = var_name[order(-loading)])) %>%
#   ungroup()
# 
# pc_plot_l = mapply(function(d,pc){
#   df = cm_pc_load_df |> filter(dataset==d & PC==pc)
#   gg = ggplot(df)+
#     geom_col(aes(x = reorder(var_name, loading, decreasing=T), y = loading))+
#     ggtitle(paste0(d,', ',pc))+
#     theme(axis.text = element_text(angle = 90))
#   return(gg)
# }, unique(cm_pc_load_df$dataset), unique(cm_pc_load_df$PC), SIMPLIFY = F)
# 
# #plot PC loadings
# library(patchwork)
# pc_plot_l = lapply(unique(cm_pc_load_df$dataset), function(d){
#   pc_l = lapply(unique(cm_pc_load_df$PC), function(pc){
#     df = cm_pc_load_df |> filter(dataset==d & PC==pc)
#     gg = ggplot(df)+
#       geom_point(aes(x = reorder(var_name, loading, decreasing=T), y = loading))+
#       ggtitle(paste0(d,', ',pc))+
#       theme(axis.text = element_text(angle = 90))
#     return(gg)
#   })
#   # names(pc_l) = unique(cm_pc_load_df$PC)
#   # return(pc_l)
#   combined_plot = patchwork::wrap_plots(pc_l, ncol=length(pc_l))
#   return(combined_plot)
# })
# pc_combined_plot = wrap_plots(pc_plot_l, nrow=length(pc_plot_l))
# pc_combined_plot
# 
# #check correlation matrix
# corrplot(cor(values(change_mosaics[[1]],na.rm=T)), order='hclust')
# 
# #get mean AUC values
# auc_mean_coefs_byfeat_df = auc_coefs |>
#   group_by(var_name, model_terms) |>
#   summarise(mean_Estimate = mean(Estimate))
