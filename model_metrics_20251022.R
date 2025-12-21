source('derive lidar metrics_202501020.R')
library(terra)
library(tidyverse)
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

pre_harv_id = "20210726_191611_86_2307"
post_harv_id = "20241011_192907_10_251a"

#----load lidar metrics----

mets_to_test = c('z_above2', 'z_p95', 'z_sd')

#load pre-harvest lidar change metrics
print('Load pre-harvest lidar change metrics')
pre_lids = lapply(block_ids, function(b){
  fl = list.files(lid_mets_pre_cropped_dir, full.names = T)
  f = fl[str_detect(fl,b)]
  r = rast(f)
  r = r[[mets_to_test]]
  return(r)
})
names(pre_lids) = block_ids

#generate a pre-harvest lidar change mosaic
lid_pre_mosaic = sprc(pre_lids) |> merge()

#load post-harvest lidar change metrics
print('Load post-harvest lidar change metrics')
post_lids = lapply(block_ids, function(b){
  fl = list.files(lid_mets_post_dec_cropped_dir, full.names = T)
  f = fl[str_detect(fl,b)]
  r = rast(f)
  r = r[[mets_to_test]]
  return(r)
})
names(post_lids) = block_ids

#generate a post-harvest lidar change mosaic
lid_post_mosaic = sprc(post_lids) |> merge()

#----load planetscope data----

pre_ps_byblock_l = pblapply(unique(results_meta_df$dataset), function(d){
  bl = pblapply(unique(results_meta_df$block_id), function(b){
    dat = results_meta_df |> 
      filter(id == pre_harv_id, block_id == b, dataset == d) |>
      select(id, file_path, block_id, dataset) |>
      distinct()
    r = rast(dat$file_path)
    r = r[[unique(results_meta_df$var_name)]] #only keep features in results_meta_df (help with multicollinearity)
    return(r)
  })
  names(bl) = unique(results_meta_df$block_id)
  return(bl)
})
names(pre_ps_byblock_l) = unique(results_meta_df$dataset)

pre_ps_mosaics_l = pblapply(pre_ps_byblock_l,function(x)sprc(x)|>merge())


post_ps_byblock_l = pblapply(unique(results_meta_df$dataset), function(d){
  bl = pblapply(unique(results_meta_df$block_id), function(b){
    dat = results_meta_df |> 
      filter(id == post_harv_id, block_id == b, dataset == d) |>
      select(id, file_path, block_id, dataset) |>
      distinct()
    r = rast(dat$file_path)
    r = r[[unique(results_meta_df$var_name)]] #only keep features in results_meta_df (help with multicollinearity)
    return(r)
  })
  names(bl) = unique(results_meta_df$block_id)
  return(bl)
})
names(post_ps_byblock_l) = unique(results_meta_df$dataset)

post_ps_mosaics_l = pblapply(post_ps_byblock_l,function(x)sprc(x)|>merge())



#----generate stratified samples for each lidar metric----


#define function to stratify a raster using histogram bins
strat_histogram = function(mraster, breaks = 10){
  #stratify a raster using percentiles
  h = hist(mraster,breaks=breaks,plot=F)
  breaks_ = h$breaks
  
  sraster = terra::classify(mraster,breaks_,include.lowest=T)
  names(sraster) = 'strata'
  return(sraster)
}
# 
# library(sgsR)
# 
# #sample lid metrics by block and metric
# lid_strat_samps_sf = pblapply(pre_lids, function(b){
#   met_l = lapply(mets_to_test, function(m){
#     r = b[[m]]
#     #stratify raster using histogram bins
#     srast = strat_histogram(r,100)
#     
#     #get total number of nonNA pixels in the scene
#     n_tot = global(r, 'notNA', na.rm=T)$notNA
#     samp_size = ceiling(0.1*n_tot)
#     
#     #get sample
#     set.seed(123)
#     samp = sample_strat(srast,nSamp = samp_size)
#   })
#   names(met_l) = mets_to_test
#   return(met_l)
# })
# names(lid_strat_samps_sf) = names(pre_lids)
# 
# #extract pre and post-harvest sample data for each metric
# pre_lid_spectral_l = pblapply(names(pre_ps_byblock_l), function(d){
#   bl = lapply(names(lid_strat_samps_sf), function(b){
#     ml = lapply(mets_to_test, function(m){
#       
#       cat(paste0('\rd="', d, '"; b="', b, '"; m="', m, '"'))
#       
#       p = pre_ps_byblock_l[[d]]
#       ps = p[[b]]
#       met = pre_lids[[b]][[m]]
#       if(!ext(met)==ext(ps)){met = extend(met,ps)} #extend lidar to match planetscope if necessary
#       r = c(met,ps) #stack lidar and planetscope
#       e = extract(y=lid_strat_samps_sf[[b]][[m]], x=r)
#       return(e)
#       
#     })
#     names(ml) = mets_to_test
#     return(ml)
#   })
#   names(bl) = names(lid_strat_samps_sf)
#   return(bl)
# })
# names(pre_lid_spectral_l) = names(pre_ps_byblock_l)

#sample by mosaic (doesn't work since z score and SM values are calculated for cropped data)
{
  # #identify sample points for each metric based on pre-harvest lidar
  # lid_strat_samps_sf = lapply(names(lid_pre_mosaic), function(x){
  #   #get metric raster
  #   r = lid_pre_mosaic[[x]]
  #   #stratify raster using histogram bins
  #   srast = strat_histogram(r,100)
  #   
  #   #get total number of nonNA pixels in the scene
  #   n_tot = global(r, 'notNA', na.rm=T)$notNA
  #   samp_size = ceiling(0.1*n_tot)
  #   
  #   #get sample
  #   set.seed(123)
  #   samp = sample_strat(srast,nSamp = samp_size)
  # })
  # names(lid_strat_samps) = names(lid_pre_mosaic)
  # 
  # #extract pre and post-harvest sample data for each metric
  # pre_lid_spectral_l = pblapply(names(lid_strat_samps), function(x){
  #   ds_l=lapply(names(pre_ps_mosaics_l), function(p){
  #     
  #     met = lid_pre_mosaic[[x]] #get lidar metric
  #     ps = pre_ps_mosaics_l[[p]] #get planetscope mosaic for a given dataset
  #     
  #     r = c(met,ps) #stack lidar and planetscope
  #     e = extract(y=lid_strat_samps[[x]], x=r)
  #     return(e)
  #   })
  #   names(ds_l) = names(pre_ps_mosaics_l)
  #   return(ds_l)
  # })
  # names(pre_lid_spectral_l) = names(lid_strat_samps)
}

#----look at correlation matrices for each dataset----

corrplotter = function(x){
  v = as.data.frame(x, na.rm=T) |> na.omit()
  if(length(v[colnames(v)=='ID'])>0){v = v |> select(-ID)}
  c = cor(v)
  corrplot::corrplot(c, type='lower',order = 'hclust',addCoef.col = 'black',addCoefasPercent = T)
}
# 
# #differences between metrics?
# corrplotter(pre_lid_spectral_l$`Non-normalized`$`12N_T3`$z_above2)
# corrplotter(pre_lid_spectral_l$`Non-normalized`$`12N_T3`$z_p95)
# corrplotter(pre_lid_spectral_l$`Non-normalized`$`12N_T3`$z_sd)
# 
# #difference between blocks?
# corrplotter(pre_lid_spectral_l$`Non-normalized`$`12N_T3`$z_above2)
# corrplotter(pre_lid_spectral_l$`Non-normalized`$`12L_C5`$z_above2)
# 
# #difference between datasets?
# corrplotter(pre_lid_spectral_l$`Non-normalized`$`12N_T3`$z_above2)
# corrplotter(pre_lid_spectral_l$`SM`$`12N_T3`$z_above2)

#----define base features based on correlation analysis----
base_feats = c('coastal_blue', 'blue', 'green_i', 'green', 'yellow', 'red', 'rededge', 'nir'
               ,'GNDVI','Hue', 'GCC','SAVI','NDVI','YNDVI')

#----make dataframes of pre and post harvest datasets with only base features----

pre_lid_base_l = pblapply(names(pre_ps_byblock_l), function(d){
  bl = lapply(block_ids,function(b){
    ml = lapply(mets_to_test, function(met){
      
      cat(paste0('\rd="', d, '"; b="', b, '"; m="', met, '"'))
      
      ps = pre_ps_byblock_l[[d]][[b]]
      ps = ps[[base_feats]]
      lid = pre_lids[[b]][[met]]
      if(!ext(lid)==ext(ps)){lid = extend(lid,ps)}
      r = c(lid,ps)
      v = as.data.frame(r,na.rm=T)
      return(v)
    })
    names(ml) = mets_to_test
    return(ml)
  })
  names(bl) = block_ids
  return(bl)
})
names(pre_lid_base_l) = names(pre_ps_byblock_l)

post_lid_base_l = pblapply(names(post_ps_byblock_l), function(d){
  bl = lapply(block_ids,function(b){
    ml = lapply(mets_to_test, function(met){
      
      cat(paste0('\rd="', d, '"; b="', b, '"; m="', met, '"'))
      
      ps = post_ps_byblock_l[[d]][[b]]
      ps = ps[[base_feats]]
      lid = post_lids[[b]][[met]]
      if(!ext(lid)==ext(ps)){lid = extend(lid,ps)}
      r = c(lid,ps)
      v = as.data.frame(r,na.rm=T)
      return(v)
    })
    names(ml) = mets_to_test
    return(ml)
  })
  names(bl) = block_ids
  return(bl)
})
names(post_lid_base_l) = names(post_ps_byblock_l)



#----initial modelling (just spectral features)----

#regularized linear models
{
  #   
  # #test regularised linear regression models on stratified sample
  # {
  # library(glmnet)
  # mod_grid = expand.grid(
  #   dataset = unique(results_meta_df$dataset),
  #   block_id = unique(results_meta_df$block_id),
  #   lid_met = mets_to_test,
  #   stringsAsFactors = F,
  #   alpha = c(0,0.5,1)
  # )
  # 
  # simple_lm__samp_l = pblapply(1:nrow(mod_grid), function(i){
  #   d = mod_grid$dataset[i]
  #   b = mod_grid$block_id[i]
  #   met = mod_grid$lid_met[i]
  #   a = mod_grid$alpha[i]
  #   
  #   cat('\r',i,'/',nrow(mod_grid))
  #   
  #   ds = pre_lid_spectral_l[[d]][[b]][[met]] |> 
  #     st_drop_geometry() |>
  #     select(all_of(c(met,base_feats)))|>
  #     # select(-ID) |>
  #     na.omit() |> 
  #     as.matrix() |>
  #     scale()
  #   
  #   cvmod = glmnet::cv.glmnet(x = ds[,2:ncol(ds)], y = ds[,1], alpha=a, family='gaussian', nfolds=10)
  #   
  #   # --- Cross-validation performance ---
  #   y_true <- ds[, 1]
  #   y_pred <- predict(cvmod, newx = ds[, -1], s = cvmod$lambda.min)
  #   
  #   rss_cv  <- sum((y_true - y_pred)^2)                     # residual sum of squares
  #   tss_cv  <- sum((y_true - mean(y_true))^2)               # total sum of squares
  #   r2_cv   <- 1 - rss_cv / tss_cv                          # coefficient of determination
  #   rmse_cv <- sqrt(mean((y_true - y_pred)^2))              # root mean square error
  #   
  #   
  #   # --- Pre-harvest full dataset ---
  #   pre_full <- na.omit(pre_lid_base_l[[d]][[b]][[met]]) |> as.matrix()
  #   pre_y_true <- pre_full[, 1]
  #   pre_y_pred <- predict(cvmod, newx = pre_full[, -1], s = cvmod$lambda.min)
  #   
  #   rss_pre  <- sum((pre_y_true - pre_y_pred)^2)
  #   tss_pre  <- sum((pre_y_true - mean(pre_y_true))^2)
  #   r2_pre   <- 1 - rss_pre / tss_pre
  #   rmse_pre <- sqrt(mean((pre_y_true - pre_y_pred)^2))
  #   
  #   
  #   # --- Post-harvest full dataset ---
  #   post_full <- na.omit(post_lid_base_l[[d]][[b]][[met]]) |> as.matrix()
  #   post_y_true <- post_full[, 1]
  #   post_y_pred <- predict(cvmod, newx = post_full[, -1], s = cvmod$lambda.min)
  #   
  #   rss_post  <- sum((post_y_true - post_y_pred)^2)
  #   tss_post  <- sum((post_y_true - mean(post_y_true))^2)
  #   r2_post   <- 1 - rss_post / tss_post
  #   rmse_post <- sqrt(mean((post_y_true - post_y_pred)^2))
  #   
  #   
  #   # --- Combine results ---
  #   res_tbl <- tibble(
  #     rss_cv = rss_cv,
  #     tss_cv = tss_cv,
  #     r2_cv = r2_cv,
  #     rmse_cv = rmse_cv,
  #     
  #     rss_pre = rss_pre,
  #     tss_pre = tss_pre,
  #     r2_pre = r2_pre,
  #     rmse_pre = rmse_pre,
  #     
  #     rss_post = rss_post,
  #     tss_post = tss_post,
  #     r2_post = r2_post,
  #     rmse_post = rmse_post
  #   )
  #   
  #   res_tbl$model[1] = list(cvmod)
  #   
  #   return(res_tbl)
  # })
  # 
  # #combine results into dataframe
  # lm_tot__samp_df = cbind(mod_grid, bind_rows(simple_lm_l))
  # }
  # 
  # #regularized linear regression models on full data
  # {
  #   simple_lm_full_l = pblapply(1:nrow(mod_grid), function(i){
  #     d = mod_grid$dataset[i]
  #     b = mod_grid$block_id[i]
  #     met = mod_grid$lid_met[i]
  #     a = mod_grid$alpha[i]
  #     
  #     cat('\r',i,'/',nrow(mod_grid),'---',paste0('d ="', d, '"; b ="', b, '"; m ="', met, '"'))
  #     
  #     ds = pre_lid_base_l[[d]][[b]][[met]] |> 
  #       st_drop_geometry() |>
  #       select(all_of(c(met,base_feats)))|>
  #       # select(-ID) |>
  #       na.omit() |> 
  #       as.matrix() |>
  #       scale()
  #     
  #     cvmod = glmnet::cv.glmnet(x = ds[,2:ncol(ds)], y = ds[,1], alpha=a, family='gaussian', nfolds=10)
  #     
  #     # # --- Cross-validation performance ---
  #     # y_true <- ds[, 1]
  #     # y_pred <- predict(cvmod, newx = ds[, -1], s = cvmod$lambda.min)
  #     # 
  #     # rss_cv  <- sum((y_true - y_pred)^2)                     # residual sum of squares
  #     # tss_cv  <- sum((y_true - mean(y_true))^2)               # total sum of squares
  #     # r2_cv   <- 1 - rss_cv / tss_cv                          # coefficient of determination
  #     # rmse_cv <- sqrt(mean((y_true - y_pred)^2))              # root mean square error
  #     
  #     
  #     # --- Pre-harvest full dataset ---
  #     pre_full <- na.omit(pre_lid_base_l[[d]][[b]][[met]]) |> as.matrix()
  #     pre_y_true <- pre_full[, 1]
  #     pre_y_pred <- predict(cvmod, newx = pre_full[, -1], s = cvmod$lambda.min)
  #     
  #     # a = assess.glmnet(object=cvmod, newx = pre_full[, -1], newy=pre_y_true)
  #     
  #     rss_pre  <- sum((pre_y_true - pre_y_pred)^2)
  #     tss_pre  <- sum((pre_y_true - mean(pre_y_true))^2)
  #     r2_pre   <- 1 - rss_pre / tss_pre
  #     rmse_pre <- sqrt(mean((pre_y_true - pre_y_pred)^2))
  #     
  #     
  #     # --- Post-harvest full dataset ---
  #     post_full <- na.omit(post_lid_base_l[[d]][[b]][[met]]) |> as.matrix()
  #     post_y_true <- post_full[, 1]
  #     post_y_pred <- predict(cvmod, newx = post_full[, -1], s = cvmod$lambda.min)
  #     
  #     rss_post  <- sum((post_y_true - post_y_pred)^2)
  #     tss_post  <- sum((post_y_true - mean(post_y_true))^2)
  #     r2_post   <- 1 - rss_post / tss_post
  #     rmse_post <- sqrt(mean((post_y_true - post_y_pred)^2))
  #     
  #     
  #     # --- Combine results ---
  #     res_tbl <- tibble(
  #       # rss_cv = rss_cv,
  #       # tss_cv = tss_cv,
  #       # r2_cv = r2_cv,
  #       # rmse_cv = rmse_cv,
  #       
  #       rss_pre = rss_pre,
  #       tss_pre = tss_pre,
  #       r2_pre = r2_pre,
  #       rmse_pre = rmse_pre,
  #       
  #       rss_post = rss_post,
  #       tss_post = tss_post,
  #       r2_post = r2_post,
  #       rmse_post = rmse_post
  #     )
  #     
  #     res_tbl$model[1] = list(cvmod)
  #     
  #     return(res_tbl)
  #   })
  #   
  #   #combine results into dataframe
  #   lm_tot_full_df = cbind(mod_grid, bind_rows(simple_lm_full_l))
  # }
  #   
}

#random forest on full data
{
  library(ranger)
  
  #compare different values of num.trees (RESULT: 1000 seems good)
  {
    rf_grid = expand.grid(
      # dataset = unique(results_meta_df$dataset),
      # block_id = unique(results_meta_df$block_id),
      dataset = 'Z',
      block_id = '12L_D345',
      lid_met = mets_to_test,
      stringsAsFactors = F,
      ntree = seq(500,3000,500)
    )
    
    rf_ntreetune_l = pblapply(1:nrow(rf_grid), function(i){
      
      d = rf_grid$dataset[i]
      b = rf_grid$block_id[i]
      met = rf_grid$lid_met[i]
      nt = rf_grid$ntree[i]
      
      cat('\r',i,'/',nrow(rf_grid),'---',paste0('d ="', d, '"; b ="', b, '"; m ="', met, '"; nt = ',nt,''))
      
      ds = pre_lid_base_l[[d]][[b]][[met]] |> 
        st_drop_geometry() |>
        select(all_of(c(met,base_feats)))|>
        # select(-ID) |>
        na.omit() |> 
        as.matrix()
      
      mod = ranger(y=ds[,1], x = ds[,2:ncol(ds)], data=ds, importance='permutation', num.trees=nt)
      
      #get performance metrics for pre-harvest dataset
      pre_met_true = ds[,1]
      pre_met_pred = predictions(mod)
      pre_rmse = mean((pre_met_true - mean(pre_met_true))^2)
      pre_Rrmse = pre_rmse/mean(pre_met_true)
      pre_bias = mean(pre_met_true-pre_met_pred)
      
      pre_dat = list(pre_met_true, pre_met_pred)
      names(pre_dat) = c('true', 'predicted')
      
      #get performance metrics for post-harvest dataset
      post_full <- na.omit(post_lid_base_l[[d]][[b]][[met]]) |> as.matrix()
      post_met_true = post_full[,1]
      post_met_pred = predict(mod, post_full[,2:ncol(post_full)])$predictions
      
      post_dat = list(post_met_true, post_met_pred)
      names(post_dat) = c('true', 'predicted')
      
      post_rss= sum((post_met_true - post_met_pred)^2)
      post_tss = sum((post_met_true - mean(post_met_true))^2)
      post_r2  = 1 - (post_rss / post_tss)
      
      post_rmse = mean((post_met_true - mean(post_met_true))^2)
      post_Rrmse = post_rmse/mean(post_met_true)
      post_bias = mean(post_met_true-post_met_pred)
      
      # post_r2 = cor(post_met_true, predictions(post_met_pred))^2
      
      res_tbl = tibble(
        pre_oob_r2 = mod$r.squared,
        post_r2,
        pre_rmse,
        post_rmse,
        pre_Rrmse,
        post_Rrmse,
        pre_bias,
        post_bias
        
        
        # ,
        # model = list(mod),
        # post_true_Vs_pred = post_dat
      )
      
      results_ = list(res_tbl, mod, pre_dat, post_dat)
      names(results_) = c('Summary_stats', 'Model', 'Pre_harv_data', 'Post_harv_data')
      
      return(results_)
    })
    
    rf_tot_df = lapply(rf_ntreetune_l, function(x)x$Summary_stats) |> bind_rows()
    
    rf_res_full_df = cbind(rf_grid,rf_tot_df)
    ggplot(rf_res_full_df)+geom_line(aes(x =ntree, y = pre_oob_r2, color = lid_met))
  }
  
  #test case_weights parameter
  {
    rf_grid = expand.grid(
      # dataset = unique(results_meta_df$dataset),
      # block_id = unique(results_meta_df$block_id),
      dataset = 'Z',
      block_id = c('12L_C7','12L_C5','12L_B8C3','12L_C4'),
      lid_met = mets_to_test,
      stringsAsFactors = F,
      weight_type = c('none', 'dist_from_mean', 'dist_from_median', 'exp_dist_from_mean', 'exp_dist_from_median', 'log_dist_from_mean', 'log_dist_from_median')
    )
    
    rf_weighttune_l = pblapply(1:nrow(rf_grid), function(i){
      
      d = rf_grid$dataset[i]
      b = rf_grid$block_id[i]
      met = rf_grid$lid_met[i]
      weight_type = rf_grid$weight_type[i]
      
      cat('\r',i,'/',nrow(rf_grid),'---',paste0('d ="', d, '"; b ="', b, '"; m ="', met, '"; weight_type = "',weight_type,'"'))
      
      ds = pre_lid_base_l[[d]][[b]][[met]] |> 
        st_drop_geometry() |>
        dplyr::select(all_of(c(met,base_feats)))|>
        # select(-ID) |>
        na.omit() |> 
        as.matrix()
      
      cw = case_when(
        weight_type == 'dist_from_mean' ~ abs(mean(ds[,1])-ds[,1])+0.000001,
        weight_type == 'dist_from_median' ~ abs(median(ds[,1])-ds[,1])+0.000001,
        weight_type == 'exp_dist_from_mean' ~ exp(abs(mean(ds[,1])-ds[,1]))+0.000001,
        weight_type == 'exp_dist_from_median' ~ exp(abs(median(ds[,1])-ds[,1]))+0.000001,
        weight_type == 'log_dist_from_mean' ~ log(abs(mean(ds[,1])-ds[,1]),0.0001),
        weight_type == 'log_dist_from_median' ~ log(abs(median(ds[,1])-ds[,1]),0.0001),
        .default=NULL
      )
      
      mod = ranger(y=ds[,1], x = ds[,2:ncol(ds)]
                   # , importance='gini'
                   , num.trees=500, case.weights = cw)
      
      #get performance metrics for pre-harvest dataset
      pre_met_true = ds[,1]
      pre_met_pred = predictions(mod)
      pre_rmse = mean((pre_met_true - mean(pre_met_true))^2)
      pre_Rrmse = pre_rmse/mean(pre_met_true)
      pre_bias = mean(pre_met_true-pre_met_pred,na.rm=T)
      
      pre_dat = list(pre_met_true, pre_met_pred)
      names(pre_dat) = c('true', 'predicted')
      
      #get performance metrics for post-harvest dataset
      post_full <- na.omit(post_lid_base_l[[d]][[b]][[met]]) |> as.matrix()
      post_met_true = post_full[,1]
      post_met_pred = predict(mod, post_full[,2:ncol(post_full)])$predictions
      
      post_dat = list(post_met_true, post_met_pred)
      names(post_dat) = c('true', 'predicted')
      
      post_rss= sum((post_met_true - post_met_pred)^2)
      post_tss = sum((post_met_true - mean(post_met_true))^2)
      post_r2  = 1 - (post_rss / post_tss)
      
      post_rmse = mean((post_met_true - mean(post_met_true))^2)
      post_Rrmse = post_rmse/mean(post_met_true)
      post_bias = mean(post_met_true-post_met_pred)
      
      # post_r2 = cor(post_met_true, predictions(post_met_pred))^2
      
      res_tbl = tibble(
        pre_oob_r2 = mod$r.squared,
        post_r2,
        pre_rmse,
        post_rmse,
        pre_Rrmse,
        post_Rrmse,
        pre_bias,
        post_bias
      )
      
      results_ = list(res_tbl, mod, pre_dat, post_dat)
      names(results_) = c('Summary_stats', 'Model', 'Pre_harv_data', 'Post_harv_data')
      
      return(results_)
    })
    
    rf_tot_df = lapply(rf_weighttune_l, function(x)x$Summary_stats) |> bind_rows()
    
    rf_res_full_df = cbind(rf_grid,rf_tot_df)
    ggplot(rf_res_full_df)+geom_col(aes(x =weight_type, y = pre_oob_r2))+facet_grid(cols=vars(lid_met), rows=vars(block_id))+theme(axis.text.x=element_text(angle=90))
    # ggplot(rf_res_full_df)+geom_col(aes(x =weight_type, y = post_Rrmse))+facet_wrap(vars(lid_met))+theme(axis.text.x=element_text(angle=90))
    
  }
}

#test MARS with and without caseweights
{
  library(earth)
  
  mars_grid = expand.grid(
    # dataset = unique(results_meta_df$dataset),
    # block_id = unique(results_meta_df$block_id),
    dataset = 'Non-normalized',
    block_id = '12L_D345',
    lid_met = mets_to_test,
    stringsAsFactors = F,
    weighted = c(1,0)
  )
  
  mars_weighttune_l = pblapply(1:nrow(rf_grid), function(i){
    
    d = rf_grid$dataset[i]
    b = rf_grid$block_id[i]
    met = rf_grid$lid_met[i]
    weight_type = rf_grid$weight_type[i]
    
    cat('\r',i,'/',nrow(rf_grid),'---',paste0('d ="', d, '"; b ="', b, '"; m ="', met, '"; weight_type = "',weight_type,'"'))
    
    ds = pre_lid_base_l[[d]][[b]][[met]] |> 
      st_drop_geometry() |>
      dplyr::select(all_of(c(met,base_feats)))|>
      # select(-ID) |>
      na.omit() |> 
      as.matrix()
    
    cw = case_when(
      weight_type == 'dist_from_mean' ~ abs(mean(ds[,1])-ds[,1])+0.0001,
      weight_type == 'dist_from_median' ~ abs(median(ds[,1])-ds[,1])+0.0001,
      weight_type == 'exp_dist_from_mean' ~ exp(abs(mean(ds[,1])-ds[,1]))+0.0001,
      weight_type == 'exp_dist_from_median' ~ exp(abs(median(ds[,1])-ds[,1]))+0.0001,
      weight_type == 'log_dist_from_mean' ~ log1p(abs(mean(ds[,1])-ds[,1])),
      weight_type == 'log_dist_from_median' ~ log1p(abs(median(ds[,1])-ds[,1])),
      .default=NULL
    )
    
    mod = ranger(y=ds[,1], x = ds[,2:ncol(ds)]
                 # , importance='gini'
                 , num.trees=500, case.weights = cw)
    
    #get performance metrics for pre-harvest dataset
    pre_met_true = ds[,1]
    pre_met_pred = predictions(mod)
    pre_rmse = mean((pre_met_true - mean(pre_met_true))^2)
    pre_Rrmse = pre_rmse/mean(pre_met_true)
    pre_bias = mean(pre_met_true-pre_met_pred,na.rm=T)
    
    pre_dat = list(pre_met_true, pre_met_pred)
    names(pre_dat) = c('true', 'predicted')
    
    #get performance metrics for post-harvest dataset
    post_full <- na.omit(post_lid_base_l[[d]][[b]][[met]]) |> as.matrix()
    post_met_true = post_full[,1]
    post_met_pred = predict(mod, post_full[,2:ncol(post_full)])$predictions
    
    post_dat = list(post_met_true, post_met_pred)
    names(post_dat) = c('true', 'predicted')
    
    post_rss= sum((post_met_true - post_met_pred)^2)
    post_tss = sum((post_met_true - mean(post_met_true))^2)
    post_r2  = 1 - (post_rss / post_tss)
    
    post_rmse = mean((post_met_true - mean(post_met_true))^2)
    post_Rrmse = post_rmse/mean(post_met_true)
    post_bias = mean(post_met_true-post_met_pred)
    
    # post_r2 = cor(post_met_true, predictions(post_met_pred))^2
    
    res_tbl = tibble(
      pre_oob_r2 = mod$r.squared,
      post_r2,
      pre_rmse,
      post_rmse,
      pre_Rrmse,
      post_Rrmse,
      pre_bias,
      post_bias
    )
    
    results_ = list(res_tbl, mod, pre_dat, post_dat)
    names(results_) = c('Summary_stats', 'Model', 'Pre_harv_data', 'Post_harv_data')
    
    return(results_)
  })
  
  rf_tot_df = lapply(rf_weighttune_l, function(x)x$Summary_stats) |> bind_rows()
  
  rf_res_full_df = cbind(rf_grid,rf_tot_df) 
}

#----generate GLCM rasters----

#define function to process GLCM rasters without edge effects

library(GLCMTextures)

#define glcm_quantization
n_levels = 32
quant_method = 'range'
glcm_dir = paste0(ps_dir,'/GLCMTextures_','QuantLevels=',n_levels,'_QuantMethod=',quant_method)
dir.check(glcm_dir)

glcmerizer_max = function(fid
                          , glcm_dir
                          , windows = c(3,7,11,15)
                          , feats = unique(results_meta_df$var_name)
                          ,textures = c(
                            "glcm_contrast", "glcm_dissimilarity", "glcm_homogeneity"
                            , "glcm_ASM", "glcm_entropy", "glcm_mean"
                            , "glcm_variance", "glcm_correlation"
                          )){
  #function for making glcm metrics using the unprocessed planetscope data in order to avoid edge effects with larger window sizes.
  #this function will make combinations of all input window sizes and features
  #... is additional arguments to glcm_textures()
  #glcm_dir is the parent directory
  
  #get scene id
  fid
  id = paste(str_split_1(basename(fid),'_')[1:4],collapse='_')
  
  #reprocess and save full scene file
  raw_fixed = paste0(ps_dir,'/raw_withindices_fixedgeom')
  dir.check(raw_fixed)
  
  raw_f = paste0(raw_fixed,'/',id,'.tif')
  if(!file.exists(raw_f)){
    #get raster of spectral features
    {
      #load raster
      r = rast(fid)
      #rescale raster
      r = r/scale_factor
      
      #calculate vegetation indices
      si = spectralIndices(img = r, blue = 'blue', green = 'green', red = 'red', nir = 'nir', redEdge1 = 'rededge'
                           # , skipRefCheck = T
      )
      #remove EVI2 (since spectralIndices doesn't calculate it well for some reason)
      si = subset(si, 'EVI2', negate=T)
      
      #calculate additional indices with custom functions
      vi = rast.batch.functions(r, include_input = F,
                                #band mapping
                                coastal_blue=1, blue=2, green_i=3, green=4, yellow=5, red=6, rededge=7, nir=8, 
                                #list of indices to calculate
                                fl = c(
                                  hue.vi,
                                  gcc.vi,
                                  ndgr.vi,
                                  bi.vi,
                                  ci.vi,
                                  cri550.vi,
                                  gli.vi,
                                  # tvi.vi, #no TVI since spectralIndices handles this
                                  varig.vi,
                                  evi2.vi,
                                  yndvi.vi
                                ))
      
      r = c(r, si, vi)
    }
    #reproject and resample raw spectral raster
    {
      rpext = terra::project(ext(r),from=crs(r),to=target_crs)
      rp_ext_snapped = snap_extent(rpext, target_res)
      template_rast = rast(rp_ext_snapped,res=target_res,crs=target_crs)
      r = project(r,template_rast)
    }
    writeRaster(r,raw_f)
  }
  
  #generate and save glcm layers
  {
    
    blocks_p_wrapped = wrap(blocks_p)
    
    grid = expand.grid(windows = windows, feats = feats, textures = textures, stringsAsFactors = F)
    
    future_lapply(1:nrow(grid), function(i){
      
      cat('\rProcessing ',i,'/',nrow(grid))
      
      feat = grid$feats[i]
      window = grid$windows[i]
      texture = grid$textures[i]
      
      
      #set parent directory for raw glcm data
      raw_g_dir = paste0(glcm_dir,'/Non-normalized_non-cropped')
      dir.check(raw_g_dir)
      
      #set parent directory based on scene id
      raw_id_dir = paste0(raw_g_dir,'/',id)
      dir.check(raw_id_dir)
      
      #check if raw texture file exists, make texture if no
      g_string = paste0(texture,'_',feat,'_',window)
      raw_g_f = paste0(raw_id_dir,'/',g_string,'.tif')
      if(!file.exists(raw_g_f)){
        rf = rast(raw_f,lyrs=feat)
        g = glcm_textures(rf,w=window,n_levels = n_levels,quant_method = quant_method, metrics = textures)
        
        #fix names
        ng = paste0(names(g),'_',feat,'_',window)
        names(g) = ng
        
        #save individual layers (will make it less annoying to load them later)
        lapply(ng, function(x){
          gf = paste0(raw_id_dir,'/',x,'.tif')
          if(!file.exists(gf)){writeRaster(g[[x]],gf)}
        })
      }
      
      #save glcm files transformed for datasets and cropped to thinning blocks
      nn_g_dir = file.path(glcm_dir,'Non-normalized')
      dir.check(nn_g_dir)
      sm_g_dir = file.path(glcm_dir, 'SM')
      dir.check(sm_g_dir)
      z_g_dir = file.path(glcm_dir,'Z')
      dir.check(z_g_dir)
      zr_g_dir = file.path(glcm_dir,'Zrobust')
      dir.check(zr_g_dir)
      
      lapply(block_ids, function(b){
        
        #save non-normalized data
        nn_blockdir = file.path(nn_g_dir,b)
        dir.check(nn_blockdir)
        
        nn_id_dir = file.path(nn_blockdir,id)
        dir.check(nn_id_dir)
        
        nn_g_f = paste0(nn_id_dir,'/',g_string,'.tif')
        if(!file.exists(nn_g_f)){
          r = rast(raw_g_f)
          .blocks_p = unwrap(blocks_p_wrapped)
          block = .blocks_p[.blocks_p$BLOCKNUM==b,]
          rc = crop(r, block, mask=T)
          writeRaster(rc,nn_g_f)
        }
        
        #save softmax data
        sm_blockdir = file.path(sm_g_dir,b)
        dir.check(sm_blockdir)
        
        sm_id_dir = file.path(sm_blockdir,id)
        dir.check(sm_id_dir)
        
        sm_g_f = paste0(sm_id_dir,'/',g_string,'.tif')
        if(!file.exists(sm_g_f)){
          r = rast(nn_g_f)
          rs = softmax(r, append_name = F)
          writeRaster(rs, sm_g_f)
        }
        
        #save Z-normalized data
        z_blockdir = file.path(z_g_dir,b)
        dir.check(z_blockdir)
        
        z_id_dir = file.path(z_blockdir,id)
        dir.check(z_id_dir)
        
        z_g_f = paste0(z_id_dir,'/',g_string,'.tif')
        if(!file.exists(z_g_f)){
          r = rast(nn_g_f)
          zr = z_rast(r)
          writeRaster(zr,z_g_f)
        }
        
        #save robust z-normalized data
        zr_blockdir = file.path(zr_g_dir,b)
        dir.check(zr_blockdir)
        
        zr_id_dir = file.path(zr_blockdir,id)
        dir.check(zr_id_dir)
        
        zr_g_f = paste0(zr_id_dir,'/',g_string,'.tif')
        if(!file.exists(zr_g_f)){
          r = rast(nn_g_f)
          zr = z_rast(r,robust = T)
          writeRaster(zr,zr_g_f)
        }
        
      })
      
      
      
    }
    # ,cl = 'future'
    ,future.seed=T
    )
    
  }
  
}



pre_raw_f = raw_rasters[str_detect(raw_rasters, pre_harv_id) & !str_detect(raw_rasters,'udm')][1]

pre_raw_f2 = raw_rasters[str_detect(raw_rasters,"20210729_191439_38_2274") & !str_detect(raw_rasters,'udm')][1]

post_raw_f = raw_rasters[str_detect(raw_rasters, post_harv_id) & !str_detect(raw_rasters,'udm')][1]

# clust =makeCluster(4)
# plan('cluster',workers=clust)

glcmerizer_max(pre_raw_f, glcm_dir, windows= c(3,5,7,11,15))
# gc()
glcmerizer_max(post_raw_f, glcm_dir, windows= c(3,5,7,11,15))

# stopCluster(clust)
# plan('sequential')

#----load GLCM rasters for desired spectral features----

#get list of all possible glcm features
glcm.feat.grid = expand.grid(texture = c("glcm_contrast", "glcm_dissimilarity", "glcm_homogeneity", "glcm_ASM",
                                         "glcm_entropy", "glcm_mean", "glcm_variance"
                                         # , "glcm_correlation"
)
# , block_id = unique(results_meta_df$block_id)
# , dataset = unique(results_meta_df$dataset)
, window = c(3
             ,5
             ,7
             ,11
             # ,15
)
, spec_feat = base_feats[9:length(base_feats)]
) |> mutate(glcm_lyr = paste0(texture,'_',spec_feat,'_',window))

#define function for loading GLCM rasters
glcmloader = function(id, 
                      glcm_lyrs, 
                      datasets = unique(results_meta_df$dataset), 
                      block_ids = unique(results_meta_df$block_id), 
                      dir = glcm_dir,
                      as.mat=F){
  #function to load glcm layers for datasets and block_ids. as.mat=T will output a matrix instead of a spatraster 
  
  #get full list of glcm_files in the target directory
  fl = list.files(dir, recursive = T, full.names = T)
  
  #restrict fl to scene id
  fl = fl[str_detect(fl,id)]
  
  dl = lapply(datasets, function(d){
    bl = lapply(block_ids, function(b){
      # fl_ = fl[str_detect(fl, d) 
      #          & str_detect(fl,b)
      #          & basename(file_path_sans_ext(fl)) %in% glcm_lyrs]
      fl_ = fl[
        basename(dirname(dirname(dirname(fl)))) == d
        & basename(dirname(dirname(fl))) == b
        & basename(file_path_sans_ext(fl)) %in% glcm_lyrs]
      r = pblapply(fl_, rast) |> rast()
      if(as.mat){r = values(r,na.rm=T)}
      return(r)
    })
    names(bl) = block_ids
    return(bl)
  })
  names(dl) = datasets
  return(dl)
}

# #load full set of glcm rasters for each block and dataset
# pre_glcm_r_l = glcmloader(id = pre_harv_id, glcm_lyrs = glcm.feat.grid$glcm_lyr)
# post_glcm_r_l = glcmloader(id = post_harv_id, glcm_lyrs = glcm.feat.grid$glcm_lyr)


#dig into GLCM feature importance
{
  # #compare post-harvest textures between mosaics for different datasets (Results: Mean, Contrast, SAVI, and w(7,11) dominate)
  # {
  #   
  #   post_vs_res_l = pblapply(names(post_ps_mosaics_l), function(d){
  #     
  #     #define output directory for VSURF models
  #     glcm_test_dir = paste0(ps_dir,'/GLCM_feat_importance')
  #     dir.check(glcm_test_dir)
  #     
  #     id_dir = paste0(glcm_test_dir, '/',post_harv_id)
  #     dir.check(id_dir)
  #     
  #     npix = 2000
  #     
  #     samp_dir = paste0(id_dir,'/Mosaic_pixperblock_w(3,5,7,11)_=',npix)
  #     dir.check(samp_dir)
  #     
  #     #get glcm
  #     glcm_r = post_glcm_r_l[[d]] |> 
  #       sprc() |> 
  #       merge()
  #     
  #     #get lidar mosaic
  #     lid = lid_post_mosaic
  #     if(ext(lid) != ext(glcm_r)){lid = extend(lid,glcm_r)}
  #     
  #     #stack_all
  #     combined_scene = c(lid, glcm_r)
  #     
  #     #get random sample of pixels based on size of smallest block 
  #     set.seed(123)
  #     combined_values_sample = lapply(block_ids, function(b){
  #       block = blocks_p[blocks_p$BLOCKNUM==b,]
  #       r_b = crop(combined_scene, block, mask=T)
  #       s = spatSample(x=r_b, size=npix, method='random', na.rm=T, exhaustive=T)
  #       return(s)
  #     }) |> bind_rows()
  #     
  #     ml = lapply(mets_to_test, function(m){
  #       
  #       cat(paste0('\nd = "',d,'"--- m = "',m,'"'))
  #       cat('\n------Start:',as.character(Sys.time()),'----')
  #       f = paste0(samp_dir,'/',m,'.rds')
  #       if(!file.exists(f)){
  #         
  #         dat = combined_values_sample |> select(all_of(c(m,names(glcm_r))))
  #         set.seed(123)
  #         vsurf = VSURF(y=dat[,1], x = dat[,2:ncol(dat)], RFimplem='ranger', parallel=T)
  #         
  #         vsurf_importance_df <- data.frame(
  #           variable = colnames(dat[, -1])[vsurf$imp.mean.dec.ind],  # exclude response column
  #           importance = vsurf$imp.mean.dec,
  #           lid_met = m,
  #           dataset = d
  #         )
  #         
  #         vs_res = list(vsurf_importance_df, vsurf, dat)
  #         names(vs_res) = c('Results_df', 'VSURF_object', 'Input_data')
  #         
  #         write_rds(vs_res,f)
  #       } else {
  #         vs_res = read_rds(f)
  #       }
  #       cat('\n------End:',as.character(Sys.time()),'----')
  #       return(vs_res)
  #     })
  #     names(ml) = mets_to_test
  #     return(ml)
  #   })
  #   names(post_vs_res_l) = names(post_ps_mosaics_l)
  #   
  #   varimp_df = lapply(vs_res_l, function(d){
  #     ml = lapply(d, function(m)m$Results_df) |> bind_rows()
  #   }) |> bind_rows()
  #   
  #   post_pred_var_l = lapply(names(post_vs_res_l), function(d){
  #     obj = post_vs_res_l[[d]]
  #     ml = lapply(names(obj), function(m){
  #       a = obj[[m]]
  #       v = a$VSURF_object
  #       v_idx = v$varselect.pred
  #       if(!is.null(v_idx)){
  #         return(names(a$Input_data)[v_idx])
  #       } else {
  #           return(NULL)
  #         }
  #     })
  #     names(ml) = names(obj)
  #     return(ml)
  #   })
  #   names(post_pred_var_l) = names(post_vs_res_l)
  #   
  # }
  # 
  # #compare pre-harvest textures between mosaics for different datasets
  # {
  #   
  #   pre_vs_res_l = pblapply(names(pre_ps_mosaics_l), function(d){
  #     
  #     #define output directory for VSURF models
  #     glcm_test_dir = paste0(ps_dir,'/GLCM_feat_importance')
  #     dir.check(glcm_test_dir)
  #     
  #     id_dir = paste0(glcm_test_dir, '/',pre_harv_id)
  #     dir.check(id_dir)
  #     
  #     npix = 2000
  #     
  #     samp_dir = paste0(id_dir,'/Mosaic_pixperblock_w(3,5,7,11)_=',npix)
  #     dir.check(samp_dir)
  #     
  #     #get glcm
  #     glcm_r = pre_glcm_r_l[[d]] |> 
  #       sprc() |> 
  #       merge()
  #     
  #     #get lidar mosaic
  #     lid = lid_pre_mosaic
  #     if(ext(lid) != ext(glcm_r)){lid = extend(lid,glcm_r)}
  #     
  #     #stack_all
  #     combined_scene = c(lid, glcm_r)
  #     
  #     #get random sample of pixels based on size of smallest block 
  #     set.seed(123)
  #     combined_values_sample = lapply(block_ids, function(b){
  #       block = blocks_p[blocks_p$BLOCKNUM==b,]
  #       r_b = crop(combined_scene, block, mask=T)
  #       s = spatSample(x=r_b, size=npix, method='random', na.rm=T, exhaustive=T)
  #       return(s)
  #     }) |> bind_rows()
  #     
  #     ml = lapply(mets_to_test, function(m){
  #       
  #       cat(paste0('\nd = "',d,'"--- m = "',m,'"'))
  #       cat('\n------Start:',as.character(Sys.time()),'----')
  #       f = paste0(samp_dir,'/',m,'.rds')
  #       if(!file.exists(f)){
  #         
  #         dat = combined_values_sample |> select(all_of(c(m,names(glcm_r))))
  #         set.seed(123)
  #         vsurf = VSURF(y=dat[,1], x = dat[,2:ncol(dat)], RFimplem='ranger', parallel=T)
  #         
  #         vsurf_importance_df <- data.frame(
  #           variable = colnames(dat[, -1])[vsurf$imp.mean.dec.ind],  # exclude response column
  #           importance = vsurf$imp.mean.dec,
  #           lid_met = m,
  #           dataset = d
  #         )
  #         
  #         vs_res = list(vsurf_importance_df, vsurf, dat)
  #         names(vs_res) = c('Results_df', 'VSURF_object', 'Input_data')
  #         
  #         write_rds(vs_res,f)
  #       } else {
  #         vs_res = read_rds(f)
  #       }
  #       cat('\n------End:',as.character(Sys.time()),'----')
  #       return(vs_res)
  #     })
  #     names(ml) = mets_to_test
  #     return(ml)
  #   })
  #   names(pre_vs_res_l) = names(pre_ps_mosaics_l)
  #   
  #   varimp_df = lapply(pre_vs_res_l, function(d){
  #     ml = lapply(d, function(m)m$Results_df) |> bind_rows()
  #   }) |> bind_rows()
  #   
  #   pre_pred_var_l = lapply(names(pre_vs_res_l), function(d){
  #     obj = pre_vs_res_l[[d]]
  #     ml = lapply(names(obj), function(m){
  #       a = obj[[m]]
  #       v = a$VSURF_object
  #       v_idx = v$varselect.pred
  #       if(!is.null(v_idx)){
  #         return(names(a$Input_data)[v_idx])
  #       } else {
  #         return(NULL)
  #       }
  #     })
  #     names(ml) = names(obj)
  #     return(ml)
  #   })
  #   names(pre_pred_var_l) = names(pre_vs_res_l)
  #   print(pre_pred_var_l)
  #   
  #   
  # }
  
}

feats_to_use = base_feats

#load subset of glcm features to use for modelling
glcm.feat.grid.sub = expand.grid(
  texture = c(
    "glcm_contrast"
    # , "glcm_dissimilarity"
    # , "glcm_homogeneity"
    # , "glcm_ASM"
    # ,"glcm_entropy"
    , "glcm_mean"
    # , "glcm_variance"
    # , "glcm_correlation"
  )
  # , block_id = unique(results_meta_df$block_id)
  # , dataset = unique(results_meta_df$dataset)
  , window = c(
    3
    # ,5
    ,7
    ,11
    # ,15
  )
  , spec_feat = feats_to_use[9:length(base_feats)]
) |> mutate(glcm_lyr = paste0(texture,'_',spec_feat,'_',window)) 

pre_glcm_r_l_sub = glcmloader(pre_harv_id, glcm.feat.grid.sub$glcm_lyr)
post_glcm_r_l_sub = glcmloader(post_harv_id, glcm.feat.grid.sub$glcm_lyr)


#----SO MUCH MODELLING----

#define function that takes a model, uses it to make predictions using another dataset, and 

#make output directories
temp_transfer_mod_dir = paste0(ps_dir,'/temporal_transfer_modelling')
dir.check(temp_transfer_mod_dir)

featureset_dir = paste0(temp_transfer_mod_dir,'/allbands_8inds_glcmMeanConW(3,7,11)')
dir.check(featureset_dir)

dates_dir = paste0(featureset_dir, '/PreHarv=',pre_harv_id,'_PostHarv=',post_harv_id)
dir.check(dates_dir)

#get combinations of dataset, lidar metric, and block_id to test
grid = expand.grid(
  dataset = unique(results_meta_df$dataset),
  lid_met = mets_to_test,
  block_id = unique(results_meta_df$block_id),
  stringsAsFactors = F
)

#define function for calculating model performance given predicted values and observed values
perf_mets_calculator = function(pred, obs){
  
  residuals = obs - pred
  
  rss = sum(residuals^2)
  tss = sum((obs-mean(obs))^2)
  
  R2 = 1-(rss/tss)
  RMSE = sqrt(mean((obs-pred)^2))
  RRMSE = RMSE/mean(obs)
  Bias = mean(residuals)
  
  results = tibble(
    Metric = c('R2', 'RMSE', 'RRMSE', 'Bias'),
    Value = c(R2, RMSE, RRMSE, Bias)
  )
  return(results)
}

library(caret)
library(VSURF)
library(ranger) #random forest
library(earth) #MARS
library(yaImpute) #kNN

pblapply(1:nrow(grid), function(i){
  
  #define variables in loop iteration
  {
    d = grid$dataset[i]
    b = grid$block_id[i]
    met = grid$lid_met[i]
    
    cat('\n',i,'/',nrow(grid),'---',paste0('d ="', d, '"; b ="', b, '"; m ="', met, '" --- ',as.character(Sys.time()),'\n'))
  }
  #data loading, variable selection
  
  #load pre-harvest data
  {
    # tic()
    pre_lid = pre_lids[[b]][[met]]
    pre_ps = pre_ps_byblock_l[[d]][[b]] |> subset(feats_to_use)
    pre_glcm = pre_glcm_r_l_sub[[d]][[b]]
    
    #make lidar match satellite if necessary (since sometimes the blocks are cut off)
    if(ext(pre_lid)!=ext(pre_ps)){pre_lid = extend(pre_lid,pre_ps)}
    
    pre_r = c(pre_lid,pre_ps,pre_glcm)
    pre_dat = values(pre_r,na.rm=T)
    # toc()
  }
  
  #VSURF variable selection on pre_harvest dataset
  {
    cat('\nVSURF variable selection\n')
    
    vs_dir = paste0(dates_dir,'/VSURF')
    dir.check(vs_dir)
    
    vs_f = paste0(vs_dir,'/',d,'_',b,'_',met,'.rds')
    if(!file.exists(vs_f)){
      set.seed(123)
      pre_vsurf = VSURF(y=pre_dat[,1], x = pre_dat[,2:ncol(pre_dat)], RFimplem='ranger', parallel=T)
      write_rds(pre_vsurf,vs_f)
    } else {
      pre_vsurf = read_rds(vs_f)
    }
  }
  
  #extract variables from vsurf
  pre_vsurf.selected_indices <- pre_vsurf$varselect.pred
  
  selected_feats <- colnames(pre_dat)[pre_vsurf.selected_indices]
  
  #subset pre_harvest data to only include selected variables
  final_pre_dat = pre_dat[,c(met, selected_feats)]
  
  #load post-harvest data
  {
    
    # tic()
    post_lid = post_lids[[b]][[met]]
    post_ps = post_ps_byblock_l[[d]][[b]] |> subset(feats_to_use)
    post_glcm = post_glcm_r_l_sub[[d]][[b]]
    
    #make lidar match satellite if necessary (since sometimes the blocks are cut off)
    if(ext(post_lid)!=ext(post_ps)){post_lid = extend(post_lid,post_ps)}
    
    post_r = c(post_lid,post_ps,post_glcm)[[c(met, selected_feats)]] #subset to only have selected variables
    final_post_dat = values(post_r,na.rm=T)
    # toc()
    
  }
  
  #random forest modelling
  {
    
    cat('\nRandom Forest modelling\n')
    
    rf_dir = paste0(dates_dir,'/RandomForest')
    dir.check(rf_dir)
    
    rf_f = paste0(rf_dir,'/',d,'_',b,'_',met,'.rds')
    if(!file.exists(rf_f)){
      
      #model on pre-harvest data
      set.seed(123)
      mod = ranger(y = final_pre_dat[,1], x = final_pre_dat[,2:ncol(final_pre_dat)], importance='permutation', num.trees=1000)
      
      #predict on post-harvest data
      p = predict(mod, final_post_dat[,2:ncol(final_post_dat)])
      
      #get pre-harvest model performance metrics
      # perf_pre = metrica::metrics_summary(obs = as.vector(final_pre_dat[,1])
      #                                     , pred = mod$predictions
      #                                     , metrics_list=c('R2', 'RMSE', 'RRMSE', 'MBE')
      #                                     , type='regression')
      perf_pre = perf_mets_calculator(pred = mod$predictions, obs = as.vector(final_pre_dat[,1])) |> 
        mutate(dataset = d,
               block_id = b,
               lid_met = met,
               temporality = 'pre_harvest')
      #get post-harvest model performance metrics
      perf_post = perf_mets_calculator(obs = final_post_dat[,1], pred = p$predictions)|> 
        mutate(dataset = d,
               block_id = b,
               lid_met = met,
               temporality = 'post_harvest')
      
      #combine pre and post-harvest metrics
      perf_df = rbind(perf_pre,perf_post)
      
      rf_res_l = list(perf_df, mod, final_pre_dat, final_post_dat)
      names(rf_res_l) = c('Model_performance_summary', 'Model', 'Pre-harvest_data', 'Post_harvest_data')
      write_rds(rf_res_l,rf_f)
    }
    
  }
  
  #MARS modelling
  # {
  #   
  #   cat('\nMARS modelling\n')
  #   
  #   mars_dir = paste0(dates_dir,'/MARS')
  #   dir.check(mars_dir)
  #   
  #   mars_f = paste0(mars_dir,'/',d,'_',b,'_',met,'.rds')
  #   if(!file.exists(mars_f)){
  #     
  #     #tune model hyperparameters using the post-harvest data for cross validation
  #     tune_grid = expand.grid(
  #       degree = c(1,2,3),
  #       case.weights = c(F
  #                        # ,T
  #                        )
  #     )
  #     
  #     #test hyperparameters in tune_grid
  #     tune_l = pblapply(1:nrow(tune_grid), function(j){
  #       cat(j)
  #       
  #       #construct case weights
  #       if(tune_grid$case.weights[j]==F){
  #         cw = NULL
  #       } else {
  #         cw = abs(final_pre_dat[,1] - mean(final_pre_dat[,1]))
  #         zro_idx = which(cw==0)
  #         cw[zro_idx] = cw[zro_idx]+1e-8 #make it so that no weights are equal to zero
  #       }
  #       
  #       #model on pre-harvest data
  #       set.seed(123)
  #       mod = earth(y = final_pre_dat[,1], 
  #                   x = final_pre_dat[,2:ncol(final_pre_dat)],
  #                   degree = tune_grid$degree[j],
  #                   weights = cw)
  #       
  #       #predict on post-harvest data
  #       p = predict(mod, final_post_dat[,2:ncol(final_post_dat)])
  #       
  #       #get pre-harvest model performance metrics
  #       # perf_pre = metrica::metrics_summary(obs = as.vector(final_pre_dat[,1])
  #       #                                     , pred = mod$predictions
  #       #                                     , metrics_list=c('R2', 'RMSE', 'RRMSE', 'MBE')
  #       #                                     , type='regression')
  #       perf_pre = perf_mets_calculator(pred = mod$fitted.values, obs = as.vector(final_pre_dat[,1])) |> 
  #         mutate(dataset = d,
  #                block_id = b,
  #                lid_met = met,
  #                temporality = 'pre_harvest')
  #       #get post-harvest model performance metrics
  #       perf_post = perf_mets_calculator(obs = final_post_dat[,1], pred = p)|> 
  #         mutate(dataset = d,
  #                block_id = b,
  #                lid_met = met,
  #                temporality = 'post_harvest')
  #       
  #       #combine pre and post-harvest metrics
  #       perf_df = rbind(perf_pre,perf_post) |>
  #         mutate(
  #           tune.degree = tune_grid$degree[j],
  #           tune.case.weights = tune_grid$case.weights[j]
  #         )
  #       
  #       tune_res_l = list(perf_df, mod)
  #       names(tune_res_l) = c('Tuning_performance_dataframe', 'Model')
  #       return(tune_res_l)
  #     })
  #     
  #     return(perf_df)
  #     
  #     res_l = list(perf_df, mod, final_pre_dat, final_post_dat)
  #     names(rf_res_l) = c('Model_performance_summary', 'Model', 'Pre-harvest_data', 'Post_harvest_data')
  #     write_rds(res_l,mars_f)
  #   }
  # }
  
  #LASSO ridge regression
  {
    
    
  }
  
  
  
  
  #Spline modelling
  # {
  #   
  #   cat('\nspline modelling\n')
  #   
  #   spline_dir = paste0(dates_dir,'/MARS')
  #   dir.check(mars_dir)
  #   
  #   rf_f = paste0(mars_dir,'/',d,'_',b,'_',met,'.rds')
  #   if(!file.exists(mars_f)){
  #     
  #     #model on pre-harvest data
  #     set.seed(123)
  #     mod = earth(y = final_pre_dat[,1], x = final_pre_dat[,2:ncol(final_pre_dat)])
  #     
  #     #predict on post-harvest data
  #     p = predict(mod, final_post_dat[,2:ncol(final_post_dat)])
  #     
  #     #get pre-harvest model performance metrics
  #     # perf_pre = metrica::metrics_summary(obs = as.vector(final_pre_dat[,1])
  #     #                                     , pred = mod$predictions
  #     #                                     , metrics_list=c('R2', 'RMSE', 'RRMSE', 'MBE')
  #     #                                     , type='regression')
  #     perf_pre = perf_mets_calculator(pred = mod$fitted.values, obs = as.vector(final_pre_dat[,1])) |> 
  #       mutate(dataset = d,
  #              block_id = b,
  #              lid_met = met,
  #              temporality = 'pre_harvest')
  #     #get post-harvest model performance metrics
  #     perf_post = perf_mets_calculator(obs = final_post_dat[,1], pred = p)|> 
  #       mutate(dataset = d,
  #              block_id = b,
  #              lid_met = met,
  #              temporality = 'post_harvest')
  #     
  #     #combine pre and post-harvest metrics
  #     perf_df = rbind(perf_pre,perf_post)
  #     
  #     rf_res_l = list(perf_df, mod, final_pre_dat, final_post_dat)
  #     names(rf_res_l) = c('Model_performance_summary', 'Model', 'Pre-harvest_data', 'Post_harvest_data')
  #     write_rds(rf_res_l,rf_f)
  #   }
  # }
  
})

