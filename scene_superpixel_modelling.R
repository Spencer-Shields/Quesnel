# source('derive lidar metrics_202501020.R')
source('scene_global_stats_filtering_20251106.R')


#define pre-harvest and post-harvest ids
pre_harv_id = "20210713_191755_66_2414"
post_harv_id = "20241011_192907_10_251a"

desired_feats = c('coastal_blue', 
                  'blue', 
                  'green_i', 
                  'green',
                  'yellow',
                  'red',
                  'rededge',
                  'nir',
                  'NDVI')
library(supercells)
library(ranger)

#----function for assessing model performance----
perf_mets_calculator = function(pred=NULL, obs=NULL){
  
  #where "obs" and "pred" are vectors of observed and predicted values
  
  #vector of performance metrics names
  Metric = c('R2', 'MSE', 'RMSE', 'RRMSE', 'ME', 'RME')
  
  #calculate performance metrics if both predictions and observations are provided
  if(!is.null(pred)&!is.null(obs)){
    
    residuals = obs - pred
    
    rss = sum(residuals^2) #residual sum of squares
    tss = sum((obs-mean(obs))^2) #total sum of squares
    denom = abs(mean(obs))
    
    R2 = ifelse(tss==0, NA, 1-(rss/tss)) #coefficient of determination
    MSE = mean((obs-pred)^2)
    RMSE = sqrt(mean((obs-pred)^2)) #root mean squared error
    RRMSE = ifelse(denom==0, NA, 100*RMSE/abs(mean(obs))) #relative rmse / percent rmse
    ME = ifelse(denom == 0, NA, mean(residuals)) #mean error / bias
    RME = ifelse(denom == 0, NA, 100*ME/abs(mean(obs))) #relative mean error
    
    results = tibble(
      Metric = Metric,
      Value = c(R2, MSE, RMSE, RRMSE, ME, RME)
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

#----generate principal component change rasters----
pc_l = pblapply(unique(results_meta_df$dataset), function(d){
  
  block_l = lapply(block_ids, function(b){
    
    #----load pre-harvest and post-harvest data----
    pre_harv_f = unique(results_meta_df$file_path[results_meta_df$id == pre_harv_id & 
                                                    results_meta_df$block_id==b & 
                                                    results_meta_df$dataset==d])
    pre_harv_r = rast(pre_harv_f)
    pre_harv_r = pre_harv_r[[desired_feats]]
    
    post_harv_f = unique(results_meta_df$file_path[results_meta_df$id == post_harv_id & 
                                                     results_meta_df$block_id==b & 
                                                     results_meta_df$dataset==d])
    post_harv_r = rast(post_harv_f)
    post_harv_r = post_harv_r[[desired_feats]]
    
    #---get the difference raster----
    
    delta_r = post_harv_r - pre_harv_r
    
    #----scale difference raster----
    
    delta_r_scaled = scale(delta_r)
    
    #----generate principal component rasters----
    
    d_pc = prcomp(delta_r_scaled)
    d_pc_r = predict(delta_r, d_pc)
    
    #----save prcomp object and raster----
    
    res_l = list(d_pc, d_pc_r)
    names(res_l) = c('prcomp_object', 'prcomp_raster')
    
    return(res_l)
    
  })
  
  names(block_l) = block_ids
  return(block_l)
})
names(pc_l) = unique(results_meta_df$dataset)

#----plot cumulative proportion of variance explained by principal components----
{
  cum_var_df = lapply(names(pc_l), function(d){
    b_df = lapply(names(pc_l[[d]]), function(b){
      
      pc = pc_l[[d]][[b]][['prcomp_object']]
      summ = summary(pc)
      summ_imp_df = summ$importance |> 
        as.data.frame()
      summ_imp_df_long = summ_imp_df |>
        mutate(metric = rownames(summ_imp_df),
               dataset = d,
               block_id = b) |>
        pivot_longer(cols = colnames(summ_imp_df), names_to = 'PC_char')
    })
    return(b_df)
  }) |>
    bind_rows() |>
    mutate(PC = as.numeric(str_replace(PC_char,'PC','')))
  
  ggplot(cum_var_df |> filter(str_detect(metric, 'Cumulative'))) +
    geom_point(aes(x = PC, y = value)) +
    geom_line(aes(x = PC, y = value)) +
    facet_grid(cols = vars(dataset), rows = vars(block_id))+
    scale_x_continuous(breaks = c(unique(cum_var_df$PC)))+
    theme_classic()
}

#----plot PC feature rotations----
{
  rotations_df = lapply(names(pc_l), function(d){
    b_df = lapply(names(pc_l[[d]]), function(b){
      
      pc = pc_l[[d]][[b]][['prcomp_object']]
      rt_df = pc$rotation |>
        as.data.frame()
      rt_df_long = rt_df |>
        mutate(var_name = rownames(rt_df),
               dataset = d,
               block_id = b) |>
        pivot_longer(cols = colnames(rt_df), names_to = 'PC_chr') |>
        mutate(PC = as.numeric(str_replace(PC_chr,'PC','')))
    })
    return(b_df)
  }) |>
    bind_rows()
  
  feat_cols <- c(
    coastal_blue = "lightblue3",
    blue         = "blue",
    green_i      = "green4",
    green        = "lightgreen",
    yellow       = "#ffd92f",
    red          = "#e31a1c",
    rededge      = "#fb9a99",
    nir          = "#6a3d9a",
    NDVI         = "#000000"
  )
  
  rotations_df$var_name <- factor(
    rotations_df$var_name,
    levels = desired_feats
  )
  
  ggplot(rotations_df |> filter(PC %in% c(1:4)))+
    geom_point(aes(x = PC, y = value, color=var_name), alpha=0.7)+
    geom_line(aes(x = PC, y = value, color=var_name), alpha=0.7)+
    scale_x_continuous(breaks = c(1:3))+
    scale_color_manual(values = feat_cols) +
    facet_grid(cols=vars(dataset), rows=vars(block_id))+
    labs(x = 'PC', y = 'Rotation')+
    geom_hline(yintercept = 0, linetype='dotted')+
    theme_classic()
}

#----superpixel modelling----

library(supercells)
library(ranger)
library(exactextractr)

#compactness values for SLIC
cmpcts = c(0.1,1,5,10,15,20)

#step values (distance between initial cluster centres) for SLIC
stps = c(2, 4, 6, 8, 10)

pcs = 1:3

#get dataframe of datasets and hyperparameters to iterate over
sp_grid = rbind(
  expand.grid(cmpct = 0, stp = 0, dataset = names(pc_l)
              # , block_id = block_ids
              , stringsAsFactors = F), #pixel-based versions
  expand.grid(cmpct=cmpcts, stp=stps, dataset = names(pc_l)
              # , block_id = block_ids
              , stringsAsFactors = F) #superpixel iterations
) |>
  arrange(stp, dataset)

#do modelling, save results

sp_dir = file.path(ps_dir,paste0('scene_superpixel_modelling_PC_results_',lid_met))
dir.check(sp_dir)
sp_mods_dir = file.path(sp_dir,'RF_models')
dir.check(sp_mods_dir)
sp_results_dir = file.path(sp_dir, 'Results_CSVs')
dir.check(sp_results_dir)

pblapply(1:nrow(sp_grid), function(i){
  
  print(paste('Processing',i,'/',nrow(sp_grid)))
  
  stp = sp_grid$stp[i]
  cmpct = sp_grid$cmpct[i]
  dataset = sp_grid$dataset[i]
  # holdout_block = sp_grid$block_id[i] #holdout data used to test random forest model in cross-validation
  
  # res_string = paste0(dataset,'_testblock=',holdout_block,'_stp=',stp,'_cmpct=',cmpct,'_pcs=',max(pcs))
  res_strings = paste0(dataset,'_testblock=',block_ids,'_stp=',stp,'_cmpct=',cmpct,'_pcs=',max(pcs),'_BalanceBlocks')
  
  csv_folds_f = paste0(sp_results_dir,'/',res_strings,'.csv')
  
  if(any(!file.exists(csv_folds_f))){
    
    
    #generate superpixels for each block
    pix_vals_df = pblapply(block_ids, function(b){
      
      # print(b)
      
      #get planetscope PC raster
      r = pc_l[[dataset]][[b]]$prcomp_raster[[pcs]]
      
      #load lidar metric
      lid_f = list.files(lid_mets_change_cropped_dir, full.names = T)[str_detect(list.files(lid_mets_change_cropped_dir, full.names = T),b)]
      lid_change_met = rast(lid_f, lyrs = lid_met) #load lidar metric from file
      
      #crop PS to match lidar
      r = crop(r, lid_change_met, mask=T)
      
      #superpixel segmentation and stacking with lidar metric
      if(stp==0){ #if st==0 skip superpixel segmentation...
        
        vals_df = r |>
          c(lid_change_met)|>
          as.data.frame() |>
          mutate(block_id = b)
        names(vals_df)[ncol(vals_df)-1] = 'lid_met'
        
      } else { #else do superpixel segmentation
        
        sp = supercells(x = r, compactness = cmpct, step = stp)
        sp_lid = exact_extract(x = lid_change_met, y = sp, 'mean')
        sp = cbind(sp, sp_lid) |>
          rename(lid_met=sp_lid)
        vals_df = sp |>
          st_drop_geometry() |>
          select(-supercells, -x, -y) |>
          mutate(block_id = b)
      }
      
      return(vals_df)
    }) |>
      bind_rows() |>
      na.omit()
    
    #random forest cross validation
    
    blocks_npix = pix_vals_df |> #get number of pixels in each block to balance classes
      group_by(block_id) |> 
      summarise(n = n())
    
    pblapply(block_ids, function(b){
      
      res_string = res_strings[str_detect(res_strings,b)]
      
      csv_f = paste0(sp_results_dir,'/', res_string,'.csv')
      if(!file.exists(csv_f)){
        
        #check if RF file exists
        rds_f = paste0(sp_mods_dir,'/',res_string,'.rds')
        if(!file.exists(rds_f)){
          
          #sample size to get even sample from each block
          samp_size = min(blocks_npix$n[blocks_npix$block_id != b])
          
          
          #block to hold out
          training_data = pix_vals_df |>
            filter(block_id != b) |>
            #get same number of observations for each block
            group_by(block_id) |>
            slice_sample(n = samp_size) |>
            ungroup() |>
            #drop block_id column
            select(-block_id)
          
          testing_data = pix_vals_df |>
            filter(block_id == b) |>
            select(-block_id)
          
          #train random forest model
          rf = ranger(lid_met ~ ., data = training_data, importance = 'permutation'
                      , num.threads = 30)
          
          #save random forest model with relevant other data
          mod_l = list(rf, training_data, testing_data, pix_vals_df)
          names(mod_l) = c('ranger_object', 'training_data', 'testing_data', 'full_data')
          write_rds(mod_l, rds_f)
          
        } else {
          
          #load model if it already exists
          mod_l = read_rds(rds_f)
          rf = mod_l$ranger_object
          training_data = mod_l$training_data
          testing_data = mod_l$testing_data
          
        }
        
        #evaluate model with testing data
        test_res = predict(data = testing_data |> select(-lid_met), 
                           object = rf)
        test_preds = test_res$predictions
        
        #get model performance metrics for the training and testing data
        
        mod_perf_df = bind_rows(
          test_perf = rbind(
            perf_mets_calculator(obs = testing_data$lid_met, pred = test_preds),
            tribble( ~Metric, ~Value,
                     'OOB_R2', rf$r.squared,
                     'OOB_MSE', rf$prediction.error)
          ) |>
            mutate(test_type = 'holdout_block'),
          train_perf = perf_mets_calculator(pred = rf$predictions, obs = training_data$lid_met) |>
            mutate(test_type = 'training_blocks'),
          v.imp_df = tibble(Value = rf$variable.importance,
                            Metric = paste0(names(rf$variable.importance))
          )
        ) |>
          mutate(holdout_block = b,
                 dataset = dataset,
                 cmpct = cmpct,
                 stp = stp,
                 PCs = max(pcs))
        
        
        write_csv(mod_perf_df, csv_f)
      }
      
    })
  }
})


#----load and visualize results----

library(vroom)
# sp_rf_results_df = pblapply(list.files(sp_results_dir, full.names = T), read_csv) |> bind_rows()
sp_rf_results_df = vroom(list.files(sp_results_dir, full.names = T), delim=',') |>
  mutate(test_type = ifelse(str_detect(Metric, 'OOB'), 'training_blocks', test_type)) #switch OOB stats from holdout data to training data

r2_cv_summ_df = sp_rf_results_df |>
  filter(Metric == 'OOB_R2') |>
  group_by(test_type, dataset, cmpct, stp) |>
  summarise(mean_R2 = mean(Value))

#summary facetted geom_tile
ggplot(data = r2_cv_summ_df |> filter(test_type == 'training_blocks'),
       aes(x = as.factor(stp), y = as.factor(cmpct)))+
  geom_tile(aes(fill=mean_R2))+
  geom_text(aes(label=round(mean_R2,3)))+
  scale_fill_viridis_c()+
  # scale_y_continuous(limits = c(-1,0.25))+
  facet_wrap(vars(dataset))+
  theme_classic()

# #summary plot facetted scatter plot
# ggplot(data = r2_cv_summ_df |> filter(test_type == 'holdout_block'),
#        aes(x = stp, y = mean_R2))+
#   geom_point()+
#   geom_line()+
#   scale_y_continuous(limits = c(-1,0.25))+
#   facet_grid(cols = vars(dataset), rows = vars(cmpct))

# #plot of all cross validation vals
# ggplot(data = sp_rf_results_df |> filter(Metric == 'R2', test_type == 'holdout_block'))+
#   geom_boxplot(aes(x = stp, y = Value, group = stp))+
#   geom_point(aes(x = stp, y = Value, color = holdout_block))+
#   facet_grid(cols = vars(dataset), rows = vars(cmpct))
