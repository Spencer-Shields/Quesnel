source('scene_global_stats_filtering_20250925.R')
# library(mgcv)
# library(dbscan)
# library(cluster)
# library(splitstackshape)
library(Rfast)

# #get raster timeseries files to work with
# d = paste0(dirs[1],'/12L_C7')
# fl = list.files(d, recursive = T, full.names = T, pattern = '\\.tif$')
# fl = fl[!str_detect(fl, 'NULL')]
# 
# feat = 'NDGR'

#----define function to perform GAM clustering----

ITS_pixel_modeller = function(fl, feat, 
                         ){
  #dir is a directory containing a timeseries of PlanetScope scenes, feat is the feature to be extracted (string representing the name of a layer in the raster dataset)
  #minPts is the tuning parameter for the hdbscan algorithm, NULL means that it will be optimized
  #min_dates is the minimum number of dates that a pixel must have in the timeseries to have a gam fit
  #smooth type is the type of basis function used in gam fitting (see mgcv documentation)
  #check_minPts means that you want to run hdbscan over a range of minPts values (defined by min and max_minPts)
  #gs months is a vector of months to fit the time series to
  #notNA_prop is the threshold proportion of notNA cells that a raster must have to be included 
  
  
  r_l = lapply(1:length(fl), function(i){
    
    if(verbose==T){
      if(i==1){
        cat('\nGetting raster timeseries layers',i,'/',length(fl))
      } else {
        cat('\rGetting raster timeseries layers',i,'/',length(fl))
      }}
    
    #load raster
    r = rast(fl[i])
    
    #check if feature exists, proceed if yes
    g = global(r, c('isNA', 'notNA'), na.rm=T)
    
    if(any(str_detect(names(r),feat))){
      if(notNA_prop*ncell(r)<g$notNA[rownames(g)==feat]){
        r_f = r[[feat]]
        
        #get raster acquisition date/time, make layer name the date/time
        ind = which(sapply(meta_df$id, function(x) grepl(x, fl[i]))) #find index of scene id in filename
        ac_date = meta_df$acquired[ind]
        names(r_f) = ac_date
        
        return(r_f)
      }
    }
  })
  
  names(r_l) = sapply(r_l, names)
  
  #remove Nulls from list
  r_l = Filter(Negate(is.null),r_l)
  
  #make list of rasters into a single spatraster object
  {
    #get matrix of raster extents
    ext_matrix <- t(sapply(r_l, function(r) {
      e <- ext(r)
      c(xmin(e), xmax(e), ymin(e), ymax(e))
    }))
    colnames(ext_matrix) <- c("xmin", "xmax", "ymin", "ymax")
    
    #find max and min extent values for largest possible extent, make spatextent spatraster object
    overall_extent <- c(
      xmin = min(ext_matrix[, "xmin"]),
      xmax = max(ext_matrix[, "xmax"]),
      ymin = min(ext_matrix[, "ymin"]),
      ymax = max(ext_matrix[, "ymax"])
    ) |>
      ext()
    #make template raster that the others should match
    template_r = rast(overall_extent, 
                      resolution = c(res(r_l[[1]])[1], res(r_l[[1]])[2]),
                      crs = crs(r_l[[1]])
    )
    #resample rasters to have consistent geometry
    r_l_resampled = lapply(1:length(r_l), function(i){
      if(verbose==T){cat('\rMaking geometry of layers consistent',i,'/',length(r_l))}
      
      r = r_l[[i]]
      if(!compareGeom(r, template_r, stopOnError=F)){
        r = resample(r, template_r)
      }
      return(r)
    })
    
    #combine rasters into a single spatraster object
    r = rast(r_l_resampled)
  }
  
  #get matrix of values in spatraster
  
  v_df = as.data.frame(r, xy=T) |>
    mutate(xy = paste(x,y))|>
    select(-x,-y) 
  v_df_long = v_df |>
    pivot_longer(cols = names(v_df)[!str_detect(names(v_df),'xy')], values_to = 'val', names_to = 'acquired') |>
    mutate(acquired = ymd_hms(acquired)) |>
    mutate(acquired_numeric = as.numeric(acquired))|>
    mutate(acquisition_month = month(acquired)) |>
    na.omit()
  
  #subset time series by gs_months
  v_df_long_gs = v_df_long |>
    filter(acquisition_month %in% gs_months)
  
  # ggplot(v_df_long)+
  #   geom_line(aes(x = acquired, y = val, group=xy))
  # 
  # ggplot(v_df_long_gs)+
  #   geom_line(aes(x = acquired, y = val, group=xy))
  
  pixel_ids = unique(v_df_long_gs$xy)
  
  #fit a gam to the time series for each pixel
  gam_l = lapply(1:length(pixel_ids), function(i){
    
    pixel_id = pixel_ids[i]
    
    if(verbose==T){cat('\rFitting GAMs',i,'/',length(pixel_ids),' ')}
    
    vals = v_df_long_gs |> filter(xy == pixel_ids[i])
    
    if(nrow(vals) > min_dates){
      mod = gam(val ~ s(acquired_numeric, bs=smooth_type), data = vals)
    } else {
      mod = NULL
    }
    return(mod)
  })
  names(gam_l) = pixel_ids
  
  gam_l = Filter(Negate(is.null),gam_l)
  
  # basis_n = sapply(gam_l, function(m){
  #   m$smooth[[1]]$bs.dim
  #   }) #all gams have 10 basis functions
  
  #extract gam coefficients, scale them
  coefs_df <- gam_l |>
    lapply(function(a) a$coefficients) |>
    do.call(rbind, args = _) |>
    `colnames<-`(c('Intercept', paste0('Coef_', 1:(length(gam_l[[1]]$coefficients)-1)))) |>
    scale() |>
    as.data.frame()
  
  #try hdbscan with different values for minPts
  
  if(check_minPts==T|is.null(minPts)){
    
    # if(verbose==T){cat('Running HDBSCAN for range of minPts')}
    minPts_to_test = min_minPts:max_minPts
    db_clust_l = lapply(minPts_to_test, function(i){
      if(verbose==T){cat('\rTesting HDBSCAN for minPts ==',i)}
      hdbscan(x = coefs_df, minPts = i)
    })
    
    #extract number of clusters, number of noise points, and average silhouette scores for each hdbscan iteration
    clust_eval_df = lapply(1:length(db_clust_l), function(i){
      if(verbose==T){cat('\rExtracting HDBSCAN cluster quality',i,'/',length(db_clust_l))}
      
      cl = db_clust_l[[i]]
      
      #get number of noise points
      noise_pts = sum(cl$cluster==0)
      
      #get number of clusters
      n_clusts = length(unique(cl$cluster[cl$cluster!=0]))
      
      #get mean cluster scores
      mean_clust_score = mean(cl$cluster_scores)
      
      #mean outlier scores
      mean_outlier_score = mean(cl$outlier_scores)
      
      # Calculate silhouette scores (excluding noise points)
      valid_clusters <- cl$cluster != 0
      
      if(sum(valid_clusters) > 0 && length(unique(cl$cluster[valid_clusters])) > 1) {
        sil <- silhouette(cl$cluster[valid_clusters], 
                          dist(coefs_df[valid_clusters, ]))
        
        # Overall mean silhouette score
        mean_sil <- mean(sil[, "sil_width"])
        # print(paste("Mean Silhouette Score:", round(mean_sil, 4)))
      } else {
        # print("Cannot compute silhouette: insufficient clusters or data")
        mean_sil = NA
      }
      
      # Return as data.frame row
      data.frame(
        Noise_points = noise_pts,
        Clusters = n_clusts, 
        Cluster_scores = mean_clust_score,
        Outlier_scores = mean_outlier_score,
        Silhouette = mean_sil
      )
    }) %>% 
      bind_rows() %>%
      as.matrix() |>
      scale() %>%
      as.data.frame()
    
    clust_eval_df$minPts = minPts_to_test
    
    clust_eval_df = clust_eval_df |> na.omit()
    
    clust_eval_df_long = clust_eval_df |>
      pivot_longer(cols = names(clust_eval_df)[names(clust_eval_df)!='minPts'],
                   names_to = 'metric',
                   values_to = 'val')
    
    
    clust_opt_plot = ggplot(clust_eval_df_long, aes(x = minPts, y = val, color = metric))+
      geom_point()+
      geom_line()
  }
  
  #set minPts value, use optimized value by default
  if(is.null(minPts)){
    mp = clust_eval_df$minPts[max(clust_eval_df$Cluster_scores)]
  } else {
    mp = minPts
  }
  
  if(verbose==T){cat('\rPerforming HDBSCAN clustering for minPts==',mp)}
  #run hdbscan
  db_clust = hdbscan(x = coefs_df
                     , minPts = mp
  )
  
  #join cluster to v_df_long
  v_df_long_clust = v_df_long_gs |>
    left_join(tibble(cluster = db_clust$cluster, xy = rownames(coefs_df)))
  
  clustered_time_plot = ggplot(v_df_long_clust, aes(x = acquired, y = val, group = xy, color=as.factor(cluster)))+
    geom_line(data=v_df_long_clust|>filter(cluster==0),alpha = 0.05)+
    geom_line(data = v_df_long_clust|> filter(cluster!=0), alpha=0.1)
  
  if(verbose==T){cat('\rMaking cluster raster')}
  #make cluster raster
  r_clust = v_df_long_clust |>
    separate_wider_delim(xy, names = c('x', 'y'), delim = ' ') |>
    select(x,y,cluster)|>
    distinct() |>
    rast(type='xyz')
  
  if(verbose==T){'\rGenerating output'}
  
  #save results to an output nested list structure
  input_parameters = list(fl, feat, minPts, min_dates, smooth_type)
  names(input_parameters) = c('file_list', 'feature', 'minPts', 'min_dates', 'smooth_type')
  
  if(check_minPts==T){
    clust_info = list(clust_eval_df_long, clust_opt_plot)
    names(clust_info) = c('HDBSCAN_model_performance', 'plot_HDBSCAN_optimization')
  } else {clust_info = NULL}
  
  results_list = list(r_clust, r, gam_l, coefs_df, db_clust, v_df_long_clust, clustered_time_plot, clust_info, input_parameters)
  names(results_list) = c('cluster_raster', 'timeseries_raster', 'GAM_list', 'GAM_coefs_df', 'HDBSCAN_model', 'pixel_cluster_timeseries_df', 'plot_pixel_cluster_timeseries', 'HDBSCAN_optimization_results', 'input_parameters')
  
  return(results_list)
}

# tc_ndvi_12lc7 = trend_cluster(fl = fl, feat='NDVI')
# 
# tc_hue_12lc7 = trend_cluster(fl = fl, feat='Hue')

#test NDGR timeseries for each combination of dataset and block id

#----generate per-pixel interrupted time series models for each feature, block, and dataset----

#filter results dataframe to remove ids that occur during harvesting
gs_months = c(5:10)

results_meta_df_prepostharv_gs = results_meta_df |>
  mutate(harvest_in_progress = ifelse(acquisition_date >= harvest_start_date & acquisition_date < harvest_finish_date,
                                      1, 0)) |>
  filter(harvest_in_progress == 0) |>
  # filter(acquisition_month %in% c(5,6,7,8,9)) |> #only model for growing season months
  #define dummy variable that measures time before harvest start and after harvest end
  mutate(time_since_harvest = ifelse(acquisition_date < harvest_start_date,
                                     as.numeric(acquisition_date - harvest_start_date),
                                     as.numeric(acquisition_date - harvest_finish_date))) |>
  filter(acquisition_month %in% gs_months)

ts_feats = unique(results_meta_df_prepostharv_gs$var_name)

ts_clust_l = pblapply(ts_feats, function(f){
  cat('Processing',f)
  ds_l = pblapply(dirs, function(d){
    
    #get list of block directories
    bd = list.dirs(d)[list.dirs(d)!=d]
    
    # bd = bd[str_detect(bd, '12L_D345')
    #         # |str_detect(bd, '12N_T3')
    # ]
    
    cat('---Processing',d)
    
    
    b_l = pblapply(bd, function(b){
      
      #get list of files in blockdir, remove files marked as invalid
      cat('------Processing',basename(b))
      
      # fl = list.files(b, recursive = T, full.names = T, pattern = '\\.tif$')
      # # fl = fl[!str_detect(fl, 'NULL')]
      # fl = fl[fl %in% results_meta_df_prepostharv_gs]
      
      fl = unique(results_meta_df_prepostharv_gs$file_path[str_detect(results_meta_df_prepostharv_gs$file_path,b)])
      
      dat = results_meta_df_prepostharv_gs[str_detect(results_meta_df_prepostharv_gs$file_path,b),] |>
        select(file_path, time_since_harvest, harvested) |> 
        distinct()
      
      #stack_rasters
      r = pblapply(dat$file_path, function(x)rast(x,lyrs=feat)) |> rast()
      names(r) = dat$time_since_harvest
      
      # v = values(r,na.rm=T)
      # vt = t(v)
      
      
      
      v = as.data.frame(r, xy=T)
      
      rownames(v) = paste(v$x, v$y)
      v = v |> select(-x,-y)
      
      # t = trend_cluster(fl=fl, feat=f, smooth_type = 'gp', notNA_prop = 0.1, verbose=T)
      return(t)
    })
    names(b_l) = basename(bd)
    return(b_l)
  })
  names(ds_l) = basename(dirs)
  return(ds_l)
})
names(ts_clust_l) = ts_feats