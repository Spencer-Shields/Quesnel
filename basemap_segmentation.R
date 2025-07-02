#general basemap preprocessing
source('basemap_setup_preprocessing.R')

#glcm info and general preprocessing
source('make_basemap_GLCM.R')

#----define time period for each study block----

#fixed image dates
pre_harvest_date = '2021_07'
post_harvest_date = '2024_07'

#change image dates per block
{
  #define growing season
  growing_season_months = c('05', '06', '07', '08', '09')
  
  #define function to round up to next growing season image
  round_to_growing_season <- function(date_vec) {
    date_vec <- as.Date(date_vec)
    year_vec <- year(date_vec)
    month_vec <- month(date_vec)
    
    result <- ifelse(
      month_vec %in% 5:9,
      ceiling_date(date_vec, "month") - days(1),
      as.Date(paste0(ifelse(month_vec >= 10, year_vec + 1, year_vec), "-05-31"))
    )
    
    return(as.Date(result, origin = "1970-01-01"))
  }
  
  #get images to use for each block
  harvest_dates_df = harvest_dates_df |>
    mutate(post_harvest_image_date = round_to_growing_season(harvest_finish_date))%>%
    mutate(
      # Try same month-day as post-harvest, but in current year
      pre_harvest_image_date_candidate = as.Date(paste0(year(harvest_start_date), "-", 
                                                        month(post_harvest_image_date), "-", 
                                                        day(post_harvest_image_date))),
      
      # If that date is after harvest_start_date, then use previous year
      pre_harvest_image_date = if_else(
        pre_harvest_image_date_candidate < harvest_start_date,
        pre_harvest_image_date_candidate,
        as.Date(paste0(year(harvest_start_date) - 1, "-", 
                       month(post_harvest_image_date), "-", 
                       day(post_harvest_image_date)))
      )
    ) %>%
    select(-pre_harvest_image_date_candidate)|>  #remove helper column
    mutate(post_harvest_image_string = format(post_harvest_image_date, '%Y_%m'),
           pre_harvest_image_string = format(pre_harvest_image_date, '%Y_%m'))
}

#----define wrapper function for processing data ----

#-----define wrapper function for processing data with GAM-----
planetscope_segmentation_fn = function(
    sp_size = 5 #target size of superpixels
    # sp_step = 3
    , sp_compactness = 0.1 #superpixel compactness
    , sp_avg_fn = 'mean' #aggregation function for superpixels
    , km_centers = 2 #number of kmeans clusters
){
  
  gam_results_df_l = pblapply(c(nonorm_dir, z_dir, zr_dir),function(d){
    #-----loading, pre-processing, raster differencing-----
    {
      #----load rasters----
      cat('\n----Loading and preparing PlanetScope data')
      ps_directory = d
      ps_flist = list.files(ps_directory, full.names = T, recursive = T, pattern = '\\.tif$')
      
      block_ids = blocks$BLOCKNUM[!str_detect(blocks$BLOCKNUM, 'NoChange')]
      
      pre_harvest_rasts = lapply(block_ids, function(x){ #single reference date
        block_paths = ps_flist[str_detect(ps_flist, x)]
        r = rast(block_paths[str_detect(block_paths, pre_harvest_date)])
        return(r)
      })
      # pre_harvest_rasts = lapply(block_ids, function(x){ #dynamic reference date
      #   block_paths = ps_flist[str_detect(ps_flist, x)]
      #   r = rast(block_paths[str_detect(block_paths, 
      #                                   harvest_dates_df$pre_harvest_image_string[harvest_dates_df$block_id==x])])
      #   return(r)
      # })
      names(pre_harvest_rasts) = block_ids
      
      post_harvest_rasts = lapply(block_ids, function(x){ #single reference date
        block_paths = ps_flist[str_detect(ps_flist, x)]
        r = rast(block_paths[str_detect(block_paths, post_harvest_date)])
        return(r)
      })
      # post_harvest_rasts = lapply(block_ids, function(x){ #dynamic post-harvest image date
      #   block_paths = ps_flist[str_detect(ps_flist, x)]
      #   r = rast(block_paths[str_detect(block_paths, 
      #                                   harvest_dates_df$post_harvest_image_string[harvest_dates_df$block_id==x])])
      #   return(r)
      # })
      names(post_harvest_rasts) = block_ids
      
      #----add GLCM layers to spectral layers----
      
      window_size = '3x3' #define window size of raster desired
      
      glcm_dataset = case_when(
        str_detect(ps_directory, 'Zrobust') ~ 'Z_robust',
        str_detect(ps_directory, 'Z') ~ 'Z',
        .default = 'Non-normalized'
      )
      glcm_dataset_dir = file.path(glcm_dir, glcm_dataset)
      glcm_flist = list.files(file.path(glcm_dataset_dir, window_size), recursive = T, pattern = '\\.tif$', full.names = T)
      
      pre_harvest_rasts = pblapply(block_ids, function(b){
        f = glcm_flist[str_detect(ps_flist, b) & str_detect(ps_flist, pre_harvest_date)]
        gr = rast(f)
        r = c(pre_harvest_rasts[[b]], gr)
        return(r)
      })
      names(pre_harvest_rasts) = block_ids
      
      post_harvest_rasts = pblapply(block_ids, function(b){
        f = glcm_flist[str_detect(ps_flist, b) & str_detect(ps_flist, post_harvest_date)]
        gr = rast(f)
        r = c(post_harvest_rasts[[b]], gr)
        return(r)
      })
      names(post_harvest_rasts) = block_ids
      
      #----calculate difference rasters----
      difference_rasts = mapply(function(x,y){x-y}, post_harvest_rasts, pre_harvest_rasts)
      
    }
    
    #----extract validation data----
    {
    cat('\n-----Loading and extracting validation data')
    #load and stack CHM
    chm_dir = 'data/Quesnel_thinning/chm_change'
    chm_files = list.files(chm_dir, full.names = T, pattern = '\\.tif$')
    
    #load canopy cover
    cc_dir = "data/quesnel_thinning_las/canopy_cover_change"
    cc_files = list.files(cc_dir, full.names = T, pattern = '\\.tif$')
    
    #stack validation with other raster data
    difference_rasts = pblapply(names(difference_rasts), function(n){
      
      #get satellite block raster
      dr = difference_rasts[[n]]
      
      #make mask for fixing geometry
      m = ifel(is.na(dr[[1]]),NA,1)
      
      #prep and stack chm
      chm_f = chm_files[str_detect(chm_files, n)]
      chm_r = rast(chm_f)
      
      if(crs(chm_r) != crs(m)){ #reproject
        chm_r = project(chm_r, m)
      }
      chm_r = mask(chm_r,m)
      names(chm_r) = 'CHM_validation'
      dr = c(dr, chm_r)
      
      #prep and stack cc
      cc_f = cc_files[str_detect(cc_files, n)]
      cc_r = rast(cc_f)
      
      if(crs(cc_r) != crs(m)){ #reproject
        cc_r = project(cc_r, m)
      }
      cc_r = mask(cc_r,m)
      names(cc_r) = 'CC_validation'
      dr = c(dr, cc_r)
      
      
      return(dr)
    })
    names(difference_rasts) = block_ids
    
    
    }
    
    #----feature extraction----
    {
      #----correlations----
      r = difference_rasts$`12N_T3`
      
      
      library(corrplot)
      
      
      cor_mat = r |> values(na.rm=T) |> cor() |> corrplot()
      
      cor_mat = r[[feats]] |> values(na.rm=T) |> cor() |> corrplot(
        # order='hclust'
        )
      cor_mat

      cor_mat = r[[c('blue','green','red','BI','ND')]] |> values(na.rm=T) |> cor() |> corrplot(
          # order='hclust'
          )

      
      #----get subset of features for processing----
      full_feats = names(difference_rasts[[1]])
      spectral_feats = c('blue', 'green', 'red'
                # , 'BI'
                )
      glcm_feats = full_feats[str_detect(full_feats, 'glcm_')]
      glcm_feats
      
      
      pattern = str_c(c(spectral_feats
                        # , paste0('_',bands)
                        ), collapse = '|') #all original bands, brightness index, and GLCM layers
      
      feats = full_feats[str_detect(full_feats, pattern)]
      feats = feats[!str_detect(feats, 'VARIgreen')]
      
      difference_rasts = pblapply(difference_rasts, function(r)r[[feats]])
      
      #----PCA----
      
      pc_l = pblapply(difference_rasts, function(r){
        #run PCA
        pc_r = r |>
          scale() |>
          prcomp()
        return(pc_r)
      })
      names(pc_l) = names(difference_rasts)
      
      
      #get dataframe of pc standard deviations
      sdev_df <- map_dfr(pc_l, ~ .x$sdev, .id = "ID") %>%
        mutate(PC = 1:20) |>
        pivot_longer(-PC, names_to = "block_id", values_to = "SD") %>%
        mutate(PC = as.integer(str_remove(PC, "^V"))) |> # optional: turn PC names into integers
        group_by(block_id) %>%
        mutate(
          variance = SD^2,
          prop_variance = variance / sum(variance)
        ) %>%
        mutate(cum_variance = cumsum(prop_variance)) %>%
        ungroup()
      
      #plot cumulative standard deviations
      ggplot(sdev_df, aes(x = PC, y = cum_variance)) +
        geom_col() +
        facet_wrap(~block_id) +
        scale_y_continuous(breaks = seq(0, 1, by = 0.1)) +   # ticks every 0.1
        labs(
          title = "Cumulative Proportion of Variance Explained",
          x = "Principal Component",
          y = "Cumulative Variance Explained"
        ) +
        theme_minimal()
      
      
      
      library(tidytext)
      # Extract first 9 PCs
      rot <- pc$rotation[, 1:9]
      
      # Convert to long format
      rotation_long <- as.data.frame(rot) %>%
        rownames_to_column("Feature") %>%
        pivot_longer(-Feature, names_to = "PC", values_to = "Loading") %>%
        mutate(AbsLoading = abs(Loading))
      
      # Reorder features within each facet
      rotation_long <- rotation_long %>%
        mutate(Feature = reorder_within(Feature, AbsLoading, PC))
      
      # Plot
      ggplot(rotation_long, aes(x = Feature, y = AbsLoading)) +
        geom_col() +
        facet_wrap(~ PC, scales = "free_x") +
        scale_x_reordered() +
        theme_minimal() +
        theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 7)) +
        labs(x = "Feature", y = "Absolute Loading", title = "Feature Contributions to Principal Components")
      
      
      pc_rasts = pblapply(block_ids, function(b){
        dif_rast = difference_rasts[[b]]
        pc = pc_l[[b]]
        pc_rast = predict(dif_rast, pc)
        return(pc_rast)
        })
      names(pc_rasts) = block_ids
  
      
      
      plot(pc_rasts$`12N_T3`[[1:5]])
      
      
      ggplot(pc_rasts$`12N_T3` |> values(), aes(x = PC1, y = PC2))+
        geom_point(alpha = 0.3)
      
      ggplot(difference_rasts2$`12N_T3` |> values())+
        geom_histogram(aes(x = red), alpha = 0.2, fill = 'red')+
        geom_histogram(aes(x = blue), alpha = 0.2, fill = 'blue')+
        geom_histogram(aes(x = green), alpha = 0.2, fill = 'green')+
        geom_histogram(aes(x = CHM_validation), alpha = 0.2, fill = 'black')
      
      
      ggplot(difference_rasts2$`12N_T3` |> values())+
        geom_density(aes(x = red), alpha = 0.2, fill = 'red', color = 'red')+
        geom_density(aes(x = blue), alpha = 0.2, fill = 'blue', color = 'blue')+
        geom_density(aes(x = green), alpha = 0.2, fill = 'green', color = 'green')+
        geom_density(aes(x = CHM_validation), alpha = 0.2, fill = 'black', color = 'black')
      
      
      ggplot(difference_rasts$`12N_T3` |> values())+
        geom_point(aes(x = red, y = green), alpha = 0.1)
      
    }
    
    #----subset each raster to only include desired bands----
    difference_rasts = lapply(difference_rasts, function(r)r[[c('blue', 'green', 'red','Hue', 'GCC', 'NDGR', 'BI', 'CI', 'CRI550', 'GLI')]])
    
    bands_to_drop = c('max_DN','TVI','VARIgreen')
    difference_rasts = lapply(difference_rasts, function(r)r[[!names(r) %in% bands_to_drop]])
    
    feats = names(difference_rasts[[1]])
    
    #----superpixel clustering----
    cat('\n-----Generating superpixels')
    if(!'package:supercells' %in% search()){library(supercells)}
    
    #define superpixel clustering hyperparameters
    
    # sp_step = sp_step
    compactness = sp_compactness #how square/regular superpixels are
    avg_fun = sp_avg_fn
    
    #define number of superpixels there should be in each block to have desired superpixel size
    sp_size = sp_size #total number of cells per supercell
    n_pixels_blocks = sapply(difference_rasts, ncell) #total pixels in each raster
    k_blocks = 
      # ceiling(
      n_pixels_blocks/sp_size
    # )
    
    #generate superpixels for each block
    superpixels_l = pblapply(1:length(difference_rasts), function(i){
      r = difference_rasts[[i]]
      k = k_blocks[i] |> unname()
      # step = sp_step
      sp = supercells(x=r
                      , k=k
                      # ,step = step
                      , compactness=compactness, avg_fun = avg_fun
                      # , verbose =1
      )
      sp$block_id = block_ids[i]
      return(sp)
    })
    names(superpixels_l) = block_ids
    
    # #----PCA----
    # pc = sp |>
    #   st_drop_geometry() |>
    #   select(all_of(feats)) |>
    #   scale() |>
    #   prcomp()
    # 
    # pc_var=pc$sdev^2
    # pc_var_prop = pc_var/sum(pc_var)
    # 
    # pc_var_tbl = tibble(PC = 1:length(pc_var),
    #                     proportion_of_variance = pc_var_prop,
    #                     cumulative_proportion_of_variance = cumsum(pc_var_prop))
    
    #----k-means segmentation----
    
    library(cluster)
    
    superpixels_km_l = pblapply(superpixels_l, function(sp){
      km = sp |>
        st_drop_geometry() |>
        select(all_of(feats)) |>
        scale() |>
        kmeans(centers = km_centers, iter.max = 50)
      km
      sp$cluster = km$cluster
      sp$
      return(sp)
    })
    
    #----extract validation data----
    cat('\n-----Loading and extracting validation data')
    #load CHM
    chm_dir = 'data/Quesnel_thinning/chm_change'
    chm_files = list.files(chm_dir, full.names = T, pattern = '\\.tif$')
    
    # chm_l = lapply(block_ids, function(x){
    #   chm_file = chm_files[str_detect(chm_files, x)]
    #   rast(chm_file)
    # })
    # names(chm_l) = block_ids
    
    #load canopy cover
    cc_dir = "data/quesnel_thinning_las/canopy_cover_change"
    cc_files = list.files(cc_dir, full.names = T, pattern = '\\.tif$')
    
    #prepare superpixels for processing
    superpixels_l_prepped = pblapply(superpixels_l, function(x){
      # sp = x
      sp = vect(x)
      project(sp, crs(rast(chm_files[[1]])))
      sp = wrap(sp)
      return(sp)
    })
    
    #extract validation info
    # plan('multisession', workers = 8)
    if(!'package:exactextractr' %in% search()){library(exactextractr)}
    
    superpixels_validation_l = pblapply(block_ids, function(b){
      
      #load superpixels block
      sp = unwrap(superpixels_l_prepped[[b]])
      
      #load canopy height
      # chm = chm_l[[b]]
      chm_f = chm_files[which(str_detect(chm_files, b))]
      chm = rast(chm_f)
      names(chm) = 'canopy_height'
      
      #load canopy cover
      cc_f = cc_files[which(str_detect(cc_files, b))]
      cc = rast(cc_f)
      names(cc) = 'canopy_cover'
      
      #stack layers
      if(ext(cc)!=ext(chm)|crs(cc)!=crs(chm)){cc = project(cc, chm)}
      validation_lyrs = c(chm,cc)
      
      #get weighted mean values for each superpixel
      # cat('Extracting',b,'validation data\n')
      extracted_vals = exact_extract(x = validation_lyrs, y = st_as_sf(sp), 'mean'
                                     ,progress = F
      )
      sp_validation = cbind(sp, extracted_vals)
      # cat(b,'finished\n')
      
      return(sp_validation)
    })
    names(superpixels_validation_l) = block_ids
    
    #make single sf object with all superpixels objects
    superpixels_validation_combined = bind_rows(lapply(superpixels_validation_l,st_as_sf))
    
    # plan('sequential')
    
    
    
    #----fit models, extract model fit metrics----
    cat('\n-----Fitting GAMs with cross-validation')
    if(!'package:mgcv' %in% search()){library(mgcv)}
    
    
    
    #write function to do k-fold cross validation with generalized additive models
    gam_cv = function(data, target_var, predictors = colnames(data)[colnames(data)!=target_var], kfolds=5, ...){
      
      #partition data into k folds
      folds = createFolds(y = data[[target_var]], k = kfolds)
      
      #generate list of gams using each fold
      gam_l = pblapply(folds, function(fold){
        
        #generate model
        data_train = data[-fold,]
        smooth_terms = smooth_terms <- paste0("s(", predictors, ")")
        form = as.formula(
          paste(target_var,'~',paste(smooth_terms, collapse = ' + ')))
        
        #fit model on training data
        gam_mod = bam(form, data=data_train, method='REML',...)
        
        #get model performance metrics based on training data
        s_train = summary(gam_mod) #r2
        r2_train = s_train$r.sq
        train_predictions = gam_mod$fitted.values
        train_actual = data_train[[target_var]]
        bias_train = mean(train_actual-train_predictions)#bias
        rmse_train = sqrt(mean((train_actual-train_predictions)^2)) #rmse
        nrmse_train = rmse_train/IQR(train_actual) #nrmse
        
        #get r2 of model applied to testing data
        data_test = data[fold,]
        test_predictions = predict.bam(gam_mod, newdata = data_test) |> as.numeric()
        test_actual = data_test[[target_var]]
        
        ss_res <- sum((test_actual - test_predictions)^2) 
        ss_tot <- sum((test_actual - mean(test_actual))^2)
        r2_test <- 1 - ss_res / ss_tot
        
        bias_test = mean(test_actual-test_predictions)#bias
        rmse_test = sqrt(mean((test_actual-test_predictions)^2)) #rmse
        nrmse_test = rmse_test/IQR(test_actual) #nrmse
        
        
        #return results and model performance stats
        gm_l = list(gam_mod, r2_test, rmse_test, nrmse_test, bias_test, r2_train, rmse_train, nrmse_train, bias_train)
        names(gm_l) =  c('gam.model', 'r.sqr_test', 'rmse_test', 'nrmse_test', 'bias_test', 'r.sqr_train', 'rmse_train', 'nrmse_train', 'bias_train')
        
        return(gm_l)
      })
      
      metrics <- c('r.sqr_test', 'rmse_test', 'nrmse_test', 'bias_test', 
                   'r.sqr_train', 'rmse_train', 'nrmse_train', 'bias_train')
      mean_metrics <- setNames(
        lapply(metrics, function(metric) {
          mean(sapply(gam_l, function(x) x[[metric]]))
        }),
        metrics
      )
      
      cv_results <- c(mean_metrics, list(kfold_results = gam_l))
      
      return(cv_results)
    }
    
    #get gam models and model performance for canopy height and canopy cover for each block
    set.seed(123)
    gam_results_l = pblapply(block_ids, function(b){
      # cat('Modelling for',b,'\n')
      sp_df = superpixels_validation_combined |>
        filter(block_id == b) |>
        st_drop_geometry() |>
        select(-c(supercells,x,y,block_id))|>
        na.omit()
      # cat('-----canopy cover\n')
      cc_results = gam_cv(data = sp_df|>select(-mean.canopy_height), target_var = 'mean.canopy_cover', kfolds = 5)
      # cat('-----canopy height\n')
      ch_results = gam_cv(data = sp_df |>select(-mean.canopy_cover), target_var = 'mean.canopy_height', kfolds = 5)
      
      results_l = list(cc_results, ch_results)
      names(results_l) = c('canopy_cover', 'canopy_height')
      return(results_l)
    })
    names(gam_results_l) = block_ids
    
    #extract results from gams
    extract_gam_metrics <- function(gam_results_l) {
      results_list <- list()
      
      for (site_id in names(gam_results_l)) {
        site_data <- gam_results_l[[site_id]]
        
        for (target_var in names(site_data)) {
          model_data <- site_data[[target_var]]
          
          # Extract overall train/test metrics
          overall_metrics <- data.frame(
            site = site_id,
            target = target_var,
            fold = "overall",
            r2_test = model_data$r.sqr_test,
            rmse_test = model_data$rmse_test,
            nrmse_test = model_data$nrmse_test,
            bias_test = model_data$bias_test,
            r2_train = model_data$r.sqr_train,
            rmse_train = model_data$rmse_train,
            nrmse_train = model_data$nrmse_train,
            bias_train = model_data$bias_train,
            stringsAsFactors = FALSE
          )
          
          results_list[[length(results_list) + 1]] <- overall_metrics
          
          # Extract k-fold results
          kfolds <- model_data$kfold_results
          if (!is.null(kfolds)) {
            for (fold_name in names(kfolds)) {
              fold_data <- kfolds[[fold_name]]
              fold_metrics <- data.frame(
                site = site_id,
                target = target_var,
                fold = fold_name,
                r2_test = fold_data$r.sqr_test,
                rmse_test = fold_data$rmse_test,
                nrmse_test = fold_data$nrmse_test,
                bias_test = fold_data$bias_test,
                r2_train = fold_data$r.sqr_train,
                rmse_train = fold_data$rmse_train,
                nrmse_train = fold_data$nrmse_train,
                bias_train = fold_data$bias_train,
                stringsAsFactors = FALSE
              )
              results_list[[length(results_list) + 1]] <- fold_metrics
            }
          }
        }
      }
      
      # Combine all results into a single dataframe
      results_df <- do.call(rbind, results_list)
      return(results_df)
    }
    
    gam_results_df = extract_gam_metrics(gam_results_l) |> 
      rename(block_id=site) |>
      left_join(harvest_dates_df) |>
      mutate(ps_directory = ps_directory)
    
    cat('\n',basename(d),'finished')
    return(gam_results_df)
  })
  
  #---combine list of model results into dataframe----
  gam_results_df = bind_rows(gam_results_df_l) |>
    mutate(dataset = basename(ps_directory)) |>
    mutate(dataset = case_when(
      str_detect(basename(ps_directory),'Zrobust') ~ 'Zrobust',
      str_detect(basename(ps_directory),'Z') ~ 'Z',
      .default = 'Non-normalized'
    ))
  
  return(gam_results_df)
  
}