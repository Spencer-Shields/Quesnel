#----load packages, set up environment----
source('basemap_setup_preprocessing.R')

#----superpixel analysis
{
  
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
  
  superpixel_modelling_fn = function(sp_size, sp_compactness){
  
  gam_results_df_l = pblapply(c(nonorm_dir, z_dir, zr_dir),function(d){
    #----load rasters----
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
    
    #----calculate difference rasters----
    difference_rasts = mapply(function(x,y){x-y}, post_harvest_rasts, pre_harvest_rasts)
    
    #----subset each raster to only include desired bands----
    difference_rasts = lapply(difference_rasts, function(r)r[[c('blue', 'green', 'red','Hue', 'GCC', 'NDGR', 'BI', 'CI', 'CRI550', 'GLI')]])
    
    #----superpixel clustering----
    library(supercells)
    
    #define total number of cells wanted per superpixel
    sp_size = 4 #total number of cells
    compactness = 0.1 #how square/regular superpixels are
    avg_fun = 'median'
    
    #define number of superpixels there should be in each block to have desired superpixel size
    n_pixels_blocks = sapply(difference_rasts, ncell) #total pixels in each raster
    k_blocks = ceiling(n_pixels_blocks/sp_size)
    
    #generate superpixels for each block
    superpixels_l = pblapply(1:length(difference_rasts), function(i){
      r = difference_rasts[[i]]
      k = k_blocks[i] |> unname()
      sp = supercells(x=r, k=k, compactness=compactness, avg_fun = avg_fun, verbose =1)
      sp$block_id = block_ids[i]
      return(sp)
    })
    names(superpixels_l) = block_ids
    
    #plot
    {
      # # Select only columns of interest and geometry
      # super_sf <- superpixels_l$`12N_V1` %>%
      #   # select(GCC, NDGR, BI, GLI, geometry)
      #   select(-x,-y)
      # 
      # # Convert to long format
      # super_long <- pivot_longer(super_sf,
      #                            cols = -c(geometry, supercells, block_id),
      #                            names_to = "feature",
      #                            values_to = "value")
      # 
      # 
      # # Separate plots by feature
      # plots <- super_long %>%
      #   group_by(feature) %>%
      #   group_split() %>%
      #   map(~ {
      #     ggplot(.x) +
      #       geom_sf(aes(fill = value), color = NA) +
      #       scale_fill_viridis_c(option = "viridis") +
      #       ggtitle(unique(.x$feature)) +
      #       theme_minimal()
      #   })
      # 
      # # Combine with patchwork
      # wrap_plots(plots, nrow = 2)
    }
    
    #----binary classification----
    
    # #remove extraneous columns
    # classified_l = pblapply(superpixels_l, function(sp){
    #   #remove extraneous columns
    #   sp_data = sp |> st_drop_geometry() |> select(-c(x,y,supercells,block_id))
    #   
    #   #kmeans_clustering
    #   km = kmeans(sp_data, centers=2)
    #   
    #   #add clusters back to original dataframe
    #   sp$cluster_id = km$cluster
    #   sp$cluster_id = as.factor(sp$cluster_id)
    #   
    #   return(sp)
    # })
    # 
    # classified_combined = do.call(rbind, classified_l)
    # 
    # # Split your sf object by block_id
    # plots <- classified_combined %>%
    #   split(.$block_id) %>%
    #   lapply(function(df) {
    #     ggplot(df) +
    #       geom_sf(aes(fill = cluster_id), color=NA) +
    #       ggtitle(paste("Block", unique(df$block_id))) +
    #       theme_minimal()
    #   })
    # 
    # # Combine using patchwork
    # wrap_plots(plots)
    # 
    # 
    # 
    #----extract validation data----
    
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
    library(exactextractr)
    
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
      
      # #get mean canopy height or all chm pixels that fall within a superpixel, weighted by proportion of chm pixel covered by superpixel
      # tic()
      # cat('Extracting',b,'canopy height\n')
      # sp_validation = terra::extract(x = chm, y = sp, weights=T, fun = 'mean', bind=T, na.rm=T)
      # 
      # #get weighted mean canopy cover for each superpixel
      # cat('Extracting',b,'canopy cover\n')
      # sp_validation = terra::extract(x = cc, y = sp_validation, weights=T, fun = 'mean', bind=T, na.rm=T)
      # toc()
      
      #get weighted mean values for each superpixel
      tic()
      cat('Extracting',b,'validation data\n')
      extracted_vals = exact_extract(x = validation_lyrs, y = st_as_sf(sp), 'mean')
      sp_validation = cbind(sp, extracted_vals)
      cat(b,'finished\n')
      toc()
      
      return(sp_validation)
    })
    names(superpixels_validation_l) = block_ids
    
    # plan('sequential')
    
    #save sf object with validation data
    superpixels_validation_combined = do.call(rbind, lapply(superpixels_validation_l,st_as_sf))
    sp_validation_filename = paste0(bm_dir, '/superpixels_withvalidation_compactness=',compactness,'_size=',sp_size,'_fn=',avg_fun,'.geojson')
    if(!file.exists(sp_validation_filename)){st_write(superpixels_validation_combined, sp_validation_filename)}
    
    #----elastic net regression----
    # library(glmnet)
    # library(caret)
    # 
    # # Your model-fitting function
    # fit_caret_models <- function(sf_obj) {
    #   # Drop unwanted columns
    #   df <- sf_obj |>
    #     st_drop_geometry() |>
    #     subset(select = -c(supercells, x, y)) |>
    #     na.omit()
    #   
    #   # Ensure required columns exist
    #   if (!all(c("mean.canopy_height", "mean.canopy_cover") %in% names(df))) {
    #     stop("Missing canopy_height or canopy_cover columns.")
    #   }
    #   
    #   # Define predictors
    #   predictors <- setdiff(names(df), c("mean.canopy_height", "mean.canopy_cover"))
    #   
    #   # Set up cross-validation
    #   train_control <- trainControl(method = "cv", number = 5)
    #   
    #   # Tuning grid for alpha (lasso to ridge) and lambda (penalty)
    #   tune_grid <- expand.grid(
    #     alpha = seq(0, 1, length = 5),
    #     lambda = 10^seq(-4, 1, length = 10)
    #   )
    #   
    #   # Train canopy height model
    #   model_height <- train(
    #     x = df[, predictors],
    #     y = df$mean.canopy_height,
    #     method = "glmnet",
    #     trControl = train_control,
    #     tuneGrid = tune_grid,
    #     standardize = TRUE
    #   )
    #   
    #   # Train canopy cover model
    #   model_cover <- train(
    #     x = df[, predictors],
    #     y = df$mean.canopy_cover,
    #     method = "glmnet",
    #     trControl = train_control,
    #     tuneGrid = tune_grid,
    #     standardize = TRUE
    #   )
    #   cm_l = list(model_height, model_cover)
    #   names(cm_l) = c('canopy_height_model', 'canopy_cover_model')
    #   return(list(height_model = model_height, cover_model = model_cover))
    # }
    # 
    # ## Fit canopy height and canopy cover models
    # caret_models = lapply(block_ids, function(b){
    #   cat('Model fitting',b,'\n')
    #   sp_sf = superpixels_validation_combined |> filter(block_id==b)
    #   cm = fit_caret_models(sp_sf)
    #   return(cm)
    # })
    # names(caret_models) = names(superpixels_validation_l)
    # 
    # #extract info about the best models used for each block
    # # Function to extract model info from a single caret model
    # extract_model_info <- function(model, block_id, model_type) {
    #   # Best tuning parameters
    #   best <- model$bestTune
    #   
    #   # Find the corresponding row in the results
    #   res <- model$results %>%
    #     filter(alpha == best$alpha, lambda == best$lambda) %>%
    #     slice(1)  # In case of duplicates
    #   
    #   # Extract coefficients
    #   coefs <- as.matrix(coef(model$finalModel, s = best$lambda))
    #   coefs_df <- as.data.frame(t(coefs))  # Transpose so variables become columns
    #   names(coefs_df) <- rownames(coefs)
    #   
    #   # Combine into one row
    #   tibble(
    #     block_id = block_id,
    #     model_type = model_type,
    #     alpha = best$alpha,
    #     lambda = best$lambda,
    #     r_squared = res$Rsquared
    #   ) %>%
    #     bind_cols(coefs_df)
    # }
    # 
    # en_model_summary <- map_dfr(names(caret_models), function(block_id) {
    #   block_models <- caret_models[[block_id]]
    #   
    #   bind_rows(
    #     extract_model_info(block_models$height_model, block_id, "canopy_height"),
    #     extract_model_info(block_models$cover_model, block_id, "canopy_cover")
    #   )
    # })
    # 
    # x11()
    # ggplot(en_model_summary, aes(x = model_type, y = block_id...1))+
    #   geom_tile(aes(fill = r_squared))+
    #   geom_text(aes(label = round(r_squared,3)))
    
    # view(model_summary)
    
    #----GAMs with VSURF feature selection----
    
    library(mgcv)
    {
      # #vsurf feature selection for each block
      # library(VSURF)
      # 
      # set.seed(123)
      # vs_l = pblapply(block_ids, function(b){ #run vsurf
      #   sp_df = superpixels_validation_combined |>
      #     filter(block_id == b) |>
      #     st_drop_geometry() |>
      #     select(-c(supercells,x,y,block_id))|>
      #     na.omit()
      #   
      #   cat('VSURFing',b,'\n')
      #   #vsurf canopy cover
      #   cat('-----canopy cover\n')
      #   vs_cc = VSURF(x = sp_df |> select(-c(mean.canopy_cover, mean.canopy_height)),
      #                 y = sp_df$mean.canopy_cover,
      #                 RFimplem = 'ranger',
      #                 verbose = F,
      #                 parallel=T)
      #   #vsurf canopy cover
      #   cat('-----canopy height\n')
      #   vs_chm = VSURF(x = sp_df |> select(-c(mean.canopy_cover, mean.canopy_height)),
      #                 y = sp_df$mean.canopy_height,
      #                 RFimplem = 'ranger',
      #                 verbose = F,
      #                 parallel=T)
      #   #combine results
      #   vs_results = list(vs_cc, vs_chm)
      #   names(vs_results) = c('canopy_cover', 'canopy_height')
      #   return(vs_results)
      # })
    }
    
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
      cat('Modelling for',b,'\n')
      sp_df = superpixels_validation_combined |>
        filter(block_id == b) |>
        st_drop_geometry() |>
        select(-c(supercells,x,y,block_id))|>
        na.omit()
      cat('-----canopy cover\n')
      cc_results = gam_cv(data = sp_df|>select(-mean.canopy_height), target_var = 'mean.canopy_cover', kfolds = 5)
      cat('-----canopy height\n')
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
    
    cat('\n',d,'finished')
    return(gam_results_df)
  })
  
  gam_results_df = bind_rows(gam_results_df_l) |>
    mutate(dataset = basename(ps_directory))
  
  return(gam_results_df)
  }
  
  #plot GAM cross-validation results
  x11()
  ggplot(gam_results_df |> filter(fold!='overall'), aes(x = block_id, y = r2_test)) +
    geom_point() +
    facet_grid(rows = vars(target), cols=vars(dataset)) +
    theme_bw() +
    labs(
      x = "Block ID (site)", 
      y = "Out-of-sample R²", 
      title = "GAM 5-fold cross validation results for each site and target variable"
    ) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  #plot GAM results for each dataset
  x11()
  ggplot(gam_results_df |> 
           filter(fold=='overall')
         , aes(x = harvest_finish_date, y = r2_test))+
    geom_point(aes(shape = 'Harvest finish date'))+
    geom_text(aes(label = block_id), nudge_x = 35, nudge_y=-0.005, size = 3) +
    facet_grid(rows = vars(target), cols = vars(dataset)) +
    theme_bw()+
    xlab('Date')+
    ylab('Mean out-of-sample R²')+
    ggtitle('Generalized Additive Model results')+
    scale_x_date(limits = as.Date(c("2021-01-01", "2024-12-31"))) +
    # geom_vline(xintercept = as.Date(post_harvest_date), aes(xintercept = as.Date(post_harvest_date), linetype = 'Post-harvest date'))+
    # geom_vline(xintercept = as.Date(pre_harvest_date), aes(xintercept = as.Date(pre_harvest_date),linetype = 'Pre-harvest date'))
    # geom_vline(aes(xintercept = as.Date(paste0(str_replace_all(post_harvest_date,'_','-'),'-28')), linetype = 'Post-harvest image'), color = 'red')+
    # geom_vline(aes(xintercept = as.Date(paste0(str_replace_all(pre_harvest_date,'_','-'),'-28')),linetype = 'Pre-harvest image'),color ='blue')+
    scale_shape_manual(name = "Point Type", values = c('Harvest finish date' = 16)) +
    # scale_linetype_manual(name = "Line Type", 
    #                       values = c('Pre-harvest image' = 'dotted', 
    #                                  'Post-harvest image' = 'dashed'))+
    guides(
      shape = guide_legend(title = NULL)
      # ,linetype = guide_legend(title = NULL)
    )
  
  
}
