#___________Predicting forest thinning by applying superpixel clustering to PlanetScope data___________

#----load packages, set up environment----
source('basemap_setup_preprocessing.R')

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


#-----define wrapper function for processing data with GAM-----
superpixel_gam_modelling_fn = function(
    sp_size = 5
    # sp_step = 3
    , sp_compactness = 0.1
    , sp_avg_fn = 'mean'
){
  
  gam_results_df_l = pblapply(c(nonorm_dir, z_dir, zr_dir),function(d){
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
    
    #----calculate difference rasters----
    difference_rasts = mapply(function(x,y){x-y}, post_harvest_rasts, pre_harvest_rasts)
    
    #----subset each raster to only include desired bands----
    difference_rasts = lapply(difference_rasts, function(r)r[[c('blue', 'green', 'red','Hue', 'GCC', 'NDGR', 'BI', 'CI', 'CRI550', 'GLI')]])
    
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

#-----define wrapper function for processing data with random forest-----
superpixel_rf_modelling_fn = function(
    sp_size = 5
    # sp_step = 3
    , sp_compactness = 0.1
    , sp_avg_fn = 'mean'
){
  
  model_results_df_l = pblapply(c(nonorm_dir, z_dir, zr_dir),function(d){
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
    
    #----calculate difference rasters----
    difference_rasts = mapply(function(x,y){x-y}, post_harvest_rasts, pre_harvest_rasts)
    
    #----subset each raster to only include desired bands----
    difference_rasts = lapply(difference_rasts, function(r)r[[c('blue', 'green', 'red','Hue', 'GCC', 'NDGR', 'BI', 'CI', 'CRI550', 'GLI')]])
    
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
    cat('\n-----Fitting random forest models')
    if(!'package:mgcv' %in% search()){library(mgcv)}
    
    
    
    #write function to do k-fold cross validation with generalized additive models
    model_cv = function(data, target_var, predictors = colnames(data)[colnames(data)!=target_var], kfolds=5, ...){
      
      #partition data into k folds
      folds = createFolds(y = data[[target_var]], k = kfolds)
      
      #generate list of models using each fold
      model_l = pblapply(folds, function(fold){
        
        #generate model
        data_train = data[-fold,]
        smooth_terms = smooth_terms <- paste0("s(", predictors, ")")
        form = as.formula(
          paste(target_var,'~',paste(smooth_terms, collapse = ' + ')))
        
        #fit model on training data
        model_mod = bam(form, data=data_train, method='REML',...)
        
        #get model performance metrics based on training data
        s_train = summary(model_mod) #r2
        r2_train = s_train$r.sq
        train_predictions = model_mod$fitted.values
        train_actual = data_train[[target_var]]
        bias_train = mean(train_actual-train_predictions)#bias
        rmse_train = sqrt(mean((train_actual-train_predictions)^2)) #rmse
        nrmse_train = rmse_train/IQR(train_actual) #nrmse
        
        #get r2 of model applied to testing data
        data_test = data[fold,]
        test_predictions = predict.bam(model_mod, newdata = data_test) |> as.numeric()
        test_actual = data_test[[target_var]]
        
        ss_res <- sum((test_actual - test_predictions)^2) 
        ss_tot <- sum((test_actual - mean(test_actual))^2)
        r2_test <- 1 - ss_res / ss_tot
        
        bias_test = mean(test_actual-test_predictions)#bias
        rmse_test = sqrt(mean((test_actual-test_predictions)^2)) #rmse
        nrmse_test = rmse_test/IQR(test_actual) #nrmse
        
        
        #return results and model performance stats
        gm_l = list(model_mod, r2_test, rmse_test, nrmse_test, bias_test, r2_train, rmse_train, nrmse_train, bias_train)
        names(gm_l) =  c('model.model', 'r.sqr_test', 'rmse_test', 'nrmse_test', 'bias_test', 'r.sqr_train', 'rmse_train', 'nrmse_train', 'bias_train')
        
        return(gm_l)
      })
      
      metrics <- c('r.sqr_test', 'rmse_test', 'nrmse_test', 'bias_test', 
                   'r.sqr_train', 'rmse_train', 'nrmse_train', 'bias_train')
      mean_metrics <- setNames(
        lapply(metrics, function(metric) {
          mean(sapply(model_l, function(x) x[[metric]]))
        }),
        metrics
      )
      
      cv_results <- c(mean_metrics, list(kfold_results = model_l))
      
      return(cv_results)
    }
    
    #get random forest models and model performance for canopy height and canopy cover for each block
    set.seed(123)
    model_results_l = pblapply(block_ids, function(b){
      # cat('Modelling for',b,'\n')
      sp_df = superpixels_validation_combined |>
        filter(block_id == b) |>
        st_drop_geometry() |>
        select(-c(supercells,x,y,block_id))|>
        na.omit()
      # cat('-----canopy cover\n')
      cc_results = model_cv(data = sp_df|>select(-mean.canopy_height), target_var = 'mean.canopy_cover', kfolds = 5)
      # cat('-----canopy height\n')
      ch_results = model_cv(data = sp_df |>select(-mean.canopy_cover), target_var = 'mean.canopy_height', kfolds = 5)
      
      results_l = list(cc_results, ch_results)
      names(results_l) = c('canopy_cover', 'canopy_height')
      return(results_l)
    })
    names(model_results_l) = block_ids
    
    #extract results from models
    extract_model_metrics <- function(model_results_l) {
      results_list <- list()
      
      for (site_id in names(model_results_l)) {
        site_data <- model_results_l[[site_id]]
        
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
    
    model_results_df = extract_model_metrics(model_results_l) |> 
      rename(block_id=site) |>
      left_join(harvest_dates_df) |>
      mutate(ps_directory = ps_directory)
    
    cat('\n',basename(d),'finished')
    return(model_results_df)
  })
  
  #---combine list of model results into dataframe----
  model_results_df = bind_rows(model_results_df_l) |>
    mutate(dataset = basename(ps_directory)) |>
    mutate(dataset = case_when(
      str_detect(basename(ps_directory),'Zrobust') ~ 'Zrobust',
      str_detect(basename(ps_directory),'Z') ~ 'Z',
      .default = 'Non-normalized'
    ))
  
  return(model_results_df)
  
}


#-----tune superpixel hyperparameters with GAM-----
sp_params = expand.grid(
  # sp_size = c(1,2,3,4,5,9,16,32),
  sp_size = seq(1,15,by=1),
  compactness = c(0.1,0.5,1,2,seq(5,10,by=0.5)),
  avg_fn = c('mean'
             # , 'median'
             )
) |>
# |>
#   filter(
#     (sp_size != 1) |
#       (sp_size == 1 & compactness == 10 & avg_fn == "mean")
#   ) |>
  mutate(avg_fn = as.character(avg_fn))

sp_opt_dir = paste0(bm_dir, '/superpixel_optimization_results_gam2')
dir.check(sp_opt_dir)

cl = makeCluster(16)
plan(future::cluster, workers = cl)
future_lapply(1:nrow(sp_params), function(i){
  cat('\nTesting combination',i,'out of',nrow(sp_params))
  sp_size = sp_params$sp_size[i]
  compactness = sp_params$compactness[i]
  avg_fn = sp_params$avg_fn[i]
  
  filename = paste0(sp_opt_dir,'/spsize=',sp_size,'_compactness=',compactness,'_fn=',avg_fn,'.arrow')
  if(!file.exists(filename)){
    
    tic()
    sp_mod_results = superpixel_gam_modelling_fn(sp_size=sp_size
                                             # sp_step = sp_step
                                             , sp_compactness = compactness
                                             , sp_avg_fn = avg_fn)
    # sp_mod_results$sp_step = sp_step
    sp_mod_results$sp_size = sp_size
    sp_mod_results$sp_compactness = compactness
    sp_mod_results$sp_avg_fn = avg_fn
    toc()
    write_feather(sp_mod_results, filename)
    return(sp_mod_results)
  }
})
stopCluster(cl)
plan('sequential')



