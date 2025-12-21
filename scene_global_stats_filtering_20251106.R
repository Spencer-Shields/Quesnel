# source('scene_setup_preprocessing_20250909.R')
source('derive lidar metrics_202501020.R')
library(varSel)
# library(pROC)
library(arrow)
library(data.table)

#----load validation data----

#load Otsu thresholds
lid_met = 'z_p95' #define lidar metric to use for generating binary thinning layers

otsu_f = paste0("data/Quesnel_thinning/Otsu_change_thresholds_Res=",target_res,'lidmet=',lid_met,".csv")
if(!file.exists(otsu_f)){
  lcmfl = list.files(lid_mets_change_cropped_dir, full.names = T)
  otsu_df = pblapply(block_ids, function(b){
    f = lcmfl[str_detect(lcmfl,b)]
    r = rast(f, lyrs=lid_met)
    o = values(r, na.rm=T) |>
      otsuThresholdCpp(bins=256)
    tbl = tibble(
      block_id = b,
      chm_change_otsu_threshold = o)
    return(tbl)
  }) |> bind_rows()
  write_csv(otsu_df,otsu_f)
} else{
  otsu_df = read_csv(otsu_f)
}

#load LiDAR data
# lidar_dir = paste0('data/Quesnel_thinning/chm_change_scenes_Res=',target_res)
# lidar_files = list.files(lidar_dir, full.names = T, recursive = T)
lidar_files = list.files(lid_mets_change_cropped_dir,full.names = T)
lidar_ids = basename(file_path_sans_ext(lidar_files))
names(lidar_files) = lidar_ids



#----extract stats from each raster and save----
harvest_threshold = paste0("Otsu_3.0m_",lid_met)

data_string = paste0('Stats_percentiles_logAUC_JMD_FDR_HarvestThreshold=',harvest_threshold,'_Res=',target_res,'_EqualClasses_NoQualChecks_20250910')
data_filename = paste0(ps_dir,'/',data_string,'.arrow')

#run if summary data file does not exist
print('Checking if summary stats file exists')
if(!file.exists(data_filename)){
  
  #----set up for processing----
  
  all_files <- unlist(sapply(dirs, list.files, pattern = "\\.tif$", full.names = TRUE, recursive = TRUE))
  all_files = all_files[!str_detect(all_files, '_NULL')] #remove dummy files
  
  #get vector of all possible dates in timeseries
  start_date <- as.Date("2020-01-01")
  end_date <- as.Date("2025-01-01")
  date_vector <- seq.Date(from = start_date, to = end_date, by = "month")
  date_vector <- format(date_vector, "%Y-%m-%d")
  date_vector = as.character(date_vector)
  date_vector = str_replace_all(date_vector,'-','_')
  
  #all possible datasets (strings to look for in filepaths to tell what dataset a raster belongs to)
  dataset_strings = paste0(basename(dirs),'/') #include / to prevent confusion between Z and Zrobust
  
  #block ids
  block_ids = thinning_block_ids
  
  #remove NoChange rasters from list of all files
  all_files_thinning = all_files[!str_detect(all_files, 'NoChange')]
  
  # all_files_thinning = rev(all_files_thinning) #reverse order of list (for processing in a separate session)
  
  #----create thinning masks, get spatial sample for each study area to make thinning/non-thinning classes same size----
  
  #load lidar CHM change layers
  lid_masks_equalclasses = pblapply(lidar_files, function(x){
    
    ##load lidar raster
    # print(paste('Processing',x))
    r = rast(x, lyrs = lid_met)
    
    #get appropriate harvest threshold
    if(str_detect(harvest_threshold, 'Otsu')){
      lid_id = basename(file_path_sans_ext(x))
      ht = otsu_df$chm_change_otsu_threshold[otsu_df$block_id == lid_id]
    } else {
      ht = harvest_threshold
    }
    
    #create raster stack with thinning/nonthinning mask layers
    m_t = ifel(r <= ht, 1, NA) #create thinning mask layer
    m_nt = ifel(r > ht, 0, NA) #create non-thinning mask layer
    m = c(m_t, m_nt) #stack layers
    names(m) = c('Thinning', 'Non_thinning') #name layers
    
    #get number of non-NA pixels in each layer of the raster stack
    gvals = global(m, 'notNA')
    
    n_t = gvals$notNA[rownames(gvals)=='Thinning'] #number of non-NA thinning pixels
    n_nt = gvals$notNA[rownames(gvals)=='Non_thinning'] #number of non-NA non-thinning pixels
    
    #if one layer has more pixels than the other, do a subsample so that they have the same number and replace the raster layer with the subsampled layer
    if(n_t > n_nt){
      target_size = n_nt
      #sample thinning mask
      m_s_xy = terra::spatSample(m_t, size = target_size, xy=T, method='random', na.rm=T, replace=F, exhaustive = T)
      m_s = rast(m_s_xy, type = 'xyz')
      m[[1]] = m_s |> resample(m_t)
    }
    if(n_t < n_nt){
      target_size = n_t
      #sample thinning mask
      m_s_xy = terra::spatSample(m_nt, size = target_size, xy=T, method='random', na.rm=T, replace=F, exhaustive = T)
      m_s = rast(m_s_xy, type = 'xyz')
      # m_s = ifel(m_s == 1, 0,NA) #reset values to zero since they become one for some reason
      m[[2]] = m_s |> resample(m_nt)
    }
    
    names(m) = c('Thinning', 'Non-thinning')
    m = wrap(m) #wrap for parallel processing
    return(m)
  })
  
  #----mask PS rasters, save global values in dataframe, calculate stats----
  
  #parallelize
  clust = makeCluster(8)
  plan('cluster', workers = clust)
  
  #make directory to store stats tables for each scene
  global_tables_dir = paste0(ps_dir,'/',data_string)
  dir.check(global_tables_dir)
  
  # df_l = pblapply(1:10, function(i){
  # pblapply(rev(1:length(all_files_thinning)), function(i){
  future_lapply(rev(1:length(all_files_thinning)), function(i){
    
    x = all_files_thinning[i]
    print(paste0('Processing ', x,', ',i,'/',length(all_files_thinning))) #uncomment to identify files where processing fails
    
    id = basename(file_path_sans_ext(x))
    dataset = find_substring(x,dataset_strings) |> str_replace('/','')
    block = find_substring(x, lidar_ids)
    
    
    filename = paste0(global_tables_dir,'/',dataset,'_',block,'_',id,'.arrow')
    if(!file.exists(filename)){
      
      #get raster
      r = rast(x)
      
      # #drop extra TVI if it's there
      # if(names(r)[35] == 'TVI'){r = r[[-35]]}
      
      #get lidar mask
      block = find_substring(x, lidar_ids)
      lid_mask_wrapped = lid_masks_equalclasses[[block]]
      m_ = unwrap(lid_mask_wrapped)
      m = resample(m_,r)
      
      #mask planetscope raster using thinning and non-thinning masks
      r_t = mask(r, m[[1]])
      r_nt = mask(r, m[[2]])
      
      #calculate global stats for thinning pixels
      d1 = global(r_t, fun = c("max", "min", "mean", "sum", "range", "rms", "sd", "std", "isNA", "notNA"), na.rm=T)
      d1 = cbind(d1, 
                 global(r_t, median, na.rm=T) |> rename(median = global)
                 ,rast_quantiles(r_t, seq(0.05,0.95,0.05))
      )
      d1[['block_pixel_stratum']] = 'Thinned'
      d1[['v1']] = rownames(d1)
      
      #calculate global stats for non-thinning pixels
      d2 = global(r_nt, fun = c("max", "min", "mean", "sum", "range", "rms", "sd", "std", "isNA", "notNA"), na.rm=T)
      d2 = cbind(d2, 
                 global(r_nt, median, na.rm=T) |> rename(median = global)
                 ,rast_quantiles(r_nt, seq(0.05,0.95,0.05))
      )
      d2[['block_pixel_stratum']] = 'Not_thinned'
      d2[['v1']] = rownames(d2)
      
      #total block stats
      d3 = global(r, fun = c("max", "min", "mean", "sum", "range", "rms", "sd", "std", "isNA", "notNA"), na.rm=T)
      d3 = cbind(d3, 
                 global(r, median, na.rm=T) |> rename(median = global)
                 ,rast_quantiles(r, seq(0.05,0.95,0.05))
      )
      d3[['block_pixel_stratum']] = 'Total'
      d3[['v1']] = rownames(d3)
      
      #combine thinning, non-thinning, and global stats
      d = rbind(d1,d2,d3)
      
      #get dataset
      d[['dataset']] = str_remove(dataset,'/')
      
      #get block
      d[['block_id']] = block
      
      #get date
      date_ = find_substring(x, date_vector)
      d[['acquisition_date']] = date_
      
      #filepath
      d[['file_path']] = x
      
      ###statistical measures between groups
      
      v_t = values(r_t, dataframe=T) |> mutate(thinning = 'Thinned')
      v_nt = values(r_nt, dataframe=T) |> mutate(thinning = 'Not_thinned') 
      v_df = rbind(v_t, v_nt)
      
      feats = unique(d$v1)
      feats = feats[!feats %in% c('max_DN'
                                  # , 'VARIgreen'
      )] #don't run tests on max_DN
      
      tests_l = pblapply(feats, function(y){
        
        print(paste('Processing',y))
        #subset values to test from dataframe
        v = v_df |> select(all_of(c(y, 'thinning')))
        v = v[!is.na(v[[y]]),] #remove NA values
        v$thinning = as.factor(v$thinning)
        
        #initialize tibble to store results
        test_df = tibble(
          v1 = y,
          logistic_p = NA,
          logistic_aic = NA,
          logistic_AUC = NA,
          fisher_discriminant_ratio = NA,
          cohen_d = NA,
          jeffries_matusita_dist = NA
          # ,
          # bhattacharya_dist = NA
        )
        
        # #get proportion not na
        # proportion_na = d1$notNA[d1$v1==y]/d1$isNA[d1$v1==y]
        # 
        # if(proportion_na > 0.05 & length(unique(v[[y]])) > 10){ #proceed if proportion na pixels is greater than threshold and there are multiple cell values
        
        ##perform tests, store results in test_df
        
        # #initialize tibble to store results
        # test_df = tibble(
        #   v1 = y,
        #   logistic_p = NA,
        #   logistic_aic = NA,
        #   logistic_AUC = NA,
        #   fisher_discriminant_ratio = NA,
        #   cohen_d = NA,
        #   jeffries_matusita_dist = NA
        #   # ,
        #   # bhattacharya_dist = NA
        # )
        
        # #get proportion not na
        # proportion_na = d1$notNA[d1$v1==y]/d1$isNA[d1$v1==y]
        # 
        # if(proportion_na > 0.05 & length(unique(v[[y]])) > 10){ #proceed if proportion na pixels is greater than threshold and there are multiple unique cell values
        ##perform tests, store results in test_df
        
        #logistic regression model with tryCatch
        log_form = paste('thinning ~', y)
        
        tryCatch({
          log_m = glm(log_form, v, family = 'binomial')
          log_m_summ = summary(log_m)
          
          test_df$logistic_p = ifelse(is.na(log_m$coefficients[[y]]), 
                                      NA,
                                      log_m_summ$coefficients[,4][[y]])
          test_df$logistic_aic = log_m_summ$aic
          
          log_p = v |> mutate(predictions = predict(log_m, v, type = 'response'))
          log_roc = pROC::roc(thinning ~ predictions, log_p, quiet=T)
          test_df$logistic_AUC = as.numeric(auc(log_roc))
        }, 
        error = function(e) {
          message("GLM error encountered: ", e$message)
          # All logistic regression values will remain NA as initialized
        })
        
        #Fisher discriminant ratio
        fdr = ((d1$mean[d1$v1==y]-d2$mean[d2$v1==y])^2)/(d1$sd[d1$v1==y]^2 + d2$sd[d2$v1==y]^2)
        test_df$fisher_discriminant_ratio = fdr
        
        #Cohen's d
        pooled_sd = sqrt(((d1$notNA[d1$v1 == y] - 1) * d1$sd[d1$v1 == y]^2 +
                            (d2$notNA[d2$v1 == y] - 1) * d2$sd[d2$v1 == y]^2) /
                           (d1$notNA[d1$v1 == y] + d2$notNA[d2$v1 == y] - 2))
        cohen_d = (d1$mean[d1$v1==y] - d2$mean[d2$v1==y])/pooled_sd
        test_df$cohen_d = cohen_d
        
        #Jeffries-Matusita distance
        tryCatch({
          jmd = varSel::JMdist(g = v$thinning, X = tibble(v[[y]]))
          test_df$jeffries_matusita_dist = jmd$jmdist
        },
        error = function(e) {
          message("Jeffries-Matusita distance error: ", e$message)
          # jeffries_matusita_dist will remain NA as initialized
        })
        
        # #Bhattacharya distance
        # tryCatch({
        #   bhat = BHATdist(g = v$thinning, X = tibble(v[[y]]))
        #   test_df$bhattacharya_dist = bhat$bhatdist
        # },
        # error = function(e) {
        #   message("Bhattacharya distance error: ", e$message)
        #   # bhattacharya_dist will remain NA as initialized
        # })
        
        # }
        
        
        # }
        
        #return results dataframe
        return(test_df)
      })
      test_results = bind_rows(tests_l)
      
      ###attach statistical test results to the summary stats dataframe
      
      d = d %>% left_join(test_results, by = 'v1')
      
      ####save the final result
      write_feather(d, filename)
      # return(d)
    }
  }
  # ,cl = 'future'
  # ,future.seed=T
  )
  
  # results_df1 = bind_rows(df_l)
  # results_df = setDT(results_df1)
  # 
  # fwrite(results_df, data_filename)
  global_tables_files = list.files(global_tables_dir, full.names = T, pattern = '\\.arrow$')
  # tables_l = pblapply(all_global_tables, function(x)setDT(read_feather(x)))
  
  #load tables into a list, save indices of tables that do not load properly
  global_tables_l = list()
  bad_indices = c()
  
  cat('Double checking for corrupted ARROW files')
  for (i in seq_along(global_tables_files)) {
    x = global_tables_files[i]
    # print(paste('Reading', x,', ',i,'out of',length(global_tables_files)))
    result <- tryCatch({
      df <- read_feather(x) |> setDT()
      global_tables_l[[length(global_tables_l) + 1]] <- df  # Add to list
      NULL  # No error
    }, error = function(e) {
      message(sprintf("Failed at index %d: %s", i, e$message))
      i  # Return the index that failed
    })
    if (!is.null(result)) {
      bad_indices <- c(bad_indices, result)
    }
    if(i %% 1000 == 0){ #keep track of progress
      cat('\rFinished',i,'of',length(global_tables_files))
    }
  }
  cat('Bad indices:',bad_indices)
  
  # #delete bad files, manually go back and recreate
  # lapply(bad_indices, function(i){file.remove(global_tables_files[i])})
  
  results_df1 = rbindlist(global_tables_l)
  results_df = results_df1 |>
    #get scene ID
    mutate(id = basename(file_path)) |>
    #fix acquisition_date
    mutate(acquisition_date = as.Date(paste0(
      substr(id,1,4),'-', #year
      substr(id,5,6),'-', #month
      substr(id,7,8) #day
    ))) |>
    rename(var_name = v1)
  
  write_feather(results_df,data_filename)
  
  plan('sequential')
  stopCluster(clust)
  
  
}


#----load stats table, clean and process dataframe----
{
  print('Loading and cleaning summary stats')
  tic()
  results_df_load = read_feather(data_filename) |> setDT()
  
  results_df_ = results_df_load |>
    # #remove unwanted indices
    # filter(!var_name %in% c('WDVI','CI','NDWI','RVI','TTVI', 'VARIgreen')) |> #indices to remove
    #only use NDVI
    filter(var_name=='NDVI') |>
    #add harvest dates
    left_join(harvest_dates_df, by = 'block_id') |>
    #define variable for time before and after harvest_start_date
    mutate(time_since_harvest_start = as.numeric(acquisition_date - harvest_start_date),
           time_since_harvest_finish = as.numeric(acquisition_date - harvest_finish_date)) |>
    #make harvested dummy variable
    mutate(harvested = ifelse(acquisition_date < harvest_start_date,0,1)) |>
    left_join(tibble(acquisition_date = sort(unique(results_df_load$acquisition_date)),
                     acquisition_date_ordinal = as.numeric(
                       1:length(unique(results_df_load$acquisition_date))))) |>
    #make other dummy variables for more refined segmentation of time series
    mutate(mid_harvest = ifelse(acquisition_date >= harvest_start_date & acquisition_date < harvest_finish_date, 1,0)
           ,post_harvest = ifelse(acquisition_date >= harvest_finish_date, 1, 0)
    )|>
    #add acquisition month
    mutate(month = month(acquisition_date)) |>
    #fix the id column since it currently has file paths
    mutate(id = file_path_sans_ext(id))|>
    # #rescale JM distance to give it the range [0-1] so that beta regression can be applied to both AUC and JMD
    mutate(jmd_scaled = jeffries_matusita_dist/sqrt(2)) |>
    #reformat datasets
    mutate(dataset = case_when(
      str_detect(dataset, 'SM') ~ 'SM',
      str_detect(dataset, 'Non-normalized') ~ 'Non-normalized', 
      str_detect(dataset, 'Zrobust') ~ 'Zrobust',
      str_detect(dataset, 'Z') ~ 'Z',
      TRUE ~ dataset  # keep original value if no match
    ))
  toc()
}

#remove scenes that should have been removed during previous steps
results_df = results_df_ |> na.omit()

#----join metadata table with results----

meta_df = setDT(meta_df) |> mutate(acquisition_date = as.Date(acquisition_date))
results_meta_df = results_df |> left_join(meta_df, by=c('id','acquisition_date'))


#----filter images based on percent valid pixels----

#get maximum number of valid pixels for each combination of block_id, variable, stratum (thinning and nonthinng), and dataset
max_valid_pixels_byfeat_byblock_bypixstrat = results_meta_df |>
  group_by(block_id,var_name,block_pixel_stratum,dataset)|>
  summarise(max_valid_pixels = max(notNA))

#join maximum number of valid pixels to the results dataframe, get the proportion of valid pixels for each row
results_meta_df = results_meta_df |> 
  left_join(max_valid_pixels_byfeat_byblock_bypixstrat) |>
  mutate(valid_proportion = notNA/max_valid_pixels) |>
  mutate(valid_proportion_rounded = round(valid_proportion, 2))


#extract quality thresholds for each block, attach to results dataframe
qual_threshs = sapply(block_ids, function(b){
  d = results_meta_df |> filter(block_pixel_stratum=='Total', block_id==b)
  h = hist(d$valid_proportion_rounded, breaks = seq(0, 1, 0.01), plot=F)
  tib = tibble(min_break = h$breaks[h$breaks!=1], count = h$counts)
  thresh = tib$min_break[tib$count==max(tib$count)]
  return(thresh)
})

qual_thresh_df = tibble(
  valid_proportion_threshold = qual_threshs, block_id = block_ids
)

results_meta_df = results_meta_df |> left_join(qual_thresh_df)


#view histograms for proportion valid pixels
{
  hdf = results_meta_df |> select(dataset,block_id,id,valid_proportion_rounded,var_name,valid_proportion_threshold) |> distinct()
  
  ggplot(hdf)+
    geom_histogram(aes(valid_proportion_rounded),bins=99)+
    facet_grid(rows=vars(block_id), cols=vars(dataset))
  
  hdf_summ = hdf |> group_by(dataset,block_id,var_name)|>summarise(n = n()) |> pivot_wider(names_from = dataset, values_from = n)
  hdf_summ
}


#get ids of scenes where the quality threshold is met for each feature
hq_scenes_df_ = results_meta_df |>
  filter(block_pixel_stratum=='Total',
         valid_proportion_rounded >= valid_proportion_threshold) |>
  pivot_wider(id_cols = c('id','dataset','block_id','acquisition_date')
              , names_from = var_name, values_from = valid_proportion)
hq_scenes_df = hq_scenes_df_ |> na.omit()

#filter out scenes where quality threshold is not met for each feature
results_meta_df = results_meta_df |>
  semi_join(hq_scenes_df |> select(id, dataset, block_id, acquisition_date),
            by = c("id", "dataset", "block_id", "acquisition_date"))



message('Final dataframe with valid scenes: results_meta_df')

#evaluate which scenes are missing (NOTE: NOT AN ISSUE WHEN ONLY NDVI IS ANALYZED)
{
  # #tally scenes per block and dataset (RESULT: softmax is consistently missing pixels)
  # hq_scenes_byblock_sum = results_meta_df |>
  #   filter(block_pixel_stratum=='Total',
  #          var_name=='NDVI') |>
  #   group_by(block_id,dataset) |>
  #   summarise(n = n()) |>
  #   pivot_wider(names_from = 'dataset', values_from = 'n')
  # 
  # #find ids of images which are missing in softmax
  # a = results_meta_df |>
  #   select(block_id,dataset,id) |>
  #   distinct() |>
  #   pivot_wider(names_from = dataset, values_from = id)
  # 
  # missing_ids_df_ = lapply(block_ids, function(b){
  #   a_ = a|>filter(block_id==b)
  #   missing_ids = setdiff(unlist(a_$`Non-normalized`), unlist(a_$SM)) #|> unlist()
  #   t = tibble(
  #     id = missing_ids,
  #     block_id = b
  #   )
  #   return(t)
  # }) |> bind_rows()
  # 
  # fl = list.files(dirs, recursive = T, full.names = T)
  # fl_sm = fl[str_detect(fl, 'SM')]
  # fl_sm_missing = c()
  # 
  # for(i in 1:nrow(missing_ids_df_)){
  #   
  #   block_id = missing_ids_df_$block_id[i]
  #   id = missing_ids_df_$id[i]
  #   
  #   f = fl_sm[which(str_detect(fl_sm,id) & str_detect(fl_sm, block_id))]
  #   
  #   fl_sm_missing = c(fl_sm_missing, f)
  # }
  # 
}

# #tally how many dates there are in each month and day for each block
# hq_scenes_permonth_summ = results_meta_df |>
#   mutate(year_mon = format(acquisition_date, "%Y-%m")) |>
#   group_by(dataset, block_id, year_mon) |>
#   summarise(count = n()) |>
#   mutate(acquisition_month = as.numeric(str_remove_all(year_mon, '.*-')))
# 
# jan_june_labels <- hq_scenes_permonth_summ$year_mon[hq_scenes_permonth_summ$acquisition_month %in% c(1, 6)]
# 
# ggplot(hq_scenes_permonth_summ |> filter(str_detect(dataset, 'SM'))
#        , aes(x = year_mon, y = count))+
#   geom_point()+
#   # geom_line()+
#   facet_grid(rows=vars(block_id)
#              # , cols = vars(block_id)
#   )+
#   scale_x_discrete(breaks = jan_june_labels)


#get ids which are common to each thinning block
# ds_ids = sapply(unique(results_meta_df$dataset), function(d)results_meta_df$id[results_meta_df$dataset==d])
# common_ids = Reduce(intersect, ds_ids)

common_ids <- results_meta_df %>%
  group_by(block_id, dataset) %>%
  summarise(ids = list(unique(id)), .groups = "drop") %>%
  summarise(common_ids = list(Reduce(intersect, ids)), .by = block_id)
ids_in_all_blocks <- Reduce(intersect, common_ids$common_ids)
