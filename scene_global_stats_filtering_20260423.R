# source('scene_setup_preprocessing_20250909.R')
source('derive lidar metrics_20260423.R')
# library(varSel)
# library(pROC)
library(arrow)
library(data.table)
library(sgsR)

#----load lidar data----

#load validation data
lid_met = 'z_above2' #define lidar metric to use for thematic thinning layers

#load LiDAR data
lidar_files = list.files(lid_mets_change_cropped_dir,full.names = T)
lidar_ids = basename(file_path_sans_ext(lidar_files))
names(lidar_files) = lidar_ids

# sbreaks = c(-0.90, -0.50, -0.25, -0.05)

#----set breaks for cc change bins----

##manual breaks
# sbreaks = c(-0.90, -0.50, -0.25, -0.05)

#jenks breaks
combined_lid_df = lapply(lidar_files, function(f){
  lr = rast(f, lyrs = lid_met)
  ld = as.data.frame(lr)
}) |> bind_rows()

b = classIntervals(combined_lid_df[[lid_met]], 4, style = 'jenks')

sbreaks = b$brks[2:(length(b$brks)-1)] |> round(2)

#----extract stats from each raster and save----
sbreaks_str = paste0(sbreaks,collapse = ' ')

#planetscope layers to include
included_lyrs = c('NDVI') |>
  paste0(collapse = ' ')

data_string = paste0('Stats_percentiles_logAUC_JMD_FDR_ClassBreaks=',sbreaks_str,'_Res=',target_res,'_IncludedLyrs=',included_lyrs)
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
  
  #make lidar canopy cover change masks
  lid_masks = pblapply(lidar_files, function(x){
    
    ##load lidar raster
    # print(paste('Processing',x))
    r = rast(x, lyrs = lid_met)
    
    #stratify raster
    rs = strat_breaks(r, sbreaks, map=T)
    
    #make each strata its own raster layer
    classes = unique(rs)[[1]]
    
    m = lapply(classes, function(i){
      r_ = ifel(rs==i,1,NA)
    }) |> rast()
    names(m) = classes
    
    return(m)
  })
  
  lid_masks = lapply(lid_masks, wrap)
  
  #----mask PS rasters, save global values in dataframe, calculate stats----
  
  # parallelize
  library(futurize)
  library(progressify)
  
  plan(multisession, workers = 8)
  
  on.exit({
    plan(sequential)
  }, add = TRUE)
  
  # clust = makeCluster(8)
  # plan('cluster', workers = clust)
  
  #make directory to store stats tables for each scene
  global_tables_dir = paste0(ps_dir,'/',data_string)
  dir.check(global_tables_dir)
  
  # df_l = pblapply(1:10, function(i){
  # pblapply(rev(1:length(all_files_thinning)), function(i){
  lapply(rev(1:length(all_files_thinning)), function(i){
    
    x = all_files_thinning[i]
    print(paste0('Processing ', x,', ',i,'/',length(all_files_thinning))) #uncomment to identify files where processing fails
    
    id = basename(file_path_sans_ext(x))
    dataset = find_substring(x,dataset_strings) |> str_replace('/','')
    block = find_substring(x, lidar_ids)
    
    
    filename = paste0(global_tables_dir,'/',dataset,'_',block,'_',id,'.arrow')
    if(!file.exists(filename)){
      
      #get raster
      r = rast(x)
      
      if(!included_lyrs %in% names(r)){
        cat('Invalid file, skipping',x)
        return(NULL)
      } else {
        r = r[[included_lyrs]]
      }
      # #drop extra TVI if it's there
      # if(names(r)[35] == 'TVI'){r = r[[-35]]}
      
      #get lidar mask
      block = find_substring(x, lidar_ids)
      lid_mask_wrapped = lid_masks[[block]]
      m_ = unwrap(lid_mask_wrapped)
      m = resample(m_,r)
      
      #mask planetscope raster using thinning and non-thinning masks
      r_t = mask(r, m[[1]])
      r_nt = mask(r, m[[2]])
      
      dx = lapply(1:nlyr(m), function(j){
        mj = m[[j]]
        r_t = mask(r,mj)
        
        #calculate global stats for thinning pixels
        d1 = global(r_t, fun = c("max", "min", "mean", "sum", "range", "rms", "sd", "std", "isNA", "notNA"), na.rm=T)
        d1 = cbind(d1, 
                   global(r_t, median, na.rm=T) |> rename(median = global)
                   ,rast_quantiles(r_t, seq(0.05,0.95,0.05))
        )
        d1[['block_pixel_stratum']] = j
        d1[['v1']] = rownames(d1)
        
        return(d1)
      }) |> bind_rows()
      
      #total block stats
      d3 = global(r, fun = c("max", "min", "mean", "sum", "range", "rms", "sd", "std", "isNA", "notNA"), na.rm=T)
      d3 = cbind(d3, 
                 global(r, median, na.rm=T) |> rename(median = global)
                 ,rast_quantiles(r, seq(0.05,0.95,0.05))
      )
      d3[['block_pixel_stratum']] = 'Total'
      d3[['v1']] = rownames(d3)
      
      #combine stats per stratum and global stats
      d = rbind(dx,d3)
      
      #get dataset
      d[['dataset']] = str_remove(dataset,'/')
      
      #get block
      d[['block_id']] = block
      
      #get date
      date_ = find_substring(x, date_vector)
      d[['acquisition_date']] = date_
      
      #filepath
      d[['file_path']] = x
      
      ####save the final result
      write_feather(d, filename)
    }
  }
  ) |> 
    progressify() |> 
    futurize() 
  
  plan(sequential)
  
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
}


#----load stats table, clean and process dataframe----
{
  print('Loading and cleaning summary stats')
  
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
    #reformat datasets
    mutate(dataset = case_when(
      str_detect(dataset, 'SM') ~ 'SM',
      str_detect(dataset, 'Non-normalized') ~ 'Non-normalized', 
      str_detect(dataset, 'Zrobust') ~ 'Zrobust',
      str_detect(dataset, 'Z') ~ 'Z',
      TRUE ~ dataset  # keep original value if no match
    )) 
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


#get ids which are common to each thinning block
# ds_ids = sapply(unique(results_meta_df$dataset), function(d)results_meta_df$id[results_meta_df$dataset==d])
# common_ids = Reduce(intersect, ds_ids)

common_ids <- results_meta_df %>%
  group_by(block_id, dataset) %>%
  summarise(ids = list(unique(id)), .groups = "drop") %>%
  summarise(common_ids = list(Reduce(intersect, ids)), .by = block_id)
ids_in_all_blocks <- Reduce(intersect, common_ids$common_ids)
