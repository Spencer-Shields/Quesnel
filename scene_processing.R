#----packages----
library(tidyverse)
library(terra)
library(RStoolbox)
library(tools)
library(sf)
library(pbapply)
library(parallel)
library(future)
library(future.apply)
library(tidyterra)
library(Rcpp)
library(tictoc)
library(arrow)
library(data.table)
# library(car)
# library(nortest)
library(pROC)
library(varSel)
library(tseries)
library(lme4)
library(patchwork)
# library(glmmTMB)
#extra functions
source('helper_functions.R')
#functions for checking or visualizing radiometric consistency
source('https://raw.githubusercontent.com/Spencer-Shields/planetscope_radiometric_correction/refs/heads/main/check_radiometric_consistency.R')


#----load data----

#set directories
list.dirs('data')
ps_dir = 'data/planet_scenes'
raw_dir = paste0(ps_dir,'/raw')

#load thinning block vector data
blocks = st_read('data/Quesnel_thinning/12l_12n_bdy.geojson')

thinning_block_ids = blocks$BLOCKNUM #get ids of thinning blocks

#load harvest dates
harvest_dates_file = 'data/Quesnel_thinning/harvest_dates.csv'
harvest_dates_df = read_csv(harvest_dates_file)

#load Otsu thresholds
otsu_df = read_csv("data/Quesnel_thinning/Otsu_3m_change_thresholds.csv")

#load LiDAR data
lidar_dir = 'data/Quesnel_thinning/chm_change_scenes'
lidar_files = list.files(lidar_dir, full.names = T, recursive = T)
lidar_ids = basename(file_path_sans_ext(lidar_files))
names(lidar_files) = lidar_ids

#load non-vegetation mask data
nonveg_mask_dir = 'data/Quesnel_thinning/nonveg_mask_pre=0.5_post=0.3'
nonveg_mask_files = list.files(nonveg_mask_dir, full.names = T)

#get dataframe of scene metadata
meta_df = ps_meta(raw_dir) |> distinct()
meta_df[['acquisition_date']] = as.Date(sapply(meta_df$acquired, function(x){
  s = str_split_1(x, 'T')[1]
  # d = as.Date(s, format = "%Y-%m-%d")
  return(s)
}))
head(meta_df |> select(acquired, acquisition_date))
meta_df = meta_df |>
  mutate(acquisition_year = year(acquisition_date),
         acquisition_month = month(acquisition_date),
         acquisition_day = day(acquisition_date)) %>%
  mutate(acquisition_month_day = paste0(acquisition_month,'-',acquisition_day))


#filter dataframe to remove undesirable scenes

meta_df_filtered = meta_df %>%
  # filter(clear_percent >= clear_threshold) %>% #get data where the percent of clear pixels is greater than or equal to a threshold
  filter(cloud_cover <= 60) |>
  filter(quality_category == 'standard') |> #get standard quality data
  filter(acquisition_year > 2020)


# #see temporal distribution of filtered images
# filtered_dates_bymonth = meta_df_filtered %>%
#   group_by(acquisition_year, acquisition_month) %>%
#   summarise(Number_of_scenes = n())
# 
# ggplot()+
#   geom_col(data = filtered_dates_bymonth, aes(x = acquisition_month, y = Number_of_scenes))+
#   facet_wrap(vars(acquisition_year), ncol = 1) +
#   ggtitle(paste0('Scenes with >=',clear_threshold,'% clear pixels'))
# 
# #see temporal distribution of ALL PlanetScope scenes
# meta_df_dates_bymonth = meta_df %>%
#   group_by(acquisition_year, acquisition_month) %>%
#   summarise(Number_of_scenes = n())
# 
# ggplot()+
#   geom_col(data = meta_df_dates_bymonth, aes(x = acquisition_month, y = Number_of_scenes))+
#   facet_wrap(vars(acquisition_year), ncol = 1) +
#   ggtitle('All scenes')

#----preprocess scenes (avoid saving intermediate data)----

ids = meta_df_filtered$id
raw_rasters = list.files(raw_dir, recursive = T, full.names = T, pattern = '\\.tif$')

if(crs(blocks) != crs(rast(raw_rasters[1]))){ #reproject blocks if necessary to match raster crs
  blocks = st_transform(blocks, crs = crs(rast(raw_rasters[1])))
}

scale_factor = 65535 #scale factor (max value of 16-bit pixel)

raw_bands = names(rast(raw_rasters[!str_detect(raw_rasters,'udm')][1])) #get names of bands in scene rasters (not indices)
indices_to_use = c('BI', 'MSAVI','SR','Hue','NDVI') #indices to test, determined by corr matrix analysis and PCA (below)
feats_to_use = c(raw_bands, indices_to_use)


nonnorm_string = 'Non-normalized'
nonnorm_dir = paste0(ps_dir,'/', nonnorm_string)
dir.check(nonnorm_dir)

z_string = 'Z'
z_dir = paste0(ps_dir,'/', z_string)
dir.check(z_dir)

zr_string = 'Zrobust'
zr_dir = paste0(ps_dir, '/', zr_string)
dir.check(zr_dir)

sm_string = 'SM'
sm_dir = paste0(ps_dir,'/',sm_string)
dir.check(sm_dir)

dirs = c(z_dir, zr_dir, sm_dir, nonnorm_dir)

# dir_strings = c('Z', 'Zrobust', 'SM', 'Non-normalized')
# dirs = paste0(ps_dir,'/',dir_strings)
# lapply(dirs,dir.check)

pblapply(dirs, function(x){ #make subdirectories to store scenes by block
  pblapply(thinning_block_ids, function(y){
    d = paste0(x,'/',y)
    dir.check(d)
  })
})

if(length(list.files(dirs,recursive = T,pattern='\\.tif$')) < length(dirs)*length(thinning_block_ids)*length(ids)){
  #process rasters
  # plan('multisession', workers = 10)
  
  pblapply(1:length(ids), function(i){
    
    id = ids[i]
    
    pblapply(1:length(thinning_block_ids), function(j){
      
      block_dir = thinning_block_ids[j]
      
      #Non-normalized raster
      nonnorm_file = paste0(nonnorm_dir,'/',block_dir,'/',id,'.tif')
      # print(paste('Processing', nonnorm_file))
      if(!file.exists(nonnorm_file)){
        
        #load raster
        f_r = raw_rasters[str_detect(raw_rasters, id) & !str_detect(raw_rasters, 'udm')][1]
        r = rast(f_r)
        
        #mask using UDM
        f_udm = raw_rasters[str_detect(raw_rasters, id) & str_detect(raw_rasters, 'udm')][1]
        udm = rast(f_udm)
        
        r = mask(r, udm[['cloud']]) #remove cloud
        r = mask(r, udm[['shadow']]) #remove shadow
        r = mask(r, udm[['haze_heavy']]) #remove heavy haze
        
        #crop to block
        block = blocks |> filter(BLOCKNUM == block_dir)
        r = crop(r, block, mask = T)
        
        #use nonveg mask
        nvf = nonveg_mask_files[which(str_detect(nonveg_mask_files,block_dir))]
        nvm = rast(nvf)
        nvm_rs = resample(nvm,r)
        r = mask(r, nvm_rs)
        
        
        #calculate vegetation indices
        
        si = spectralIndices(img = r, blue = 'blue', green = 'green', red = 'red', nir = 'nir', redEdge1 = 'rededge'
                             ,scaleFactor = scale_factor, skipRefCheck = T)
        vi = rast.batch.functions(r, include_input = F, fl = c(
          hue.vi,
          gcc.vi,
          ndgr.vi,
          bi.vi,
          ci.vi,
          cri550.vi,
          gli.vi,
          # tvi.vi,
          varig.vi
        ))
        
        r = c(r, si, vi)
        
        #save file
        terra::writeRaster(r, nonnorm_file, datatype='FLT8S')
        
      } else {
        r = rast(nonnorm_file)
      }
      
      #Z-score raster
      z_file = paste0(z_dir,'/',block_dir,'/',id,'.tif')
      # print(paste('Processing',z_file))
      if(!file.exists(z_file)){
        z_r = rast(lapply(r, z_rast))
        writeRaster(z_r, z_file)
      }
      
      #robust Z-score raster
      zr_file = paste0(zr_dir,'/',block_dir,'/',id,'.tif')
      # print(paste('Processing',zr_file))
      if(!file.exists(zr_file)){
        zr_r = rast(lapply(r, function(x)z_rast(x,robust = T)))
        writeRaster(zr_r, zr_file)
      }
      
      #softmax raster
      sm_file = paste0(sm_dir,'/',block_dir,'/',id,'.tif')
      # print(paste('Processing',sm_file))
      if(!file.exists(sm_file)){
        
        #scale bands for use with softmax
        rs = rast(lapply(1:nlyr(r), function(k){
          b = r[[k]]
          nom = names(r)[k]
          if(nom %in% raw_bands){
            b = b/scale_factor
          } else {
            b 
          }
          return(b)
        }))
        
        sm_r = softmax(rs, append_name = F)
        writeRaster(sm_r, sm_file, datatype='FLT8S')
      }
    })
    
  }
  # ,cl = 'future'
  # ,future.seed=T
  )
  
  # plan('sequential')
}


#----timeseries analysis by thinned vs not-thinned pixels----
{
  
  #extract and clean data
  {
    #----extract stats from remote sensing data----
    harvest_threshold = "Otsu3m"
    
    data_string = paste0('Stats_logAUC_JMD_HarvestThreshold=',harvest_threshold,'_EqualClasses')
    data_filename = paste0(ps_dir,'/',data_string,'.arrow')
    
    #run if summary data file does not exist
    if(!file.exists(data_filename)){
      
      #----set up for processing----
      
      #get list of all files
      dirs = c(
        z_dir
        ,zr_dir
        # ,
        # sm_dir
        ,
        nonnorm_dir
      )
      
      all_files <- unlist(sapply(dirs, list.files, pattern = "\\.tif$", full.names = TRUE, recursive = TRUE))
      
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
        r = rast(x)
        
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
      
      # clust = makeCluster(10)
      # plan('cluster', workers = clust)
      
      global_tables_dir = paste0(ps_dir,'/',data_string)
      dir.check(global_tables_dir)
      
      # df_l = pblapply(1:10, function(i){
      pblapply(rev(1:length(all_files_thinning)), function(i){
        
        x = all_files_thinning[i]
        print(paste0('Processing ', x,', ',i,'/',length(all_files_thinning))) #uncomment to identify files where processing fails
        
        id = basename(file_path_sans_ext(x))
        dataset = find_substring(x,dataset_strings) |> str_replace('/','')
        block = find_substring(x, lidar_ids)
        
        
        filename = paste0(global_tables_dir,'/',dataset,'_',block,'_',id,'.arrow')
        if(!file.exists(filename)){
          
          #get raster
          r = rast(x)
          
          #drop extra TVI if it's there
          if(names(r)[35] == 'TVI'){r = r[[-35]]}
          
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
          )
          d1[['block_pixel_stratum']] = 'Thinned'
          d1[['v1']] = rownames(d1)
          
          #calculate global stats for non-thinning pixels
          d2 = global(r_nt, fun = c("max", "min", "mean", "sum", "range", "rms", "sd", "std", "isNA", "notNA"), na.rm=T)
          d2 = cbind(d2, 
                     global(r_nt, median, na.rm=T) |> rename(median = global)
          )
          d2[['block_pixel_stratum']] = 'Not_thinned'
          d2[['v1']] = rownames(d2)
          
          #total block stats
          d3 = global(r, fun = c("max", "min", "mean", "sum", "range", "rms", "sd", "std", "isNA", "notNA"), na.rm=T)
          d3 = cbind(d3, 
                     global(r, median, na.rm=T) |> rename(median = global)
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
          feats = feats[!feats %in% c('max_DN', 'VARIgreen')] #don't run tests on max_DN
          
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
              jeffries_matusita_dist = NA,
              bhattacharya_dist = NA
            )
            
            #get proportion not na
            proportion_na = d1$notNA[d1$v1==y]/d1$isNA[d1$v1==y]
            
            if(proportion_na > 0.05 & length(unique(v[[y]])) > 10){ #proceed if proportion na pixels is greater than threshold and there are multiple cell values
              ##perform tests, store results in test_df
              
              #initialize tibble to store results
              test_df = tibble(
                v1 = y,
                logistic_p = NA,
                logistic_aic = NA,
                logistic_AUC = NA,
                fisher_discriminant_ratio = NA,
                cohen_d = NA,
                jeffries_matusita_dist = NA,
                bhattacharya_dist = NA
              )
              
              #get proportion not na
              proportion_na = d1$notNA[d1$v1==y]/d1$isNA[d1$v1==y]
              
              if(proportion_na > 0.05 & length(unique(v[[y]])) > 10){ #proceed if proportion na pixels is greater than threshold and there are multiple cell values
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
                  log_roc = roc(thinning ~ predictions, log_p, quiet=T)
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
                  jmd = JMdist(g = v$thinning, X = tibble(v[[y]]))
                  test_df$jeffries_matusita_dist = jmd$jmdist
                },
                error = function(e) {
                  message("Jeffries-Matusita distance error: ", e$message)
                  # jeffries_matusita_dist will remain NA as initialized
                })
                
                #Bhattacharya distance
                tryCatch({
                  bhat = BHATdist(g = v$thinning, X = tibble(v[[y]]))
                  test_df$bhattacharya_dist = bhat$bhatdist
                },
                error = function(e) {
                  message("Bhattacharya distance error: ", e$message)
                  # bhattacharya_dist will remain NA as initialized
                })
                
              }
              
              
            }
            
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
      
      # plan('sequential')
      # stopCluster(clust)
      
      
    }
    
    
    #----clean and process dataframe----
    
    results_df = read_feather(data_filename)
    
    # #define harvest dates
    # harvest_dates_df = tribble(
    #   ~block_id, ~harvest_start_date, ~harvest_finish_date,
    #   '12L_C5', '2022-07-31', '2022-11-28',
    #   '12L_D345','2022-07-06', '2022-07-20',
    #   '12L_C4', '2022-10-30', '2023-01-28',
    #   '12L_C7', '2022-09-01', '2022-11-28',
    #   '12L_B8C3', '2022-09-13', '2022-11-02',
    #   '12N_T3', '2023-01-28', '2024-02-29',
    #   '12N_1X1W', '2024-02-29', '2024-03-22',
    #   '12N_V1', '2024-03-04', '2024-04-17'
    # ) |> 
    #   mutate(harvest_start_date = as.Date(harvest_start_date),
    #          harvest_finish_date = as.Date(harvest_finish_date))
    
    results_df = results_df |>
      #add harvest dates
      left_join(harvest_dates_df, by = 'block_id') |>
      #define variable for time before and after harvest_start_date
      mutate(time_since_harvest_start = as.numeric(acquisition_date - harvest_start_date),
             time_since_harvest_finish = as.numeric(acquisition_date - harvest_finish_date)) |>
      #make harvested dummy variable
      mutate(harvested = ifelse(acquisition_date < harvest_start_date,0,1)) |>
      left_join(tibble(acquisition_date = sort(unique(results_df$acquisition_date)),
                       acquisition_date_ordinal = as.numeric(
                         1:length(unique(results_df$acquisition_date))))) |>
      #add acquisition month
      mutate(month = month(acquisition_date)) |>
      #fix the id column since it currently has file paths
      mutate(id = file_path_sans_ext(id))|>
      # #rescale JM distance to give it the range [0-1] so that beta regression can be applied to both AUC and JMD
      mutate(jmd_scaled = jeffries_matusita_dist/sqrt(2)) |>
      filter(!str_detect(dataset,'SM'))
    
    #make dataframe that only has the results from the tests (i.e. remove unnecessary rows)
    test_results_df = results_df |> 
      filter(block_pixel_stratum == 'Total') |> dplyr:: select(-block_pixel_stratum)
    
    #----join metadata table with results to find high quality images----
    
    meta_df = setDT(meta_df) |> mutate(acquisition_date = as.Date(acquisition_date))
    
    results_meta_df = results_df |> left_join(meta_df, by=c('id','acquisition_date'))
    
    max_valid_pixels_byfeat_byblock_bypixstrat = results_meta_df |>
      group_by(block_id,var_name,block_pixel_stratum)|>
      summarise(max_valid_pixels = max(notNA))
    results_meta_df = results_meta_df |> left_join(max_valid_pixels_byfeat_byblock_bypixstrat)
    
    
    notNA_threshold = 0.55
    {
      hq_scenes_df = results_meta_df |>
        # left_join(max_valid_pixels_byfeat_byblock_bypixstrat) |>
        mutate(notNA_proportion = notNA/(notNA+isNA))|>
        filter(var_name %in% feats_to_use,
               block_pixel_stratum=='Total',
               notNA >= notNA_threshold*max_valid_pixels) |>
        pivot_wider(id_cols = c('id','dataset','block_id','acquisition_date'), names_from = var_name, values_from = notNA_proportion) |>
        na.omit() |>
        mutate(idblock = paste(id,block_id))
      
      # ggplot(hq_scenes_df, aes(x = acquisition_date, y=0))+
      #   geom_point()+
      #   facet_wrap(vars(dataset))
      
      hq_scenes_summ = hq_scenes_df |>
        mutate(year_mon = format(acquisition_date, "%Y-%m")) |>
        group_by(dataset, block_id, year_mon) |>
        summarise(count = n())
      ggplot(hq_scenes_summ |> filter(dataset == 'Z')
             , aes(x = year_mon, y = count))+
        geom_point()+
        # geom_line()+
        facet_grid(rows = vars(block_id)
                   # , cols = vars(block_id)
        )
    }
    
    
  }
  
  #random forest over time
  {
    #----random forest analysis over time for all blocks----
    
    results_df_forRF = results_meta_df |>
      mutate(idblock = paste(id,block_id))|>
      filter(
        idblock %in% hq_scenes_df$idblock,
        var_name %in% feats_to_use
        ,block_pixel_stratum!='Total'
      )
    
    #full list of files to test
    all_files_forRF = unique(results_df_forRF$file_path)
    
    library(ranger)
    
    #set up directories
    rf_dir = paste0(ps_dir,'/random_forest_models_nopredict'
                    ,'_',harvest_threshold
    )
    dir.check(rf_dir)
    
    pblapply(all_files_forRF, function(f){
      
      #get dataset from input filepath, create output directory if necessary
      data_dir = str_replace(dirname(dirname(f)),ps_dir,rf_dir)
      dir.check(data_dir)
      
      #get block id from input filepath, create output directory if necessary
      block_dir = str_replace(dirname(f), ps_dir,rf_dir)
      block_id = basename(block_dir)
      dir.check(block_dir)
      
      #get id of scene from filepath, make output filepath
      id = basename(file_path_sans_ext(f))
      rf_filename = paste0(block_dir,'/',id,'.rds')
      rf_filename_null = paste0(block_dir,'/',id,'_NULL.rds')
      
      #run random forest model if output file does not exist
      # cat('\rProcessing ',rf_filename)
      cat('\nProcessing',f)
      if(!file.exists(rf_filename) && !file.exists(rf_filename_null)){
        
        #load and subset raster
        r = rast(f)
        r = r[[feats_to_use]]
        
        #load chm_change file
        chm = lidar_files[str_detect(lidar_files, block_id)] |> rast() #load height change raster
        ht = otsu_df$chm_change_otsu_threshold[otsu_df$block_id==block_id] #get height change threshold
        thinning_binary = ifel(chm <= ht,1,0) #binarize raster with threshold
        names(thinning_binary) = 'thinned' #name layer
        
        #prepare lidar data
        t_r = resample(thinning_binary,r)
        r_ = c(t_r,r)
        v = values(r_,dataframe=T) |> 
          na.omit() |> 
          mutate(thinned = ifelse(thinned == 1, 'Thinned','Not_thinned'))
        v$thinned = as.factor(v$thinned)
        
        #balance classes
        min_class_size <- v %>%
          count(thinned) %>%
          pull(n) %>%
          min()
        v_b <- v %>%
          group_by(thinned) %>%
          slice_sample(n = min_class_size) %>%
          ungroup()
        
        #fit random forest model, message if error thrown
        rf <- tryCatch({
          ranger(thinned ~ ., data = v_b, importance = 'impurity'
                 ,num.threads = floor(detectCores()/2)
          )
        }, error = function(e) {
          message('\nError fitting RF for file: ', f)
          message('Error: ', e$message)
          return(NULL)
        })
        
        # save model if successful
        if(!is.null(rf)){
          saveRDS(rf, rf_filename)
        } else { #save null objects to rds files so that you know that these scenes have been processed
          saveRDS(NULL, rf_filename_null)
        }
        
      }
      cat('\nDone ',match(f,all_files_forRF),'/',length(all_files_forRF))
      # else{
      #     rf = readRDS(rf_filename)
      #     return(rf)
      # }
    })
    
  }
  
  #----focus on analyzing 12N_T3----
  {
    #----preprocess block dataframe----
    id_df_12N_T3 = hq_scenes_df |> filter(block_id == '12N_T3')

    results_12N_T3 = results_meta_df |>
      # filter(notNA >= 0.8*max_valid_pixels) |>
      filter(block_id == '12N_T3'
             ,id %in% id_df_12N_T3$id,
             var_name %in% feats_to_use
             # ,month %in% c(6,7,8,9)
             # ,dataset!='SM'
             ,block_pixel_stratum!='Total'
             # ,acquisition_year!=2021
      ) |>
      filter(!acquisition_date %in% c(
        '2024-04-01'
        ,'2022-11-17'
        ,'2022-11-19'
        ,'2021-12-06'
      ))|>
      filter(!id %in% c(
        '20240401_182741_00_24b2'
      ))

    #get full list of files to test
    all_files_12nt3 = unique(results_12N_T3$file_path)

    #for plotting
    {
      results_nt3_plotting = results_12N_T3 |>
        mutate(var_name = str_replace(var_name, 'coastal_blue', 'cost.blu'))|>
        mutate(var_name = factor(var_name, levels = c('cost.blu', 'blue', 'green_i', 'green',
                                                      'yellow', 'red', 'rededge', 'nir', 
                                                      'BI'
                                                      , 'Hue',
                                                      'MSAVI', 'NDVI', 'SR')))|>
        mutate(block_pixel_stratum = str_replace(block_pixel_stratum,'_',' '))|>
        mutate(block_pixel_stratum = factor(block_pixel_stratum, levels = c('Thinned', 'Not thinned'))) |>
        filter(var_name == 'BI')


    #pixels plot
    {
      z_pix = ggplot(data = results_nt3_plotting |> filter(dataset=='Z'),
                     aes(x = acquisition_date,y=mean))+
        # geom_bin_2d(aes(x = acquisition_date, y = ean, fill = block_pixel_stratum),alpha=0.5)+
        geom_point(aes(color = block_pixel_stratum),alpha=0.2)+
        geom_ribbon(aes(ymin = mean - sd, ymax = mean + sd, fill = block_pixel_stratum), alpha = 0.2, color = NA) +
        geom_vline(aes(xintercept = harvest_start_date, linetype = 'Harvest start'))+
        # geom_vline(aes(xintercept = harvest_finish_date, linetype = 'Harvest end'))+
        facet_grid(rows = vars(var_name), scales = 'free')+
        ggtitle('Z-score')+
        xlab('Acquisition date')+
        theme_classic()+
        ylab('Mean pixel value')+
        theme(strip.text = element_blank(),
              legend.title = element_blank(),
              axis.title.y = element_blank())
      zr_pix = ggplot(data = results_nt3_plotting |>
                        filter(dataset=='Zrobust'),
                      aes(x = acquisition_date,y=mean))+
        # geom_bin_2d(aes(x = acquisition_date, y = mean, fill = block_pixel_stratum),alpha=0.5)+
        geom_point(aes(color = block_pixel_stratum),alpha=0.2)+
        geom_ribbon(aes(ymin = mean - sd, ymax = mean + sd, fill = block_pixel_stratum), alpha = 0.2, color = NA) +
        geom_vline(aes(xintercept = harvest_start_date, linetype = 'Harvest start'))+
        # geom_vline(aes(xintercept = harvest_finish_date, linetype = 'Harvest end'))+
        facet_grid(rows = vars(var_name), scales = 'free')+
        ggtitle('Robust Z-score')+
        ylab(NULL)+
        theme_classic()+
        theme(legend.title = element_blank(),
              axis.title.x = element_blank())
      nn_pix = ggplot(data = results_nt3_plotting |> filter(dataset=='Non-normalized'),
                      aes(x = acquisition_date,y=mean))+
        # geom_bin_2d(aes(x = acquisition_date, y = mean, fill = block_pixel_stratum),alpha=0.5)+
        geom_point(aes(color = block_pixel_stratum),alpha=0.2)+
        geom_ribbon(aes(ymin = mean - sd, ymax = mean + sd, fill = block_pixel_stratum), alpha = 0.2, color = NA) +
        geom_vline(aes(xintercept = harvest_start_date, linetype = 'Harvest start'))+
        # geom_vline(aes(xintercept = harvest_finish_date, linetype = 'Harvest end'))+
        facet_grid(rows = vars(var_name), scales = 'free')+
        ggtitle('Non-normalized')+
        ylab('Mean pixel value')+
        xlab(NULL)+
        theme_classic()+
        theme(strip.text = element_blank(),
              legend.title = element_blank())
      # sm_pix = ggplot(data = results_12N_T3 |> filter(dataset=='SM'),
      #                 aes(x = acquisition_date,y=median))+
      #   # geom_bin_2d(aes(x = acquisition_date, y = mean, fill = block_pixel_stratum),alpha=0.5)+
      #   geom_point(aes(color = block_pixel_stratum),alpha=0.5)+
      #   #geom_ribbon(aes(ymin = mean - sd, ymax = mean + sd, fill = block_pixel_stratum), alpha = 0.2, color = NA) +
      #   geom_vline(aes(xintercept = harvest_start_date, linetype = 'Harvest start'))+
      #   geom_vline(aes(xintercept = harvest_finish_date, linetype = 'Harvest end'))+
      #   facet_grid(rows = vars(var_name), scales = 'free')+
      #   ggtitle('Softmax')+
      #   theme_classic()
      library(patchwork)
      nn_pix+z_pix+zr_pix+plot_layout(guides = 'collect',ncol=3)
      }

    #just BI plot
    {
      bi_df = results_nt3_plotting |> 
        filter(var_name=='BI') |>
        mutate(var_name = ifelse(var_name=='BI', 'Brightness Index', var_name))

      z_pix = ggplot(data = bi_df |> filter(dataset=='Z'),
                     aes(x = acquisition_date,y=mean))+
        # geom_bin_2d(aes(x = acquisition_date, y = ean, fill = block_pixel_stratum),alpha=0.5)+
        geom_point(aes(color = block_pixel_stratum),alpha=0.2)+
        geom_ribbon(aes(ymin = mean - sd, ymax = mean + sd, fill = block_pixel_stratum), alpha = 0.2, color = NA) +
        geom_vline(aes(xintercept = harvest_start_date, linetype = 'Harvest start'))+
        geom_vline(aes(xintercept = harvest_finish_date, linetype = 'Harvest end'))+
        facet_grid(rows = vars(var_name), scales = 'free')+
        ggtitle('Z-score')+
        xlab('Acquisition date')+
        theme_classic()+
        ylab('Mean pixel value')+
        theme(strip.text = element_blank(),
              legend.title = element_blank(),
              axis.title.y = element_blank())
      zr_pix = ggplot(data = bi_df |>
                        filter(dataset=='Zrobust'),
                      aes(x = acquisition_date,y=mean))+
        # geom_bin_2d(aes(x = acquisition_date, y = mean, fill = block_pixel_stratum),alpha=0.5)+
        geom_point(aes(color = block_pixel_stratum),alpha=0.2)+
        geom_ribbon(aes(ymin = mean - sd, ymax = mean + sd, fill = block_pixel_stratum), alpha = 0.2, color = NA) +
        geom_vline(aes(xintercept = harvest_start_date, linetype = 'Harvest start'))+
        geom_vline(aes(xintercept = harvest_finish_date, linetype = 'Harvest end'))+
        facet_grid(rows = vars(var_name), scales = 'free')+
        ggtitle('Robust Z-score')+
        ylab(NULL)+
        theme_classic()+
        theme(legend.title = element_blank(),
              axis.title.x = element_blank())
      nn_pix = ggplot(data = bi_df |> filter(dataset=='Non-normalized'),
                      aes(x = acquisition_date,y=mean))+
        # geom_bin_2d(aes(x = acquisition_date, y = mean, fill = block_pixel_stratum),alpha=0.5)+
        geom_point(aes(color = block_pixel_stratum),alpha=0.2)+
        geom_ribbon(aes(ymin = mean - sd, ymax = mean + sd, fill = block_pixel_stratum), alpha = 0.2, color = NA) +
        geom_vline(aes(xintercept = harvest_start_date, linetype = 'Harvest start'))+
        geom_vline(aes(xintercept = harvest_finish_date, linetype = 'Harvest end'))+
        facet_grid(rows = vars(var_name), scales = 'free')+
        ggtitle('Non-normalized')+
        ylab('Mean pixel value')+
        xlab(NULL)+
        theme_classic()+
        theme(strip.text = element_blank(),
              legend.title = element_blank())
      # sm_pix = ggplot(data = results_12N_T3 |> filter(dataset=='SM'),
      #                 aes(x = acquisition_date,y=median))+
      #   # geom_bin_2d(aes(x = acquisition_date, y = mean, fill = block_pixel_stratum),alpha=0.5)+
      #   geom_point(aes(color = block_pixel_stratum),alpha=0.5)+
      #   #geom_ribbon(aes(ymin = mean - sd, ymax = mean + sd, fill = block_pixel_stratum), alpha = 0.2, color = NA) +
      #   geom_vline(aes(xintercept = harvest_start_date, linetype = 'Harvest start'))+
      #   geom_vline(aes(xintercept = harvest_finish_date, linetype = 'Harvest end'))+
      #   facet_grid(rows = vars(var_name), scales = 'free')+
      #   ggtitle('Softmax')+
      #   theme_classic()
      library(patchwork)
      nn_pix+z_pix+zr_pix+plot_layout(guides = 'collect',ncol=3)
      
      nn_pix/z_pix/zr_pix+plot_layout(guides = 'collect',nrow=3)
      
    }
    }
    
    #gifs
    {
      #---make raster stack of layers to animate
      datasets = c(nonnorm_dir, z_dir, zr_dir) |> basename()
      
      
      gif_files_12nt3 = results_12N_T3 |>
        filter(
          # month %in% c(5,6,7,8,9,10),
          var_name=='BI'
          ,acquisition_year>=2022
          ,!acquisition_date %in% c(
            '2022-07-17',
            '2022-10-06',
            '2023-05-18',
            '2023-07-05',
            '2023-07-11',
            '2023-08-26',
            '2024-05-09',
            '2024-07-09',
            '2024-09-07'
          )
          )
      
      n_pixels_summ = gif_files_12nt3 |> group_by(dataset, notNA) |> summarise(num = n())
      nonnorm_target = n_pixels_summ |> 
        filter(dataset == "Non-normalized") %>%
        slice_max(`num`, n = 1, with_ties = FALSE) %>%
        pull(notNA)
      
      gif_ids = gif_files_12nt3$id[gif_files_12nt3$notNA == nonnorm_target]
      
      gif_files_12nt3 = gif_files_12nt3 |>
        filter(id %in% gif_ids) |>
        select(dataset, id, file_path, acquired) |> distinct()
      
      vi_stacks = pblapply(datasets, function(d){
        dataset_df = gif_files_12nt3 |> 
          filter(dataset==d) |>
          arrange(acquired)
        files = dataset_df$file_path
        stack = pblapply(1:length(files), function(i){
          f = files[i]
          r = rast(f, lyrs ='BI')
        })
        stack = rast(stack)
        names(stack) = dataset_df$acquired
        return(stack)
      })
      names(vi_stacks) = datasets
      
      road_nonveg_mask_files = list.files('data/Quesnel_thinning/nonveg_mask_pre=0.5_post=0.3', full.names = T)
      rnv_f = road_nonveg_mask_files[str_detect(road_nonveg_mask_files, '12N_T3')]
      rnv_m = rast(rnv_f)
      
      vi_stacks = pblapply(vi_stacks, function(r){
        m = resample(rnv_m,r)
        r_m = mask(r,m)
        return(r_m)
        })
      
      library(gifski)
      
      ##save gifs with terra::animate and gifski
      #
      # pause = 0.1
      # nn_gif = terra::animate(vi_stacks$`Non-normalized`, pause)
      # z_gif= terra::animate(vi_stacks$Z, pause)
      # zr_gif= terra::animate(vi_stacks$Zrobust, pause)
      
      # save_gif(expr = animate(vi_stacks$`Non-normalized`, pause),
      #          gif_file = paste0('figures/NN_',pause,'s_22to24.gif'))
      # save_gif(expr = animate(vi_stacks$Z, pause),
      #          gif_file = paste0('figures/Z_',pause,'s_22to24.gif'))
      # save_gif(expr = animate(vi_stacks$Zrobust, pause),
      #          gif_file = paste0('figures/Zrobust_',pause,'s_22to24.gif'))
      
      
      #try again with gganimate
      library(tidyterra)
      library(gganimate)
      
      
      vi_stacks_dfs = pblapply(vi_stacks, function(r){
        d = as.data.frame(r, xy=T) |>
          pivot_longer(
            cols = -c(x, y),
            names_to = "lyr",
            values_to = "value")
      })
      
      gg_anims_l = pblapply(vi_stacks, function(r){
        ggplot()+
          geom_spatraster(data = r)+
          scale_fill_viridis_c(option='viridis', na.value = 'white', name = 'Value')+
          transition_manual(lyr) +
          labs(title = '{current_frame}')+
          theme_void()
      })
      names(gg_anims_l) = datasets
      
      vi_stacks_dfs = pblapply(vi_stacks, function(r){
        d = as.data.frame(r, xy=T) |>
          pivot_longer(
            cols = -c(x, y),
            names_to = "lyr",
            values_to = "value")
      })
      
      gg_anims_l = pblapply(vi_stacks_dfs, function(df){
        ggplot(df) +
          geom_raster(aes(x = x, y = y, fill = value)) +
          scale_fill_viridis_c(option = 'viridis', na.value = 'white', name = 'Value') +
          # transition_states(lyr, transition_length = 0, state_length = 0) +
          transition_manual(lyr)+
          labs(title = "{current_frame}") +
          theme_void()
      })
      names(gg_anims_l) = datasets
      
      
      nn_anim=gganimate::animate(gg_anims_l$`Non-normalized`, nframes = nlyr(vi_stacks$`Non-normalized`), fps = 20)
      # nn_anim1=gganimate::animate(gg_anims_l$`Non-normalized`, fps=1)
      z_anim = gganimate::animate(gg_anims_l$Z, nframes = nlyr(vi_stacks$`Non-normalized`), fps=20)
      zr_anim = gganimate::animate(gg_anims_l$Zrobust, nframes = nlyr(vi_stacks$`Non-normalized`), fps=20)
      
      
      pblapply(1:length(list(nn_anim, z_anim, zr_anim)), function(i){
        an = 
        anim_save(filename = paste0('figures/',names))
      })
      
      anim_save(filename='figures/NN_22to24_allscenes_12NT3.gif',animation = nn_anim)
      anim_save(filename='figures/Z_22to24_allscenes_12NT3.gif',animation = z_anim)
      anim_save(filename='figures/Zrobust_22to24_allscenes_12NT3.gif',animation = zr_anim)
      
    }
    
    #----do random forest for each date for 12N_T3----
    
    library(ranger)

    #set up directories
    rf_dir = paste0(ps_dir,'/random_forest_models_nopredict'
                    ,'_',harvest_threshold
                    )
    dir.check(rf_dir)

    #load chm_change file
    chm = lidar_files[str_detect(lidar_files, '12N_T3')] |> rast()
    ht = otsu_df$chm_change_otsu_threshold[otsu_df$block_id=='12N_T3']
    thinning_binary = ifel(chm <= ht,1,0)
    names(thinning_binary) = 'thinned'

    pblapply(1:length(all_files_12nt3), function(i){
      f = all_files_12nt3[i]
      data_dir = str_replace(dirname(dirname(f)),ps_dir,rf_dir)
      dir.check(data_dir)
      block_dir = str_replace(dirname(f), ps_dir,rf_dir)
      dir.check(block_dir)

      id = basename(file_path_sans_ext(f))
      rf_filename = paste0(block_dir,'/',id,'.rds')
      cat('\rProcessing ',rf_filename,' ',i,'/',length(all_files_12nt3))
      if(!file.exists(rf_filename)){

        r = rast(f)
        r = r[[feats_to_use]]
        t_r = resample(thinning_binary,r)
        r_ = c(t_r,r)
        v = values(r_,dataframe=T) |>
          na.omit() |>
          mutate(thinned = ifelse(thinned == 1, 'Thinned','Not_thinned'))
        v$thinned = as.factor(v$thinned)

        #balance classes
        min_class_size <- v %>%
          count(thinned) %>%
          pull(n) %>%
          min()
        v_b <- v %>%
          group_by(thinned) %>%
          slice_sample(n = min_class_size) %>%
          ungroup()

        #fit random forest model, message if error thrown
        rf <- tryCatch({
          ranger(thinned ~ ., data = v_b
                 , importance = 'impurity'
                 ,write.forest = F
                 ,verbose = F
                 # ,num.trees = 300
                 ,num.threads = floor(detectCores/2)
          )
        }, error = function(e) {
          message('\nError fitting RF for file: ', f)
          message('Error: ', e$message)
          return(NULL)
        })

        # save model if successful
        if (!is.null(rf)) saveRDS(rf, rf_filename)

      }
    })
    
    
    #----analyze random forest results----
    rffl = list.files(rf_dir, full.names=T, recursive=T, pattern = '\\.rds$')
    rf_l = lapply(rffl, readRDS)
    
    #make tibble to store info
    rf_ids = basename(file_path_sans_ext(rffl))
    
    rf_df = tibble(
      id = rf_ids,
      dataset = basename(dirname(dirname(rffl))),
      block_id = basename(dirname(rffl)),
      OOBE = sapply(rf_l, function(x)x$prediction.error)
    ) |>
      mutate(acquisition_date = as.Date(paste0( #get acquisition date
        substr(id, 1, 4), "-",  # year
        substr(id, 5, 6), "-",  # month
        substr(id, 7, 8)        # day
      ))) |>
      mutate(overall_accuracy = 1 - OOBE) |>
      left_join(results_meta_df, by = c('id', 'dataset', 'block_id', 'acquisition_date')) |>
      filter(
        # month %in% c(5,6,7,8,9),
        dataset != 'SM')
    
    #train interrupted time series model for each RF dataset
    datasets = str_remove(dataset_strings, '/')
    datasets = datasets[!str_detect(datasets, 'SM')]
    datasets=c('Z','Zrobust','Non-normalized')
    
    rf_lm_l = pblapply(datasets, function(x){
      rf_df_sub = rf_df |>
        filter(dataset == x)
      lin_mod = lm(overall_accuracy ~ time_since_harvest_start * harvested, data = rf_df_sub)
    })
    names(rf_lm_l) = datasets
    
    #function for extracting model coefficients
    extract_coefs_lm = function(a){ #x is a linear model created by lm
      summ = summary(a)
      m = summ$coefficients
      df = as.data.frame(m)
      df$term = row.names(summ$coefficients)
      df$r.squared = summ$r.squared
      return(df)
    }
    
    rf_lm_coeffs = lapply(datasets, function(x){
      m = rf_lm_l[[x]]
      d = extract_coefs_lm(m)
      d$dataset = x
      return(d)
    }) |> 
      bind_rows() |>
      mutate(Estimate_rounded = round(Estimate,3))
    
    rf_lm_coeffs_long = rf_lm_coeffs |>
      pivot_longer(cols = c('Estimate'))
    
    ggplot(rf_lm_coeffs)+
      geom_tile(aes(x = dataset, y = term, fill = Estimate))+
      geom_text(aes(x = dataset, y = term, label = Estimate))
    
    #make predicted values (to visualize model)
    rf_pred = pblapply(datasets, function(x){
      rf_df_sub = rf_df %>% 
        filter(dataset == x)
      model = rf_lm_l[[x]]
      rf_df_sub$predicted_OA = predict(model, newdata = rf_df_sub)
      return(rf_df_sub)
    }) |> bind_rows()
    
    #make counterfactual predictions (to visualize model)
    rf_counterfactual_lm = pblapply(datasets, function(x){
      rf_df_sub = rf_df |>
        filter(acquisition_date<harvest_start_date)|>
        filter(dataset == x)
      lin_mod = lm(overall_accuracy ~ time_since_harvest_start * harvested, data = rf_df_sub)
    })
    names(rf_counterfactual_lm) = datasets
    rf_counterfactual_pred = pblapply(datasets, function(x){
      rf_df_sub = rf_df %>% 
        filter(dataset == x)
      model = rf_counterfactual_lm[[x]]
      rf_df_sub$counterfactual_OA = predict(model, newdata = rf_df_sub)
      return(rf_df_sub)
    }) |> bind_rows()
    
    rf_df = left_join(rf_df, rf_pred) |> left_join(rf_counterfactual_pred)
    
    rf_df_plotting = rf_df |>
      mutate(dataset = case_match(dataset,
                                  'Z' ~ 'Z-score',
                                  'Zrobust' ~ 'Robust Z-score',
                                  .default='Non-normalized'))
    
    ggplot(rf_df_plotting, aes(x = acquisition_date, y = overall_accuracy))+
      geom_point(shape = 1, alpha = 0.5)+
      geom_line(aes(y = counterfactual_OA, color='Counterfactual prediction'))+
      geom_line(aes(y = predicted_OA, color = 'Linear prediction'))+
      geom_vline(aes(xintercept = harvest_finish_date, linetype = 'Harvest end'))+
      geom_vline(aes(xintercept = harvest_start_date, linetype = 'Harvest start'))+
      ylim(0.5,1)+
      facet_grid(rows = vars(dataset))+
      guides(linetype = guide_legend(title = NULL))+
      xlab('Acquisition date')+
      ylab('Overall accuracy')+
      theme_classic()
    
    
    
    
    
    #----get RMSE between adjacent scenes----
    
    delta_dir = paste0(ps_dir,'/delta_scenes')
    dir.check(delta_dir)
    
    delta_globs_dir = paste0(ps_dir,'/delta_global_vals')
    dir.check(delta_globs_dir)
    
    # cl = makeCluster(4)
    # plan('cluster', workers = cl)
    plan('multisession', workers=4)
    
    future_lapply(dirs, function(d){
      dat_dir = paste0(delta_dir,'/',basename(d))
      dir.check(dat_dir)
      
      block_indir_l = list.dirs(d)
      block_indir_l <- block_indir_l[block_indir_l != d]               # remove top-level dir
      block_indir_l <- block_indir_l[str_detect(block_indir_l, "12N_T3")]
      
      
      pblapply(block_indir_l, function(b){
        block_dir = paste0(dat_dir,'/',basename(b))
        dir.check(block_dir)
        
        #get list of scenes in dir
        scenes_l = list.files(b, full.names = T, pattern = '\\.tif$')
        
        #subset dataframe to only include items which appear in this list
        df_sub = results_12N_T3 |> 
          filter(file_path %in% scenes_l) |>
          distinct(file_path,id,acquisition_date) |>
          #sort in increasing order of 
          arrange(acquisition_date, file_path)
        
        pblapply(2:nrow(df_sub), function(i){
          filename = paste0(block_dir,'/',df_sub$id[i],'.tif')
          vals_file = paste0(delta_globs_dir,'/',basename(d),'_',basename(b),'_',df_sub$id[i],'.arrow')
          
          cat('\rProcessing',df_sub$file_path[i-1], '&', df_sub$file_path[i])
          # if(!file.exists(filename)){
          r = rast(df_sub$file_path[i])
          r = r[[feats_to_use]]
          
          r_past = rast(df_sub$file_path[i-1])
          r_past = r_past[[feats_to_use]]
          
          r_past = resample(r_past,r)
          r_c = crop(r,r_past,mask=T)
          r_past_c = crop(r_past,r,mask=T)
          
          if(!file.exists(filename)){
            dr = r_c - r_past_c
            writeRaster(dr, filename)
          } else {
            dr = rast(filename)
          }
          
          if(!file.exists(vals_file)){
            
            
            dv = global(dr, fun = c("max", "min", "mean", "sum", "range", "rms", "sd", "std", "isNA", "notNA"), na.rm=T)
            dv = cbind(dv, 
                       global(dr, median, na.rm=T) |> rename(median = global))
            dv$dataset = d
            dv$block_id = b
            dv$acquisition_date = df_sub$acquisition_date[i]
            dv$id = df_sub$id[i]
            dv$file_path = df_sub$file_path[i]
            dv$var_name = rownames(dv)
            
            dv = cbind(dv, rmse(r_c, r_past_c, normalise_rmse = T, method = 'iqr'))
            
            
            if(!file.exists(vals_file)){write_feather(dv,vals_file)}
          }
        })
      })
    })
    
    plan('sequential')
    # stopCluster(cl)
    
  }
  
  
  #----analyze feature importance and correlations between features----
  
  #get vector of high-quality images from each block
  max_valid_pixels_byblock = results_meta_df |>
    group_by(block_id)|>
    summarise(max_valid_pixels = max(notNA))
  
  qual_df = results_meta_df |>
    left_join(max_valid_pixels_byblock, by=c('block_id')) #|>
    filter(notNA >= 0.95*max_valid_pixels,
           # month %in% c(6,7,8,9)
    ) |>
    distinct(id, block_id, month, acquisition_year, keep_all=T)
  
  #load high-quality rasters
  hq_rl = lapply(1:nrow(qual_df), function(i){
    id = high_qual_df$id[i]
    block = high_qual_df$block_id[i]
    path = paste0(sm_dir,'/',block,'/',id,'.tif')
    r = rast(path)
    if(names(r)[35]=='TVI'){r = r[[-35]]}
    if(names(r)[33]=='TVI'){r = r[[-33]]}
    return(r)
  })
  
  #inspect plots
  il = lapply(hq_rl, function(x)x[['RVI']])
  names(il) = high_qual_df$id
  list_plot(il)
  
  #get values from each scene
  hq_vl = pblapply(hq_rl, function(x)values(x,dataframe=T))
  hq_v = rbindlist(hq_vl) |>
    select(-EVI2,-NRVI,-VARIgreen) |> na.omit()
  nrow(hq_v)
  
  #principal component analysis
  p = prcomp(hq_v|>
               select(-coastal_blue,-blue,-green_i,-green,-yellow,-red,-rededge,-nir)#exclude raw bands
  )
  summary(p) #most variance explained by first two PCs
  rotations = p$rotation#[,1:2]
  rotations_df = as.data.frame(rotations) |> 
    mutate(var_name = rownames(rotations))
  rotations_df_long = rotations_df |>
    pivot_longer(cols = colnames(rotations),names_to = 'PC',values_to = 'rotation')
  
  ggplot(rotations_df_long |> filter(PC %in% c('PC1','PC2')))+
    geom_tile(aes(x = PC,y=var_name,fill=rotation))+
    scale_fill_viridis_c()
  
  library(ggrepel)
  ggplot(rotations_df, aes(x = PC1, y = PC2, label=var_name))+
    geom_point()+
    geom_label_repel(max.overlaps = Inf
                     ,box.padding = 0.7
                     ,min.segment.length = 0.2
    )
  
  #do correlation matrix
  library(corrplot)
  
  cm = cor(hq_v, method = 'pearson')
  corrplot(cm)
  
  cm_sub=cm[9:nrow(cm),9:ncol(cm)]
  corrplot(cm_sub)
  
  
  cm_sub_summ = tibble(var_name = colnames(cm_sub),
                       mean_cor = rowMeans(cm_sub),
                       tot_cor = rowSums(cm_sub)) |>
    mutate(abs_mean_cor = abs(mean_cor))
  cm_sub_summ_long =cm_sub_summ |> 
    pivot_longer(cols = c('mean_cor', 'tot_cor'), names_to = 'val_type', values_to = 'score')
  ggplot()+
    geom_tile(data=cm_summ_long |> filter(val_type=='mean_cor'), aes(x = val_type, y = var_name, fill = score))+
    scale_fill_viridis_c()
  
  
  indices_to_use = c('BI', 'MSAVI','SR','Hue','NDVI')
  #----plotting vegetation indices over time----
  {
    indices = c(
      # "blue", "green", "red"
      # ,
      # "Hue"
      # ,
      # "GCC"
      # ,
      # "NDGR"
      # ,
      "BI"
      # ,
      # "CI", "CRI550",
      # "GLI", "TVI"
      # , "VARIgreen"
    ) #vector of indices to look at

    blocks_ = c(
      "12L_C5"
      ,
      "12L_D345"
      ,  "12L_C4",    "12L_C7",    "12L_B8C3",  "12N_T3",    "12N_1X1W",  "12N_V1"
      # ,"NoChange1.1_conifer", "NoChange1.2_conifer",
      # "NoChange4_conifer",   "NoChange5_conifer",
      # "NoChange6_conifer",   "NoChange7_conifer"
    )

    tests_ = c(
      'levene_p',
      'anderson_darling_p_thinning',
      'anderson_darling_p_notthinning',
      'anderson_darling_p_pooled',
      'student_t_p',
      'welch_t_p',
      'mann_whitney_p',
      'anova_p',
      'kruskal_wallis_p',
      'logistic_p',
      'logistic_aic',
      'logistic_AUC',
      'fisher_discriminant_ratio',
      'cohen_d',
      'jeffries_matusita_dist',
      'divergence'
    )

    datasets_ = c(
      'Z'
      # ,
      # 'Z_delta'
      # ,
      # 'Zrobust'
      # , 'Zrobust_delta'
      # ,
      # 'SM'
      # ,
      # 'SM_delta'
      # 'Non-normalized'
    )

    months_ = c(
      '01','02','03',
      '04',
      '05','06','07','08','09','10','11'
      ,'12'
    )

    # harvest_dates_df = tibble(
    #   block_id = blocks_p$BLOCKNUM[!str_detect(blocks_p$BLOCKNUM, 'NoChange')],
    #   harvest_month = c('2023-02'
    #                     ,'2022-07'
    #                     ,'2023-02'
    #                     ,'2023-02'
    #                     ,'2023-02'
    #                     ,'2023-02'
    #                     ,'2024-03'
    #                     ,'2024-03')
    #   ,harvest_note = c('-','-','-','-','-','partial, complete 2023-03','-','partial, complete 2024-04')
    # )

    pixel_stratum = c('Not_thinned', 'Thinned'
                      # , 'Total'
    )

    subset_df = results_df %>%
      filter(var_name %in% indices
             , str_detect(block_id, paste(blocks_, collapse = "|"))
             , dataset %in% datasets_
             , month %in% months_
             , block_pixel_stratum %in% pixel_stratum)

    #harvest pixel values
    pixels_p = ggplot(subset_df, aes(x = acquisition_date, y = mean)) +
      geom_point(aes(color = block_pixel_stratum))+
      geom_line(aes(color = block_pixel_stratum))+
      geom_ribbon(aes(ymin = mean - sd, ymax = mean + sd, fill = block_pixel_stratum), alpha = 0.2, color = NA) +
      # geom_ribbon(aes(x = acquisition_date, ymin = mean-sds, ymax = mean+sds))+
      scale_x_date(
        date_breaks = "1 year",       # Major tick marks every year
        date_labels = "%Y %m",           # Year labels on major ticks
        date_minor_breaks = "1 month" # Minor tick marks every month
      )+
      facet_grid(rows = vars(block_id), cols = vars(var_name), scales = 'free') +
      # geom_vline(xintercept = as.Date('2022-07-01')) + #12L_D345
      geom_vline(aes(xintercept = harvest_finish_date, linetype = 'Harvest end date')) + #get harvest date for each plot
      geom_vline(aes(xintercept = harvest_start_date, linetype = 'Harvest start date')) +
      ggtitle(unique(subset_df$dataset))+
      theme_classic()
    
    #logistic_AUC results
    alpha = 0.5
    log_p = ggplot(subset_df |> filter(block_pixel_stratum == 'Thinned') #remove redundant info (since thinning, not_thinning, and total have same stats info)
                   , aes(x = acquisition_date, y = logistic_AUC)) +
      geom_point()+
      # geom_line()+
      geom_point(data = subset_df|>
                   mutate(logistic_sig = ifelse(logistic_p < alpha,1,0)) |>
                   filter(logistic_sig == 1)|>
                   mutate(label = paste0('p<',alpha)),
                 aes(x = acquisition_date, y =1, color=label), shape = 8)+ #significance of model
      # scale_shape_manual(name = "Model Significance", values = c("Significant Model" = 8))+
      scale_x_date(
        date_breaks = "1 year",       # Major tick marks every year
        date_labels = "%Y %m",           # Year labels on major ticks
        date_minor_breaks = "1 month" # Minor tick marks every month
      )+
      facet_grid(rows = vars(block_id), cols = vars(var_name)
                 , scales = 'free'
      ) +
      ylim(0.3,1)+
      geom_vline(aes(xintercept = harvest_finish_date, linetype = 'Harvest end date')) + #get harvest date for each plot
      geom_vline(aes(xintercept = harvest_start_date, linetype = 'Harvest start date')) +
      ggtitle(unique(subset_df$dataset))+
      theme_classic()
    
    #jeffries-matusita results
    jm_p = ggplot(subset_df |> filter(block_pixel_stratum == 'Thinned') #remove redundant info (since thinning, not_thinning, and total have same stats info)
                  , aes(x = acquisition_date, y = jmd_scaled)) +
      geom_point()+
      # geom_line()+
      # scale_shape_manual(name = "Model Significance", values = c("Significant Model" = 8))+
      scale_x_date(
        date_breaks = "1 year",       # Major tick marks every year
        date_labels = "%Y %m",           # Year labels on major ticks
        date_minor_breaks = "1 month" # Minor tick marks every month
      )+
      facet_grid(rows = vars(block_id), cols = vars(var_name)
                 , scales = 'free'
      ) +
      ylim(0,0.3)+
      geom_vline(aes(xintercept = harvest_finish_date, linetype = 'Harvest end date')) + #get harvest date for each plot
      geom_vline(aes(xintercept = harvest_start_date, linetype = 'Harvest start date')) +
      ggtitle(unique(subset_df$dataset))+
      theme_classic()
    
    #fisher discriminant ratio
    fdr_p = ggplot(subset_df |> filter(block_pixel_stratum == 'Thinned') #remove redundant info (since thinning, not_thinning, and total have same stats info)
                   , aes(x = acquisition_date, y = fisher_discriminant_ratio)) +
      geom_point()+
      # geom_line()+
      # scale_shape_manual(name = "Model Significance", values = c("Significant Model" = 8))+
      scale_x_date(
        date_breaks = "1 year",       # Major tick marks every year
        date_labels = "%Y %m",           # Year labels on major ticks
        date_minor_breaks = "1 month" # Minor tick marks every month
      )+
      facet_grid(rows = vars(block_id), cols = vars(var_name)
                 , scales = 'free'
      ) +
      ylim(0,0.35)+
      geom_vline(aes(xintercept = harvest_finish_date, linetype = 'Harvest end date')) + #get harvest date for each plot
      geom_vline(aes(xintercept = harvest_start_date, linetype = 'Harvest start date')) +
      ggtitle(unique(subset_df$dataset))+
      theme_classic()
    
    library(patchwork)
    pixels_p+log_p+jm_p+plot_layout(guides = 'collect')
    
  }

  #----statistical modelling to assess effect of harvesting on class separability----
  #define features and datasets that should be tested
  ps_feats = unique(test_results_df$var_name)
  ps_feats = ps_feats[!ps_feats %in% c('VARIgreen', 'max_DN', 'max')]

  datasets = unique(test_results_df$dataset)
  datasets = datasets[datasets %in% c('Z', 'Non-normalized', 'SM', 'Zrobust')]

  block_ids = unique(test_results_df$block_id)
  
  #----LMMs for AUC and JM distance----
  
  
  
  # plan('multisession', workers = 8)
  jmd_lm_list = pblapply(1:length(datasets), function(i){

    dsi = datasets[i]

    lm_l = pblapply(1:length(ps_feats), function(j){

      varj = ps_feats[j]

      #print for debugging
      print(paste0(
        'Processing ',i,'-',j,', ',dsi,'-',varj
      ))

      dat = test_results_df |> filter(dataset == dsi,
                                      var_name == varj)


      # jmd_lm = glmmTMB(jmd_scaled ~ acquisition_date_ordinal + harvested + (acquisition_date_ordinal * harvested | block_id)
      #                  , data = dat
      #                  ,family=ordbeta(link="logit")) #ordered beta regression since the dataset includes zeroes and ones

      jmd_lm = lmer(jmd_scaled ~ 1 + months_since_harvest * harvested +
                      (1 + months_since_harvest * harvested | block_id)
                    , data = dat)

      return(jmd_lm)
    })
    names(lm_l) = ps_feats
    return(lm_l)
  }
  # ,cl = 'future'
  )
  names(jmd_lm_list) = datasets

  auc_lm_list = pblapply(1:length(datasets), function(i){

    dsi = datasets[i]

    lm_l = pblapply(1:length(ps_feats), function(j){

      varj = ps_feats[j]

      #print for debugging
      print(paste0(
        'Processing ',i,'-',j,', ',dsi,'-',varj
      ))

      dat = test_results_df |> filter(dataset == dsi,
                                      var_name == varj)

      # auc_lm = glmmTMB(logistic_AUC ~ acquisition_date_ordinal + harvested + (acquisition_date_ordinal * harvested | block_id)
      #                  , data = dat
      #                  ,family=beta_family(link="logit"))

      auc_lm = lmer(logistic_AUC ~ 1 + months_since_harvest * harvested +
                      (1 + months_since_harvest * harvested | block_id)
                    , data = dat)

      return(auc_lm)
    })
    names(lm_l) = ps_feats
    return(lm_l)
  }
  # ,cl = 'future'
  )
  names(auc_lm_list) = datasets
  # plan('sequential')

  #extract model coefficients and other info

  extract_coefs_lme4 = function(L){ #L is a nested structure with lme4 models in lists according to feature in lists according to dataset
    coefs_df_l = pblapply(datasets, function(x){
      df_l = pblapply(ps_feats, function(y){
        m = L[[x]][[y]]
        s = summary(m)
        coefs_df = s$coefficients |> as.data.frame()
        coefs_df$model_terms = rownames(coefs_df)
        coefs_df$AIC = AIC(m)
        coefs_df$REML = REMLcrit(m)
        coefs_df$dataset = x
        coefs_df$var_name = y

        return(coefs_df)
      })
      df = bind_rows(df_l)
      return(df)
    })
    df = bind_rows(coefs_df_l)
    rownames(df) = NULL
    return(df)
  }

  jmd_coefs = extract_coefs_lme4(jmd_lm_list) |>
    mutate(Estimate_rounded = round(Estimate, 3)) |>
    mutate(metric = 'Area Under Curve')
  auc_coefs = extract_coefs_lme4(auc_lm_list) |>
    mutate(Estimate_rounded = round(Estimate, 3)) |>
    mutate(metric = 'Jeffries-Matusita distance')

  coefs_combined_df = bind_rows(jmd_coefs, auc_coefs)

  #plot results
  {
    ggplot(data = auc_coefs |> filter(!model_terms %in% c('months_since_harvest', '(Intercept)')), aes(x = dataset, y = var_name))+
      geom_tile(aes(fill = Estimate_rounded))+
      geom_text(aes(label = Estimate_rounded))+
      scale_fill_viridis_b()+
      facet_grid(vars(model_terms))+
      theme_classic()

    # ggplot(data = auc_coefs |> filter(model_terms %in% c('(Intercept)')), aes(x = dataset, y = var_name))+
    #   geom_tile(aes(fill = Estimate_rounded))+
    #   geom_text(aes(label = Estimate_rounded))+
    #   scale_fill_viridis_b()+
    #   facet_grid(vars(model_terms))+
    #   theme_classic()
  }
  
  
  
  
  
  
  
  #----LMMs for AUC, JM distance, and fisher discriminant ratio----
  ps_feats = unique(test_results_df$var_name)
  ps_feats = ps_feats[!ps_feats %in% c('VARIgreen', 'max_DN', 'max')]
  
  datasets = unique(test_results_df$dataset)
  datasets = datasets[datasets %in% c('Z', 'Non-normalized', 'SM', 'Zrobust')]
  
  # models from whole timeseries
  { 
    # datasets = c(nonnorm_dir, z_dir, zr_dir) |> basename()
    
    
    # plan('multisession', workers = 8)
    jmd_lm_list = pblapply(1:length(datasets), function(i){
      
      dsi = datasets[i]
      
      lm_l = pblapply(1:length(ps_feats), function(j){
        
        varj = ps_feats[j]
        
        #print for debugging
        print(paste0(
          'Processing ',i,'-',j,', ',dsi,'-',varj
        ))
        
        dat = test_results_df |> filter(dataset == dsi,
                                        var_name == varj,
                                        !is.na(jeffries_matusita_dist))
        
        
        # jmd_lm = glmmTMB(jmd_scaled ~ acquisition_date_ordinal + harvested + (acquisition_date_ordinal * harvested | block_id)
        #                  , data = dat
        #                  ,family=ordbeta(link="logit")) #ordered beta regression since the dataset includes zeroes and ones
        
        if(nrow(dat)>4*length(unique(dat$block_id))){
          jmd_lm = lmer(jmd_scaled ~ 1 + time_since_harvest_start * harvested +
                          (1 + time_since_harvest_start * harvested | block_id)
                        , data = dat)
          
          return(jmd_lm)
        }
      })
      names(lm_l) = ps_feats
      return(lm_l)
    }
    # ,cl = 'future'
    )
    names(jmd_lm_list) = datasets
    
    auc_lm_list = pblapply(1:length(datasets), function(i){
      
      dsi = datasets[i]
      
      lm_l = pblapply(1:length(ps_feats), function(j){
        
        varj = ps_feats[j]
        
        #print for debugging
        print(paste0(
          'Processing ',i,'-',j,', ',dsi,'-',varj
        ))
        
        dat = test_results_df |> filter(dataset == dsi,
                                        var_name == varj,
                                        !is.na(logistic_AUC))
        
        # auc_lm = glmmTMB(logistic_AUC ~ acquisition_date_ordinal + harvested + (acquisition_date_ordinal * harvested | block_id)
        #                  , data = dat
        #                  ,family=beta_family(link="logit"))
        
        if(nrow(dat)>4*length(unique(dat$block_id))){
          auc_lm = lmer(logistic_AUC ~ 1 + time_since_harvest_start * harvested +
                          (1 + time_since_harvest_start * harvested | block_id)
                        , data = dat)
          
          return(auc_lm)
        }
        
      })
      names(lm_l) = ps_feats
      return(lm_l)
    }
    # ,cl = 'future'
    )
    names(auc_lm_list) = datasets
    # plan('sequential')
    
    fdr_lm_list = pblapply(1:length(datasets), function(i){
      
      dsi = datasets[i]
      
      lm_l = pblapply(1:length(ps_feats), function(j){
        
        varj = ps_feats[j]
        
        #print for debugging
        print(paste0(
          'Processing ',i,'-',j,', ',dsi,'-',varj
        ))
        
        dat = test_results_df |> filter(dataset == dsi,
                                        var_name == varj,
                                        !is.na(fisher_discriminant_ratio))
        
        # auc_lm = glmmTMB(logistic_AUC ~ acquisition_date_ordinal + harvested + (acquisition_date_ordinal * harvested | block_id)
        #                  , data = dat
        #                  ,family=beta_family(link="logit"))
        
        if(nrow(dat)>4*length(unique(dat$block_id))){
          auc_lm = lmer(fisher_discriminant_ratio ~ 1 + time_since_harvest_start * harvested +
                          (1 + time_since_harvest_start * harvested | block_id)
                        , data = dat)
          
          return(auc_lm)
        }
        
      })
      names(lm_l) = ps_feats
      return(lm_l)
    }
    # ,cl = 'future'
    )
    names(fdr_lm_list) = datasets
    
    #extract model coefficients and other info
    
    extract_coefs_lme4 = function(L){ #L is a nested structure with lme4 models in lists according to feature in lists according to dataset
      coefs_df_l = pblapply(datasets, function(x){
        df_l = pblapply(ps_feats, function(y){
          m = L[[x]][[y]]
          if(!is.null(m)){
            s = summary(m)
            coefs_df = s$coefficients |> as.data.frame()
            coefs_df$model_terms = rownames(coefs_df)
            coefs_df$AIC = AIC(m)
            coefs_df$REML = REMLcrit(m)
            coefs_df$dataset = x
            coefs_df$var_name = y
            
            return(coefs_df)
          }
        })
        df = bind_rows(df_l)
        return(df)
      })
      df = bind_rows(coefs_df_l)
      rownames(df) = NULL
      return(df)
    }
    
    jmd_coefs = extract_coefs_lme4(jmd_lm_list) |>
      mutate(Estimate_rounded = round(Estimate, 3)) |>
      mutate(metric = 'Area Under Curve')
    auc_coefs = extract_coefs_lme4(auc_lm_list) |>
      mutate(Estimate_rounded = round(Estimate, 3)) |>
      mutate(metric = 'Jeffries-Matusita distance')
    fdr_coefs = extract_coefs_lme4(fdr_lm_list) |>
      mutate(Estimate_rounded = round(Estimate, 3)) |>
      mutate(metric = 'Fisher discriminant ratio')
    
    coefs_combined_df = bind_rows(jmd_coefs, auc_coefs, fdr_coefs)
    
    #plot results
    {
      ggplot(data = coefs_combined_df |> 
               filter(!model_terms %in% c('months_since_harvest', '(Intercept)')
                      ,var_name %in% feats_to_use
               )
             , aes(
               x = dataset,
               y = var_name))+
        geom_tile(aes(fill = Estimate_rounded))+
        geom_text(aes(label = Estimate_rounded))+
        scale_fill_viridis_b()+
        # facet_grid(rows= vars(model_terms), cols = vars(dataset), scales = 'free_x')+
        facet_grid(rows = vars(model_terms), cols = vars(metric), scales = 'free')+
        ggtitle('ITS model coefficients (whole time series)')+
        theme_classic()
      
      # ggplot(data = auc_coefs |> filter(model_terms %in% c('(Intercept)')), aes(x = dataset, y = var_name))+
      #   geom_tile(aes(fill = Estimate_rounded))+
      #   geom_text(aes(label = Estimate_rounded))+
      #   scale_fill_viridis_b()+
      #   facet_grid(vars(model_terms))+
      #   theme_classic()
    }
  }
  
  # models from growing season
  { 
    gs_months = c(5,6,7,8,9,10)
    test_results_df_gs = test_results_df |> filter(month %in% gs_months)
    
    # plan('multisession', workers = 8)
    jmd_lm_list = pblapply(1:length(datasets), function(i){
      
      dsi = datasets[i]
      
      lm_l = pblapply(1:length(ps_feats), function(j){
        
        varj = ps_feats[j]
        
        #print for debugging
        print(paste0(
          'Processing ',i,'-',j,', ',dsi,'-',varj
        ))
        
        dat = test_results_df_gs |> filter(dataset == dsi,
                                        var_name == varj,
                                        !is.na(jeffries_matusita_dist))
        
        
        # jmd_lm = glmmTMB(jmd_scaled ~ acquisition_date_ordinal + harvested + (acquisition_date_ordinal * harvested | block_id)
        #                  , data = dat
        #                  ,family=ordbeta(link="logit")) #ordered beta regression since the dataset includes zeroes and ones
        
        if(nrow(dat)>4*length(unique(dat$block_id))){
          jmd_lm = lmer(jmd_scaled ~ 1 + time_since_harvest_start * harvested +
                          (1 + time_since_harvest_start * harvested | block_id)
                        , data = dat)
          
          return(jmd_lm)
        }
      })
      names(lm_l) = ps_feats
      return(lm_l)
    }
    # ,cl = 'future'
    )
    names(jmd_lm_list) = datasets
    
    auc_lm_list = pblapply(1:length(datasets), function(i){
      
      dsi = datasets[i]
      
      lm_l = pblapply(1:length(ps_feats), function(j){
        
        varj = ps_feats[j]
        
        #print for debugging
        print(paste0(
          'Processing ',i,'-',j,', ',dsi,'-',varj
        ))
        
        dat = test_results_df_gs |> filter(dataset == dsi,
                                        var_name == varj,
                                        !is.na(logistic_AUC))
        
        # auc_lm = glmmTMB(logistic_AUC ~ acquisition_date_ordinal + harvested + (acquisition_date_ordinal * harvested | block_id)
        #                  , data = dat
        #                  ,family=beta_family(link="logit"))
        
        if(nrow(dat)>4*length(unique(dat$block_id))){
          auc_lm = lmer(logistic_AUC ~ 1 + time_since_harvest_start * harvested +
                          (1 + time_since_harvest_start * harvested | block_id)
                        , data = dat)
          
          return(auc_lm)
        }
        
      })
      names(lm_l) = ps_feats
      return(lm_l)
    }
    # ,cl = 'future'
    )
    names(auc_lm_list) = datasets
    # plan('sequential')
    
    fdr_lm_list = pblapply(1:length(datasets), function(i){
      
      dsi = datasets[i]
      
      lm_l = pblapply(1:length(ps_feats), function(j){
        
        varj = ps_feats[j]
        
        #print for debugging
        print(paste0(
          'Processing ',i,'-',j,', ',dsi,'-',varj
        ))
        
        dat = test_results_df_gs |> filter(dataset == dsi,
                                        var_name == varj,
                                        !is.na(fisher_discriminant_ratio))
        
        # auc_lm = glmmTMB(logistic_AUC ~ acquisition_date_ordinal + harvested + (acquisition_date_ordinal * harvested | block_id)
        #                  , data = dat
        #                  ,family=beta_family(link="logit"))
        
        if(nrow(dat)>4*length(unique(dat$block_id))){
          auc_lm = lmer(fisher_discriminant_ratio ~ 1 + time_since_harvest_start * harvested +
                          (1 + time_since_harvest_start * harvested | block_id)
                        , data = dat)
          
          return(auc_lm)
        }
        
      })
      names(lm_l) = ps_feats
      return(lm_l)
    }
    # ,cl = 'future'
    )
    names(fdr_lm_list) = datasets
    
    #extract model coefficients and other info
    
    extract_coefs_lme4 = function(L){ #L is a nested structure with lme4 models in lists according to feature in lists according to dataset
      coefs_df_l = pblapply(datasets, function(x){
        df_l = pblapply(ps_feats, function(y){
          m = L[[x]][[y]]
          if(!is.null(m)){
            s = summary(m)
            coefs_df = s$coefficients |> as.data.frame()
            coefs_df$model_terms = rownames(coefs_df)
            coefs_df$AIC = AIC(m)
            coefs_df$REML = REMLcrit(m)
            coefs_df$dataset = x
            coefs_df$var_name = y
            
            return(coefs_df)
          }
        })
        df = bind_rows(df_l)
        return(df)
      })
      df = bind_rows(coefs_df_l)
      rownames(df) = NULL
      return(df)
    }
    
    jmd_coefs = extract_coefs_lme4(jmd_lm_list) |>
      mutate(Estimate_rounded = round(Estimate, 3)) |>
      mutate(metric = 'Area Under Curve')
    auc_coefs = extract_coefs_lme4(auc_lm_list) |>
      mutate(Estimate_rounded = round(Estimate, 3)) |>
      mutate(metric = 'Jeffries-Matusita distance')
    fdr_coefs = extract_coefs_lme4(fdr_lm_list) |>
      mutate(Estimate_rounded = round(Estimate, 3)) |>
      mutate(metric = 'Fisher discriminant ratio')
    
    coefs_combined_df = bind_rows(jmd_coefs, auc_coefs, fdr_coefs)
    
    #tile plot for coefficients
    {
      ggplot(data = coefs_combined_df |> 
               filter(!model_terms %in% c('months_since_harvest', '(Intercept)')
                      ,var_name %in% feats_to_use
               )
             , aes(
               x = dataset,
               y = var_name))+
        geom_tile(aes(fill = Estimate_rounded))+
        geom_text(aes(label = Estimate_rounded))+
        scale_fill_viridis_b()+
        # facet_grid(rows= vars(model_terms), cols = vars(dataset), scales = 'free_x')+
        facet_grid(rows = vars(model_terms), cols = vars(metric), scales = 'free')+
        ggtitle('ITS model coefficients (growing season months)') +
        theme_classic()
      
      # ggplot(data = auc_coefs |> filter(model_terms %in% c('(Intercept)')), aes(x = dataset, y = var_name))+
      #   geom_tile(aes(fill = Estimate_rounded))+
      #   geom_text(aes(label = Estimate_rounded))+
      #   scale_fill_viridis_b()+
      #   facet_grid(vars(model_terms))+
      #   theme_classic()
    }
    
    
    #plot ITS models with counterfactual
    
    desired_feats = c('SR')
    desired_datasets = c('Z')
    
    plotting_df = test_results_df |>
      filter(var_name %in% desired_feats
             , dataset %in% desired_datasets) |>
      pivot_longer(cols = c('logistic_AUC', 'jeffries_matusita_dist', 'fisher_discriminant_ratio')
                   ,names_to = 'separability_metric'
                   ,values_to = 'separability_value') |>
      left_join(harvest_dates_df)
      
    
    ggplot(plotting_df
           )+ 
      geom_point(aes(x = acquisition_date, y = separability_value), alpha = 0.2)+
      facet_grid(cols = vars(separability_metric), rows = vars(block_id)) +
      geom_vline(aes(xintercept = harvest_finish_date, linetype = 'Harvest end date')) + #get harvest date for each plot
      geom_vline(aes(xintercept = harvest_start_date, linetype = 'Harvest start date'))
    
  }
  
  # models from growing season, only using pre and post harvest dates
  { 
    gs_months = c(5,6,7,8,9,10)
    test_results_df_gs = test_results_df |> 
      filter(month %in% gs_months) |> #limit to growing season months
      filter(acquisition_date < harvest_start_date|acquisition_date >= harvest_finish_date)
    
    # plan('multisession', workers = 8)
    jmd_lm_list = pblapply(1:length(datasets), function(i){
      
      dsi = datasets[i]
      
      lm_l = pblapply(1:length(ps_feats), function(j){
        
        varj = ps_feats[j]
        
        #print for debugging
        print(paste0(
          'Processing ',i,'-',j,', ',dsi,'-',varj
        ))
        
        dat = test_results_df_gs |> filter(dataset == dsi,
                                           var_name == varj,
                                           !is.na(jeffries_matusita_dist))
        
        
        # jmd_lm = glmmTMB(jmd_scaled ~ acquisition_date_ordinal + harvested + (acquisition_date_ordinal * harvested | block_id)
        #                  , data = dat
        #                  ,family=ordbeta(link="logit")) #ordered beta regression since the dataset includes zeroes and ones
        
        if(nrow(dat)>4*length(unique(dat$block_id))){
          jmd_lm = lmer(jmd_scaled ~ 1 + time_since_harvest_start * harvested +
                          (1 + time_since_harvest_start * harvested | block_id)
                        , data = dat)
          return(jmd_lm)
        }
      })
      names(lm_l) = ps_feats
      return(lm_l)
    }
    # ,cl = 'future'
    )
    names(jmd_lm_list) = datasets
    
    auc_lm_list = pblapply(1:length(datasets), function(i){
      
      dsi = datasets[i]
      
      lm_l = pblapply(1:length(ps_feats), function(j){
        
        varj = ps_feats[j]
        
        #print for debugging
        print(paste0(
          'Processing ',i,'-',j,', ',dsi,'-',varj
        ))
        
        dat = test_results_df_gs |> filter(dataset == dsi,
                                           var_name == varj,
                                           !is.na(logistic_AUC))
        
        # auc_lm = glmmTMB(logistic_AUC ~ acquisition_date_ordinal + harvested + (acquisition_date_ordinal * harvested | block_id)
        #                  , data = dat
        #                  ,family=beta_family(link="logit"))
        
        if(nrow(dat)>4*length(unique(dat$block_id))){
          auc_lm = lmer(logistic_AUC ~ 1 + time_since_harvest_start * harvested +
                          (1 + time_since_harvest_start * harvested | block_id)
                        , data = dat)
          
          return(auc_lm)
        }
        
      })
      names(lm_l) = ps_feats
      return(lm_l)
    }
    # ,cl = 'future'
    )
    names(auc_lm_list) = datasets
    # plan('sequential')
    
    fdr_lm_list = pblapply(1:length(datasets), function(i){
      
      dsi = datasets[i]
      
      lm_l = pblapply(1:length(ps_feats), function(j){
        
        varj = ps_feats[j]
        
        #print for debugging
        print(paste0(
          'Processing ',i,'-',j,', ',dsi,'-',varj
        ))
        
        dat = test_results_df_gs |> filter(dataset == dsi,
                                           var_name == varj,
                                           !is.na(fisher_discriminant_ratio))
        
        # auc_lm = glmmTMB(logistic_AUC ~ acquisition_date_ordinal + harvested + (acquisition_date_ordinal * harvested | block_id)
        #                  , data = dat
        #                  ,family=beta_family(link="logit"))
        
        if(nrow(dat)>4*length(unique(dat$block_id))){
          auc_lm = lmer(fisher_discriminant_ratio ~ 1 + time_since_harvest_start * harvested +
                          (1 + time_since_harvest_start * harvested | block_id)
                        , data = dat)
          
          return(auc_lm)
        }
        
      })
      names(lm_l) = ps_feats
      return(lm_l)
    }
    # ,cl = 'future'
    )
    names(fdr_lm_list) = datasets
    
    #extract model coefficients and other info
    
    extract_coefs_lme4 = function(L){ #L is a nested structure with lme4 models in lists according to feature in lists according to dataset
      coefs_df_l = pblapply(datasets, function(x){
        df_l = pblapply(ps_feats, function(y){
          m = L[[x]][[y]]
          if(!is.null(m)){
            s = summary(m)
            coefs_df = s$coefficients |> as.data.frame()
            coefs_df$model_terms = rownames(coefs_df)
            coefs_df$AIC = AIC(m)
            coefs_df$REML = REMLcrit(m)
            coefs_df$dataset = x
            coefs_df$var_name = y
            
            return(coefs_df)
          }
        })
        df = bind_rows(df_l)
        return(df)
      })
      df = bind_rows(coefs_df_l)
      rownames(df) = NULL
      return(df)
    }
    
    jmd_coefs = extract_coefs_lme4(jmd_lm_list) |>
      mutate(Estimate_rounded = round(Estimate, 3)) |>
      mutate(metric = 'Area Under Curve')
    auc_coefs = extract_coefs_lme4(auc_lm_list) |>
      mutate(Estimate_rounded = round(Estimate, 3)) |>
      mutate(metric = 'Jeffries-Matusita distance')
    fdr_coefs = extract_coefs_lme4(fdr_lm_list) |>
      mutate(Estimate_rounded = round(Estimate, 3)) |>
      mutate(metric = 'Fisher discriminant ratio')
    
    coefs_combined_df = bind_rows(jmd_coefs, auc_coefs, fdr_coefs)
    
    #tile plot for coefficients
    {
      ggplot(data = coefs_combined_df |> 
               filter(!model_terms %in% c('months_since_harvest', '(Intercept)')
                      ,var_name %in% feats_to_use
               )
             , aes(
               x = dataset,
               y = var_name))+
        geom_tile(aes(fill = Estimate_rounded))+
        geom_text(aes(label = Estimate_rounded))+
        scale_fill_viridis_b()+
        # facet_grid(rows= vars(model_terms), cols = vars(dataset), scales = 'free_x')+
        facet_grid(rows = vars(model_terms), cols = vars(metric), scales = 'free')+
        ggtitle('ITS model coefficients (growing season months)') +
        theme_classic()
      
      # ggplot(data = auc_coefs |> filter(model_terms %in% c('(Intercept)')), aes(x = dataset, y = var_name))+
      #   geom_tile(aes(fill = Estimate_rounded))+
      #   geom_text(aes(label = Estimate_rounded))+
      #   scale_fill_viridis_b()+
      #   facet_grid(vars(model_terms))+
      #   theme_classic()
    }
    
    
    #plot ITS models with counterfactual
    
    desired_feats = c('SR')
    desired_datasets = c('Z')
    
    plotting_df = test_results_df |>
      filter(var_name %in% desired_feats
             , dataset %in% desired_datasets) |>
      pivot_longer(cols = c('logistic_AUC', 'jeffries_matusita_dist', 'fisher_discriminant_ratio')
                   ,names_to = 'separability_metric'
                   ,values_to = 'separability_value') |>
      left_join(harvest_dates_df)
    
    
    ggplot(plotting_df
    )+ 
      geom_point(aes(x = acquisition_date, y = separability_value), alpha = 0.2)+
      facet_grid(cols = vars(separability_metric), rows = vars(block_id)) +
      geom_vline(aes(xintercept = harvest_finish_date, linetype = 'Harvest end date')) + #get harvest date for each plot
      geom_vline(aes(xintercept = harvest_start_date, linetype = 'Harvest start date'))
    
  }
  

  #----separate linear models for AUC, JM distance, and FDR for each block and to compare pre-harvest with post-harvest----
  
  ps_feats = unique(test_results_df$var_name)
  ps_feats = ps_feats[!ps_feats %in% c('VARIgreen', 'max_DN', 'max')]
  
  datasets = unique(test_results_df$dataset)
  datasets = datasets[datasets %in% c('Z', 'Non-normalized', 'SM', 'Zrobust')]
  
  #before vs after harvesting, just growing season
  {
    jmd_lm_list = pblapply(1:length(datasets), function(i){
      
      dsi = datasets[i]
      
      feat_l = pblapply(1:length(ps_feats), function(j){
        
        varj = ps_feats[j]
        
        blocks_l = pblapply(1:length(block_ids), function(k){
          
          block_id = block_ids[k]
          
          #print for debugging
          print(paste0(
            'Processing ',i,'---',j,'---',k,'; ',dsi,'---',varj,'---',block_id
          ))
          
          dat = test_results_df_gs |> filter(dataset == dsi,
                                             var_name == varj,
                                             !is.na(jeffries_matusita_dist))
          
          
          # jmd_lm = glmmTMB(jmd_scaled ~ acquisition_date_ordinal + harvested + (acquisition_date_ordinal * harvested | block_id)
          #                  , data = dat
          #                  ,family=ordbeta(link="logit")) #ordered beta regression since the dataset includes zeroes and ones
          
          if(nrow(dat)>4*length(unique(dat$block_id))){
            jmd_lm = lmer(jmd_scaled ~ 1 + time_since_harvest_start * harvested +
                            (1 + time_since_harvest_start * harvested | block_id)
                          , data = dat)
            
            return(jmd_lm)
          }
          
        })
          
          #print for debugging
          print(paste0(
            'Processing ',i,'-',j,', ',dsi,'-',varj
          ))
        
        dat = test_results_df_gs |> filter(dataset == dsi,
                                           var_name == varj,
                                           !is.na(jeffries_matusita_dist))
        
        
        # jmd_lm = glmmTMB(jmd_scaled ~ acquisition_date_ordinal + harvested + (acquisition_date_ordinal * harvested | block_id)
        #                  , data = dat
        #                  ,family=ordbeta(link="logit")) #ordered beta regression since the dataset includes zeroes and ones
        
        if(nrow(dat)>4*length(unique(dat$block_id))){
          jmd_lm = lmer(jmd_scaled ~ 1 + time_since_harvest_start * harvested +
                          (1 + time_since_harvest_start * harvested | block_id)
                        , data = dat)
          
          return(jmd_lm)
        }
      })
      names(lm_l) = ps_feats
      return(lm_l)
    }
    # ,cl = 'future'
    )
    names(jmd_lm_list) = datasets
    
    
  }
  }
#----plotting----

# nt3_r_l = pblapply(1:nrow(results_12N_T3 |> filter(dataset=='Non-normalized')), function(i){
#   d = results_12N_T3 |> filter(dataset == 'Non-normalized',
#                                block_pixel_stratum=='Total')
#   f = d$file_path[i]
#   r = rast(f)
#   return(r)
# })

id_ext = '20220901_184746_98_24a4.tif'

nn_id_ext = paste0(ps_dir,'/Non-normalized_REDO/12N_T3/',id_ext)

sample_nn_rast_preharv = rast(nn_id_ext)
sample_nn_rast_preharv = sample_nn_rast_preharv[[feats_to_use]]


layer_plots = pblapply(names(sample_nn_rast_preharv), function(x){
  
  layr = sample_nn_rast_preharv[[x]]
  
  p = ggplot()+
    geom_spatraster(data = layr)+
    facet_wrap(~lyr)+
    scale_fill_viridis_c(option='viridis', na.value = 'white', name = 'Value')+
    theme_classic()+
    # coord_equal(expand = FALSE) +
    theme(
      axis.title = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank()
    )+
    theme(legend.position = 'none')
  return(p)
  
})

combined_layers_plot = wrap_plots(layer_plots, nrow=1)
# combined_layers_plot

z_id_ext = paste0(ps_dir,'/Z_REDO/12N_T3/',id_ext)
zr_id_ext = paste0(ps_dir,'/Zrobust_REDO/12N_T3/',id_ext)

sample_z_rast_preharv = rast(z_id_ext)
sample_zr_rast_preharv = rast(zr_id_ext)

combined_ndvi = c(
  sample_nn_rast_preharv[['NDVI']],
  sample_z_rast_preharv[['NDVI']],
  sample_zr_rast_preharv[['NDVI']]
)
names(combined_ndvi) = c('Non-normalized NDVI', 'Z NDVI', 'Zr NDVI')

ggplot()+
  geom_spatraster(data = combined_ndvi)+
  facet_wrap(~lyr)+
  scale_fill_viridis_c(option='viridis', na.value = 'white', name = 'Value')+
  theme_classic()+
  # theme(legend.position = 'none')+
  theme(
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank()
  )

ndvi_plots = pblapply(names(combined_ndvi), function(x){
  
  layr = combined_ndvi[[x]]
  
  p = ggplot()+
    geom_spatraster(data = layr)+
    facet_wrap(~lyr)+
    scale_fill_viridis_c(option='viridis', na.value = 'white', name = 'Value')+
    theme_classic()+
    # coord_equal(expand = FALSE) +
    theme(
      axis.title = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank()
    )+
    theme(legend.position = 'none')
  return(p)
})
ndvi_combined_plot = wrap_plots(ndvi_plots)
ndvi_combined_plot


#plot temporal resolution

temp_res_nt3 = results_12N_T3 |>
  distinct(id, acquisition_date, acquisition_year, acquisition_month_day) |>
  mutate(month_day_dummy = as.Date(paste0('2000-',acquisition_month_day)))

ggplot(temp_res_nt3, aes(x = month_day_dummy, y = 0))+
  geom_point(alpha = 0.3, color = 'purple4')+
  scale_x_date(date_labels = "%b", date_breaks = "1 month") +
  facet_grid(rows = vars(acquisition_year))+
  theme_classic()+
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank(),
        panel.grid.major.x = element_line(color = "grey80"),
        panel.grid.minor.x = element_blank())+
  xlab('Acquisition date')

dates = unique(temp_res_nt3$acquisition_date) |> sort()
date_gaps = sapply(2:length(dates), function(i){
  d = dates[i]
  d_past = dates[i-1]
  gap = d - d_past
  return(gap)
})
mean(date_gaps)

#----plot sample images

redo_datasets = c('Z_REDO', 'Zrobust_REDO', 'Non-normalized_REDO')
redo_dirs = paste0(ps_dir,'/',redo_datasets,'/12N_T3')

redo_files_lists = lapply(redo_dirs, function(x)list.files(x,full.names = T))
names(redo_files_lists) = redo_datasets

sample_ids = c("20210626_182716_53_2235"
               , "20230123_184930_42_248f"
               # ,  "20230128_181143_14_241b"
               ,'20230128_185503_32_2473'
               ,"20230129_181440_00_2442"
               ,'20241011_192907_10_251a')


sample_timeseries_rasts = pblapply(redo_dirs, function(x){
  r_l = pblapply(sample_ids, function(y){
    path = paste0(x,'/',y,'.tif')
    r = rast(path)
    return(r)
  })
  return(r_l)
})
names(sample_timeseries_rasts) = str_remove(redo_datasets, '_REDO')

a = sample_timeseries_rasts$Z[[3]]
ggRGB(img = a, r = 6, g = 4, b = 2)


ndvi_timeseries_rasts = pblapply(redo_dirs, function(x){
  ndvi_l = pblapply(sample_ids, function(y){
    path = paste0(x,'/',y,'.tif')
    r = rast(path)
    ndvi = r[['NDVI']]
    return(ndvi)
  })
  ndvi_rast = rast(ndvi_l)
  
  acquisition_dates <- as.Date(substr(sample_ids, 1, 8), format = "%Y%m%d")
  
  names(ndvi_rast)= acquisition_dates
  return(ndvi_rast)
})
names(ndvi_timeseries_rasts) = str_remove(redo_datasets, '_REDO')

z_ndvi_plot = ggplot()+
  geom_spatraster(data = ndvi_timeseries_rasts$Z)+
  facet_wrap(~lyr)
z_ndvi_plot



#----save raw values to datafiles----

data_tables_dir = paste0(ps_dir,'/Raw_value_data_tables')
dir.check(data_tables_dir)

library(arrow)

# plan('multisession', workers = 10)

pblapply(dirs, function(x){
  
  dataset = basename(x)
  dataset_dir = paste0(data_tables_dir,'/',dataset) #create subdirectory
  dir.check(dataset_dir)
  
  block_subdirs = list.dirs(x, full.names = T)[x != list.dirs(x, full.names = T)]
  bock_subdirs = block_subdirs[str_detect(block_subdirs,'12N_T3')]
  
  pblapply(block_subdirs, function(y){
    
    block = basename(y)
    block_dir = paste0(dataset_dir, '/', block)
    dir.check(block_dir)
    
    rasters = list.files(y, full.names = T, pattern = '\\.tif$') #get list of rasters
    
    pblapply(rasters, function(z){
      
      id = file_path_sans_ext(basename(z))
      
      filename = paste0(block_dir,'/',id,'.parquet')
      print(paste('Processing',filename))
      if(!file.exists(filename)){
        
        #load raster
        r = rast(z)
        
        if(names(r)[35] == 'TVI'){r = r[[-35]]} #get rid of extra TVI layer
        
        #load lidar change raster
        lid = rast(lidar_files[[block]])
        names(lid) = 'CHM_change'
        
        #resample lidar, combine with scene data
        if(ext(r) != ext(lid) | !identical(res(r), res(lid))){
          lid = resample(lid, r)
        }
        
        r = c(r, lid)
        
        #get values
        v = values(r, dataframe = T) |> setDT()
        
        v[, X := xFromCell(r,1:ncell(r))] #get X values of cells
        
        v[, Y := yFromCell(r,1:ncell(r))] # get y values of cells
        
        v[, block_id := block] #add block_id
        
        v[, dataset := dataset] #add dataset
        
        v[, cellnum := 1:ncell(r)] #assign number to each cell
        
        v[, id := id] #add id
        
        v[, acquisition_date := as.Date(paste0( #get acquisition date
          substr(id, 1, 4), "-",  # year
          substr(id, 5, 6), "-",  # month
          substr(id, 7, 8)        # day
        ))]
        
        
        #write data table
        write_parquet(v, filename)
      }
    }
    # ,cl='future'
    )
  })
}
# ,cl='future'
)
# plan('sequential')
