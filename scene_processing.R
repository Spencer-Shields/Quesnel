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
library(car)
library(nortest)
library(pROC)
library(varSel)
library(tseries)
library(lme4)
library(glmmTMB)
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

#load Otsu thresholds
otsu_df = read_csv("data/Quesnel_thinning/Otsu_change_thresholds.csv")

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
  #----extract stats from remote sensing data----
  harvest_threshold = "Otsu"
  
  data_string = paste0('Stats_logAUC_JMD_HarvestThreshold=',harvest_threshold,'_EqualClasses')
  data_filename = paste0(ps_dir,'/',data_string,'.arrow')
  
  #run if summary data file does not exist
  if(!file.exists(data_filename)){
    
    #----set up for processing----
    
    #get list of all files
    dirs = c(
      z_dir
      ,zr_dir
      ,
      sm_dir
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
    
    #----create thinning masks, get spatial sample for each study area to make thinning/non-thinning classes same size----
    
    #load lidar CHM change layers
    lid_masks_equalclasses = pblapply(lidar_files, function(x){
      
      ##load lidar raster
      # print(paste('Processing',x))
      r = rast(x)
      
      #get appropriate harvest threshold
      if(harvest_threshold == 'Otsu'){
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
    
    plan('multisession', workers = 10)
    
    global_tables_dir = paste0(ps_dir,'/',data_string)
    dir.check(global_tables_dir)
    
    # df_l = pblapply(1:10, function(i){
    future_lapply(1:length(all_files_thinning), function(i){
      
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
            
            #logistic regression model
            log_form = paste('thinning ~', y)
            log_m = glm(log_form, v, family = 'binomial')
            log_m_summ = summary(log_m)
            
            test_df$logistic_p = ifelse(is.na(log_m$coefficients[[y]]), 
                                        NA,
                                        log_m_summ$coefficients[,4][[y]])
            test_df$logistic_aic = log_m_summ$aic
            
            log_p = v |> mutate(predictions = predict(log_m, v, type = 'response'))
            log_roc = roc(thinning ~ predictions, log_p, quiet=T)
            test_df$logistic_AUC = as.numeric(auc(log_roc))
            
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
            jmd = JMdist(g = v$thinning, X = tibble(v[[y]]))
            test_df$jeffries_matusita_dist = jmd$jmdist
            
            #Bhattacharya distance
            bhat = BHATdist(g = v$thinning, X = tibble(v[[y]]))
            test_df$bhattacharya_dist = bhat$bhatdist
            
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
      all_global_tables = list.files(global_tables_dir, full.names = T, pattern = '\\.tif$')
      tables_l = pblapply(all_global_tables, function(x)setDT(read_feather(x)))
      
      results_df = rbindlist(tables_l)
      write_feather(results_df,data_filename)
      
      plan('sequential')
      }
  
  #----clean and process dataframe----
  
  # results_df = fread(data_filename)
  # 
  # results_df[, acquisition_date2 := as.Date(paste0( #get acquisition date
  #   substr(basename(file_path_sans_ext(file_path)), 1, 4), "-",  # year
  #   substr(basename(file_path_sans_ext(file_path)), 5, 6), "-",  # month
  #   substr(basename(file_path_sans_ext(file_path)), 7, 8)        # day
  # ))]
  # results_df[, julian_day := yday(acquisition_date2)]
  # results_df[, acquisition_year := year(acquisition_date2)]
  # results_df[, acquisition_month := month(acquisition_date2)]
  # results_df[, dataset := str_remove(dataset,'/')]
  # results_df[, var_name := v1]
  # 
  # 
  # 
  #   #reformat dates
  #   mutate(
  #   ) %>%
  #   mutate(julian_day = yday(acquisition_date2)
  #   ) %>%
  #   #modify variable names
  #   mutate(var_name = str_remove(v1, '_.*')) %>%
  #   # mutate(var_name = str_remove(var_name, "\\.\\.\\..*")) %>%
  #   filter(var_name != 'max_DN') %>%
  #   # fix dataset names
  #   mutate(block_type = ifelse(str_detect(block_id, 'NoChange'), 'Control', 'Thinning')) %>%
  #   mutate(dataset = str_replace(dataset, 'BlockClipped','')) %>%
  #   mutate(dataset = ifelse(str_detect(dataset, 'Z|SM'), dataset, paste0('Non-normalized', dataset))) %>%
  #   mutate(dataset = str_remove(dataset, "^_")) %>%
  #   mutate(dataset = str_remove(dataset, '/$'))
  # 
  #add harvest dates
  harvest_dates_df = tribble(
    ~block_id, ~harvest_start_date, ~harvest_finish_date,
    '12L_C5', '2022-07-31', '2022-11-28',
    '12L_D345','2022-07-06', '2022-07-20',
    '12L_C4', '2022-10-30', '2023-01-28',
    '12L_C7', '2022-09-01', '2022-11-28',
    '12L_B8C3', '2022-09-13', '2022-11-02',
    '12N_T3', '2023-01-28', '2024-02-29',
    '12N_1X1W', '2024-02-29', '2024-03-22',
    '12N_V1', '2024-03-04', '2024-04-17'
  )
  # |>
  # results_df = results_df |>
  #   left_join(harvest_dates_df, by = 'block_id') |>
  #   mutate(harvest_date = as.Date(paste0(harvest_month, '-01')))
  # #calculate time since harvest
  # results_df = results_df |>
  #   mutate(months_since_harvest = as.numeric(acquisition_date2 - harvest_date)/30.44) |>
  #   mutate(harvested = ifelse(acquisition_date2 < harvest_date,0,1)) |>
  #   left_join(tibble(acquisition_date2 = sort(unique(results_df$acquisition_date2)),
  #                    acquisition_date_ordinal = as.numeric(1:length(unique(results_df$acquisition_date2)))))
  # #rescale JM distance to give it the range [0-1] so that beta regression can be applied to both AUC and JMD
  # results_df = results_df |>
  #   mutate(jmd_scaled = jeffries_matusita_dist/sqrt(2))
  # 
  # #make dataframe that only has the results from the tests (i.e. remove unnecessary rows)
  # test_results_df = results_df |> filter(block_pixel_stratum == 'Total') |> dplyr:: select(-block_pixel_stratum)
  # 
  
  
  
  #----plotting vegetation indices over time----
  # {
  #   indices = c(
  #     # "blue", "green", "red"
  #     # ,
  #     # "Hue"
  #     # ,
  #     # "GCC"
  #     # ,
  #     # "NDGR"
  #     # ,
  #     "BI"
  #     # ,
  #     # "CI", "CRI550",
  #     # "GLI", "TVI"
  #     # , "VARIgreen"
  #   ) #vector of indices to look at
  #   
  #   blocks_ = c(
  #     "12L_C5"
  #     ,
  #     "12L_D345"
  #     ,  "12L_C4",    "12L_C7",    "12L_B8C3",  "12N_T3",    "12N_1X1W",  "12N_V1"
  #     # ,"NoChange1.1_conifer", "NoChange1.2_conifer",
  #     # "NoChange4_conifer",   "NoChange5_conifer",  
  #     # "NoChange6_conifer",   "NoChange7_conifer"
  #   )
  #   
  #   tests_ = c(
  #     'levene_p',
  #     'anderson_darling_p_thinning',
  #     'anderson_darling_p_notthinning',
  #     'anderson_darling_p_pooled',
  #     'student_t_p',
  #     'welch_t_p',
  #     'mann_whitney_p',
  #     'anova_p',
  #     'kruskal_wallis_p',
  #     'logistic_p',
  #     'logistic_aic',
  #     'logistic_AUC',
  #     'fisher_discriminant_ratio',
  #     'cohen_d',
  #     'jeffries_matusita_dist',
  #     'divergence'
  #   )
  #   
  #   datasets_ = c(
  #     'Z'
  #     # ,
  #     # 'Z_delta'
  #     # ,
  #     # 'Zrobust'
  #     # , 'Zrobust_delta'
  #     # ,
  #     # 'SM'
  #     # ,
  #     # 'SM_delta'
  #     # 'Non-normalized'
  #   )
  #   
  #   months_ = c(
  #     '01','02','03',
  #     '04',
  #     '05','06','07','08','09','10','11'
  #     ,'12'
  #   )
  #   
  #   # harvest_dates_df = tibble(
  #   #   block_id = blocks_p$BLOCKNUM[!str_detect(blocks_p$BLOCKNUM, 'NoChange')],
  #   #   harvest_month = c('2023-02'
  #   #                     ,'2022-07'
  #   #                     ,'2023-02'
  #   #                     ,'2023-02'
  #   #                     ,'2023-02'
  #   #                     ,'2023-02'
  #   #                     ,'2024-03'
  #   #                     ,'2024-03')
  #   #   ,harvest_note = c('-','-','-','-','-','partial, complete 2023-03','-','partial, complete 2024-04')
  #   # )
  #   
  #   pixel_stratum = c('Not_thinned', 'Thinned'
  #                     # , 'Total'
  #   )
  #   
  #   subset_df = results_df %>% 
  #     filter(var_name %in% indices
  #            , str_detect(block_id, paste(blocks_, collapse = "|"))
  #            , dataset %in% datasets_
  #            , month %in% months_
  #            , block_pixel_stratum %in% pixel_stratum)
  #   
  #   #harvest pixel values
  #   pixels_p = ggplot(subset_df, aes(x = acquisition_date2, y = mean)) +
  #     geom_point(aes(color = block_pixel_stratum))+
  #     # geom_line(aes(color = block_pixel_stratum))+
  #     geom_ribbon(aes(ymin = mean - sd, ymax = mean + sd, fill = block_pixel_stratum), alpha = 0.2, color = NA) +
  #     # geom_ribbon(aes(x = acquisition_date, ymin = mean-sds, ymax = mean+sds))+
  #     scale_x_date(
  #       date_breaks = "1 year",       # Major tick marks every year
  #       date_labels = "%Y %m",           # Year labels on major ticks
  #       date_minor_breaks = "1 month" # Minor tick marks every month
  #     )+
  #     facet_grid(rows = vars(block_id), cols = vars(var_name), scales = 'free') +
  #     # geom_vline(xintercept = as.Date('2022-07-01')) + #12L_D345
  #     geom_vline(aes(xintercept = harvest_date, linetype = 'Harvest date')) + #get harvest date for each plot
  #     ggtitle(unique(subset_df$dataset))+
  #     theme_classic()
  #   
  #   #logistic_AUC results
  #   log_p = ggplot(subset_df |> filter(block_pixel_stratum == 'Thinned') #remove redundant info (since thinning, not_thinning, and total have same stats info)
  #                  , aes(x = acquisition_date2, y = logistic_AUC)) +
  #     geom_point()+
  #     # geom_line()+
  #     geom_point(data = subset_df|>
  #                  filter(logistic_sig == 1)|>
  #                  mutate(label = paste0('p<',alpha)), 
  #                aes(x = acquisition_date2, y =1, color=label), shape = 8)+ #significance of model
  #     # scale_shape_manual(name = "Model Significance", values = c("Significant Model" = 8))+
  #     scale_x_date(
  #       date_breaks = "1 year",       # Major tick marks every year
  #       date_labels = "%Y %m",           # Year labels on major ticks
  #       date_minor_breaks = "1 month" # Minor tick marks every month
  #     )+
  #     facet_grid(rows = vars(block_id), cols = vars(var_name)
  #                , scales = 'free'
  #     ) +
  #     ylim(0.3,1)+
  #     geom_vline(aes(xintercept = harvest_date, linetype = 'Harvest date')) + #get harvest date for each plot
  #     ggtitle(unique(subset_df$dataset))+
  #     theme_classic()
  #   
  #   #jeffries-matusita results
  #   jm_p = ggplot(subset_df |> filter(block_pixel_stratum == 'Thinned') #remove redundant info (since thinning, not_thinning, and total have same stats info)
  #                 , aes(x = acquisition_date2, y = jmd_scaled)) +
  #     geom_point()+
  #     # geom_line()+
  #     # scale_shape_manual(name = "Model Significance", values = c("Significant Model" = 8))+
  #     scale_x_date(
  #       date_breaks = "1 year",       # Major tick marks every year
  #       date_labels = "%Y %m",           # Year labels on major ticks
  #       date_minor_breaks = "1 month" # Minor tick marks every month
  #     )+
  #     facet_grid(rows = vars(block_id), cols = vars(var_name)
  #                , scales = 'free'
  #     ) +
  #     ylim(0,0.3)+
  #     geom_vline(aes(xintercept = harvest_date, linetype = 'Harvest date')) + #get harvest date for each plot
  #     ggtitle(unique(subset_df$dataset))+
  #     theme_classic()
  #   
  #   library(patchwork)
  #   pixels_p+log_p+jm_p+plot_layout(guides = 'collect')
  # }
  # 
  #----statistical modelling to assess effect of harvesting on class separability----
  # #define features and datasets that should be tested
  # ps_feats = unique(test_results_df$var_name)
  # ps_feats = ps_feats[!ps_feats %in% c('VARIgreen', 'max_DN', 'max')]
  # 
  # datasets = unique(test_results_df$dataset)
  # datasets = datasets[datasets %in% c('Z', 'Non-normalized', 'SM', 'Zrobust')]
  # 
  # block_ids = unique(test_results_df$block_id)
  
  #----LMMs for AUC and JM distance----
  
  # plan('multisession', workers = 8)
  # jmd_lm_list = pblapply(1:length(datasets), function(i){
  #   
  #   dsi = datasets[i]
  #   
  #   lm_l = pblapply(1:length(ps_feats), function(j){
  #     
  #     varj = ps_feats[j]
  #     
  #     #print for debugging
  #     print(paste0(
  #       'Processing ',i,'-',j,', ',dsi,'-',varj
  #     ))
  #     
  #     dat = test_results_df |> filter(dataset == dsi,
  #                                     var_name == varj)
  #     
  #     
  #     # jmd_lm = glmmTMB(jmd_scaled ~ acquisition_date_ordinal + harvested + (acquisition_date_ordinal * harvested | block_id)
  #     #                  , data = dat
  #     #                  ,family=ordbeta(link="logit")) #ordered beta regression since the dataset includes zeroes and ones
  #     
  #     jmd_lm = lmer(jmd_scaled ~ 1 + months_since_harvest * harvested + 
  #                     (1 + months_since_harvest * harvested | block_id)
  #                   , data = dat)
  #     
  #     return(jmd_lm)
  #   })
  #   names(lm_l) = ps_feats
  #   return(lm_l)
  # }
  # ,cl = 'future'
  # )
  # names(jmd_lm_list) = datasets
  # 
  # auc_lm_list = pblapply(1:length(datasets), function(i){
  #   
  #   dsi = datasets[i]
  #   
  #   lm_l = pblapply(1:length(ps_feats), function(j){
  #     
  #     varj = ps_feats[j]
  #     
  #     #print for debugging
  #     print(paste0(
  #       'Processing ',i,'-',j,', ',dsi,'-',varj
  #     ))
  #     
  #     dat = test_results_df |> filter(dataset == dsi,
  #                                     var_name == varj)
  #     
  #     # auc_lm = glmmTMB(logistic_AUC ~ acquisition_date_ordinal + harvested + (acquisition_date_ordinal * harvested | block_id)
  #     #                  , data = dat
  #     #                  ,family=beta_family(link="logit"))
  #     
  #     auc_lm = lmer(logistic_AUC ~ 1 + months_since_harvest * harvested + 
  #                     (1 + months_since_harvest * harvested | block_id)
  #                   , data = dat)
  #     
  #     return(auc_lm)
  #   })
  #   names(lm_l) = ps_feats
  #   return(lm_l)
  # }
  # ,cl = 'future'
  # )
  # names(auc_lm_list) = datasets
  # plan('sequential')
  # 
  # #extract model coefficients and other info
  # 
  # extract_coefs_lme4 = function(L){ #L is a nested structure with lme4 models in lists according to feature in lists according to dataset
  #   coefs_df_l = pblapply(datasets, function(x){
  #     df_l = pblapply(ps_feats, function(y){
  #       m = L[[x]][[y]]
  #       s = summary(m)
  #       coefs_df = s$coefficients |> as.data.frame()
  #       coefs_df$model_terms = rownames(coefs_df)
  #       coefs_df$AIC = AIC(m)
  #       coefs_df$REML = REMLcrit(m)
  #       coefs_df$dataset = x
  #       coefs_df$var_name = y
  #       
  #       return(coefs_df)
  #     })
  #     df = bind_rows(df_l)
  #     return(df)
  #   })
  #   df = bind_rows(coefs_df_l)
  #   rownames(df) = NULL
  #   return(df)
  # }
  # 
  # jmd_coefs = extract_coefs_lme4(jmd_lm_list) |> 
  #   mutate(Estimate_rounded = round(Estimate, 3)) |>
  #   mutate(metric = 'Area Under Curve')
  # auc_coefs = extract_coefs_lme4(auc_lm_list) |> 
  #   mutate(Estimate_rounded = round(Estimate, 3)) |>
  #   mutate(metric = 'Jeffries-Matusita distance')
  # 
  # coefs_combined_df = bind_rows(jmd_coefs, auc_coefs)
  # 
  # #plot results
  # {
  #   ggplot(data = auc_coefs |> filter(!model_terms %in% c('months_since_harvest', '(Intercept)')), aes(x = dataset, y = var_name))+
  #     geom_tile(aes(fill = Estimate_rounded))+
  #     geom_text(aes(label = Estimate_rounded))+
  #     scale_fill_viridis_b()+
  #     facet_grid(vars(model_terms))+
  #     theme_classic()
  #   
  #   # ggplot(data = auc_coefs |> filter(model_terms %in% c('(Intercept)')), aes(x = dataset, y = var_name))+
  #   #   geom_tile(aes(fill = Estimate_rounded))+
  #   #   geom_text(aes(label = Estimate_rounded))+
  #   #   scale_fill_viridis_b()+
  #   #   facet_grid(vars(model_terms))+
  #   #   theme_classic()
  # }
  
  
  }

#----get RMSE between adjacent scenes in each block----

#----save raw values to datafiles----

data_tables_dir = paste0(ps_dir,'/Raw_value_data_tables')
dir.check(data_tables_dir)

library(arrow)

plan('multisession', workers = 4)

pblapply(dirs, function(x){
  
  dataset = basename(x)
  dataset_dir = paste0(data_tables_dir,'/',dataset) #create subdirectory
  dir.check(dataset_dir)
  
  block_subdirs = list.dirs(x, full.names = T)[x != list.dirs(x, full.names = T)]
  
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
    ,cl='future'
    )
  })
}
# ,cl='future'
)
plan('sequential')
