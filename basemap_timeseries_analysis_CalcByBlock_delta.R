library(tidyverse)
library(terra)
library(sf)
source('helper_functions.R')
library(future)
library(future.apply)
library(tools)
library(data.table)
library(pbapply)
library(RStoolbox)

#----load directories----

bm_dir = 'data/planet_basemaps/global_monthly'
subdir_strings = c(
  'nonorm-blocks_delta',
  'Z-scoreCalcByBlock_delta',
  'SoftmaxCalcByBlock_delta'
)

dirs <- paste0(bm_dir,'/', subdir_strings)

# List and combine .tif files from all directories
all_files <- unlist(sapply(dirs, list.files, pattern = "\\.tif$", full.names = TRUE, recursive = TRUE))

#get vector of block names
blocks = st_read('data/Quesnel_thinning/12l_12n_bdy.gpkg')
no_change = st_read('Arc project/NoChangeStands.shp')
no_change = no_change %>% rename(geom = geometry)
blocks = rbind(blocks, no_change)
block_ids = blocks$BLOCKNUM

#all possible dates in timeseries
start_date <- as.Date("2021-01-01")
end_date <- as.Date("2024-12-31")
date_vector <- seq.Date(from = start_date, to = end_date, by = "month")
date_vector <- format(date_vector, "%Y-%m")
date_vector = as.character(date_vector)
date_vector = str_replace_all(date_vector,'-','_')

#helper functions
# find_substring <- function(filepath, v) {
#   match <- v[grep(paste(v, collapse = "|"), filepath)]
#   if (length(match) > 0) return(match) else return(NA)  # Return NA if no match is found
# }

find_substring = function(filepath, v){
  g = sapply(v, function(x)str_detect(filepath, x))
  if(!T %in% g){
    return(NA)
  } else {
    i = which(g)
    a = v[i]
    return(a)
  }
}

#----Turn each raster into a dataframe with global statistics----

data_string = 'GlobalStats_ByBlock_delta'
data_filename = paste0(bm_dir,'/',data_string,'.csv')

if(!file.exists(data_filename)){
  
  # plan('multisession', workers = 8)
  
  df_l = pblapply(all_files, function(x){
    
    print(paste('Processing', x))
    
    r = rast(x)
    d = global(r, fun = c("max", "min", "mean", "sum", "range", "rms", "sd", "std", "isNA", "notNA"), na.rm=T)
    
    #get dataset
    dataset = find_substring(x,subdir_strings)
    d[['dataset']] = dataset
    
    #get block
    block = find_substring(x,block_ids)
    d[['block_id']] = block
    
    #get date
    date_ = find_substring(x, date_vector)
    d[['acquisition_date']] = date_
    
    #filepath
    d[['path']] = x
    
    #rownames
    d[['v1']] = rownames(d)
    
    return(d)
  })
  
  results_df = bind_rows(df_l)
  
  write.csv(results_df, data_filename)
  
  # plan('sequential')
}




#----process combined data----
results_df = read_csv(data_filename) %>%
  mutate(
    year = str_remove(acquisition_date, '_.*'), 
    month = str_remove(acquisition_date, '.*_')
  ) %>%
  mutate(acquisition_date2 = as.Date(paste0(year,'-',month,'-01'))
  ) %>%
  mutate(julian_day = yday(acquisition_date2)
  ) %>%
  mutate(var_name = str_remove(v1, '_.*')) %>%
  # mutate(var_name = str_remove(var_name, "\\.\\.\\..*")) %>%
  filter(var_name != 'max_DN')

var_summ = results_df %>%
  group_by(var_name, dataset) %>%
  summarise(count = n())

#---- correlation matrices----

# var_df = results_df %>% pivot_wider(names_from = var_name, values_from = mean)

#----plotting vegetation indices over time----

{
  indices = c(
    # "blue", "green", "red"
    # ,
    "Hue"
    # ,
    # "GCC"
    # ,       
    # "NDGR"
    # ,
    # "BI",        
    # "CI", "CRI550",   
    # "GLI", "TVI"
    # , "VARIgreen"
  ) #vector of indices to look at
  
  blocks_ = c(
    "12L_C5",
    "12L_D345"
    ,  "12L_C4",    "12L_C7",    "12L_B8C3",  "12N_T3",    "12N_1X1W",  "12N_V1"
    ,"NoChange1"
    ,'NoChange2'
    ,'NoChange3')
  
  datasets_ = c(
    'Z-scoreCalcByBlock_delta'
    , 'SoftmaxCalcByBlock_delta', 'nonorm-blocks_delta'
  )
  
  months_ = c(
    # '01','02','03',
    # '04',
    '05','06','07','08','09','10','11'
    # ,'12'
  )
  
  harvest_dates_df = tibble(
    block = blocks_[!str_detect(blocks_, 'NoChange')],
    harvest_month = c('2023-02'
                      ,'2022-07'
                      ,'2023-02'
                      ,'2023-02'
                      ,'2023-02'
                      ,'2023-02'
                      ,'2024-03'
                      ,'2024-03')
    ,note = c('-','-','-','-','-','partial, complete 2023-03','-','partial, complete 2024-04')
  )
  
  
  
  subset_df = results_df %>% 
    filter(var_name %in% indices
           , str_detect(block_id, paste(blocks_, collapse = "|"))
           , dataset %in% datasets_ 
           , month %in% months_)
  
  ggplot(subset_df, aes(x = acquisition_date2, y = rms, color = block_id)) +
    geom_point()+
    geom_line()+
    # geom_ribbon(aes(x = acquisition_date, ymin = mean-sds, ymax = mean+sds))+
    scale_x_date(
      date_breaks = "1 year",       # Major tick marks every year
      date_labels = "%Y %m",           # Year labels on major ticks
      date_minor_breaks = "1 month" # Minor tick marks every month
    )+
    facet_grid(rows = vars(dataset), cols = vars(var_name), scales = 'free') +
    geom_vline(xintercept = as.Date('2022-07-01')) +
    ggtitle('Change from previous image')
}

#----Plot rasters to examine----
# library(patchwork)
# 
# zl = list.files(paste0(dirs[2],'/12L_D345'), full.names = T)
# names(zl) = basename(file_path_sans_ext(zl))
# 
# plot(rast(zl[['global_monthly_2022_06_mosaic.tif_12L_D345_delta']])[[7]],main = '2022-06') 
# plot(rast(zl[['global_monthly_2022_07_mosaic.tif_12L_D345_delta']])[[7]],main= '2022-07')
# plot(rast(zl[['global_monthly_2022_04_mosaic.tif_12L_D345_delta']])[[7]],main = '2022-04')
# 
# p3+p1+p2

wanted_dates = c('2022-09')


