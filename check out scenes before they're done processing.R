
#c

global_tables_dir = "D:/Quesnel/data/planet_scenes/Stats_logAUC_JMD_HarvestThreshold=Otsu_EqualClasses"
global_tables_files = list.files(global_tables_dir, full.names = T, pattern = '\\.arrow$')

global_tables_l = list()
bad_indices = c()

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
}

results_df1 = rbindlist(global_tables_l)
results_df = results_df1 |>
  #get scene ID
  mutate(id = basename(file_path)) |>
  # fix acquisition_date
  mutate(acquisition_date =as.Date(
           paste0(
    substr(id,1,4),'-', #year
    substr(id,5,6),'-', #month
    substr(id,7,8) #day
  )
           )) |>
  rename(var_name = v1) |>
  mutate(acquisition_date2 = acquisition_date,
         month = month(acquisition_date))

results_summ = results_df |>
  filter(block_pixel_stratum=='Total') |>
  group_by(dataset,block_id,var_name) |>
  summarize(count = n())

ggplot(data = results_summ, aes(x = block_id, y = var_name))+
  geom_point(aes(size=count))+
  facet_grid(vars(dataset))

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
) |> mutate(harvest_start_date = as.Date(harvest_start_date),
            harvest_finish_date = as.Date(harvest_finish_date))
# |>
results_df = results_df |>
  left_join(harvest_dates_df, by = 'block_id') 
# |>
#   mutate(harvest_date = as.Date(paste0(harvest_month, '-01')))
#calculate time since harvest
results_df = results_df |>
  mutate(months_since_harvest = as.numeric(acquisition_date2 - harvest_start_date)) |>
  mutate(harvested = ifelse(acquisition_date2 < harvest_start_date,0,1)) |>
  left_join(tibble(acquisition_date2 = sort(unique(results_df$acquisition_date2)),
                   acquisition_date_ordinal = as.numeric(1:length(unique(results_df$acquisition_date2)))))
#rescale JM distance to give it the range [0-1] so that beta regression can be applied to both AUC and JMD
results_df = results_df |>
  mutate(jmd_scaled = jeffries_matusita_dist/sqrt(2))

#make dataframe that only has the results from the tests (i.e. remove unnecessary rows)
test_results_df = results_df |> filter(block_pixel_stratum == 'Total') |> dplyr:: select(-block_pixel_stratum)




# ----plotting vegetation indices over time----
{
  indices = c(
    # "blue", "green", "red",
    'nir'
    # ,
    # "Hue"
    # ,
    # "GCC"
    # ,
    # "NDGR"
    # ,
    # "BI"
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
  
  subsumm = subset_df |> 
    group_by(block_id, block_pixel_stratum)|>
    summarise(n_dates = n())

  #harvest pixel values
  pixels_p = ggplot(subset_df, aes(x = acquisition_date2, y = mean)) +
    geom_point(aes(color = block_pixel_stratum))+
    # geom_line(aes(color = block_pixel_stratum))+
    geom_ribbon(aes(ymin = mean - sd, ymax = mean + sd, fill = block_pixel_stratum), alpha = 0.2, color = NA) +
    # geom_ribbon(aes(x = acquisition_date, ymin = mean-sds, ymax = mean+sds))+
    scale_x_date(
      date_breaks = "1 year",       # Major tick marks every year
      date_labels = "%Y %m",           # Year labels on major ticks
      date_minor_breaks = "1 month" # Minor tick marks every month
    )+
    facet_grid(rows = vars(block_id), cols = vars(var_name), scales = 'free') +
    # geom_vline(xintercept = as.Date('2022-07-01')) + #12L_D345
    geom_vline(aes(xintercept = harvest_start_date, linetype = 'Harvest date')) + #get harvest date for each plot
    ggtitle(unique(subset_df$dataset))+
    theme_classic()

  #logistic_AUC results
  log_p = ggplot(subset_df |> filter(block_pixel_stratum == 'Thinned') #remove redundant info (since thinning, not_thinning, and total have same stats info)
                 , aes(x = acquisition_date2, y = logistic_AUC)) +
    geom_point()+
    # geom_line()+
    # geom_point(data = subset_df|>
    #              filter(logistic_sig == 1)|>
    #              mutate(label = paste0('p<',alpha)),
    #            aes(x = acquisition_date2, y =1, color=label), shape = 8)+ #significance of model
    # scale_shape_manual(name = "Model Significance", values = c("Significant Model" = 8))+
    scale_x_date(
      date_breaks = "1 year",       # Major tick marks every year
      date_labels = "%Y %m",           # Year labels on major ticks
      date_minor_breaks = "1 month" # Minor tick marks every month
    )+
    facet_grid(rows = vars(block_id), cols = vars(var_name)
               , scales = 'free'
    ) +
    # ylim(0.3,1)+
    geom_vline(aes(xintercept = harvest_start_date, linetype = 'Harvest date')) + #get harvest date for each plot
    ggtitle(unique(subset_df$dataset))+
    theme_classic()

  #jeffries-matusita results
  jm_p = ggplot(subset_df |> filter(block_pixel_stratum == 'Thinned') #remove redundant info (since thinning, not_thinning, and total have same stats info)
                , aes(x = acquisition_date2, y = jmd_scaled)) +
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
    # ylim(0,0.3)+
    geom_vline(aes(xintercept = harvest_start_date, linetype = 'Harvest date')) + #get harvest date for each plot
    ggtitle(unique(subset_df$dataset))+
    theme_classic()

  library(patchwork)
  pixels_p+log_p+jm_p+plot_layout(guides = 'collect')
}

#get temporal distribution of scenes
ggplot()+
  geom_point(data = results_df |> filter(var_name=='NDVI',block_pixel_stratum=='Total'), aes(x = acquisition_date, y =1))+
  facet_wrap(vars(block_id,dataset))

