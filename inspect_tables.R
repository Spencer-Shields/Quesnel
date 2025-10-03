tbls_dir = "D:/Quesnel/data/planet_scenes/Stats_percentiles_logAUC_JMD_FDR_HarvestThreshold=Otsu3m_EqualClasses_NoQualChecks_20250910"

tbls_files_l = list.files(tbls_dir, full.names = T)

df_l = pblapply(tbls_files_l, read_feather)

d = rbindlist(df_l) |>
  #get scene ID
  mutate(id = basename(file_path)) |>
  #fix acquisition_date
  mutate(acquisition_date = as.Date(paste0(
    substr(id,1,4),'-', #year
    substr(id,5,6),'-', #month
    substr(id,7,8) #day
  ))) |>
  rename(var_name = v1) |>
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
  mutate(jmd_scaled = jeffries_matusita_dist/sqrt(2))


#visualize jmd
ggplot(data = d |> filter(var_name=='Hue'), aes(x = acquisition_date, y = jmd_scaled)) +
  geom_point(alpha=0.05)+
  facet_grid(cols=vars(dataset), rows=vars(block_id))

#visualise distribution of pixel values over time
ggplot(data = d |> filter(var_name=='NDVI', block_pixel_stratum!='Total'), aes(x = acquisition_date)) +
  geom_point(aes(y = median, color = block_pixel_stratum), alpha = 0.1)+
  geom_ribbon(aes(ymin = X25., ymax=X75., fill = block_pixel_stratum), alpha=0.1)+
  facet_grid(cols=vars(dataset), rows=vars(block_id), scales='free')



#visualize the percent NA values with each stat
d_summ = d |>
  group_by(block_id
           # ,dataset
           # ,v1
           ) |>
  summarise(
    percent_logAUC_na = 100*sum(is.na(logistic_AUC)/n())
    , percent_JMD_na = 100*sum(is.na(jeffries_matusita_dist)/n())
    , percent_FDR_na = 100*sum(is.na(fisher_discriminant_ratio)/n())
    )

ggplot(d_summ |>
         pivot_longer(cols = c('percent_logAUC_na', 'percent_JMD_na', 'percent_FDR_na'), names_to = 'stat', values_to = 'perc')) +
  geom_col(aes(x = block_id, y = perc, fill = stat))
