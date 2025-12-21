source('scene_global_stats_filtering_20251106.R')

library(ggbreak)
# library(dtplyr)
# library(sgsR)

#----plot from single dataframe----

# #function for stratifying sample
# strat_histogram = function(mraster, breaks = 10){
#   #stratify a raster using percentiles
#   h = hist(mraster,breaks=breaks,plot=F)
#   breaks_ = h$breaks
#   
#   sraster = terra::classify(mraster,breaks_,include.lowest=T)
#   names(sraster) = 'strata'
#   return(sraster)
# }

sample_size = 1500

plot_dfl = pblapply(block_ids, function(b){
  
  b_df = results_meta_df |>
    filter(!acquisition_month %in% gs_months) |>
    select(block_id, dataset, acquired, file_path) |>
    filter(block_id == b) |>
    distinct()
  
  #load canopy height change data
  lidfl = list.files(lid_mets_change_cropped_dir, full.names = T)
  lidf = lidfl[basename(file_path_sans_ext(lidfl))==b]
  lid_r = rast(lidf, lyrs='z_p95')
  
  # #get sampling points using delta-CHM data
  # samp_size = 1500
  # 
  # lid_strat = strat_histogram(mraster = lid_r)
  # samp = sample_strat()
  # ss = spatSample(x = lid_strat, size = samp_size, method = 'stratified', as.points=T)
  
  #extract data for first dataset separately in order to get points to get sample points
  datasets = b_df$dataset |> unique()
  
  {
    # d = datasets[1]
    # 
    # bd_df = b_df |>
    #   filter(dataset==d)
    # 
    # #load PlanetScope timeseries data
    # rl = lapply(bd_df$file_path, function(f){
    #   # cat('\r',f); 
    #   rast(f,lyrs='NDVI')})
    # r = rl |> rast()
    # names(r) = bd_df$acquired
    # 
    # #stack with lidar change data
    # r = c(r, lid_r)
    # 
    # #get cell values as dataframe
    # df = as.data.frame(r, xy=T) |>
    #   # setDT()|>
    #   mutate(xy = paste(x,y)) |>
    #   select(-x, -y)
    # 
    # #extract sample for plotting
    # samp_pix = df$xy |> sample(1000)
    # df = df[df$xy %in% samp_pix,]
    # 
    # #make long dataframe for plotting time series
    # df_long = df |>
    #   pivot_longer(cols = bd_df$acquired, names_to = 'acquired', values_to = 'NDVI') |>
    #   mutate(acquired = as.POSIXct(acquired, format = "%Y-%m-%dT%H:%M:%OSZ", tz = "UTC"),
    #          block_id = b,
    #          dataset = d
    #          # , z_p95_binned = case_when(
    #          #   z_p95 >= 0
    #          # )
    #   )
  }
  
  
  dsl = pblapply(datasets, function(d){
    
    bd_df = b_df |>
      filter(dataset==d)
    
    #load PlanetScope timeseries data
    rl = lapply(bd_df$file_path, function(f){
      # cat('\r',f)
      rast(f,lyrs='NDVI')})
    r = rl |> rast()
    names(r) = bd_df$acquired
    
    #stack with lidar change data
    # if(ext(lid_r) < ext(r)){lid_r = extend(lid_r,r)}
    # if(ext(r) < ext(lid_r)){r = extend(r, lid_r)}
    if(any(as.vector(ext(lid_r)) > as.vector(ext(r)))){lid_r = extend(x=lid_r,y=r)}
    if(any(as.vector(ext(lid_r)) < as.vector(ext(r)))){r = extend(r, lid_r)}
    
    r = c(r, lid_r)
    
    #get cell values as dataframe
    tic()
    df = as.data.frame(r, xy=T) |>
      # setDT()|>
      mutate(xy = paste(x,y)) |>
      select(-x, -y)
    
    #extract sample for plotting
    set.seed(123)
    samp_pix = df$xy |> sample(sample_size)
    df = df[df$xy %in% samp_pix,]
    
    #make long dataframe for plotting time series
    df_long = df |>
      pivot_longer(cols = bd_df$acquired, names_to = 'acquired', values_to = 'NDVI') |>
      mutate(acquired = as.POSIXct(acquired, format = "%Y-%m-%dT%H:%M:%OSZ", tz = "UTC"),
             block_id = b,
             dataset = d
             # , z_p95_binned = case_when(
             #   z_p95 >= 0
             # )
      )
    toc()
    
    return(df_long)
  })
  ds_df = rbindlist(dsl)
  
  return(ds_df)
  
})
plot_df = rbindlist(plot_dfl) |>
  mutate(dataset = case_match(dataset,
                              'Non-normalized' ~ 'NN',
                              'Zrobust' ~ 'Z',
                              .default = dataset
  )
  )


#generate facet_grid plot for each dataset (since scale_break_x doesn't play well with facetting by column)

col_plots_l = pblapply(unique(plot_df$dataset), function(d){

p = ggplot(plot_df |>
             filter(dataset == d)|>
             mutate(dataset = case_match(dataset,
                                         'Non-normalized' ~ 'NN',
                                         'Zrobust' ~ 'Z',
                                         .default = dataset
             )
             )
)+
  geom_line(aes(x = acquired, y = NDVI, color = z_p95, group = xy), alpha=0.1)+
  scale_color_viridis_c(option='turbo', direction= -1)+
  #gaps to show only growing season
  scale_x_break(c(as.POSIXct('2020-11-01 00:00:00 UTC'),as.POSIXct('2021-04-30 23:59:59 UTC')))+
  scale_x_break(c(as.POSIXct('2021-11-01 00:00:00 UTC'),as.POSIXct('2022-04-30 23:59:59 UTC')))+
  scale_x_break(c(as.POSIXct('2022-11-01 00:00:00 UTC'),as.POSIXct('2023-04-30 23:59:59 UTC')))+
  scale_x_break(c(as.POSIXct('2023-11-01 00:00:00 UTC'),as.POSIXct('2024-04-30 23:59:59 UTC')))+
  #format x axis date/time info
  scale_x_datetime(
    date_labels = "%b %Y",       # e.g., "Jan 2021"
    # date_breaks = "3 months",     # choose interval for ticks
    limits = c(as.POSIXct('2020-05-01 00:00:00 UTC'),as.POSIXct('2024-10-31 23:59:59 UTC'))
  )+
  #show vertical lines for start and end of harvesting
  geom_vline(data = harvest_dates_df, aes(xintercept = harvest_finish_date, linetype = 'Harvest end date'))+
  geom_vline(data = harvest_dates_df ,aes(xintercept = harvest_start_date, linetype = 'Harvest start date'))+
  scale_linetype_manual(
    name = "Harvest date",
    values = c(
      "Harvest start date" = "dotted",
      "Harvest end date" = "dashed"
    )
  )+
  theme_minimal()+
  theme(
    #remove x axis text on top of plot (put there by ggbreaks)
    axis.text.x.top = element_blank(),
    axis.ticks.x.top = element_blank(),
    #format lower x axis text
    axis.text.x = element_text(angle=90)
  ) +
  facet_grid(rows = vars(block_id)
             # , cols = vars(dataset)
             , scales = 'free')

return(p)
})


#----make individual plots----
# 
# plot_list = pblapply(block_ids[block_ids=='12L_D345'], function(b){
#   b_df = results_meta_df |>
#     select(block_id, dataset, acquired, file_path) |>
#     filter(block_id == b) |>
#     distinct()
#   
#   datasets = b_df$dataset |> unique()
#   
#   dsl = pblapply(datasets, function(d){
#     bd_df = b_df |>
#       filter(dataset==d)
#     
#     lidfl = list.files(lid_mets_change_cropped_dir, full.names = T)
#     lidf = lidfl[basename(file_path_sans_ext(lidfl))==b]
#     lid_r = rast(lidf, lyrs='z_p95')
#     
#     #load PlanetScope timeseries data
#     rl = lapply(bd_df$file_path, function(f){
#       # cat('\r',f); 
#       rast(f,lyrs='NDVI')})
#     r = rl |> rast()
#     names(r) = bd_df$acquired
#     
#     #stack with lidar change data
#     r = c(r, lid_r)
#     
#     #get cell values as dataframe
#     df = as.data.frame(r, xy=T) |>
#       mutate(xy = paste(x,y)) |>
#       select(-x, -y)
#     
#     #extract sample for plotting
#     samp_pix = df$xy |> sample(1000)
#     df = df[df$xy %in% samp_pix,]
#     
#     #make long dataframe for plotting time series
#     df_long = df |>
#       pivot_longer(cols = bd_df$acquired, names_to = 'acquired', values_to = 'NDVI') |>
#       mutate(acquired = as.POSIXct(acquired, format = "%Y-%m-%dT%H:%M:%OSZ", tz = "UTC")
#              # , z_p95_binned = case_when(
#              #   z_p95 >= 0
#              # )
#       )
#     
#     #plot time series
#     p = ggplot(df_long)+
#       geom_line(aes(x = acquired, y = NDVI, color = z_p95, group = xy), alpha=0.1)+
#       scale_color_viridis_c(option='turbo', direction= -1)+
#       #gaps to show only growing season
#       scale_x_break(c(as.POSIXct('2020-11-01 00:00:00 UTC'),as.POSIXct('2021-04-30 23:59:59 UTC')))+
#       scale_x_break(c(as.POSIXct('2021-11-01 00:00:00 UTC'),as.POSIXct('2022-04-30 23:59:59 UTC')))+
#       scale_x_break(c(as.POSIXct('2022-11-01 00:00:00 UTC'),as.POSIXct('2023-04-30 23:59:59 UTC')))+
#       scale_x_break(c(as.POSIXct('2023-11-01 00:00:00 UTC'),as.POSIXct('2024-04-30 23:59:59 UTC')))+
#       #format x axis date/time info
#       scale_x_datetime(
#         date_labels = "%b %Y",       # e.g., "Jan 2021"
#         # date_breaks = "3 months",     # choose interval for ticks
#         limits = c(as.POSIXct('2020-05-01 00:00:00 UTC'),as.POSIXct('2024-10-31 23:59:59 UTC'))
#       )+
#       #show vertical lines for start and end of harvesting
#       geom_vline(data = harvest_dates_df[harvest_dates_df$block_id ==b,], aes(xintercept = harvest_finish_date, linetype = 'Harvest end date'))+
#       geom_vline(data = harvest_dates_df[harvest_dates_df$block_id ==b,] ,aes(xintercept = harvest_start_date, linetype = 'Harvest start date'))+
#       scale_linetype_manual(
#         name = "Harvest date",
#         values = c(
#           "Harvest start date" = "dotted",
#           "Harvest end date" = "dashed"
#         )
#       )+
#       theme_minimal()+
#       theme(
#         #remove x axis text on top of plot (put there by ggbreaks)
#         axis.text.x.top = element_blank(),
#         axis.ticks.x.top = element_blank(),
#         #format lower x axis text
#         axis.text.x = element_text(angle=90)
#       )
#     
#     return(p)
#   })
#   names(dsl) = datasets
#   return(dsl)
# })
# names(plot_list) = block_ids

