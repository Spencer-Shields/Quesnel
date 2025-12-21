source('scene_SeparabilityAnalysis_20251103.R')

#----plot pixel values before and after thinning----
#setup data
{
  feats_to_plot = c(
    'NDVI'
    # ,
    # 'VARIgreen'
  )
  
  blocks_to_plot = c(
    # '12N_T3'
    # ,
    # '12L_D345'
    # ,
    # '12N_T3'
    # ,
    '12L_C7'
  )
  
  subset_df = results_meta_df |>
    filter(
      block_pixel_stratum != 'Total'
      ,var_name %in% feats_to_plot
      ,block_id %in% blocks_to_plot
    ) |>
    mutate(dataset = case_match(dataset, 
                                'Non-normalized' ~ 'NN', 
                                'Zrobust' ~ 'Zr',
                                .default = dataset)) |>
    mutate(block_pixel_stratum = str_replace_all(block_pixel_stratum, '_', ' '))|>
    mutate(block_pixel_stratum = factor(block_pixel_stratum, levels = c("Not thinned", "Thinned")  # desired order
    ))
  
  subset_df$acquired = ymd_hms(subset_df$acquired)
  
  harvest_dates_df <- harvest_dates_df |>
    mutate(
      harvest_start_date = as.POSIXct(harvest_start_date),
      harvest_finish_date = as.POSIXct(harvest_finish_date)
    )
  
  library(ggbreak)
  
  pixel_strata_color_map = c('Not thinned' = "#00BFC4", 'Thinned' = "#F8766D")
}

#growing season plot
{
  pix_plot_gs = ggplot(subset_df |> filter(month %in% gs_months))+
    #show 5th to 95th percentiles for pixel strata
    # geom_ribbon(aes(x = acquired, ymin = X5., ymax = X95., fill = block_pixel_stratum, colour = block_pixel_stratum), alpha = 0.1)+
    #show interquartile range of pixel strata
    geom_ribbon(aes(x = acquired, ymin = X25., ymax = X75., fill = block_pixel_stratum, color=block_pixel_stratum), alpha = 0.3)+
    #show median pixel values for the different strata
    geom_point(aes(x = acquired, y = median, color=block_pixel_stratum),alpha=0.5)+
    #format pixel strata
    # Custom colors and legend title
    scale_color_manual(
      name = "Pixel stratum"
      ,values = pixel_strata_color_map
    ) +
    scale_fill_manual(
      name = "Pixel stratum"
      ,values = pixel_strata_color_map
    ) +
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
    geom_vline(data = harvest_dates_df[harvest_dates_df$block_id %in% unique(subset_df$block_id),], aes(xintercept = harvest_finish_date, linetype = 'Harvest end date'))+
    geom_vline(data = harvest_dates_df[harvest_dates_df$block_id %in% unique(subset_df$block_id),] ,aes(xintercept = harvest_start_date, linetype = 'Harvest start date'))+
    scale_linetype_manual(
      name = "Harvest date",
      values = c(
        "Harvest start date" = "dotted",
        "Harvest end date" = "dashed"
      )
    )+
    #facet by dataset (NOTE: using ncols messes up the broken x axis)
    facet_grid(rows=vars(dataset)
               # , cols = vars(var_name)
               , scales='free_y')+
    labs(x = 'Acquisition date', y = unique(subset_df$var_name))+
    # ggtitle(unique(subset_df$block_id))+
    theme_classic()+
    theme(
      #remove x axis text on top of plot (put there by ggbreaks)
      axis.text.x.top = element_blank(),
      axis.ticks.x.top = element_blank(),
      #format lower x axis text
      axis.text.x = element_text(angle=90)
    )
  # pix_plot_gs
}

#full timeseries plot
{
  pix_plot_full = ggplot(subset_df)+
    #show 5th to 95th percentiles for pixel strata
    # geom_ribbon(aes(x = acquired, ymin = X5., ymax = X95., fill = block_pixel_stratum, colour = block_pixel_stratum), alpha = 0.1)+
    #show interquartile range of pixel strata
    geom_ribbon(aes(x = acquired, ymin = X25., ymax = X75., fill = block_pixel_stratum, color=block_pixel_stratum), alpha = 0.3)+
    #show median pixel values for the different strata
    geom_point(aes(x = acquired, y = mean, color=block_pixel_stratum),alpha=0.3)+
    #format pixel strata
    # Custom colors and legend title
    scale_color_manual(
      name = "Pixel stratum"
      ,values = pixel_strata_color_map
    ) +
    scale_fill_manual(
      name = "Pixel stratum"
      ,values = pixel_strata_color_map
    ) +
    # #gaps to show only growing season
    # scale_x_break(c(as.POSIXct('2020-11-01 00:00:00 UTC'),as.POSIXct('2021-04-30 23:59:59 UTC')))+
    # scale_x_break(c(as.POSIXct('2021-11-01 00:00:00 UTC'),as.POSIXct('2022-04-30 23:59:59 UTC')))+
    # scale_x_break(c(as.POSIXct('2022-11-01 00:00:00 UTC'),as.POSIXct('2023-04-30 23:59:59 UTC')))+
    # scale_x_break(c(as.POSIXct('2023-11-01 00:00:00 UTC'),as.POSIXct('2024-04-30 23:59:59 UTC')))+
    #format x axis date/time info
    scale_x_datetime(
      date_labels = "%b %Y",       # e.g., "Jan 2021"
      date_breaks = "3 months",     # choose interval for ticks
      limits = c(as.POSIXct('2020-01-01 00:00:00 UTC'),as.POSIXct('2024-10-31 23:59:59 UTC'))
    )+
    #show vertical lines for start and end of harvesting
    geom_vline(data = harvest_dates_df[harvest_dates_df$block_id %in% unique(subset_df$block_id),], aes(xintercept = harvest_finish_date, linetype = 'Harvest end date'))+
    geom_vline(data = harvest_dates_df[harvest_dates_df$block_id %in% unique(subset_df$block_id),] ,aes(xintercept = harvest_start_date, linetype = 'Harvest start date'))+
    scale_linetype_manual(
      name = "Harvest date",
      values = c(
        "Harvest start date" = "dotted",
        "Harvest end date" = "dashed"
      )
    )+
    #facet by dataset (NOTE: using ncols messes up the broken x axis)
    facet_grid(rows=vars(dataset)
               # , cols = vars(var_name)
               , scales='free_y')+
    labs(x = 'Acquisition date', y = unique(subset_df$var_name))+
    # ggtitle(unique(subset_df$block_id))+
    theme_classic()+
    theme(
      #remove x axis text on top of plot (put there by ggbreaks)
      axis.text.x.top = element_blank(),
      axis.ticks.x.top = element_blank(),
      #format lower x axis text
      axis.text.x = element_text(angle=90)
      #remove legend
      # ,legend.position = 'none'
    )
  pix_plot_full
}

#combined plot
{
  library(patchwork)
  (pix_plot_full+ggtitle('A')) + 
    (pix_plot_gs+ggtitle('B')) + 
    plot_layout(guides='collect', axis_titles = 'collect'
                , ncol=1
    )
}

#----plot jeffries matusita distance timeseries----


#jeffries matusita scatter full timeseries
{
  subset_df_jm = subset_df |> 
    select(acquired,dataset,var_name,block_id,jmd_scaled,harvest_start_date,harvest_finish_date) |> 
    distinct()
  
  
  jm_plot_full = ggplot(subset_df_jm)+
    #show 5th to 95th percentiles for pixel strata
    # geom_ribbon(aes(x = acquired, ymin = X5., ymax = X95., fill = block_pixel_stratum, colour = block_pixel_stratum), alpha = 0.1)+
    #show interquartile range of pixel strata
    # geom_ribbon(aes(x = acquired, ymin = X25., ymax = X75., fill = block_pixel_stratum, color=block_pixel_stratum), alpha = 0.3)+
    #show median pixel values for the different strata
    geom_point(aes(x = acquired, y = jmd_scaled),alpha=0.3)+
    #format pixel strata
    # Custom colors and legend title
    # scale_color_manual(
    #   name = "Pixel stratum"
    #   ,values = pixel_strata_color_map
    # ) +
    # #gaps to show only growing season
    # scale_x_break(c(as.POSIXct('2020-11-01 00:00:00 UTC'),as.POSIXct('2021-04-30 23:59:59 UTC')))+
    # scale_x_break(c(as.POSIXct('2021-11-01 00:00:00 UTC'),as.POSIXct('2022-04-30 23:59:59 UTC')))+
    # scale_x_break(c(as.POSIXct('2022-11-01 00:00:00 UTC'),as.POSIXct('2023-04-30 23:59:59 UTC')))+
    # scale_x_break(c(as.POSIXct('2023-11-01 00:00:00 UTC'),as.POSIXct('2024-04-30 23:59:59 UTC')))+
    #format x axis date/time info
    scale_x_datetime(
      date_labels = "%b %Y",       # e.g., "Jan 2021"
      date_breaks = "3 months",     # choose interval for ticks
      limits = c(as.POSIXct('2020-01-01 00:00:00 UTC'),as.POSIXct('2024-10-31 23:59:59 UTC'))
    )+
    #show vertical lines for start and end of harvesting
    geom_vline(data = harvest_dates_df[harvest_dates_df$block_id %in% unique(subset_df$block_id),], aes(xintercept = harvest_finish_date, linetype = 'Harvest end date'))+
    geom_vline(data = harvest_dates_df[harvest_dates_df$block_id %in% unique(subset_df$block_id),] ,aes(xintercept = harvest_start_date, linetype = 'Harvest start date'))+
    scale_linetype_manual(
      name = "Harvest date",
      values = c(
        "Harvest start date" = "dotted",
        "Harvest end date" = "dashed"
      )
    )+
    #facet by dataset (NOTE: using ncols messes up the broken x axis)
    facet_grid(rows=vars(dataset)
               # , cols = vars(var_name)
               , scales='free_y')+
    labs(x = 'Acquisition date', y = unique(subset_df$var_name))+
    # ggtitle(unique(subset_df$block_id))+
    theme_classic()+
    theme(
      #remove x axis text on top of plot (put there by ggbreaks)
      axis.text.x.top = element_blank(),
      axis.ticks.x.top = element_blank(),
      #format lower x axis text
      axis.text.x = element_text(angle=90)
      #remove legend
      # ,legend.position = 'none'
    )
  jm_plot_full
}

#set up data for all blocks and timeseries
{
  jm_allblocks_allds_df = results_meta_df |>
    filter(var_name=='Hue',
           block_pixel_stratum=='Total'
           # , block_id %in% c('12L_D345')
    ) |>
    # mutate(block_pixel_stratum = str_replace_all(block_pixel_stratum, '_', ' '))|>
    # mutate(block_pixel_stratum = factor(block_pixel_stratum, levels = c("Not thinned", "Thinned")  # desired order
    # )) |>
   mutate( time_since_harvest = ifelse(acquisition_date < harvest_start_date,
             as.numeric(acquisition_date - harvest_start_date),
             as.numeric(acquisition_date - harvest_finish_date)))
  
  jm_allblocks_allds_df$acquired = ymd_hms(jm_allblocks_allds_df$acquired)
  
  #get model predictions
  {
    
    #fit predictions
    its_preds_df = lapply(unique(jm_allblocks_allds_df$dataset), function(d){
      mod = jmd_lm_list_full[[d]][[unique(jm_allblocks_allds_df$var_name)]]
      dat = jm_allblocks_allds_df |> 
        filter(dataset==d) |>
        select(harvested, time_since_harvest, block_id, jmd_scaled, id)
      pred = predict(mod,dat)
      t = tibble(
        time_since_harvest = dat$time_since_harvest,
        block_id = dat$block_id,
        dataset = d,
        jmd_scaled = dat$jmd_scaled,
        jmd_predicted = pred
      )
      return(t)
    }) |> 
      bind_rows() |>
      left_join(jm_allblocks_allds_df) |>
      select(acquired, block_id, dataset, jmd_predicted, time_since_harvest) |>
      mutate(acquired = as.POSIXct(acquired))|>
      mutate(dataset = case_match(dataset, 
                                  'Non-normalized' ~ 'NN', 
                                  'Zrobust' ~ 'Zr',
                                  .default = dataset))
  }
  
}

#facet grid all blocks and full timeseries
{
  ggplot(jm_allblocks_allds_df |>
           mutate(dataset = case_match(dataset, 
                                       'Non-normalized' ~ 'NN', 
                                       'Zrobust' ~ 'Zr',
                                       .default = dataset)))+
    geom_point(aes(x = acquired, y = jmd_scaled), alpha=0.1)+
    scale_x_datetime(
      date_labels = "%b %Y",       # e.g., "Jan 2021"
      date_breaks = "6 months",     # choose interval for ticks
      limits = c(as.POSIXct('2020-01-01 00:00:00 UTC'),as.POSIXct('2024-10-31 23:59:59 UTC'))
    )+
    #show vertical lines for start and end of harvesting
    geom_vline(data = harvest_dates_df[harvest_dates_df$block_id %in% unique(jm_allblocks_allds_df$block_id),], aes(xintercept = as.POSIXct(harvest_finish_date), linetype = 'Harvest end date'),color='red')+
    geom_vline(data = harvest_dates_df[harvest_dates_df$block_id %in% unique(jm_allblocks_allds_df$block_id),] ,aes(xintercept = as.POSIXct(harvest_start_date), linetype = 'Harvest start date'),color='red')+
    scale_linetype_manual(
      # name = "Harvest date",
      values = c(
        "Harvest start date" = "dotted",
        "Harvest end date" = "dashed"
      )
    )+
    facet_grid(cols = vars(dataset), rows = vars(block_id))+
    labs(x = 'Acquisition date', y = paste0(unique(jm_allblocks_allds_df$var_name), ' scaled JMD'))+
    # ggtitle(unique(jm_allblocks_allds_df$block_id))+
    theme_classic()+
    theme(
      #remove x axis text on top of plot (put there by ggbreaks)
      axis.text.x.top = element_blank(),
      axis.ticks.x.top = element_blank(),
      #format lower x axis text
      axis.text.x = element_text(angle=90)
      #remove legend
      # ,legend.position = 'none'
    )+
    #add regression lines
    geom_line(data=its_preds_df, aes(x = acquired, y = jmd_predicted, color='ITS model fit'), linewidth=1)+
    scale_color_manual(values = c('ITS model fit' = 'skyblue3'), name = NULL)
}

#facet wrap all datasets for one block
{
  jmd_ts_wrap_plot = ggplot(jm_allblocks_allds_df |>
                              mutate(dataset = case_match(dataset, 
                                                          'Non-normalized' ~ 'NN', 
                                                          'Zrobust' ~ 'Zr',
                                                          .default = dataset))
  )+
    geom_point(aes(x = acquired, y = jmd_scaled), alpha = 0.3)+
    scale_x_datetime(
      date_labels = "%b %Y",       # e.g., "Jan 2021"
      date_breaks = "6 months",     # choose interval for ticks
      limits = c(as.POSIXct('2020-01-01 00:00:00 UTC'),as.POSIXct('2024-10-31 23:59:59 UTC'))
    )+
    #show vertical lines for start and end of harvesting
    geom_vline(data = harvest_dates_df[harvest_dates_df$block_id %in% unique(jm_allblocks_allds_df$block_id),], aes(xintercept = as.POSIXct(harvest_finish_date), linetype = 'Harvest end date'),color='red')+
    geom_vline(data = harvest_dates_df[harvest_dates_df$block_id %in% unique(jm_allblocks_allds_df$block_id),] ,aes(xintercept = as.POSIXct(harvest_start_date), linetype = 'Harvest start date'),color='red')+
    scale_linetype_manual(
      # name = "Harvest date",
      name = NULL,
      values = c(
        "Harvest start date" = "dotted",
        "Harvest end date" = "dashed"
      )
    )+
    facet_wrap(vars(dataset),ncol=1)+
    labs(x = 'Acquisition date', y = paste0(unique(jm_allblocks_allds_df$var_name), ' scaled JMD'))+
    # ggtitle(unique(jm_allblocks_allds_df$block_id))+
    theme_classic()+
    theme(
      #remove x axis text on top of plot (put there by ggbreaks)
      axis.text.x.top = element_blank(),
      axis.ticks.x.top = element_blank(),
      #format lower x axis text
      axis.text.x = element_text(angle=90)
      #remove legend
      # ,legend.position = 'none'
    )
  
  #get regresson predictions for plotting:
  its_preds_df = lapply(unique(jm_allblocks_allds_df$dataset), function(d){
    mod = jmd_lm_list_full[[d]][[unique(jm_allblocks_allds_df$var_name)]]
    dat = jm_allblocks_allds_df |> 
      filter(dataset==d) |>
      select(harvested, time_since_harvest, block_id, jmd_scaled, id)
    pred = predict(mod,dat)
    t = tibble(
      time_since_harvest = dat$time_since_harvest,
      block_id = dat$block_id,
      dataset = d,
      jmd_scaled = dat$jmd_scaled,
      jmd_predicted = pred
    )
    return(t)
  }) |> 
    bind_rows() |>
    left_join(jm_allblocks_allds_df) |>
    select(acquired, block_id, dataset, jmd_predicted, time_since_harvest) |>
    mutate(acquired = as.POSIXct(acquired))|>
    mutate(dataset = case_match(dataset, 
                                'Non-normalized' ~ 'NN', 
                                'Zrobust' ~ 'Zr',
                                .default = dataset))
  
  #plot regression line on top of scaled JMD values
  jmd_ts_wrap_plot+
    geom_line(data=its_preds_df, aes(x = acquired, y = jmd_predicted, color='ITS model fit'))+
    scale_color_manual(values = c('ITS model fit' = 'navy'), name = NULL)
}
