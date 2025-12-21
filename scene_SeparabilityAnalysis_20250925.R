source('scene_global_stats_filtering_20250925.R')

#----plot temporal distribution of final planetscope dataset for each block----

ggplot(results_meta_df |>
         filter(block_pixel_stratum == 'Total'
                ,dataset == 'Z'
                ,var_name=='blue'
         ) |>
         mutate(acquisition_julian = yday(acquisition_date))
       # mutate(acquisition_year = year(acquisition_date)
       #        ,acquisition_month = month(acquisition_date)
       #        ,acquisition_day = day(acquisition_date)
       #        )
)+
  geom_point(aes(x = acquisition_julian, y=acquisition_year), alpha=0.3)+
  facet_grid(rows=vars(block_id))


#----plot thinning vs non-thinning values for specified features and plots over time----

desired_blocks = c(
  "12L_B8C3"
  ,
  "12L_C4"
  ,
  "12L_C5"
  ,
  "12L_C7"
  ,
  "12L_D345"
  ,
  "12N_1X1W"
  ,
  "12N_T3"
  ,
  "12N_V1"
)
# desired_datasets = 

#----interrupted time series analysis----

library(lme4)
library(lmerTest)
library(performance)

#filter out dates that occur between start and end of harvesting
results_meta_df_prepostharv = results_meta_df |>
  mutate(harvest_in_progress = ifelse(acquisition_date >= harvest_start_date & acquisition_date < harvest_finish_date,
                                      1, 0)) |>
  filter(harvest_in_progress == 0) |>
  # filter(acquisition_month %in% c(5,6,7,8,9)) |> #only model for growing season months
  #define dummy variable that measures time before harvest start and after harvest end
  mutate(time_since_harvest = ifelse(acquisition_date < harvest_start_date,
                                     as.numeric(acquisition_date - harvest_start_date),
                                     as.numeric(acquisition_date - harvest_finish_date)))

#define functions for generating lists of linear mixed models for every combination of dataset and variable
mod_lister = function(data, 
                      formula, 
                      datasets = unique(data[['dataset']]), 
                      feats = unique(data[['var_name']])
){
  #produce lists of linear mixed models. Output will be a nested list structure
  #"formula" should be an object defined using the as.formula() function, data should be a data frame
  
  mod_list = pblapply(1:length(datasets), function(i){
    
    dsi = datasets[i]
    
    mod_l = pblapply(1:length(feats), function(j){
      
      varj = feats[j]
      
      #print for debugging
      print(paste0(
        'Processing ',i,'-',j,', ',dsi,'-',varj
      ))
      
      dat = data |> filter(dataset == dsi,
                           var_name == varj)
      
      
      # jmd_lm = glmmTMB(jmd_scaled ~ acquisition_date_ordinal + harvested + (acquisition_date_ordinal * harvested | block_id)
      #                  , data = dat
      #                  ,family=ordbeta(link="logit")) #ordered beta regression since the dataset includes zeroes and ones
      
      mod = lmer(formula, data = dat)
      
      return(mod)
    })
    names(mod_l) = feats
    return(mod_l)
  }
  )
  
  
  names(mod_list) = datasets
  return(mod_list)
} 

#define ITS base formula
# ITS_base_formula = '~ 1 + time_since_harvest * harvested +
#   (1 + time_since_harvest * harvested | block_id)' #fixed and random slopes and intercepts (doesn't work as well)

ITS_base_formula = '~ time_since_harvest * harvested + (1 | block_id)' #random intercepts, fixed slopes


#define function for extracting coefficients from nested lists of mixed models

datasets = unique(results_meta_df_prepostharv[['dataset']]) 
ps_feats = unique(results_meta_df_prepostharv[['var_name']])

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
      
      #get additional performance metrics using the performance function
      perf = performance(m)
      
      coefs_df$R2_conditional = perf$R2_conditional
      coefs_df$R2_marginal = perf$R2_marginal
      coefs_df$RMSE = perf$RMSE
      
      
      return(coefs_df)
    })
    df = bind_rows(df_l)
    return(df)
  })
  df = bind_rows(coefs_df_l)
  rownames(df) = NULL
  return(df)
}

#generate interrupted time series models for full time series
print('Generating ITS models for full time series')
{
  jmd_lm_list = mod_lister(data = results_meta_df_prepostharv,
                           formula = as.formula(paste('jmd_scaled',ITS_base_formula)))
  
  # auc_lm_list = mod_lister(data = results_meta_df_prepostharv,
  #                          formula = as.formula(paste('logistic_AUC',ITS_base_formula)))
  
  #extract model coefficients and other info
  
  jmd_coefs = extract_coefs_lme4(jmd_lm_list) |>
    mutate(Estimate_rounded = round(Estimate, 3)) |>
    mutate(metric = 'Area Under Curve')
  # auc_coefs = extract_coefs_lme4(auc_lm_list) |>
  #   mutate(Estimate_rounded = round(Estimate, 3)) |>
  #   mutate(metric = 'Jeffries-Matusita distance')
  
  # coefs_combined_df = bind_rows(jmd_coefs
  #                               # , auc_coefs
  #                               )
  
  #tile plot of coefficients results
  {
    # auc_coefs_plot = ggplot(data = auc_coefs |> 
    #                           filter(!model_terms %in% c('months_since_harvest', '(Intercept)')) |>
    #                           filter(model_terms == 'harvested')
    #                         , aes(x = dataset, y = var_name))+
    #   geom_tile(aes(fill = Estimate_rounded))+
    #   geom_text(aes(label = Estimate_rounded))+
    #   scale_fill_viridis_c()+
    #   facet_grid(cols=vars(model_terms))+
    #   ggtitle('AUC ITS coefficients')+
    #   theme_classic()
    
    jmd_full_coefs_plot = ggplot(data = jmd_coefs |> 
                                   filter(!model_terms %in% c('months_since_harvest', '(Intercept)')) |>
                                   filter(model_terms == 'harvested')
                                 , aes(x = dataset, y = var_name))+
      geom_tile(aes(fill = Estimate_rounded))+
      geom_text(aes(label = Estimate_rounded))+
      scale_fill_viridis_c()+
      facet_grid(cols=vars(model_terms))+
      ggtitle('JMD ITS coefficients')+
      theme_classic()
    
    # x11() 
    # library(patchwork) 
    # its_full_coefs_plot = auc_coefs_plot+jmd_coefs_plot
    # # +plot_layout(guides='collect')
    # 
    # its_full_coefs_plot
    }
}

#generate interrupted time series models just based on growing season data
print('Generating ITS models for growing season data')
{
  gs_months = c(5:10)
  
  results_meta_df_prepostharv_gs = results_meta_df_prepostharv |>
    filter(acquisition_month %in% gs_months)
  
  jmd_lm_list = mod_lister(data = results_meta_df_prepostharv_gs,
                           formula = as.formula(paste('jmd_scaled',ITS_base_formula)))
  
  # auc_lm_list = mod_lister(data = results_meta_df_prepostharv_gs,
  #                          formula = as.formula(paste('logistic_AUC',ITS_base_formula)))
  
  #extract model coefficients and other info
  jmd_coefs = extract_coefs_lme4(jmd_lm_list) |>
    mutate(Estimate_rounded = round(Estimate, 2)) |>
    mutate(metric = 'Scaled JMD')
  # auc_coefs = extract_coefs_lme4(auc_lm_list) |>
  #   mutate(Estimate_rounded = round(Estimate, 3)) |>
  #   mutate(metric = 'Jeffries-Matusita distance')
  
  # coefs_combined_df = bind_rows(jmd_coefs, auc_coefs)
  
  #plot results
  {
    # auc_coefs_plot = ggplot(data = auc_coefs |> 
    #                           filter(!model_terms %in% c('months_since_harvest', '(Intercept)')) |>
    #                           filter(model_terms == 'harvested')
    #                         , aes(x = dataset, y = var_name))+
    #   geom_tile(aes(fill = Estimate_rounded))+
    #   geom_text(aes(label = Estimate_rounded))+
    #   scale_fill_viridis_c()+
    #   facet_grid(cols=vars(model_terms))+
    #   ggtitle('AUC ITS coefficients')+
    #   theme_classic()
    
    # var_order <- rev(unique(results_meta_df$var_name))
    var_order = rev(
      c(unique(results_meta_df$var_name)[1:8],
        sort(unique(results_meta_df$var_name)[9:length(unique(results_meta_df$var_name))])
      )
    )
    
    jmd_coefs_gs_plot <- ggplot(
      data = jmd_coefs |>
        # filter(!model_terms %in% c('months_since_harvest', '(Intercept)')) |>
        filter(model_terms == 'harvested') |>
        mutate(var_name = factor(var_name, levels = var_order)),  # ⬅️ set order
      aes(x = dataset, y = var_name)
    ) +
      geom_tile(aes(fill = Estimate_rounded)) +
      geom_text(aes(label = Estimate_rounded)) +
      scale_fill_viridis_c(name = 'JMD step change') +
      labs(y = 'Spectral feature', x = 'Dataset') +
      theme_classic()
    
    x11()
    library(patchwork)
    # its_gs_coefs_plot = auc_coefs_plot + jmd_coefs_gs_plot
    # +plot_layout(guides='collect')
    
    # its_gs_coefs_plot
    jmd_coefs_gs_plot
    
    jmd_r2_gs_plot = ggplot(
      data = jmd_coefs |>
        select(dataset, var_name, R2_conditional) |>
        distinct() |>
        mutate(var_name = factor(var_name, levels = var_order),
               R2_conditional = round(R2_conditional,2)
        )
      ,aes(x = dataset, y = var_name)
    ) +
      geom_tile(aes(fill = R2_conditional)) +
      geom_text(aes(label = R2_conditional)) +
      scale_fill_viridis_c(name = 'R2') +
      labs(y = 'Spectral feature', x = 'Dataset') +
      theme_classic()
    
    x11()
    jmd_r2_gs_plot
    
    jmd_coefs_gs_plot + (jmd_r2_gs_plot+theme(axis.title.y = element_blank()))
    }
}

#----plot median pixel values before and after thinning----
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
    '12L_D345'
  )
  
  subset_df = results_meta_df |>
    filter(
      block_pixel_stratum != 'Total'
      ,var_name %in% feats_to_plot
      ,block_id %in% blocks_to_plot
    ) |>
    mutate(dataset = case_match(dataset, 'Non-normalized' ~ 'NN', .default = dataset)) |>
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
    ggtitle(paste0('Growing season ',unique(subset_df$var_name),' (Block ',unique(subset_df$block_id),')'))+
    theme_classic()+
    theme(
      #remove x axis text on top of plot (put there by ggbreaks)
      axis.text.x.top = element_blank(),
      axis.ticks.x.top = element_blank(),
      #format lower x axis text
      axis.text.x = element_text(angle=90)
    )
  pix_plot_gs
}

#full timeseries plot
{
  pix_plot_full = ggplot(subset_df)+
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
    # #gaps to show only growing season
    # scale_x_break(c(as.POSIXct('2020-11-01 00:00:00 UTC'),as.POSIXct('2021-04-30 23:59:59 UTC')))+
    # scale_x_break(c(as.POSIXct('2021-11-01 00:00:00 UTC'),as.POSIXct('2022-04-30 23:59:59 UTC')))+
    # scale_x_break(c(as.POSIXct('2022-11-01 00:00:00 UTC'),as.POSIXct('2023-04-30 23:59:59 UTC')))+
    # scale_x_break(c(as.POSIXct('2023-11-01 00:00:00 UTC'),as.POSIXct('2024-04-30 23:59:59 UTC')))+
    #format x axis date/time info
    scale_x_datetime(
      date_labels = "%b %Y",       # e.g., "Jan 2021"
      date_breaks = "3 months",     # choose interval for ticks
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
    ggtitle(paste0('Full time series ',unique(subset_df$var_name),' (Block ',unique(subset_df$block_id),')'))+
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
  pix_plot_full + pix_plot_gs + plot_layout(guides='collect', axis_titles = 'collect'
                                            , ncol=1
                                            )
}
