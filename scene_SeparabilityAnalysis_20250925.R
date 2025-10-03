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
  # ,cl = 'future'
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
{
  jmd_lm_list = mod_lister(data = results_meta_df_prepostharv,
                           formula = as.formula(paste('jmd_scaled',ITS_base_formula)))
  
  auc_lm_list = mod_lister(data = results_meta_df_prepostharv,
                           formula = as.formula(paste('logistic_AUC',ITS_base_formula)))
  
  #extract model coefficients and other info
  
  jmd_coefs = extract_coefs_lme4(jmd_lm_list) |>
    mutate(Estimate_rounded = round(Estimate, 3)) |>
    mutate(metric = 'Area Under Curve')
  auc_coefs = extract_coefs_lme4(auc_lm_list) |>
    mutate(Estimate_rounded = round(Estimate, 3)) |>
    mutate(metric = 'Jeffries-Matusita distance')
  
  coefs_combined_df = bind_rows(jmd_coefs, auc_coefs)
  
  #tile plot of coefficients results
  {
    auc_coefs_plot = ggplot(data = auc_coefs |> 
                              filter(!model_terms %in% c('months_since_harvest', '(Intercept)')) |>
                              filter(model_terms == 'harvested')
                            , aes(x = dataset, y = var_name))+
      geom_tile(aes(fill = Estimate_rounded))+
      geom_text(aes(label = Estimate_rounded))+
      scale_fill_viridis_c()+
      facet_grid(cols=vars(model_terms))+
      ggtitle('AUC ITS coefficients')+
      theme_classic()
    
    jmd_coefs_plot = ggplot(data = jmd_coefs |> 
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
{
  gs_months = c(5:9)
  
  results_meta_df_prepostharv_gs = results_meta_df_prepostharv |>
    filter(acquisition_month %in% gs_months)
  
  jmd_lm_list = mod_lister(data = results_meta_df_prepostharv_gs,
                           formula = as.formula(paste('jmd_scaled',ITS_base_formula)))
  
  auc_lm_list = mod_lister(data = results_meta_df_prepostharv_gs,
                           formula = as.formula(paste('logistic_AUC',ITS_base_formula)))
  
  #extract model coefficients and other info
  jmd_coefs = extract_coefs_lme4(jmd_lm_list) |>
    mutate(Estimate_rounded = round(Estimate, 3)) |>
    mutate(metric = 'Area Under Curve')
  auc_coefs = extract_coefs_lme4(auc_lm_list) |>
    mutate(Estimate_rounded = round(Estimate, 3)) |>
    mutate(metric = 'Jeffries-Matusita distance')
  
  coefs_combined_df = bind_rows(jmd_coefs, auc_coefs)
  
  #plot results
  {
    auc_coefs_plot = ggplot(data = auc_coefs |> 
                              filter(!model_terms %in% c('months_since_harvest', '(Intercept)')) |>
                              filter(model_terms == 'harvested')
                            , aes(x = dataset, y = var_name))+
      geom_tile(aes(fill = Estimate_rounded))+
      geom_text(aes(label = Estimate_rounded))+
      scale_fill_viridis_c()+
      facet_grid(cols=vars(model_terms))+
      ggtitle('AUC ITS coefficients')+
      theme_classic()
    
    jmd_coefs_plot = ggplot(data = jmd_coefs |> 
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
    # its_gs_coefs_plot = auc_coefs_plot+jmd_coefs_plot
    # # +plot_layout(guides='collect')
    
    # its_gs_coefs_plot
    
  }
}