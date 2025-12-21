source('scene_global_stats_filtering_20251106.R')


#----set up data and helper functions for ITS----

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
  # mutate(time_since_harvest = ifelse(acquisition_date < harvest_start_date,
  #                                    as.numeric(acquisition_date - harvest_start_date),
  #                                    as.numeric(acquisition_date - harvest_finish_date)))
  mutate(time_since_harvest = case_when(
    acquisition_date < harvest_start_date ~ as.numeric(acquisition_date - harvest_start_date),
    acquisition_date > harvest_finish_date ~ as.numeric(acquisition_date - harvest_finish_date),
    .default = NA))


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
      
      if(any(str_detect(as.character(formula), '\\|'))){ #check for random effects
      
      mod = lmer(formula, data = dat) #lmer if random effects
      
      } else {
        
        mod = lm(formula, data = dat) #lm if no random effects
      }
      
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

# ITS_base_formula = '~ time_since_harvest * harvested + (1 | block_id)' #random intercepts, fixed slopes

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
      if(class(m)=='lm'){
        coefs_df$REML = NULL
      } else { 
        coefs_df$REML = REMLcrit(m)}
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

#define dataframe of model structures to test
its_models_df = tribble(
  ~type, ~model,
  'only_fixed', '~ time_since_harvest * harvested', # Only fixed effects - all blocks assumed identical
  'random_intercepts', '~ time_since_harvest * harvested + (1 | block_id)', # Random intercepts - blocks have different baselines
  'random_slopes_time', '~ time_since_harvest * harvested + (1 + time_since_harvest | block_id)', # Random intercepts + random slopes for time - blocks differ in baseline and time trends
  'random_slopes_full', '~ time_since_harvest * harvested + (1 + time_since_harvest + harvested | block_id)', # Random intercepts + random slopes for time and harvested - blocks differ in baseline, time effect, and harvesting effect
  'random_ints_&_slopes', '~ time_since_harvest * harvested + (1 + time_since_harvest * harvested | block_id)' # Maximal random effects - blocks differ in all effects including the time×harvested interaction
)

#----fit ITS models----

#fit mixed effects models with each structure,
response_var = 'jmd_scaled'

its_dir = file.path(ps_dir, paste0(response_var,'_linear_ITS')) 
dir.check(its_dir)

#define function for fitting models given input data, response variable, and a dataframe of candidate model structures
modfitter_max = function(response = 'jmd_scaled', dat, mods = its_models_df, outdir){
  
  #check if outdir exists
  dir.check(outdir)
  
  pblapply(1:nrow(mods), function(i){
    
    mod_type = mods$type[i]
    ITS_base_formula = mods$model[i]
    
    #make directory for model type
    mod_dir = file.path(outdir,mod_type)
    dir.check(mod_dir)
    
    #check for model summary stats file
    summ_f = file.path(mod_dir,paste0(mod_type,'_summary.csv'))
    if(!file.exists(summ_f)){
      #generate nested lists of models
      mod_list = mod_lister(data = dat, formula = as.formula(paste0(response,ITS_base_formula)))
      
      #save nested lists of models
      pblapply(names(mod_list), function(n){
        mod_ds_dir = file.path(mod_dir, n)
        dir.check(mod_ds_dir)
        
        mod_ds_list = mod_list[[n]]
        
        pblapply(names(mod_ds_list), function(x){
          mod_f = file.path(mod_ds_dir, paste0(x,'.rds'))
          if(!file.exists(mod_f)){
            mod = mod_ds_list[[x]]
            write_rds(mod,mod_f)
            }
        })
      })
      
      #extract coefficients and performance metrics
      coefs = extract_coefs_lme4(mod_list)
      coefs$metric = response
      coefs_model_type = mod_type
      write_csv(coefs, summ_f)
    }
    
    })
}

#test models for full timeseries
jmd_ft_dir = file.path(its_dir,'full_timeseries')
dir.check(jmd_ft_dir)

modfitter_max(dat = results_meta_df_prepostharv, outdir = jmd_ft_dir)

#compare AIC values for different model structures for full time series (random intercepts model has lowest AIC while also avoiding singular fits and convergence failure)
jmd_summ_stats_full_fl = list.files(jmd_ft_dir, recursive = T, full.names = T, pattern='\\.csv$')
jmd_summ_stats_full_df = lapply(jmd_summ_stats_full_fl, function(x){
  d = read_csv(x)
  mod = basename(file_path_sans_ext(x))
  d$model_type = mod
  return(d)
  }) |> 
  bind_rows() |>
  select(dataset,var_name, AIC, model_type) |>
  distinct()

aic_cols <- paste0(its_models_df$type,'_summary')

jmd_best <- jmd_summ_stats_full_df %>%
  select(dataset,var_name, AIC, model_type) |>
  distinct() |>
  pivot_wider(names_from = 'model_type', values_from = 'AIC')
  rowwise() %>%
  mutate(best_model = aic_cols[ which.min(c_across(all_of(aic_cols))) ]) %>%
  ungroup()

ggplot(jmd_summ_stats_full_df)+
  # geom_boxplot(aes(x = model_type, y = AIC))+
  geom_point(aes(x = model_type, y = AIC, color=var_name))+
  geom_line(aes(x = model_type, y = AIC, colour = var_name))+
  facet_wrap(vars(dataset))+
  theme(axis.text.x = element_text(angle=90))
  
# #get desired results and models for final model type for full time series
# jmd_final_mod_dir = list.files(jmd_ft_dir, full.names = T) |> purrr::keep(~ str_detect(.x, "random_intercepts"))
# jmd_lm_list = lapply(list.dirs(jmd_final_mod_dir, recursive = F), function(d){
#   dfl = list.files(d, full.names = T, pattern = '\\.rds$')
#   d_l = pblapply(dfl, function(x){cat(bname(d),'---',bname(x)); read_rds(x)})
#   names(d_l) = bname(dfl)
#   return(d_l)
# })

ITS_base_formula = its_models_df$model[its_models_df$type == 'random_intercepts']

#generate interrupted time series models for full time series
print('Generating ITS models for full time series')
{
  jmd_lm_list_full = mod_lister(data = results_meta_df_prepostharv,
                                formula = as.formula(paste('jmd_scaled',ITS_base_formula)))
  
  #extract model coefficients and other info
  jmd_coefs_full = extract_coefs_lme4(jmd_lm_list_full) |>
    mutate(Estimate_rounded = round(Estimate, 3)) |>
    mutate(metric = 'Scaled JMD')
}

#generate interrupted time series models just based on growing season data
print('Generating ITS models for growing season data')
{
  
  results_meta_df_prepostharv_gs = results_meta_df_prepostharv |>
    filter(acquisition_month %in% gs_months)
  
  jmd_lm_list_gs = mod_lister(data = results_meta_df_prepostharv_gs,
                              formula = as.formula(paste('jmd_scaled',ITS_base_formula)))
  
  #extract model coefficients and other info
  jmd_coefs_gs = extract_coefs_lme4(jmd_lm_list_gs) |>
    mutate(Estimate_rounded = round(Estimate, 3)) |>
    mutate(metric = 'Scaled JMD')
  
}

#combine full ts and gs ts into single dataframe
{
  jmd_coefs_combined = rbind(
    jmd_coefs_full |> mutate(mod_time_series = 'Full time series'),
    jmd_coefs_gs |> mutate(mod_time_series = 'May - October')
  )
  jmd_coefs_combined = jmd_coefs_combined |> 
    mutate()
  
}
