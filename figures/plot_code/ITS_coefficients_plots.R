source('scene_SeparabilityAnalysis_20251103.R')
library(ggtext)

jmd_coefs_combined_ = jmd_coefs_combined |>
  mutate(dataset = case_match(dataset,
                              'Non-normalized' ~ 'NN',
                              'Zrobust' ~ 'Zr',
                              .default=dataset)) |>
  # Add significance stars
  mutate(significance = case_when(
    `Pr(>|t|)` < 0.001 ~ "a",
    `Pr(>|t|)` < 0.01 ~ "b",
    `Pr(>|t|)` < 0.05 ~ "c",
    TRUE ~ ""
  )) |>
  #text label
  mutate(plot_label = as.character(Estimate_rounded))|>
  mutate(plot_label = case_match(plot_label,
                                 '0.1' ~ '0.10',
                                 '0.2' ~ '0.20',
                                 '0.3' ~ '0.30',
                                 .default=plot_label))|>
  # Combine estimate with significance
  mutate(label_with_sig = ifelse(
    significance != "",
    paste0(plot_label, "^", significance),
    plot_label
  ))


#define spectral variable order for y axis
var_order = rev(
  c(unique(results_meta_df$var_name)[1:8],
    sort(unique(results_meta_df$var_name)[9:length(unique(results_meta_df$var_name))])
  )
)

#----plot step-change coefficient for full and may-oct time series----

jmd_stepchange_plot <- ggplot(
  data = jmd_coefs_combined_ |>
    filter(model_terms == 'harvested') |>
    mutate(var_name = factor(var_name, levels = var_order)),  # ⬅️ set order
  aes(x = dataset, y = var_name)
) +
  geom_tile(aes(fill = Estimate_rounded)) +
  geom_text(aes(label = label_with_sig),parse = TRUE) +
  scale_fill_viridis_c(name = 'JMD step change') +
  labs(y = 'Spectral feature', x = 'Dataset') +
  facet_wrap(vars(mod_time_series),ncol=2)+
  # ggtitle('May - October')+
  theme_classic()

#----plot intercepts for full and may-oct time series----

jmd_intercept_plot <- ggplot(
  data = jmd_coefs_combined_ |>
    filter(model_terms == '(Intercept)') |>
    mutate(var_name = factor(var_name, levels = var_order)),  # ⬅️ set order
  aes(x = dataset, y = var_name)
) +
  geom_tile(aes(fill = Estimate_rounded)) +
  geom_text(aes(label = label_with_sig),parse = TRUE) +
  scale_fill_viridis_c(name = 'Intercept') +
  labs(y = 'Spectral feature', x = 'Dataset') +
  facet_wrap(vars(mod_time_series),ncol=2)+
  # ggtitle('May - October')+
  theme_classic()

#----plot pre-harvest slope for full and may-oct time series (no significant slopes----

jmd_slope_plot <- ggplot(
  data = jmd_coefs_combined_ |>
    filter(model_terms == 'time_since_harvest') |>
    mutate(var_name = factor(var_name, levels = var_order)),  # ⬅️ set order
  aes(x = dataset, y = var_name)
) +
  geom_tile(aes(fill = Estimate_rounded)) +
  geom_text(aes(label = label_with_sig),parse = TRUE) +
  scale_fill_viridis_c(name = 'Slope') +
  labs(y = 'Spectral feature', x = 'Dataset') +
  facet_wrap(vars(mod_time_series),ncol=2)+
  # ggtitle('May - October')+
  theme_classic()

#----plot jmd slope change (there is none)----
jmd_slopechange_plot <- ggplot(
  data = jmd_coefs_combined_ |>
    filter(model_terms == 'time_since_harvest:harvested') |>
    mutate(var_name = factor(var_name, levels = var_order)),  # ⬅️ set order
  aes(x = dataset, y = var_name)
) +
  geom_tile(aes(fill = Estimate_rounded)) +
  geom_text(aes(label = label_with_sig),parse = TRUE) +
  scale_fill_viridis_c(name = 'Slope change') +
  labs(y = 'Spectral feature', x = 'Dataset') +
  facet_wrap(vars(mod_time_series),ncol=2)+
  # ggtitle('May - October')+
  theme_classic()

#----get median coefficients----

cat('median ITS intercept full TS:',
    median(jmd_coefs_combined_$Estimate[jmd_coefs_combined_$model_terms=='(Intercept)']))

cat('median ITS intercept snow_free TS:',
    median(jmd_coefs_combined_$Estimate[jmd_coefs_combined_$model_terms=='(Intercept)' & jmd_coefs_combined_$mod_time_series!='Full time series']))

cat('median ITS step-change full TS:',
    median(jmd_coefs_combined_$Estimate[jmd_coefs_combined_$model_terms=='harvested']))

cat('median ITS step-change snow_free TS:',
    median(jmd_coefs_combined_$Estimate[jmd_coefs_combined_$model_terms=='harvested' & jmd_coefs_combined_$mod_time_series!='Full time series']))

#----plot RMSE values----

jmd_rmse_plot <- ggplot(
  data = jmd_coefs_combined_ |>
    select(dataset, var_name, RMSE, mod_time_series) |>
    mutate(var_name = factor(var_name, levels = var_order)),  # ⬅️ set order
  aes(x = dataset, y = var_name)
) +
  geom_tile(aes(fill = RMSE)) +
  geom_text(aes(label = round(RMSE,2)),parse = TRUE) +
  scale_fill_viridis_c(name = 'RMSE') +
  labs(y = 'Spectral feature', x = 'Dataset') +
  facet_wrap(vars(mod_time_series),ncol=2)+
  # ggtitle('May - October')+
  theme_classic()
