mv_form = bf(
  mvbind(jmd_scaled, logistic_AUC) ~ 
    months_since_harvest + harvested + months_since_harvest * harvested + (1 | block_id) 
  # + cor_ar(~ acquisition_date_ordinal | block_id, p = 1)
  # + ar(time = acquisition_date_ordinal, gr = block_id, p = 1)
)


m = ordbetareg(
  # formula = mvbind(jmd_scaled, logistic_AUC) ~ months_since_harvest + harvested + (1 + months_since_harvest + harvested | block_id)
  mv_form,
  # prior = c(prior(normal(0,1), class = Intercept),
  #           prior(normal(0,1), class = b)
  #           ),
  data = dat_ds_var,
  true_bounds = c(0,1),
  chains = 4, cores = 4, iter = 4000,
  control = list(adapt_delta = 0.99)  # helps with convergence in complex models
)

m_s = summary(m)

p = predict(m, dat_ds_var)

logistic_auc_predictions = p[,, 'logisticAUC'] |>
  as.data.frame()
names(logistic_auc_predictions) = paste0('logAUC_',names(logistic_auc_predictions))

jmd_predictions = p[,, 'jmdscaled'] |>
  as.data.frame()
names(jmd_predictions) = paste0('jmd_',names(jmd_predictions))

dat_ds_var2 = cbind(dat_ds_var, logistic_auc_predictions, jmd_predictions)

plog_p = ggplot(data = dat_ds_var2)+
  geom_ribbon(aes(ymin = logAUC_Q2.5, ymax = logAUC_Q97.5, y = logAUC_Estimate, x = acquisition_date2)
              ,fill = 'seagreen'
              , alpha = 0.2)+
  geom_point(aes(x = acquisition_date2, y = logistic_AUC))+
  geom_line(aes(x = acquisition_date2, y = logAUC_Estimate)) +
  facet_grid(rows = vars(block_id), cols = vars(var_name)
             , scales = 'free'
  )
pjmd_p = ggplot(data = dat_ds_var2)+
  geom_ribbon(aes(ymin = jmd_Q2.5, ymax = jmd_Q97.5, y = jmd_Estimate, x = acquisition_date2)
              ,fill = 'seagreen'
              , alpha = 0.2)+
  geom_point(aes(x = acquisition_date2, y = jmd_scaled))+
  geom_line(aes(x = acquisition_date2, y = jmd_Estimate)) +
  facet_grid(rows = vars(block_id), cols = vars(var_name)
             , scales = 'free'
  )


plog_p+pjmd_p
