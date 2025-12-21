source('scene_global_stats_filtering_20251106.R')

#12N_T3
lid_f = list.files(lid_mets_change_cropped_dir, full.names = T, pattern = '12N_T3.tif')
lid_r = rast(lid_f, lyrs='z_p95')
v = as.data.frame(lid_r,na.rm=T)

ggplot(v) +
  geom_histogram(aes(z_p95), bins = 256) +
  geom_vline(
    aes(xintercept = otsu_df$chm_change_otsu_threshold[otsu_df$block_id == '12N_T3']),
    linetype = "dashed",
    color = "red"
  ) +
  annotate("text", 
           x = otsu_df$chm_change_otsu_threshold[otsu_df$block_id == '12N_T3'] - 1, 
           y = Inf, 
           label = "Thinned", 
           vjust = 2, hjust = 1, size = 4) +
  annotate("text", 
           x = otsu_df$chm_change_otsu_threshold[otsu_df$block_id == '12N_T3'] + 1, 
           y = Inf, 
           label = "Not thinned", 
           vjust = 2, hjust = 0, size = 4) +
  labs(x = 'dCHM',
       y = 'Num. pixels')+
  theme_classic()

#12L_C5
lid_f = list.files(lid_mets_change_cropped_dir, full.names = T, pattern = '12L_C5.tif')
lid_r = rast(lid_f, lyrs='z_p95')
v = as.data.frame(lid_r,na.rm=T)

ggplot(v) +
  geom_histogram(aes(z_p95), bins = 256) +
  geom_vline(
    aes(xintercept = otsu_df$chm_change_otsu_threshold[otsu_df$block_id == '12L_C5']),
    linetype = "dashed",
    color = "red"
  ) +
  annotate("text", 
           x = otsu_df$chm_change_otsu_threshold[otsu_df$block_id == '12L_C5'] - 1, 
           y = Inf, 
           label = "Thinned", 
           vjust = 2, hjust = 1, size = 4) +
  annotate("text", 
           x = otsu_df$chm_change_otsu_threshold[otsu_df$block_id == '12L_C5'] + 1, 
           y = Inf, 
           label = "Not thinned", 
           vjust = 2, hjust = 0, size = 4) +
  labs(x = 'dCHM',
       y = 'Num. pixels')+
  theme_classic()