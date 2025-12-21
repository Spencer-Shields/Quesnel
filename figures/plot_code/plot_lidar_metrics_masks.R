source('scene_global_stats_filtering_20251106.R')

#load lidar change metrics
lid_mets_change_l = lapply(list.files(lid_mets_change_cropped_dir,full.names = T), function(x)rast(x,lyrs='z_p95'))
names(lid_mets_change_l) = basename(file_path_sans_ext(list.files(lid_mets_change_cropped_dir)))

#segment change metrics into thematic layers using otsu thresholds
otsu_thinning_lyrs = pblapply(block_ids, function(b){
  rc = lid_mets_change_l[[b]]
  o = otsu_df$chm_change_otsu_threshold[otsu_df$block_id==b]
  m = ifel(rc <= o ,1,0)
  levels(m) = tibble(z_p95=c(1,0), status = c('Thinned', 'Not thinned'))
  return(m)
})
list_plot(otsu_thinning_lyrs)

#get number of pixels in each category in each block
otsu_npix_summ = pblapply(block_ids, function(b){
  r = otsu_thinning_lyrs[[b]]
  v = as.data.frame(r, na.rm=T)|>
    mutate(block_id=b)
  return(v)
}) |> 
  bind_rows()|>
  group_by(status, block_id) |>
  summarise(n = n()) |>
  pivot_wider(names_from = status, values_from = n) 
otsu_npix_summ$total = otsu_npix_summ$Thinned + otsu_npix_summ$`Not thinned`
otsu_npix_summ$proportion_thinned = round(otsu_npix_summ$Thinned/otsu_npix_summ$total, 3)  
print(otsu_npix_summ)  
#
library(patchwork)

# mets_p_l = lapply(block_ids[1:4], function(b){
#   r = lid_mets_change_l[[b]]
#   names(r) = paste0(names(r),' change')
#   p = ggplot()+
#     geom_spatraster(data=lid_mets_change_l[[b]])+
#     scale_fill_viridis_c(na.value='white')+
#     theme_void()
#   return(p)
# })
# mets_p = wrap_plots(mets_p_l, ncol=1)
# mets_p

ots_mets_p_l = lapply(block_ids#[1:4]
                      , function(b){
                        r = lid_mets_change_l[[b]]
                        or = otsu_thinning_lyrs[[b]]
                        names(r) = paste0(names(r),' change')
                        rp = ggplot()+
                          geom_spatraster(data=r)+
                          scale_fill_viridis_c(na.value='white'
                                               ,name = 'dCHM'
                                               )+
                          theme_void()#+
                          # ggtitle('z_p95 change')
                        orp = ggplot()+
                          geom_spatraster(data=or)+
                          scale_fill_manual(na.value='white',
                                               labels = c('Thinned', 'Not thinned'),
                                               values = c('Thinned' = 'red3',
                                                          'Not thinned' = 'seagreen')
                                            ,name='Class'
                                            ,breaks = c("Thinned", "Not thinned"),  # removes NA from legend
                                               )+
                          theme_void()#+
                          # ggtitle('Thematic thinning map')
                        p = (rp+orp)
                        return(p)
                      })

# ots_mets_p_l = lapply(block_ids#[1:4]
#                       , function(b){
#                         
#                         r = lid_mets_change_l[[b]]
#                         or = otsu_thinning_lyrs[[b]]
#                         r = c(r,or)
#                         names(r) = c('z_p95 change','Thinning stratum')
#                         
#                         rp = ggplot()+
#                           geom_spatraster(data=r)+
#                           # scale_fill_viridis_c(na.value='white')+
#                           facet_wrap(~lyr)+
#                           theme_void()
#                         # orp = ggplot()+
#                         #   geom_spatraster(data=or)+
#                         #   scale_fill_discrete(na.value='white')+
#                         #   theme_void()
#                         # 
#                         # p = (rp+orp)
#                         p=rp
#                         return(p)
#                       })

ots_mets_p_l$`12L_C7`/ots_mets_p_l$`12N_T3`
