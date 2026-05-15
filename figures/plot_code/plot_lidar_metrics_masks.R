source('scene_global_stats_filtering_20251106.R')

print(paste('Loading metric:', lid_met))

#load lidar change metrics
lid_mets_change_l = lapply(list.files(lid_mets_change_cropped_dir,full.names = T), function(x){
  r = rast(x,lyrs=lid_met)
  })
names(lid_mets_change_l) = basename(file_path_sans_ext(list.files(lid_mets_change_cropped_dir)))

#segment change metrics into thematic layers using otsu thresholds
otsu_thinning_lyrs = pblapply(block_ids, function(b){
  rc = lid_mets_change_l[[b]]
  o = otsu_df$chm_change_otsu_threshold[otsu_df$block_id==b]
  m = ifel(rc <= o ,1,0)
  levels(m) = tibble(z_above2=c(1,0), status = c('Thinned', 'Not thinned'))
  return(m)
})
# list_plot(otsu_thinning_lyrs)

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

####----plot spatrasters----

#make scale bar to add to plots
library(ggspatial)

scalebar_layer =
  annotation_scale(
    location = "bl",
    width_hint = 0.25,
    height = unit(0.15, "cm"),
    text_cex = 0.6,
    pad_x = unit(0.2, "cm"),
    pad_y = unit(0.05, "cm"), #distance from bottom of panel
    line_width = 0.5
  )

#assemble plots of lidar metrics and thematic maps
ots_mets_p_l = lapply(sort(block_ids), function(b){
  
                        r = lid_mets_change_l[[b]]*100 #multiply by 100 to get percent dCC
                        or = otsu_thinning_lyrs[[b]]
                        names(r) = paste0(names(r),' change')
                        rp = ggplot()+
                          geom_spatraster(data=r)+
                          scale_fill_viridis_c(na.value='white'
                                               ,name = 'dCC(%)'
                                               # , limits = c(-1,0.75)
                                               ,limits = c(-100,75)
                                               )+
                          ggtitle(b)+
                          scalebar_layer+
                          theme_void()+
                          theme(
                            plot.margin = margin(t = 0, r = 0, b = 10, l = 0) #pad plot margin so scale bar doesn't overlap plot
                            ,title = element_text(size=8)
                          
                            )
                  
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

#plot all blocks stacked together with patchwork
final_plot = wrap_plots(ots_mets_p_l, ncol=2, guides = 'collect')
final_plot

#plot all blocks in row for presentation
final_plot_wide = wrap_plots(ots_mets_p_l, nrow = 2, guides = 'collect')
final_plot_wide

####----plot histograms----

#get dataframe of values for lidar change metric for each block

lid_change_df = pblapply(sort(block_ids), function(b){
  r = lid_mets_change_l[[b]]
  df = r |> 
    as.data.frame(na.rm=T) |>
    mutate(block_id = b)
    
}) |>
  bind_rows() |>
  mutate(dcc_percent = z_above2*100)

otsu_df_ = otsu_df |>
  mutate(otsu_threshold = round(chm_change_otsu_threshold*100, digits = 1))
  # mutate(otsu_threshold = sprintf('%.1f', chm_change_otsu_threshold*100))

ggplot()+
  #color plot area
  #color plot area
  geom_rect(data = otsu_df_, aes(xmin = -Inf, xmax = otsu_threshold, ymin = -Inf, ymax = Inf,
                                       fill = "Thinned"), alpha = 1)+
  geom_rect(data = otsu_df_,  aes(xmin = otsu_threshold, xmax = Inf, ymin = -Inf, ymax = Inf,
                                        fill = "Not thinned"), alpha = 0.5)+
  scale_fill_manual(values = c("Thinned" = "lightyellow", "Not thinned" = "lightblue"),
                    name = NULL)+
  #add density plot
  geom_density(data = lid_change_df, aes(dcc_percent))+
  xlab('dCC (%)')+
  facet_wrap(vars(block_id), scales='free')+
  #add vertical line and text label for thinning threshold
  geom_vline(data = otsu_df_, aes(xintercept = otsu_threshold
                                        , linetype = 'Thinning threshold'
  )
  ,color='darkred'
  )+
  scale_linetype_manual(values = c("Thinning threshold" = "solid"),
                        name = NULL)+
  geom_text(data = otsu_df_, aes(x = otsu_threshold-20
                                 # , label = round(otsu_threshold, 1)
                                 ,label = sprintf('%.1f', otsu_threshold)
                                 ),
            y = Inf, vjust = 1.2, hjust = 0.5, color = 'darkred',
            # angle = 90,
            size = 3)+
  ylab('Density')+
  theme_classic()

#----plot for slide----

ots_mets_p_l$`12L_C7`

