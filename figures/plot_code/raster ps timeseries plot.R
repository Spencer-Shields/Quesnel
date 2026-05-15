source('scene_SeparabilityAnalysis_20251103.R')
library(patchwork)

subset_df = results_meta_df |>
  filter(block_id == '12L_C7', dataset == 'Non-normalized', clear_confidence_percent>90) |>
  select(acquisition_date, id, file_path, acquired, harvest_start_date, harvest_finish_date) |>
  distinct()

desired_ids = c(
  '20210424_191331_92_227a'
  ,'20220522_184820_70_2461'
  ,'20230513_182107_37_24a1'
  # ,'20240509_192241_42_24f2'
  ,'20240505_183451_59_2447'
  
)  

plotxrgb = function(x,newwin=T,falsecol=F){
  
  r = rast(x)
  
  max_val = 0.09
  
  if(falsecol){
    x11()
    terra::plotRGB(x=r, r = 8, g = 6, b = 4, scale = max_val)
  } else {
    x11()
    terra::plotRGB(x=r, r = 6, g = 4, b = 2
                   # , scale = max_val
                   , stretch = 'lin'
    )
  }
}

#load each desired raster file
rl = lapply(desired_ids, function(i){
  f = subset_df$file_path[subset_df$id==i]
  r = rast(f) 
})
names(rl) = desired_ids

#get min and max values of desired bands for each raster (to stretch values for display)
# desired_lyrs = c('red','green','blue') #true color image
desired_lyrs = c('red', 'green', 'blue') #false color image
minmax_df = lapply(1:length(rl), function(i){
  r = rl[[i]]
  gdf = global(r[[desired_lyrs]], c('min','max'), na.rm=T)
  
  gdf$lyr = rownames(gdf)
  gdf$id = names(rl)[i]
  
  return(gdf)
}) |> bind_rows()

glob_max = max(minmax_df$max)
glob_min = min(minmax_df$min)


#make rgb plots for each raster
plot_l = lapply(1:length(rl), function(i){
  r = rl[[i]]
  r = r[[desired_lyrs]]
  
  #stretch
  r = 255*(r-glob_min)/(glob_max-glob_min)
  
  # p = 
  p = ggplot()+
    geom_spatraster_rgb(data = r, stretch = 'lin')+
    # scale_fill_continuous(na.value = 'white')+
    ggtitle(subset_df$acquisition_date[subset_df$id==names(rl)[i]])+
    theme_void()
})

wrap_plots(plot_l,nrow=1)


#add scale bar
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

plot_l[[1]] = plot_l[[1]] +
  annotation_scale(
    location = "bl",
    width_hint = 0.25,
    height = unit(0.15, "cm"),
    text_cex = 0.6,
    pad_x = unit(0.2, "cm"),
    pad_y = unit(0.05, "cm"),
    line_width = 0.5
  )

wrap_plots(plot_l, nrow = 1)
