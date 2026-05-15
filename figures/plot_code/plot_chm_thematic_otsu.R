# source('make_chm_change_validation_20250929.R')

source('derive lidar metrics_202501020.R')
library(patchwork)

#----load chm layers----

#metric to load
met = lid_met
print(met)

#function to load metric
load_mets = function(d){
  fl = list.files(d, full.names = T)
  rl = lapply(fl, function(f){
    r = rast(f, lyrs=met)
    return(r)
  })
  names(rl) = file_path_sans_ext(basename(fl))
  return(rl)
}

list_plot_mets = function(d){
  mets = load_mets(d)
  pl = lapply(mets, function(m){
    p = ggplot()+
      geom_spatraster(data=m)+
      scale_fill_viridis_c(na.value='white', guide=F)+
      theme_void()
    return(p)
  })
  
  wp = wrap_plots(pl,ncol=1)
}

#pre-thinning lidar

pre_p = list_plot_mets(lid_mets_pre_cropped_dir)

#post_thinning_lidar

post_p = list_plot_mets(lid_mets_post_dec_cropped_dir)

#chm change

change_p = list_plot_mets(lid_mets_change_cropped_dir)

(pre_p | post_p | change_p) #+ plot_layout(nrow=1)


