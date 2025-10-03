source('scene_setup_preprocessing_20250909.R')
library(data.table)
# library(supercells)
# library(dbscan)

#----Load data----
dir = 'data/Quesnel_thinning'
d_l = list.dirs(dir, full.names = T)

chm_post_dir = d_l[str_detect(d_l, 'chm_post')]
chm_pre_dir = d_l[str_detect(d_l, 'chm_pre')]

chm_post_files = list.files(chm_post_dir, full.names = T)
chm_post_files = chm_post_files[str_detect(chm_post_files, 'tif$')]
chm_pre_files = list.files(chm_pre_dir, full.names = T)
chm_pre_files = chm_pre_files[str_detect(chm_pre_files, 'tif$')]

#get names of thinning blocks that have pre and post-thinning lidar data
chm_post_ids = file_path_sans_ext(basename(chm_post_files)) %>% str_extract(".*(?=_chm_)")
chm_pre_ids = file_path_sans_ext(basename(chm_pre_files)) %>% str_extract(".*(?=_chm_)")
chm_common_ids = intersect(chm_pre_ids, chm_post_ids)

post_rasts_l = lapply(chm_post_files, rast)
names(post_rasts_l) = chm_post_ids

pre_rasts_l = lapply(chm_pre_files, rast)
names(pre_rasts_l) = chm_pre_ids

#rasters which cover the same blocks
post_rasts_common_l = post_rasts_l[chm_common_ids]
pre_rasts_common_l = pre_rasts_l[chm_common_ids]

#load data for plotting
aoi = vect('data/Quesnel_thinning/AOI_fullsite.geojson')
if(crs(aoi) != crs(pre_rasts_common_l[[1]])){
  aoi = project(aoi, crs(pre_rasts_common_l[[1]]))}


blocks = vect('data/Quesnel_thinning/12l_12n_bdy.gpkg')
if(crs(blocks) != crs(pre_rasts_common_l[[1]])){
  blocks = project(blocks, crs(pre_rasts_common_l[[1]]))}

areas = expanse(blocks)
total_area_ha = sum(areas)/10000

#----Calculate difference rasters----

chm_change_dir = 'data/Quesnel_thinning/chm_change'
dir.check(chm_change_dir)

change_rasts_l = pblapply(X = 1:length(post_rasts_common_l), function(i){
  post_rast = post_rasts_common_l[[i]]
  pre_rast = pre_rasts_common_l[[i]]
  
  #align rasters
  if(crs(post_rast)!=crs(pre_rast)){
    post_rast = project(post_rast, pre_rast)
  }
  if(ext(post_rast) != ext(pre_rast)){
    post_rast = resample(post_rast, pre_rast)
  }
  #calculate change
  change_rast = post_rast - pre_rast
  
  #return change
  change_rast
})

names(change_rasts_l) = chm_common_ids
# list_plot(change_rasts_l)

#----Export difference rasters----

lapply(1:length(change_rasts_l), function(i){
  filename = paste0(chm_change_dir,'/',chm_common_ids[i],'.tif')
  if(!file.exists(filename)){
    writeRaster(change_rasts_l[[i]], filename)
  }
})

#----Look at histograms of difference rasters to establish harvest threshold----

change_rasts_clipped = pblapply(1:length(change_rasts_l), function(i){
  rast_id = names(change_rasts_l)[i]
  block = st_as_sf(blocks) %>% filter(BLOCKNUM == rast_id)
  clipped_rast = crop(change_rasts_l[[i]], block, mask = T)
  return(clipped_rast)
})
names(change_rasts_clipped) = names(change_rasts_l)

change_df_l = pblapply(1:length(change_rasts_clipped), function(i){
  d = data.table(
    vals = values(change_rasts_clipped[[i]])
    ,
    block = names(change_rasts_clipped)[i]
  )
})
change_df = data.table::rbindlist(change_df_l)

#run otsu on each raster
otsu_l = pbsapply(names(change_rasts_clipped), function(x){
  d = change_df %>%
    filter(block == x) 
  t = otsuThresholdCpp(values = d[['vals.Z']], bins = 256)
})

#plot otsu thresholds on density plots of change rasters
change_df = change_df %>% 
  left_join(data.frame(block = names(otsu_l), otsu_threshold = otsu_l), by = 'block')

# otsu_t = round(otsu_l[['12N_T3']],2)

# #plot all
# density_otsu_plot_0.25m=ggplot(data = change_df)+
#   # annotate('rect', xmin=-Inf, xmax=otsu_t, ymin=-Inf, ymax = Inf,
#   #          fill = 'lightblue', alpha = 0.3)+
#   # annotate('text', x = -20, y = 0.16, label = 'Thinned', color = 'blue')+
#   # annotate('rect', xmin= otsu_t, xmax = Inf, ymin=-Inf, ymax = Inf,
#   #          fill = 'lightyellow', alpha = 0.5)+
#   # annotate('text', x = 20, y = 0.16, label = 'Not thinned', color = 'yellow4')+
#   geom_density(aes(x = vals.Z),fill='seagreen')+
#   geom_vline(aes(xintercept = otsu_threshold
#                  # , linetype = 'Otsu threshold'
#   )
#   ,color='darkred'
#   )+
#   # annotate("text", x = otsu_t+1, y = 0.16, label = paste0("Thinning threshold: ",otsu_t)
#   #          , angle = 0, vjust = -0.5, hjust = 0, color = 'darkred')+
#   ylab('Density')+
#   theme_classic()+
#   xlab("Height change (m)")+
#   theme(legend.position = NULL)+
#   facet_grid(vars(block))

#save otsu thresholds
otsu_tbl = tibble(
  block_id = names(otsu_l),
  chm_change_otsu_threshold = otsu_l
)

otsu_file = paste0(dir,'/','Otsu_change_thresholds.csv')
if(!file.exists(otsu_file)){write_csv(otsu_tbl,otsu_file)}

#----save change validation rasters which are cropped to the blocks, resampled and reprojected to match basemap data----

#get first file from each subdirectory so that there's one planetscope raster for each thinning block
bm_dir = 'data/planet_basemaps/global_monthly/CroppedMosaics_Indices_BlockClipped'
bm_subdirs = list.dirs(bm_dir)
bm_subdirs = bm_subdirs[bm_subdirs != bm_dir]
bm_subdirs = bm_subdirs[!str_detect(bm_subdirs, 'NoChange')]
bm_rasts = lapply(bm_subdirs, function(x){
  fl = list.files(x, full.names = T)
  r = rast(fl[1])
  return(r)
})
bm_ids = basename(bm_subdirs)
names(bm_rasts) = bm_ids

chm_change_dir_bm = paste0(chm_change_dir, '_basemap') #make output directory
dir.check(chm_change_dir_bm)

pblapply(1:length(change_rasts_clipped), function(i){
  n = names(change_rasts_clipped)[i]
  filename = paste0(chm_change_dir_bm,'/',n, '.tif')
  if(!file.exists(filename)){
    chm = change_rasts_clipped[[n]]
    ps = bm_rasts[[n]]
    chm_p = project(chm, ps)
    writeRaster(chm_p, filename)
  }
})

#----save change validation rasters for blocks which are resampled to match scenes----

#get first file from each subdirectory so that there's one planetscope raster for each thinning block
scene_dir = dirs[str_detect(dirs,'Z_')]
scene_subdirs = list.dirs(scene_dir)
scene_subdirs = scene_subdirs[scene_subdirs != scene_dir]
scene_subdirs = scene_subdirs[!str_detect(scene_subdirs, 'NoChange')]

# #resample based on a single scene which fully contains the AOI
# scene_r = rast(list.files(scene_subdirs, recursive = T, full.names = T, pattern='\\.tif$')[1])
# if(terra::relate(ext(blocks), ext(scene_r), 'within')[1,1]){
#   
#   #make output directory
#   chm_change_dir_scene = paste0(chm_change_dir, '_scenes_',target_res)
#   dir.check(chm_change_dir_scene)
#   
#   #resample chm_change rasters to match scenes, save
#   pblapply(1:length(change_rasts_clipped), function(i){
#     n = names(change_rasts_clipped)[i]
#     filename = paste0(chm_change_dir_scene,'/',n, '.tif')
#     if(!file.exists(filename)){
#       chm = change_rasts_clipped[[n]]
#       ps = scene_r
#       chm_p = project(chm, ps)
#       writeRaster(chm_p, filename)
#     }
#   })
# }

scene_rasts = lapply(scene_subdirs, function(x){
  fl = list.files(x, full.names = T)
  f = fl[1]
  print(f)
  r = rast(f)
  return(r)
})
scene_ids = basename(scene_subdirs)
names(scene_rasts) = scene_ids

#make output directory
chm_change_dir_scene = paste0(chm_change_dir, '_scenes_Res=',target_res)
dir.check(chm_change_dir_scene)

#resample chm_change rasters to match scenes, save
pblapply(1:length(change_rasts_clipped), function(i){
  n = names(change_rasts_clipped)[i]
  filename = paste0(chm_change_dir_scene,'/',n, '.tif')
  if(!file.exists(filename)){
    chm = change_rasts_clipped[[n]]
    ps = scene_rasts[[n]]
    chm_p = project(chm, ps)
    writeRaster(chm_p, filename)
  }
})


#reload resampled change rasters
chm_changes_scene = lapply(list.files(chm_change_dir_scene, full.names = T), rast)
names(chm_changes_scene) = file_path_sans_ext(list.files(chm_change_dir_scene))

#create dataframe of values
change_scene_df_l = pblapply(1:length(chm_changes_scene), function(i){
  d = data.table(
    vals = values(chm_changes_scene[[i]])
    ,
    block = names(chm_changes_scene)[i]
  )
})
change_scene_df = data.table::rbindlist(change_scene_df_l)

#get otsu thresholds
otsu_scene_l = pbsapply(names(chm_changes_scene), function(x){
    d = change_scene_df %>%
      filter(block == x) 
    t = otsuThresholdCpp(values = d[['vals.Z']], bins = 256)
  })

#export otsu thresholds for planetscope data with target_resolution
otsu_scene_tbl = tibble(
  block_id = names(otsu_scene_l),
  chm_change_otsu_threshold = otsu_scene_l
)

otsu_scene_file = paste0(dir,'/','Otsu_change_thresholds_Res=',target_res,'.csv')
if(!file.exists(otsu_scene_file)){write_csv(otsu_tbl,otsu_scene_file)}

#add otsu thresholds to dataframe
change_scene_df = change_scene_df %>% 
  left_join(data.frame(block = names(otsu_scene_l), otsu_threshold = otsu_scene_l), by = 'block')

# #plot density graphs with otsu thresholds
# density_otsu_plot_3m = ggplot(data = change_scene_df)+
#   geom_density(aes(x = vals.Z),fill='seagreen')+
#   geom_vline(aes(xintercept = otsu_threshold
#                  # , linetype = 'Otsu threshold'
#   )
#   ,color='darkred'
#   )+
#   # annotate("text", x = otsu_t+1, y = 0.16, label = paste0("Thinning threshold: ",otsu_t)
#   #          , angle = 0, vjust = -0.5, hjust = 0, color = 'darkred')+
#   ylab('Density')+
#   theme_classic()+
#   xlab("Height change (m)")+
#   theme(legend.position = NULL)+
#   facet_grid(vars(block))

#combine 0.25m and 3.76m data, make combined plot
change_df_combined = rbind(
  change_df |> mutate(resolution = '0.25m'),
  change_scene_df |> mutate(resolution = '3.76m')
)

#remove values which are many standard deviations removed from the rest of the data to make plotting and plot generation better
change_df_combined_ = change_df_combined |> filter(vals.Z <= 10 & vals.Z >= -20)

#get subset of data to make generating plots much faster
# set.seed(123)
# change_df_combined_idx = sample(1:nrow(change_df_combined_), size = 0.05*nrow(change_df_combined_))
# change_df_combined_ = change_df_combined_[change_df_combined_idx,]

#make smaller dataset for plotting rectangles and text
rect_text_data = change_df_combined_ |>
  distinct(block, resolution, otsu_threshold)

#make density plots for each block at 0.25m and planetscope scene resolution
{
# ggplot(data = change_df_combined_ #|>filter(resolution=='3.76m') #filter used to tune plot parameters, comment out for final plot
#        )+
#   #color plot area
#   geom_rect(data = rect_text_data, aes(xmin = -Inf, xmax = otsu_threshold, ymin = -Inf, ymax = Inf, 
#                 fill = "Thinned"), alpha = 1)+
#   geom_rect(data = rect_text_data,  aes(xmin = otsu_threshold, xmax = Inf, ymin = -Inf, ymax = Inf, 
#                 fill = "Not thinned"), alpha = 0.5)+
#   scale_fill_manual(values = c("Thinned" = "lightyellow", "Not thinned" = "lightblue"),
#                     name = NULL)+
#   #add density plot
#   geom_density(aes(x = vals.Z),fill='seagreen')+
#   #add vertical line and text label for thinning threshold
#   geom_vline(data = rect_text_data, aes(xintercept = otsu_threshold
#                  , linetype = 'Thinning threshold'
#   )
#   ,color='darkred'
#   )+
#   scale_linetype_manual(values = c("Thinning threshold" = "solid"), 
#                         name = NULL)+
#   geom_text(data = rect_text_data, aes(x = otsu_threshold-2.5, label = round(otsu_threshold, 2)), 
#             y = Inf, vjust = 1.2, hjust = 0.5, color = 'darkred', 
#             # angle = 90, 
#             size = 3)+
#   ylab('Density')+
#   theme_classic()+
#   xlab("Height change (m)")+
#   # theme(legend.position = NULL)+
#   facet_grid(rows=vars(block), cols = vars(resolution))
# 
# pfile = 'figures/density_plots_0.25m_3.76m.png'
# if(!file.exists(pfile)){ggsave(filename = pfile, dpi=600)}
  }


#----make road and non-vegetated ground mask----

blocks_scene_p = project(blocks, crs(scene_rasts[[1]])) #reproject blocks to match scene crs
blocks_bm_p = project(blocks, crs(bm_rasts[[1]]))

pre_threshold = 0.5
post_threshold = 0.3
change_threshold = 0.1

#combine scene rasters
combined_lid_rasts_l = pblapply(blocks_, function(x){
  
  block = blocks_scene_p[blocks_p$BLOCKNUM==x,]
  # scene = scene_rasts[[x]]
  
  # print(paste('Processing', block$BLOCKNUM,'combining and cropping lidar rasters'))
  pre = pre_rasts_l[[x]]
  post = post_rasts_l[[x]]
  change = change_rasts_l[[x]]
  
  pl = list(pre,post,change)
  names(pl) = c('pre','post','change')
  
  #align rasters and project to match blocks
  r_l = lapply(names(pl), function(y){
    r = project(pl[[y]], crs(block))
    # print(paste('Projected',y))
    r = crop(r,block,mask=T)
    # print(paste('Cropped',y))
    return(r)
  })
  names(r_l) = names(pl)
  
  r_l[['post']] = resample(r_l[['post']], r_l[['change']])
  # print(paste('Resample post'))
  r_l[['pre']] = resample(r_l[['pre']], r_l[['change']])
  # print(paste('Resample pre'))
  
  c_r = rast(r_l)
  # print(paste('Created combined raster'))
  
  names(c_r) = c('pre', 'post', 'change')
  return(c_r)
})
names(combined_lid_rasts_l) = blocks_

combined_lid_dir = paste0(dir,'/chm_combined_cropped') #save
dir.check(combined_lid_dir)
pblapply(names(combined_lid_rasts_l), function(n){
  filename=paste0(combined_lid_dir,'/',n,'.tif')
  if(!file.exists(filename)){
    writeRaster(combined_lid_rasts_l[[n]],filename)
  }
})

#make masks for scenes
nonveg_masks_l = pblapply(combined_lid_rasts_l, function(r){
  m = ifel(r[['pre']] < pre_threshold & r[['post']] < post_threshold
           ,NA
           ,1)
  return(m)
})
names(nonveg_masks_l) = blocks_
# list_plot(scene_nonveg_masks_l)

nonveg_mask_dir = paste0(dir,'/nonveg_mask_pre=',pre_threshold,'_post=',post_threshold) #save masks
dir.check(nonveg_mask_dir)
pblapply(names(nonveg_masks_l), function(n){
  filename=paste0(nonveg_mask_dir,'/',n,'.tif')
  if(!file.exists(filename)){
    writeRaster(nonveg_masks_l[[n]],filename)
  }
})

#resample to match scene resolution
resample_mask = pblapply(blocks_, function(x){
  scene = scene_rasts[[x]]
  m = scene_nonveg_masks_l[[x]]
  
  mrs = resample(m, scene)
  return(mrs)
})

#----make map of chm and PS data----

library(tmap)
library(tidyterra)

#load combined lidar data
combined_lid_files = list.files(combined_lid_dir, full.names=T)

blocks_v = vect(blocks)
block_12nt3 = blocks_v[blocks_v$BLOCKNUM=='12N_T3']

bs = '12N_T3'


combined_lid_files = list.files(combined_lid_dir, full.names=T)

comb_lid_12nt3 = rast(combined_lid_files[str_detect(combined_lid_files,bs)])

names(comb_lid_12nt3) = c('July 2021 canopy height', 'October 2024 canopy height', 'Height change')

ggplot()+
  geom_spatraster(data = comb_lid_12nt3)+
  facet_wrap(~lyr)+
  scale_fill_viridis_c(option='viridis', na.value = 'white', name = 'Value')+
  theme_classic()+
  # coord_equal(expand = FALSE) +
  theme(
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank()
  )

# scene_nonveg_mask_dir = paste0(dir,'/nonveg_mask_scene_pre=',pre_threshold,'_post=',post_threshold)
# dir.check(scene_lid_dir)
# pblapply(names(scene_combined_rasts_l), function(n){
#   filename=paste0(scene_lid_dir,'/',n,'.tif')
#   if(!file.exists(filename)){
#     writeRaster(scene_combined_rasts_l[[n]],filename)
#   }
# })

# #check results
# scene_masked_rasts_l = pblapply(1:length(scene_nonveg_masks_l), function(i){
#   nvm = scene_nonveg_masks_l[[i]]
#   r = scene_combined_rasts_l[[i]]
#   r_m = mask(r, nvm)
#   return(r_m)
# })
# names(scene_masked_rasts_l) = blocks_
# pre_l = lapply(scene_masked_rasts_l, function(x)x[[1]])
# names(pre_l) = names(scene_masked_rasts_l)
# list_plot(pre_l)
# # list_plot(masked_rasts_l, single_plot = T)
# 
# #combine basemap rasters
# bm_combined_rasts_l = pblapply(blocks_, function(x){
#   
#   block = blocks_bm_p[blocks_p$BLOCKNUM==x,]
#   print(paste('Processing', block$BLOCKNUM,'combining and cropping lidar rasters'))
#   pre = pre_rasts_l[[x]]
#   post = post_rasts_l[[x]]
#   change = change_rasts_l[[x]]
#   
#   pl = list(pre,post,change)
#   names(pl) = c('pre','post','change')
#   
#   #align rasters and project to match blocks
#   r_l = pblapply(names(pl), function(y){
#     r = project(pl[[y]], crs(block))
#     print(paste('Projected',y))
#     r = crop(r,block,mask=T)
#     print(paste('Cropped',y))
#     return(r)
#   })
#   names(r_l) = names(pl)
#   
#   r_l[['post']] = resample(r_l[['post']], r_l[['change']])
#   print(paste('Resample post'))
#   r_l[['pre']] = resample(r_l[['pre']], r_l[['change']])
#   print(paste('Resample pre'))
#   
#   c_r = rast(r_l)
#   print(paste('Created combined raster'))
#   
#   names(c_r) = c('pre', 'post', 'change')
#   return(c_r)
# })
# names(bm_combined_rasts_l) = blocks_
# 
# #make masks for basemaps
# bm_nonveg_masks_l = pblapply(bm_combined_rasts_l, function(r){
#   m = ifel(r[['pre']] < pre_threshold & r[['post']] < post_threshold
#            ,NA
#            ,1)
#   return(m)
# })
# names(bm_nonveg_masks_l) = blocks_
# list_plot(bm_nonveg_masks_l)
