library(tidyverse)
library(terra)
library(tools)
library(pbapply)
library(Rcpp)
source('helper_functions.R')
library(tidyterra)
source('https://raw.githubusercontent.com/Spencer-Shields/planetscope_radiometric_correction/refs/heads/main/check_radiometric_consistency.R')
library(data.table)
# library(supercells)
# library(dbscan)

#----Load data----
dir = 'data/Quesnel_thinning'
d_l = list.dirs(dir, full.names = T)

chm_post_dir = d_l[str_detect(d_l, 'post')]
chm_pre_dir = d_l[str_detect(d_l, 'pre')]

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

chm_change_dir = 'data/Quesnel_thinning/chm_change'
dir.check(chm_change_dir)

lapply(1:length(change_rasts_l), function(i){
  filename = paste0(chm_change_dir,'/',chm_common_ids[i],'.tif')
  if(!file.exists(filename)){
    writeRaster(change_rasts_l[[i]], filename)
  }
})


#----Look at histograms of difference rasters to establish harvest threshold----

change_rasts_clipped = pblapply(1:length(change_rasts_l), function(i){
  rast_id = names(change_rasts_l)[i]
  block = blocks %>% filter(BLOCKNUM == rast_id)
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

change_df = change_df %>% 
  left_join(data.frame(block = names(otsu_l), otsu_threshold = otsu_l), by = 'block')

ggplot(data = change_df)+
  geom_density(aes(x = vals.Z))+
  geom_vline(aes(xintercept = otsu_threshold, linetype = 'Otsu threshold'))+
  xlab("Height change (m)")+
  facet_grid(vars(block))

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


# x11()
# plot(pre_rasts_common_l[[1]], main = 'pre')
# x11()
# plot(post_rasts_common_l[[1]], main = 'post')
# x11()
# plot(change_rasts_l[[1]], main = 'pre minus post')

# ch = change_rasts_l[[1]]
# ch_df = as.data.frame(values(ch))
# names(ch_df) = 'val'
# 
# 
# h=hist(ch, breaks=30, maxcell = 3000000, plot=F)
# h_d = data.frame(
#   # breaks = h$breaks, 
#                  counts = h$counts, 
#                  mids = h$mids)
# 
# ggplot()+
#   geom_histogram(data = ch_df, aes(val), binwidth = 2)
# 
# ggplot(ch_df, aes(val))+
#   geom_bar()+
#   scale_x_binned()

# #----Make difference raster masks----
# 
# height_change_threshold = 3
# 
# harvested_l = pblapply(change_rasts_l, function(x){
#   msk = ifel(x >= height_change_threshold, 1, NA)
#   masked = mask(x, msk)
#   masked
# })
# names(harvested_l) = chm_common_ids
# 
# list_plot(harvested_l)
# list_plot(pre_rasts_common_l)
# list_plot(post_rasts_common_l)
# # list_plot(harvested_l, single_plot = F)

# #----segmentation----
# 
# a = change_rasts_l[[1]]
# a
# density(a)
# 
# library(future)
# plan(multisession, workers = 8)
# # slic_a = supercells(a
# #                     # , k = round(ncell(a)/32)
# #                     ,1000000
# #                     , compactness = 1
# #                     # ,chunks = 150
# #                     , future = T
# #                     , verbose = 1)
# # 
# # ggplot()+
# #   geom_sf(data = slic_a, aes(fill = Z), color = NA)+
# #   scale_fill_viridis_b()
# plan(sequential)
# 
# 
# 
# #kmeans
# min_clusts = 2
# max_clusts = 50
# 
# # k = terra::k_means(a, 2)
# a_df = as.data.frame(na.omit(values(a)))
# 
# # k_2 = kmeans(ps_di_df, centers = 2, iter.max = 20)
# # k_2
# # 
# plan(multisession, workers = 8)
# km_l = pblapply(min_clusts:max_clusts, function(i){
#   kmeans(a_df, i)
# }
# ,cl = 'future')
# plan(sequential)
# 
# k_results = data.frame(
#   wcss = sapply(km_l, function(k)k[['tot.withinss']]),
#   n_clusts = min_clusts:max_clusts)
# 
# ggplot(k_results, aes(x = n_clusts, y = wcss))+
#   geom_point()
# 
# k = terra::k_means(a,20)
# plot(k)
# # ps_di_v = as.numeric(na.omit(values(ps_di)))
# # j = getJenksBreaks(ps_di_v, k = 1)
# # #classify image based on jenks threshold
# # ps_di_c = ifel(ps_di >= j, 1, 0)
# # 
# # ggplot() +
# #   # geom_histogram(data = di_v, aes(x = DI), bins = 100) +
# #   geom_density(data = di_v, aes(DI))+
# #   geom_vline(aes(xintercept = j), linetype = "dashed", color = "red") +
# #   annotate("text", x = j, y = Inf, label = "Threshold", vjust = 1.5, color = "red")
