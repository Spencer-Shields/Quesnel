library(tidyverse)
library(terra)
library(RStoolbox)
library(tools)
source('helper_functions.R')
library(BAMMtools)
library(pbapply)
library(microbenchmark)
library(future)
library(future.apply)
library(tidyterra)
library(sf)

#---- Load data ----
#load landsat-8 scene
l8_dir = 'data/landsat'
l8_scene = landstack(l8_dir)[[1]]

#load near-coincident Planetscope scene
ps_dir = list.dirs('data/planet_scenes/raw')[3]
ps_scene = rast(paste0(ps_dir,"/20240831_192428_37_24f6_3B_AnalyticMS_SR_8b_harmonized_clip_reproject.tif"))

#---- calculate L8 tasseled cap indices ----

l8_tc_coefs = data.frame(
  bands = c('B2','B3','B4','B5','B6','B7'),
  brightness = c(0.3029, 0.2786, 0.4733, 0.5599, 0.508, 0.1872),
  greenness = c(-0.2941, -0.243, -0.5424, 0.7276, 0.0713, -0.1608),
  wetness = c(0.1511, 0.1973, 0.3283, 0.3407, -0.7117, -0.4559)
)

l8_scene_sub = l8_scene[[2:7]]

tc_band_names = c('brightness', 'greenness', 'wetness')

tc_list = pblapply(tc_band_names, function(n){
  coefs = l8_tc_coefs[[n]]
  l = pblapply(1:length(coefs), function(i){
    r = l8_scene_sub[[i]]*coefs[i]
    r
  })
  rtc = app(rast(l), 'sum')
  rtc
})

l8_tc = rast(tc_list)
names(l8_tc) = tc_band_names

#---- model tasseled cap indices using planetscope data ----

#get site bounding box, reproject to match l8, crop l8
ex = ext(ps_scene)
aoi_v = as.polygons(ex)
crs(aoi_v) = crs(ps_scene)
aoi_p = project(aoi_v, crs(l8_tc))
l8_tc_crop = crop(l8_tc, aoi_p)


#reproject and resample ps to l8 geometry
if(crs(l8_tc_crop)!=crs(ps_scene)){
  ps_r = project(ps_scene, l8_tc_crop)
}


#stack all, generate models
combined_r = rast(list(l8_tc_crop, ps_r))
names(combined_r) #check band names

combined_v = as.data.frame(na.omit(values(combined_r)))

lm_l = pblapply(tc_band_names, function(n){
  d = combined_v %>% select(all_of(c(n, names(ps_r))))
  m = lm(as.formula(paste(n, "~ .")), data = d)
})

lm_summs_l = lapply(lm_l, summary) #model summary stats

#----export data----
ps_tc_band_names = paste0('ps_', tc_band_names)

ps_tc_dir = 'data/ps_TasseledCap_models'
dir.check(ps_tc_dir)

pblapply(1:length(lm_l), function(i){
  filename = paste0(ps_tc_dir,'/',ps_tc_band_names[i],'_model.rds')
  if(!file.exists(filename)){
    write_rds(x = lm_l[[i]], file = filename)
  }
})

coefs_df = rbind(lm_l[[1]]$coefficients, lm_l[[2]]$coefficients, lm_l[[3]]$coefficients) %>% 
  as.data.frame()
coefs_df = cbind(ps_tc_band_names, coefs_df)
filename = paste0(ps_tc_dir, '/PS_TC_lm_coefficients.csv')
if(!file.exists(filename))(write.csv(coefs_df, filename))

#apply models to other rasters
ps_tc_l = pblapply(lm_l, function(m){predict(ps_scene, m)})

ps_tc = rast(ps_tc_l)
names(ps_tc) = paste0('ps_', tc_band_names)


ps_di = di_fn(ps_tc) #calculate disturbance index


di_v = as.data.frame(na.omit(values(ps_di))) %>% mutate(znormed = 'Yes')
x11()
ggplot()+
  # geom_histogram(data = di_v, aes(x = DI), bins = 100)+
  geom_density(data = di_v, aes(DI))+
  xlab('Disturbance index')

blocks_p = st_transform(blocks, crs = crs(ps_di))

x11()
ggplot() +
  geom_spatraster(data = ps_di) +
  scale_fill_viridis_c(name = "Disturbance Index") +  # Rename raster legend
  geom_sf(data = blocks_p, aes(color = BLOCKNUM), fill = NA) +
  scale_color_discrete(name = "Block ID Number") +  # Rename BLOCKNUM legend
  theme_minimal()

#test against non-normalized disturbance index

ps_di_nonz = di_fn(ps_tc, znormed = F)
di_v_nonz = as.data.frame(na.omit(values(ps_di_nonz))) %>% mutate(znormed = 'No')

di_df_combined = rbind(di_v, di_v_nonz)

ggplot(data = di_df_combined, aes(DI))+
  geom_density()+
  facet_wrap(facets = vars(znormed), scales = 'free')

vars = di_df_combined %>%
  group_by(znormed)%>%
  summarise(coefficient_of_variation = sd(DI)/mean(DI),
            sd = sd(DI))

#test orthogonality by calculating pairwise correlation coefficients
{
  # ps_cor_matrix = layerCor(ps_tc, fun = 'pearson') #indices are not really orthogonal
  # l8_cor_matrix = layerCor(l8_tc, fun = 'pearson') #neither are these
  }
#--- try another approach, use PCA ----

# ps_v_pc = terra::prcomp(ps_scene, scale. = T)
# summary(ps_pc)
# 
# ps_pc = predict(ps_scene, ps_v_pc)
# 
# pspc_cor_mat = layerCor(ps_pc, fun='pearson')
# 
# x11()
# plot(ps_pc[[1:3]])




#---- clustering ----

# min_clusts = 2
# max_clusts = 50
# 
# ps_di_df = as.data.frame(na.omit(values(ps_di)))
# 
# k_2 = kmeans(ps_di_df, centers = 2, iter.max = 20)
# k_2

# km_l = pblapply(min_clusts:max_clusts, function(i){
#   kmeans(ps_di_df, i)
# })
# 
# k_results = data.frame(
#   wcss = sapply(km_l, function(k)k[['tot.withinss']]),
#   n_clusts = min_clusts:max_clusts)
# 
# ggplot(k_results, aes(x = n_clusts, y = wcss))+
#   geom_point()

# ps_di_v = as.numeric(na.omit(values(ps_di)))
# j = getJenksBreaks(ps_di_v, k = 1)
# #classify image based on jenks threshold
# ps_di_c = ifel(ps_di >= j, 1, 0)

# ggplot() +
#   # geom_histogram(data = di_v, aes(x = DI), bins = 100) +
#   geom_density(data = di_v, aes(DI))+
#   geom_vline(aes(xintercept = j), linetype = "dashed", color = "red") +
#   annotate("text", x = j, y = Inf, label = "Threshold", vjust = 1.5, color = "red")


# #write a function for running Otsu's method on a numerical vector
# otsu = function(x, threads = NULL, threshold_increment = NULL){
#   #x is a numerical vector.
#   #threshold_increment is an optional parameter to set the increment for what threshold values will be tested. Making it larger will reduce computational requirements
#   
#   #sort the vector, generate sequence of threshold values
#   x = sort(x)
#   min_x = x[1]
#   max_x = x[length(x)]
#   
#   #generate sequence of threshold values
#   if(is.null(threshold_increment)){
#     threshold_increment = min(diff(x))
#     
#     seq_length = (max_x-min_x)/threshold_increment
#     max_seq_length = 10*length(x)
#     
#     while(seq_length > max_seq_length){
#       threshold_increment = 10*threshold_increment
#       seq_length = (max_x-min_x)/threshold_increment
#     }
#   }
#   
#   t_seq = seq(from = min_x, to = max_x, by = threshold_increment)
#   
#   if(!is.null(threads)){future::plan(multisession)}
#   
#   pbsapply(t_seq, function(n){
#     
#     #subset vector using threshold
#     class1 = x[x<=n]
#     class2 = x[x>n]
#     
#     #calculate variance of each class
#     var1 = var(class1)
#   })
#   
#   
# }
# 
# microbenchmark(sort(x),min(x),x[1],length(x),seq(from = min_x, to = max_x))

#calculate NDWI

ps_ndwi = (ps_scene[['green']] - ps_scene[['nir']])/
  (ps_scene[['green']] + ps_scene[['nir']])
plot(ps_ndwi)
density(ps_ndwi)
