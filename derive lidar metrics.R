# source('scene_setup_preprocessing_20250909.R')
source('scene_global_stats_filtering_20250925.R')
library(lidR)
library(lasR)


parent_dir = 'data/quesnel_thinning_las'

pre_dir = paste0(parent_dir,'/pre')
pre_blockdirs = list.dirs(pre_dir, recursive=F)

post_dir = paste0(parent_dir,'/post')
post_blockdirs = list.dirs(post_dir, recursive=F)

#----inspect point clouds----

pre_harvest_las = pblapply(pre_blockdirs, function(x)readLAScatalog(file.path(x, "input/las/norm")))
names(pre_harvest_las) = basename(pre_blockdirs)

post_harvest_las = pblapply(post_blockdirs, function(x)readLAScatalog(file.path(x, "input/las/norm")))
names(post_harvest_las) = basename(post_blockdirs)

#----process point clouds, generate lidar metrics----

#set parameters for building pipelines
res = target_res #defined in the initial pre-processing script

metrics = c(
  #height
  'z_max', 'z_min', paste0('z_p',seq(5,95,5)), 'z_mean'
  #canopy cover
  , paste0('z_above',seq(1,5,1))
  #height variability
  ,'z_sd', 'z_cv'
  #number of points in each pixel
  ,'count'
  )

# #numer of first returns for pre_harvest data (used to get target sample sizes for post-harvest rasters)
# pre_harv_firsts_counts = pblapply(pre_harvest_las, function(lc){
#   set_parallel_strategy(concurrent_files(6))
#   r = lasR::exec(pipeline = reader()+rasterize(res,'count',filter=keep_first()), on = lc, with = list(progress = T))
#   return(r)
# })
# 
# sample_size_rasts = lapply(pre_harv_firsts_counts, function(r)global(r, c('mean', 'sd'), na.rm=T)) |> bind_rows() |> mutate(block = names(pre_harv_firsts_counts))


#generate pre_harvest lidar metrics
# pre_harv_pipe = reader() + normalize() + rasterize(res, operators = metrics, filter=keep_first()) #don't need to normalize since that's already been done
pre_harv_pipe = reader() + rasterize(res, operators = metrics, filter=keep_first())


# pre_harvest_metrics = pblapply(pre_harvest_las[1], function(lc){
#   set_parallel_strategy(concurrent_files(6))
#   r = lasR::exec(pre_harv_pipe, on = lc, with = list(progress = T))
#   return(r)
# })

tic()
pre_harvest_metrics = lapply(pre_blockdirs, function(lc){
  lc_ = paste0(lc,'/input/las/norm')
  set_parallel_strategy(concurrent_files(6))
  r = lasR::exec(pre_harv_pipe
                 , on = lc_, with = list(progress = T))
  return(r)
})
names(pre_harvest_metrics) = basename(pre_blockdirs)
toc()


#see how post-harvest metrics look without decimating the pointclouds
# post_harvest_metrics = pblapply(post_harvest_las, function(lc){
#   set_parallel_strategy(concurrent_files(12))
#   r = lasR::exec(pre_harv_pipe, on = lc, with = list(progress = T))
#   return(r)
# })

# tic()
post_harvest_metrics = pblapply(post_blockdirs, function(lc){
  lc_ = paste0(lc,'/input/las/norm')
  set_parallel_strategy(concurrent_files(6))
  r = lasR::exec(pre_harv_pipe
                 , on = lc_, with = list(progress = T))
  return(r)
})
names(post_harvest_metrics) = basename(post_blockdirs)
# toc()


