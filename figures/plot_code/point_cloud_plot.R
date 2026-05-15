source('scene_global_stats_filtering_20251106.R')
library(lidR)


#load post-harvest points

post_harvest_las = pblapply(post_blockdirs, function(x)readLAScatalog(file.path(x, "input/las/norm"), filter = '-keep_first'))
names(post_harvest_las) = basename(post_blockdirs)

target_block = '12L_D345'

lc = post_harvest_las[target_block]
plot(lc)
lc_bb = st_bbox(lc)
lc_x_center <- mean(c(lc_bb["xmin"], lc_bb["xmax"]))
lc_y_center <- mean(c(lc_bb["ymin"], lc_bb["ymax"]))

lc_sample = clip_circle(lc, lc_x_center, lc_y_center, 60)
plot(lc_sample)

#load pre-harvest points
pre_harvest_las = pblapply(pre_blockdirs, function(x)readLAScatalog(file.path(x, "input/las/norm"), filter='-keep_first'))
names(pre_harvest_las) = basename(pre_blockdirs)

lc = pre_harvest_las[[target_block]]
lc_sample = clip_circle(lc, lc_x_center, lc_y_center, 60)
plot(lc_sample)