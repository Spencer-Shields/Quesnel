library(tidyverse)
library(terra)
library(RStoolbox)
library(tools)
library(sf)
source('helper_functions.R')
library(pbapply)
library(parallel)
library(future)
library(future.apply)

#set up number of cores to be used for parallel processing
cores = detectCores()

#----load and examine basemap data----
bm_dir = 'data/planet_basemaps/global_monthly'

raw_dir = paste0(bm_dir, '/raw')
# indir = 'data/planet_basemaps/ps_monthly_sen2_normalized_analytic_8b_sr/raw'
# indir = 'data/planet_basemaps/global_quarterly/raw'

rast_files = list.files(raw_dir, full.names = T, recursive = T
                        , pattern = '\\.tif$'
)
r_l = lapply(rast_files,rast) #get rasters


#----get thinning blocks and aoi----
ps_dir = 'data/planet_scenes'
blocks = st_read('data/Quesnel_thinning/12l_12n_bdy.gpkg')
no_change = st_read('Arc project/NoChangeStands_conifer.shp')
no_change = no_change %>% rename(geom = geometry)
blocks = rbind(blocks, no_change)
blocks_p = st_transform(blocks, crs = crs(r_l[[1]]))

aoi = st_read('data/AOI_fullsite_wgs84.kml') %>%
  st_transform(crs = crs(r_l[[1]]))

# ggRGB(r_l[[1]])+
#   geom_sf(data = blocks_p, aes(color = BLOCKNUM))
# 
# ggRGB(r_l[[2]]) +
#   geom_sf(data = blocks_p, aes(color = BLOCKNUM))

#----mosaic basemap quads from each month, crop to aoi----

#get subdirectories that contain monthly basemap quads
month_dirs = list.dirs(raw_dir, full.names = T)
month_dirs = month_dirs[month_dirs != raw_dir] #exclude parent directory

desired_months = c('01','02','03','04', '05','06', '07', '08', '09', '10', '11','12')
month_str = paste0('_',desired_months,'_')
month_pattern = paste0("_(", paste(desired_months, collapse = "|"), ")_")

desired_month_dirs = month_dirs[str_detect(month_dirs, month_pattern)]

mosaic_string = 'cropped_month_mosaics'
outdir = paste0(bm_dir, '/', mosaic_string)
dir.check(outdir)


# cl = makeCluster(ceiling(cores/3))
# plan('cluster', workers = cl)

future_lapply(desired_month_dirs, function(d){
  output_file = paste0(outdir, '/', basename(d),'.tif')
  if(!file.exists(output_file)){
    f = list.files(d, full.names = T, recursive = T, pattern = '\\.tif$')
    rl = lapply(f,rast)
    sc = sprc(rl)
    m = mosaic(sc)
    mc = crop(m,aoi)
    writeRaster(mc, output_file)
  }
})
print('Mosaic crop done')
# stopCluster(cl)
# plan('sequential')

#----calculate vegetation indices----

indir = paste0(bm_dir,'/',mosaic_string)
files_list = list.files(indir, full.names = T, recursive = T, pattern = '\\.tif$')

vi_string = 'VegIndices'
outdir = paste0(bm_dir,'/', vi_string)
dir.check(outdir)

# cl = makeCluster(ceiling(cores/4))
# plan('cluster', workers = cl)

future_lapply(files_list, function(x){
  output_file = paste0(outdir, '/', basename(x),'.tif')
  if(!file.exists(output_file)){
    
    #load raster
    r = rast(x)
    
    #calculate vegetation indices, include original raster bands in the output
    vi = rast.batch.functions(r, include_input = T, fl = c(
      hue.vi,
      gcc.vi,
      ndgr.vi,
      bi.vi,
      ci.vi,
      cri550.vi,
      gli.vi,
      tvi.vi,
      varig.vi
    ))
    names(vi[[1:4]]) = c('blue', 'green', 'red', 'max_DN') #rename original bands
    
    #save raster
    writeRaster(vi, output_file)
    gc()
  }
})

print('Calculate VegIndices done')
# stopCluster(cl)
# plan('sequential')

#----crop data to blocks----

indir = outdir
rasters = list.files(indir, pattern = '\\.tif$', full.names = T)

clipped_string = 'VegIndices_clipped'
clipped_dir = paste0(bm_dir, '/', clipped_string)
dir.check(clipped_dir)

outdir = clipped_dir

lapply(1:nrow(blocks_p), function(i){ #create subdirectories for every block
  d = paste0(outdir,'/',blocks_p$BLOCKNUM[i])
  dir.check(d)
})

pblapply(1:length(rasters), function(i){
  r_id = basename(file_path_sans_ext(rasters[i]))
  future_lapply(1:nrow(blocks_p), function(j){
    block_id = blocks_p$BLOCKNUM[j]
    filename = paste0(outdir,'/',block_id,'/',r_id)
    print(paste0('Processing ',filename))
    if(!file.exists(filename)){
      r = rast(rasters[i])
      b = blocks_p[j,]
      r_c = crop(r,b,mask=T)
      terra::writeRaster(r_c, filename)
    }
  })
})





#----Calculate z_scores and softmax using cropped data----

z_dir_string = 'Z-scoreCalcByBlock'
z_dir = paste0(bm_dir, '/', z_dir_string)
dir.check(z_dir)

lapply(1:nrow(blocks_p), function(i){ #create subdirectories for every block
  d = paste0(z_dir,'/',blocks_p$BLOCKNUM[i])
  dir.check(d)
})

sm_dir_string = 'SoftmaxCalcByBlock'
sm_dir = paste0(bm_dir, '/', sm_dir_string)
dir.check(sm_dir)

lapply(1:nrow(blocks_p), function(i){ #create subdirectories for every block
  d = paste0(sm_dir,'/',blocks_p$BLOCKNUM[i])
  dir.check(d)
})

instring = 'nonorm-blocks'
indir = paste0(bm_dir, '/',instring) #this directory and data came from first running the basemap_timeseries script
rasters = list.files(indir, full.names = T, recursive = T, pattern = '\\.tif$')

sf = 255 #scale factor to use for calculating softmax (i.e. maximum possible band pixel value)

z_outdir = z_dir
sm_outdir = sm_dir

bands = c('blue', 'green', 'red')

# cl = makeCluster(ceiling(cores/4))
# plan('cluster', workers = cl)

#calculate z-scores on cropped data
future_lapply(1:length(rasters), function(i){
  filename = str_replace(rasters[i], instring, z_dir_string) 
  if(!file.exists(filename)){
    #load scene
    scene = rast(rasters[i])
    
    #calculate and save z-score
    scene_z_l = lapply(scene, z_rast)
    scene_z = rast(scene_z_l)
    names(scene_z) = paste0(names(scene),'_',z_dir_string)
    terra::writeRaster(scene_z, filename, overwrite = T)
  }
})

#calculate softmax on cropped data
future_lapply(1:length(rasters), function(i){
  filename = str_replace(rasters[i], instring, sm_dir_string) 
  if(!file.exists(filename)){
    #load scene
    scene = rast(rasters[i])
    
    #calculate and save softmax
    
    s_l = lapply(1:nlyr(scene), function(i){ #scale the raw bands for use with softmax
      b = scene[[i]]
      nom = names(scene)[i]
      if(nom %in% bands){
        b = b/sf
      } else {
        b 
      }
      return(b)
    })
    scene_s = rast(s_l)
    
    scene_sm = softmax(scene_s, append_name = T)
    terra::writeRaster(scene_sm, filename, overwrite = T)
  }
})

# stopCluster(cl)
# plan('sequential')

#----Calculate difference between time t and time t-1 for whole time series----

#----z_score delta

#create output dir
dz_string = paste0(z_dir_string,'_delta')
dz_dir = paste0(bm_dir,'/',dz_string)
dir.check(dz_dir)

lapply(1:nrow(blocks_p), function(i){ #create subdirectories for every block
  d = paste0(dz_dir,'/',blocks_p$BLOCKNUM[i])
  dir.check(d)
})

indir = z_dir #input files

in_subdirs = list.dirs(indir)
in_subdirs = in_subdirs[in_subdirs != indir]

#calculate difference rasters
lapply(in_subdirs, function(d){
  #load input block subdirectory
  d_id = basename(d)
  files = list.files(d, full.names = T, pattern = '\\.tif$')
  df = tibble(
    files
  ) %>%
    #get dates from filepaths
    mutate(year_month = str_extract(files, "\\d{4}_\\d{2}")) %>%
    mutate(year_month_day = paste0(year_month,'_01')) %>%
    mutate(date = as.Date(str_replace_all(year_month_day,'_','-'))) %>%
    #sort filepaths by increasing order of date
    arrange(date)
  
  files_sorted = df$files
  
  #calculate rasters from files in subdirectory
  pblapply(2:length(files_sorted), function(i){ #start at 
    
    f_id = file_path_sans_ext(basename(files[i]))
    filename = paste0(dz_dir,'/',d_id,'/',f_id,'_delta.tif')
    if(!file.exists(filename)){
      rd = rast(files_sorted[i]) - rast(files_sorted[i-1])
      terra::writeRaster(rd, filename = filename)
    }
  })
})

#----not normalized delta

#create output dir
dn_string = paste0('nonorm-blocks','_delta')
dn_dir = paste0(bm_dir,'/',dn_string)
dir.check(dn_dir)

lapply(1:nrow(blocks_p), function(i){ #create subdirectories for every block
  d = paste0(dn_dir,'/',blocks_p$BLOCKNUM[i])
  dir.check(d)
})

indir = paste0(bm_dir,'/nonorm-blocks') #input files

in_subdirs = list.dirs(indir)
in_subdirs = in_subdirs[in_subdirs != indir]

#calculate difference rasters
lapply(in_subdirs, function(d){
  #load input block subdirectory
  d_id = basename(d)
  files = list.files(d, full.names = T, pattern = '\\.tif$')
  df = tibble(
    files
  ) %>%
    #get dates from filepaths
    mutate(year_month = str_extract(files, "\\d{4}_\\d{2}")) %>%
    mutate(year_month_day = paste0(year_month,'_01')) %>%
    mutate(date = as.Date(str_replace_all(year_month_day,'_','-'))) %>%
    #sort filepaths by increasing order of date
    arrange(date)
  
  files_sorted = df$files
  
  #calculate rasters from files in subdirectory
  pblapply(2:length(files_sorted), function(i){ #start at 
    
    f_id = file_path_sans_ext(basename(files[i]))
    filename = paste0(dn_dir,'/',d_id,'/',f_id,'_delta.tif')
    if(!file.exists(filename)){
      rd = rast(files_sorted[i]) - rast(files_sorted[i-1])
      terra::writeRaster(rd, filename = filename)
    }
  })
})

#----softmax delta

#create output dir
dsm_string = paste0(sm_dir_string,'_delta')
dsm_dir = paste0(bm_dir,'/',dsm_string)
dir.check(dsm_dir)

lapply(1:nrow(blocks_p), function(i){ #create subdirectories for every block
  d = paste0(dsm_dir,'/',blocks_p$BLOCKNUM[i])
  dir.check(d)
})

indir = sm_dir #input files

in_subdirs = list.dirs(indir)
in_subdirs = in_subdirs[in_subdirs != indir]

#calculate difference rasters
lapply(in_subdirs, function(d){
  #load input block subdirectory
  d_id = basename(d)
  files = list.files(d, full.names = T, pattern = '\\.tif$')
  df = tibble(
    files
  ) %>%
    #get dates from filepaths
    mutate(year_month = str_extract(files, "\\d{4}_\\d{2}")) %>%
    mutate(year_month_day = paste0(year_month,'_01')) %>%
    mutate(date = as.Date(str_replace_all(year_month_day,'_','-'))) %>%
    #sort filepaths by increasing order of date
    arrange(date)
  
  files_sorted = df$files
  
  #calculate rasters from files in subdirectory
  pblapply(2:length(files_sorted), function(i){ #start at 
    
    f_id = file_path_sans_ext(basename(files[i]))
    filename = paste0(dsm_dir,'/',d_id,'/',f_id,'_delta.tif')
    if(!file.exists(filename)){
      rd = rast(files_sorted[i]) - rast(files_sorted[i-1])
      terra::writeRaster(rd, filename = filename)
    }
  })
})
