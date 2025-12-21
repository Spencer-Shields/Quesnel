library(tidyverse)
library(terra)
library(Rcpp)
library(tools)
library(future)
library(future.apply)
library(parallel)
library(RStoolbox)
source('helper_functions.R')
source('https://raw.githubusercontent.com/Spencer-Shields/planetscope_radiometric_correction/refs/heads/main/check_radiometric_consistency.R')

#set directories
ps_dir = 'data/planet_scenes'
# raw_dir = file.path(ps_dir,'raw_additional')
raw_dir = file.path(ps_dir, "Quesnel_Clip_Harm_Coreg_Product=analytic_sr_udm2_Time=2013-01-01to2024-12-31_MaxCloudProp=0.99")

#load thinning block vector data
blocks = terra::vect('data/Quesnel_thinning/12l_12n_bdy.geojson')
block_ids = blocks$BLOCKNUM
names(block_ids) = block_ids

#tweak ps_meta function to handle cases where coregistration metadata is also present
ps_meta = function(dir, coreg_qa=F, recursive = T){
  
  #Returns a tibble of PlanetScope meta data based on metadata.json files in a directory.
  files = list.files(dir, recursive = recursive, full.names = T)
  jsons = files[str_detect(files, 'metadata\\.json$')]
  
  #exclude coregistration quality assurance metadata
    jsons = jsons[!str_detect(jsons, 'coreg_qa')]
  
  if(length(jsons)==0){
    print('No metadata files found in directory.')
  } else {
    jsons_l = lapply(jsons, jsonify::from_json)
    
    json_df <- map_dfr(jsons_l, function(x) {
      tibble(id = x$id, !!!x$properties)
    })
    return(json_df)
  }
}

#get dataframe of scene metadata, get dates
meta_df = ps_meta(raw_dir) |> distinct()
meta_df[['acquisition_date']] = as.Date(sapply(meta_df$acquired, function(x){
  s = str_split_1(x, 'T')[1]
  return(s)
}))
meta_df = meta_df |>
  mutate(acquisition_year = year(acquisition_date),
         acquisition_month = month(acquisition_date),
         acquisition_day = day(acquisition_date)) %>%
  mutate(acquisition_month_day = paste0(acquisition_month,'-',acquisition_day))

# coreg_qa_fl = list.files(raw_dir, pattern='qa_metadata.json$', full.names = T, recursive = T) 
# coreg_qa_df = pblapply(coreg_qa_fl, function(f){
#   j = jsonify::from_json(f)
#   jdf <- data.frame(
#     coreg_method               = j$coreg_method,
#     coreg_method_version       = j$coreg_method_version,
#     anchor_item                = j$anchor_item,
#     target_item                = j$target_item,
#     projection_epsg_before     = j$projection_epsg_before,
#     projection_epsg_after      = j$projection_epsg_after,
#     pixel_shift_x              = j$pixel_shift[1],
#     pixel_shift_y              = j$pixel_shift[2],
#     matching_score_before      = j$matching_score_before,
#     matching_score_after       = j$matching_score_after,
#     matching_score_improvement = j$matching_score_improvement,
#     details                    = j$details,
#     output_item_scene          = j$output_item[1],
#     output_item_udm            = j$output_item[2]
#   )
# }) |> bind_rows()

#get files list
all_files = list.files(raw_dir, full.names = T, recursive = T)
tif_fl = list.files(raw_dir, recursive = T, full.names = T, pattern = '\\.tif$')
udm_fl = tif_fl[str_detect(basename(tif_fl), 'udm')]
scene_fl = tif_fl[str_detect(basename(tif_fl), '_SR_')]

#filter meta_df to only contain items that have both a scene and a UDM
# meta_df_ = meta_df |>
#   filter(str_detect(scene_fl,id) & str_detect(udm_fl,id))

meta_df <- meta_df |>
  filter(
    sapply(id, function(x) any(str_detect(scene_fl, fixed(x)))) &
      sapply(id, function(x) any(str_detect(udm_fl, fixed(x))))
  )

#reproject blocks
blocks = project(blocks, crs(rast(scene_fl[1])))
blocks_wrapped = wrap(blocks)

#check resolutions and projections (result: all equal, resolution 3m crs wgs84)

# projections_l = list()
# resolutions_l = list()
# 
# for(i in 1:length(scene_fl)){
#   cat('\r',i,'/',length(scene_fl))
#   f = scene_fl[i]
#   r = rast(f)
#   
#   reso = res(r)[1]
#   resolutions_l[i] = reso
#   
#   proj = crs(r, describe=T)$name
#   projections_l[i] = proj
# }

#make output directory
scenes_dir = file.path(ps_dir, 'additional_scenes')
dir.check(scenes_dir)

#get ids of raw rasters and list of raw raster files
ids = meta_df$id

#set scale factor
scale_factor = 10000 #what the reflectance values get divided by to facilitate calculating indices and softmax values

#define growing season months (May through October)
gs_months = c(5:10)

#create directories for each type of dataset to store the processed rasters
nonnorm_string = 'Non-normalized'
nonnorm_dir = paste0(scenes_dir,'/', nonnorm_string,'_20250910')
dir.check(nonnorm_dir)

z_string = 'Z'
z_dir = paste0(scenes_dir,'/', z_string,'_20250910')
dir.check(z_dir)

zr_string = 'Zrobust'
zr_dir = paste0(scenes_dir, '/', zr_string,'_20250910')
dir.check(zr_dir)

sm_string = 'SM'
sm_dir = paste0(scenes_dir,'/',sm_string,'_20250910')
dir.check(sm_dir)

dirs = c(z_dir, zr_dir, sm_dir, nonnorm_dir)

#make subdirectories within dataset directories to store scenes for each thinning block
pblapply(dirs, function(x){
  pblapply(block_ids, function(y){
    d = paste0(x,'/',y)
    dir.check(d)
  })
})

#check if there are fewer files in the output folders than expected, engage in preprocessing if yes
if(length(list.files(dirs,recursive = T,pattern='\\.tif$')) < length(dirs)*length(block_ids)*length(ids)){
  
  #process rasters
  ps_preprocessor_wrapper = function(){
    
    # #define dummy raster to save if the actual raster ends up being blank or invalid (saves time reloading data later)
    # r_dummy = rast(nrow=1,ncol=1)
    # values(r_dummy) = NA #r_dummy has a single cell, the value is NA
    # r_dummy_wrapped = wrap(r_dummy)
    
    #set up parallel processing
    clust = makeCluster(8)
    plan('cluster', workers = clust)
    on.exit(stopCluster(clust))
    on.exit(plan('sequential'), add=T)
    
    start_idx = 1
    
    future_lapply(1:length(ids), function(i){
    # pblapply(start_idx:length(ids), function(i){
      
      #get file id
      id = ids[i]
      
      #make dummy raster to save
      r_dummy = rast(nrow=1,ncol=1)
      values(r_dummy) = NA #r_dummy has a single cell, the value is NA
      
      lapply(1:length(block_ids), function(j){
        
        block_dir = block_ids[j]
        
        #Non-normalized raster
        nonnorm_file = paste0(nonnorm_dir,'/',block_dir,'/',id,'.tif')
        nonnorm_dummy_file = paste0(file_path_sans_ext(nonnorm_file),'_NULL.tif')
        cat('\nProcessing',i,'/',length(ids),'---',j,'/',length(block_ids),'; ', nonnorm_file, '\n')
        if(!file.exists(nonnorm_file) & !file.exists(nonnorm_dummy_file)){
          
          # #unwrap r_dummy for workers
          # r_dummy = unwrap(r_dummy_wrapped)
          
          #load raster
          f_r = scene_fl[str_detect(scene_fl,id)]
          if(length(f_r)>1){f_r = f_r[which.max(str_length(f_r))]}
          r = rast(f_r)
          # r = r[[feats_to_use]] #subset raster to use features described above
          
          #mask using UDM
          f_udm = udm_fl[str_detect(udm_fl,id)]
          udm = rast(f_udm)
          
          r = mask(r, udm[['cloud']]) #remove cloud
          r = mask(r, udm[['shadow']]) #remove shadow
          r = mask(r, udm[['haze_heavy']]) #remove heavy haze
          
          #load thinning block boundary
          blocks = unwrap(blocks_wrapped)
          block = blocks[blocks$BLOCKNUM==block_dir,]
          
          
          #check if thinning block and raster overlap do not overlap OR if overlap only contains NAs, save dummy file, proceed otherwise
          if(
            !any(relate(r,block,'intersects')) | 
            all(is.na(extract(r,block)|>select(-ID)))
          ){
            
            terra::writeRaster(r_dummy, nonnorm_dummy_file)
            
          } else {
            
            #crop to block
            r = crop(r, block
                     , mask = T
            )
            
            #check if raster only contains NA values, save dummy file if yes, proceed with preprocessing otherwise
            g = global(r, 'isNA', na.rm=T)
            
            if(g[1,1]==ncell(r)){
              
              terra::writeRaster(r_dummy, nonnorm_dummy_file)
              
            } else {
              
              #calculate vegetation indices, combine all features into single spatraster
              {
                #rescale raster
                r = r/scale_factor
                
                #calculate vegetation indices
                si = spectralIndices(img = r, blue = 'blue', green = 'green', red = 'red', nir = 'nir'
                                     # , skipRefCheck = T
                )
                #remove EVI2 (since spectralIndices doesn't calculate it well for some reason)
                if('EVI2' %in% names(si)){si = subset(si, 'EVI2', negate=T)}
                
                #calculate additional indices with custom functions
                vi = rast.batch.functions(r, include_input = F,
                                          #band mapping
                                          blue=1, green=2, red=3, nir=4, 
                                          #list of indices to calculate
                                          fl = c(
                                            hue.vi,
                                            gcc.vi,
                                            ndgr.vi,
                                            bi.vi,
                                            ci.vi,
                                            cri550.vi,
                                            gli.vi,
                                            # tvi.vi, #no TVI since spectralIndices handles this
                                            varig.vi,
                                            evi2.vi
                                          ))
                
                r = c(r, si, vi)
              }
              
              #save file
              terra::writeRaster(r, nonnorm_file, datatype='FLT8S')
            }
          }
          
        } else { #if the file doesn't exist make it, if it does exist load it for further processing checks
          
          #load the dummy file if nonnorm file does not exist
          if(file.exists(nonnorm_file)){
            r = rast(nonnorm_file)
          } else {
            r = rast(nonnorm_dummy_file)
          }
        }
        
        #Z-score raster
        z_file = paste0(z_dir,'/',block_dir,'/',id,'.tif')
        z_dummy_file = paste0(file_path_sans_ext(z_file),'_NULL.tif')
        
        # print(paste('Processing',z_file))
        if(!file.exists(z_file)&!file.exists(z_dummy_file)){
          
          if(ncell(r) == ncell(r_dummy)){
            writeRaster(r, z_dummy_file)
          } else {
            z_r = z_rast(r)
            writeRaster(z_r, z_file)
          }
        }
        
        #robust Z-score raster
        zr_file = paste0(zr_dir,'/',block_dir,'/',id,'.tif')
        zr_dummy_file = paste0(file_path_sans_ext(zr_file),'_NULL.tif')
        
        # print(paste('Processing',zr_file))
        if(!file.exists(zr_file)&!file.exists(zr_dummy_file)){
          
          if(ncell(r) == ncell(r_dummy)){
            writeRaster(r, zr_dummy_file)
          } else {
            zr_r = z_rast(r,robust = T)
            writeRaster(zr_r, zr_file)
          }
        }
        
        #softmax raster
        sm_file = paste0(sm_dir,'/',block_dir,'/',id,'.tif')
        sm_dummy_file = paste0(file_path_sans_ext(sm_file),'_NULL.tif')
        
        # print(paste('Processing',sm_file))
        if(!file.exists(sm_file)&!file.exists(sm_dummy_file)){
          
          if(ncell(r) == ncell(r_dummy)){
            writeRaster(r, sm_dummy_file)
          } else{
            
            sm_r = softmax(r, append_name = F)
            # sm_r = sm_r * scale_up_factor
            writeRaster(sm_r, sm_file, datatype='FLT8S')
          }
        }
        
      })
      
    }
    # ,cl = 'future'
    # ,future.seed=T
    )
    
    stopCluster(clust)
    plan('sequential')
  }
  
  ps_preprocessor_wrapper()
}




