source('scene_setup_preprocessing_20250909.R')

hls_dir = 'data/HLS'
dir.check(hls_dir)
landsat_dir = 'data/landsat'
dir.check(landsat_dir)
hlss30_raw_dir = file.path(hls_dir,'raw_hlss30')
hlsl30_raw_dir = file.path(hls_dir,'raw_hlsl30')
landsat_raw_dir = file.path(landsat_dir,'raw')

#----Pre-process HLS data----

#get list of all HLS files
hls_fl_raw = c(
  list.files(hlsl30_raw_dir, full.names=T),
  list.files(hlss30_raw_dir, full.names = T)
)

#get list of HLS scene IDs
hls_raw_ids = unique( str_replace(
  basename(hls_fl_raw), "v2\\.0.*$", "v2.0"))

#subset ids to only include growing season
hls_meta = tibble(id = hls_raw_ids) |>
  mutate(acquired_dt = as.POSIXct(substr(id,16,29), format = "%Y%jT%H%M%S", tz = 'UTC')) |>
  mutate(acquisition_month = month(as.Date(acquired_dt))) |>
  mutate(acquisition_day = day(as.Date(acquired_dt))) |>
  mutate(acquired = as.Date(acquired_dt))
# |>
#   filter(acquisition_month >= min(gs_months) & acquisition_month < (max(gs_months)+1))
  

#reproject blocks to match HLS
blocks_hlscrs = project(blocks_p,
                        rast(hls_fl_raw[1]))

#make NBR directory
hls_nbr_dir = file.path(hls_dir,'NBR')
dir.check(hls_nbr_dir)

#loop to preprocess HLS files
pblapply(1:length(hls_raw_ids), function(i){
  
  cat('\r',i,'/',length(hls_raw_ids))
  
  #get scene id to process
  id = hls_raw_ids[i]
  
  #stack, mask
  hls_masked_dir = file.path(hls_dir,'masked')
  dir.check(hls_masked_dir)
  
  masked_f = paste0(hls_masked_dir,'/',id,'.tif')
  if(!file.exists(masked_f)){
    #get list of files which correspond with the target id
    idl = hls_fl_raw[str_detect(hls_fl_raw,id)]
    
    #stack surface reflectance bands
    idl = hls_fl_raw[str_detect(hls_fl_raw,id)]
    r = lapply(idl[str_detect(idl,'2.0.B')], rast) |> rast()
    
    #load fmask
    fmask = idl[str_detect(idl,'Fmask')] |> rast()
    
    #convert fmask from integer into multiband QA raster
    {
      bits = 0:5
      fv = values(fmask)
      qa = lapply(bits, function(b){
        m = as.integer(bitwAnd(fv, bitwShiftL(1, b)) != 0)
        r = fmask
        values(r) = m
        return(r)
      }) |> rast()
      names(qa) <- c(
        "cirrus",
        "cloud",
        "adjacent_to_cloud_shadow",
        "cloud_shadow",
        "snow_ice",
        "water"
      )
    }
    
    #mask HLS using all QA layers
    {
      qa_sum = app(qa, sum)
      qa_m = ifel(qa_sum > 0,NA,1)
      r = mask(r, qa_m)
    }
    
    writeRaster(r, masked_f, overwrite=T)
  }
  
  
  #calculate NBR
  {
    nbr_notcropped_dir = file.path(hls_dir,'NBR_notcropped')
    dir.check(nbr_notcropped_dir)
    
    nbr_nc_f = paste0(nbr_notcropped_dir,'/',id,'.tif')
    if(!file.exists(nbr_nc_f)){
      r = rast(masked_f)
      nir = ifelse( #get appropriate name for the NIR band
        any(str_detect(names(r), 'NIR_Broad'))
        , 'NIR_Broad'
        , 'NIR')
      nbr = (r[[nir]] - r[['SWIR1']])/(r[['SWIR1']] + r[[nir]])
      names(nbr) = 'NBR'
      writeRaster(nbr,nbr_nc_f)
    }
  }
  
  #reproject and crop to thinning blocks
  {
    nbr_fl = paste0(hls_nbr_dir,'/',block_ids,'/',id,'.tif')
    if(any(!file.exists(nbr_fl))){
      
      #load data
      r = rast(nbr_nc_f)
      
      #reproject to match PlanetScope
      r = project(r, target_crs)
      
      #crop to thinning blocks
      lapply(block_ids, function(b){
        blockdir = file.path(hls_nbr_dir,b)
        dir.check(blockdir)
        f = paste0(blockdir,'/',id,'.tif')
        if(!file.exists(f)){
          block = blocks_p[blocks_p$BLOCKNUM==b]
          rc = crop(r,block
                    , mask=T)
          writeRaster(rc,f)
        }
      })
    }
  }
  
})

# ggplot()+
#   geom_spatraster(data=w[['NIR_Broad']])+
#   geom_spatvector(data=blocks_hlscrs, aes(color=BLOCKNUM))
# 
# ggplot()+
#   geom_spatraster(data = rc)+
#   geom_spatvector(data = block, fill=NA)

#----Pre-process Landsat data----

#----BEAST timeseries modelling----
beast_dir = file.path('data','BEAST')
dir.check(beast_dir)

beast_data_dir = file.path(beast_dir,'HLS_2013-2024')
dir.check(beast_data_dir)

beast_feat_dir = file.path(beast_data_dir,'NBR')
dir.check(beast_feat_dir)

indir = hls_nbr_dir #directory where all the input rasters are stored

fl_in = list.files(indir, full.names = T, recursive=T)

#develop BEAST model for each thinning block
library(Rbeast)
library(raster)
pblapply(block_ids, function(b){
  
  beast_f = paste0(beast_feat_dir,'/',b,'.rds')
  if(!file.exists(beast_f)){
    
    # #example data
    # zipfile = system.file("extdata/ndvi.zip", package="Rbeast") 
    # tmpfld  = file.path(tempdir(), "RbeastTmpDir" )                    # get a temp folder named "RbeastTmpDir"
    # dir.create(tmpfld)                                                 # create the tmp folder
    # unzip( zipfile, exdir= normalizePath(tmpfld ))                     # unzip the images to the tmp folder
    # flb=list.files(path=tmpfld, pattern='*.tif', full.names=TRUE) # the list of all  the tiff filenames
    # b_dates = basename(file_path_sans_ext(flb))
    
    # my data
    flb = fl_in[str_detect(fl_in,b)]
    b_dates = hls_meta$acquired[hls_meta$id %in% file_path_sans_ext(basename(flb))]
    
    # r = lapply(flb,rast)|>rast()
    r = raster::stack(flb)
    names(r) = b_dates
    inMemory(r)
    r = raster::readAll(r)
    r = raster::brick(r)
    raster::inMemory(r)
    # dim(rs) = c(ncol(r),nrow(r),nlyr(r))
    # names(r) = b_dates
    
    v = raster::values(r)
    dim(v) = c(ncol(r),nrow(r),nlayers(r))
    
    #set up BEAST run
    metadata                  = list()        
    metadata$isRegularOrdered = FALSE    # IRREGULAR input
    metadata$whichDimIsTime   = 3        # 437 is the ts length, so set it to '3' here.
    metadata$time$datestr     = b_dates  # date info is contained in the file names
    # metadata$time$strfmt      = 'LT05_018032_20080311.yyyy-mm-dd'
    metadata$deltaTime        = '1d'     # Aggregate Y into a monthly ts: 1/12 year
    metadata$period           = 1.0      # The period is 1 year: deltaTime*freq=1/12*12=1.0
    
    extra = list()
    # extra$numThreadsPerCPU = 3
    # extra$numParThreads = 30
    extra$numThreadsPerCPU = 1
    extra$numParThreads = 1
    
    # tic()
    BEAST = beast123(Y = v
                     , metadata = metadata
                     ,extra = extra
                     )
    # toc()
    out = list(BEAST, r)
    names(out) = c('BEAST_object', 'input_raster')
    write_rds(out, beast_f)
    
  }
})
