source('planet_scenes_BEAST_preprocessing.R')

#define target feature to use for timeseries decomposition
target_feat = 'NDVI'

#set up directories for beast output
beast_dir = file.path('data','BEAST')
dir.check(beast_dir)

beast_data_dir = file.path(beast_dir,'ps_scenes_harmonized_coreg_2013-2024')
dir.check(beast_data_dir)

beast_feat_dir = file.path(beast_data_dir,target_feat)
dir.check(beast_feat_dir)

indir = scenes_dir #directory where all the input rasters are stored

#get files list
fl_in = list.files(indir, full.names = T, recursive=T)
fl_in = fl_in[!str_detect(fl_in,'_NULL')] #exclude dummy files

# #check layers in raster files
# lyr_l = pblapply(1:length(fl_in), function(i){
#   f = fl_in[i]
#   # cat('\n',
#   #     # f,' --- ',
#   #     i,'/',length(fl_in))
#   r = rast(f)
#   ns = names(r)
#   bool = target_feat %in% ns
#   # cat('\n',bool)
#   return(bool)
# })
# 
# lyr_l_vect = unlist(lyr_l)
# 
# tibble(missing = lyr_l_vect) |> group_by(missing) |> summarise(n = n())
# 
# missing_lyr_l = fl_in[!lyr_l_vect]
# 
# pblapply(missing_lyr_l, function(f){ #remove files with missing layers
#   if(file.exists(f)){file.remove(f)}
#   cat('\nRemoved',f)
# })

#develop BEAST model for each thinning block
library(Rbeast)
# library(raster)

datasets = basename(dirs)

lapply(datasets, function(d){dir.check(file.path(beast_feat_dir,d))}) #initialize dataset directories for beast models

dbgrid = expand.grid(datasets=datasets,block_ids=block_ids)

pblapply(1:nrow(dbgrid), function(i){
  
  b = dbgrid$block_ids[i]
  d = dbgrid$datasets[i]
  
  beast_f = paste0(beast_feat_dir,'/',d,'/',b,'.rds')
  
  cat('\nProcessing',beast_f)
  if(!file.exists(beast_f)){
    
    # #example data
    # zipfile = system.file("extdata/ndvi.zip", package="Rbeast") 
    # tmpfld  = file.path(tempdir(), "RbeastTmpDir" )                    # get a temp folder named "RbeastTmpDir"
    # dir.create(tmpfld)                                                 # create the tmp folder
    # unzip( zipfile, exdir= normalizePath(tmpfld ))                     # unzip the images to the tmp folder
    # flb=list.files(path=tmpfld, pattern='*.tif', full.names=TRUE) # the list of all  the tiff filenames
    # b_dates = basename(file_path_sans_ext(flb))
    
    # my data
    flb = fl_in[str_detect(fl_in,paste0(d,'/',b))]
    b_dates = meta_df$acquisition_date[meta_df$id %in% file_path_sans_ext(basename(flb))]
    
    cat('\nLoading rasters\n')
    # # r = lapply(flb,rast)|>rast()
    # r = raster::stack(flb, layers = target_feat)
    # names(r) = b_dates
    # inMemory(r)
    # r = raster::readAll(r)
    # r = raster::brick(r)
    # raster::inMemory(r)
    # # dim(rs) = c(ncol(r),nrow(r),nlyr(r))
    # # names(r) = b_dates
    
    # #get extents of all rasters, get maximum extent of raster
    # el = pblapply(flb, function(x){
    #   cat('\n',x)
    #   ext(rast(x,lyrs=target_feat))
    #   })
    
    r_geom_tbl = tibble(
      flb = flb,
      nrows = pbsapply(flb, function(x)nrow(rast(x,lyrs=target_feat))),
      ncols = pbsapply(flb, function(x)ncol(rast(x,lyrs=target_feat))),
      origin_x = pbsapply(flb, function(x)origin(rast(x,lyrs=target_feat))[1]),
      origin_y = pbsapply(flb, function(x)origin(rast(x,lyrs=target_feat))[2])
    )
    r_geom_summ = r_geom_tbl |>
      group_by(nrows, ncols, origin_x, origin_y) |>
      summarise(count = n())
    r_geom_tbl = r_geom_tbl |> left_join(r_geom_summ)
    r_geom_candidates = r_geom_tbl |>
      filter(count == max(count))
    template_r = rast(r_geom_candidates$flb[1],lyrs=target_feat)
    
    #load rasters, extend to max extent, stack
    rl = pblapply(flb, function(f){
      # cat('\nLoading',f)
      r = rast(f, lyrs=target_feat)
      # re = extend(r, max_ext, snap='out')
      rr = resample(r, template_r, method='cubic')
      return(rr)
    })
    
    r = rast(rl)
    names(r) = b_dates
    
    v = values(r)
    dim(v) = c(ncol(r),nrow(r),nlyr(r))
    
    #set up BEAST run
    metadata                  = list()        
    metadata$isRegularOrdered = FALSE    # IRREGULAR input
    metadata$whichDimIsTime   = 3        # 3rd dimension is the ts length
    metadata$time$datestr     = b_dates  # date info
    # metadata$time$strfmt      = 'LT05_018032_20080311.yyyy-mm-dd'
    metadata$deltaTime        = '1d'     # Aggregate Y into a daily ts
    metadata$period           = 1.0      # The period is 1 year: deltaTime*freq=1/12*12=1.0
    
    extra = list()
    extra$numThreadsPerCPU = 3
    extra$numParThreads = 30
    
    # tic()
    cat('\nRunning BEAST\n')
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
