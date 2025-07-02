#generate basemap GLCM textures

source('basemap_setup_preprocessing.R')
library(GLCMTextures)

#setup parent directory for storing textures
glcm_dir = file.path(bm_dir,'basemap_GLCM_textures')
dir.check(glcm_dir)


#function for calculating glcm_textures for a multiband raster
glcm_textures_multiband = function(r, ...){
  #calculate glcm textures for a raster, combine output into a single multiband raster
  bands = names(r)
  glcm_l = pblapply(bands, function(b){
    r_sub = r[[b]]
    g = glcm_textures(r=r_sub, ...)
    names(g) = paste0(names(g),'_',b)
    return(g)
  })
  glcm_combined = rast(glcm_l)
  return(glcm_combined)
}
tic()
n = glcm_textures_multiband(r[[1:4]], n_levels = 32, quant_method='prob')
toc()

#function for calculating glcm metrics for every block file in a direcotry
make_glcm_textures_fn = function(indir, outdir, w = c(3,3), layers=NULL, ...){
  #Calculate glcm textures for each basemap in each block. 
  #indir: the input directory 
  #layers: vector of layer names for which textures will be calculated (NULL will use all)
  #... other arguments to be passed to glcm_textures
  
  #check that output directory exists, create if it doesn't
  dir.check(outdir)
  
  #make directory in outdir based on window dimensions
  w_dir = file.path(outdir, paste0(w[1],'x',w[2]))
  dir.check(w_dir)
  
  #get list of input directories that contain basemaps for each block
  in_block_dirs = list.dirs(indir, recursive=F)
  
  #iterate over input block directories
  pblapply(1:length(in_block_dirs), function(i){
    
    #get list of input rasters
    in_block_dir = in_block_dirs[i]
    in_rasts = list.files(in_block_dir, full.names = T, pattern = '\\.tif$')
    
    #make output directory
    block_id = basename(in_block_dir)
    block_dir = file.path(w_dir, block_id)
    dir.check(block_dir)
    
    #iterate over input rasters
    pblapply(in_rasts, function(f){
      f_id = basename(f)
      filename = file.path(block_dir, f_id)
      if(!file.exists(filename)){
        
        #load raster
        r = rast(f)
        
        #subset raster using specified layers
        if(!is.null(layers)){
          r = r[[layers]]
        }
        
        #get multiband glcm metrics
        
        glcm = glcm_textures_multiband(r, w, ...)
        
        writeRaster(glcm, filename)
      }
    })
  })
}

dirs_tbl = tibble(
  indirs = c(nonorm_dir, z_dir, zr_dir),
  outdirs = file.path(glcm_dir, c('Non-normalized', 'Z', 'Z_robust'))
)

# plan('multisession', workers=4)
pblapply(1:nrow(dirs_tbl), function(i){
  make_glcm_textures_fn(indir = dirs_tbl$indirs[i],
  outdir = dirs_tbl$outdirs[i],
  w = c(3,3),
  quant_method = 'prob',
  n_levels = 32,
  na.rm=T,
  metrics = c('glcm_mean', 'glcm_contrast', 'glcm_entropy', 'glcm_variance'),
  layers = c('blue', 'green', 'red', 'Hue', 'GCC', 'NDGR', 'BI', 'CI'))
}
# ,cl='future'
)
# plan('sequential')
