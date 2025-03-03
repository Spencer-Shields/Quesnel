#helper functions

#----define function for plotting lists of rasters----

list_plot = function(l, single_plot = T){
  #function for plotting list of named single band rasters
  #rasters must be named
  #single_plot = T makes plots faceted in a grid, otherwise plots will appear in individual windows.
  if(single_plot == T){
    x11()
    par(mfrow = c(ceiling(sqrt(length(l))), ceiling(sqrt(length(l)))))
    
    for(i in 1:
        length(l)
    )
    {
      
      plot(l[[i]]
           , main = names(l)[i]
           , sub = as.character(substitute(l))
      )
    }
    par(mfrow = c(1,1))
    
  } else {
    
    lapply(X = 1:length(l), function(i){
      x11()
      plot(l[[i]], main = names(l)[i])
    })
  }
}

#----define function for loading separate landsat bands into a spatraster stack----

landstack = function(dir, what= 'bands'){
  #Stack individual Landsat bands into SpatRaster objects using metadata in the ML.txt file. Replacement for the RStoolbox function since that doesn't seem to be working. 
  #Input is a directory and output is a LIST of multiband spatrasters.
  #dir: the path to a directory containing landsat data and metadata.
  #what: the type of data to stack. Can be 'bands' or 'qa' for quality assurance layers.
  
  lsat_files = list.files(dir, full.names = T)
  meta_files = lsat_files[stringr::str_detect(lsat_files, 'MTL.txt')]
  
  raster_files = switch(what,
                        'bands' = lsat_files[stringr::str_detect(lsat_files, 'B\\d.TIF')],
                        'qa' = lsat_files[stringr::str_detect(lsat_files, '_QA_[A-Za-z]+\\.TIF')],
                        stop('Value provided to `what` not recognized.'))
  
  stack_list = pbapply::pblapply(X = 1:length(meta_files), FUN = function(i){
    
    #get product id
    lsatMeta = RStoolbox::readMeta(meta_files[i], raw=T)
    product_id = lsatMeta$PRODUCT_CONTENTS[rownames(lsatMeta$PRODUCT_CONTENTS)=='LANDSAT_PRODUCT_ID',]
    
    id_rast_files = raster_files[stringr::str_detect(raster_files, product_id)]
    rast_list = lapply(id_rast_files, rast)
    stacked_rast = terra::rast(rast_list)
    
    #apply band names
    band_names = sapply(1:length(rast_list), function(j){
      n = names(rast_list[[j]])
      sl = stringr::str_split_1(n, '_')
      sl[length(sl)]
    })
    names(stacked_rast) = band_names
    
    return(stacked_rast)
  })
  print('Your list of Landsat rasters is ready!')
  return(stack_list)
}

#----check if a directory exists, create it if it doesn't----
dir.check = function(dir){
  #check if a directory exists, create it if it doesn't
  if(!dir.exists(dir)){dir.create(dir)}
}

#----Calculate disturbance index using Tasseled Cap indices----
z_rast = function(r, maxcell = Inf){ 
  #function for calculating z-scores for a single-band raster
  g = global(r, c('mean','std'), maxcell = maxcell, na.rm=T)
  zr = (r-g[['mean']])/g[['std']]
} 

di_fn <- function(r, znormed = T) {
  #Calculate disturbance index based on Tasseled Cap indices.
  #Make sure that raster stack R has brightness, greenness, and wetness as the first, second, and third bands.
  #The znormed parameter gives the option of not calculating the z-score for each band prior to calculating the z-score.
  if(znormed == T){
  z_layers <- sapply(1:nlyr(r), function(i) z_rast(r[[i]]))  # Apply calculate z-scores for each layer
  } else{
    z_layers = r
  }
  d <- z_layers[[1]] - (z_layers[[2]] + z_layers[[3]])  # Compute disturbance index
  names(d) <- "DI" #rename band
  d
}

#----Extract PlanetScope metadata from a collection of json files----

ps_meta = function(dir, recursive = T){
  #Returns a tibble of PlanetScope meta data based on metadata.json files in a directory.
  files = list.files(dir, recursive = recursive, full.names = T)
  jsons = files[str_detect(files, 'metadata\\.json$')]
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
