#helper functions

#----define function for plotting lists of rasters----

list_plot = function(l, single_plot = T){
  #function for plotting list of named single-band rasters
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

list_plot_ggRGB = function(l, r = 3, g = 2, b = 1,...){
  x11()
  par(mfrow = c(ceiling(sqrt(length(l))), ceiling(sqrt(length(l)))))
  
  for(i in 1:length(l)){
    ggRGB(l[[i]])
  }
  par(mfrow = c(1,1))
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


#----Define softmax function----

softmax = function(raster, append_name = FALSE){
  #Calculate softmax function on a multiband spatraster.
  #Input and output == multiband spatraster.
  #append_names will add '_softmax' on the end of hte original band names.
  #NOTE: it is recommended to rescale data to low values (e.g. [0,1]) before using this function.
  sums = global(exp(raster), 'sum', na.rm = T)
  
  sm_l = pblapply(X = 1:nlyr(raster), FUN = function(i){
    num = exp(raster[[i]])
    denom = sums$sum[i]
    sm = num/denom
    return(sm)
  })
  
  sm = rast(sm_l)
  
  if(append_name == T){
    names(sm) = paste0(names(raster),'_softmax')
  }
  
  return(sm)
  
}

#----Calculate Hue Index----

hue.vi = function(r, red = 3, green = 2, blue = 1, type = 'Escadafal'){
  #calculate the Hue Index using RGB bands
  #bands can be specified using numerical indexing or character vectors representing band names
  #the 'type' argument specifies whether the equation given by Escadafal et al. (1994) or Mandal (2016)
  
  b_r = r[[red]]
  b_b = r[[blue]]
  b_g = r[[green]]
  
  if(type == 'Escadafal'){
    
    h = atan((b_g-b_b)*(2*b_r-b_g-b_b)/30.5)
    
  } else {
    
    if(type == 'Mandal'){
      
      h = (2*b_r-b_g-b_b)/(b_g-b_b)
      
    } else {
      print('Error: type of Hue unrecognized')
    }
  }
  
  names(h) = 'Hue'
  return(h)
}

#----Calculate Normalized Difference Green-Red index----

ndgr.vi = function(r, red = 3, green = 2){
 # Calculate Normalized Difference Green-Red index
  b_r = r[[red]]
  b_g = r[[green]]
  
  ndgr = (b_g-b_r)/(b_g+b_r)
  names(ndgr) = 'NDGR'
  return(ndgr)
}

#----Calculate Green Chromatic Coordinate----

gcc.vi = function(r, red = 3, green = 2, blue = 1){
  #Calculate Green Chromatic Coordinate
  b_r = r[[red]]
  b_g = r[[green]]
  b_b = r[[blue]]
  
  gcc = b_g/(b_g+b_r+b_b)
  names(gcc) = 'GCC'
  return(gcc)
}

#----Calculate Carotenoid Reflectance Index 550 ----

cri550.vi = function(r, red = 3, green = 2){
  #Calculate Carotenoid Reflectance Index 550
  b_r = r[[red]] 
  b_g = r[[green]]
  
  b_r = b_r+1 #add one so that results are not undefined
  b_g = b_g+1 #add one so that results are not undefined
  
  vi = (1/b_r)-(1/b_g)
  names(vi) = 'CRI550'
  return(vi)
}

#----Calculate Green Leaf Index----

gli.vi = function(r, red = 3, green = 2, blue = 1){
  #Calculate Green Leaf Index
  b_r = r[[red]]
  b_g = r[[green]]
  b_b = r[[blue]]
  
  vi = (2*b_g-b_r-b_b)/(2*b_g+b_r+b_b)
  names(vi) = 'GLI'
  return(vi)
}

#----Calculate Coloration Index----

ci.vi = function(r, red = 3, green = 2){
  #Calculate Coloration Index
  b_r = r[[red]]
  b_g = r[[green]]
  
  vi = (b_r-b_g)/(b_r+b_g)
  names(vi) = 'CI'
  return(vi)
}

#----Calculate Brightness Index----

bi.vi = function(r, red = 3, green = 2){
  #Calculate Brightness Index
  b_r = r[[red]]
  b_g = r[[green]]
  
  vi = sqrt((b_r^2+b_g^2)/2)
  names(vi) = 'BI'
  return(vi)
}

#----Calculate Transformed Vegetation Index----

tvi.vi = function(r, red = 3, green = 2){
  #Calculate Transformed Vegetation Index
  b_r = r[[red]]
  b_g = r[[green]]
  
  vi = sqrt(0.5+((b_r-b_g)/(b_r+b_g)))
  names(vi) = 'TVI'
  return(vi)
}

#----Calculate Visible Atmospherically Resistant Index Green----

varig.vi = function(r, red = 3, green = 2, blue = 1){
  #Calculate Visible Atmospherically Resistant Index Green
  b_r = r[[red]]
  b_g = r[[green]]
  b_b = r[[blue]]
  
  vi = (b_r-b_b)/(b_g+b_r-b_b)
  names(vi) = 'VARIgreen'
  return(vi)
}

#----Calculate multiple vegetation indices for a spatraster by providing a list of functions ----
rast.batch.functions = function(r, fl, rast_out = T, include_input = T, ...){
  #apply a list of functions - including custom functions - which take the same arguments to a spatraster
  #the purpose is to calculate multiple vegetation indices with one compact call using custom functions
  #if rast_out ==T, result will be a spatraster, else result will be a list of spatrasters
  result = lapply(fl, function(f)f(r,...))
  
  if(rast_out==T){ #option to make result a spatraster
    result = rast(result)
  }
  
  if(include_input == T){
    result = c(r, result)
  }
  
  return(result)
}
