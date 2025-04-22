# #helper functions
# library(Rcpp)
# library(RcppArmadillo)

#----define function for plotting lists of rasters----

list_plot = function(l, single_plot = T, new_window=T, col = viridis::viridis(101)){
  #function for plotting list of named single-band rasters
  #rasters must be named
  #single_plot = T makes plots faceted in a grid, otherwise plots will appear in individual windows.
  if(single_plot == T){
    if(new_window){x11()}
    par(mfrow = c(ceiling(sqrt(length(l))), ceiling(sqrt(length(l)))))
    
    for(i in 1:
        length(l)
    )
    {
      
      plot(l[[i]]
           , main = names(l)[i]
           , sub = as.character(substitute(l))
           , col = col
      )
    }
    par(mfrow = c(1,1))
    
  } else {
    
    lapply(X = 1:length(l), function(i){
      if(new_window){x11()}
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

#----Calculate mean absolute deviation----
mad = function(x, na.rm = FALSE) {
  if (na.rm) {
    x = x[!is.na(x)]  # Remove NA values
  }
  mean(abs(x - mean(x)))
}
#----Calculate z-score for a single band raster----
z_rast = function(r, maxcell = Inf, robust = F){ 
  #function for calculating z-scores for a single-band raster
  #robust = T will do the calculation using the median and mean absolute deviation instead of mean and std
  if(robust == F){
    g = global(r, c('mean','std'), maxcell = maxcell, na.rm=T)
    z = (r-g[['mean']])/g[['std']]
  } else {
    g1 = global(r, median, maxcell = maxcell, na.rm=T)
    g2 = global(r, mad, maxcell = maxcell, na.rm=T)
    g = cbind(g1,g2)
    names(g) = c('median', 'mad')
    z = (r - g[['median']])/g[['mad']]
  }
  return(z)
} 


#----Calculate disturbance index using Tasseled Cap indices----
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

#----function for choosing which substring in a vector appears in a filepath or other long string----
find_substring = function(filepath, v, first = T){
  #given a filepath (or other long string) and a character vector, return the element of the vector which appears in the path
  #first = T will return the first result if a vector is returned
  g = sapply(v, function(x)str_detect(filepath, x))
  # if(!T %in% g){
  #   return(NA)
  # } else {
  #   i = which(g)
  #   a = v[i]
  #   return(a)
  # }
  
  if (!any(g)) {
    return(NA)
  } else {
    a <- v[which(g)]
    if (first) {
      return(a[1])  # Return only the first match
    } else {
      return(a)     # Return all matches
    }
  }
}


#----Implement Otsu's method for SpatRasters----
# C++ implementation of Otsu's algorithm
cppFunction('
double otsuThresholdCpp(NumericVector values, int bins) {
    // Remove NA values
    LogicalVector na_mask = is_na(values);
    NumericVector validValues = values[!na_mask];
    int n = validValues.size();
    
    if (n == 0) {
        return NA_REAL;
    }
    
    // Find min and max values
    double minVal = min(validValues);
    double maxVal = max(validValues);
    double range = maxVal - minVal;
    
    if (range <= 0) {
        return minVal; // No variation in the data
    }
    
    // Create histogram
    IntegerVector histogram(bins, 0);
    for (int i = 0; i < n; i++) {
        int bin = std::min(bins - 1, (int)((validValues[i] - minVal) / range * bins));
        histogram[bin]++;
    }
    
    // Calculate cumulative sums
    std::vector<double> p(bins);
    std::vector<double> cumsum_p(bins);
    std::vector<double> cumsum_mean(bins);
    
    double total = n;
    double sum_intensity = 0;
    
    // Calculate probabilities and intensity sums
    for (int i = 0; i < bins; i++) {
        double intensity = minVal + (i + 0.5) * range / bins; // bin center
        p[i] = histogram[i] / total;
        sum_intensity += p[i] * intensity;
    }
    
    // Compute cumulative values
    cumsum_p[0] = p[0];
    cumsum_mean[0] = p[0] * (minVal + 0.5 * range / bins);
    
    for (int i = 1; i < bins; i++) {
        double intensity = minVal + (i + 0.5) * range / bins;
        cumsum_p[i] = cumsum_p[i-1] + p[i];
        cumsum_mean[i] = cumsum_mean[i-1] + p[i] * intensity;
    }
    
    // Find optimal threshold
    double max_variance = -INFINITY;
    int optimal_k = 0;
    
    for (int k = 0; k < bins-1; k++) {
        double w0 = cumsum_p[k];
        if (w0 <= 0 || w0 >= 1) continue;
        
        double w1 = 1.0 - w0;
        double mu0 = cumsum_mean[k] / w0;
        double mu1 = (sum_intensity - cumsum_mean[k]) / w1;
        
        double variance = w0 * w1 * (mu0 - mu1) * (mu0 - mu1);
        
        if (variance > max_variance) {
            max_variance = variance;
            optimal_k = k;
        }
    }
    
    // Calculate threshold from bin index
    return minVal + (optimal_k + 1) * range / bins;
}

NumericVector otsuMultipleCpp(NumericMatrix values, int bins) {
    int ncols = values.ncol();
    NumericVector thresholds(ncols);
    
    for (int col = 0; col < ncols; col++) {
        NumericVector col_values = values(_, col);
        thresholds[col] = otsuThresholdCpp(col_values, bins);
    }
    
    return thresholds;
}
')

#' Apply fast Otsu's thresholding algorithm to a SpatRaster
#'
#' @param rast A SpatRaster object (single layer/band)
#' @param bins Number of bins for the histogram (default: 256)
#' @param return_binary Logical; if TRUE, returns a binary SpatRaster (default: FALSE)
#' @param sample_size Number of pixels to sample for large rasters (default: NULL - use all pixels)
#' @return If return_binary=TRUE, returns a binary SpatRaster. Otherwise, returns the threshold value.
#' @export
otsu_spatraster_fast <- function(rast, bins=256, return_binary=FALSE, sample_size=NULL) {
  # Check if input is a SpatRaster
  if (!inherits(rast, "SpatRaster")) {
    stop("Input must be a SpatRaster object")
  }
  
  # Check if raster has multiple layers
  if (nlyr(rast) > 1) {
    warning("Input has multiple layers. Using only the first layer.")
    rast <- rast[[1]]
  }
  
  # Sample values if needed (for large rasters)
  if (!is.null(sample_size) && sample_size < ncell(rast)) {
    set.seed(42)  # For reproducibility
    cell_indices <- sample(ncell(rast), sample_size)
    values <- terra::extract(rast, cell_indices)[,1]
  } else {
    values <- terra::values(rast, mat=FALSE)
  }
  
  # Calculate threshold using C++ function
  threshold <- otsuThresholdCpp(values, bins)
  
  # Return binary raster or threshold value
  if (return_binary) {
    return(rast > threshold)
  } else {
    return(threshold)
  }
}

#' Apply fast Otsu's thresholding algorithm to a multi-band SpatRaster
#'
#' @param rast A SpatRaster object (can be multi-band)
#' @param method Method to combine bands: "all" (logical AND), "any" (logical OR), 
#'               "majority" (majority vote), or "individual" (return list of thresholds)
#' @param bins Number of bins for the histogram (default: 256)
#' @param sample_size Number of pixels to sample for large rasters (default: NULL - use all pixels)
#' @return If method is "all", "any", or "majority", returns a binary SpatRaster.
#'         If method is "individual", returns a list of threshold values.
#' @export
otsu_multiband_fast <- function(rast, method="majority", bins=256, sample_size=NULL) {
  # Check if input is a SpatRaster
  if (!inherits(rast, "SpatRaster")) {
    stop("Input must be a SpatRaster object")
  }
  
  # Get number of layers
  n_layers <- nlyr(rast)
  
  # If single layer, use standard function
  if (n_layers == 1) {
    return(otsu_spatraster_fast(rast, bins=bins, return_binary=method != "individual", 
                                sample_size=sample_size))
  }
  
  # Sample values if needed (for large rasters)
  if (!is.null(sample_size) && sample_size < ncell(rast)) {
    set.seed(42)  # For reproducibility
    cell_indices <- sample(ncell(rast), sample_size)
    values <- terra::extract(rast, cell_indices)
  } else {
    values <- terra::values(rast)
  }
  
  # Calculate thresholds using C++ function
  thresholds <- otsuMultipleCpp(values, bins)
  
  # Return result based on method
  if (method == "individual") {
    return(thresholds)
  } else {
    # Create binary rasters
    binary_rasters <- vector("list", n_layers)
    for (i in 1:n_layers) {
      binary_rasters[[i]] <- rast[[i]] > thresholds[i]
    }
    
    if (method == "all") {
      result <- binary_rasters[[1]]
      for (i in 2:n_layers) {
        result <- result & binary_rasters[[i]]
      }
      return(result)
    } else if (method == "any") {
      result <- binary_rasters[[1]]
      for (i in 2:n_layers) {
        result <- result | binary_rasters[[i]]
      }
      return(result)
    } else if (method == "majority") {
      # Stack binary rasters
      binary_stack <- do.call(c, binary_rasters)
      # Sum them up
      sum_rast <- sum(binary_stack)
      # Majority vote
      return(sum_rast > (n_layers/2))
    } else {
      stop("Invalid method. Use 'all', 'any', 'majority', or 'individual'")
    }
  }
}

#----Implementation of DBSCAN for spatrasters----
# Define the optimized C++ functions for DBSCAN with spatial indexing
sourceCpp(code='
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

// Structure for grid-based spatial indexing
struct SpatialIndex {
  std::vector<std::vector<int>> cells;
  int rows;
  int cols;
  double cell_size;
  double min_x;
  double min_y;
  double max_x;
  double max_y;
};

// Function to create spatial index
// [[Rcpp::export]]
List create_spatial_index(NumericMatrix data, double cell_size) {
  int n = data.nrow();
  int d = data.ncol();
  
  // Find min/max values for each dimension
  NumericVector min_vals(d, R_PosInf);
  NumericVector max_vals(d, R_NegInf);
  
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < d; j++) {
      if (data(i, j) < min_vals[j]) min_vals[j] = data(i, j);
      if (data(i, j) > max_vals[j]) max_vals[j] = data(i, j);
    }
  }
  
  // We only use the first two dimensions (typically x,y) for spatial indexing
  int rows = std::ceil((max_vals[0] - min_vals[0]) / cell_size) + 1;
  int cols = std::ceil((max_vals[1] - min_vals[1]) / cell_size) + 1;
  
  // Create grid cells
  std::vector<std::vector<int>> cells(rows * cols);
  
  // Assign points to cells
  for (int i = 0; i < n; i++) {
    int row = std::floor((data(i, 0) - min_vals[0]) / cell_size);
    int col = std::floor((data(i, 1) - min_vals[1]) / cell_size);
    int cell_idx = row * cols + col;
    cells[cell_idx].push_back(i);
  }
  
  return List::create(
    Named("cells") = wrap(cells),
    Named("rows") = rows,
    Named("cols") = cols,
    Named("cell_size") = cell_size,
    Named("min_x") = min_vals[0],
    Named("min_y") = min_vals[1],
    Named("max_x") = max_vals[0],
    Named("max_y") = max_vals[1]
  );
}

// Calculate squared Euclidean distance between two points
double dist_squared(const NumericMatrix& data, int i, int j, int dims) {
  double sum = 0.0;
  for (int d = 0; d < dims; d++) {
    double diff = data(i, d) - data(j, d);
    sum += diff * diff;
  }
  return sum;
}

// DBSCAN implementation with spatial indexing
// [[Rcpp::export]]
List dbscan_cpp_optimized(NumericMatrix data, double eps, int minPts, List spatial_idx) {
  int n = data.nrow();
  int dims = data.ncol();
  double eps_squared = eps * eps;
  
  // Convert spatial index from R list to C++ struct
  std::vector<std::vector<int>> cells = as<std::vector<std::vector<int>>>(spatial_idx["cells"]);
  int rows = spatial_idx["rows"];
  int cols = spatial_idx["cols"];
  double cell_size = spatial_idx["cell_size"];
  double min_x = spatial_idx["min_x"];
  double min_y = spatial_idx["min_y"];
  
  IntegerVector cluster_ids(n, 0);  // 0 = unclassified
  int cluster_count = 0;
  
  // Function to get cell indices to search in
  auto get_neighbor_cells = [&](int row, int col) {
    std::vector<int> neighbor_cells;
    
    for (int r = std::max(0, row-1); r <= std::min(rows-1, row+1); r++) {
      for (int c = std::max(0, col-1); c <= std::min(cols-1, col+1); c++) {
        neighbor_cells.push_back(r * cols + c);
      }
    }
    
    return neighbor_cells;
  };
  
  // Function to find neighbors using spatial index
  auto get_neighbors = [&](int point_idx) {
    std::vector<int> neighbors;
    
    // Get the cell of the current point
    int row = std::floor((data(point_idx, 0) - min_x) / cell_size);
    int col = std::floor((data(point_idx, 1) - min_y) / cell_size);
    
    // Get neighboring cells to search
    std::vector<int> cell_indices = get_neighbor_cells(row, col);
    
    // Check points in all neighboring cells
    for (int cell_idx : cell_indices) {
      for (int candidate : cells[cell_idx]) {
        if (dist_squared(data, point_idx, candidate, dims) <= eps_squared) {
          neighbors.push_back(candidate);
        }
      }
    }
    
    return neighbors;
  };
  
  // Main DBSCAN algorithm
  for (int i = 0; i < n; i++) {
    if (cluster_ids[i] != 0) continue;  // Skip already processed points
    
    std::vector<int> neighbors = get_neighbors(i);
    
    if (neighbors.size() < (size_t)minPts) {
      cluster_ids[i] = -1;  // Mark as noise
      continue;
    }
    
    // Create new cluster
    cluster_count++;
    cluster_ids[i] = cluster_count;
    
    // Process neighbors (expand cluster)
    size_t j = 0;
    while (j < neighbors.size()) {
      int current = neighbors[j];
      
      // Update noise points to border points
      if (cluster_ids[current] == -1) {
        cluster_ids[current] = cluster_count;
      }
      
      // Skip already classified points
      if (cluster_ids[current] != 0) {
        j++;
        continue;
      }
      
      // Add to current cluster
      cluster_ids[current] = cluster_count;
      
      // Find neighbors of current point
      std::vector<int> current_neighbors = get_neighbors(current);
      
      // If core point, add its neighbors
      if (current_neighbors.size() >= (size_t)minPts) {
        for (int neighbor : current_neighbors) {
          // Check if this neighbor is already in our list
          bool found = false;
          for (size_t k = 0; k < neighbors.size(); k++) {
            if (neighbors[k] == neighbor) {
              found = true;
              break;
            }
          }
          if (!found) {
            neighbors.push_back(neighbor);
          }
        }
      }
      
      j++;
    }
  }
  
  return List::create(
    Named("cluster_ids") = wrap(cluster_ids),
    Named("n_clusters") = cluster_count
  );
}
')

# R function that uses the optimized C++ implementation
dbscan_raster_rcpp_optimized <- function(rast, eps, minPts, variables = 1:nlyr(rast)) {
  # Extract coordinates and values
  pts <- terra::as.data.frame(rast, xy = TRUE)
  coords <- pts[, c("x", "y")]
  
  # Select the variables to use for clustering
  if (is.numeric(variables)) {
    variables <- names(pts)[variables + 2]  # +2 to account for x,y columns
  }
  values <- pts[, variables, drop = FALSE]
  
  # Standardize feature values
  values_scaled <- scale(values)
  
  # Create matrix with spatial coordinates and feature values
  cluster_data <- cbind(coords, values_scaled)
  
  # Create spatial index (use cell size proportional to eps)
  cell_size <- eps * 1.5  # slightly larger than eps for efficiency
  spatial_idx <- create_spatial_index(as.matrix(cluster_data), cell_size)
  
  # Run the optimized C++ DBSCAN implementation
  result <- dbscan_cpp_optimized(as.matrix(cluster_data), eps, minPts, spatial_idx)
  
  # Add cluster IDs back to the original data
  pts$cluster <- result$cluster_ids
  
  # Create a new raster with cluster IDs
  result_raster <- terra::rast(rast, nlyrs=1)
  terra::values(result_raster) <- result$cluster_ids
  names(result_raster) <- "cluster"
  
  return(list(
    raster = result_raster,
    points = pts,
    n_clusters = result$n_clusters
  ))
}