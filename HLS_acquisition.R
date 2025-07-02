

#----install and load packages, configure session----
packages <- c(
  'tidyverse', 'earthdatalogin', 'rstac','imager','lubridate',
  'xts','dygraphs','leaflet','terra','pbapply'
  )

# Identify missing (not installed) packages
new.packages = packages[!(packages %in% installed.packages()[,"Package"])]

# Install new (not installed) packages
if(length(new.packages)) install.packages(new.packages, repos='http://cran.rstudio.com/', dependencies = TRUE) else print('All required packages are installed.')

#load packages
invisible(lapply(packages, library, character.only = TRUE))

#set up earthdata login
earthdatalogin::edl_netrc(
  username = 'spenshi', password = '3@rth D@t@ Plz'
  )

#define CMR-STAC endpoint
s = stac("https://cmr.earthdata.nasa.gov/stac/LPCLOUD/")


#----query data collections----
#more info: https://git.earthdata.nasa.gov/projects/LPDUR/repos/data-discovery---cmr-stac-api/browse

#search landsat 8/9 and sentinel 2a/b
HLS_col <- list("HLSS30_2.0", "HLSL30_2.0")

#load region of interest
roi_path = 'data/Quesnel_thinning/AOI_fullsite_wgs84.geojson'

roi = vect(roi_path)


# #visualize roi
# leaflet() %>% 
#   addPolygons(data = roi, fill = FALSE) %>% 
#   addProviderTiles(providers$Esri.WorldImagery) %>% 
#   addMiniMap(zoomLevelFixed = 5)

#get bounding box coordinates around ROI
roi_extent <- terra::ext(roi)
bbox <- c(roi_extent$xmin, roi_extent$ymin, roi_extent$xmax, roi_extent$ymax)

#define temporal range

start_dt = '2000-01-01T00:00:00Z'
end_dt = '2024-12-31T23:59:59Z'

roi_datetime <- paste0(start_dt,'/',end_dt)   # YYYY-MM-DDTHH:MM:SSZ/YYYY-MM-DDTHH:MM:SSZ

#submit query, get items
items <- s %>%
  stac_search(collections = HLS_col,
              bbox = bbox,
              datetime = roi_datetime,
              limit = 100) %>%
  post_request()
print(items)

#view item assets
assets = items_assets(items)
print(assets)

#view images
browse_image_url <- items$features[[1]]$assets$browse$href

browse <-load.image(browse_image_url)
plot(browse)


#put items into sf
sf_items <- items_as_sf(items)
View(sf_items)

# Retrieve Granule ID for each feature
granule_id <- sapply(items$features, function(feature) feature$id)
# Add as first column in sf_items
sf_items <- cbind(granule = granule_id, sf_items)
View(sf_items)

#mapping of band numbers to bands for landsat and sentinel
bands_df = tribble(
  ~land_asset, ~sent_asset, ~band,
  'Fmask', 'Fmask', 'Fmask',
  # 'metadata', 'metadata',
  'B01','B01', 'coastal_aerosol',
  'B02', 'B02', 'blue',
  'B03','B03', 'green',
  'B04','B04', 'red',
  NA, 'B05', 'rededge_1',
  NA, 'B06', 'rededge_2',
  NA, 'B07', 'rededge_3',
  NA, 'B08', 'nir_broad',
  'B05', 'B08A', 'nir_narrow',
  'B06', 'B11', 'swir_1',
  'B07', 'B12', 'swir_2',
  NA, 'B09', 'water_vapor',
  'B09', 'B10', 'cirrus',
  'B10', NA, 'tir_1',
  'B11', NA, 'tir_2'
)

#choose what bands you want to load
desired_bands = c('Fmask', 'nir_narrow', 'swir_1')

desired_land_bands = bands_df$land_asset[bands_df$band %in% desired_bands]
desired_sent_bands = bands_df$sent_asset[bands_df$band %in% desired_bands]

# Define a function to extract asset urls for selected bands
# # This also includes a check to ensure the correct bands are extracted
# # # depending on the collection (HLSL30 or HLSS30)

extract_asset_urls <- function(feature, land_bands, sent_bands) {
  collection_id <- feature$collection
  if (collection_id == "HLSS30_2.0") {
    bands = sent_bands
  } else if (collection_id == "HLSL30_2.0") {
    bands = land_bands}
  sapply(bands, function(band) feature$assets[[band]]$href)
}
# Retrieve Asset URLs for each feature using our extract_asset_urls function and transpose them to columns
asset_urls <- t(sapply(items$features, function(x)extract_asset_urls(x, land_bands = desired_land_bands, sent_bands = desired_sent_bands)))
View(asset_urls)

colnames(asset_urls) <- c('fmask', 'nir', 'swir')
sf_items <- cbind(sf_items, asset_urls)
# View(sf_items)
head(sf_items)


## Filter based on cloud cover
# sf_items1 <- sf_items[sf_items$eo.cloud_cover < 70,]

## Reset Row Indices
row.names(sf_items) <- NULL
# View(sf_items)
head(sf_items)


#----subset HLS COGs spatially and stack data layers----

#gdal configuration
setGDALconfig("GDAL_HTTP_UNSAFESSL", value = "YES")
setGDALconfig("GDAL_HTTP_COOKIEFILE", value = ".rcookies")
setGDALconfig("GDAL_HTTP_COOKIEJAR", value = ".rcookies")
setGDALconfig("GDAL_DISABLE_READDIR_ON_OPEN", value = "EMPTY_DIR")
setGDALconfig("CPL_VSIL_CURL_ALLOWED_EXTENSIONS", value = "TIF")

# This function reads an HLS scene from a URL, applies the scale factor if necessary, and optionally crops and
# masks the scene based on a polygon. It requries the above GDAL configurations and a .netrc file. A .netrc
# can be created by running `earthdatalogin::edl_netrc()`.
open_hls <- function(url, roi = NULL) {
  # Add VSICURL prefix
  url <- paste0('/vsicurl/', url)
  # Retrieve metadata
  meta <- describe(url)
  # Check if dataset is Quality Layer (Fmask) - no scaling this asset (int8 datatype)
  is_fmask <- any(grep("Fmask", meta))
  # Check if Scale is present in band metadata
  will_autoscale <- any(grep("Scale:", meta))
  # Read the raster
  r <- rast(url)
  # Apply Scale Factor if necessary
  if (!will_autoscale && !is_fmask){
    print(paste("No scale factor found in band metadata. Applying scale factor of 0.0001 to", basename(url)))
    r <- r * 0.0001
  }
  # Crop if roi specified
  if (!is.null(roi)){
    # Reproject roi to match crs of r
    roi_reproj <- project(roi, crs(r))
    r <- mask(crop(r, roi_reproj), roi_reproj)
  }
  return(r)
}

# Test opening and cropping a single raster
nir <- open_hls(sf_items$nir[1]
                , roi
                )
#make sure roi is same crs as raster data
if(crs(roi) != crs(nir)){
  roi = project(roi, nir)
  }
# plot(nir)
# plot(roi, add=T)

library(tidyterra)
ggplot()+
  geom_spatraster(data = nir)+
  geom_spatvector(data=roi, fill = NA, color = 'red')

# leaflet() %>%
#   addPolygons(data = roi, fill = FALSE, color = 'red') %>%
#   # addRasterImage(nir) |>
#   addProviderTiles(providers$Esri.WorldImagery) %>%
#   addRasterImage(nir) |>
#   # addLegend()|>
#   addMiniMap(zoomLevelFixed = 5)

#stack bands by type
{
stacks_l = pblapply(colnames(asset_urls), function(x){
  col = sf_items[[x]]
  stack = pblapply(col, function(y){
    cat('\rProcessing',y)
    r = open_hls(y, roi = roi)
    return(y)
  } 
  # open_hls, roi=roi
  )
  cat(x,'done')
  return(stack)
})
names(stacks_l) = colnames(asset_urls)
  }

combined_rasts_l = pblapply(1:length(stacks_l[[1]]), function(i){
  rl = lapply(stacks_l, function(x)x[[i]])
  r = rast(rl)
  names(r) = names(stacks_l)
  return(r)
})
names(combined_rasts_l) = sf_items$datetime

#----Process HLS data----

#define NBR function
nbr.vi = function(swir, nir){
  nbr = (swir - nir)/(swir + nir)
  return(nbr)
}

#calculate nbr for all scenes
nbr_l = pblapply(1:nrow(sf_items), function(i){
  nbr = nbr.vi(swir = stacks_l$swir[[i]], nir = stacks_l$nir[[i]])
  return(nbr)
  })
names(nbr_l) = sf_items$datetime #make date layernames

for(i in 2:length(nbr_l)){ #make crs consistent if it's not
  if(crs(nbr_l[[i]]) != crs(nbr_l[[1]])){
    nbr_fixed = project(nbr_l[[i]], nbr_l[[1]])
    nbr_l[[i]] = nbr_fixed
  }
}

nbr_stack = rast(nbr_l) #list to spatraster


