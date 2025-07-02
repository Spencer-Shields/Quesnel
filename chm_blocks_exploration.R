library(tidyverse)
library(terra)
library(sf)
library(tools)

#load thinning blocks
blocks_path = 'data/Quesnel_thinning/12l_12n_bdy.gpkg'
st_layers(blocks_path) #1 layer

blocks = st_read(blocks_path)
ggplot()+
  geom_sf(data = blocks, aes(color = BLOCKNUM))

#load chm data

chm_pre_dir = 'data/Quesnel_thinning/chm_pre_Jul2021'
chm_pre_files = list.files(chm_pre_dir, full.names = T)

chm_pre_names = list.files(chm_pre_dir) %>% file_path_sans_ext()

chm_pre_l = lapply(chm_pre_files, rast)
names(chm_pre_l) = list.files(chm_pre_dir) %>% file_path_sans_ext()

chm_pre_vrt = vrt(chm_pre_files)

#reproject blocks

blocks = st_transform(blocks, crs = crs(chm_pre_vrt))

plot(blocks)
lapply(chm_pre_l, function(x)plot(x, add = T))

#get bounding box for all data
b_blocks = st_bbox(blocks) %>% 
  st_as_sfc() %>%
  st_sf()

b_chms_pre_l = lapply(chm_pre_l, function(x){
  st_bbox(x) %>% 
    st_as_sfc() %>%
    st_sf()})
b_chms_pre = bind_rows(b_chms_pre_l) %>% 
  st_union() %>%
  st_sf()


b = bind_rows(list(b_chms_pre, b_blocks)) %>% 
  st_union() %>%
  st_bbox() %>%
  st_as_sfc()%>%
  st_sf()

ggplot()+
  geom_sf(data = b)+
  geom_sf(data = blocks, aes(fill = BLOCKNUM)) +
  geom_sf(data = b_chms_pre)

b = st_transform(b, crs= 'WGS84')

#export bbox
b_filename = 'data/Quesnel_thinning/AOI_fullsite_wgs84.geojson'
if(!file.exists(b_filename)){st_write(b, b_filename)}
st_read(b_filename)

#----make map of study area with chm data----

library(tmap)
library(tidyterra)

blocks_v = vect(blocks)
block_12nt3 = blocks_v[blocks_v$BLOCKNUM=='12N_T3']

i = 7

chm = c(chm_pre_l[[i]], chm_post_l[[i]])

ggplot() +
  geom_spatraster(data = chm_pre_l[[7]], aes(fill = Z)) +
  scale_fill_viridis_c(option='turbo', aesthetics = 'fill')+
  geom_spatvector(data = block_12nt3, color = 'red', fill=NA)
  


