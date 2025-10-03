source('scene_setup_preprocessing.R')

#---site map plotting----

#----base plot of thinning blocks
library(RColorBrewer)

blockplot =  #geom for the thinning blocks
  # ggplot()+
  list(
  geom_sf(data = blocks, aes(color = BLOCKNUM), fill = NA, linewidth = 1),
  scale_color_manual(values = brewer.pal(8, 'Set1')),
  labs(color = "Block ID"),
  theme_minimal()
  )

blocks |> st_area() |> sum() |> units::set_units("hectares") #get total area of thinning blocks


#----full area plot

aoi = vect('data/Quesnel_thinning/AOI_fullsite_wgs84.geojson') |> project(vect(blocks))

library(rnaturalearth)
library(rnaturalearthdata)
library(rnaturalearthhires)
library(osmdata)
library(ggspatial)
library(sf)

# Get Canada boundary
# canada <- ne_countries(country = "canada", returnclass = "sf", scale = "large")

# Download 10m resolution (highest available) admin boundaries
provinces_hires <- ne_download(scale = "large", 
                               type = "admin_1_states_provinces", 
                               category = "cultural") %>%
  filter(admin == "Canada") |>
  filter(name == 'British Columbia')

canada <- list(
  geom_sf(data = provinces_hires, 
          fill = "lightgrey", 
          color = "white", 
          size = 0.7) 
  # ,geom_sf_text(data = provinces_hires, 
  #              aes(label = name), 
  #              size = 3.5, 
  #              fontface = "bold",
  #              check_overlap = TRUE) ,
  # labs(title = "Canada - High Resolution Provincial Boundaries") ,
  # ,theme_minimal()
  )

# canada

# Convert your AOI to sf
aoi_sf <- st_as_sf(aoi) |> st_centroid()
aoi_sf$label = 'Study site location'

# Create a buffer around your AOI
aoi_buffer <- st_buffer(aoi_sf, dist = 450000)  # 100km buffer for detailed map
bbox <- st_bbox(aoi_buffer)

# # Get roads
# roads <- opq(bbox = bbox) %>%
#   add_osm_feature(key = "highway", 
#                   value = c("motorway"
#                             , "trunk"
#                             #, "primary", "secondary", "tertiary"
#                             )
#                   ) %>%
#   osmdata_sf()

# # Get water bodies
# water <- opq(bbox = bbox) %>%
#   add_osm_feature(key = "natural", value = c("water")) %>%
#   osmdata_sf()
# 
# lakes <- opq(bbox = bbox) %>%
#   add_osm_feature(key = "water", value = c("lake", "reservoir")) %>%
#   osmdata_sf()
# 
# # Get rivers
# rivers <- opq(bbox = bbox) %>%
#   add_osm_feature(key = "waterway", value = c("river", "stream")) %>%
#   osmdata_sf()

# Get populated places
towns <- opq(bbox = bbox) %>%
  add_osm_feature(key = "place", value = c("city"
                                           # , "town"
                                           )
                  ) %>%
  osmdata_sf()


bb = st_bbox(provinces_hires)

# Create detailed map
detailed_map <- ggplot() +
  # geom_sf(data = canada, fill = "lightgray", color = "black", size = 0.3) +
  canada+
  # # Add water bodies
  # geom_sf(data = water$osm_polygons, fill = "lightblue", color = "blue", size = 0.2) +
  # geom_sf(data = lakes$osm_polygons, fill = "lightblue", color = "blue", size = 0.2) +
  # geom_sf(data = rivers$osm_lines, color = "blue", size = 0.5) +
  # # Add roads
  # geom_sf(data = roads$osm_lines, color = "darkred", size = 0.5, alpha = 0.5) +
  # Add towns
  geom_sf(data = towns$osm_points |> filter(name %in% c('Vancouver', 'Quesnel', 'Prince George')), color = "black", size = 2) +
  geom_sf_text(data = towns$osm_points |> filter(name %in% c('Vancouver', 'Quesnel', 'Prince George')), aes(label = name), 
               size = 3, nudge_y = 0.3, nudge_x = 1.5, color = "black", fontface = "bold") +
  # Add your AOI
  geom_sf(data = aoi_sf, fill='red', color = "darkred", size = 2, alpha = 0.7) +
  geom_sf_text(data = aoi_sf, color = 'darkred', label = aoi_sf$label, size = 3, nudge_x = -2.5)+
  # scale_fill_manual(c(red')+
  # annotation_scale()+
  # annotation_north_arrow()+
  # coord_sf(xlim = c(bbox[1], bbox[3]), ylim = c(bbox[2], bbox[4])) +
  # labs(title = "Detailed Study Area Location near Quesnel, BC") +
  theme_void()+
  labs(title = 'A')+
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank()
  )
detailed_map






#----planetscope plot

ps_raw_list = list.files(raw_dir, recursive = T, full.names = T, pattern = '\\.tif$')
ps_raw_scenes_list = ps_raw_list[!str_detect(ps_raw_list, 'udm2')]
ps_raw_ids = ps_raw_scenes_list |> basename()


# ps_2021_scene_file = ps_raw_scenes_list[str_detect(ps_raw_scenes_list,"20210712_182629_31_2460")]
ps_2021_scene_file = ps_raw_scenes_list[3]
ps_2021_scene = rast(ps_2021_scene_file) |> project(y = crs(aoi)) |> crop(aoi)

mcv = 255 #max color value (default 255, adjust to moderate saturation vs darkness)

ps2021_p = ggplot()+
  geom_spatraster_rgb(data = ps_2021_scene, r = 6, g = 4, b = 2
                      , stretch = 'lin'
                      ,max_col_value = 300
                      ) +
  blockplot +
  labs(title='B', subtitle =  str_remove(basename(ps_2021_scene_file), "_3B.*"))
# ps2021_p

ps_2024_scene_file = ps_raw_scenes_list[828]
ps_2024_scene = rast(ps_2024_scene_file) |> project(y = crs(aoi)) |> crop(aoi)

ps2024_p = ggplot()+
  geom_spatraster_rgb(data = ps_2024_scene, r = 6, g = 4, b = 2
                      , stretch = 'lin'
                      ,max_col_value = 300
  ) +
  blockplot +
  labs(title = 'B', subtitle= str_remove(basename(ps_2024_scene_file), "_3B.*"))
# ps2024_p

library(patchwork)
ps2021_p+ps2024_p+plot_layout(guides = 'collect',nrow=2)

#----planetscope next to location

library(patchwork)
(detailed_map+theme(plot.background = element_rect(color = "black", fill = NA, size = 1))) + (ps2021_p+theme(plot.background = element_rect(color = "black", fill = NA, size = 1)))


#----gifs----
{
#----make raster stack for non-normalized, z_normalized, and zrobust
  
  
}