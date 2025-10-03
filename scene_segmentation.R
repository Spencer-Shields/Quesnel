source('scene_setup_preprocessing.R')

#----load lidar data----

chm_l = pblapply(lidar_files, rast)
names(chm_l) = block_ids

pblapply(1:length(chm_l), function(i){
  v = values(chm_l[[i]], dataframe=T) |> mutate(block_id = names(chm_l)[i])
  return(v)
}) |> 
  bind_rows() |>
  ggplot(aes(x = Z))+
  geom_density(alpha = 1)+
  facet_grid(rows = vars(block_id), scales = 'free_y')

#----otsu thresholding for lidar data resampled to 3m----

change_df = pblapply(1:length(chm_l), function(i){
  d = data.table(
    vals = values(chm_l[[i]])
    ,
    block = names(chm_l)[i]
  )
}) |>
data.table::rbindlist()

otsu3m_l = pbsapply(names(chm_l), function(x){
  d = change_df %>%
    filter(block == x) 
  t = otsuThresholdCpp(values = d[['vals.Z']], bins = 256)
})

otsu3m_df = tibble(block_id = names(otsu3m_l),
                   chm_change_otsu_threshold = otsu3m_l)

otsu3m_path = "data/Quesnel_thinning/Otsu_3m_change_thresholds.csv"
if(!file.exists(otsu3m_path)){write_csv(otsu3m_df, otsu3m_path)}


#----plot Otsu segmented height change rasters----

#thresholds calculated for 3m raster
chm_binarized_Otsu3m_l = pblapply(block_ids, function(b){
  chm = chm_l[[b]]
  otsu_t = otsu3m_df$chm_change_otsu_threshold[otsu3m_df$block_id==b]
  chm_bin = ifel(chm < otsu_t, 1, 0)
  return(chm_bin)
})
names(chm_binarized_Otsu3m_l) = block_ids
list_plot(chm_binarized_Otsu3m_l)

#thresholds calculated for 0.25m raster
chm_binarized_Otsu0.25m_l = pblapply(block_ids, function(b){
  chm = chm_l[[b]]
  otsu_t = otsu_df$chm_change_otsu_threshold[otsu3m_df$block_id==b]
  chm_bin = ifel(chm < otsu_t, 1, 0)
  return(chm_bin)
})
names(chm_binarized_Otsu0.25m_l) = block_ids
list_plot(chm_binarized_Otsu0.25m_l)


#----plot density plots for lidar height change----

pblapply(1:length(chm_l), function(i){
  v = values(chm_l[[i]], dataframe=T) |> mutate(block_id = names(chm_l)[i])
  return(v)
}) |> 
  bind_rows() |>
  left_join(otsu3m_df) |>
  mutate(otsu_t = round(chm_change_otsu_threshold,2))|>
  
  group_by(block_id) |> #get density for height aesthetics
  mutate(max_density = max(density(Z,na.rm=T)$y)) |>
  ungroup() |>
  
  ggplot(aes(x = Z))+
  geom_density(alpha = 1)+
  geom_vline(aes(xintercept = chm_change_otsu_threshold
                 # , linetype = 'Otsu threshold'
  )
  ,color = 'darkred'
  )+
  geom_text(
    aes(x = chm_change_otsu_threshold + 1, y = 0.9 * max_density,
        label = paste0(
          # "Thinning threshold: ", 
          otsu_t)),
    color = 'darkred',
    hjust = 0
    # ,vjust = -0.5
    ,size=2.5
  ) +
  facet_grid(rows = vars(block_id), scales = 'free_y')+
  theme_minimal()

#----get stats for timeseries of 


#----dbscan segmentation----

library(dbscan)

r = chm_l$`12N_T3`
d = r |>
  values(dataframe=T)
d = cbind(d, xyFromCell(r, 1:ncell(r))) |>
  scale() |>
  as.data.frame()
d[['cell_id']] = 1:nrow(d)

d_noNA = d |> na.omit()

# db = hdbscan(d_noNA |> select(Z),
#             # eps = 0.0001,
#             minPts = 10
#             # *global(r, 'notNA')[1,1],
#             ,verbose = T
#             )
# 
# d_noNA[['db_cluster']] = db$cluster
# 
# d = left_join(d, d_noNA)
# 
# r_clust = r
# values(r_clust) = NA
# values(r_clust) = d$db_cluster
# plot(r_clust)


sc = supercells::supercells(r, step = 2, compactness = 0.1, avg_fun = 'median')

ggplot() +
  geom_sf(data = sc, aes(fill = Z), color = NA) +  # remove borders
  # scale_fill_viridis_c() +  # continuous color scale
  theme_minimal()
  
db = dbscan(sc |> 
               st_drop_geometry()|> 
               select(Z) |>
               as.matrix()
            , eps = 0.1
             , minPts = 5
            # , verbose=T
            )

sc$db_cluster_id = db$cluster
ggplot() +
  geom_sf(data = sc, aes(fill = db_cluster_id), color = NA) +  # remove borders
  # scale_fill_viridis_c() +  # continuous color scale
  theme_minimal()


km = kmeans(sc |> 
              st_drop_geometry()|> 
              select(Z) 
            , centers = 2)
sc$km_cluster_id = km$cluster
ggplot() +
  geom_sf(data = sc, aes(fill = km_cluster_id), color = NA) +  # remove borders
  # scale_fill_viridis_c() +  # continuous color scale
  theme_minimal()

ggplot(sc |> 
         st_drop_geometry() |> 
         mutate(km_cluster_id = as.factor(km_cluster_id)))+
  geom_boxplot(aes(fill = km_cluster_id, y = Z))
  

# sc$height_change=exactextractr::exact_extract(r, sc, fun = 'weighted_mean')

r2 = rast(list.files("data/planet_scenes/Z_REDO/12N_T3"
                , recursive = T, full.names = T)[562])
sc2 = supercells::supercells(x = r2 |> scale() |> select(BI), step = 2, compactness = 0.1, avg_fun = 'median')
ggplot() +
  geom_sf(data = sc2, aes(fill = BI), color = NA) +  # remove borders
  # scale_fill_viridis_c() +  # continuous color scale
  theme_minimal()

library(units)
sc2$area = st_area(sc2)
sc2$n_pixels = as.numeric(sc2$area/(3^2))
hist(sc2$n_pixels, nclass=30)

ggplot(data=sc2)+
  geom_sf(aes(fill=n_pixels), color=NA)+
  scale_fill_viridis_c()

boxplot(sc2$n_pixels)

km2 = kmeans(sc2 |> 
              st_drop_geometry()|> 
              select(BI) 
            , centers = 2)
sc2$km_cluster_id = km2$cluster
ggplot() +
  geom_sf(data = sc2, aes(fill = km_cluster_id), color = NA) +  # remove borders
  # scale_fill_viridis_c() +  # continuous color scale
  theme_minimal()

plot(r2$nir)
ggR(img=r2, layer = 31)
hist(r2[['nir']])
