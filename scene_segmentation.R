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

otsu_l = pbsapply(names(chm_l), function(x){
  d = change_df %>%
    filter(block == x) 
  t = otsuThresholdCpp(values = d[['vals.Z']], bins = 256)
})

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
