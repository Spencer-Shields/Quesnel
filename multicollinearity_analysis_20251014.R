source('scene_SeparabilityAnalysis_20250925.R')
source('scene_global_stats_filtering_20250925.R')

#----get corrlation matrices for each dataset, block_id, and scene id----

gs_months = c(5,6,7,8,9,10)
results_meta_df_gs = results_meta_df |> 
  filter(month %in% gs_months) |>
  select(acquired, id, dataset, block_id, file_path) |>
  distinct()

clust = makeCluster(8)
plan('cluster', workers=clust)

corrs_l = future_lapply(unique(results_meta_df_gs$dataset), function(ds){
  df = results_meta_df_gs |> filter(dataset == ds)
  
  #generate correlation matrices for each scene id
  cm_l = lapply(1:length(unique(df$id)), function(i){
  # cm_l = lapply(1:10, function(i){
    cat(i,'/',length(unique(df$id)))
    x = unique(df$id)[i]
    df_ = df |> filter(id == x)
    
    #load raster for each block, merge into a mosaic, produce matrix
    m = lapply(df_$file_path, rast) |>
      sprc() |>
      merge() |>
      values(na.rm=T)
    # m = m[,colnames(m) %in% unique(results_meta_df$var_name)]
    cm = cor(m)
    return(cm)
  })
  # ca = abind(cm_l, along=3)
  ca = simplify2array(cm_l)
  mean_mat <- apply(ca, c(1, 2), mean, na.rm = TRUE)
  sd_mat = apply(ca, c(1,2), sd, na.rm=T)
  stat_array = simplify2array(list(mean_mat, sd_mat))
  dimnames(stat_array)[[3]] = c('mean', 'sd')
  return(stat_array)
})

stopCluster(clust)
plan('sequential')

names(corrs_l) = unique(results_meta_df_gs$dataset)

#remove unwanted features
corrs_l_sub <- lapply(corrs_l, function(a) {
  b = a[intersect(dimnames(a)[[1]], unique(results_meta_df$var_name)),
    intersect(dimnames(a)[[2]], unique(results_meta_df$var_name)),
    , drop = FALSE]
  b
})

library(patchwork)
library(ggcorrplot)

corplot_mean_l = lapply(1:length(corrs_l_sub), function(i){
  x = corrs_l_sub[[i]]
  y = x[,,'mean']
  y = round(y,2)
  ggcorrplot(y,
             hc.order = TRUE,
             type = "lower",
             lab = TRUE,
             title=names(corrs_l_sub)[i])
})

corrplot_mean = wrap_plots(corplot_mean_l, nrow=1)


ggcorrplot(corr,
           hc.order = TRUE,
           type = "lower",
           lab = TRUE)