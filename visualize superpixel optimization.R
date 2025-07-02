library(tidyverse)
library(arrow)
library(GGally)
library(viridis)
# library(Polychrome)
library(paletteer)

pth = "D:/Quesnel/data/planet_basemaps/global_monthly/superpixel_optimization_results1"

fl = list.files(pth, full.names = T)

df = lapply(fl, read_feather) |> bind_rows()

df_s = df |> filter(fold == 'overall')

ggplot(df |> filter(fold =='overall'
                    , target=='canopy_cover'
                    , block_id=='12N_V1'
                    ))+
  geom_point(aes(x = sp_compactness, y = sp_size, shape = sp_avg_fn, color = r2_test)) +
  scale_color_viridis_c()+
  # facet_grid(rows=vars(dataset), cols=vars(block_id))+
  scale_x_log10()

library(scatterplot3d)

d_v1 = df |> filter(fold =='overall'
                    , target=='canopy_cover'
                    , dataset=='Non-normalized'
                    , sp_avg_fn=='mean'
                    , block_id=='12N_V1'
)

scatterplot3d(x = d_v1$sp_compactness, y = d_v1$sp_size, z = d_v1$r2_test)

library(plotly)
plot_ly(data = d_v1, x = ~sp_compactness, y = ~sp_size, z = ~r2_test,
        type = "scatter3d", mode = "markers")

across_blocks_summ = df_s |>
  group_by(target, sp_avg_fn, sp_compactness, sp_size, dataset) |>
  summarise(mean_r2_test = mean(r2_test)
            , mean_rmse_test = mean(rmse_test)
            , mean_r2_train = mean(r2_train)
            , mean_rmse_train = mean(rmse_train))

plot_ly(data = across_blocks_summ, x = ~sp_compactness, y = ~sp_size, z = ~mean_r2_test,
        type = "scatter3d", mode = "markers")

plot_ly(
  data = across_blocks_summ |> filter(target=='canopy_cover'),
  x = ~sp_compactness,
  y = ~sp_size,
  z = ~mean_r2_test,
  type = "scatter3d",
  mode = "markers",
  color = ~sp_avg_fn,  # works with categorical variables
  colors = c("blue", "red")  # choose any two colors
)

#plot showing average r2 values for each dataset
ggplot(across_blocks_summ |> filter(sp_avg_fn!='median'))+
  # geom_point(aes(x = sp_size, y = mean_r2_test))+
  geom_line(aes(x = sp_size, y = mean_r2_test, linetype=dataset))+
  facet_grid(cols = vars(sp_compactness), rows=vars(target))

#plot showing block_id spread


block_ids <- unique(df_s$block_id) #get polychrome color palette
n_colors <- length(block_ids)
block_colors <- createPalette(n_colors
                              # , seedcolors = c("#000000", "#FFFFFF")
                              , seedcolors = c('#eb3434','#3d34eb')
                              )
names(block_colors) <- block_ids

ggplot()+
  geom_line(data = df_s |> 
               filter(sp_avg_fn!='median'
                      ,sp_size <= 15
                      )
             , aes(x = sp_size, y = r2_test, color=block_id, linetype=dataset)
             ,alpha=0.5
            )+
  geom_line(data = across_blocks_summ |> 
              filter(sp_avg_fn!='median'
                     ,sp_size <= 15
              )
            , aes(x = sp_size, y = mean_r2_test,linetype = dataset)
            ,linewidth = 0.8
            )+
  # scale_color_manual(values=block_colors)+
  scale_color_manual(values = as.character(paletteer_d("ggsci::default_nejm")[1:length(unique(df_s$block_id))])) +
  ylim(c(-0.05,0.81))+
  facet_grid(cols = vars(sp_compactness), rows=vars(target))+
  theme_minimal()
