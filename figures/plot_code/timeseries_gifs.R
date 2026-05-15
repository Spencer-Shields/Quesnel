source('scene_SeparabilityAnalysis_20251103.R')
library(gganimate)
library(ggtext)
library(ggbreak)
library(magick)


#set parameters
feat = 'NDVI'
block = '12N_T3'
datasets = c('NN', 'SM', 'Z', 'Zr')

#set up data
{
  feats_to_plot = feat
  
  blocks_to_plot = block
  
  datasets_to_plot = datasets
  
  subset_df_ = results_meta_df |>
    filter(
      # block_pixel_stratum == 'Total',
      var_name %in% feats_to_plot
      ,block_id %in% blocks_to_plot
      # ,id %in% ids_in_all_blocks
    ) |>
    mutate(dataset = case_match(dataset, 
                                'Non-normalized' ~ 'NN', 
                                'Zrobust' ~ 'Zr',
                                .default = dataset)) |>
    mutate(block_pixel_stratum = str_replace_all(block_pixel_stratum, '_', ' '))|>
    mutate(block_pixel_stratum = factor(block_pixel_stratum, levels = c("Not thinned", "Thinned")  # desired order
    )) |>
    filter(dataset %in% datasets)
  
  subset_df$acquired = ymd_hms(subset_df$acquired)
  
  harvest_dates_df <- harvest_dates_df |>
    mutate(
      harvest_start_date = as.POSIXct(harvest_start_date),
      harvest_finish_date = as.POSIXct(harvest_finish_date)
    )
}

#make directory
gif_dir = file.path('figures','gifs')
dir.check(gif_dir)

#----raster time series gifs----

subset_df = subset_df_ |> 
  filter(block_pixel_stratum == 'Total') |>
  select(acquired, month)

#individual rasters, full time series
{
  #set length of gif
  duration_full = 10 #seconds
  
  #define individual gif files, make if they don't exist
  rast_gif_files = file.path(gif_dir, paste0('rast_',datasets,'_',block,'_',feat,'_',duration_full,'.gif'))
  lapply(rast_gif_files, function(gf){
    if(!file.exists(gf)){
      
      d = datasets[which(str_detect(gf, paste0(datasets,'_')))] #add underscore to differentiate between Z and Zr
      
      cat('\nMaking full time series raster gif,',d)
      cat('\n---Loading rasters',d)
      df = subset_df |> filter(dataset==d) |> distinct()
      r = lapply(df$file_path, function(f){
        # cat('\r',f); 
        rast(f,lyrs=feat)}) |> rast()
      # names(r) = df$acquired
      names(r) = paste0(
        "<b>", format(df$acquired, "%Y-%m-%d"), "</b> ",
        format(df$acquired, "%H:%M:%S")
      )
      
      cat('\n---Making gif',d)
      # gifl = pblapply(rl,function(r){
      p = ggplot()+
        geom_spatraster(data = r)+
        scale_fill_viridis_c(na.value = 'transparent')+
        transition_manual(lyr)+
        theme_void()+
        labs(subtitle = '{current_frame}'
             ,title=paste0('<b>',d,'</b>')
        )+
        theme(plot.subtitle = element_markdown()
              ,plot.title = element_markdown())
      gif = gganimate::animate(p,duration=duration_full)
      cat('\n---Saving gif',d)
      anim_save(gf, gif)
    }
  })
}

#individual rasters, growing season
{
  #set length of gif
  duration_gs = 10 #seconds
  
  #filter data for only gs months
  subset_df_gs = subset_df |> filter(month %in% gs_months)
  
  #define individual gif files, make if they don't exist
  rast_gif_files_gs = file.path(gif_dir, paste0('rast_GS_',datasets,'_',block,'_',feat,'_',duration_gs,'.gif'))
  lapply(rast_gif_files_gs, function(gf){
    if(!file.exists(gf)){
      
      d = datasets[which(str_detect(gf, paste0(datasets,'_')))] #add underscore to differentiate between Z and Zr
      
      cat('\nMaking growing season raster gif',d)
      cat('\n---Loading raster',d)
      df = subset_df_gs |> filter(dataset==d) |> distinct()
      r = pblapply(df$file_path, function(f){
        cat('\r',f);
        rast(f,lyrs=feat)}) |> rast()
      # names(r) = df$acquired
      names(r) = paste0(
        "<b>", format(df$acquired, "%Y-%m-%d"), "</b> ",
        format(df$acquired, "%H:%M:%S")
      )
      
      cat('\n---Making gif',d)
      # gifl = pblapply(rl,function(r){
      p = ggplot()+
        geom_spatraster(data = r)+
        scale_fill_viridis_c(na.value = 'transparent')+
        transition_manual(lyr)+
        theme_void()+
        labs(subtitle = '{current_frame}'
             ,title=paste0('<b>',d,'</b>')
        )+
        theme(plot.subtitle = element_markdown()
              ,plot.title = element_markdown())
      gif = gganimate::animate(p,duration=duration_gs)
      
      cat('\n---Saving gif',d)
      anim_save(gf, gif)
    }
  })
}

#reload growing season raster gifs as magick objects
rast_l_magick = lapply(rast_gif_files_gs, image_read)
names(rast_l_magick) = datasets

#----scatter plot time series gifs----

subset_df = subset_df_ |> 
  filter(block_pixel_stratum != 'Total') |>
  mutate(acquired = as.POSIXct(acquired))

#raster parameters to set plot dimensions
ii = image_info(rast_l_magick[[1]])
rwidth = unique(ii$width)
rheight = unique(ii$height)*(2/3)

pixel_strata_color_map = c('Not thinned' = "#00BFC4", 'Thinned' = "#F8766D")

#growing season plot
{
  
  scat_plot_l_gs = pblapply(unique(subset_df$dataset), function(d){
    
    # 1. Build the static base plot (no transition)
    static_base <- ggplot(
      subset_df |> filter(month %in% gs_months, dataset == d)
    ) +
      geom_ribbon(
        aes(x = acquired, ymin = X25., ymax = X75.,
            fill = block_pixel_stratum, color = block_pixel_stratum),
        alpha = 0.3
      ) +
      geom_point(
        aes(x = acquired, y = median, color = block_pixel_stratum),
        alpha = 0.5
      ) +
      scale_color_manual(name = "Pixel stratum", values = pixel_strata_color_map) +
      scale_fill_manual(name = "Pixel stratum",  values = pixel_strata_color_map) +
      scale_x_break(c(as.POSIXct('2020-11-01 UTC'), as.POSIXct('2021-04-30 UTC'))) +
      scale_x_break(c(as.POSIXct('2021-11-01 UTC'), as.POSIXct('2022-04-30 UTC'))) +
      scale_x_break(c(as.POSIXct('2022-11-01 UTC'), as.POSIXct('2023-04-30 UTC'))) +
      scale_x_break(c(as.POSIXct('2023-11-01 UTC'), as.POSIXct('2024-04-30 UTC'))) +
      scale_x_datetime(
        date_labels = "%b %Y",
        limits = c(as.POSIXct('2020-05-01 UTC'), as.POSIXct('2024-10-31 UTC'))
      ) +
      geom_vline(
        data = harvest_dates_df[harvest_dates_df$block_id %in% unique(subset_df$block_id), ],
        aes(xintercept = harvest_finish_date, linetype = 'Harvest end date')
      ) +
      geom_vline(
        data = harvest_dates_df[harvest_dates_df$block_id %in% unique(subset_df$block_id), ],
        aes(xintercept = harvest_start_date, linetype = 'Harvest start date')
      ) +
      scale_linetype_manual(
        name = "Harvest date",
        values = c("Harvest start date" = "dotted", "Harvest end date" = "dashed")
      ) +
      labs(x = 'Acquisition date', y = unique(subset_df$var_name)) +
      theme_classic() +
      theme(
        axis.text.x.top  = element_blank(),
        axis.ticks.x.top = element_blank(),
        axis.text.x = element_text(angle = 90)
        ,legend.position = 'none' #remove legend
      )
    
    # 2. Build a separate data frame for the moving line — one row per frame
    frame_times <- subset_df |>
      filter(month %in% gs_months, dataset == d) |>
      distinct(acquired) |>
      arrange(acquired) |>
      mutate(frame = row_number())   # <-- key: explicit frame index
    
    # 3. Add the animated layer using the SEPARATE data frame + transition_manual
    animated_plot <- static_base +
      geom_vline(
        data    = frame_times,          # own data, not inherited
        mapping = aes(xintercept = acquired, group = frame),  # group = frame prevents filtering bleed
        color     = "blue",
        linewidth = 1
      ) +
      transition_manual(frame)          # transition on integer frame, not a shared column
    
    # 4. Render
    animate(
      animated_plot,
      nframes  = nrow(frame_times),
      duration = duration_gs
      ,height = rheight,
      width=rwidth
    )
  })
  names(scat_plot_l_gs) = unique(subset_df$dataset)
  
}

#----arrange rasters and scatter plots----

#convert scatter plot gganimate gifs to magick objects
scat_plot_l_magick = pblapply(scat_plot_l_gs, image_read)

#arrange rasters with scatterplots

rast_scats = pblapply(datasets, function(d){
  r = rast_l_magick[[d]]
  s = scat_plot_l_magick[[d]]
  combined <- image_append(c(r, s), stack = TRUE)
})









#growing season
{
  #combine raster gifs with magick
  {
    library(magick)
    
    # Read all GIFs as magick objects
    gif_list <- lapply(rast_gif_files_gs, image_read)
    
    # Check number of frames in each (should be the same)
    sapply(gif_list, length)
    
    # Combine all 4 GIFs horizontally for each frame
    n_frames <- length(gif_list[[1]])
    
    combined_frames <- lapply(1:n_frames, function(i) {
      # Get frame i from each gif
      frame_list <- lapply(gif_list, function(gif) gif[i])
      # Append horizontally
      image_append(image_join(frame_list), stack = FALSE)
    })
    
    # Optional: Add labels above each GIF
    combined_frames_labeled <- lapply(1:n_frames, function(i) {
      frame_list <- lapply(seq_along(gif_list), function(j) {
        # Add text label
        img <- gif_list[[j]][i]
        img = image_trim(img)
        image_annotate(img, names(gifl_gs)[j], 
                       size = 20, color = "white", 
                       location = "+10+10")
      })
      image_append(image_join(frame_list), stack = FALSE)
    })
    combined_gif_labeled <- image_join(combined_frames_labeled)
    
    # Join all frames into a single animation
    combined_gif <- image_join(combined_frames)
    combined_f = paste0(gif_dir,'/combined_',unique(subset_df_gs$block_id),'_',unique(subset_df_gs$var_name),'.gif')
    
  }
  
}
