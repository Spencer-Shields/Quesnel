source('scene_SeparabilityAnalysis_20251103.R')
library(gganimate)
library(ggtext)
library(ggbreak)
library(magick)


#set parameters
feat = 'NDVI'
block = '12L_C7'
datasets = c('NN', 'SM', 'Z', 'Zr')

color_palette = 'turbo'

#make a video for selected blocks
blocks_to_animate = c(
  block_ids
)

for(block in blocks_to_animate){
  
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
    
    subset_df_$acquired = ymd_hms(subset_df_$acquired)
    
    harvest_dates_df <- harvest_dates_df |>
      mutate(
        harvest_start_date = as.POSIXct(harvest_start_date),
        harvest_finish_date = as.POSIXct(harvest_finish_date)
      )
    
    subset_df_gs = subset_df_ |> filter(month %in% gs_months)
    
    gs_dates = unique(subset_df_gs$acquired)
  }
  
  #check that the same dates exist for all datasets
  all(sapply(datasets, function(d){
    print(length(unique(subset_df_ |> filter(dataset==d))))
  }))
  
  
  #make directory to save frames
  gif_dir = file.path('figures','gifs')
  dir.check(gif_dir)
  
  
  #----generate frames with raster and scatterplot----
  target_fps = 20
  combined_animation_f = file.path(gif_dir, paste0('combinedanimation_gs_',block,'_',feat,'_',target_fps,'fps_',color_palette,'yaxlabel.mp4'))
  if(!file.exists(combined_animation_f)){
    cat('Plotting, reshaping frames')
    ds_frames = pblapply(datasets, function(ds){
      
      #get dataframe
      df = subset_df_gs |> 
        filter(dataset==ds
               , !is.na(block_pixel_stratum)
        )
      
      #get raster stack
      rs = pblapply(unique(df$file_path), function(f){
        cat('\r',f);
        rast(f,lyrs=feat)}) |> rast()
      # names(rs) = format(gs_dates, "%Y-%m-%d %H:%M:%S", tz = "UTC")
      names(rs) = gs_dates
      # names(rs) = paste0(
      #   "<b>", format(df$acquired, "%Y-%m-%d"), "</b> ",
      #   format(df$acquired, "%H:%M:%S")
      # )  
      
      #get min and max values across raster stack to set global color ramp
      global_min = minmax(rs)[1,] |> min()
      global_max = minmax(rs)[2,] |> max()
      
      #generate animation frames
      frames = pblapply(gs_dates, function(d){
        
        #1. Make scatterplot with vline at current acquisition date
        scat_p = ggplot(
          df
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
          # geom_vline(
          #   data = harvest_dates_df[harvest_dates_df$block_id %in% unique(subset_df$block_id), ],
          #   aes(xintercept = harvest_finish_date, linetype = 'Harvest end date')
          # ) +
          # geom_vline(
          #   data = harvest_dates_df[harvest_dates_df$block_id %in% unique(subset_df$block_id), ],
          #   aes(xintercept = harvest_start_date, linetype = 'Harvest start date')
          # ) +
          # scale_linetype_manual(
          #   name = "Harvest date",
          #   values = c("Harvest start date" = "dotted", "Harvest end date" = "dashed")
          # ) +
          labs(x = 'Acquisition date', y = unique(subset_df$var_name)) +
          theme_classic() +
          theme(
            axis.text.x.top  = element_blank(),
            axis.ticks.x.top = element_blank(),
            axis.text.x = element_text(angle = 90)
            ,legend.position = 'none' #remove legend
          ) +
          geom_vline(
            xintercept = d,
            color     = "blue",
            linewidth = 1
          )
        
        #2. plot spatraster
        
        r = rs[[str_replace(d,' UTC','')]]
        
        r_p = ggplot()+
          geom_spatraster(data = r)+
          coord_sf(expand = FALSE) + # Removes expansion
          scale_fill_viridis_c(
            limits   = c(global_min, global_max),
            na.value = "transparent",
            oob      = scales::squish   # handles any values outside limits gracefully
            ,option = color_palette
          ) +
          theme_void()+
          labs(
            # subtitle = d
            subtitle = paste0('<b>', format(d, "%Y-%m-%d"), '</b> ', format(d, "%H:%M:%S"))
            ,
            title=paste0('<b>',ds,'</b>')
          )+
          theme(
            plot.subtitle = element_markdown(margin = margin(b = 5.5)),  # default ggplot spacing
            plot.title    = element_markdown(margin = margin(b = 5.5))
            ,plot.margin = margin(0, 0, 0, 0, "cm")
          )
        
        #3. knit together with patchwork
        
        library(patchwork)
        combined_p = r_p+scat_p+plot_layout(ncol=1, widths=c(1,1))
        
        
        return(combined_p)
      })
      names(frames) = gs_dates
      return(frames)
    })
    names(ds_frames) = datasets
    
    
    #merge frames so all datasets are next to each other
    cat('Merging frames')
    frames_wide_l <- pblapply(1:length(gs_dates), function(i){
      
      imgs <- lapply(datasets, function(ds) {
        tmp <- tempfile(fileext = ".png")
        ggsave(tmp, plot = ds_frames[[ds]][[i]], width = 6, height = 8, dpi = 200)
        image_read(tmp)
      })
      
      image_append(do.call(c, imgs), stack = FALSE)  # horizontal
    })
    
    #save as mp4
    wide_dir <- tempfile("wide_frames_")
    dir.create(wide_dir)
    
    lapply(seq_along(frames_wide_l), function(i) {
      image_write(
        frames_wide_l[[i]],
        file.path(wide_dir, sprintf("frame_%04d.png", i))
      )
    })
    
    av::av_encode_video(
      sort(list.files(wide_dir, full.names = TRUE)),
      output    = combined_animation_f,
      framerate = target_fps
    )
    
    # #make gif
    # cat('Animating')
    # animation <- image_animate(
    #   do.call(c, frames_wide_l),
    #   fps = target_fps
    # )
  }
  
  print(paste(block,'done'))
}
