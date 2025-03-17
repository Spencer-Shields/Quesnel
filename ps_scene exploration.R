library(tidyverse)
library(terra)
library(devtools)
source('https://raw.githubusercontent.com/Spencer-Shields/planetscope_radiometric_correction/refs/heads/main/check_radiometric_consistency.R')

base_dir = paste0(getwd(),'/data/planet_scenes')
raw_dir = paste0(base_dir,'/raw/a2c2dc75-8eab-428f-894e-6ff0baaea934/PSScene')

raw_files = list.files(raw_dir, full.names = T)
raw_rasts = raw_files[str_detect(raw_files, 'tif$')]
raw_scenes = raw_rasts[!str_detect(raw_rasts, '_udm2_')]

scene_l = lapply(raw_scenes, rast)
x11()
plot(scene_l[[1]])

x11()
plot(scene_l[[length(scene_l)]])

rmse(scene_l[[1]], scene_l[[2]])
rmse(scene_l[[1]], scene_l[[length(scene_l)]])

