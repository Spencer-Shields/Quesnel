#since the python script dumped all the basemaps in the code folder for some reason

library(tidyverse)

from_dir = 'D:/Quesnel/code'
dir.exists(to_dir)

to_dir = 'D:/Quesnel/data/planet_basemaps/monthly_normalized_sr/raw'
dir.exists(to_dir)

from_l = list.files(from_dir, recursive = T, full.names = T)
from_l_base = from_l[str_detect(from_l, 'global_monthly_20')]

lapply(from_l_base, function(x){
  new_path = str_replace(x, from_dir, to_dir)
  new_path = sub("(.*)/", "\\1_", new_path)
  file.rename(x, new_path)
})

#check if it worked
moved_l = list.files(to_dir, full.names = T)
moved_l
rast(moved_l[1])
plot(rast(moved_l[2]))
