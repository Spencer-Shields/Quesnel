
#example data
f1 = "D:/Quesnel/data/planet_scenes/SM_REDO/12L_D345/20220722_185151_49_2479.tif"
f2 = "D:/Quesnel/data/planet_scenes/Non-normalized_REDO/12L_D345/20220722_185151_49_2479.tif"
raw_list = list.files(raw_dir, pattern='\\.tif',recursive = T, full.names = T)
f3 = raw_list[str_detect(raw_list, '20220722_185151_49_2479') & !str_detect(raw_list, 'udm')]
r = rast(f)
r2 = rast(f2)

raw_r = rast(f3)

#scale raster, calculate indices

# scale_factor = 65535
scale_factor = 10000
raw_scaled = raw_r / scale_factor
raw_inds = spectralIndices(raw_scaled, blue = 2, green = 4, red = 6, nir=8)

r = raw_scaled

si = spectralIndices(img = r, blue = 'blue', green = 'green', red = 'red', nir = 'nir', redEdge1 = 'rededge'
                     ,scaleFactor = scale_factor, skipRefCheck = T)
si = subset(si, c('EVI2','TVI'), negate=T)
vi = rast.batch.functions(r, include_input = F, red = 6, rededge = 7, nir = 8, yellow = 5, green = 4, blue = 2
                          , fl = c(
  hue.vi,
  gcc.vi,
  ndgr.vi,
  bi.vi,
  ci.vi,
  cri550.vi,
  gli.vi,
  tvi.vi,
  varig.vi,
  yndvi.vi,
  evi2.vi
))
plot(vi)

nn_r = c(r, si, vi)

histeraster = function(x, title = NULL){
  v = values(x, dataframe=T)
  v_long = pivot_longer(data = v, cols = names(x))
  
  x11()
  ggplot(v_long |> na.omit())+ #that appeared to work, every band is displaying correctly
    geom_histogram(aes(value)) +
    facet_wrap(vars(name), scales = 'free') +
    ggtitle(title)
}

histeraster(nn_r, title = 'Scaled non-normalized')
# corrplot::corrplot(cor(v |> na.omit()), order = 'hclust')

#calculate softmax

sm_r = softmax(nn_r)
histeraster(x = sm_r, 'Softmax')

plot(sm_r) #everything is fine now except for SR

v = values(sm_r, dataframe=T) |> na.omit()
n_vals = sapply(names(v), function(n)length(unique(v[[n]])))
t = tibble(n_vals, names = names(v))

# ggplot()+
#   geom_spatraster(data = sm_r)+
#   facet_wrap(~lyr, scales='free')

#test different scale factors to use with softmax SR
test_scale_facts = seq(from=1, to=100, by=1)

scaled_sr = pblapply(test_scale_facts, function(n){
  nn_r[['SR']]/n
})
scaled_sr = rast(scaled_sr)
names(scaled_sr) = test_scale_facts
sm_sr = softmax(scaled_sr)

v = values(sm_sr, dataframe=T) |> na.omit()

n_uniques = sapply(names(v), function(n)length(unique(v[[n]])))

ggplot(tibble(scale_facts = test_scale_facts, n_unique = n_uniques) |> filter(scale_facts<25))+
  geom_line(aes(x = scale_facts, y = n_unique))
  

#calculate z-score
# z_r = z_rast(raw_scaled)

zz_r = z_rast(nn_r)
histeraster(zz_r)

ndvi = (zz_r[['nir']] - zz_r[['red']])/(zz_r[['nir']] + zz_r[['red']])


plotx(zz_r[['TVI']])

#calculate robust z_score

zr_r = z_rast(nn_r, robust=T)
histeraster(zr_r)
plot(zr_r)

 