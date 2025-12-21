source('HLS_beast.R')

#----reload HLS beast data----

hls_beast_fl = list.files(beast_data_dir, recursive = T, full.names = T)

hls_beast_l = pblapply(hls_beast_fl,read_rds)
names(hls_beast_l) = basename(file_path_sans_ext(hls_beast_fl))

a = hls_beast_l[[2]]$BEAST_object

plot(a$trend$ncp)


f = "data/BEAST/ps_scenes_harmonized_coreg_2013-2024/NDVI/Z_20250910/12L_C5.rds"
beast_out = read_rds(f)
input_rast = beast_out$input_raster
b = beast_out$BEAST_object

#number change points
b_ncp=b$trend$ncp
b_ncp = t(b_ncp)
r_ncp = rast(b_ncp)

#changepoint locations
b_cp = b$trend$cp
b_cp = aperm(b_cp, c(2,1,3))
r_cp = rast(b_cp)

#changepoint probabilities
b_prcp = b$trend$cpOccPr
b_prcp = aperm(b_prcp, c(2,1,3))
r_prcp = rast(b_prcp)

r_totprob = app(r_prcp,'sum')
