source('scene_global_stats_filtering_20251106.R')

smdf = results_meta_df |>
  filter(dataset=='SM'
         , var_name=='SR'
         , acquired<harvest_start_date
         , jmd_scaled==1
         , block_id == '12N_T3'
         )
fsm = smdf$file_path[1]
rsm = rast(fsm,lyrs='SR')
x11()
hist(rsm)

fsm_id = smdf$id[smdf$file_path==fsm] |> unique()

fnn = str_replace(fsm, 'SM', 'Non-normalized')
rnn = rast(fnn,lyrs='SR')

rnn_sm = softmax(rnn)

r = c(rnn,rnn_sm)
names(r) = c('NN', 'SM')
plotx(r)

datatype(rnn_sm) <- "FLT8S"
