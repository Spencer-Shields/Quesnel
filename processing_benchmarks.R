#spectral separability stats
benchmarks = tribble(
  ~cores, ~n_completed,
  24, 20,
  14, 27,
  19, 23,
  8, 56,
  11, 36,
  4,15
)

ggplot(benchmarks, aes(x = cores, y = n_completed))+
  geom_line()+
  geom_point()
  