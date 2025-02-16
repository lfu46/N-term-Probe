#import packages
packages_names <- c("tidyverse")
lapply(packages_names, require, character.only = TRUE)

## figure3A, workflow
# medium exchange 3h
data_0h <- tibble(
  Medium = c('Heavy', 'Light'),
  time_period = c(0/24, 24/24) 
)

medium_exchange_0h <- data_0h |> 
  ggplot() +
  geom_bar(aes(
    x = "", y = -time_period, fill = Medium
  ) , stat = "identity", width = 1, color = 'white') +
  coord_polar("y") +
  scale_fill_manual(values = c("Heavy" = color_1, "Light" = color_2)) +
  theme_void() +
  theme(
    legend.position = 'none'
  )

ggsave(
  filename = 'figures/figure3/medium_exchange_0h.eps',
  plot = medium_exchange_0h,
  height = 0.5, width = 0.5, units = 'in'
)

# medium exchange 3h
data_3h <- tibble(
  Medium = c('Heavy', 'Light'),
  time_period = c(3/24, 21/24) 
)

medium_exchange_3h <- data_3h |> 
  ggplot() +
  geom_bar(aes(
    x = "", y = -time_period, fill = Medium
  ) , stat = "identity", width = 1, color = 'white') +
  coord_polar("y") +
  scale_fill_manual(values = c("Heavy" = color_1, "Light" = color_2)) +
  theme_void() +
  theme(
    legend.position = 'none'
  )

ggsave(
  filename = 'figures/figure3/medium_exchange_3h.eps',
  plot = medium_exchange_3h,
  height = 0.5, width = 0.5, units = 'in'
)

# medium exchange 6h
data_6h <- tibble(
  Medium = c('Heavy', 'Light'),
  time_period = c(6/24, 18/24) 
)

medium_exchange_6h <- data_6h |> 
  ggplot() +
  geom_bar(aes(
    x = "", y = -time_period, fill = Medium
  ) , stat = "identity", width = 1, color = 'white') +
  coord_polar("y") +
  scale_fill_manual(values = c("Heavy" = color_1, "Light" = color_2)) +
  theme_void() +
  theme(
    legend.position = 'none'
  )

ggsave(
  filename = 'figures/figure3/medium_exchange_6h.eps',
  plot = medium_exchange_6h,
  height = 0.5, width = 0.5, units = 'in'
)

# medium exchange 9h
data_9h <- tibble(
  Medium = c('Heavy', 'Light'),
  time_period = c(9/24, 15/24) 
)

medium_exchange_9h <- data_9h |> 
  ggplot() +
  geom_bar(aes(
    x = "", y = -time_period, fill = Medium
  ) , stat = "identity", width = 1, color = 'white') +
  coord_polar("y") +
  scale_fill_manual(values = c("Heavy" = color_1, "Light" = color_2)) +
  theme_void() +
  theme(
    legend.position = 'none'
  )

ggsave(
  filename = 'figures/figure3/medium_exchange_9h.eps',
  plot = medium_exchange_9h,
  height = 0.5, width = 0.5, units = 'in'
)

# medium exchange 12h
data_12h <- tibble(
  Medium = c('Heavy', 'Light'),
  time_period = c(12/24, 12/24) 
)

medium_exchange_12h <- data_12h |> 
  ggplot() +
  geom_bar(aes(
    x = "", y = -time_period, fill = Medium
  ) , stat = "identity", width = 1, color = 'white') +
  coord_polar("y") +
  scale_fill_manual(values = c("Heavy" = color_1, "Light" = color_2)) +
  theme_void() +
  theme(
    legend.position = 'none'
  )

ggsave(
  filename = 'figures/figure3/medium_exchange_12h.eps',
  plot = medium_exchange_12h,
  height = 0.5, width = 0.5, units = 'in'
)

# medium exchange 24h
data_24h <- tibble(
  Medium = c('Heavy', 'Light'),
  time_period = c(24/24, 0/24) 
)

medium_exchange_24h <- data_24h |> 
  ggplot() +
  geom_bar(aes(
    x = "", y = -time_period, fill = Medium
  ) , stat = "identity", width = 1, color = 'white') +
  coord_polar("y") +
  scale_fill_manual(values = c("Heavy" = color_1, "Light" = color_2)) +
  theme_void() +
  theme(
    legend.position = 'none'
  )

ggsave(
  filename = 'figures/figure3/medium_exchange_24h.eps',
  plot = medium_exchange_24h,
  height = 0.5, width = 0.5, units = 'in'
)

# MS2 spectrum 
MS2_spectrum_Nt1_17_spectrum_20109 <- read_csv(
  'figures/figure3/E_LF_Nterm_Deg_HEK_Nt_1_17_01162025_MS2_20109.csv',
  skip = 7,
  col_names = TRUE,
  name_repair = 'universal'
)

MS2_spectrum <- MS2_spectrum_Nt1_17_spectrum_20109 |> 
  ggplot() +
  geom_bar(
    aes(
      x = Mass,
      y = Intensity
    ), stat = 'identity', color = 'black' 
  ) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  coord_cartesian(xlim = c(0, 1200)) +
  theme_classic() +
  theme(
    axis.text = element_blank(),
    axis.title = element_blank()
  )

ggsave(
  filename = 'figures/figure3/MS2_spectrum.eps',
  plot = MS2_spectrum,
  height = 2, width = 2, units = 'in'
)

# MS2 spectrum TMT report ion intensity
MS2_spectrum_TMT_intensity <- read_tsv(
  'data_source/raw_data/HEK_Nt_1_psm.tsv'
) |> 
  filter(Spectrum == 'E_LF_Nterm_Deg_HEK_Nt_1_17_01162025.20109.20109.3') |> 
  select(`126_0h`:`131_24h`) |> 
  pivot_longer(cols = `126_0h`:`131_24h`, names_to = 'TMT_report_ion', values_to = 'Intensity')

MS2_spectrum_TMT <- MS2_spectrum_TMT_intensity |> 
  ggplot() +
  geom_bar(
    aes(
      x = TMT_report_ion,
      y = Intensity
    ), stat = 'identity', color = 'transparent', fill = 'black', width = 0.2
  ) +
  scale_y_continuous(expand = c(0, 0)) +
  theme_classic() +
  theme(
    axis.text = element_blank(),
    axis.title = element_blank()
  )

ggsave(
  filename = 'figures/figure3/MS2_spectrum_TMT.eps',
  plot = MS2_spectrum_TMT,
  height = 2, width = 0.6, units = 'in'
)

##figure 3B, 
quantification_filter <- tibble(
  Level = c('cell doubling time normalization', 'filter criteria', 'quantification'),
  Number = c(4317, 8942-4317, 14912-8942)
)

quantification_filter$fraction <- quantification_filter$Number/sum(quantification_filter$Number)

quantification_filter$ymax <- cumsum(quantification_filter$fraction)

quantification_filter$ymin = c(0, head(quantification_filter$ymax, n=-1))

quantification_filter_plot <- quantification_filter |>
  ggplot(aes(
    x = 1, ymin = ymin, ymax = ymax, 
    fill = factor(Level, levels = rev(Level))
  )) +
  geom_rect(aes(xmin = 3, xmax = 4), color = "white") +
  coord_polar(theta = "y") +
  xlim(c(2, 4)) +
  scale_fill_manual(values = c('gray', color_2, color_1)) +
  theme_void() +
  theme(legend.position = "none")

ggsave(
  filename = 'figures/figure3/quantification_filter_plot.eps',
  plot = quantification_filter_plot,
  height = 2, width = 2, units = 'in'
)

