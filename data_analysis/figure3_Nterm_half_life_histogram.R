#import packages
packages_names <- c("tidyverse", "readxl", "writexl", 'rstatix', 'showtext')
lapply(packages_names, require, character.only = TRUE)

#histogram for Nterm half life
font_add(family = 'arial', regular = 'arial.ttf')
showtext_auto()

Nterm_half_life_median <- HEK_Nterm_Kd_half_life |> 
  get_summary_stats(half_life, type = 'median') |> 
  pull(median)

histogram_Nterm_half_life <- HEK_Nterm_Kd_half_life |> 
  ggplot() +
  geom_histogram(aes(x = half_life), color = 'black', fill = color_1) +
  geom_vline(xintercept = Nterm_half_life_median, linetype = 'dashed', color = 'gray') +
  annotate(
    'text', label = paste(round(Nterm_half_life_median, digits = 1), ' hr'),
    x = 23, y = 1200, size = 3, color = 'black', 
  ) +
  xlim(0, 100) +
  labs(x = 'Half life (hr)', y = 'Count') +
  theme(
    axis.title = element_text(size = 10),
    axis.text = element_text(size = 8)
  )

ggsave(
  'figures/figure3/histogram_Nterm_half_life.eps',
  plot = histogram_Nterm_half_life,
  height = 1.5, width = 2, units = 'in'
)