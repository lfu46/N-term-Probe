# import packages
library(tidyverse)

### figure 4A, workflow
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
  filename = 'figures/figure4/medium_exchange_0h.eps',
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
  filename = 'figures/figure4/medium_exchange_3h.eps',
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
  filename = 'figures/figure4/medium_exchange_6h.eps',
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
  filename = 'figures/figure4/medium_exchange_9h.eps',
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
  filename = 'figures/figure4/medium_exchange_12h.eps',
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
  filename = 'figures/figure4/medium_exchange_24h.eps',
  plot = medium_exchange_24h,
  height = 0.5, width = 0.5, units = 'in'
)

# MS2 spectrum 
MS2_spectrum_Nt1_15_spectrum_09260 <- read_csv(
  'figures/figure4/E_LF_Nterm_Deg_HEK_Nt_1_15_01162025_MS2_09260.csv',
  skip = 7,
  col_names = TRUE,
  name_repair = 'universal'
)

MS2_spectrum <- MS2_spectrum_Nt1_15_spectrum_09260 |> 
  mutate(
    rel_intensity = Intensity/max(Intensity)
  ) |> 
  ggplot() +
  geom_bar(
    aes(
      x = Mass,
      y = rel_intensity
    ), stat = 'identity', color = 'black' 
  ) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  coord_cartesian(xlim = c(0, 1700)) +
  theme_classic() +
  theme(
    axis.text = element_blank(),
    axis.title = element_blank()
  )

ggsave(
  filename = 'figures/figure4/MS2_spectrum.eps',
  plot = MS2_spectrum,
  height = 2, width = 2, units = 'in'
)

# MS2 spectrum TMT report ion intensity
MS2_spectrum_TMT_intensity <- read_tsv(
  'data_source/raw_data/HEK_Nt_1_psm.tsv'
) |> 
  filter(Spectrum == 'E_LF_Nterm_Deg_HEK_Nt_1_15_01162025.09260.09260.3') |> 
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
  filename = 'figures/figure4/MS2_spectrum_TMT.eps',
  plot = MS2_spectrum_TMT,
  height = 2, width = 0.6, units = 'in'
)

### figure 4B, quantification, curve fitting and filtering
## TMT report ion triplicates
Nterm_P07737_47_TMT_rep1 <- HEK_Nterm_1_deg_ratio |> 
  filter(Index == 'P07737_47') |> 
  ggplot() +
  geom_bar(
    aes(
      x = factor(timepoint),
      y = deg_ratio
    ), stat = 'identity', color = 'transparent', fill = 'black', width = 0.1
  ) +
  scale_y_continuous(expand = c(0, 0)) +
  theme(
    axis.text = element_blank(),
    axis.ticks = element_line(linewidth = 0),
    axis.title = element_blank()
  )

ggsave(
  filename = 'figures/figure4/Nterm_P07737_47_TMT_rep1.eps',
  plot = Nterm_P07737_47_TMT_rep1,
  height = 0.5, width = 1, units = 'in'
)

Nterm_P07737_47_TMT_rep2 <- HEK_Nterm_2_deg_ratio |> 
  filter(Index == 'P07737_47') |> 
  ggplot() +
  geom_bar(
    aes(
      x = factor(timepoint),
      y = deg_ratio
    ), stat = 'identity', color = 'transparent', fill = 'black', width = 0.1
  ) +
  scale_y_continuous(expand = c(0, 0)) +
  theme(
    axis.text = element_blank(),
    axis.ticks = element_line(linewidth = 0),
    axis.title = element_blank()
  )

ggsave(
  filename = 'figures/figure4/Nterm_P07737_47_TMT_rep2.eps',
  plot = Nterm_P07737_47_TMT_rep2,
  height = 0.5, width = 1, units = 'in'
)

Nterm_P07737_47_TMT_rep3 <- HEK_Nterm_3_deg_ratio |> 
  filter(Index == 'P07737_47') |> 
  ggplot() +
  geom_bar(
    aes(
      x = factor(timepoint),
      y = deg_ratio
    ), stat = 'identity', color = 'transparent', fill = 'black', width = 0.1
  ) +
  scale_y_continuous(expand = c(0, 0)) +
  theme(
    axis.text = element_blank(),
    axis.ticks = element_line(linewidth = 0),
    axis.title = element_blank()
  )

ggsave(
  filename = 'figures/figure4/Nterm_P07737_47_TMT_rep3.eps',
  plot = Nterm_P07737_47_TMT_rep3,
  height = 0.5, width = 1, units = 'in'
)

## non-lineaer curve fitting and linear curve fitting
# P07737_47, non-linear curve fitting
Nterm_P07737_47_deg_ratio <- bind_rows(
  HEK_Nterm_1_deg_ratio |> 
    filter(Index == 'P07737_47'),
  HEK_Nterm_2_deg_ratio |> 
    filter(Index == 'P07737_47'),
  HEK_Nterm_3_deg_ratio |> 
    filter(Index == 'P07737_47')
)

Nterm_P07737_47_nonlinear_model <- HEK_Nterm_curve_fitting_combined |> 
  filter(Index == 'P07737_47') |> 
  pivot_wider(names_from = parameters, values_from = values)

time_series <- seq(0, 24, length.out = 100)
Nterm_P07737_47_nonlinear_fitting <- (Nterm_P07737_47_nonlinear_model$A - Nterm_P07737_47_nonlinear_model$B)*exp(-Nterm_P07737_47_nonlinear_model$Kd*time_series) + Nterm_P07737_47_nonlinear_model$B

Nterm_P07737_47_nonlinear_fitting_line <- tibble(
  timepoint = time_series,
  deg_ratio = Nterm_P07737_47_nonlinear_fitting
)

font_add(family = 'arial', regular = 'arial.ttf')
showtext_auto()

Nterm_P07737_47_nonlinear_example <- ggplot() +
  geom_point(
    data = Nterm_P07737_47_deg_ratio,
    aes(
      x = timepoint,
      y = deg_ratio
    ),
    shape = 21, fill = color_1, color = 'transparent', size = 0.8
  ) +
  geom_line(
    data = Nterm_P07737_47_nonlinear_fitting_line,
    aes(
      x = timepoint,
      y = deg_ratio
    )
  ) +
  labs(x = NULL, y = NULL) +
  coord_cartesian(ylim = c(0, 1)) +
  theme(
    axis.text = element_text(size = 5),
    axis.ticks = element_line(linewidth = 0.2),
    axis.ticks.length = unit(0.03, 'in')
  )

ggsave(
  filename = 'figures/figure4/Nterm_P07737_47_nonlinear_example.eps',
  plot = Nterm_P07737_47_nonlinear_example,
  height = 0.8, width = 1.2, units = 'in'
)

# Q16204_29, linear curve fitting
Nterm_Q16204_29_deg_ratio <- bind_rows(
  HEK_Nterm_1_deg_ratio |> 
    filter(Index == 'Q16204_29'),
  HEK_Nterm_2_deg_ratio |> 
    filter(Index == 'Q16204_29'),
  HEK_Nterm_3_deg_ratio |> 
    filter(Index == 'Q16204_29')
)

Nterm_Q16204_29_linear_model <- HEK_Nterm_curve_fitting_combined |> 
  filter(Index == 'Q16204_29') |> 
  pivot_wider(names_from = parameters, values_from = values)

time_series <- seq(0, 24, length.out = 100)
Nterm_Q16204_29_linear_fitting <- exp(Nterm_Q16204_29_linear_model$lnA - Nterm_Q16204_29_linear_model$Kd * time_series)

Nterm_Q16204_29_linear_fitting_line <- tibble(
  timepoint = time_series,
  deg_ratio = Nterm_Q16204_29_linear_fitting
)

font_add(family = 'arial', regular = 'arial.ttf')
showtext_auto()

Nterm_Q16204_29_linear_example <- ggplot() +
  geom_point(
    data = Nterm_Q16204_29_deg_ratio,
    aes(
      x = timepoint,
      y = deg_ratio
    ),
    shape = 21, fill = color_2, color = 'transparent', size = 0.8
  ) +
  geom_line(
    data = Nterm_Q16204_29_linear_fitting_line,
    aes(
      x = timepoint,
      y = deg_ratio
    )
  ) +
  labs(x = NULL, y = NULL) +
  coord_cartesian(ylim = c(0, 1)) +
  theme(
    axis.text = element_text(size = 5),
    axis.ticks = element_line(linewidth = 0.2),
    axis.ticks.length = unit(0.03, 'in')
  )

ggsave(
  filename = 'figures/figure4/Nterm_Q16204_29_linear_example.eps',
  plot = Nterm_Q16204_29_linear_example,
  height = 0.8, width = 1.2, units = 'in'
)

# filtering criteria
quantification_filter <- tibble(
  Level = c('Filtering', 'quantification'),
  Number = c(6484, 14912-6484)
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
  scale_fill_manual(values = c('gray', color_3)) +
  theme_void() +
  theme(legend.position = "none")

ggsave(
  filename = 'figures/figure4/quantification_filter_plot.eps',
  plot = quantification_filter_plot,
  height = 2, width = 2, units = 'in'
)

### figure 4C, overlap of quantified proteins for Nterm and Whole Proteome
Nterm_number <- HEK_Nterm_Kd_half_life_LaminB_Tcomplex |> distinct(Index) |> nrow()
Nterm_protein_number <- HEK_Nterm_Kd_half_life_LaminB_Tcomplex |> distinct(UniProt_Accession) |> nrow()
total_protein_number <- HEK_WP_Kd_half_life_LaminB_Tcomplex |> distinct(UniProt_Accession) |> nrow()
overlap_Nterm_protein_total_protein_number <- intersect(
  HEK_Nterm_Kd_half_life_LaminB_Tcomplex |> distinct(UniProt_Accession),
  HEK_WP_Kd_half_life_LaminB_Tcomplex |> distinct(UniProt_Accession)
) |> nrow()

Nterm_protein_list <- intersect(
  HEK_Nterm_Kd_half_life_LaminB_Tcomplex |> distinct(UniProt_Accession),
  HEK_WP_Kd_half_life_LaminB_Tcomplex |> distinct(UniProt_Accession)
)

Nterm_overlap_number <- HEK_Nterm_Kd_half_life_LaminB_Tcomplex |> 
  filter(UniProt_Accession %in% Nterm_protein_list$UniProt_Accession) |> 
  nrow()

# generate data frame
Nterm_WP_overlap <- tibble(
  category = factor(
    c('Proteoform', 'Nterm Protein', 'Nterm Protein', 'total protein'),
    levels = c('total protein', 'Nterm Protein', 'Proteoform')
  ),
  count = c(Nterm_number, Nterm_protein_number - overlap_Nterm_protein_total_protein_number, overlap_Nterm_protein_total_protein_number, total_protein_number),
  group = c('group_1', 'group_2', 'group_3', 'group_4')
)

# bar plot
barplot_Nterm_WP_overlap <- Nterm_WP_overlap |> 
  ggplot() +
  geom_bar(
    aes(
      x = category,
      y = count,
      fill = group
    ),
    stat = 'identity', show.legend = FALSE
  ) +
  labs(x = '', y = 'Count') +
  scale_fill_manual(
    values = c(
      'group_1' = color_1,
      'group_2' = 'gray70',
      'group_3' = color_3,
      'group_4' = color_2
    )
  ) +
  coord_flip() +
  theme(
    axis.text = element_text(size = 8, family = c('arial')),
    axis.title = element_text(size = 8, family = c('arial'))
  )

ggsave(
  filename = 'figures/figure4/barplot_Nterm_WP_overlap.eps',
  plot = barplot_Nterm_WP_overlap, 
  height = 1.5, width = 2.5, units = 'in'
)

### figure 4D, half life distribution
library(rstatix)

# Nterm histogram
Nterm_half_life_median <- HEK_Nterm_Kd_half_life_LaminB_Tcomplex |> 
  get_summary_stats(half_life, type = 'median') |> 
  pull(median)

histogram_Nterm_half_life <- HEK_Nterm_Kd_half_life_LaminB_Tcomplex |> 
  ggplot() +
  geom_histogram(aes(x = half_life), color = 'black', fill = color_1, bins = 30) +
  geom_vline(xintercept = Nterm_half_life_median, linetype = 'dashed', color = 'gray') +
  annotate(
    'text', label = paste(round(Nterm_half_life_median, digits = 1), ' hr'),
    x = 60, y = 3900, size = 3, color = 'black', 
  ) +
  labs(x = 'Half-life (hr)', y = 'Count') +
  theme(
    axis.title = element_text(size = 8),
    axis.text = element_text(size = 8)
  )

ggsave(
  'figures/figure4/histogram_Nterm_half_life.eps',
  plot = histogram_Nterm_half_life,
  height = 1.5, width = 2, units = 'in'
)

# Whole proteome histogram
library(rstatix)

Whole_Proteome_half_life_median <- HEK_WP_Kd_half_life_LaminB_Tcomplex |> 
  filter(
    UniProt_Accession %in% HEK_Nterm_Kd_half_life_LaminB_Tcomplex$UniProt_Accession
  ) |> 
  get_summary_stats(half_life, type = 'median') |> 
  pull(median)

histogram_WP_half_life <- HEK_WP_Kd_half_life_LaminB_Tcomplex |> 
  filter(
    UniProt_Accession %in% HEK_Nterm_Kd_half_life_LaminB_Tcomplex$UniProt_Accession
  ) |> 
  ggplot() +
  geom_histogram(aes(x = half_life), color = 'black', fill = color_2, bins = 30) +
  geom_vline(xintercept = Whole_Proteome_half_life_median, linetype = 'dashed', color = 'gray') +
  annotate(
    'text', label = paste(round(Whole_Proteome_half_life_median, digits = 1), ' hr'),
    x = 65, y = 760, size = 3, color = 'black', 
  ) +
  labs(x = 'Half-life (hr)', y = 'Count') +
  theme(
    axis.title = element_text(size = 8),
    axis.text = element_text(size = 8)
  )

ggsave(
  'figures/figure4/histogram_WP_half_life.eps',
  plot = histogram_WP_half_life,
  height = 1.5, width = 2, units = 'in'
)

# Kolmogorov-Smirnov (KS) test for the half-life of N-terminal proteoforms and total proteins
Nterm_WP_half_life_ks_test <- ks.test(
  HEK_WP_Kd_half_life_LaminB_Tcomplex$half_life,
  HEK_Nterm_Kd_half_life_LaminB_Tcomplex$half_life
)

