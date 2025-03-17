# import packages
library(tidyverse)

### figure 6A, N-terminus secondary structure, solvent accessibility and IDR analysis
## Tip 1: NetSurfP - 3.0 does not accept protein sequence more than 5,000 residues. 
## You can manually check the sequence length and split longer sequences into shorter segments.
## Tip 2: NetSurfP - 3.0 does not accept U as an amino acid residue.
## You can use the GRAVY Calculator (https://www.gravy-calculator.de/index.php) to locate U residues in your FASTA file and replace them with C.
## Tip 3: Do not use NetSurfP - 3.0.
## Wilcoxon rank-sum test
library(rstatix)

Nterm_degradation_alphafold_half_life <- read_csv(
  'data_source/Nterm_degradation_structuremap/Nterm_degradation_alphafold_half_life.csv'
)

Nterm_degradation_alphafold_half_life_property <- Nterm_degradation_alphafold_half_life |> 
  select(Index, half_life, structure_group, accessibility, IDR) |> 
  pivot_longer(cols = c(structure_group, accessibility, IDR), names_to = 'property', values_to = 'category')

Nterm_degradation_alphafold_half_life_property |> 
  group_by(category) |> 
  get_summary_stats(half_life, type = 'median')

write_csv(
  Nterm_degradation_alphafold_half_life_property, 
  'data_source/Nterm_structuremap/Nterm_degradation_alphafold_half_life_property.csv'
)

Nterm_alphafold_half_life_wilcox_test <- Nterm_alphafold_half_life_property |> 
  group_by(property) |> 
  wilcox_test(half_life ~ category) |> 
  add_significance('p') |> 
  filter(p.signif != 'ns')

# point range plot
library(ggpubr)

point_range_Nterm_alphafold_property <- Nterm_alphafold_half_life_property |> 
  ggplot() +
  geom_point(
    aes(
      x = category, y = half_life
    ),
    position = position_jitter(width = 0.3),
    size = 0.5,
    color = 'black',
    alpha = 0.3
  ) +
  stat_summary(
    aes(
      x = category, y = half_life, colour = property
    ),
    fun.data = 'mean_cl_boot', linewidth = 0.2, size = 0.5, show.legend = FALSE
  ) +
  scale_color_manual(
    values = c(
      'structure_group' = color_1,
      'accessibility' = color_2,
      'IDR' = color_3
    )
  ) +
  labs(x = '', y = '') +
  stat_pvalue_manual(
    data = Nterm_alphafold_half_life_wilcox_test, label = 'p.signif', tip.length = 0,
    size = 6, hide.ns = 'p', y.position = c(170, 160, 180, 170)
  ) +
  facet_grid(cols = vars(factor(property, levels = c(
    'structure_group', 'accessibility', 'IDR'
  ))), scales = 'free_x', space = 'free_x') +
  theme(
    axis.text.x = element_text(family = 'arial', color = 'black', angle = 30, hjust = 1, size = 8),
    axis.text.y = element_text(family = 'arial', color = 'black', size = 8),
  )

ggsave(
  filename = 'figures/figure6/point_range_Nterm_alphafold_property.eps',
  plot = point_range_Nterm_alphafold_property,
  device = cairo_ps,
  height = 2.5, width = 4, units = 'in',
  fallback_resolution = 1200
)

### figure 6B, ubiquitination sites occurrence
# Spearman correlation test
cor.test(
  Nterm_ubiquitination_site_occurrence$half_life,
  Nterm_ubiquitination_site_occurrence$n,
  method = 'spearman'
)

# spearman correlation point range plot
plot_half_life_ubi_occurence <- Nterm_ubiquitination_site_occurrence |> 
  ggplot() +
  geom_point(
    aes(
      x = category, y = n
    ),
    position = position_jitter(width = 0.3),
    color = 'black',
    alpha = 0.1
  ) +
  stat_summary(
    aes(
      x = category, y = n
    ),
    fun.data = 'mean_cl_boot', color = color_4, linewidth = 0.2, size = 0.5
  ) +
  labs(x = '', y = '') +
  theme(
    axis.text.x = element_text(family = 'arial', color = 'black', angle = 30, hjust = 1),
    axis.text.y = element_text(family = 'arial', color = 'black')
  )

ggsave(
  filename = 'figures/figure6/plot_half_life_ubi_occurence.eps',
  device = cairo_ps,
  plot = plot_half_life_ubi_occurence,
  height = 2.2, width = 1.5, units = 'in',
  fallback_resolution = 1200
)

### figure 6C, phosphorylation sites occurrence
# Spearman correlation test
cor.test(
  Nterm_phosphorylation_site_occurrence$half_life,
  Nterm_phosphorylation_site_occurrence$n,
  method = 'spearman'
)

# spearman correlation point range plot
plot_half_life_phos_occurence <- Nterm_phosphorylation_site_occurrence |> 
  ggplot() +
  geom_point(
    aes(
      x = category, y = n
    ),
    position = position_jitter(width = 0.3),
    color = 'black',
    alpha = 0.1
  ) +
  stat_summary(
    aes(
      x = category, y = n
    ),
    fun.data = 'mean_cl_boot', color = color_1, linewidth = 0.2, size = 0.5
  ) +
  labs(x = '', y = '') +
  theme(
    axis.text.x = element_text(family = 'arial', color = 'black', angle = 30, hjust = 1),
    axis.text.y = element_text(family = 'arial', color = 'black')
  )

ggsave(
  filename = 'figures/figure6/plot_half_life_phos_occurence.eps',
  device = cairo_ps,
  plot = plot_half_life_phos_occurence,
  height = 2.2, width = 1.5, units = 'in',
  fallback_resolution = 1200
)

### figure 6D, hydropathy, NCPR and kappa of Nterm 13 mer 
Nterm_13mer_sequence_features_selected <- Nterm_13mer_sequence_features |> 
  select(Index, category, hydropathy, NCPR, kappa) |> 
  pivot_longer(cols = hydropathy:kappa, names_to = 'features', values_to = 'values')

# Wilcoxon rank-sum test
library(rstatix)

Nterm_13mer_sequence_features_wilcox_test <- Nterm_13mer_sequence_features_selected |> 
  group_by(features) |> 
  wilcox_test(values ~ category) |> 
  filter(p.adj.signif != "ns")

# boxplot
library(ggpubr)

boxplot_Nterm_13mer_sequence_features <- Nterm_13mer_sequence_features_selected |> 
  ggplot() +
  geom_boxplot(
    aes(
      x = category, y = values, color = features
    ), show.legend = FALSE
  ) +
  scale_color_manual(
    values = c(
      'hydropathy' = color_2,
      'NCPR' = color_3,
      'kappa' = color_4
    )
  ) +
  labs(x = '', y = '') +
  stat_pvalue_manual(
    data = Nterm_13mer_sequence_features_wilcox_test, label = 'p.adj.signif', tip.length = 0,
    size = 6, hide.ns = 'p', y.position = c(0.9, 0.7, 0.5, 7, 6.5, 6, 1.2, 1)
  ) +
  facet_grid(rows = vars(factor(features, levels = c(
    'hydropathy', 'NCPR', 'kappa'
  ))), scales = 'free_y', space = 'free_x') +
  theme(
    axis.text.x = element_text(family = 'arial', color = 'black', angle = 30, hjust = 1, size = 8),
    axis.text.y = element_text(family = 'arial', color = 'black', size = 8),
  )

ggsave(
  filename = 'figures/figure6/boxplot_Nterm_13mer_sequence_features.eps',
  plot = boxplot_Nterm_13mer_sequence_features,
  height = 5, width = 1.8, units = 'in'
)
