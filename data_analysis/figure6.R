# import packages
library(tidyverse)

### figure 6A, ELM N-degrons
# Wilcoxon rank-sum test
library(rstatix)

HEK_Nterm_ELM_N_degron_wilcoxon_test <- HEK_Nterm_ELM_N_degron |>
  wilcox_test(half_life ~ ELM_N_degron) |> 
  add_significance('p') |> 
  filter(p < 0.05)

# point range plot
library(ggpubr)

point_range_Nterm_degron <- HEK_Nterm_ELM_N_degron |> 
  ggplot() +
  geom_point(
    aes(
      x = fct_reorder(ELM_N_degron, half_life),
      y = half_life
    ),
    position = position_jitter(width = 0.3),
    color = 'black',
    alpha = 0.1,
    size = 0.5
  ) +
  stat_summary(
    aes(
      x = fct_reorder(ELM_N_degron, half_life),
      y = half_life 
    ),
    fun.data = 'mean_cl_boot', color = color_1, linewidth = 0.2, size = 0.5
  ) +
  stat_pvalue_manual(
    data = HEK_Nterm_ELM_N_degron_wilcoxon_test, label = 'p.signif', label.size = 6, 
    tip.length = 0, y.position = c(160, 180)
  ) +
  labs(x = '', y = '') +
  theme(
    axis.text.x = element_text(color = 'black', size = 9, angle = 30, hjust = 1),
    axis.text.y = element_text(color = 'black', size = 9)
  )

ggsave(
  filename = 'figures/figure5/point_range_Nterm_degron.eps',
  plot = point_range_Nterm_degron,
  device = cairo_ps,
  height = 2, width = 2.5, units = 'in',
  fallback_resolution = 1200
)

### figure 6B, ELM motifs
# import GSEA result from ELM motif analysis
Nterm_ELM_motif_GSEA_des_Kd <- read_csv(
  'data_source/ELM_degron/Nterm_ELM_motif_GSEA_des_Kd.csv'
)

# filter the enriched term which p < 0.05, rearrange the terms and export as a table
library(gridExtra)

Nterm_ELM_motif_GSEA_des_Kd_enriched_term <- Nterm_ELM_motif_GSEA_des_Kd |> 
  filter(pvalue < 0.05) |> 
  select(
    'ELM.Name' = Description, NES, pvalue
  ) |> 
  mutate(
    NES = round(NES, 2) , 
    'P.value' = format(pvalue, format = "e", digits = 2)
  ) |> 
  select(-pvalue) |> 
  arrange(desc(NES))

grid.table(Nterm_ELM_motif_GSEA_des_Kd_enriched_term)

# generate proteoform list for related motif
PIP_degron_list <- Nterm_ELM_motif_GSEA_des_Kd |> 
  filter(Description == 'CRL4-Cdt2 binding PIP degron') |> 
  select(core_enrichment) |> 
  separate_rows(core_enrichment, sep = '/') |> 
  pull(core_enrichment)

CASK_binding_list <- Nterm_ELM_motif_GSEA_des_Kd |> 
  filter(Description == 'CASK CaMK domain binding ligand motif') |> 
  select(core_enrichment) |> 
  separate_rows(core_enrichment, sep = '/') |> 
  pull(core_enrichment)

AGC_docking_list <- Nterm_ELM_motif_GSEA_des_Kd |> 
  filter(Description == 'AGC Kinase docking motif') |> 
  select(core_enrichment) |> 
  separate_rows(core_enrichment, sep = '/') |> 
  pull(core_enrichment)

HCF_binding_list <- Nterm_ELM_motif_GSEA_des_Kd |> 
  filter(Description == 'HCF-1 binding motif') |> 
  select(core_enrichment) |> 
  separate_rows(core_enrichment, sep = '/') |> 
  pull(core_enrichment)

# Nterm proteoforms ranking plot
Nterm_Kd_rank <- HEK_Nterm_Kd_half_life_sequence |> 
  arrange(desc(Kd)) |> 
  select(Index, Kd) |> 
  mutate(ranking = rank(Kd))

ranking_plot_Nterm_Kd <- ggplot() +
  # for all the proteoforms that are not related to selected enriched motif
  geom_point(
    data = Nterm_Kd_rank |> filter(
      ! Index %in% PIP_degron_list,
      ! Index %in% CASK_binding_list,
      ! Index %in% AGC_docking_list,
      ! Index %in% HCF_binding_list
    ),
    aes(
      x = ranking,
      y = Kd
    ),
    shape = 21, color = 'gray50', fill = 'transparent'
  ) +
  # Nterm proteoforms related to CRL4-Cdt2 binding PIP degron
  geom_point(
    data = Nterm_Kd_rank |> filter(Index %in% PIP_degron_list),
    aes(
      x = ranking,
      y = Kd
    ),
    shape = 21, color = 'black', fill = color_1, size = 3
  ) +
  # Nterm proteoforms related to CASK CaMK domain binding ligand motif
  geom_point(
    data = Nterm_Kd_rank |> filter(Index %in% CASK_binding_list),
    aes(
      x = ranking,
      y = Kd
    ),
    shape = 21, color = 'black', fill = color_2, size = 3
  ) +
  # Nterm proteoforms related to AGC Kinase docking motif
  geom_point(
    data = Nterm_Kd_rank |> filter(Index %in% AGC_docking_list),
    aes(
      x = ranking,
      y = Kd
    ),
    shape = 21, color = 'black', fill = color_3, size = 3
  ) +
  # Nterm proteoforms related to HCF-1 binding motif
  geom_point(
    data = Nterm_Kd_rank |> filter(Index %in% HCF_binding_list),
    aes(
      x = ranking,
      y = Kd
    ),
    shape = 21, color = 'black', fill = color_4, size = 3
  ) +
  labs(x = '', y = '') +
  theme(
    axis.text = element_text(color = 'black', family = 'arial', size = 9)
  )

ggsave(
  filename = 'figures/figure5/ranking_plot_Nterm_Kd.eps',
  plot = ranking_plot_Nterm_Kd,
  height = 1.5, width = 3, units = 'in'
)

## protein example plot
# without motif
protein_example_data <- tribble(
  ~ start, ~ end, ~ top, ~ bottom,
  0, 1000, 1, 2
)

protein_example_without_motif <- ggplot() +
  geom_rect(
    data = protein_example_data,
    aes(
      xmin = start, xmax = end, ymin = bottom, ymax = top
    ),
    fill = "gray80", color = "black"
  ) +
  theme_void() +
  theme(
    axis.text = element_blank()
  )

ggsave(
  filename = 'figures/figure5/protein_example_without_motif.eps',
  plot = protein_example_without_motif,
  height = 1, width = 4, units = 'in'
)

# with motif
protein_example_data <- tribble(
  ~ start, ~ end, ~ top, ~ bottom,
  0, 1000, 1, 2
)

motif_data <- tribble(
  ~ motif, ~ start_position, ~ end_position, 
  'motif_1', 100, 300,
  'motif_2', 600, 900
)

protein_example_with_motif <- ggplot() +
  geom_rect(
    data = protein_example_data,
    aes(
      xmin = start, xmax = end, ymin = bottom, ymax = top
    ),
    fill = "gray80", color = "black"
  ) +
  geom_rect(
    data = motif_data, 
    aes(
      xmin = start_position, 
      xmax = end_position, 
      ymin = protein_example_data$bottom, 
      ymax = protein_example_data$top, 
      fill = motif
    ), 
    show.legend = FALSE, color = 'transparent'
  ) +
  scale_fill_manual(
    values = c(
      "motif_1" = color_1,
      "motif_2" = color_2
    )
  ) +
  theme_void() +
  theme(
    axis.text = element_blank()
  )

ggsave(
  filename = 'figures/figure5/protein_example_with_motif.eps',
  plot = protein_example_with_motif,
  height = 1, width = 4, units = 'in'
)

### figure 6C, degron-related ELM motif
# Wilcoxon rank-sum test
library(rstatix)

HEK_Nterm_sequence_ELM_motif_wilcoxon_test <- HEK_Nterm_sequence_ELM_motif |> 
  filter(str_detect(matched_motifs, 'degron')) |> 
  filter(matched_motifs != 'N-degron') |> 
  wilcox_test(half_life ~ matched_motifs) |> 
  add_significance('p') |> 
  filter(p < 0.05)

# boxplot
font_add(family = 'arial', regular = 'arial.ttf')
showtext_auto()

boxplot_Nterm_ELM_motifs <- HEK_Nterm_sequence_ELM_motif |> 
  filter(
    matched_motifs %in% c(
      'APC/C_Apc2-docking motif', 
      'LEDGF/P75 IBD binding site', 
      'Trans-Golgi Network-Endosome-Lysosome-sorting signals',
      'Corepressor nuclear receptor box',
      'FBP Nbox motif',
      'Gamma-adaptin ear interaction motif',
      'KLHDC2 C-terminal GG degrons',
      'CRL4-Cdt2 binding PIP degron',
      'CASK CaMK domain binding ligand motif',
      'DCAF12 C-terminal diGlu degrons'
    )
  ) |> 
  ggplot() +
  geom_boxplot(
    aes(
      x = fct_reorder(matched_motifs, half_life),
      y = half_life
    ),
    color = color_2, outliers = FALSE
  ) +
  stat_pvalue_manual(
    data = HEK_Nterm_sequence_ELM_motif_wilcoxon_test, label = 'p.signif', label.size = 6, 
    tip.length = 0, y.position = c(23, 26, 29)
  ) +
  labs(x = '', y = '') +
  coord_cartesian(ylim = c(0, 30)) +
  theme_bw() +
  theme(
    panel.grid.major = element_line(color = 'gray', linewidth = 0.2),
    panel.grid.minor = element_line(color = 'gray', linewidth = 0.1),
    axis.text.x = element_blank(),
    axis.text.y = element_text(color = 'black', size = 9)
  )

ggsave(
  filename = 'figures/figure5/boxplot_Nterm_ELM_motifs.eps',
  plot = boxplot_Nterm_ELM_motifs,
  height = 1.5, width = 3.5, units = 'in'
)

### figure 5C
## get 7mer of N-terminal sequence
# FYLIW
Nterm_degron_FYLIW <- HEK_Nterm_ELM_N_degron |> 
  filter(ELM_N_degron == 'FYLIW') |> 
  mutate(
    Nterm_7mer = substr(Nterm_sequence, start = 1, stop = 7)
  ) |> 
  arrange(half_life)

write_csv(
  Nterm_degron_FYLIW |> select(Nterm_7mer),
  file = 'data_source/ELM_degron/Nterm_degron_FYLIW_7mer.csv'
)

# ED
Nterm_degron_ED <- HEK_Nterm_ELM_N_degron |> 
  filter(ELM_N_degron == 'ED') |> 
  mutate(
    Nterm_7mer = substr(Nterm_sequence, start = 1, stop = 7)
  ) |> 
  arrange(half_life)

write_csv(
  Nterm_degron_ED |> select(Nterm_7mer),
  file = 'data_source/ELM_degron/Nterm_degron_ED_7mer.csv'
)

# NQ
Nterm_degron_NQ <- HEK_Nterm_ELM_N_degron |> 
  filter(ELM_N_degron == 'NQ') |> 
  mutate(
    Nterm_7mer = substr(Nterm_sequence, start = 1, stop = 7)
  ) |> 
  arrange(half_life)

write_csv(
  Nterm_degron_NQ |> select(Nterm_7mer),
  file = 'data_source/ELM_degron/Nterm_degron_NQ_7mer.csv'
)

# FYLIW heatmap
Nterm_degron_FYLIW_top25 <- Nterm_degron_FYLIW |> 
  slice(1:25) |> 
  select(Index, half_life)

Nterm_degron_FYLIW_top25_matrix <- data.matrix(Nterm_degron_FYLIW_top25)
rownames(Nterm_degron_FYLIW_top25_matrix) <- Nterm_degron_FYLIW_top25$Index

col_mat <- colorRamp2(
  breaks = c(0.2, 0.6, 1.0),
  colors = c('blue', 'white', 'red')
)

Heatmap(
  matrix = Nterm_degron_FYLIW_top25_matrix[,2],
  col = col_mat,
  cluster_rows = FALSE
)

# ED heatmap
Nterm_degron_ED_top25 <- Nterm_degron_ED |> 
  slice(1:25) |> 
  select(Index, half_life)

Nterm_degron_ED_top25_matrix <- data.matrix(Nterm_degron_ED_top25)
rownames(Nterm_degron_ED_top25_matrix) <- Nterm_degron_ED_top25$Index

col_mat <- colorRamp2(
  breaks = c(0.2, 0.6, 1.0),
  colors = c('blue', 'white', 'red')
)

Heatmap(
  matrix = Nterm_degron_ED_top25_matrix[,2],
  col = col_mat,
  cluster_rows = FALSE
)

# NQ heatmap
Nterm_degron_NQ_top25 <- Nterm_degron_NQ |> 
  slice(1:25) |> 
  select(Index, half_life)

Nterm_degron_NQ_top25_matrix <- data.matrix(Nterm_degron_NQ_top25)
rownames(Nterm_degron_NQ_top25_matrix) <- Nterm_degron_NQ_top25$Index

col_mat <- colorRamp2(
  breaks = c(0.2, 0.6, 1.0),
  colors = c('blue', 'white', 'red')
)

Heatmap(
  matrix = Nterm_degron_NQ_top25_matrix[,2],
  col = col_mat,
  cluster_rows = FALSE
)

### figure 5D, cleaving proteases
# Wilcoxon rank-sum test
Nterm_topfiner_cleaving_proteases_wilcoxon_test <- Nterm_topfinder_cleaving_proteases |> 
  wilcox_test(half_life ~ Cleaving.proteases) |> 
  add_significance('p') |> 
  filter(p < 0.05)

# boxplot
font_add(family = 'arial', regular = 'arial.ttf')
showtext_auto()

boxplot_Nterm_cleaving_proteases <- Nterm_topfinder_cleaving_proteases |> 
  filter(Cleaving.proteases != 'HTRA2') |> 
  ggplot() +
  geom_boxplot(
    aes(x = fct_reorder(Cleaving.proteases, half_life), y = half_life),
    color = color_3, fill = 'transparent', outliers = FALSE
  ) +
  labs(x = '', y = '') +
  stat_pvalue_manual(
    data = Nterm_topfiner_cleaving_proteases_wilcoxon_test,
    label = 'p.signif',
    y.position = c(40, 32, 28, 36),
    tip.length = 0,
    coord.flip = TRUE,
    label.size = 6
  ) +
  coord_flip(ylim = c(0, 40)) +
  theme_bw() +
  theme(
    panel.grid.major = element_line(color = 'gray', linewidth = 0.2),
    panel.grid.minor = element_line(color = 'gray', linewidth = 0.1),
    axis.text.x = element_text(color = 'black', size = 10),
    axis.text.y = element_text(color = 'black', size = 9)
  )

ggsave(
  filename = 'figures/figure5/boxplot_Nterm_cleaving_proteases.eps',
  plot = boxplot_Nterm_cleaving_proteases,
  height = 4, width = 2.2, units = 'in'
)
