# import packages
packages_names <- c('tidyverse', 'rstatix', 'ggpubr', 'showtext', 'ComplexHeatmap', 'circlize')
lapply(packages_names, require, character.only = TRUE)

### figure 5A, ELM N-degrons
# Wilcoxon rank-sum test
HEK_Nterm_ELM_N_degron_wilcoxon_test <- HEK_Nterm_ELM_N_degron |>
  wilcox_test(half_life ~ ELM_N_degron) |> 
  add_significance('p') |> 
  filter(p < 0.01)

# boxplot
font_add(family = 'arial', regular = 'arial.ttf')
showtext_auto()

boxplot_Nterm_degron <- HEK_Nterm_ELM_N_degron |> 
  ggplot() +
  geom_boxplot(
    aes(
      x = fct_reorder(ELM_N_degron, half_life),
      y = half_life
    ),
    color = color_1, outliers = FALSE
  ) +
  stat_pvalue_manual(
    data = HEK_Nterm_ELM_N_degron_wilcoxon_test, label = 'p.signif', label.size = 6, 
    tip.length = 0, y.position = c(23, 26, 29)
  ) +
  labs(x = '', y = '') +
  coord_cartesian(ylim = c(0, 30)) +
  theme_bw() +
  theme(
    panel.grid.major = element_line(color = 'gray', linewidth = 0.2),
    panel.grid.minor = element_line(color = 'gray', linewidth = 0.1),
    axis.text.x = element_text(color = 'black', size = 9, angle = 30, hjust = 1),
    axis.text.y = element_text(color = 'black', size = 9)
  )

ggsave(
  filename = 'figures/figure5/boxplot_Nterm_degron.eps',
  plot = boxplot_Nterm_degron,
  height = 2, width = 2.5, units = 'in'
)

### figure 5B, ELM motifs
# Wilcoxon rank-sum test
HEK_Nterm_sequence_ELM_motif_wilcoxon_test <- HEK_Nterm_sequence_ELM_motif |> 
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
  wilcox_test(half_life ~ matched_motifs) |> 
  add_significance('p') |> 
  filter(p < 0.05) |> 
  slice(7, 13, 16)

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
