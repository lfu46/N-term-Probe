# import packages
library(tidyverse)

### figure S2A, modification-related ELM motif
# Wilcoxon rank-sum test
library(rstatix)

HEK_Nterm_sequence_ELM_motif_modification_wilcoxon_test <- HEK_Nterm_sequence_ELM_motif |> 
  filter(
    matched_motifs %in% c(
      'MOD_NMyristoyl',
      'MOD_OFUCOSY',
      'MOD_Cter_Amidation',
      'MOD_N-GLC_1',
      'MOD_N-GLC_2',
      'MOD_LOK_YxT_1',
      'MOD_NEK2_1',
      'MOD_CMANNOS',
      'MOD_TYR_CSK',
      'MOD_SUMO_for_1'
    )
  ) |> 
  wilcox_test(half_life ~ matched_motifs) |> 
  add_significance('p') |> 
  filter(p < 0.05)

# point range plot
library(ggpubr)

point_range_Nterm_ELM_motifs_modification <- HEK_Nterm_sequence_ELM_motif |> 
  filter(
    matched_motifs %in% c(
      'MOD_NMyristoyl',
      'MOD_OFUCOSY',
      'MOD_Cter_Amidation',
      'MOD_N-GLC_1',
      'MOD_N-GLC_2',
      'MOD_LOK_YxT_1',
      'MOD_NEK2_1',
      'MOD_CMANNOS',
      'MOD_TYR_CSK',
      'MOD_SUMO_for_1'
    )
  ) |>
  ggplot() +
  stat_summary(
    aes(
      x = matched_motifs,
      y = half_life
    ),
    fun.data = 'mean_cl_boot', 
    color = color_2, 
    linewidth = 0.2, 
    size = 0.5
  ) +
  # stat_pvalue_manual(
  #   data = HEK_Nterm_sequence_ELM_motif_modification_wilcoxon_test |> slice(), 
  #   label = 'p.signif', label.size = 6, 
  #   tip.length = 0, y.position = c(), coord.flip = TRUE
  # ) +
  labs(x = '', y = '') +
  coord_flip(ylim = c()) +
  theme(
    axis.text.x = element_text(color = 'black', size = 8, family = 'arial'),
    axis.text.y = element_text(color = 'black', size = 8, family = 'arial')
  )

ggsave(
  filename = 'figures/figureS2/point_range_Nterm_ELM_motifs_modification.eps',
  plot = point_range_Nterm_ELM_motifs_modification,
  height = 3, width = 2.6, units = 'in'
)

### figure S2B, N-degron motif analysis
### PSSM Search (http://slim.ucd.ie/pssmsearch/, scoring method: Log Relative Binomial)
## get 7mer of N-terminal sequence
# FYLIW
Nterm_degron_FYLIW <- HEK_Nterm_ELM_N_degron |> 
  filter(ELM_N_degron == 'FYLIW') |> 
  mutate(
    Nterm_fasta_name = paste0('>', Index, '_FYLIW'),
    Nterm_7mer = substr(Nterm_sequence, start = 1, stop = 7)
  ) |> 
  arrange(half_life) |> 
  select(Nterm_fasta_name, Nterm_7mer) |> 
  mutate(
    Nterm_fasta_format = paste(Nterm_fasta_name, Nterm_7mer, sep = '\n')
  )

writeLines(
  Nterm_degron_FYLIW$Nterm_fasta_format,
  'data_source/ELM_degron/Nterm_degron_FYLIW_7mer.fasta'
)

# ED
Nterm_degron_ED <- HEK_Nterm_ELM_N_degron |> 
  filter(ELM_N_degron == 'ED') |> 
  mutate(
    Nterm_fasta_name = paste0('>', Index, '_ED'),
    Nterm_7mer = substr(Nterm_sequence, start = 1, stop = 7)
  ) |> 
  arrange(half_life) |> 
  select(Nterm_fasta_name, Nterm_7mer) |> 
  mutate(
    Nterm_fasta_format = paste(Nterm_fasta_name, Nterm_7mer, sep = '\n')
  )

writeLines(
  Nterm_degron_ED$Nterm_fasta_format,
  'data_source/ELM_degron/Nterm_degron_ED_7mer.fasta'
)

# RK
Nterm_degron_RK <- HEK_Nterm_ELM_N_degron |> 
  filter(ELM_N_degron == 'RK') |> 
  mutate(
    Nterm_fasta_name = paste0('>', Index, '_RK'),
    Nterm_7mer = substr(Nterm_sequence, start = 1, stop = 7)
  ) |> 
  arrange(half_life) |> 
  mutate(
    Nterm_fasta_format = paste(Nterm_fasta_name, Nterm_7mer, sep = '\n')
  )

writeLines(
  Nterm_degron_RK$Nterm_fasta_format,
  'data_source/ELM_degron/Nterm_degron_RK_7mer.fasta'
)

# NQ
Nterm_degron_NQ <- HEK_Nterm_ELM_N_degron |> 
  filter(ELM_N_degron == 'NQ') |> 
  mutate(
    Nterm_fasta_name = paste0('>', Index, '_NQ'),
    Nterm_7mer = substr(Nterm_sequence, start = 1, stop = 7)
  ) |> 
  arrange(half_life) |> 
  mutate(
    Nterm_fasta_format = paste(Nterm_fasta_name, Nterm_7mer, sep = '\n')
  )

writeLines(
  Nterm_degron_NQ$Nterm_fasta_format,
  'data_source/ELM_degron/Nterm_degron_NQ_7mer.fasta'
)

# Others
Nterm_degron_Others <- HEK_Nterm_ELM_N_degron |> 
  filter(ELM_N_degron == 'Others') |> 
  mutate(
    Nterm_fasta_name = paste0('>', Index, '_Others'),
    Nterm_7mer = substr(Nterm_sequence, start = 1, stop = 7)
  ) |> 
  arrange(half_life) |> 
  mutate(
    Nterm_fasta_format = paste(Nterm_fasta_name, Nterm_7mer, sep = '\n')
  )

writeLines(
  Nterm_degron_Others$Nterm_fasta_format,
  'data_source/ELM_degron/Nterm_degron_Others_7mer.fasta'
)

# C
Nterm_degron_C <- HEK_Nterm_ELM_N_degron |> 
  filter(ELM_N_degron == 'C') |> 
  mutate(
    Nterm_fasta_name = paste0('>', Index, '_C'),
    Nterm_7mer = substr(Nterm_sequence, start = 1, stop = 7)
  ) |> 
  arrange(half_life) |> 
  mutate(
    Nterm_fasta_format = paste(Nterm_fasta_name, Nterm_7mer, sep = '\n')
  )

writeLines(
  Nterm_degron_C$Nterm_fasta_format,
  'data_source/ELM_degron/Nterm_degron_C_7mer.fasta'
)

### figure S2C, cleaving proteases
# Wilcoxon rank-sum test
library(rstatix)

Cleaving.Proteases.List <- Nterm_topfinder_cleaving_proteases |> 
  group_by(Cleaving.proteases) |> 
  get_summary_stats(half_life, type = 'median') |> 
  filter(n > 1) |> 
  pull(Cleaving.proteases)

Nterm_topfiner_cleaving_proteases_wilcoxon_test <- Nterm_topfinder_cleaving_proteases |> 
  filter(Cleaving.proteases %in% Cleaving.Proteases.List) |> 
  wilcox_test(half_life ~ Cleaving.proteases) |> 
  add_significance('p') |> 
  filter(p < 0.05)

# boxplot
library(ggpubr)

boxplot_Nterm_cleaving_proteases <- Nterm_topfinder_cleaving_proteases |> 
  ggplot() +
  stat_summary(
    aes(
      x = Cleaving.proteases, 
      y = half_life
    ),
    fun.data = 'mean_cl_boot', color = color_3, linewidth = 0.2, size = 0.5
  ) +
  labs(x = '', y = '') +
  stat_pvalue_manual(
    data = Nterm_topfiner_cleaving_proteases_wilcoxon_test,
    label = 'p.signif',
    y.position = c(160, 180, 170, 190),
    tip.length = 0,
    label.size = 6
  ) +
  theme(
    axis.text.x = element_text(color = 'black', size = 8, angle = 30, hjust = 1),
    axis.text.y = element_text(color = 'black', size = 8)
  )

ggsave(
  filename = 'figures/figureS2/boxplot_Nterm_cleaving_proteases.eps',
  plot = boxplot_Nterm_cleaving_proteases,
  height = 2, width = 6, units = 'in'
)
