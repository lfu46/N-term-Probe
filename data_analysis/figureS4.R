# import packages
library(tidyverse)

### figure S4A, N-degron motif analysis
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

# GASTC
Nterm_degron_GASTC <- HEK_Nterm_ELM_N_degron |> 
  filter(ELM_N_degron == 'GASTC') |> 
  mutate(
    Nterm_fasta_name = paste0('>', Index, '_GASTC'),
    Nterm_7mer = substr(Nterm_sequence, start = 1, stop = 7)
  ) |> 
  arrange(half_life) |> 
  mutate(
    Nterm_fasta_format = paste(Nterm_fasta_name, Nterm_7mer, sep = '\n')
  )

writeLines(
  Nterm_degron_GASTC$Nterm_fasta_format,
  'data_source/ELM_degron/Nterm_degron_GASTC_7mer.fasta'
)

### figure S4B, GO and KEGG enrichment analysis for FYLIW/N-degron and GASTC/N-degron
library(tidyverse)
library(clusterProfiler)
library(org.Hs.eg.db)

## Nterm total protein
Nterm_protein <- HEK_Nterm_ELM_N_degron |> 
  distinct(UniProt_Accession) |> 
  pull()

## FYLIW/N-degron
FYLIW_N_degron_protein <- HEK_Nterm_ELM_N_degron |> 
  filter(ELM_N_degron == 'FYLIW') |> 
  distinct(UniProt_Accession) |> 
  pull()

# GO analysis
FYLIW_N_degron_GO <- enrichGO(
  gene = FYLIW_N_degron_protein,
  OrgDb = org.Hs.eg.db,
  universe = Nterm_protein,
  keyType = 'UNIPROT',
  ont = 'ALL',
  pvalueCutoff = 1,
  qvalueCutoff = 1
)

write.csv(
  FYLIW_N_degron_GO@result,
  'data_source/ELM_degron/FYLIW_N_degron_GO.csv',
)

# KEGG analysis
FYLIW_N_degron_KEGG <- enrichKEGG(
  gene = FYLIW_N_degron_protein,
  organism = 'hsa',
  keyType = 'uniprot',
  universe = Nterm_protein,
  pvalueCutoff = 1,
  qvalueCutoff = 1
)

write.csv(
  FYLIW_N_degron_KEGG@result,
  'data_source/ELM_degron/FYLIW_N_degron_KEGG.csv',
)

## GASTC/N-degron
GASTC_N_degron_protein <- HEK_Nterm_ELM_N_degron |> 
  filter(ELM_N_degron == 'GASTC') |> 
  distinct(UniProt_Accession) |> 
  pull()

# GO analysis
GASTC_N_degron_GO <- enrichGO(
  gene = GASTC_N_degron_protein,
  OrgDb = org.Hs.eg.db,
  universe = Nterm_protein,
  keyType = 'UNIPROT',
  ont = 'ALL',
  pvalueCutoff = 1,
  qvalueCutoff = 1
)

write_csv(
  GASTC_N_degron_GO@result,
  'data_source/ELM_degron/GASTC_N_degron_GO.csv',
)

# KEGG analysis
GASTC_N_degron_KEGG <- enrichKEGG(
  gene = GASTC_N_degron_protein,
  organism = 'hsa',
  keyType = 'uniprot',
  universe = Nterm_protein,
  pvalueCutoff = 1,
  qvalueCutoff = 1
)

write_csv(
  GASTC_N_degron_KEGG@result,
  'data_source/ELM_degron/GASTC_N_degron_KEGG.csv',
)

## bar plot
# FWLIW/N-degron
FYLIW_N_degron_GO <- read_csv(
  'data_source/ELM_degron/FYLIW_N_degron_GO.csv'
)

FYLIW_N_degron_GO_KEGG <- bind_rows(
  FYLIW_N_degron_GO |> 
    filter(Description %in% c(
      'Vesicle lumen',
      'Secretory granule',
      'ATP-dependent protein folding chaperone',
      'Ubiquitin protein ligase binding',
      'Ribonucleoprotein complex'
    ))
)

barplot_FYLIW_N_degron_GO_KEGG <- FYLIW_N_degron_GO_KEGG |> 
  ggplot() +
  geom_bar(
    aes(
      x = fct_reorder(Description, -log10(p.adjust)), 
      y = -log10(p.adjust)
    ),
    fill = color_1, color = 'transparent', stat = 'identity'
  ) +
  labs(x = 'GO & KEGG', y = '-log10(adjust P value)') +
  coord_flip() +
  theme(
    axis.text = element_text(size = 8, color = 'black', family = 'arial'),
    axis.title = element_text(size = 8, color = 'black', family = 'arial')
  )

ggsave(
  filename = 'figures/figureS4/barplot_FYLIW_N_degron_GO_KEGG.eps',
  plot = barplot_FYLIW_N_degron_GO_KEGG,
  height = 1.5, width = 4, units = 'in'
)

# GASTC/N-degron
GASTC_N_degron_GO <- read_csv(
  'data_source/ELM_degron/GASTC_N_degron_GO.csv'
)

GASTC_N_degron_GO_KEGG <- bind_rows(
  GASTC_N_degron_GO |> 
    filter(Description %in% c(
      'Cell adhesion molecule binding',
      'mRNA binding',
      'Chromatin binding',
      'mRNA metabolic process',
      'Regulation of apoptotic process'
    ))
)

barplot_GASTC_N_degron_GO_KEGG <- GASTC_N_degron_GO_KEGG |> 
  ggplot() +
  geom_bar(
    aes(
      x = fct_reorder(Description, -log10(p.adjust)), 
      y = -log10(p.adjust)
    ),
    fill = color_2, color = 'transparent', stat = 'identity'
  ) +
  labs(x = 'GO & KEGG', y = '-log10(adjust P value)') +
  coord_flip() +
  theme(
    axis.text = element_text(size = 8, color = 'black', family = 'arial'),
    axis.title = element_text(size = 8, color = 'black', family = 'arial')
  )

ggsave(
  filename = 'figures/figureS4/barplot_GASTC_N_degron_GO_KEGG.eps',
  plot = barplot_GASTC_N_degron_GO_KEGG,
  height = 1.5, width = 4, units = 'in'
)

# ### figure S4A, modification-related ELM motif
# # Wilcoxon rank-sum test
# library(rstatix)
# 
# HEK_Nterm_sequence_ELM_motif_modification_wilcoxon_test <- HEK_Nterm_sequence_ELM_motif |> 
#   filter(
#     matched_motifs %in% c(
#       'MOD_NMyristoyl',
#       'MOD_OFUCOSY',
#       'MOD_Cter_Amidation',
#       'MOD_N-GLC_1',
#       'MOD_N-GLC_2',
#       'MOD_LOK_YxT_1',
#       'MOD_NEK2_1',
#       'MOD_CMANNOS',
#       'MOD_TYR_CSK',
#       'MOD_SUMO_for_1'
#     )
#   ) |> 
#   wilcox_test(half_life ~ matched_motifs) |> 
#   add_significance('p') |> 
#   filter(p < 0.05)
# 
# # point range plot
# library(ggpubr)
# 
# point_range_Nterm_ELM_motifs_modification <- HEK_Nterm_sequence_ELM_motif |> 
#   filter(
#     matched_motifs %in% c(
#       'MOD_NMyristoyl',
#       'MOD_OFUCOSY',
#       'MOD_Cter_Amidation',
#       'MOD_N-GLC_1',
#       'MOD_N-GLC_2',
#       'MOD_LOK_YxT_1',
#       'MOD_NEK2_1',
#       'MOD_CMANNOS',
#       'MOD_TYR_CSK',
#       'MOD_SUMO_for_1'
#     )
#   ) |>
#   ggplot() +
#   stat_summary(
#     aes(
#       x = matched_motifs,
#       y = half_life
#     ),
#     fun.data = 'mean_cl_boot', 
#     color = color_2, 
#     linewidth = 0.2, 
#     size = 0.5
#   ) +
#   # stat_pvalue_manual(
#   #   data = HEK_Nterm_sequence_ELM_motif_modification_wilcoxon_test |> slice(), 
#   #   label = 'p.signif', label.size = 6, 
#   #   tip.length = 0, y.position = c(), coord.flip = TRUE
#   # ) +
#   labs(x = '', y = '') +
#   coord_flip(ylim = c()) +
#   theme(
#     axis.text.x = element_text(color = 'black', size = 8, family = 'arial'),
#     axis.text.y = element_text(color = 'black', size = 8, family = 'arial')
#   )
# 
# ggsave(
#   filename = 'figures/figureS2/point_range_Nterm_ELM_motifs_modification.eps',
#   plot = point_range_Nterm_ELM_motifs_modification,
#   height = 3, width = 2.6, units = 'in'
# )
# 
# ### figure S4C,
# 
# 
# ### figure 7D, docking-related ELM motif
# # Wilcoxon rank-sum test
# library(rstatix)
# 
# HEK_Nterm_sequence_ELM_motif_docking_wilcoxon_test <- HEK_Nterm_sequence_ELM_motif |> 
#   filter(str_detect(matched_motifs, 'DOC')) |> 
#   wilcox_test(half_life ~ matched_motifs) |> 
#   add_significance('p') |> 
#   filter(p < 0.05)
# 
# # point range plot
# library(ggpubr)
# 
# point_range_Nterm_ELM_motifs_docking <- HEK_Nterm_sequence_ELM_motif |> 
#   filter(
#     matched_motifs %in% c(
#       'DOC_MIT_MIM_1',
#       'DOC_MAPK_gen_1',
#       'DOC_PP2B_LxvP_1',
#       'DOC_USP7_MATH_1',
#       'DOC_USP7_MATH_2',
#       'DOC_PP1_RVXF_1',
#       'DOC_CKS1_1',
#       'DOC_CYCLIN_RxL_1',
#       'DOC_PP4_MxPP_1',
#       'DOC_RSK_DDVF_1'
#     )
#   ) |> 
#   ggplot() +
#   stat_summary(
#     aes(
#       x = matched_motifs,
#       y = half_life
#     ),
#     fun.data = 'mean_cl_boot', 
#     color = color_3, 
#     linewidth = 0.2, 
#     size = 0.5
#   ) +
#   stat_pvalue_manual(
#     data = HEK_Nterm_sequence_ELM_motif_docking_wilcoxon_test |> slice(2, 4, 6), 
#     label = 'p.signif', label.size = 6, 
#     tip.length = 0, y.position = c(100, 110, 120), coord.flip = TRUE
#   ) +
#   labs(x = '', y = '') +
#   coord_flip(ylim = c(0, 130)) +
#   theme(
#     axis.text.x = element_text(color = 'black', size = 8, family = 'arial'),
#     axis.text.y = element_text(color = 'black', size = 8, family = 'arial')
#   )
# 
# ggsave(
#   filename = 'figures/figure7/point_range_Nterm_ELM_motifs_docking.eps',
#   plot = point_range_Nterm_ELM_motifs_docking,
#   height = 3, width = 2.5, units = 'in'
# )
# 
# ### figure 7E, targeting-related ELM motif
# # Wilcoxon rank-sum test
# library(rstatix)
# 
# HEK_Nterm_sequence_ELM_motif_targeting_wilcoxon_test <- HEK_Nterm_sequence_ELM_motif |> 
#   filter(str_detect(matched_motifs, 'TRG')) |> 
#   wilcox_test(half_life ~ matched_motifs) |> 
#   add_significance('p') |> 
#   filter(p < 0.05)
# 
# # point range plot
# library(ggpubr)
# 
# point_range_Nterm_ELM_motifs_targeting <- HEK_Nterm_sequence_ELM_motif |> 
#   filter(
#     matched_motifs %in% c(
#       'TRG_ER_diArg_1',
#       'TRG_NLS_Bipartite_1',
#       'TRG_PTS1',
#       'TRG_ENDOCYTIC_2',
#       'TRG_ER_diLys_1',
#       'TRG_ER_FFAT_1',
#       'TRG_ER_KDEL_1',
#       'TRG_NES_CRM1_1',
#       'TRG_NLS_MonoCore_2',
#       'TRG_NLS_MonoExtC_3',
#       'TRG_NLS_MonoExtN_4'
#     )
#   ) |>
#   ggplot() +
#   stat_summary(
#     aes(
#       x = matched_motifs,
#       y = half_life
#     ),
#     fun.data = 'mean_cl_boot', 
#     color = color_4, 
#     linewidth = 0.2, 
#     size = 0.5
#   ) +
#   stat_pvalue_manual(
#     data = HEK_Nterm_sequence_ELM_motif_targeting_wilcoxon_test |> slice(12, 19, 20, 21, 22), 
#     label = 'p.signif', label.size = 6, 
#     tip.length = 0, y.position = c(160, 120, 130, 140, 150), coord.flip = TRUE
#   ) +
#   labs(x = '', y = '') +
#   coord_flip(ylim = c(0, 170)) +
#   theme(
#     axis.text.x = element_text(color = 'black', size = 8, family = 'arial'),
#     axis.text.y = element_text(color = 'black', size = 8, family = 'arial')
#   )
# 
# ggsave(
#   filename = 'figures/figure7/point_range_Nterm_ELM_motifs_targeting.eps',
#   plot = point_range_Nterm_ELM_motifs_targeting,
#   height = 3.2, width = 2.6, units = 'in'
# )
# 
# ### figure 7F, binding-related ELM motif
# # Wilcoxon rank-sum test
# library(rstatix)
# 
# HEK_Nterm_sequence_ELM_motif_binding_wilcoxon_test <- HEK_Nterm_sequence_ELM_motif |> 
#   filter(
#     matched_motifs %in% c(
#       'LIG_APCC_Cbox_2',
#       'LIG_AP_GAE_1',
#       'LIG_G3BP_FGDF_1',
#       'LIG_PTB_Phospho_1',
#       'LIG_NBox_RRM_1',
#       'LIG_RRM_PRI_1',
#       'LIG_14-3-3_CterR_2',
#       'LIG_eIF4E_2',
#       'LIG_UBA3_1',
#       'LIG_CNOT1_NIM_1'
#     )
#   ) |> 
#   wilcox_test(half_life ~ matched_motifs) |> 
#   add_significance('p') |> 
#   filter(p < 0.05)
# 
# # point range plot
# library(ggpubr)
# 
# point_range_Nterm_ELM_motifs_binding <- HEK_Nterm_sequence_ELM_motif |> 
#   filter(
#     matched_motifs %in% c(
#       'LIG_APCC_Cbox_2',
#       'LIG_AP_GAE_1',
#       'LIG_G3BP_FGDF_1',
#       'LIG_PTB_Phospho_1',
#       'LIG_NBox_RRM_1',
#       'LIG_RRM_PRI_1',
#       'LIG_14-3-3_CterR_2',
#       'LIG_eIF4E_2',
#       'LIG_UBA3_1',
#       'LIG_CNOT1_NIM_1'
#     )
#   ) |>
#   ggplot() +
#   stat_summary(
#     aes(
#       x = matched_motifs,
#       y = half_life
#     ),
#     fun.data = 'mean_cl_boot', 
#     color = color_1, 
#     linewidth = 0.2, 
#     size = 0.5
#   ) +
#   stat_pvalue_manual(
#     data = HEK_Nterm_sequence_ELM_motif_binding_wilcoxon_test |> slice(3, 4, 7), 
#     label = 'p.signif', label.size = 6, 
#     tip.length = 0, y.position = c(170, 180, 190), coord.flip = TRUE
#   ) +
#   labs(x = '', y = '') +
#   coord_flip(ylim = c(0, 200)) +
#   theme(
#     axis.text.x = element_text(color = 'black', size = 8, family = 'arial'),
#     axis.text.y = element_text(color = 'black', size = 8, family = 'arial')
#   )
# 
# ggsave(
#   filename = 'figures/figure7/point_range_Nterm_ELM_motifs_binding.eps',
#   plot = point_range_Nterm_ELM_motifs_binding,
#   height = 3, width = 2.6, units = 'in'
# )
# 
# ### figure 7C, degron-related ELM motif
# # Wilcoxon rank-sum test
# library(rstatix)
# 
# HEK_Nterm_sequence_ELM_motif_degron_wilcoxon_test <- HEK_Nterm_sequence_ELM_motif |> 
#   filter(str_detect(matched_motifs, 'DEG')) |> 
#   wilcox_test(half_life ~ matched_motifs) |> 
#   add_significance('p') |> 
#   filter(p < 0.05)
# 
# # point range plot
# library(ggpubr)
# 
# point_range_Nterm_ELM_motifs_degron <- HEK_Nterm_sequence_ELM_motif |> 
#   filter(
#     matched_motifs %in% c(
#       'DEG_APCC_KENBOX_2',
#       'DEG_Cend_FEM1B_2',
#       'DEG_Cend_KLHDC2_1',
#       'DEG_APCC_DBOX_1',
#       'DEG_COP1_1',
#       'DEG_Kelch_Keap1_1',
#       'DEG_Kelch_KLHL12_1',
#       'DEG_APCC_DBOX_1',
#       'DEG_ODPH_VHL_1',
#       'DEG_CRBN_cyclicCter_1',
#       'DEG_SCF_FBW7_1'
#     )
#   ) |> 
#   ggplot() +
#   stat_summary(
#     aes(
#       x = matched_motifs,
#       y = half_life
#     ),
#     fun.data = 'mean_cl_boot', 
#     color = color_2, 
#     linewidth = 0.2, 
#     size = 0.5
#   ) +
#   stat_pvalue_manual(
#     data = HEK_Nterm_sequence_ELM_motif_degron_wilcoxon_test |> slice(3, 14, 16), 
#     label = 'p.signif', label.size = 6, 
#     tip.length = 0, y.position = c(170, 190, 160), coord.flip = TRUE
#   ) +
#   labs(x = '', y = '') +
#   coord_flip(ylim = c(0, 200)) +
#   theme(
#     axis.text.x = element_text(color = 'black', size = 8, family = 'arial'),
#     axis.text.y = element_text(color = 'black', size = 8, family = 'arial')
#   )
# 
# ggsave(
#   filename = 'figures/figure7/point_range_Nterm_ELM_motifs_degron.eps',
#   plot = point_range_Nterm_ELM_motifs_degron,
#   height = 3, width = 2.7, units = 'in'
# )
# 
