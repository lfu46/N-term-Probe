# import packages
library(tidyverse)

### figure 7A, ELM N-degrons
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
  filename = 'figures/figure7/point_range_Nterm_degron.eps',
  plot = point_range_Nterm_degron,
  device = cairo_ps,
  height = 2, width = 2.5, units = 'in',
  fallback_resolution = 1200
)

### figure 7B, ELM motifs
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
  filename = 'figures/figure7/ranking_plot_Nterm_Kd.eps',
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
  filename = 'figures/figure7/protein_example_without_motif.eps',
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
  filename = 'figures/figure7/protein_example_with_motif.eps',
  plot = protein_example_with_motif,
  height = 1, width = 4, units = 'in'
)

### figure 7C, degron-related ELM motif
# Wilcoxon rank-sum test
library(rstatix)

HEK_Nterm_sequence_ELM_motif_degron_wilcoxon_test <- HEK_Nterm_sequence_ELM_motif |> 
  filter(str_detect(matched_motifs, 'DEG')) |> 
  wilcox_test(half_life ~ matched_motifs) |> 
  add_significance('p') |> 
  filter(p < 0.05)

# point range plot
library(ggpubr)

point_range_Nterm_ELM_motifs_degron <- HEK_Nterm_sequence_ELM_motif |> 
  filter(
    matched_motifs %in% c(
      'DEG_APCC_KENBOX_2',
      'DEG_Cend_FEM1B_2',
      'DEG_Cend_KLHDC2_1',
      'DEG_APCC_DBOX_1',
      'DEG_COP1_1',
      'DEG_Kelch_Keap1_1',
      'DEG_Kelch_KLHL12_1',
      'DEG_APCC_DBOX_1',
      'DEG_ODPH_VHL_1',
      'DEG_CRBN_cyclicCter_1',
      'DEG_SCF_FBW7_1'
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
  stat_pvalue_manual(
    data = HEK_Nterm_sequence_ELM_motif_degron_wilcoxon_test |> slice(3, 14, 16), 
    label = 'p.signif', label.size = 6, 
    tip.length = 0, y.position = c(170, 190, 160), coord.flip = TRUE
  ) +
  labs(x = '', y = '') +
  coord_flip(ylim = c(0, 200)) +
  theme(
    axis.text.x = element_text(color = 'black', size = 8, family = 'arial'),
    axis.text.y = element_text(color = 'black', size = 8, family = 'arial')
  )

ggsave(
  filename = 'figures/figure7/point_range_Nterm_ELM_motifs_degron.eps',
  plot = point_range_Nterm_ELM_motifs_degron,
  height = 3, width = 2.7, units = 'in'
)

### figure 7D, docking-related ELM motif
# Wilcoxon rank-sum test
library(rstatix)

HEK_Nterm_sequence_ELM_motif_docking_wilcoxon_test <- HEK_Nterm_sequence_ELM_motif |> 
  filter(str_detect(matched_motifs, 'DOC')) |> 
  wilcox_test(half_life ~ matched_motifs) |> 
  add_significance('p') |> 
  filter(p < 0.05)

# point range plot
library(ggpubr)

point_range_Nterm_ELM_motifs_docking <- HEK_Nterm_sequence_ELM_motif |> 
  filter(
    matched_motifs %in% c(
      'DOC_MIT_MIM_1',
      'DOC_MAPK_gen_1',
      'DOC_PP2B_LxvP_1',
      'DOC_USP7_MATH_1',
      'DOC_USP7_MATH_2',
      'DOC_PP1_RVXF_1',
      'DOC_CKS1_1',
      'DOC_CYCLIN_RxL_1',
      'DOC_PP4_MxPP_1',
      'DOC_RSK_DDVF_1'
    )
  ) |> 
  ggplot() +
  stat_summary(
    aes(
      x = matched_motifs,
      y = half_life
    ),
    fun.data = 'mean_cl_boot', 
    color = color_3, 
    linewidth = 0.2, 
    size = 0.5
  ) +
  stat_pvalue_manual(
    data = HEK_Nterm_sequence_ELM_motif_docking_wilcoxon_test |> slice(2, 4, 6), 
    label = 'p.signif', label.size = 6, 
    tip.length = 0, y.position = c(100, 110, 120), coord.flip = TRUE
  ) +
  labs(x = '', y = '') +
  coord_flip(ylim = c(0, 130)) +
  theme(
    axis.text.x = element_text(color = 'black', size = 8, family = 'arial'),
    axis.text.y = element_text(color = 'black', size = 8, family = 'arial')
  )

ggsave(
  filename = 'figures/figure7/point_range_Nterm_ELM_motifs_docking.eps',
  plot = point_range_Nterm_ELM_motifs_docking,
  height = 3, width = 2.5, units = 'in'
)

### figure 7E, targeting-related ELM motif
# Wilcoxon rank-sum test
library(rstatix)

HEK_Nterm_sequence_ELM_motif_targeting_wilcoxon_test <- HEK_Nterm_sequence_ELM_motif |> 
  filter(str_detect(matched_motifs, 'TRG')) |> 
  wilcox_test(half_life ~ matched_motifs) |> 
  add_significance('p') |> 
  filter(p < 0.05)

# point range plot
library(ggpubr)

point_range_Nterm_ELM_motifs_targeting <- HEK_Nterm_sequence_ELM_motif |> 
  filter(
    matched_motifs %in% c(
      'TRG_ER_diArg_1',
      'TRG_NLS_Bipartite_1',
      'TRG_PTS1',
      'TRG_ENDOCYTIC_2',
      'TRG_ER_diLys_1',
      'TRG_ER_FFAT_1',
      'TRG_ER_KDEL_1',
      'TRG_NES_CRM1_1',
      'TRG_NLS_MonoCore_2',
      'TRG_NLS_MonoExtC_3',
      'TRG_NLS_MonoExtN_4'
    )
  ) |>
  ggplot() +
  stat_summary(
    aes(
      x = matched_motifs,
      y = half_life
    ),
    fun.data = 'mean_cl_boot', 
    color = color_4, 
    linewidth = 0.2, 
    size = 0.5
  ) +
  stat_pvalue_manual(
    data = HEK_Nterm_sequence_ELM_motif_targeting_wilcoxon_test |> slice(12, 19, 20, 21, 22), 
    label = 'p.signif', label.size = 6, 
    tip.length = 0, y.position = c(160, 120, 130, 140, 150), coord.flip = TRUE
  ) +
  labs(x = '', y = '') +
  coord_flip(ylim = c(0, 170)) +
  theme(
    axis.text.x = element_text(color = 'black', size = 8, family = 'arial'),
    axis.text.y = element_text(color = 'black', size = 8, family = 'arial')
  )

ggsave(
  filename = 'figures/figure7/point_range_Nterm_ELM_motifs_targeting.eps',
  plot = point_range_Nterm_ELM_motifs_targeting,
  height = 3.2, width = 2.6, units = 'in'
)

### figure 7F, binding-related ELM motif
# Wilcoxon rank-sum test
library(rstatix)

HEK_Nterm_sequence_ELM_motif_binding_wilcoxon_test <- HEK_Nterm_sequence_ELM_motif |> 
  filter(
    matched_motifs %in% c(
      'LIG_APCC_Cbox_2',
      'LIG_AP_GAE_1',
      'LIG_G3BP_FGDF_1',
      'LIG_PTB_Phospho_1',
      'LIG_NBox_RRM_1',
      'LIG_RRM_PRI_1',
      'LIG_14-3-3_CterR_2',
      'LIG_eIF4E_2',
      'LIG_UBA3_1',
      'LIG_CNOT1_NIM_1'
    )
  ) |> 
  wilcox_test(half_life ~ matched_motifs) |> 
  add_significance('p') |> 
  filter(p < 0.05)

# point range plot
library(ggpubr)

point_range_Nterm_ELM_motifs_binding <- HEK_Nterm_sequence_ELM_motif |> 
  filter(
    matched_motifs %in% c(
      'LIG_APCC_Cbox_2',
      'LIG_AP_GAE_1',
      'LIG_G3BP_FGDF_1',
      'LIG_PTB_Phospho_1',
      'LIG_NBox_RRM_1',
      'LIG_RRM_PRI_1',
      'LIG_14-3-3_CterR_2',
      'LIG_eIF4E_2',
      'LIG_UBA3_1',
      'LIG_CNOT1_NIM_1'
    )
  ) |>
  ggplot() +
  stat_summary(
    aes(
      x = matched_motifs,
      y = half_life
    ),
    fun.data = 'mean_cl_boot', 
    color = color_1, 
    linewidth = 0.2, 
    size = 0.5
  ) +
  stat_pvalue_manual(
    data = HEK_Nterm_sequence_ELM_motif_binding_wilcoxon_test |> slice(3, 4, 7), 
    label = 'p.signif', label.size = 6, 
    tip.length = 0, y.position = c(170, 180, 190), coord.flip = TRUE
  ) +
  labs(x = '', y = '') +
  coord_flip(ylim = c(0, 200)) +
  theme(
    axis.text.x = element_text(color = 'black', size = 8, family = 'arial'),
    axis.text.y = element_text(color = 'black', size = 8, family = 'arial')
  )

ggsave(
  filename = 'figures/figure7/point_range_Nterm_ELM_motifs_binding.eps',
  plot = point_range_Nterm_ELM_motifs_binding,
  height = 3, width = 2.6, units = 'in'
)
