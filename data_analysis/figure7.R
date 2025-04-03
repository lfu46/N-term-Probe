# import packages
library(tidyverse)

### figure 7A, ELM N-degrons
# Wilcoxon rank-sum test
library(rstatix)

HEK_Nterm_ELM_N_degron_wilcoxon_test <- HEK_Nterm_ELM_N_degron |>
  wilcox_test(half_life ~ ELM_N_degron) |> 
  add_significance('p') |> 
  filter(p.adj < 0.05) |> 
  slice(1, 2, 5, 8)

# point range plot
library(ggpubr)

violin_boxplot_Nterm_degron <- HEK_Nterm_ELM_N_degron |> 
  ggplot() +
  geom_violin(
    aes(
      x = fct_reorder(ELM_N_degron, half_life),
      y = half_life
    ),
    fill = color_1, color = 'transparent'
  ) +
  geom_boxplot(
    aes(
      x = fct_reorder(ELM_N_degron, half_life),
      y = half_life
    ), 
    color = 'black', width = 0.2, outliers = FALSE
  ) +
  stat_pvalue_manual(
    data = HEK_Nterm_ELM_N_degron_wilcoxon_test, label = 'p.adj.signif', label.size = 6, hide.ns = TRUE,
    tip.length = 0, y.position = c(48, 42, 45, 39)
  ) +
  labs(x = '', y = '') +
  coord_cartesian(ylim = c(0, 50)) +
  theme(
    panel.grid.major = element_line(color = "gray", linewidth = 0.2),
    panel.grid.minor = element_line(color = "gray", linewidth = 0.1),
    axis.text.x = element_text(color = 'black', size = 8, angle = 30, hjust = 1),
    axis.text.y = element_text(color = 'black', size = 8)
  )

ggsave(
  filename = 'figures/figure7/violin_boxplot_Nterm_degron.eps',
  plot = violin_boxplot_Nterm_degron,
  device = 'eps',
  height = 2, width = 2.2, units = 'in'
)

### figure 7B, insights from half-life, N-degron, protease, structure, motif and domain
library(tidyverse)

## combine all features
# half-life and N-degron
HEK_Nterm_ELM_N_degron <- read_csv(
  'data_source/ELM_degron/HEK_Nterm_ELM_N_degron.csv'
) |> 
  select(
    Index, UniProt_Accession, Gene, Entry.Name, 
    Kd_adj, half_life, category, Nterm_terminus, ELM_N_degron
  )

# protease
Nterm_topfinder_result_protease_feature <- read_delim(
  'data_source/Nterm_topfinder/2025_04_02_Nterm_degradation_04022025_R/2025_04_02_Nterm_degradation_04022025_R_Full_Table.txt',
  col_names = TRUE,
  name_repair = 'universal'
) |> 
  filter(Protein.found == 'YES') |> 
  mutate(Index = paste(Accession, P1..Position, sep = '_')) |> 
  select(
    Index, 
    Cleaving.proteases, 
    Distance.To.signal.peptide,
    Distance.to.last.transmembrane.domain..shed.,
    C.terminal.Features..P1..to.End.
  )

# structure
Nterm_degradation_alphafold_half_life <- read_csv(
  'data_source/Nterm_degradation_structuremap/Nterm_degradation_alphafold_half_life.csv'
) |> 
  select(Index, structure_group, IDR)

# combine all features
HEK_Nterm_N_degron_feature_comb <- HEK_Nterm_ELM_N_degron |> 
  left_join(Nterm_topfinder_result_protease_feature, by = 'Index') |> 
  left_join(Nterm_degradation_alphafold_half_life, by = 'Index')

write_csv(
  HEK_Nterm_N_degron_feature_comb,
  file = 'data_source/ELM_degron/HEK_Nterm_N_degron_feature_comb.csv'
)

## UBR-box, N-degron example
## Q13740_28, Q9Y4L1_33, Q92643_30, Q9NS69_106
# Q13740, ALCAM, Q13740_28
Q13740_database_info <- tribble(
  ~ name, ~ start, ~ end,
  'protein', 1, 583,
  'C2-set_2', 137, 230,
  'Ig_3', 254, 314
)

Q13740_result <- tibble(
  cleavage_site = c(28)
)

# example plot
Q13740_example <- ggplot() +
  geom_rect(
    data = Q13740_database_info,
    aes(
      xmin = start,
      xmax = end,
      ymin = 1,
      ymax = 2,
      fill = name, 
      color = name
    ),
    show.legend = FALSE
  ) +
  geom_segment(
    aes(
      x = Q13740_result$cleavage_site, 
      xend = Q13740_result$cleavage_site, 
      y = 2.6, 
      yend = 2
    ),
    arrow = arrow(length = unit(0.04, "in")), 
    color = "black",
    linewidth = 0.3
  ) +
  scale_fill_manual(
    values = c(
      'protein' = 'grey70',
      'C2-set_2' = color_1,
      'Ig_3' = color_2
    )
  ) +
  scale_color_manual(
    values = c(
      'protein' = 'black',
      'C2-set_2' = 'transparent',
      'Ig_3' = 'transparent'
    )
  ) +
  theme_void()

ggsave(
  filename = 'figures/figure7/Q13740_example.eps',
  height = 0.2, width = 2, units = 'in'
)

# Q9Y4L1, HYOU1, Q9Y4L1_33
Q9Y4L1_database_info <- tribble(
  ~ name, ~ start, ~ end,
  'protein', 1, 999,
  'HSP70', 36, 632
)

Q9Y4L1_result <- tibble(
  cleavage_site = c(33)
)

# example plot
Q9Y4L1_example <- ggplot() +
  geom_rect(
    data = Q9Y4L1_database_info,
    aes(
      xmin = start,
      xmax = end,
      ymin = 1,
      ymax = 2,
      fill = name, 
      color = name
    ),
    show.legend = FALSE
  ) +
  geom_segment(
    aes(
      x = Q9Y4L1_result$cleavage_site, 
      xend = Q9Y4L1_result$cleavage_site, 
      y = 2.6, 
      yend = 2
    ),
    arrow = arrow(length = unit(0.04, "in")), 
    color = "black",
    linewidth = 0.3
  ) +
  scale_fill_manual(
    values = c(
      'protein' = 'grey70',
      'HSP70' = color_3
    )
  ) +
  scale_color_manual(
    values = c(
      'protein' = 'black',
      'HSP70' = 'transparent'
    )
  ) +
  theme_void()

ggsave(
  filename = 'figures/figure7/Q9Y4L1_example.eps',
  height = 0.2, width = 2, units = 'in'
)

# Q92643, PIGK, Q92643_30
Q92643_database_info <- tribble(
  ~ name, ~ start, ~ end,
  'protein', 1, 395,
  'Peptidase_C13', 45, 277
)

Q92643_result <- tibble(
  cleavage_site = c(30)
)

# example plot
Q92643_example <- ggplot() +
  geom_rect(
    data = Q92643_database_info,
    aes(
      xmin = start,
      xmax = end,
      ymin = 1,
      ymax = 2,
      fill = name, 
      color = name
    ),
    show.legend = FALSE
  ) +
  geom_segment(
    aes(
      x = Q92643_result$cleavage_site, 
      xend = Q92643_result$cleavage_site, 
      y = 2.6, 
      yend = 2
    ),
    arrow = arrow(length = unit(0.04, "in")), 
    color = "black",
    linewidth = 0.3
  ) +
  scale_fill_manual(
    values = c(
      'protein' = 'grey70',
      'Peptidase_C13' = color_4
    )
  ) +
  scale_color_manual(
    values = c(
      'protein' = 'black',
      'Peptidase_C13' = 'transparent'
    )
  ) +
  theme_void()

ggsave(
  filename = 'figures/figure7/Q92643_example.eps',
  height = 0.2, width = 2, units = 'in'
)

# Q9NS69, TOMM22, Q9NS69_106
Q9NS69_database_info <- tribble(
  ~ name, ~ start, ~ end,
  'protein', 1, 142,
  'Tom22', 26, 117
)

Q9NS69_result <- tibble(
  cleavage_site = c(106)
)

# example plot
Q9NS69_example <- ggplot() +
  geom_rect(
    data = Q9NS69_database_info,
    aes(
      xmin = start,
      xmax = end,
      ymin = 1,
      ymax = 2,
      fill = name, 
      color = name
    ),
    show.legend = FALSE
  ) +
  geom_segment(
    aes(
      x = Q9NS69_result$cleavage_site, 
      xend = Q9NS69_result$cleavage_site, 
      y = 2.6, 
      yend = 2
    ),
    arrow = arrow(length = unit(0.04, "in")), 
    color = "black",
    linewidth = 0.3
  ) +
  scale_fill_manual(
    values = c(
      'protein' = 'grey70',
      'Tom22' = color_5
    )
  ) +
  scale_color_manual(
    values = c(
      'protein' = 'black',
      'Tom22' = 'transparent'
    )
  ) +
  theme_void()

ggsave(
  filename = 'figures/figure7/Q9NS69_example.eps',
  height = 0.2, width = 2, units = 'in'
)

## ZER1, N-degron example
# P26368_129, Q13263_686
# P26368, U2AF2, P26368_129
P26368_database_info <- tribble(
  ~ name, ~ start, ~ end,
  'protein', 1, 475,
  'RRM_1', 151, 224,
  'RRM_1', 261, 330,
  'RRM_1', 400, 459
)

P26368_result <- tibble(
  cleavage_site = c(129)
)

# example plot
P26368_example <- ggplot() +
  geom_rect(
    data = P26368_database_info,
    aes(
      xmin = start,
      xmax = end,
      ymin = 1,
      ymax = 2,
      fill = name, 
      color = name
    ),
    show.legend = FALSE
  ) +
  geom_segment(
    aes(
      x = P26368_result$cleavage_site, 
      xend = P26368_result$cleavage_site, 
      y = 2.6, 
      yend = 2
    ),
    arrow = arrow(length = unit(0.04, "in")), 
    color = "black",
    linewidth = 0.3
  ) +
  scale_fill_manual(
    values = c(
      'protein' = 'grey70',
      'RRM_1' = color_1
    )
  ) +
  scale_color_manual(
    values = c(
      'protein' = 'black',
      'RRM_1' = 'transparent'
    )
  ) +
  theme_void()

ggsave(
  filename = 'figures/figure7/P26368_example.eps',
  height = 0.2, width = 2, units = 'in'
)

# Q13263, TRIM28, Q13263_686
Q13263_database_info <- tribble(
  ~ name, ~ start, ~ end,
  'protein', 1, 835,
  'PHD', 628, 669,
  'zf-B_box', 150, 194,
  'zf-B_box', 208, 243,
  'zf-RING_5', 64, 122
)

Q13263_result <- tibble(
  cleavage_site = c(686)
)

# example plot
Q13263_example <- ggplot() +
  geom_rect(
    data = Q13263_database_info,
    aes(
      xmin = start,
      xmax = end,
      ymin = 1,
      ymax = 2,
      fill = name, 
      color = name
    ),
    show.legend = FALSE
  ) +
  geom_segment(
    aes(
      x = Q13263_result$cleavage_site, 
      xend = Q13263_result$cleavage_site, 
      y = 2.6, 
      yend = 2
    ),
    arrow = arrow(length = unit(0.04, "in")), 
    color = "black",
    linewidth = 0.3
  ) +
  scale_fill_manual(
    values = c(
      'protein' = 'grey70',
      'PHD' = color_2,
      'zf-B_box' = color_3,
      'zf-RING_5' = color_4
    )
  ) +
  scale_color_manual(
    values = c(
      'protein' = 'black',
      'PHD' = 'transparent',
      'zf-B_box' = 'transparent',
      'zf-RING_5' = 'transparent'
    )
  ) +
  theme_void()

ggsave(
  filename = 'figures/figure7/Q13263_example.eps',
  height = 0.2, width = 2, units = 'in'
)

## NCPR for UBR-box and ZER1 N-recognin
library(tidyverse)
library(reticulate)
library(Biostrings)

## get amino acid sequences for UBR-box and ZER1
# import human fasta downloaded from UniProt (https://www.uniprot.org/)
human_fasta <- readAAStringSet(
  'data_source/fasta_file/uniprotkb_reviewed_true_AND_model_organ_2025_02_15.fasta'
)

# build tibble using human fasta
human_fasta_tibble <- tibble(
  Name = names(human_fasta),
  Sequence = as.character(human_fasta),
  Length = width(human_fasta)
) |> 
  mutate(
    Name = sub(' .*', '', Name),
    Full_Protein_Length = as.numeric(Length)
  ) |> 
  separate(Name, into = c('sp', 'UniProt_Accession', 'name'), sep = '\\|') |> 
  select(UniProt_Accession, Sequence, Full_Protein_Length)

# binding region sequence
N_recognin_binding_region_sequence <- human_fasta_tibble |> 
  filter(
    UniProt_Accession %in% c(
      'Q8IWV7', # UBR1
      'Q7Z7L7' # ZER1
    )
  ) |> 
  mutate(
    start = c(401, 99),
    end = c(755, 166),
    binding_region = substr(Sequence, start = start, stop = end)
  )

write_csv(
  N_recognin_binding_region_sequence,
  file = 'data_source/ELM_degron/N_recognin_binding_region_sequence.csv'
)

# use specific vitual env
library(reticulate)

use_condaenv(
  condaenv = '/opt/anaconda3/envs/Nterm_probe',
  required = TRUE
)

# excute python code for NCPR of binding region sequence
source_python('data_analysis/N_recognin_NCPR.py')

# convert list result to tibble
library(tidyverse)

NCPR_result_tibble_list <- lapply(names(NCPR_result), function(protein_id) {
  tibble(UniProt_Accession = protein_id, NCPR = NCPR_result[[protein_id]])
})

NCPR_result_tibble <- bind_rows(NCPR_result_tibble_list) |> 
  mutate(NCPR = lapply(NCPR, function(x) as.data.frame(t(x)))) |> 
  unnest(NCPR)

colnames(NCPR_result_tibble) <- c('UniProt_Accession', 'Position', 'NCPR', 'charge')

NCPR_result_tibble_adj <- NCPR_result_tibble |> 
  mutate(
    charge = case_when(
      NCPR > 0 ~ 'positive',
      NCPR < 0 ~ 'negative',
      NCPR == 0 ~ 'none'
    )
  )

## NCPR distribution for each N-recoginin binding region
# Q8IWV7, UBR1
NCPR_UBR1 <- NCPR_result_tibble_adj |> 
  filter(UniProt_Accession == 'Q8IWV7') |> 
  ggplot() +
  geom_col(
    aes(
      x = Position,
      y = NCPR,
      color = charge
    ),
    show.legend = FALSE
  ) +
  scale_color_manual(
    values = c(
      'positive' = color_1,
      'negative' = color_2,
      'none' = 'transparent'
      )
  ) +
  labs(x = '', y = '') +
  coord_cartesian(ylim = c(-0.8, 0.6)) +
  theme(
    panel.grid.major = element_line(color = "gray", linewidth = 0.2),
    panel.grid.minor = element_line(color = "gray", linewidth = 0.1),
    axis.text.x = element_text(family = 'arial', color = 'black', size = 6),
    axis.text.y = element_text(family = 'arial', color = 'black', size = 6)
  )

ggsave(
  file = 'figures/figure7/NCPR_UBR1.eps',
  plot = NCPR_UBR1,
  height = 1, width = 1.8, units = 'in'
)

# Q7Z7L7, ZER1
NCPR_ZER1 <- NCPR_result_tibble_adj |> 
  filter(UniProt_Accession == 'Q7Z7L7') |> 
  ggplot() +
  geom_col(
    aes(
      x = Position,
      y = NCPR,
      color = charge
    ),
    show.legend = FALSE
  ) +
  scale_color_manual(
    values = c(
      'positive' = color_1,
      'negative' = color_2,
      'none' = 'transparent'
    )
  ) +
  labs(x = '', y = '') +
  coord_cartesian(ylim = c(-0.8, 0.6)) +
  theme(
    panel.grid.major = element_line(color = "gray", linewidth = 0.2),
    panel.grid.minor = element_line(color = "gray", linewidth = 0.1),
    axis.text.x = element_text(family = 'arial', color = 'black', size = 6),
    axis.text.y = element_text(family = 'arial', color = 'black', size = 6)
  )

ggsave(
  file = 'figures/figure7/NCPR_ZER1.eps',
  plot = NCPR_ZER1,
  height = 1, width = 1.8, units = 'in'
)

### figure 7C, half-life comparison of cleaving proteases
# Wilcoxon rank-sum test
library(rstatix)

Cleaving.Proteases.List <- Nterm_degron_topfinder_cleaving_proteases |> 
  group_by(Cleaving.proteases) |> 
  get_summary_stats(half_life, type = 'median') |> 
  filter(n > 1) |> 
  pull(Cleaving.proteases)

Nterm_topfiner_cleaving_proteases_wilcoxon_test <- Nterm_degron_topfinder_cleaving_proteases |> 
  filter(Cleaving.proteases %in% Cleaving.Proteases.List) |> 
  wilcox_test(half_life ~ Cleaving.proteases) |> 
  add_significance('p') |> 
  filter(p < 0.05)

Nterm_cleaving_protease_list <- c(
  'CATS',
  'CATB',
  
  'GRAB',
  
  'MEP1A',
  'MEP1B',

  'MMP11',
  'HTRA2',
  
  'CASP1',
  'CASP3'
)

# boxplot
library(ggpubr)

point_range_plot_Nterm_cleaving_proteases <- Nterm_degron_topfinder_cleaving_proteases |> 
  filter(Cleaving.proteases %in% Nterm_cleaving_protease_list) |> 
  ggplot() +
  stat_summary(
    aes(
      x = fct_reorder(Cleaving.proteases, half_life), 
      y = half_life
    ),
    fun.data = 'mean_cl_boot', color = color_3, linewidth = 0.2, size = 0.5
  ) +
  labs(x = '', y = '') +
  stat_pvalue_manual(
    data = Nterm_topfiner_cleaving_proteases_wilcoxon_test |> slice(2, 4),
    label = 'p.signif',
    y.position = c(35, 32),
    tip.length = 0,
    label.size = 6
  ) +
  coord_cartesian(ylim = c(10, 38)) +
  theme(
    panel.grid.major = element_line(color = "gray", linewidth = 0.2),
    panel.grid.minor = element_line(color = "gray", linewidth = 0.1),
    axis.text.x = element_text(color = 'black', size = 8, angle = 30, hjust = 1),
    axis.text.y = element_text(color = 'black', size = 8)
  )

ggsave(
  filename = 'figures/figure7/point_range_plot_Nterm_cleaving_proteases.eps',
  plot = point_range_plot_Nterm_cleaving_proteases,
  device = 'eps',
  height = 2, width = 2.2, units = 'in'
)

### figure 7D, ELM motifs
library(tidyverse)

# import GSEA result from ELM motif analysis
Nterm_ELM_motif_GSEA_des_Kd <- read_csv(
  'data_source/ELM_degron/Nterm_ELM_motif_GSEA_des_Kd.csv'
)

# filter the enriched term which p < 0.05, rearrange the terms and export as a table
library(gridExtra)

Nterm_ELM_motif_GSEA_des_Kd_enriched_term <- Nterm_ELM_motif_GSEA_des_Kd |> 
  filter(
    Description %in% c(
      'DOC_AGCK_PIF_1',
      'LIG_PDZ_Class_2',
      'LIG_CaMK_CASK_1',
      'LIG_HCF-1_HBM_1',
      'DEG_CRL4_CDT2_2'
    )
  ) |> 
  select(
    'ELM.Name' = Description, NES, pvalue
  ) |> 
  mutate(
    NES = round(NES, 2) , 
    'P.value' = format(pvalue, format = "e", digits = 2)
  ) |> 
  select(-pvalue) |> 
  arrange(P.value)

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
