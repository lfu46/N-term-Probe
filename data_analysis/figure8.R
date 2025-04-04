# import packages
library(tidyverse)

### figure 8A, GO and KEGG analysis for top20% and bottom 20% proteins
## import GO and KEGG enrichment results
top20_comparison_GO <- read_csv(
  'data_source/Nterm_WP_comparison/top20_comparison_GO.csv'
)

bottom20_comparison_GO <- read_csv(
  'data_source/Nterm_WP_comparison/bottom20_comparison_GO.csv'
)

top20_comparison_KEGG <- read_csv(
  'data_source/Nterm_WP_comparison/top20_comparison_KEGG.csv'
)

bottom20_comparison_KEGG <- read_csv(
  'data_source/Nterm_WP_comparison/bottom20_comparison_KEGG.csv'
)

## dot plot
top20_bottom20_comb <- bind_rows(
  top20_comparison_GO |> 
    filter(
      Description %in% c(
        'Histone modifying activity',
        'Catalytic activity, acting on RNA',
        'Kinase activity',
        'Transferase activity',
        'Macromolecule modification'
      )
    ) |> 
    mutate(
      category = 'top20'
    ),
  
  bottom20_comparison_GO |> 
    filter(
      Description %in% c(
        'Unfolded protein binding',
        'Protein folding chaperone',
        'mRNA binding',
        'Regulation of mRNA metabolic process',
        'Catalytic step 2 spliceosome'
      )
    ) |> 
    mutate(
      category = 'bottom20'
    )
)

write_csv(
  top20_bottom20_comb,
  file = 'data_source/Nterm_WP_comparison/top20_bottom20_comb.csv'
)

dotplot_top20_bottom20_GO_KEGG <- top20_bottom20_comb |> 
  mutate(
    FoldEnrichment = ifelse(category == 'bottom20', -FoldEnrichment, FoldEnrichment)
  ) |> 
  ggplot() +
  geom_point(
    aes(
      x = FoldEnrichment,
      y = -log10(p.adjust),
      size = Count,
      fill = category
    ),
    shape = 21
  ) +
  geom_vline(xintercept = 0, color = 'black', linetype = 'dashed') +
  theme(
    axis.text = element_text(size = 8, color = 'black', family = 'arial'),
    axis.title = element_text(size = 8, color = 'black', family = 'arial'),
    legend.title = element_text(size = 8, color = 'black', family = 'arial'),
    legend.text = element_text(size = 8, color = 'black', family = 'arial')
  )

ggsave(
  filename = 'figures/figure8/dotplot_top20_bottom20_GO_KEGG.eps',
  plot = dotplot_top20_bottom20_GO_KEGG,
  height = 2, width = 3.5, units = 'in'
)

### figure 8B, protein complex and domain analysis for top20% and bottom 20% proteins
## import enrichment results
top20_comparison_corum <- read_csv(
  'data_source/Nterm_WP_comparison/top20_comparison_corum.csv'
)

top20_comparison_pfam <- read_csv(
  'data_source/Nterm_WP_comparison/top20_comparison_pfam.csv'
)

bottom20_comparison_corum <- read_csv(
  'data_source/Nterm_WP_comparison/bottom20_comparison_corum.csv'
)

bottom20_comparison_pfam <- read_csv(
  'data_source/Nterm_WP_comparison/bottom20_comparison_pfam.csv'
)

## combine results from corum and pfam
# top20
top20_comparison_corum_pfam_comb <- top20_comparison_corum |> 
  filter(
    Description %in% c(
      'corum_id_308'
    )
  ) |> 
  bind_rows(
    top20_comparison_pfam |> 
      filter(
        Description %in% c(
          'PF00595', 
          'PF00069', 
          'PF00271', 
          'PF07653'
        )
      )
  )

# bottom20
bottom20_comparison_corum_pfam_comb <- bottom20_comparison_corum |> 
  filter(
    Description %in% c(
      'corum_id_8372',
      'corum_id_1181',
      'corum_id_3082',
      'corum_id_8391'
    )
  ) |> 
  bind_rows(
    bottom20_comparison_pfam |> 
      filter(
        Description %in% c(
          'PF00076'
        )
      )
  )

## bar plot
# top20
barplot_top20_corum_pfam <- top20_comparison_corum_pfam_comb |> 
  ggplot() +
  geom_bar(
    aes(
      x = fct_reorder(Description, -pvalue),
      y = -log10(pvalue)
    ),
    stat = 'identity', fill = color_3
  ) +
  labs(x = '', y = '') +
  coord_flip() +
  theme(
    axis.text = element_text(size = 8, color = 'black')
  )

ggsave(
  filename = 'figures/figure8/barplot_top20_corum_pfam.eps',
  plot = barplot_top20_corum_pfam,
  height = 2, width = 2.5, units = 'in'
)

# bottom20
barplot_bottom20_corum_pfam <- bottom20_comparison_corum_pfam_comb |> 
  ggplot() +
  geom_bar(
    aes(
      x = fct_reorder(Description, -pvalue),
      y = -log10(pvalue)
    ),
    stat = 'identity', fill = color_4
  ) +
  labs(x = '', y = '') +
  coord_flip() +
  theme(
    axis.text = element_text(size = 8, color = 'black')
  )

ggsave(
  filename = 'figures/figure8/barplot_bottom20_corum_pfam.eps',
  plot = barplot_bottom20_corum_pfam,
  height = 2, width = 2.5, units = 'in'
)

### figure 8C, kinase family
# import kinase family data
# http://www.kinase.com./web/current/human/
# Update: (Dec 07): An updated table (Excel) includes Gene IDs, refseq accessions, HGNC names and improved sequences for some kinases.
library(tidyverse)
library(readxl)

kinase_family <- read_xls(
  'data_source/Nterm_WP_comparison/Kincat_Hsap.08.02.xls',
  col_names = TRUE
) |> 
  select(Entrez_Symbol, Group, Family)

kinase_family_delta_half_life <- HEK_Nterm_WP_delta_half_life |> 
  left_join(
    kinase_family,
    by = c('Gene' = 'Entrez_Symbol'),
    relationship = "many-to-many"
  ) |> 
  filter(
    !is.na(Group)
  ) |> 
  distinct()

write_csv(
  kinase_family_delta_half_life,
  file = 'data_source/Nterm_WP_comparison/kinase_family_delta_half_life.csv'
)

# Wilcoxon rank-sum test
library(rstatix)

kinase_family_delta_half_life |> 
  group_by(Group) |>
  get_summary_stats(delta_half_life)

kinase_family_wilcoxon_test <- kinase_family_delta_half_life |> 
  wilcox_test(delta_half_life ~ Group) |> 
  filter(
    p.adj < 0.05
  )

# point range plot
library(ggpubr)

point_range_plot_kinase_family_delta_half_life <- kinase_family_delta_half_life |> 
  ggplot() +
  stat_summary(
    aes(
      x = fct_reorder(Group, delta_half_life, .fun = mean),
      y = delta_half_life
    ),
    color = color_1,
    fun.data = 'mean_cl_boot', linewidth = 0.2, size = 0.5, show.legend = FALSE
  ) +
  labs(x = '', y = '') +
  stat_pvalue_manual(
    data = kinase_family_wilcoxon_test |> slice(2, 7),
    label = 'p.adj.signif',
    y.position = c(0, -4),
    tip.length = 0,
    label.size = 6
  ) +
  coord_cartesian(
    ylim = c(-80, 5)
  ) +
  theme(
    panel.grid.major = element_line(color = "gray", linewidth = 0.2),
    panel.grid.minor = element_line(color = "gray", linewidth = 0.1),
    axis.text.x = element_text(color = 'black', size = 8, family = 'arial', angle = 60, hjust = 1),
    axis.text.y = element_text(color = 'black', size = 8, family = 'arial')
  )

ggsave(
  filename = 'figures/figure8/point_range_plot_kinase_family_delta_half_life.eps',
  device = 'eps',
  plot = point_range_plot_kinase_family_delta_half_life,
  height = 2, width = 2, units = 'in'
)

### figure 8D, Q05397_890, PTK2, Focal adhesion kinase 1
Q05397_database_info <- tribble(
  ~ name, ~ start, ~ end,
  'protein', 1, 1052,
  'FERM_M', 136, 249,
  'Focal_AT', 914, 1045,
  'PK_Tyr_Ser-Thr', 422, 676,
  'FERM_N_2', 35, 130,
  'FERM_C_FAK1', 260, 357
)

Q05397_result <- tibble(
  cleavage_site = c(890)
)

# example plot
Q05397_example <- ggplot() +
  geom_rect(
    data = Q05397_database_info,
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
      x = Q05397_result$cleavage_site, 
      xend = Q05397_result$cleavage_site, 
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
      'FERM_N_2' = color_1,
      'FERM_C_FAK1' = color_2,
      'FERM_M' = color_3,
      'Focal_AT' = color_4,
      'PK_Tyr_Ser-Thr' = color_5
    )
  ) +
  scale_color_manual(
    values = c(
      'protein' = 'black',
      'PK_Tyr_Ser-Thr' = 'transparent',
      'FERM_C_FAK1' = 'transparent',
      'FERM_N_2' = 'transparent',
      'FERM_M' = 'transparent',
      'Focal_AT' = 'transparent'
    )
  ) +
  theme_void()

ggsave(
  filename = 'figures/figure8/Q05397_example.eps',
  height = 0.2, width = 2, units = 'in'
)

## Q05397_890, PTK2, Focal adhesion kinase 1
## Nterm
Nterm_Q05397_890_deg_ratio <- HEK_Nterm_deg_ratio_adj |> 
  filter(Index == 'Q05397_890')

Nterm_Q05397_890_linear_model <- HEK_Nterm_linear_fitting |> 
  filter(Index == 'Q05397_890') |> 
  pivot_wider(names_from = parameters, values_from = values)

Nterm_Q05397_890_linear_model_adj <- HEK_Nterm_Kd_half_life_LaminB_Tcomplex |> 
  filter(Index == 'Q05397_890')

time_point <- c(0, 3, 6, 9, 12, 24)
Nterm_Q05397_890_time_point <- exp(Nterm_Q05397_890_linear_model$lnA - Nterm_Q05397_890_linear_model$Kd * time_point)
Nterm_Q05397_890_time_point_adj <- exp(Nterm_Q05397_890_linear_model$lnA - Nterm_Q05397_890_linear_model_adj$Kd_adj * time_point)

adjusted_raito <- Nterm_Q05397_890_time_point_adj / Nterm_Q05397_890_time_point
Nterm_Q05397_890_deg_ratio_adj <- Nterm_Q05397_890_deg_ratio |> 
  mutate(
    deg_ratio_adj_adj = deg_ratio_adj * adjusted_raito
  ) |> 
  mutate(
    rel_intensity = deg_ratio_adj_adj / max(deg_ratio_adj_adj)
  )

time_series <- seq(0, 24, length.out = 100)
Nterm_Q05397_890_linear_fitting_adj <- exp(Nterm_Q05397_890_linear_model$lnA - Nterm_Q05397_890_linear_model_adj$Kd_adj * time_series)

Nterm_Q05397_890_linear_fitting_line <- tibble(
  timepoint = time_series,
  deg_ratio = Nterm_Q05397_890_linear_fitting_adj
)

## WP
WP_Q05397_deg_ratio <- HEK_WP_deg_ratio_adj |> 
  filter(UniProt_Accession == 'Q05397')

WP_Q05397_linear_model <- HEK_WP_linear_fitting |>
  filter(UniProt_Accession == 'Q05397') |> 
  pivot_wider(names_from = parameters, values_from = values)

WP_Q05397_linear_model_adj <- HEK_WP_Kd_half_life_LaminB_Tcomplex |>
  filter(UniProt_Accession == 'Q05397')

time_point <- c(0, 3, 6, 9, 12, 24)
WP_Q05397_time_point <- exp(WP_Q05397_linear_model$lnA - WP_Q05397_linear_model$Kd * time_point)
WP_Q05397_time_point_adj <- exp(WP_Q05397_linear_model$lnA - WP_Q05397_linear_model_adj$Kd_adj * time_point)

adjusted_raito <- WP_Q05397_time_point_adj / WP_Q05397_time_point
WP_Q05397_deg_ratio_adj <- WP_Q05397_deg_ratio |> 
  mutate(
    deg_ratio_adj_adj = deg_ratio_adj * adjusted_raito
  ) |> 
  mutate(
    rel_intensity = deg_ratio_adj_adj / max(deg_ratio_adj_adj)
  )

time_series <- seq(0, 24, length.out = 100)
WP_Q05397_linear_fitting_adj <- exp(WP_Q05397_linear_model$lnA - WP_Q05397_linear_model_adj$Kd_adj * time_series)

WP_Q05397_linear_fitting_line <- tibble(
  timepoint = time_series,
  deg_ratio = WP_Q05397_linear_fitting_adj
)

# plot
Nterm_WP_Q05397_890_linear_example <- ggplot() +
  geom_point(
    data = Nterm_Q05397_890_deg_ratio_adj,
    aes(
      x = timepoint,
      y = rel_intensity
    ),
    shape = 21, fill = color_1, color = 'transparent', size = 1.2
  ) +
  geom_point(
    data = WP_Q05397_deg_ratio_adj,
    aes(
      x = timepoint,
      y = rel_intensity
    ),
    shape = 21, fill = color_2, color = 'transparent', size = 1.2
  ) +
  geom_line(
    data = Nterm_Q05397_890_linear_fitting_line,
    aes(
      x = timepoint,
      y = deg_ratio
    )
  ) +
  geom_line(
    data = WP_Q05397_linear_fitting_line,
    aes(
      x = timepoint,
      y = deg_ratio
    )
  ) +
  labs(x = NULL, y = NULL) +
  coord_cartesian(ylim = c(0, 1)) +
  theme(
    axis.text = element_text(size = 8, color = 'black', family = 'arial')
  )

ggsave(
  filename = 'figures/figure8/Nterm_WP_Q05397_890_linear_example.eps',
  plot = Nterm_WP_Q05397_890_linear_example,
  height = 1, width = 2, units = 'in'
)

### figure 8E, kinase family example
# TK
# P42684_543, ABL2, Tyrosine-protein kinase ABL2
P42684_database_info <- tribble(
  ~ name, ~ start, ~ end,
  'protein', 1, 1182,
  'SH2', 173, 248,
  'SH3_1', 113, 159,
  'PK_Tyr_Ser-Thr', 288, 538,
  'F_actin_bind', 1083, 1182
)

P42684_result <- tibble(
  cleavage_site = c(543)
)

# example plot
P42684_example <- ggplot() +
  geom_rect(
    data = P42684_database_info,
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
      x = P42684_result$cleavage_site, 
      xend = P42684_result$cleavage_site, 
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
      'SH2' = color_2,
      'SH3_1' = color_3,
      'F_actin_bind' = color_4,
      'PK_Tyr_Ser-Thr' = color_5
    )
  ) +
  scale_color_manual(
    values = c(
      'protein' = 'black',
      'PK_Tyr_Ser-Thr' = 'transparent',
      'SH2' = 'transparent',
      'SH3_1' = 'transparent',
      'F_actin_bind' = 'transparent'
    )
  ) +
  theme_void()

ggsave(
  filename = 'figures/figure8/P42684_example.eps',
  height = 0.2, width = 2, units = 'in'
)

## TKL
# Q13546_577, RIPK1, Receptor-interacting serine/threonine-protein kinase 1
Q13546_database_info <- tribble(
  ~ name, ~ start, ~ end,
  'protein', 1, 671,
  'Death', 585, 666,
  'PK_Tyr_Ser-Thr', 21, 284
)

Q13546_result <- tibble(
  cleavage_site = c(577)
)

# example plot
Q13546_example <- ggplot() +
  geom_rect(
    data = Q13546_database_info,
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
      x = Q13546_result$cleavage_site, 
      xend = Q13546_result$cleavage_site, 
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
      'Death' = color_4,
      'PK_Tyr_Ser-Thr' = color_5
    )
  ) +
  scale_color_manual(
    values = c(
      'protein' = 'black',
      'PK_Tyr_Ser-Thr' = 'transparent',
      'Death' = 'transparent'
    )
  ) +
  theme_void()

ggsave(
  filename = 'figures/figure8/Q13546_example.eps',
  height = 0.2, width = 2, units = 'in'
)

# O43318_467, MAP3K7, Mitogen-activated protein kinase kinase kinase 7
O43318_database_info <- tribble(
  ~ name, ~ start, ~ end,
  'protein', 1, 606,
  'PK_Tyr_Ser-Thr', 36, 283
)

O43318_result <- tibble(
  cleavage_site = c(467)
)

# example plot
O43318_example <- ggplot() +
  geom_rect(
    data = O43318_database_info,
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
      x = O43318_result$cleavage_site, 
      xend = O43318_result$cleavage_site, 
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
      'PK_Tyr_Ser-Thr' = color_5
    )
  ) +
  scale_color_manual(
    values = c(
      'protein' = 'black',
      'PK_Tyr_Ser-Thr' = 'transparent'
    )
  ) +
  theme_void()

ggsave(
  filename = 'figures/figure8/O43318_example.eps',
  height = 0.2, width = 2, units = 'in'
)

# P10398_236, ARAF, Serine/threonine-protein kinase A-Raf
P10398_database_info <- tribble(
  ~ name, ~ start, ~ end,
  'protein', 1, 606,
  'C1_1', 99, 145,
  'RBD', 20, 89,
  'PK_Tyr_Ser-Thr', 311, 565
)

P10398_result <- tibble(
  cleavage_site = c(236)
)

# example plot
P10398_example <- ggplot() +
  geom_rect(
    data = P10398_database_info,
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
      x = P10398_result$cleavage_site, 
      xend = P10398_result$cleavage_site, 
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
      'C1_1' = color_3,
      'RBD' = color_4,
      'PK_Tyr_Ser-Thr' = color_5
    )
  ) +
  scale_color_manual(
    values = c(
      'protein' = 'black',
      'C1_1' = 'transparent',
      'RBD' = 'transparent',
      'PK_Tyr_Ser-Thr' = 'transparent'
    )
  ) +
  theme_void()

ggsave(
  filename = 'figures/figure8/P10398_example.eps',
  height = 0.2, width = 2, units = 'in'
)

### figure 8F-H, spliceosome complex example
# corum_id_8391, spliceosome E complex
# corum_id_8372, spliceosome A complex
# corum_id_1181, spliceosome C complex
# import result from corum bottom20
library(tidyverse)

spliceosome_complex_unique <- read_csv(
  'data_source/Nterm_WP_comparison/bottom20_comparison_corum.csv'
) |> 
  filter(
    Description %in% c(
      'corum_id_8391',
      'corum_id_8372',
      'corum_id_1181'
    )
  ) |> 
  mutate(
    spliceosome_complex = case_when(
      Description == 'corum_id_8391' ~ 'spliceosome E complex',
      Description == 'corum_id_8372' ~ 'spliceosome A complex',
      Description == 'corum_id_1181' ~ 'spliceosome C complex'
    )
  ) |> 
  select(
    spliceosome_complex,
    geneID
  ) |> 
  separate_rows(geneID, sep = '/') |> 
  add_count(geneID) |>
  filter(n == 1) |>
  select(-n)

spliceosome_complex_unique_half_life <- spliceosome_complex_unique |> 
  left_join(
    HEK_Nterm_WP_delta_half_life, by = join_by(geneID == UniProt_Accession)
  )

write_csv(
  spliceosome_complex_unique_half_life,
  file = 'data_source/Nterm_WP_comparison/spliceosome_complex_unique_half_life.csv'
)

### figure 8F, P38919_77 example plot
## P38919, EIF4A3, Eukaryotic initiation factor 4A-III
P38919_database_info <- tribble(
  ~ name, ~ start, ~ end,
  'protein', 1, 411,
  'DEAD', 63, 225,
  'Helicase_C', 264, 372
)

P38919_result <- tibble(
  cleavage_site = c(77)
)

# example plot
P38919_example <- ggplot() +
  geom_rect(
    data = P38919_database_info,
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
      x = P38919_result$cleavage_site, 
      xend = P38919_result$cleavage_site, 
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
      'DEAD' = color_1,
      'Helicase_C' = color_2
    )
  ) +
  scale_color_manual(
    values = c(
      'protein' = 'black',
      'DEAD' = 'transparent',
      'Helicase_C' = 'transparent'
    )
  ) +
  theme_void()

ggsave(
  filename = 'figures/figure8/P38919_example.eps',
  height = 0.2, width = 2, units = 'in'
)

## degradation profile
## Nterm
Nterm_P38919_77_deg_ratio <- HEK_Nterm_deg_ratio_adj |> 
  filter(Index == 'P38919_77')

Nterm_P38919_77_linear_model <- HEK_Nterm_linear_fitting |> 
  filter(Index == 'P38919_77') |> 
  pivot_wider(names_from = parameters, values_from = values)

Nterm_P38919_77_linear_model_adj <- HEK_Nterm_Kd_half_life_LaminB_Tcomplex |> 
  filter(Index == 'P38919_77')

time_point <- c(0, 3, 6, 9, 12, 24)
Nterm_P38919_77_time_point <- exp(Nterm_P38919_77_linear_model$lnA - Nterm_P38919_77_linear_model$Kd * time_point)
Nterm_P38919_77_time_point_adj <- exp(Nterm_P38919_77_linear_model$lnA - Nterm_P38919_77_linear_model_adj$Kd_adj * time_point)

adjusted_raito_P38919_77 <- Nterm_P38919_77_time_point_adj / Nterm_P38919_77_time_point
Nterm_P38919_77_deg_ratio_adj <- Nterm_P38919_77_deg_ratio |> 
  mutate(
    deg_ratio_adj_adj = deg_ratio_adj * adjusted_raito_P38919_77
  ) |> 
  mutate(
    rel_intensity = deg_ratio_adj_adj / max(deg_ratio_adj_adj)
  )

time_series <- seq(0, 24, length.out = 100)
Nterm_P38919_77_linear_fitting_adj <- exp(Nterm_P38919_77_linear_model$lnA - Nterm_P38919_77_linear_model_adj$Kd_adj * time_series)

Nterm_P38919_77_linear_fitting_line <- tibble(
  timepoint = time_series,
  deg_ratio = Nterm_P38919_77_linear_fitting_adj
)

## WP
WP_P38919_deg_ratio <- HEK_WP_deg_ratio_adj |> 
  filter(UniProt_Accession == 'P38919')

WP_P38919_linear_model <- HEK_WP_linear_fitting |>
  filter(UniProt_Accession == 'P38919') |> 
  pivot_wider(names_from = parameters, values_from = values)

WP_P38919_linear_model_adj <- HEK_WP_Kd_half_life_LaminB_Tcomplex |>
  filter(UniProt_Accession == 'P38919')

time_point <- c(0, 3, 6, 9, 12, 24)
WP_P38919_time_point <- exp(WP_P38919_linear_model$lnA - WP_P38919_linear_model$Kd * time_point)
WP_P38919_time_point_adj <- exp(WP_P38919_linear_model$lnA - WP_P38919_linear_model_adj$Kd_adj * time_point)

adjusted_raito_P38919 <- WP_P38919_time_point_adj / WP_P38919_time_point
WP_P38919_deg_ratio_adj <- WP_P38919_deg_ratio |> 
  mutate(
    deg_ratio_adj_adj = deg_ratio_adj * adjusted_raito
  ) |> 
  mutate(
    rel_intensity = deg_ratio_adj_adj / max(deg_ratio_adj_adj)
  )

time_series <- seq(0, 24, length.out = 100)
WP_P38919_linear_fitting_adj <- exp(WP_P38919_linear_model$lnA - WP_P38919_linear_model_adj$Kd_adj * time_series)

WP_P38919_linear_fitting_line <- tibble(
  timepoint = time_series,
  deg_ratio = WP_P38919_linear_fitting_adj
)

# plot
Nterm_WP_P38919_77_linear_example <- ggplot() +
  geom_point(
    data = Nterm_P38919_77_deg_ratio_adj,
    aes(
      x = timepoint,
      y = rel_intensity
    ),
    shape = 21, fill = color_1, color = 'transparent', size = 1.2
  ) +
  geom_point(
    data = WP_P38919_deg_ratio_adj,
    aes(
      x = timepoint,
      y = rel_intensity
    ),
    shape = 21, fill = color_2, color = 'transparent', size = 1.2
  ) +
  geom_line(
    data = Nterm_P38919_77_linear_fitting_line,
    aes(
      x = timepoint,
      y = deg_ratio
    )
  ) +
  geom_line(
    data = WP_P38919_linear_fitting_line,
    aes(
      x = timepoint,
      y = deg_ratio
    )
  ) +
  labs(x = NULL, y = NULL) +
  coord_cartesian(ylim = c(0, 1)) +
  theme(
    axis.text = element_text(size = 8, color = 'black', family = 'arial')
  )

ggsave(
  filename = 'figures/figure8/Nterm_WP_P38919_77_linear_example.eps',
  plot = Nterm_WP_P38919_77_linear_example,
  height = 1, width = 2, units = 'in'
)

### figure 8G, P11142_142 example plot
## P11142, HSPA8, Heat shock cognate 71 kDa protein
P11142_database_info <- tribble(
  ~ name, ~ start, ~ end,
  'protein', 1, 646,
  'HSP70', 6, 612
)

P11142_result <- tibble(
  cleavage_site = c(142)
)

# example plot
P11142_example <- ggplot() +
  geom_rect(
    data = P11142_database_info,
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
      x = P11142_result$cleavage_site, 
      xend = P11142_result$cleavage_site, 
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
  filename = 'figures/figure8/P11142_example.eps',
  height = 0.2, width = 2, units = 'in'
)

## degradation profile
## Nterm
Nterm_P11142_142_deg_ratio <- HEK_Nterm_deg_ratio_adj |> 
  filter(Index == 'P11142_142')

Nterm_P11142_142_linear_model <- HEK_Nterm_linear_fitting |> 
  filter(Index == 'P11142_142') |> 
  pivot_wider(names_from = parameters, values_from = values)

Nterm_P11142_142_linear_model_adj <- HEK_Nterm_Kd_half_life_LaminB_Tcomplex |> 
  filter(Index == 'P11142_142')

time_point <- c(0, 3, 6, 9, 12, 24)
Nterm_P11142_142_time_point <- exp(Nterm_P11142_142_linear_model$lnA - Nterm_P11142_142_linear_model$Kd * time_point)
Nterm_P11142_142_time_point_adj <- exp(Nterm_P11142_142_linear_model$lnA - Nterm_P11142_142_linear_model_adj$Kd_adj * time_point)

adjusted_raito_P11142_142 <- Nterm_P11142_142_time_point_adj / Nterm_P11142_142_time_point
Nterm_P11142_142_deg_ratio_adj <- Nterm_P11142_142_deg_ratio |> 
  mutate(
    deg_ratio_adj_adj = deg_ratio_adj * adjusted_raito_P11142_142
  ) |> 
  mutate(
    rel_intensity = deg_ratio_adj_adj / max(deg_ratio_adj_adj)
  )

time_series <- seq(0, 24, length.out = 100)
Nterm_P11142_142_linear_fitting_adj <- exp(Nterm_P11142_142_linear_model$lnA - Nterm_P11142_142_linear_model_adj$Kd_adj * time_series)

Nterm_P11142_142_linear_fitting_line <- tibble(
  timepoint = time_series,
  deg_ratio = Nterm_P11142_142_linear_fitting_adj
)

## WP
WP_P11142_deg_ratio <- HEK_WP_deg_ratio_adj |> 
  filter(UniProt_Accession == 'P11142')

WP_P11142_linear_model <- HEK_WP_linear_fitting |>
  filter(UniProt_Accession == 'P11142') |> 
  pivot_wider(names_from = parameters, values_from = values)

WP_P11142_linear_model_adj <- HEK_WP_Kd_half_life_LaminB_Tcomplex |>
  filter(UniProt_Accession == 'P11142')

time_point <- c(0, 3, 6, 9, 12, 24)
WP_P11142_time_point <- exp(WP_P11142_linear_model$lnA - WP_P11142_linear_model$Kd * time_point)
WP_P11142_time_point_adj <- exp(WP_P11142_linear_model$lnA - WP_P11142_linear_model_adj$Kd_adj * time_point)

adjusted_raito_P11142 <- WP_P11142_time_point_adj / WP_P11142_time_point
WP_P11142_deg_ratio_adj <- WP_P11142_deg_ratio |> 
  mutate(
    deg_ratio_adj_adj = deg_ratio_adj * adjusted_raito_P11142
  ) |> 
  mutate(
    rel_intensity = deg_ratio_adj_adj / max(deg_ratio_adj_adj)
  )

time_series <- seq(0, 24, length.out = 100)
WP_P11142_linear_fitting_adj <- exp(WP_P11142_linear_model$lnA - WP_P11142_linear_model_adj$Kd_adj * time_series)

WP_P11142_linear_fitting_line <- tibble(
  timepoint = time_series,
  deg_ratio = WP_P11142_linear_fitting_adj
)

# plot
Nterm_WP_P11142_142_linear_example <- ggplot() +
  geom_point(
    data = Nterm_P11142_142_deg_ratio_adj,
    aes(
      x = timepoint,
      y = rel_intensity
    ),
    shape = 21, fill = color_1, color = 'transparent', size = 1.2
  ) +
  geom_point(
    data = WP_P11142_deg_ratio_adj,
    aes(
      x = timepoint,
      y = rel_intensity
    ),
    shape = 21, fill = color_2, color = 'transparent', size = 1.2
  ) +
  geom_line(
    data = Nterm_P11142_142_linear_fitting_line,
    aes(
      x = timepoint,
      y = deg_ratio
    )
  ) +
  geom_line(
    data = WP_P11142_linear_fitting_line,
    aes(
      x = timepoint,
      y = deg_ratio
    )
  ) +
  labs(x = NULL, y = NULL) +
  coord_cartesian(ylim = c(0, 1)) +
  theme(
    axis.text = element_text(size = 8, color = 'black', family = 'arial')
  )

ggsave(
  filename = 'figures/figure8/Nterm_WP_P11142_142_linear_example.eps',
  plot = Nterm_WP_P11142_142_linear_example,
  height = 1, width = 2, units = 'in'
)

### figure 8G,  Q14103_99 example plot
## Q14103, HNRNPD, Heterogeneous nuclear ribonucleoprotein D0
Q14103_database_info <- tribble(
  ~ name, ~ start, ~ end,
  'protein', 1, 355,
  'RRM_1', 100, 167,
  'RRM_1', 184, 243,
  'CBFNT', 1, 78
)

Q14103_result <- tibble(
  cleavage_site = c(99)
)

# example plot
Q14103_example <- ggplot() +
  geom_rect(
    data = Q14103_database_info,
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
      x = Q14103_result$cleavage_site, 
      xend = Q14103_result$cleavage_site, 
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
      'RRM_1' = color_4,
      'CBFNT' = color_5
    )
  ) +
  scale_color_manual(
    values = c(
      'protein' = 'black',
      'RRM_1' = 'transparent',
      'CBFNT' = 'transparent'
    )
  ) +
  theme_void()

ggsave(
  filename = 'figures/figure8/Q14103_example.eps',
  height = 0.2, width = 2, units = 'in'
)

## degradation profile
## Nterm
Nterm_Q14103_99_deg_ratio <- HEK_Nterm_deg_ratio_adj |> 
  filter(Index == 'Q14103_99')

Nterm_Q14103_99_linear_model <- HEK_Nterm_linear_fitting |> 
  filter(Index == 'Q14103_99') |> 
  pivot_wider(names_from = parameters, values_from = values)

Nterm_Q14103_99_linear_model_adj <- HEK_Nterm_Kd_half_life_LaminB_Tcomplex |> 
  filter(Index == 'Q14103_99')

time_point <- c(0, 3, 6, 9, 12, 24)
Nterm_Q14103_99_time_point <- exp(Nterm_Q14103_99_linear_model$lnA - Nterm_Q14103_99_linear_model$Kd * time_point)
Nterm_Q14103_99_time_point_adj <- exp(Nterm_Q14103_99_linear_model$lnA - Nterm_Q14103_99_linear_model_adj$Kd_adj * time_point)

adjusted_raito_Q14103_99 <- Nterm_Q14103_99_time_point_adj / Nterm_Q14103_99_time_point
Nterm_Q14103_99_deg_ratio_adj <- Nterm_Q14103_99_deg_ratio |> 
  mutate(
    deg_ratio_adj_adj = deg_ratio_adj * adjusted_raito_Q14103_99
  ) |> 
  mutate(
    rel_intensity = deg_ratio_adj_adj / max(deg_ratio_adj_adj)
  )

time_series <- seq(0, 24, length.out = 100)
Nterm_Q14103_99_linear_fitting_adj <- exp(Nterm_Q14103_99_linear_model$lnA - Nterm_Q14103_99_linear_model_adj$Kd_adj * time_series)

Nterm_Q14103_99_linear_fitting_line <- tibble(
  timepoint = time_series,
  deg_ratio = Nterm_Q14103_99_linear_fitting_adj
)

## WP
WP_Q14103_deg_ratio <- HEK_WP_deg_ratio_adj |> 
  filter(UniProt_Accession == 'Q14103')

WP_Q14103_linear_model <- HEK_WP_linear_fitting |>
  filter(UniProt_Accession == 'Q14103') |> 
  pivot_wider(names_from = parameters, values_from = values)

WP_Q14103_linear_model_adj <- HEK_WP_Kd_half_life_LaminB_Tcomplex |>
  filter(UniProt_Accession == 'Q14103')

time_point <- c(0, 3, 6, 9, 12, 24)
WP_Q14103_time_point <- exp(WP_Q14103_linear_model$lnA - WP_Q14103_linear_model$Kd * time_point)
WP_Q14103_time_point_adj <- exp(WP_Q14103_linear_model$lnA - WP_Q14103_linear_model_adj$Kd_adj * time_point)

adjusted_raito_Q14103 <- WP_Q14103_time_point_adj / WP_Q14103_time_point
WP_Q14103_deg_ratio_adj <- WP_Q14103_deg_ratio |> 
  mutate(
    deg_ratio_adj_adj = deg_ratio_adj * adjusted_raito_Q14103
  ) |> 
  mutate(
    rel_intensity = deg_ratio_adj_adj / max(deg_ratio_adj_adj)
  )

time_series <- seq(0, 24, length.out = 100)
WP_Q14103_linear_fitting_adj <- exp(WP_Q14103_linear_model$lnA - WP_Q14103_linear_model_adj$Kd_adj * time_series)

WP_Q14103_linear_fitting_line <- tibble(
  timepoint = time_series,
  deg_ratio = WP_Q14103_linear_fitting_adj
)

# plot
Nterm_WP_Q14103_99_linear_example <- ggplot() +
  geom_point(
    data = Nterm_Q14103_99_deg_ratio_adj,
    aes(
      x = timepoint,
      y = rel_intensity
    ),
    shape = 21, fill = color_1, color = 'transparent', size = 1.2
  ) +
  geom_point(
    data = WP_Q14103_deg_ratio_adj,
    aes(
      x = timepoint,
      y = rel_intensity
    ),
    shape = 21, fill = color_2, color = 'transparent', size = 1.2
  ) +
  geom_line(
    data = Nterm_Q14103_99_linear_fitting_line,
    aes(
      x = timepoint,
      y = deg_ratio
    )
  ) +
  geom_line(
    data = WP_Q14103_linear_fitting_line,
    aes(
      x = timepoint,
      y = deg_ratio
    )
  ) +
  labs(x = NULL, y = NULL) +
  coord_cartesian(ylim = c(0, 1)) +
  theme(
    axis.text = element_text(size = 8, color = 'black', family = 'arial')
  )

ggsave(
  filename = 'figures/figure8/Nterm_WP_Q14103_99_linear_example.eps',
  plot = Nterm_WP_Q14103_99_linear_example,
  height = 1, width = 2, units = 'in'
)

### figure 8H, top20 and bottom20 structure analysis
## import result from StructureMap
library(tidyverse)

top20_bottom20_protein_structure <- read_csv(
  'data_source/Nterm_WP_comparison/top20_bottom20_protein_structure.csv'
)

## secondary structure
top20_secondary_structure <- top20_bottom20_protein_structure |> 
  filter(top20 == 1) |> 
  count(structure_group) |> 
  mutate(
    percentage = n/sum(n),
    category = 'top20'
  )

bottom20_secondary_structure <- top20_bottom20_protein_structure |> 
  filter(bottom20 == 1) |> 
  count(structure_group) |> 
  mutate(
    percentage = n/sum(n),
    category = 'bottom20'
  )

# combine result for top20 and bottom20
top20_bottom20_secondary_structure_comb <- bind_rows(
  top20_secondary_structure,
  bottom20_secondary_structure
)

# ANOVA test/Turkey HSD
top20_bottom20_secondary_structure_anova_test <- aov(
  percentage ~ category, data = top20_bottom20_secondary_structure_comb
)

summary(top20_bottom20_secondary_structure_anova_test)

# bar plot
barplot_top20_bottom20_secondary_structure_comb <- top20_bottom20_secondary_structure_comb |> 
  ggplot() +
  geom_bar(
    aes(
      x = factor(category, levels = c('top20', 'bottom20')),
      y = percentage,
      fill = structure_group
    ), 
    stat = 'identity',
    position = 'stack'
  ) +
  labs(x = '', y = '') +
  scale_fill_manual(
    name = '',
    values = c(
      'BEND' = color_1,
      'HELX' = color_2,
      'STRN' = color_3,
      'TURN' = color_4,
      'unstructured' = 'gray70'
    )
  ) +
  ylim(0, 4/3) +
  xlim('1', '2', 'top20', 'bottom20') +
  coord_polar(theta = "y") +
  theme(
    axis.text.x = element_text(size = 8, color = 'black', family = 'arial'),
    axis.text.y = element_text(size = 8, color = 'black', family = 'arial'),
    legend.text = element_text(size = 8, color = 'black', family = 'arial')
  )

ggsave(
  filename = 'figures/figure8/barplot_top20_bottom20_secondary_structure_comb.eps',
  plot = barplot_top20_bottom20_secondary_structure_comb,
  height = 2, width = 5, units = 'in'
)

## IDR
top20_IDR <- Nterm_WP_delta_half_life_alphafold_N_terminus |> 
  as_tibble() |> 
  filter(top20 == 1) |> 
  count(IDR) |> 
  mutate(
    percentage = n/sum(n),
    category = 'top20'
  )

bottom20_IDR <- Nterm_WP_delta_half_life_alphafold_N_terminus |> 
  as_tibble() |> 
  filter(bottom20 == 1) |> 
  count(IDR) |> 
  mutate(
    percentage = n/sum(n),
    category = 'bottom20'
  )

# combine result for top20 and bottom20
top20_bottom20_IDR_comb <- bind_rows(
  top20_IDR,
  bottom20_IDR
)

# ANOVA test/Turkey HSD
top20_bottom20_IDR_anova_test <- aov(
  percentage ~ category, data = top20_bottom20_IDR_comb
)

summary(top20_bottom20_IDR_anova_test)

# bar plot
barplot_top20_bottom20_IDR_comb <- top20_bottom20_IDR_comb |> 
  mutate(
    IDR = as.character(IDR)
  ) |> 
  ggplot() +
  geom_bar(
    aes(
      x = factor(category, levels = c('top20', 'bottom20')),
      y = percentage,
      fill = IDR
    ), 
    stat = 'identity',
    position = 'stack'
  ) +
  labs(x = '', y = '') +
  scale_fill_manual(
    name = '',
    values = c(
      '1' = color_1,
      '0' = color_2
    )
  ) +
  theme(
    axis.text.x = element_text(size = 8, color = 'black', angle = 30, hjust = 1),
    axis.text.y = element_text(size = 8, color = 'black'),
    legend.text = element_text(size = 8, color = 'black')
  )

ggsave(
  filename = 'figures/figure8/barplot_top20_bottom20_IDR_comb.eps',
  plot = barplot_top20_bottom20_IDR_comb,
  height = 2, width = 2, units = 'in'
)

### figure 8I, Nucleolin example
# P06748, NPM1, Nucleophosmin
P06748_database_info <- tribble(
  ~ name, ~ start, ~ end,
  'protein', 1, 294,
  'Nucleoplasmin', 18, 117,
  'NPM1-C', 245, 293
)

P06748_result <- tibble(
  cleavage_site = c(
    46, 89, 90, 91
  )
)

# example plot
P06748_example <- ggplot() +
  geom_rect(
    data = P06748_database_info,
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
      x = P06748_result$cleavage_site, 
      xend = P06748_result$cleavage_site, 
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
      'Nucleoplasmin' = color_1,
      'NPM1-C' = color_2
    )
  ) +
  scale_color_manual(
    values = c(
      'protein' = 'black',
      'Nucleoplasmin' = 'transparent',
      'NPM1-C' = 'transparent'
    )
  ) +
  theme_void()

ggsave(
  filename = 'figures/figure8/P06748_example.eps',
  height = 0.2, width = 2, units = 'in'
)

# ## P06748 degradation profile
# ## Nterm
# # P06748_46
# Nterm_P06748_46_deg_ratio <- HEK_Nterm_deg_ratio_adj |> 
#   filter(Index == 'P06748_46')
# 
# Nterm_P06748_46_linear_model <- HEK_Nterm_linear_fitting |> 
#   filter(Index == 'P06748_46') |> 
#   pivot_wider(names_from = parameters, values_from = values)
# 
# Nterm_P06748_46_linear_model_adj <- HEK_Nterm_Kd_half_life_LaminB_Tcomplex |> 
#   filter(Index == 'P06748_46')
# 
# time_point <- c(0, 3, 6, 9, 12, 24)
# Nterm_P06748_46_time_point <- exp(Nterm_P06748_46_linear_model$lnA - Nterm_P06748_46_linear_model$Kd * time_point)
# Nterm_P06748_46_time_point_adj <- exp(Nterm_P06748_46_linear_model$lnA - Nterm_P06748_46_linear_model_adj$Kd_adj * time_point)
# 
# adjusted_raito_P06748_46 <- Nterm_P06748_46_time_point_adj / Nterm_P06748_46_time_point
# Nterm_P06748_46_deg_ratio_adj <- Nterm_P06748_46_deg_ratio |> 
#   mutate(
#     deg_ratio_adj_adj = deg_ratio_adj * adjusted_raito_P06748_46
#   ) |> 
#   mutate(
#     rel_intensity = deg_ratio_adj_adj / max(deg_ratio_adj_adj)
#   )
# 
# time_series <- seq(0, 24, length.out = 100)
# Nterm_P06748_46_linear_fitting_adj <- exp(Nterm_P06748_46_linear_model$lnA - Nterm_P06748_46_linear_model_adj$Kd_adj * time_series)
# 
# Nterm_P06748_46_linear_fitting_line <- tibble(
#   timepoint = time_series,
#   deg_ratio = Nterm_P06748_46_linear_fitting_adj
# )
# 
# # P06748_89
# Nterm_P06748_89_deg_ratio <- HEK_Nterm_deg_ratio_adj |> 
#   filter(Index == 'P06748_89')
# 
# Nterm_P06748_89_linear_model <- HEK_Nterm_linear_fitting |> 
#   filter(Index == 'P06748_89') |> 
#   pivot_wider(names_from = parameters, values_from = values)
# 
# Nterm_P06748_89_linear_model_adj <- HEK_Nterm_Kd_half_life_LaminB_Tcomplex |> 
#   filter(Index == 'P06748_89')
# 
# time_point <- c(0, 3, 6, 9, 12, 24)
# Nterm_P06748_89_time_point <- exp(Nterm_P06748_89_linear_model$lnA - Nterm_P06748_89_linear_model$Kd * time_point)
# Nterm_P06748_89_time_point_adj <- exp(Nterm_P06748_89_linear_model$lnA - Nterm_P06748_89_linear_model_adj$Kd_adj * time_point)
# 
# adjusted_raito_P06748_89 <- Nterm_P06748_89_time_point_adj / Nterm_P06748_89_time_point
# Nterm_P06748_89_deg_ratio_adj <- Nterm_P06748_89_deg_ratio |> 
#   mutate(
#     deg_ratio_adj_adj = deg_ratio_adj * adjusted_raito_P06748_89
#   ) |> 
#   mutate(
#     rel_intensity = deg_ratio_adj_adj / max(deg_ratio_adj_adj)
#   )
# 
# time_series <- seq(0, 24, length.out = 100)
# Nterm_P06748_89_linear_fitting_adj <- exp(Nterm_P06748_89_linear_model$lnA - Nterm_P06748_89_linear_model_adj$Kd_adj * time_series)
# 
# Nterm_P06748_89_linear_fitting_line <- tibble(
#   timepoint = time_series,
#   deg_ratio = Nterm_P06748_89_linear_fitting_adj
# )
# 
# # P06748_90
# Nterm_P06748_90_deg_ratio <- HEK_Nterm_deg_ratio_adj |> 
#   filter(Index == 'P06748_90')
# 
# Nterm_P06748_90_linear_model <- HEK_Nterm_linear_fitting |> 
#   filter(Index == 'P06748_90') |> 
#   pivot_wider(names_from = parameters, values_from = values)
# 
# Nterm_P06748_90_linear_model_adj <- HEK_Nterm_Kd_half_life_LaminB_Tcomplex |> 
#   filter(Index == 'P06748_90')
# 
# time_point <- c(0, 3, 6, 9, 12, 24)
# Nterm_P06748_90_time_point <- exp(Nterm_P06748_90_linear_model$lnA - Nterm_P06748_90_linear_model$Kd * time_point)
# Nterm_P06748_90_time_point_adj <- exp(Nterm_P06748_90_linear_model$lnA - Nterm_P06748_90_linear_model_adj$Kd_adj * time_point)
# 
# adjusted_raito_P06748_90 <- Nterm_P06748_90_time_point_adj / Nterm_P06748_90_time_point
# Nterm_P06748_90_deg_ratio_adj <- Nterm_P06748_90_deg_ratio |> 
#   mutate(
#     deg_ratio_adj_adj = deg_ratio_adj * adjusted_raito_P06748_90
#   ) |> 
#   mutate(
#     rel_intensity = deg_ratio_adj_adj / max(deg_ratio_adj_adj)
#   )
# 
# time_series <- seq(0, 24, length.out = 100)
# Nterm_P06748_90_linear_fitting_adj <- exp(Nterm_P06748_90_linear_model$lnA - Nterm_P06748_90_linear_model_adj$Kd_adj * time_series)
# 
# Nterm_P06748_90_linear_fitting_line <- tibble(
#   timepoint = time_series,
#   deg_ratio = Nterm_P06748_90_linear_fitting_adj
# )
# 
# # P06748_91
# Nterm_P06748_91_deg_ratio <- HEK_Nterm_deg_ratio_adj |> 
#   filter(Index == 'P06748_91')
# 
# Nterm_P06748_91_linear_model <- HEK_Nterm_linear_fitting |> 
#   filter(Index == 'P06748_91') |> 
#   pivot_wider(names_from = parameters, values_from = values)
# 
# Nterm_P06748_91_linear_model_adj <- HEK_Nterm_Kd_half_life_LaminB_Tcomplex |> 
#   filter(Index == 'P06748_91')
# 
# time_point <- c(0, 3, 6, 9, 12, 24)
# Nterm_P06748_91_time_point <- exp(Nterm_P06748_91_linear_model$lnA - Nterm_P06748_91_linear_model$Kd * time_point)
# Nterm_P06748_91_time_point_adj <- exp(Nterm_P06748_91_linear_model$lnA - Nterm_P06748_91_linear_model_adj$Kd_adj * time_point)
# 
# adjusted_raito_P06748_91 <- Nterm_P06748_91_time_point_adj / Nterm_P06748_91_time_point
# Nterm_P06748_91_deg_ratio_adj <- Nterm_P06748_91_deg_ratio |> 
#   mutate(
#     deg_ratio_adj_adj = deg_ratio_adj * adjusted_raito_P06748_91
#   ) |> 
#   mutate(
#     rel_intensity = deg_ratio_adj_adj / max(deg_ratio_adj_adj)
#   )
# 
# time_series <- seq(0, 24, length.out = 100)
# Nterm_P06748_91_linear_fitting_adj <- exp(Nterm_P06748_91_linear_model$lnA - Nterm_P06748_91_linear_model_adj$Kd_adj * time_series)
# 
# Nterm_P06748_91_linear_fitting_line <- tibble(
#   timepoint = time_series,
#   deg_ratio = Nterm_P06748_91_linear_fitting_adj
# )
# 
# ## WP
# WP_P06748_deg_ratio <- HEK_WP_deg_ratio_adj |> 
#   filter(UniProt_Accession == 'P06748')
# 
# WP_P06748_linear_model <- HEK_WP_linear_fitting |>
#   filter(UniProt_Accession == 'P06748') |> 
#   pivot_wider(names_from = parameters, values_from = values)
# 
# WP_P06748_linear_model_adj <- HEK_WP_Kd_half_life_LaminB_Tcomplex |>
#   filter(UniProt_Accession == 'P06748')
# 
# time_point <- c(0, 3, 6, 9, 12, 24)
# WP_P06748_time_point <- exp(WP_P06748_linear_model$lnA - WP_P06748_linear_model$Kd * time_point)
# WP_P06748_time_point_adj <- exp(WP_P06748_linear_model$lnA - WP_P06748_linear_model_adj$Kd_adj * time_point)
# 
# adjusted_raito <- WP_P06748_time_point_adj / WP_P06748_time_point
# WP_P06748_deg_ratio_adj <- WP_P06748_deg_ratio |> 
#   mutate(
#     deg_ratio_adj_adj = deg_ratio_adj * adjusted_raito
#   ) |> 
#   mutate(
#     rel_intensity = deg_ratio_adj_adj / max(deg_ratio_adj_adj)
#   )
# 
# time_series <- seq(0, 24, length.out = 100)
# WP_P06748_linear_fitting_adj <- exp(WP_P06748_linear_model$lnA - WP_P06748_linear_model_adj$Kd_adj * time_series)
# 
# WP_P06748_linear_fitting_line <- tibble(
#   timepoint = time_series,
#   deg_ratio = WP_P06748_linear_fitting_adj
# )
# 
# # plot
# Nterm_WP_P06748_linear_example <- ggplot() +
#   geom_point(
#     data = Nterm_P06748_46_deg_ratio_adj,
#     aes(
#       x = timepoint,
#       y = rel_intensity
#     ),
#     shape = 21, fill = color_1, color = 'transparent', size = 1.2
#   ) +
#   geom_point(
#     data = Nterm_P06748_89_deg_ratio_adj,
#     aes(
#       x = timepoint,
#       y = rel_intensity
#     ),
#     shape = 21, fill = color_3, color = 'transparent', size = 1.2
#   ) +
#   geom_point(
#     data = Nterm_P06748_90_deg_ratio_adj,
#     aes(
#       x = timepoint,
#       y = rel_intensity
#     ),
#     shape = 21, fill = color_4, color = 'transparent', size = 1.2
#   ) +
#   geom_point(
#     data = Nterm_P06748_91_deg_ratio_adj,
#     aes(
#       x = timepoint,
#       y = rel_intensity
#     ),
#     shape = 21, fill = color_5, color = 'transparent', size = 1.2
#   ) +
#   geom_point(
#     data = WP_P06748_deg_ratio_adj,
#     aes(
#       x = timepoint,
#       y = rel_intensity
#     ),
#     shape = 21, fill = color_2, color = 'transparent', size = 1.2
#   ) +
#   geom_line(
#     data = Nterm_P06748_46_linear_fitting_line,
#     aes(
#       x = timepoint,
#       y = deg_ratio
#     )
#   ) +
#   geom_line(
#     data = Nterm_P06748_89_linear_fitting_line,
#     aes(
#       x = timepoint,
#       y = deg_ratio
#     )
#   ) +
#   geom_line(
#     data = Nterm_P06748_90_linear_fitting_line,
#     aes(
#       x = timepoint,
#       y = deg_ratio
#     )
#   ) +
#   geom_line(
#     data = Nterm_P06748_91_linear_fitting_line,
#     aes(
#       x = timepoint,
#       y = deg_ratio
#     )
#   ) +
#   geom_line(
#     data = WP_P06748_linear_fitting_line,
#     aes(
#       x = timepoint,
#       y = deg_ratio
#     )
#   ) +
#   labs(x = NULL, y = NULL) +
#   coord_cartesian(ylim = c(0, 1)) +
#   theme(
#     axis.text = element_text(size = 8, color = 'black', family = 'arial')
#   )
# 
# ggsave(
#   filename = 'figures/figure8/Nterm_WP_P06748_linear_example.eps',
#   plot = Nterm_WP_P06748_linear_example,
#   height = 1, width = 2, units = 'in'
# )

# ### figure 8C, spliceosome, proteasome and ribosome
# ## GO analysis of overlap proteins
# library(clusterProfiler)
# library(org.Hs.eg.db)
# 
# Nterm_WP_overlap_GO <- enrichGO(
#   gene = HEK_Nterm_WP_delta_half_life |> distinct(UniProt_Accession) |> pull(),
#   OrgDb = org.Hs.eg.db,
#   keyType = 'UNIPROT',
#   ont = 'ALL',
#   pvalueCutoff = 1,
#   qvalueCutoff = 1
# )
# 
# write_csv(
#   Nterm_WP_overlap_GO@result,
#   file = 'data_source/Nterm_WP_comparison/Nterm_WP_overlap_GO.csv'
# )
# 
# ## import GO enrichment analysis result
# Nterm_WP_overlap_GO <- read_csv(
#   'data_source/Nterm_WP_comparison/Nterm_WP_overlap_GO.csv'
# )
# 
# ## ribosome
# # GO:0015935, small ribosomal subunit
# # GO:0022625, cytosolic large ribosomal subunit
# # GO:0000313, organellar ribosome
# ribosome_related_protein_list <- Nterm_WP_overlap_GO |> 
#   filter(
#     Description %in% c(
#       'small ribosomal subunit',
#       'cytosolic large ribosomal subunit',
#       'organellar ribosome'
#     )
#   ) |> 
#   dplyr::select(geneID) |> 
#   separate_rows(geneID, sep = '/') |> 
#   distinct() |> 
#   pull()
# 
# ## proteasome
# # GO:0000502, proteasome complex
# # GO:0022624, proteasome accessory complex
# # GO:0005838, proteasome regulatory particle
# proteasome_related_protein_list <- Nterm_WP_overlap_GO |> 
#   filter(
#     Description %in% c(
#       'proteasome complex',
#       'proteasome accessory complex', 
#       'proteasome regulatory particle'
#     )
#   ) |> 
#   dplyr::select(geneID) |> 
#   separate_rows(geneID, sep = '/') |> 
#   distinct() |> 
#   pull()
# 
# ## spliceosome
# # GO:0097525, spliceosomal snRNP complex
# # GO:0071011, precatalytic spliceosome
# # GO:0071006, U2-type catalytic step 1 spliceosome
# # GO:0071007, U2-type catalytic step 2 spliceosome
# # GO:0071010, prespliceosome
# spliceosome_related_protein_list <- Nterm_WP_overlap_GO |> 
#   filter(
#     Description %in% c(
#       'spliceosomal snRNP complex',
#       'precatalytic spliceosome',
#       'U2-type catalytic step 1 spliceosome',
#       'U2-type catalytic step 2 spliceosome',
#       'prespliceosome'
#     )
#   ) |> 
#   dplyr::select(geneID) |> 
#   separate_rows(geneID, sep = '/') |> 
#   distinct() |> 
#   pull()
# 
# # combine protein from different organelle
# Nterm_linear_fit_protein_list <- HEK_Nterm_curve_fitting_combined |> 
#   filter(model == 'linear fitting') |> 
#   distinct(Index) |> 
#   pull()
# 
# WP_linear_fit_protein_list <- HEK_WP_curve_fitting_combined |> 
#   filter(model == 'linear fitting') |> 
#   distinct(UniProt_Accession) |> 
#   pull()
# 
# ribosome_proteasome_spliceosome_comb <- bind_rows(
#   HEK_Nterm_WP_delta_half_life |> 
#     filter(UniProt_Accession %in% ribosome_related_protein_list) |> 
#     mutate(GO_category = 'ribosome'),
#   
#   HEK_Nterm_WP_delta_half_life |> 
#     filter(UniProt_Accession %in% proteasome_related_protein_list) |> 
#     mutate(GO_category = 'proteasome'),
#   
#   HEK_Nterm_WP_delta_half_life |> 
#     filter(UniProt_Accession %in% spliceosome_related_protein_list) |> 
#     mutate(GO_category = 'spliceosome')
# ) |> 
#   mutate(
#     Nterm_curve_fitting = ifelse(
#       Index %in% Nterm_linear_fit_protein_list,
#       'linear fitting',
#       'nonlinear fitting'
#     ),
#     WP_curve_fitting = ifelse(
#       UniProt_Accession %in% WP_linear_fit_protein_list,
#       'linear fitting',
#       'nonlinear fitting'
#     )
#   )
# 
# write_csv(
#   ribosome_proteasome_spliceosome_comb,
#   file = 'data_source/Nterm_WP_comparison/ribosome_proteasome_spliceosome_comb.csv'
# )
# 
# # Wilcoxon rank-sum test
# library(rstatix)
# 
# ribosome_proteasome_spliceosome_comb_wilcoxon_test <- ribosome_proteasome_spliceosome_comb |> 
#   wilcox_test(delta_half_life ~ GO_category)
# 
# # point range plot
# library(ggpubr)
# 
# point_range_plot_ribosome_proteasome_spliceosome_comb <- ribosome_proteasome_spliceosome_comb |> 
#   ggplot() +
#   geom_point(
#     aes(
#       x = GO_category,
#       y = delta_half_life
#     ),
#     position = position_jitter(width = 0.3),
#     color = 'black',
#     alpha = 0.3,
#     size = 0.5
#   ) +
#   stat_summary(
#     aes(
#       x = GO_category,
#       y = delta_half_life,
#       color = GO_category
#     ),
#     fun.data = 'mean_cl_boot', linewidth = 0.2, size = 0.5, show.legend = FALSE
#   ) +
#   labs(x = '', y = '') +
#   stat_pvalue_manual(
#     data = ribosome_proteasome_spliceosome_comb_wilcoxon_test,
#     label = 'p.adj.signif',
#     y.position = c(70, 100, 130),
#     tip.length = 0,
#     label.size = 6
#   ) +
#   scale_color_manual(
#     values = c(
#       'proteasome' = color_1,
#       'ribosome' = color_2,
#       'spliceosome' = color_3
#     )
#   ) +
#   theme(
#     panel.grid.major = element_line(color = "gray", linewidth = 0.2),
#     panel.grid.minor = element_line(color = "gray", linewidth = 0.1),
#     axis.text.x = element_text(color = 'black', size = 8, family = 'arial', angle = 30, hjust = 1),
#     axis.text.y = element_text(color = 'black', size = 8, family = 'arial')
#   )
# 
# ggsave(
#   filename = 'figures/figure8/point_range_plot_ribosome_proteasome_spliceosome_comb.eps',
#   device = cairo_ps,
#   plot = point_range_plot_ribosome_proteasome_spliceosome_comb,
#   height = 2, width = 1.8, units = 'in',
#   fallback_resolution = 1200
# )
# 
# ### figure 8D, ribosome example
# # P63244, RACK1, 188, 226, 229
# 
# 
# # P62899, RPL31
