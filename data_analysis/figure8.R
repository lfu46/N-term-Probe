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
        'protein folding',
        'endoplasmic reticulum lumen',
        'vesicle lumen',
        'secretory granule lumen',
        'extracellular organelle',
        'proteasome complex',
        'proteasome regulatory particle, base subcomplex'
      )
    ) |> 
    mutate(
      category = 'top20'
    ),
  
  top20_comparison_KEGG |> 
    filter(
      Description %in% c(
        'Proteasome',
        'Protein processing in endoplasmic reticulum'
      )
    ) |> 
    mutate(
      category = 'top20'
    ),
  
  bottom20_comparison_GO |> 
    filter(
      Description %in% c(
        'nucleotide binding',
        'anion binding',
        'ribonucleotide binding',
        'carbohydrate derivative binding',
        'ATP binding',
        'helicase activity',
        'ATP hydrolysis activity',
        'pyrophosphatase activity',
        'structural constituent of ribosome'
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

### figure 8B, protein complex analysis for top20% and bottom 20% proteins
## import corum enrichment results
top20_comparison_corum <- read_csv(
  'data_source/Nterm_WP_comparison/top20_comparison_corum.csv'
)

bottom20_comparison_corum <- read_csv(
  'data_source/Nterm_WP_comparison/bottom20_comparison_corum.csv'
)

## bar plot
# top20
barplot_top20_corum <- top20_comparison_corum |> 
  filter(
    Description %in% c(
      'corum_id_181',
      'corum_id_8372',
      'corum_id_193',
      'corum_id_5615',
      'corum_id_8391'
    )
  ) |> 
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
  filename = 'figures/figure8/barplot_top20_corum.eps',
  plot = barplot_top20_corum,
  height = 2, width = 2.5, units = 'in'
)

# bottom20
barplot_bottom20_corum <- bottom20_comparison_corum |> 
  filter(
    Description %in% c(
      'corum_id_3055',
      'corum_id_306',
      'corum_id_305',
      'corum_id_3040',
      'corum_id_7627'
    )
  ) |> 
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
  filename = 'figures/figure8/barplot_bottom20_corum.eps',
  plot = barplot_bottom20_corum,
  height = 2, width = 2.5, units = 'in'
)

### figure 8C, spliceosome, proteasome and ribosome
## GO analysis of overlap proteins
library(clusterProfiler)
library(org.Hs.eg.db)

Nterm_WP_overlap_GO <- enrichGO(
  gene = HEK_Nterm_WP_delta_half_life |> distinct(UniProt_Accession) |> pull(),
  OrgDb = org.Hs.eg.db,
  keyType = 'UNIPROT',
  ont = 'ALL',
  pvalueCutoff = 1,
  qvalueCutoff = 1
)

write_csv(
  Nterm_WP_overlap_GO@result,
  file = 'data_source/Nterm_WP_comparison/Nterm_WP_overlap_GO.csv'
)

## import GO enrichment analysis result
Nterm_WP_overlap_GO <- read_csv(
  'data_source/Nterm_WP_comparison/Nterm_WP_overlap_GO.csv'
)

## ribosome
# GO:0015935, small ribosomal subunit
# GO:0022625, cytosolic large ribosomal subunit
# GO:0000313, organellar ribosome
ribosome_related_protein_list <- Nterm_WP_overlap_GO |> 
  filter(
    Description %in% c(
      'small ribosomal subunit',
      'cytosolic large ribosomal subunit',
      'organellar ribosome'
    )
  ) |> 
  dplyr::select(geneID) |> 
  separate_rows(geneID, sep = '/') |> 
  distinct() |> 
  pull()

## proteasome
# GO:0000502, proteasome complex
# GO:0022624, proteasome accessory complex
# GO:0005838, proteasome regulatory particle
proteasome_related_protein_list <- Nterm_WP_overlap_GO |> 
  filter(
    Description %in% c(
      'proteasome complex',
      'proteasome accessory complex', 
      'proteasome regulatory particle'
    )
  ) |> 
  dplyr::select(geneID) |> 
  separate_rows(geneID, sep = '/') |> 
  distinct() |> 
  pull()

## spliceosome
# GO:0097525, spliceosomal snRNP complex
# GO:0071011, precatalytic spliceosome
# GO:0071006, U2-type catalytic step 1 spliceosome
# GO:0071007, U2-type catalytic step 2 spliceosome
# GO:0071010, prespliceosome
spliceosome_related_protein_list <- Nterm_WP_overlap_GO |> 
  filter(
    Description %in% c(
      'spliceosomal snRNP complex',
      'precatalytic spliceosome',
      'U2-type catalytic step 1 spliceosome',
      'U2-type catalytic step 2 spliceosome',
      'prespliceosome'
    )
  ) |> 
  dplyr::select(geneID) |> 
  separate_rows(geneID, sep = '/') |> 
  distinct() |> 
  pull()

# combine protein from different organelle
Nterm_linear_fit_protein_list <- HEK_Nterm_curve_fitting_combined |> 
  filter(model == 'linear fitting') |> 
  distinct(Index) |> 
  pull()

WP_linear_fit_protein_list <- HEK_WP_curve_fitting_combined |> 
  filter(model == 'linear fitting') |> 
  distinct(UniProt_Accession) |> 
  pull()

ribosome_proteasome_spliceosome_comb <- bind_rows(
  HEK_Nterm_WP_delta_half_life |> 
    filter(UniProt_Accession %in% ribosome_related_protein_list) |> 
    mutate(GO_category = 'ribosome'),
  
  HEK_Nterm_WP_delta_half_life |> 
    filter(UniProt_Accession %in% proteasome_related_protein_list) |> 
    mutate(GO_category = 'proteasome'),
  
  HEK_Nterm_WP_delta_half_life |> 
    filter(UniProt_Accession %in% spliceosome_related_protein_list) |> 
    mutate(GO_category = 'spliceosome')
) |> 
  mutate(
    Nterm_curve_fitting = ifelse(
      Index %in% Nterm_linear_fit_protein_list,
      'linear fitting',
      'nonlinear fitting'
    ),
    WP_curve_fitting = ifelse(
      UniProt_Accession %in% WP_linear_fit_protein_list,
      'linear fitting',
      'nonlinear fitting'
    )
  )

write_csv(
  ribosome_proteasome_spliceosome_comb,
  file = 'data_source/Nterm_WP_comparison/ribosome_proteasome_spliceosome_comb.csv'
)

# Wilcoxon rank-sum test
library(rstatix)

ribosome_proteasome_spliceosome_comb_wilcoxon_test <- ribosome_proteasome_spliceosome_comb |> 
  wilcox_test(delta_half_life ~ GO_category)

# point range plot
library(ggpubr)

point_range_plot_ribosome_proteasome_spliceosome_comb <- ribosome_proteasome_spliceosome_comb |> 
  ggplot() +
  geom_point(
    aes(
      x = GO_category,
      y = delta_half_life
    ),
    position = position_jitter(width = 0.3),
    color = 'black',
    alpha = 0.3,
    size = 0.5
  ) +
  stat_summary(
    aes(
      x = GO_category,
      y = delta_half_life,
      color = GO_category
    ),
    fun.data = 'mean_cl_boot', linewidth = 0.2, size = 0.5, show.legend = FALSE
  ) +
  labs(x = '', y = '') +
  stat_pvalue_manual(
    data = ribosome_proteasome_spliceosome_comb_wilcoxon_test,
    label = 'p.adj.signif',
    y.position = c(70, 100, 130),
    tip.length = 0,
    label.size = 6
  ) +
  scale_color_manual(
    values = c(
      'proteasome' = color_1,
      'ribosome' = color_2,
      'spliceosome' = color_3
    )
  ) +
  theme(
    panel.grid.major = element_line(color = "gray", linewidth = 0.2),
    panel.grid.minor = element_line(color = "gray", linewidth = 0.1),
    axis.text.x = element_text(color = 'black', size = 8, family = 'arial', angle = 30, hjust = 1),
    axis.text.y = element_text(color = 'black', size = 8, family = 'arial')
  )

ggsave(
  filename = 'figures/figure8/point_range_plot_ribosome_proteasome_spliceosome_comb.eps',
  device = cairo_ps,
  plot = point_range_plot_ribosome_proteasome_spliceosome_comb,
  height = 2, width = 1.8, units = 'in',
  fallback_resolution = 1200
)

### figure 8D, ribosome example
# P63244, RACK1, 188, 226, 229


# P62899, RPL31


### figure 8E, top20 and bottom20 structure analysis
## import result from StructureMap
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

### figure 8F, Nucleolin example
# P19338, NCL, Nucleolin
P19338_database_info <- tribble(
  ~ name, ~ start, ~ end,
  'protein', 1, 710,
  'RRM_1', 309, 377,
  'RRM_1', 398, 459,
  'RRM_1', 488, 554,
  'RRM_1', 574, 640
)

P19338_result <- tibble(
  cleavage_site = c(
    308, 311, 325, 332, 363, 
    406, 411, 414, 430, 461,
    491, 494, 526, 530, 576,
    580
  )
)

# example plot
P19338_example <- ggplot() +
  geom_rect(
    data = P19338_database_info,
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
      x = P19338_result$cleavage_site, 
      xend = P19338_result$cleavage_site, 
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
  filename = 'figures/figure8/P19338_example.eps',
  height = 0.2, width = 2, units = 'in'
)

# heatmap for different protein N-termini
library(ComplexHeatmap)
library(circlize)

P19338_N_termini <- str_c(
  'P19338', 
  c(
    308, 311, 325, 332, 363, 
    406, 411, 414, 430, 461,
    491, 494, 526, 530, 576,
    580
  ),
  sep = "_"
)

df <- HEK_Nterm_WP_delta_half_life |> 
  filter(
    Index %in% P19338_N_termini
  ) |> 
  select(Index, half_life.Nt)

mat <- data.matrix(df |> select(half_life.Nt))
rownames(mat) <- df$Index

mat_color = colorRamp2(
  breaks = c(5, 30, 55),
  colors = c('blue', 'yellow', 'red')
)

Heatmap(
  matrix = mat,
  col = mat_color,
  cluster_rows = FALSE
)
