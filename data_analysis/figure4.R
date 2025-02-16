# import packages
packages_names <- c('tidyverse', 'showtext', 'rstatix', 'ggpubr', 'ComplexHeatmap', 'circlize')
lapply(packages_names, require, character.only = TRUE)

### figure 4A
## import GO and KEGG results
# Nterm short half life, top 25%
Nterm_short_half_life_GO <- read_csv(
  'data_source/GO_KEGG_analysis/Nterm_short_half_life_GO.csv'
)
# Nterm long half life, bottom 25%
Nterm_long_half_life_GO <- read_csv(
  'data_source/GO_KEGG_analysis/Nterm_long_half_life_GO.csv'
)

# combine selected terms
Nterm_GO_KEGG <- bind_rows(
  Nterm_short_half_life_GO |> 
    select(Description, Count, p.adjust) |> 
    filter(
      Description %in% c(
        'Extracellular space',
        'DNA biosynthetic process',
        'DNA helicase activity',
        'ATP-dependent protein folding chaperone',
        'ATP-dependent activity',
        'ATP hydrolysis activity',
        'Small molecule metabolic process',
        'Unfolded protein binding'
      )
    ) |> 
    mutate(category = 'Top 25%'),
  Nterm_long_half_life_GO |> 
    select(Description, Count, p.adjust) |> 
    filter(
      Description %in% c(
        'Cytoplasmic stress granule',
        'Ribonucleoprotein granule',
        'Cell surface',
        'mRNA binding',
        'miRNA binding',
        'Regulation of protein import into nucleus',
        'Regulation of cytoplasmic translation',
        'Positive regulation of translation',
        'Regulation of protein metabolic process'
      )
    ) |> 
    mutate(category = 'Bottom 25%')
)

# dot plot
font_add(family = 'arial', regular = 'arial.ttf')
showtext_auto()

dot_plot_Nterm_GO_KEGG <- Nterm_GO_KEGG |> 
  ggplot(aes(x = factor(category, levels = c('Top 25%', 'Bottom 25%')), 
             y = factor(Description, levels = Nterm_GO_KEGG |> pull(Description)))) +
  geom_point(aes(fill = category, size = Count, alpha = p.adjust), shape = 21) +
  scale_y_discrete(labels = function(x) str_wrap(x, width = 40), position = "left") +
  scale_size(range = c(3, 7), breaks = c(50, 125, 200)) +
  scale_alpha_binned(range = c(1, 0.2), breaks = c(0.01, 0.02, 0.03, 0.04)) +
  scale_fill_manual(
    values = c(
      'Top 25%' = color_1,
      'Bottom 25%' = color_2
    )
  ) +
  theme_bw() +
  theme(
    panel.grid.major = element_line(color = "gray", linewidth = 0.2),
    panel.grid.minor = element_line(color = "gray", linewidth = 0.1),
    axis.title = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_text(size = 90, color = "black", lineheight = 0.05),
    legend.text = element_text(size = 90),
    legend.title = element_text(size = 90, color = "black"),
    legend.frame = element_rect(color = "black", linewidth = 0.2),
    legend.ticks = element_line(color = "black", linewidth = 0.2),
    legend.key.size = unit(0.1, "in")
  )

ggsave(
  'figures/figure4/dot_plot_Nterm_GO_KEGG.tiff',
  plot = dot_plot_Nterm_GO_KEGG,
  height = 3.5, width = 4, units = 'in',
  dpi = 1200
)

### figure 4B
# selected sub-cellular location
HEK_Nterm_Kd_half_life_subcellular_adj <- HEK_Nterm_Kd_half_life_subcellular |> 
  filter(Main.location %in% c(
    'Focal adhesion sites', 
    'Nucleoplasm', 'Nuclear speckles', 'Nucleoli', 'Nucleoli rim', 'Nucleoli fibrillar center', 
    'Peroxisomes',  'Endosomes', 'Lysosomes',
    'Golgi apparatus', 'Endoplasmic reticulum', 'Mitochondria', 'Cytosol',
    'Cytoplasmic bodies'
  ))

# Wilcoxon rank-sum test
subcellular_half_life_wilcox_test <- HEK_Nterm_Kd_half_life_subcellular_adj |> 
  wilcox_test(half_life ~ Main.location) |> 
  add_significance('p') |> 
  filter(p.signif != 'ns')

# point plot
font_add(family = 'arial', regular = 'arial.ttf')
showtext_auto()

point_boxplot_Nterm_subcellular <- HEK_Nterm_Kd_half_life_subcellular_adj |> 
  ggplot() +
  geom_point(
    aes(x = fct_reorder(Main.location, half_life), y = half_life), 
    position = position_jitter(width = 0.3),
    color = 'black',
    alpha = 0.3
  ) +
  geom_boxplot(
    aes(x = fct_reorder(Main.location, half_life), y = half_life),
    color = color_1, fill = 'transparent', outliers = FALSE
  ) +
  labs(x = '', y = '') +
  # stat_pvalue_manual(
  #   data = subcellular_half_life_wilcox_test,
  #   label = 'p.signif',
  #   y.position = c(42, 46, 50, 38), 
  #   tip.length = 0, 
  #   coord.flip = TRUE,
  #   label.size = 60
  # ) +
  coord_flip(ylim = c(0, 40)) +
  theme_bw() +
  theme(
    panel.grid.major = element_line(color = 'gray', linewidth = 0.2),
    panel.grid.minor = element_line(color = 'gray', linewidth = 0.1),
    axis.text.x = element_text(color = 'black', size = 100),
    axis.text.y = element_text(color = 'black', size = 90)
  )

ggsave(
  filename = 'figures/figure4/point_boxplot_Nterm_subcellular.tiff',
  plot = point_boxplot_Nterm_subcellular,
  height = 3.5, width = 3, units = 'in',
  dpi = 1200
)

### figure 4C
## import GSEA enrichment results
# descending Kd
Nterm_protein_functional_enrichment_des_Kd <- read_csv(
  'data_source/Enzyme_Motif_Domain_Complex_analysis/Nterm_protein_functional_property_GSEA_des_Kd.csv'
) |> 
  filter(pvalue < 0.05) |> 
  mutate(
    NES = 0-NES
  ) |> 
  select(Description, setSize, NES, pvalue)

# descending half life
Nterm_protein_functional_enrichment_des_half_life <- read_csv(
  'data_source/Enzyme_Motif_Domain_Complex_analysis/Nterm_protein_functional_property_GSEA_des_half_life.csv'
) |> 
  filter(pvalue < 0.05) |> 
  select(Description, setSize, NES, pvalue)

# combine descending Kd and descending half life result
Nterm_protein_functional_enrichment_combined <- bind_rows(
  Nterm_protein_functional_enrichment_des_Kd,
  Nterm_protein_functional_enrichment_des_half_life
)

# dot plot
font_add(family = 'arial', regular = 'arial.ttf')
showtext_auto()

dot_plot_Nterm_protein_functional_enrichment <- Nterm_protein_functional_enrichment_combined |> 
  mutate(
    log_pvalue = -log10(pvalue),
    category = case_when(
      str_detect(Description, 'PF') ~ 'PFAM',
      str_detect(Description, 'PS') ~ 'PROSITE',
      str_detect(Description, 'ENZYME') ~ 'ENZYME'
    )
  ) |> 
  ggplot() +
  geom_point(
    aes(x = NES, y = log_pvalue, size = setSize, fill = category), shape = 21
  ) +
  scale_fill_manual(
    name = '',
    values = c(
      'PFAM' = color_1, 
      'PROSITE' = color_2,
      'ENZYME' = color_3
    ),
    guide = 'none'
  ) +
  scale_size(range = c(4, 7)) +
  labs(x = '', y = '') +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    axis.text = element_text(color = 'black', size = 9)
  )

ggsave(
  filename = 'figures/figure4/dot_plot_Nterm_protein_functional_enrichment.eps',
  plot = dot_plot_Nterm_protein_functional_enrichment,
  height = 2, width = 3, units = 'in'
)


### figure 4D, mitochondrial complexes
# Wilcoxon rank-sum test
Nterm_mito_complex_wilcoxon_test <- Nterm_mito_half_life |> 
  filter(
    str_detect(complex_name, 'MICOS|ribosom|TIM23|Respiratory')
  ) |> 
  mutate(
    complex_name = case_when(
      str_detect(complex_name, 'MICOS') ~ 'MICOS complex related',
      str_detect(complex_name, 'ribosom') ~ 'Mitochondrial ribosome related',
      str_detect(complex_name, 'TIM23') ~ 'TIM23 complex related',
      str_detect(complex_name, 'Respiratory chain complex') ~ 'Respiratory chain complex I'
    )
  ) |> 
  wilcox_test(half_life ~ complex_name) |> 
  add_significance('p') |> 
  filter(p < 0.05)

# violin boxplot
font_add(family = 'arial', regular = 'arial.ttf')
showtext_auto()

violin_boxplot_mito_complex <- Nterm_mito_half_life |> 
  filter(
    str_detect(complex_name, 'MICOS|ribosom|TIM23|Respiratory')
  ) |> 
  mutate(
    complex_name = case_when(
      str_detect(complex_name, 'MICOS') ~ 'MICOS complex related',
      str_detect(complex_name, 'ribosom') ~ 'Mitochondrial ribosome related',
      str_detect(complex_name, 'TIM23') ~ 'TIM23 complex related',
      str_detect(complex_name, 'Respiratory chain complex') ~ 'Respiratory chain complex I'
    )
  ) |> 
  ggplot() +
  geom_violin(
    aes(x = complex_name, y = half_life), fill = color_1, color = 'transparent'
  ) +
  geom_boxplot(
    aes(x = complex_name, y = half_life), width = 0.15
  ) +
  stat_pvalue_manual(
    data = Nterm_mito_complex_wilcoxon_test, label = 'p.signif',
    tip.length = 0, y.position = c(42, 46), label.size = 6
  ) +
  coord_cartesian(ylim = c(0, 50)) +
  labs(x = '', y = '') +
  theme_bw() +
  theme(
    panel.grid.major = element_line(color = 'gray', linewidth = 0.2),
    panel.grid.minor = element_line(color = 'gray', linewidth = 0.1),
    axis.text.x = element_blank(),
    axis.text.y = element_text(color = 'black', size = 9)
  )

ggsave(
  filename = 'figures/figure4/violin_boxplot_mito_complex.eps',
  plot = violin_boxplot_mito_complex,
  height = 2, width = 2, units = 'in'
)

### figure 4E, nucleolus complexes
# Wilcoxon rank-sum test
Nterm_nucleolus_localization_wilcoxon_test <- Nterm_nucleolus_localization_half_life |> 
  filter(Localization %in% c(
    'DFC', 'PDFC', 'GC', 'GC (aggregation)', 'NR'
  )) |> 
  wilcox_test(half_life ~ Localization) |> 
  add_significance('p') |> 
  filter(p < 0.05)

# violin boxplot
font_add(family = 'arial', regular = 'arial.ttf')
showtext_auto()

violin_boxplot_nucleolus_complex <- Nterm_nucleolus_localization_half_life |> 
  filter(Localization %in% c(
    'DFC', 'PDFC', 'GC', 'GC (aggregation)', 'NR'
  )) |> 
  ggplot() +
  geom_violin(
    aes(x = Localization, y = half_life), fill = color_2, color = 'transparent'
  ) +
  geom_boxplot(
    aes(x = Localization, y = half_life), width = 0.15
  ) +
  stat_pvalue_manual(
    data = Nterm_nucleolus_localization_wilcoxon_test, label = 'p.signif',
    tip.length = 0, y.position = c(38, 42, 46, 50), label.size = 6
  ) +
  coord_cartesian(ylim = c(0, 50)) +
  labs(x = '', y = '') +
  theme_bw() +
  theme(
    panel.grid.major = element_line(color = 'gray', linewidth = 0.2),
    panel.grid.minor = element_line(color = 'gray', linewidth = 0.1),
    axis.text.x = element_blank(),
    axis.text.y = element_text(color = 'black', size = 9)
  )

ggsave(
  filename = 'figures/figure4/violin_boxplot_nucleolus_complex.eps',
  plot = violin_boxplot_nucleolus_complex,
  height = 2, width = 2, units = 'in'
)

### figure 4F, Nterm nucleolus proteoforms heatmap
# PDFC, 35 proteoforms
Nterm_PDFC <- Nterm_nucleolus_localization_half_life |> 
  filter(Localization == 'PDFC') |> 
  select(Index, half_life) |> 
  arrange(desc(half_life))

Nterm_PDFC_matrix <- data.matrix(Nterm_PDFC)
rownames(Nterm_PDFC_matrix) <- Nterm_PDFC$Index

mat_col <- colorRamp2(
  breaks = c(0, 25, 50), colors = c('blue', 'white', 'red')
)

Heatmap(
  matrix = Nterm_PDFC_matrix[,2],
  col = mat_col,
  cluster_rows = FALSE
)

# DFC, 10 proteoforms
Nterm_DFC <- Nterm_nucleolus_localization_half_life |> 
  filter(Localization == 'DFC') |> 
  select(Index, half_life) |> 
  arrange(desc(half_life))

Nterm_DFC_matrix <- data.matrix(Nterm_DFC)
rownames(Nterm_DFC_matrix) <- Nterm_DFC$Index

mat_col <- colorRamp2(
  breaks = c(0, 25, 50), colors = c('blue', 'white', 'red')
)

Heatmap(
  matrix = Nterm_DFC_matrix[,2],
  col = mat_col,
  cluster_rows = FALSE
)

# NR, 23 proteoforms
Nterm_NR <- Nterm_nucleolus_localization_half_life |> 
  filter(Localization == 'NR') |> 
  select(Index, half_life) |> 
  arrange(desc(half_life))

Nterm_NR_matrix <- data.matrix(Nterm_NR)
rownames(Nterm_NR_matrix) <- Nterm_NR$Index

mat_col <- colorRamp2(
  breaks = c(0, 25, 50), colors = c('blue', 'white', 'red')
)

Heatmap(
  matrix = Nterm_NR_matrix[,2],
  col = mat_col,
  cluster_rows = FALSE
)
