# import packages
packages_names <- c('tidyverse', 'showtext', 'rstatix', 'ggpubr', 'ComplexHeatmap', 'circlize')
lapply(packages_names, require, character.only = TRUE)

### figure 4A, GO and KEGG analysis
## import GO and KEGG results
# Nterm, Fast turnover, half-life < 7 h
Nterm_fast_turnover_GO <- read_csv(
  'data_source/GO_KEGG_analysis/Nterm_fast_turnover_GO.csv'
)
# Nterm, Stable, half-life = 200 h
Nterm_stable_GO <- read_csv(
  'data_source/GO_KEGG_analysis/Nterm_stable_GO.csv'
)

# combine selected terms
Nterm_GO_KEGG <- bind_rows(
  Nterm_fast_turnover_GO |> 
    select(Description, Count, pvalue) |> 
    filter(
      Description %in% c(
        'Spindle organization',
        'Sister chromatid segregation',
        'Nuclear chromosome segregation',
        'Unfolded protein binding',
        'Non-membrane-bounded organelle assembly',
        'Lyase activity',
        'ATP-dependent protein folding chaperone',
        'ATP binding',
        'Aminoacyl-tRNA ligase activity',
        'Intramolecular oxidoreductase activity',
        'DNA helicase activity',
        'ATP hydrolysis activity'
      )
    ) |> 
    mutate(Category = 'Fast turnover'),
  Nterm_stable_GO |> 
    select(Description, Count, pvalue) |> 
    filter(
      Description %in% c(
        'Cell junction',
        'Anchoring junction',
        'Focal adhesion',
        'Cell-substrate junction'
      )
    ) |> 
    mutate(Category = 'Stable')
)

# dot plot
font_add(family = 'arial', regular = 'arial.ttf')
showtext_auto()

dot_plot_Nterm_GO_KEGG <- Nterm_GO_KEGG |> 
  ggplot(aes(x = factor(Category, levels = c('Fast turnover', 'Stable')), 
             y = factor(Description, levels = Nterm_GO_KEGG |> pull(Description)))) +
  geom_point(aes(fill = Category, size = Count, alpha = pvalue), shape = 21) +
  scale_y_discrete(labels = function(x) str_wrap(x, width = 40), position = "left") +
  scale_size(range = c(3, 7), breaks = c(50, 150, 250)) +
  scale_alpha_binned(range = c(1, 0.2), breaks = c(0.001, 0.002, 0.003, 0.004)) +
  scale_fill_manual(
    values = c(
      'Fast turnover' = color_1,
      'Stable' = color_2
    )
  ) +
  theme_bw() +
  theme(
    panel.grid.major = element_line(color = "gray", linewidth = 0.2),
    panel.grid.minor = element_line(color = "gray", linewidth = 0.1),
    axis.title = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_text(size = 9, color = "black", lineheight = 0.05),
    legend.text = element_text(size = 9),
    legend.title = element_text(size = 9, color = "black"),
    legend.frame = element_rect(color = "black", linewidth = 0.2),
    legend.ticks = element_line(color = "black", linewidth = 0.2),
    legend.key.size = unit(0.1, "in")
  )

ggsave(
  'figures/figure4/dot_plot_Nterm_GO_KEGG.eps',
  device = cairo_ps,
  plot = dot_plot_Nterm_GO_KEGG,
  height = 3.5, width = 4.5, units = 'in',
  fallback_resolution = 1200
)

### figure 4B, sub-cellular location analysis
# selected sub-cellular location
HEK_Nterm_Kd_half_life_subcellular_adj <- HEK_Nterm_Kd_half_life_subcellular |> 
  filter(Main.location %in% c(
    'Midbody', 'Midbody ring', 'Centrosome',
    'Nuclear speckles', 'Nucleoli', 'Nucleoli rim', 'Nucleoli fibrillar center', 'Nucleoplasm',
    'Golgi apparatus', 'Endoplasmic reticulum', 'Mitochondrion', 'Cell Junctions', 'Cytoplasmic bodies',
    'Peroxisomes',  'Endosomes', 'Lysosomes'
  ))

# Wilcoxon rank-sum test
subcellular_half_life_wilcox_test <- HEK_Nterm_Kd_half_life_subcellular_adj |> 
  wilcox_test(half_life ~ Main.location) |> 
  add_significance('p') |> 
  filter(p.signif != 'ns')

# point boxplot
font_add(family = 'arial', regular = 'arial.ttf')
showtext_auto()

point_boxplot_Nterm_subcellular <- HEK_Nterm_Kd_half_life_subcellular_adj |> 
  ggplot() +
  geom_point(
    aes(
      x = factor(Main.location, levels = c(
        'Midbody', 'Midbody ring', 'Centrosome',
        'Nuclear speckles', 'Nucleoli', 'Nucleoli rim', 'Nucleoli fibrillar center', 'Nucleoplasm',
        'Golgi apparatus', 'Endoplasmic reticulum', 'Mitochondrion', 'Cell Junctions', 'Cytoplasmic bodies',
        'Peroxisomes',  'Endosomes', 'Lysosomes'
      )), 
      y = half_life
    ), 
    position = position_jitter(width = 0.3),
    color = 'black',
    alpha = 0.3
  ) +
  stat_summary(
    aes(
      x = factor(Main.location, levels = c(
        'Midbody', 'Midbody ring', 'Centrosome',
        'Nuclear speckles', 'Nucleoli', 'Nucleoli rim', 'Nucleoli fibrillar center', 'Nucleoplasm',
        'Golgi apparatus', 'Endoplasmic reticulum', 'Mitochondrion', 'Cell Junctions', 'Cytoplasmic bodies',
        'Peroxisomes',  'Endosomes', 'Lysosomes'
      )), 
      y = half_life
    ),
    fun.data = 'mean_cl_boot', color = color_3, linewidth = 0.2, size = 0.5
  ) +
  labs(x = '', y = '') +
  stat_pvalue_manual(
    data = subcellular_half_life_wilcox_test,
    label = 'p.signif',
    y.position = c(170, 180, 160),
    tip.length = 0,
    coord.flip = TRUE,
    label.size = 6
  ) +
  coord_flip() +
  theme(
    panel.grid.major = element_line(color = 'gray', linewidth = 0.2),
    panel.grid.minor = element_line(color = 'gray', linewidth = 0.1),
    axis.text.x = element_text(color = 'black', size = 10),
    axis.text.y = element_text(color = 'black', size = 9)
  )

ggsave(
  filename = 'figures/figure4/point_boxplot_Nterm_subcellular.eps',
  device = cairo_ps,
  plot = point_boxplot_Nterm_subcellular,
  height = 3.5, width = 3, units = 'in',
  fallback_resolution = 1200
)

### figure 4C, cycle cell heatmap
# cell cycle diagram
cell_cycle_df = data.frame(
  phase = factor(c("G1", "S", "G2", "M"), levels = c("G1", "S", "G2", "M")),
  hour = c(11, 8, 4, 1)
)

cell_cycle_color = c(color_1, color_2, color_3, color_4)

circos.par(start.degree = 90)

circos.initialize(cell_cycle_df$phase, xlim = cbind(rep(0, 4), cell_cycle_df$hour))

circos.track(ylim = c(0, 1), panel.fun = function(x, y) {
  circos.arrow(CELL_META$xlim[1], CELL_META$xlim[2], 
               arrow.head.width = CELL_META$yrange*0.8, arrow.head.length = cm_x(0.5),
               col = cell_cycle_color[CELL_META$sector.numeric.index])
  circos.text(CELL_META$xcenter, CELL_META$ycenter, CELL_META$sector.index,
              cex = 1, col = 'white', facing = "downward")
}, bg.border = NA, track.height = 0.6)

circos.clear()

# heatmap
Nterm_cell_cycle_example <- HEK_Nterm_deg_ratio |> 
  filter(Index %in% c(
    'Q9NXR7_207',
    'Q8IYL3_28',
    'P20248_174',
    'P20248_23',
    'P20248_20',
    'P20248_21',
    'P14635_153',
    'P14635_155',
    'Q02224_230',
    'Q02224_627',
    'Q02224_612',
    'Q9HAW4_873',
    'O00273_105',
    'O00273_108',
    'O75496_90',
    'P62805_68',
    'P62805_61',
    'P62805_69',
    'P17096_1',
    'P17096_8',
    'Q96EA4_535',
    'Q9Y6A5_7',
    'Q9Y6A5_479',
    'Q9UNY4_893',
    'Q16763_138',
    'Q16763_123',
    'P17028_128'
  )) |> 
  select(Index, timepoint, deg_ratio_avg) |> 
  pivot_wider(names_from = timepoint, values_from = deg_ratio_avg)

Nterm_cell_cycle_example_matrix <- data.matrix(Nterm_cell_cycle_example)
rownames(Nterm_cell_cycle_example_matrix) <- Nterm_cell_cycle_example$Index

mat_col <- colorRamp2(
  breaks = c(1.0, 0.75, 0.5),
  colors = c('red', 'yellow', 'blue')
)

Heatmap(
  matrix = Nterm_cell_cycle_example_matrix[,-1],
  col = mat_col,
  cluster_rows = FALSE, 
  cluster_columns = FALSE
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
