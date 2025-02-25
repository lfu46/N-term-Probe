# import packages
library(tidyverse)

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
library(rstatix)
subcellular_half_life_wilcox_test <- HEK_Nterm_Kd_half_life_subcellular_adj |> 
  wilcox_test(half_life ~ Main.location) |> 
  add_significance('p') |> 
  filter(p.signif != 'ns')

# point range plot
library(showtext)
library(ggpubr)

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
library(circlize)

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
library(ComplexHeatmap)

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

### figure 4D, nucleolus complexes
# Wilcoxon rank-sum test
library(rstatix)

Nterm_nucleolus_localization_wilcoxon_test <- Nterm_nucleolus_localization_half_life |> 
  filter(Localization %in% c(
    'DFC', 'PDFC', 'GC', 'GC (aggregation)', 'NR'
  )) |> 
  wilcox_test(half_life ~ Localization) |> 
  add_significance('p') |> 
  filter(p < 0.05)

# point range plot
library(showtext)
library(ggpubr)

font_add(family = 'arial', regular = 'arial.ttf')
showtext_auto()

violin_boxplot_nucleolus_complex <- Nterm_nucleolus_localization_half_life |> 
  filter(Localization %in% c(
    'DFC', 'PDFC', 'GC', 'GC (aggregation)', 'NR'
  )) |> 
  ggplot() +
  geom_point(
    aes(
      x = Localization,
      y = half_life
    ),
    position = position_jitter(width = 0.3),
    color = 'black',
    alpha = 0.3
  ) +
  stat_summary(
    aes(
      x = Localization, 
      y = half_life
    ),
    fun.data = 'mean_cl_boot', color = color_4, linewidth = 0.2, size = 0.5
  ) +
  stat_pvalue_manual(
    data = Nterm_nucleolus_localization_wilcoxon_test, label = 'p.signif',
    tip.length = 0, y.position = c(180, 170), label.size = 6
  ) +
  labs(x = '', y = '') +
  theme(
    panel.grid.major = element_line(color = 'gray', linewidth = 0.2),
    panel.grid.minor = element_line(color = 'gray', linewidth = 0.1),
    axis.text.x = element_blank(),
    axis.text.y = element_text(color = 'black', size = 9)
  )

ggsave(
  filename = 'figures/figure4/violin_boxplot_nucleolus_complex.eps',
  device = cairo_ps,
  plot = violin_boxplot_nucleolus_complex,
  height = 2, width = 2, units = 'in',
  fallback_resolution = 1200
)

### figure 4E, protein interaction network of nucleoli, stress granule, processing body and cajal body
# combine nucleolus and cytoplasmic body proteins
nucleolus_cytoplasmic_body_protein <- bind_rows(
  Nterm_nucleolus_localization_half_life |> 
    select(UniProt_Accession, category, half_life) |> 
    distinct(),
  
  Nterm_cytoplasmic_body_half_life |> 
    select(UniProt_Accession, category, half_life) |> 
    distinct()
)

# only keep the proteins with unique localization
nucleolus_cytoplasmic_body_protein_unique_localization <- nucleolus_cytoplasmic_body_protein |> 
  group_by(UniProt_Accession) |> 
  filter(n() == 1) |> 
  ungroup()

write_csv(
  nucleolus_cytoplasmic_body_protein_unique_localization,
  file = 'figures/figure4/nucleolus_cytoplasmic_body_protein_unique_localization.csv'
)

# Wilcoxon rank-sum test
library(rstatix)

nucleolus_cytoplasmic_body_protein_unique_localization_wilcoxon_test <- nucleolus_cytoplasmic_body_protein_unique_localization |> 
  wilcox_test(half_life ~ category) |> 
  add_significance('p') |> 
  filter(p < 0.05)

# point range plot
library(showtext)
library(ggpubr)

font_add(family = 'arial', regular = 'arial.ttf')
showtext_auto()

point_range_plot_Nterm_nucleolus_cytoplasmic_body <- nucleolus_cytoplasmic_body_protein_unique_localization |> 
  ggplot() +
  geom_point(
    aes(
      x = category,
      y = half_life
    ),
    position = position_jitter(width = 0.3),
    color = 'black',
    alpha = 0.3
  ) +
  stat_summary(
    aes(
      x = category, 
      y = half_life
    ),
    fun.data = 'mean_cl_boot', color = color_1, linewidth = 0.2, size = 0.5
  ) +
  stat_pvalue_manual(
    data = nucleolus_cytoplasmic_body_protein_unique_localization_wilcoxon_test, label = 'p.signif',
    tip.length = 0, y.position = c(180, 170), label.size = 6, coord.flip = TRUE
  ) +
  labs(x = '', y = '') +
  coord_flip() +
  theme(
    panel.grid.major = element_line(color = 'gray', linewidth = 0.2),
    panel.grid.minor = element_line(color = 'gray', linewidth = 0.1),
    axis.text.x = element_text(color = 'black', size = 9),
    axis.text.y = element_blank()
  )

ggsave(
  filename = 'figures/figure4/point_range_plot_Nterm_nucleolus_cytoplasmic_body.eps',
  device = cairo_ps,
  height = 2, width = 2, units = 'in',
  fallback_resolution = 1200
)

# extract unique protein list
library(tidyverse)

nucleolus_cytoplasmic_body_protein_unique_localization <- read_csv(
  'figures/figure4/nucleolus_cytoplasmic_body_protein_unique_localization.csv'
)

unique_protein_list <- nucleolus_cytoplasmic_body_protein_unique_localization |> 
  distinct(UniProt_Accession) |> 
  pull()

# convert to STRING compatible format
protein_query <- paste(unique_protein_list, collapse = ",")

# ensure the connection with Cytoscape
library(RCy3)

cytoscapePing()

# query STRING database via Cytoscape
# Tip 1: run 'help string' in Cytoscape command line or commandsHelp('string') in R to check all the available commands
# Tip 2: run 'help string protein query' in Cytoscape command line or run commandsHelp('string protein query') in R to check the available parameters
commandsRun(paste(
  'string protein query query="', 
  protein_query, 
  '" species="Homo sapiens" limit=0 newNetName="nucleolus_cytoplasmic_body_network"', 
  sep = ""
))

# annotate Node Table with the protein category
loadTableData(
  nucleolus_cytoplasmic_body_protein_unique_localization, 
  data.key.column = "UniProt_Accession",
  table = "node",
  table.key.column = "query term"
)

# STRING functional enrichment
commandsRun('string retrieve enrichment background="genome" selectedNodesOnly=false')

# delte Cytoscape network
deleteAllNetworks()

# heatmap of example Nterm proteoforms half-life
library(ComplexHeatmap)
library(circlize)

protein_list <- c(
  # nucleolus
  'Q14137',
  'O15160',
  'P0DPB6',
  'Q3B726',
  'Q9P1U0',
  # stress granule
  'Q9UN86',
  'Q14671',
  'Q13347',
  'P78344',
  'Q9Y6M1',
  # P-body
  'O15116',
  'Q8NDV7',
  'O00165',
  'Q9UKM9',
  'Q92540',
  # cajal body
  'Q8WWY3',
  'P08621',
  'P09012',
  'Q01081'
)

Nterm_example_half_life <- HEK_Nterm_Kd_half_life_LaminB_Tcomplex |> 
  filter(UniProt_Accession %in% protein_list) |> 
  select(Index, half_life)

Nterm_example_half_life_matrix <- data.matrix(Nterm_example_half_life)
rownames(Nterm_example_half_life_matrix) <- Nterm_example_half_life$Index

mat_col <- colorRamp2(
  breaks = c(0, 25, 50),
  colors = c('blue', 'yellow', 'red')
)

Heatmap(
  matrix = Nterm_example_half_life_matrix[,-1],
  col = mat_col,
  cluster_rows = FALSE
)
