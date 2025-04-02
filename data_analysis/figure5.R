# import packages
library(tidyverse)

### figure 5A, GO and KEGG analysis
## import GO and KEGG results
# Nterm, Fast turnover, Top 20%
Nterm_fast_turnover_GO <- read_csv(
  'data_source/GO_KEGG_analysis/Nterm_fast_turnover_GO.csv'
)
Nterm_fast_turnover_KEGG <- read_csv(
  'data_source/GO_KEGG_analysis/Nterm_fast_turnover_KEGG.csv'
)
# Nterm, Stable, Bottom 20%
Nterm_stable_GO <- read_csv(
  'data_source/GO_KEGG_analysis/Nterm_stable_GO.csv'
)
Nterm_stable_KEGG <- read_csv(
  'data_source/GO_KEGG_analysis/Nterm_stable_KEGG.csv'
)

# combine selected terms
Nterm_GO_KEGG <- bind_rows(
  Nterm_fast_turnover_GO |> 
    select(Description, Count, pvalue) |> 
    filter(
      Description %in% c(
        'Regulation of telomere maintenance',
        'mRNA binding',
        'Telomere organization',
        'mRNA processing',
        'ATP-dependent protein folding chaperone',
        'Spliceosomal complex',
        'Chromosome organization',
        'Cytoplasmic stress granule',
        'Protein localization to nucleus',
        'RNA splicing'
      )
    ) |> 
    mutate(Category = 'Fast turnover'),
  
  Nterm_fast_turnover_KEGG |> 
    select(Description, Count, pvalue) |> 
    filter(
      Description %in% c(
        'Cell cycle',
        'Spliceosome'
      )
    ) |> 
    mutate(Category = 'Fast turnover'),
  
  Nterm_stable_GO |> 
    select(Description, Count, pvalue) |> 
    filter(
      Description %in% c(
        'ATP-dependent protein folding chaperone',
        'Cytoplasmic translation',
        'Translation',
        'Nuclear periphery',
        'Structural molecule activity',
        'Hydrolase activity',
        'mRNA stabilization',
        'Carbohydrate derivative binding',
        'Endonuclease activity'
      )
    ) |> 
    mutate(Category = 'Stable'),
  
  Nterm_stable_KEGG |> 
    select(Description, Count, pvalue) |> 
    filter(
      Description %in% c(
        'Cell cycle',
        'Spliceosome'
      )
    ) |> 
    mutate(Category = 'Stable')
)

# dot plot
dot_plot_Nterm_GO_KEGG <- Nterm_GO_KEGG |> 
  ggplot(aes(x = factor(Category, levels = c('Fast turnover', 'Stable')), 
             y = factor(Description, levels = Nterm_GO_KEGG |> distinct(Description) |> pull()))) +
  geom_point(aes(fill = Category, size = Count, alpha = pvalue), shape = 21) +
  scale_y_discrete(labels = function(x) str_wrap(x, width = 40), position = "left") +
  scale_size(range = c(3, 7), breaks = c(50, 150, 250)) +
  scale_alpha_binned(range = c(1, 0.2), breaks = c(0.002, 0.005, 0.008)) +
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
    axis.text.y = element_text(size = 8, color = "black", lineheight = 0.05, family = 'arial'),
    legend.text = element_text(size = 8),
    legend.title = element_text(size = 8, color = "black"),
    legend.frame = element_rect(color = "black", linewidth = 0.2),
    legend.ticks = element_line(color = "black", linewidth = 0.2),
    legend.key.size = unit(0.1, "in")
  )

ggsave(
  'figures/figure5/dot_plot_Nterm_GO_KEGG.eps',
  device = cairo_ps,
  plot = dot_plot_Nterm_GO_KEGG,
  height = 3.5, width = 4, units = 'in',
  fallback_resolution = 1200
)

### figure 5B, sub-cellular location analysis
# selected sub-cellular location
HEK_Nterm_Kd_half_life_subcellular_adj <- HEK_Nterm_Kd_half_life_subcellular |> 
  filter(Main.location %in% c(
    'Centrosome', 'Centriolar satellite',
    'Cytoplasmic bodies',
    'Nuclear speckles', 'Nucleoli', 'Nucleoli rim', 'Nuclear bodies','Nucleoli fibrillar center', 'Nucleoplasm', 
    'Golgi apparatus', 'Cytosol', 'Endoplasmic reticulum', 'Plasma membrane','Mitochondria', 
    'Peroxisomes'
  ))

write_csv(
  HEK_Nterm_Kd_half_life_subcellular_adj,
  file = 'data_source/HPA_subcellular_location/HEK_Nterm_Kd_half_life_subcellular_adj.csv'
)

# Wilcoxon rank-sum test
library(rstatix)

HEK_Nterm_Kd_half_life_subcellular_adj |> 
  count(Main.location)

subcellular_half_life_wilcoxon_test <- HEK_Nterm_Kd_half_life_subcellular_adj |> 
  wilcox_test(half_life ~ Main.location) |> 
  add_significance('p') |> 
  filter(p.signif != 'ns')

# point range plot
library(ggpubr)

boxplot_Nterm_subcellular <- HEK_Nterm_Kd_half_life_subcellular_adj |> 
  ggplot() +
  geom_boxplot(
    aes(
      x = factor(Main.location, levels = c(
        'Centrosome', 'Centriolar satellite',
        'Cytoplasmic bodies',
        'Nuclear speckles', 'Nucleoli', 'Nucleoli rim', 'Nuclear bodies','Nucleoli fibrillar center', 'Nucleoplasm', 
        'Golgi apparatus', 'Cytosol', 'Endoplasmic reticulum', 'Plasma membrane','Mitochondria', 
        'Peroxisomes'
      )), 
      y = half_life
    ), 
    color = color_3, outliers = FALSE
  ) +
  labs(x = '', y = '') +
  stat_pvalue_manual(
    data = subcellular_half_life_wilcoxon_test,
    label = 'p.signif',
    y.position = c(43, 40, 37, 34, 31, 28),
    tip.length = 0,
    coord.flip = TRUE,
    label.size = 6
  ) +
  coord_flip(ylim = c(0, 45)) +
  theme(
    panel.grid.major = element_line(color = 'gray', linewidth = 0.2),
    panel.grid.minor = element_line(color = 'gray', linewidth = 0.1),
    axis.text.x = element_text(color = 'black', size = 8, family = 'arial'),
    axis.text.y = element_text(color = 'black', size = 8, family = 'arial')
  )

ggsave(
  filename = 'figures/figure5/boxplot_Nterm_subcellular.eps',
  device = 'eps',
  plot = boxplot_Nterm_subcellular,
  height = 3.5, width = 3, units = 'in'
)

### figure 5C, cycle cell heatmap
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

# cell-cycle related proteoform examples
Nterm_cell_cycle_example <- HEK_Nterm_linear_fitting |> 
  filter(Index %in% c(
    # G2/M
    'Q02224_230',
    'Q02224_612',
    'P62805_25',
    'P62805_46',
    'P62805_47',
    'P62805_61',
    'P62805_68',
    'Q96EA4_30',
    'Q96EA4_535',
    'P11388_1177',
    'P11388_1392',
    'Q16763_45',
    'Q16763_123',
    'Q16763_140',
    # G1 + S
    'P17096_1',
    'P17096_8',
    'O00273_105',
    'O00273_108',
    'O00273_114',
    # S + G2/M
    'Q9NXR7_204',
    'Q9NXR7_205',
    'Q9NXR7_207',
    # S
    'P17028_28',
    'P17028_128'
  ))

## generate ratio using the fitted model
# linear model
linear_model <- function(lnA, Kd_adj) {
  
  timepoint <- c(3, 6, 9, 12, 24)
  
  exp(lnA - Kd_adj * timepoint)
}

Nterm_cell_cycle_linear_fitting <- Nterm_cell_cycle_example |> 
  pivot_wider(
    names_from = parameters, values_from = values
  ) |> 
  left_join(HEK_Nterm_Kd_half_life_LaminB_Tcomplex, by = c(
    'Index', 'UniProt_Accession', 'Protein.Start', 'Gene', 'Entry.Name'
  )) |> 
  mutate(ratio = map2(lnA, Kd_adj, linear_model)) |> 
  unnest_wider(ratio, names_sep = "_")

write_csv(
  Nterm_cell_cycle_linear_fitting,
  file = 'data_source/cell_cycle/Nterm_cell_cycle_linear_fitting.csv'
)

## import result of the raito of example N-terminal proteoforms
library(tidyverse)

# linear fitting
Nterm_cell_cycle_linear_fitting <- read_csv(
  'data_source/cell_cycle/Nterm_cell_cycle_linear_fitting.csv'
)

# heatmap
library(ComplexHeatmap)
library(circlize)

Nterm_cell_cycle_curve_fitting_adj <- Nterm_cell_cycle_linear_fitting |> 
  mutate(ratio_0 = 1)

Nterm_cell_cycle_example_matrix <- data.matrix(
  Nterm_cell_cycle_curve_fitting_adj |> 
    select(ratio_0, ratio_1:ratio_5)
)
rownames(Nterm_cell_cycle_example_matrix) <- Nterm_cell_cycle_curve_fitting_adj$Index

mat_col <- colorRamp2(
  breaks = c(1.0, 0.65, 0.3),
  colors = c('red', 'yellow', 'blue')
)

Heatmap(
  matrix = Nterm_cell_cycle_example_matrix,
  col = mat_col,
  cluster_rows = FALSE, 
  cluster_columns = FALSE
)

### figure 5D, protein interaction network of nucleoli, stress granule, processing body and cajal body
# combine nucleolus and cytoplasmic body proteins
nucleolus_cytoplasmic_body_protein <- bind_rows(
  Nterm_nucleolus_localization_half_life |> 
    select(UniProt_Accession, category) |> 
    distinct(),
  
  Nterm_cytoplasmic_body_half_life |> 
    select(UniProt_Accession, category) |> 
    distinct()
)

# only keep the proteins with unique localization
nucleolus_cytoplasmic_body_protein_unique_localization <- nucleolus_cytoplasmic_body_protein |> 
  count(UniProt_Accession) |> 
  filter(n == 1)

# get half-life for those proteoforms with unique parent protein localization
nucleolus_cytoplasmic_body_Nterm_half_life <- bind_rows(
  Nterm_nucleolus_localization_half_life |> 
    select(Index, UniProt_Accession, Gene = Gene.Symbol, category, half_life) |> 
    distinct() |> 
    filter(
      UniProt_Accession %in% nucleolus_cytoplasmic_body_protein_unique_localization$UniProt_Accession
    ),
  
  Nterm_cytoplasmic_body_half_life |> 
    select(Index, UniProt_Accession, Gene, category, half_life)
) |> 
  filter(
    UniProt_Accession %in% nucleolus_cytoplasmic_body_protein_unique_localization$UniProt_Accession
  )

write_csv(
  nucleolus_cytoplasmic_body_Nterm_half_life,
  file = 'data_source/nucleolus_cytoplasmic_body/nucleolus_cytoplasmic_body_Nterm_half_life.csv'
)

# Wilcoxon rank-sum test
library(tidyverse)

nucleolus_cytoplasmic_body_Nterm_half_life <- read_csv(
  'data_source/nucleolus_cytoplasmic_body/nucleolus_cytoplasmic_body_Nterm_half_life.csv'
)

library(rstatix)

nucleolus_cytoplasmic_body_protein_unique_localization_wilcoxon_test <- nucleolus_cytoplasmic_body_Nterm_half_life |> 
  count(category)
  # group_by(category) |>
  # get_summary_stats(half_life)
  wilcox_test(half_life ~ category) |> 
  add_significance('p') |>
  filter(p.signif != 'ns')

# point range plot
library(ggpubr)

point_range_plot_Nterm_nucleolus_cytoplasmic_body <- nucleolus_cytoplasmic_body_Nterm_half_life |> 
  ggplot() +
  # geom_point(
  #   aes(
  #     x = category,
  #     y = half_life
  #   ),
  #   position = position_jitter(width = 0.3),
  #   color = 'black',
  #   alpha = 0.3
  # ) +
  stat_summary(
    aes(
      x = category, 
      y = half_life
    ),
    fun.data = 'mean_cl_boot', color = color_1, linewidth = 0.2,
  ) +
  stat_pvalue_manual(
    data = nucleolus_cytoplasmic_body_protein_unique_localization_wilcoxon_test, label = 'p.signif',
    tip.length = 0, y.position = c(24, 23), label.size = 6, coord.flip = TRUE
  ) +
  labs(x = '', y = '') +
  coord_flip(ylim = c(15, 25)) +
  theme(
    panel.grid.major = element_line(color = 'gray', linewidth = 0.2),
    panel.grid.minor = element_line(color = 'gray', linewidth = 0.1),
    axis.text.x = element_text(color = 'black', size = 8),
    axis.text.y = element_blank()
  )

ggsave(
  filename = 'figures/figure5/point_range_plot_Nterm_nucleolus_cytoplasmic_body.eps',
  device = 'eps',
  height = 2, width = 2, units = 'in'
)

# extract unique protein list
library(tidyverse)

nucleolus_cytoplasmic_body_Nterm_half_life <- read_csv(
  'data_source/nucleolus_cytoplasmic_body/nucleolus_cytoplasmic_body_Nterm_half_life.csv'
)

unique_protein_list <- nucleolus_cytoplasmic_body_Nterm_half_life |> 
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
  nucleolus_cytoplasmic_body_Nterm_half_life, 
  data.key.column = "UniProt_Accession",
  table = "node",
  table.key.column = "query term"
)

# Layout, and then Group Attributes Layout, then select category
# Use the filter in sidebar to select specific group of proteins 

# STRING functional enrichment
commandsRun('string retrieve enrichment background="genome" selectedNodesOnly=false')

# delte Cytoscape network
deleteAllNetworks()

# heatmap of example Nterm proteoforms half-life
library(ComplexHeatmap)
library(circlize)

proteoform_list <- c(
  ## nucleolus
  # P06748, NPM1
  'P06748_89',
  'P06748_90',
  'P06748_91',
  'P06748_211',
  'P06748_213',
  ## stress granule
  # P60842, EIF4A1
  'P60842_147',
  'P60842_149',
  'P60842_212',
  'P60842_213',
  'P60842_214',
  ## P-body
  # Q9NPI6, DCP1A
  'Q9NPI6_383',
  'Q9NPI6_420',
  'Q9NPI6_502',
  # Q9NRA8, EIF4ENIF1
  'Q9NRA8_282',
  'Q9NRA8_434',
  'Q9NRA8_923',
  ## cajal body
  # Q9UQ35, SRRM2
  'Q9UQ35_2308',
  'Q9UQ35_2310',
  'Q9UQ35_2311',
  'Q9UQ35_2313',
  'Q9UQ35_2359',
  'Q9UQ35_2360',
  'Q9UQ35_2361'
)

Nterm_example_half_life <- HEK_Nterm_Kd_half_life_LaminB_Tcomplex |> 
  filter(Index %in% proteoform_list) |> 
  select(Index, half_life)

Nterm_example_half_life_matrix <- data.matrix(Nterm_example_half_life |> select(half_life))
rownames(Nterm_example_half_life_matrix) <- Nterm_example_half_life$Index

mat_col <- colorRamp2(
  breaks = c(10, 23, 36),
  colors = c('blue', 'yellow', 'red')
)

Heatmap(
  matrix = Nterm_example_half_life_matrix,
  col = mat_col,
  cluster_rows = FALSE
)

### figure 5E, nucleolus complexes
# Wilcoxon rank-sum test
library(rstatix)

Nterm_nucleolus_localization_wilcoxon_test <- Nterm_nucleolus_localization_half_life |> 
  # count(Localization)
  filter(Localization %in% c(
    'DFC', 'FC', 'PDFC', 'PNC', 'GC', 'GC (aggregation)', 'NR'
  )) |> 
  # group_by(Localization) |> 
  # get_summary_stats(half_life, type = 'mean')
  wilcox_test(half_life ~ Localization) |> 
  add_significance('p') |> 
  filter(p < 0.05)

# point range plot
library(ggpubr)

boxplot_nucleolus_complex <- Nterm_nucleolus_localization_half_life |> 
  filter(Localization %in% c(
    'DFC', 'FC', 'PDFC', 'PNC', 'GC', 'GC (aggregation)', 'NR'
  )) |> 
  ggplot() +
  geom_boxplot(
    aes(
      x = Localization,
      y = half_life
    ),
    color = color_4, outliers = FALSE
  ) +
  stat_pvalue_manual(
    data = Nterm_nucleolus_localization_wilcoxon_test, label = 'p.signif',
    tip.length = 0, y.position = c(33), label.size = 6
  ) +
  labs(x = '', y = '') +
  coord_cartesian(ylim = c(10, 35)) +
  theme(
    panel.grid.major = element_line(color = 'gray', linewidth = 0.2),
    panel.grid.minor = element_line(color = 'gray', linewidth = 0.1),
    axis.text.x = element_text(color = 'black', size = 8, angle = 90, hjust = 1),
    axis.text.y = element_text(color = 'black', size = 8)
  )

ggsave(
  filename = 'figures/figure5/boxplot_nucleolus_complex.eps',
  device = 'eps',
  plot = boxplot_nucleolus_complex,
  height = 3, width = 2, units = 'in'
)
