# import packages
library(tidyverse)

### figure 6A, N-terminus secondary structure, solvent accessibility and IDR analysis
## Tip 1: NetSurfP - 3.0 does not accept protein sequence more than 5,000 residues. 
## You can manually check the sequence length and split longer sequences into shorter segments.
## Tip 2: NetSurfP - 3.0 does not accept U as an amino acid residue.
## You can use the GRAVY Calculator (https://www.gravy-calculator.de/index.php) to locate U residues in your FASTA file and replace them with C.
## Tip 3: Do not use NetSurfP - 3.0.
## Wilcoxon rank-sum test
library(rstatix)

Nterm_degradation_alphafold_half_life <- read_csv(
  'data_source/Nterm_degradation_structuremap/Nterm_degradation_alphafold_half_life.csv'
)

Nterm_degradation_alphafold_half_life_property <- Nterm_degradation_alphafold_half_life |> 
  select(Index, half_life, structure_group, accessibility, IDR) |> 
  pivot_longer(cols = c(structure_group, accessibility, IDR), names_to = 'property', values_to = 'category')

Nterm_degradation_alphafold_half_life_property |> 
  group_by(category) |> 
  get_summary_stats(half_life, type = 'median')

write_csv(
  Nterm_degradation_alphafold_half_life_property, 
  'data_source/Nterm_degradation_structuremap/Nterm_degradation_alphafold_half_life_property.csv'
)

Nterm_degradation_alphafold_half_life_wilcox_test <- Nterm_degradation_alphafold_half_life_property |> 
  group_by(property) |> 
  wilcox_test(half_life ~ category) |> 
  add_significance('p') |> 
  filter(p.adj.signif != 'ns')

# point range plot
library(ggpubr)

point_range_Nterm_degradation_alphafold_property <- Nterm_degradation_alphafold_half_life_property |> 
  ggplot() +
  stat_summary(
    aes(
      x = category, y = half_life, colour = property
    ),
    fun.data = 'mean_cl_boot', linewidth = 0.2, size = 0.5, show.legend = FALSE
  ) +
  scale_color_manual(
    values = c(
      'structure_group' = color_1,
      'accessibility' = color_2,
      'IDR' = color_3
    )
  ) +
  labs(x = '', y = '') +
  stat_pvalue_manual(
    data = Nterm_degradation_alphafold_half_life_wilcox_test, label = 'p.signif', tip.length = 0,
    size = 6, hide.ns = 'p', y.position = c(20.8, 21.2, 21.6, 21.2, 20.8)
  ) +
  facet_grid(cols = vars(factor(property, levels = c(
    'structure_group', 'accessibility', 'IDR'
  ))), scales = 'free_x', space = 'free_x') +
  theme(
    panel.grid.major = element_line(color = "gray", linewidth = 0.2),
    panel.grid.minor = element_line(color = "gray", linewidth = 0.1),
    axis.text.x = element_text(family = 'arial', color = 'black', angle = 30, hjust = 1, size = 8),
    axis.text.y = element_text(family = 'arial', color = 'black', size = 8),
  )

ggsave(
  filename = 'figures/figure6/point_range_Nterm_degradation_alphafold_property.eps',
  plot = point_range_Nterm_degradation_alphafold_property,
  device = 'eps',
  height = 2.5, width = 4, units = 'in'
)

### figure 6B, hydropathy, NCPR and kappa of Nterm 13-mer 
# import the result from localCIDER (data_analysis/Nterm_sequence_features.py)
library(tidyverse)

HEK_Nterm_13mer_sequence_features <- read_csv(
  'data_source/Nterm_13mer_sequence_features/HEK_Nterm_13mer_Kd_half_life_sequence_features.csv'
)

HEK_Nterm_13mer_sequence_features_selected <- HEK_Nterm_13mer_sequence_features |> 
  select(Index, category, hydropathy, NCPR, kappa) |> 
  pivot_longer(cols = hydropathy:kappa, names_to = 'features', values_to = 'values')

# Wilcoxon rank-sum test
library(rstatix)

HEK_Nterm_13mer_sequence_features_wilcoxon_test <- HEK_Nterm_13mer_sequence_features_selected |> 
  # group_by(category, features) |> 
  # get_summary_stats(values)
  group_by(features) |>
  wilcox_test(values ~ category) |> 
  filter(p.adj.signif != "ns")

# boxplot
library(ggpubr)

boxplot_HEK_Nterm_13mer_sequence_features <- HEK_Nterm_13mer_sequence_features_selected |> 
  ggplot() +
  geom_boxplot(
    aes(
      x = category, y = values, color = features
    ), show.legend = FALSE
  ) +
  scale_color_manual(
    values = c(
      'hydropathy' = color_2,
      'NCPR' = color_3,
      'kappa' = color_4
    )
  ) +
  labs(x = '', y = '') +
  stat_pvalue_manual(
    data = HEK_Nterm_13mer_sequence_features_wilcoxon_test, label = 'p.adj.signif', tip.length = 0,
    size = 6, hide.ns = 'p', y.position = c(0.9, 0.7, 0.5, 7, 6.5, 6, 1.2, 1)
  ) +
  facet_grid(rows = vars(factor(features, levels = c(
    'hydropathy', 'NCPR', 'kappa'
  ))), scales = 'free_y', space = 'free_x') +
  theme(
    panel.grid.major = element_line(color = "gray", linewidth = 0.2),
    panel.grid.minor = element_line(color = "gray", linewidth = 0.1),
    axis.text.x = element_text(family = 'arial', color = 'black', angle = 30, hjust = 1, size = 8),
    axis.text.y = element_text(family = 'arial', color = 'black', size = 8),
  )

ggsave(
  filename = 'figures/figure6/boxplot_HEK_Nterm_13mer_sequence_features.eps',
  plot = boxplot_HEK_Nterm_13mer_sequence_features,
  height = 5, width = 1.8, units = 'in'
)

### figure 6C, D and E, heatmap of 15-mer around N-terminus in different categories
## import the result from plogo 
## (https://plogo.uconn.edu/#, foreground: 13-mer sequence, background: protein-human)
library(tidyverse)
library(ComplexHeatmap)
library(circlize)

# Fast-turnover
fast_turnover_15mer_matrix <- read_delim(
  'figures/figure6/fast_turnover_15mer_matrix.txt',
  col_names = FALSE, delim = '\n'
) |> 
  mutate(extracted = str_extract(X1, "\\[(.*?)\\]")) %>%
  mutate(extracted = str_replace_all(extracted, "\\[|\\]", "")) %>%
  separate(extracted, into = c("amino_acid_position", "value"), sep = "=") %>%
  mutate(value = as.numeric(value)) |> 
  separate(amino_acid_position, into = c("amino_acid", "position"), sep = "_") |> 
  mutate(
    position = as.numeric(position)
  ) |>
  filter(!amino_acid %in% c('U', 'Z', 'B', 'X')) |> 
  select(!X1) |> 
  pivot_wider(names_from = amino_acid, values_from = value)

fast_turnover_15mer_matrix_data <- fast_turnover_15mer_matrix |> 
  column_to_rownames("position") |> 
  as.matrix()

# Median
median_15mer_matrix <- read_delim(
  'figures/figure6/median_15mer_matrix.txt',
  col_names = FALSE, delim = '\n'
) |> 
  mutate(extracted = str_extract(X1, "\\[(.*?)\\]")) %>%
  mutate(extracted = str_replace_all(extracted, "\\[|\\]", "")) %>%
  separate(extracted, into = c("amino_acid_position", "value"), sep = "=") %>%
  mutate(value = as.numeric(value)) |> 
  separate(amino_acid_position, into = c("amino_acid", "position"), sep = "_") |> 
  mutate(
    position = as.numeric(position)
  ) |>
  filter(!amino_acid %in% c('U', 'Z', 'B', 'X')) |> 
  select(!X1) |> 
  pivot_wider(names_from = amino_acid, values_from = value)

median_15mer_matrix_data <- median_15mer_matrix |> 
  column_to_rownames("position") |> 
  as.matrix()

# Stable
stable_15mer_matrix <- read_delim(
  'figures/figure6/stable_15mer_matrix.txt',
  col_names = FALSE, delim = '\n'
) |> 
  mutate(extracted = str_extract(X1, "\\[(.*?)\\]")) %>%
  mutate(extracted = str_replace_all(extracted, "\\[|\\]", "")) %>%
  separate(extracted, into = c("amino_acid_position", "value"), sep = "=") %>%
  mutate(value = as.numeric(value)) |> 
  separate(amino_acid_position, into = c("amino_acid", "position"), sep = "_") |> 
  mutate(
    position = as.numeric(position)
  ) |>
  filter(!amino_acid %in% c('U', 'Z', 'B', 'X')) |> 
  select(!X1) |> 
  pivot_wider(names_from = amino_acid, values_from = value)

stable_15mer_matrix_data <- stable_15mer_matrix |> 
  column_to_rownames("position") |> 
  as.matrix()

# define the amino acid groups for clustering
aa_groups <- list(
  hydrophobic = c("L", "V", "I", "M", "C"),
  small = c("A", "S", "G", "T", "P"),
  aromatic = c("F", "Y", "W"),
  acidic_amide = c("E", "D", "N", "Q"),
  basic = c("K", "R", "H")
)

# flatten the list to get the order of all amino acids
custom_col_order <- unlist(aa_groups)

# subset and reorder the matrix based on our custom order
stable_15mer_matrix_data <- stable_15mer_matrix_data[, custom_col_order]

# create a named vector to map each amino acid to its group
aa_to_group <- rep(names(aa_groups), sapply(aa_groups, length))
names(aa_to_group) <- unlist(aa_groups)

# create split information for columns to add borders around groups
column_split <- factor(aa_to_group[colnames(stable_15mer_matrix_data)], 
                       levels = names(aa_groups))

# create a color scale from blue (negative values) to white (zero) to red (positive values)
col_fun <- colorRamp2(
  breaks = c(-50, 0, 100),
  colors = c("blue", "white", "red")
)

# heatmap
Heatmap(
  stable_15mer_matrix_data,
  name = "Z-score",
  
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  
  show_row_names = TRUE,
  show_column_names = TRUE,
  
  column_split = column_split,
  border = TRUE,
  column_gap = unit(2, "mm"),
  border_gp = gpar(col = "black", lwd = 1),
  column_title_gp = gpar(fontsize = 0),
  
  col = col_fun,
  rect_gp = gpar(col = "white", lwd = 0.5)
)

