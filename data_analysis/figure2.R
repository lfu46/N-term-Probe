# import packages
library(tidyverse)

### figure 2A, UpSet plot for HEK293T, Jurkat and THP-1
library(ComplexHeatmap)

# generate the list which contains all identified unique N-term in each cell line
Nterm_list_comb <- list(
  'HEK293T' = c(
    overlap_HEK293T, 
    HEK293T_1_07102024_N_term_unique_filtered$Index,
    HEK293T_2_07102024_N_term_unique_filtered$Index
  ),
  'Jurkat' = Jurkat_07102024_N_term$Index,
  'THP-1' = THP1_07102024_N_term$Index
)

# make the combination matrix
comb_mat <- make_comb_mat(Nterm_list_comb)

# make UpSet plot
UpSet(
  comb_mat,
  set_order = c('HEK293T', 'Jurkat', 'THP-1'),
  comb_order = order(comb_size(comb_mat)),
  comb_col = c(color_1, color_2, color_3)[comb_degree(comb_mat)],
  pt_size = unit(20, 'pt'), 
  lwd = 3
)

### figure 2B, GO analysis of commonly identified N-term in HEK293T, Jurkat and THP-1 cells
# generate the list which contains all identified unique N-term in each cell line
Nterm_list_comb <- list(
  'HEK293T' = c(
    overlap_HEK293T, 
    HEK293T_1_07102024_N_term_unique_filtered$Index,
    HEK293T_2_07102024_N_term_unique_filtered$Index
  ),
  'Jurkat' = Jurkat_07102024_N_term$Index,
  'THP-1' = THP1_07102024_N_term$Index
)

# generate the list which contains overlapping N-term from three cell lines
common_Nterm <- tibble(Reduce(intersect, Nterm_list_comb))
colnames(common_Nterm) <- 'common_Nterm'

# generate the list which contains proteins of N-terminal proteoforms from three cell lines
common_Nterm_protein <- common_Nterm |> 
  separate(common_Nterm, into = c('UniProt_Accession', 'start.position'), sep = '_', convert = TRUE) |> 
  distinct(UniProt_Accession)

# GO analysis
library(clusterProfiler)
library(org.Hs.eg.db)

common_Nterm_protein_GO <- enrichGO(
  gene = common_Nterm_protein$UniProt_Accession,
  OrgDb = org.Hs.eg.db,
  keyType = 'UNIPROT',
  ont = 'ALL',
  pvalueCutoff = 1,
  qvalueCutoff = 1
)

write_csv(
  common_Nterm_protein_GO@result,
  file = 'data_source/common_Nterm/GO_common_Nterm.csv'
)

# import result from GO enrichment analysis
GO_common_Nterm <- read_csv(
  'data_source/common_Nterm/GO_common_Nterm.csv'
)

### figure 2B
## bar plot
# BP
barplot_GO_BP_common_Nterm <- GO_common_Nterm |> 
  filter(Description %in% c(
    'chaperone-mediated protein folding',
    "'de novo' protein folding",
    'positive regulation of RNA splicing',
    'positive regulation of mRNA processing',
    'telomere maintenance via telomerase',
    'RNA localization to Cajal body',
    'regulation of protein localization to Cajal body',
    'cell redox homeostasis',
    'U2-type prespliceosome assembly',
    'protein localization to nuclear body'
  )) |> 
  ggplot() +
  geom_bar(
    aes(
      x = fct_reorder(Description, -log10(p.adjust)), 
      y = -log10(p.adjust)
    ),
    fill = color_1, color = 'transparent', stat = 'identity'
  ) +
  labs(x = 'Biological Process', y = '-log10(adjust P value)') +
  coord_flip() +
  theme(
    axis.text = element_text(size = 8, color = 'black')
  )

ggsave(
  filename = 'figures/figure2/barplot_GO_BP_common_Nterm.eps',
  plot = barplot_GO_BP_common_Nterm,
  height = 2, width = 4, units = 'in'
)

### figure 2C, KEGG analysis of commonly identified N-term in HEK293T, Jurkat and THP-1 cells
# generate the list which contains overlapping N-term from three cell lines
common_Nterm <- tibble(Reduce(intersect, Nterm_list_comb))
colnames(common_Nterm) <- 'common_Nterm'

# generate the list which contains proteins of N-terminal proteoforms from three cell lines
common_Nterm_protein <- common_Nterm |> 
  separate(common_Nterm, into = c('UniProt_Accession', 'start.position'), sep = '_', convert = TRUE) |> 
  distinct(UniProt_Accession)

# KEGG analysis
library(clusterProfiler)

common_Nterm_protein_KEGG <- enrichKEGG(
  gene = common_Nterm_protein$UniProt_Accession,
  organism = "hsa",
  keyType = 'uniprot',
  pvalueCutoff = 1,
  qvalueCutoff = 1
)

write_csv(
  common_Nterm_protein_KEGG@result,
  'data_source/common_Nterm/KEGG_common_Nterm.csv'
)

# import result from KEGG enrichment analysis
KEGG_common_Nterm <- read_csv(
  'data_source/common_Nterm/KEGG_common_Nterm.csv'
)

# bar plot
barplot_KEGG_common_Nterm <- KEGG_common_Nterm |> 
  filter(Description %in% c(
    'Spliceosome',
    'Ribosome',
    'Carbon metabolism',
    'Protein processing in endoplasmic reticulum',
    'Glycolysis / Gluconeogenesis',
    'Oxidative phosphorylation',
    'Pentose phosphate pathway',
    'RNA degradation',
    'Proteasome',
    'mRNA surveillance pathway'
  )) |> 
  ggplot() +
  geom_bar(
    aes(
      x = fct_reorder(Description, -log10(p.adjust)), 
      y = -log10(p.adjust)
    ),
    fill = color_4, color = 'transparent', stat = 'identity'
  ) +
  labs(x = 'KEGG pathway', y = '-log10(adjust P value)') +
  coord_flip() +
  theme(
    axis.text = element_text(size = 8, color = 'black')
  )

ggsave(
  filename = 'figures/figure2/barplot_KEGG_common_Nterm.eps',
  plot = barplot_KEGG_common_Nterm,
  height = 2, width = 4, units = 'in'
)

### figure 2D, N-termial proteoform topFinder information
# TopFIND ExploreR retrieves general and position specific information 
# for a list of protein termini as well as protease specific analysis tools. 
# (https://topfind.clip.msl.ubc.ca/topfinder)
# generate the list which contains all identified unique N-term in each cell line
Nterm_list_comb <- list(
  'HEK293T' = c(
    overlap_HEK293T, 
    HEK293T_1_07102024_N_term_unique_filtered$Index,
    HEK293T_2_07102024_N_term_unique_filtered$Index
  ),
  'Jurkat' = Jurkat_07102024_N_term$Index,
  'THP-1' = THP1_07102024_N_term$Index
)

# generate the list which contains overlapping N-term from three cell lines
common_Nterm <- tibble(Reduce(intersect, Nterm_list_comb))
colnames(common_Nterm) <- 'common_Nterm'

# generate the list which contains proteins of N-terminal proteoforms from three cell lines
common_Nterm_protein <- common_Nterm |> 
  separate(
    common_Nterm, 
    into = c('UniProt_Accession', 'start.position'), 
    sep = '_', 
    convert = TRUE,
    remove = FALSE
  )

write_csv(
  common_Nterm_protein,
  'data_source/common_Nterm/common_Nterm_protein.csv'
)

## generate Nterm sequence
library(Biostrings)

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

# generate topfinder id
topfinder_id_common_Nterm <- common_Nterm_protein |> 
  left_join(human_fasta_tibble, by = 'UniProt_Accession') |> 
  mutate(
    Nterm_sequence = substr(Sequence, start = start.position, stop = Full_Protein_Length)
  ) |> 
  mutate(
    UniProt_Accession_Nterm_sequence = paste(UniProt_Accession, Nterm_sequence, sep = ' ')
  ) 

writeLines(
  topfinder_id_common_Nterm$UniProt_Accession_Nterm_sequence,
  'data_source/common_Nterm/topfinder_id_common_Nterm.txt'
)

## import topfinder result
library(tidyverse)

common_Nterm_topfinder <- read_delim(
  'data_source/common_Nterm/2025_02_26_common_Nterm_02262025/2025_02_26_common_Nterm_02262025_Full_Table.txt',
  col_names = TRUE,
  name_repair = 'universal'
)

# database feature
common_Nterm_database_feature <- common_Nterm_topfinder |> 
  select(
    UniProt.curated.start, 
    Alternative.Spliced.Start,
    Cleaving.proteases,
    Other.experimental.terminus.evidences,
    Alternative.Translation.Start
) |> 
  pivot_longer(UniProt.curated.start:Alternative.Translation.Start, names_to = 'category', values_to = 'value') |> 
  filter(!is.na(value)) |> 
  count(category) |> 
  mutate(
    category = factor(category, levels = c(
      'UniProt.curated.start',
      'Alternative.Spliced.Start',
      'Cleaving.proteases',
      'Other.experimental.terminus.evidences',
      'Alternative.Translation.Start'
    ))
  )

# bar plot
barplot_common_Nterm_database_feature <- common_Nterm_database_feature |> 
  ggplot() +
  geom_bar(
    aes(
      x = category, 
      y = n
    ), 
    stat = 'identity', 
    fill = color_1
  ) + 
  labs(x= '', y = '# of N-terminal Proteoform') +
  theme(
    axis.text.x = element_blank(),
    axis.text.y = element_text(size = 8, color = "black"),
  )

ggsave(
  filename = 'figures/figure2/barplot_common_Nterm_database_feature.eps',
  plot = barplot_common_Nterm_database_feature,
  height = 2, width = 2, units = 'in'
)

### figure 2E, enriched cleaving proteases
# Protease_histogram.svg, q < 0.05
# MEP1A, CATL1, MEP1B, CATS, MMP11
common_Nterm_cleaving_proteases <- common_Nterm_topfinder |> 
  select(
    Cleaving.proteases
  ) |> 
  separate_rows(Cleaving.proteases) |> 
  filter(!is.na(Cleaving.proteases)) |> 
  count(Cleaving.proteases) |> 
  filter(
    Cleaving.proteases %in% c(
      'MEP1A', 'MEP1B', 'CATL1', 'CATS', 'MMP11'
    )
  )

# generate data frame for the five enriched proteases
common_Nterm_enriched_cleaving_proteases <- common_Nterm_topfinder |> 
  separate_rows(
    Cleaving.proteases, sep = ';'
  ) |> 
  filter(
    Cleaving.proteases %in% c(
      'MEP1A', 'MEP1B', 'CATL1', 'CATS', 'MMP11'
    )
  ) |> 
  select(
    Accession, Recommended.Protein.Name, Other.Names.and.IDs, Cleaving.proteases
  )

write_csv(
  common_Nterm_enriched_cleaving_proteases,
  file = 'data_source/common_Nterm/common_Nterm_enriched_cleaving_proteases.csv'
)

# bar plot
barplot_common_Nterm_cleaving_proteases <- common_Nterm_cleaving_proteases |> 
  ggplot() +
  geom_bar(
    aes(
      x = fct_reorder(Cleaving.proteases, n), 
      y = n
    ), 
    fill = color_2, 
    stat = 'identity'
  ) +
  labs(x= '', y = '# of N-terminal Proteoform') +
  theme(
    axis.text.x = element_text(size = 8, color = "black", angle = 30, hjust = 1),
    axis.text.y = element_text(size = 8, color = "black"),
  )

ggsave(
  filename = 'figures/figure2/barplot_common_Nterm_cleaving_proteases.eps',
  plot = barplot_common_Nterm_cleaving_proteases,
  height = 2, width = 2, units = 'in'
)

### figure 2F, enriched cleaving proteases example
library(tidyverse)

## MEP1A/MEP1B
# P07900, HSP90AA1, Heat shock protein HSP 90-alpha
P07900_database_info <- tribble(
  ~ name, ~ start, ~ end,
  'protein', 1, 732,
  'HSP90', 196, 714,
  'HATPase_c_3', 43, 159
)

P07900_result <- tibble(
  cleavage_site = c(422, 632, 633, 634, 636, 637, 638)
)

# example plot
P07900_example <- ggplot() +
  geom_rect(
    data = P07900_database_info,
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
      x = P07900_result$cleavage_site, 
      xend = P07900_result$cleavage_site, 
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
      'HSP90' = color_1,
      'HATPase_c_3' = color_2
    )
  ) +
  scale_color_manual(
    values = c(
      'protein' = 'black',
      'HSP90' = 'transparent',
      'HATPase_c_3' = 'transparent'
    )
  ) +
  theme_void()

ggsave(
  filename = 'figures/figure2/P07900_example.eps',
  height = 0.2, width = 2, units = 'in'
)

# P08238, HSP90AB1, Heat shock protein HSP 90-beta
library(tidyverse)

P08238_database_info <- tribble(
  ~ name, ~ start, ~ end,
  'protein', 1, 724,
  'HSP90', 191, 702,
  'HATPase_c_3', 39, 154
)

P08238_result <- tibble(
  cleavage_site = c(292, 492, 624, 625, 626, 627, 628, 629, 630)
)

# example plot
P08238_example <- ggplot() +
  geom_rect(
    data = P08238_database_info,
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
      x = P08238_result$cleavage_site, 
      xend = P08238_result$cleavage_site, 
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
      'HSP90' = color_3,
      'HATPase_c_3' = color_4
    )
  ) +
  scale_color_manual(
    values = c(
      'protein' = 'black',
      'HSP90' = 'transparent',
      'HATPase_c_3' = 'transparent'
    )
  ) +
  theme_void()

ggsave(
  filename = 'figures/figure2/P08238_example.eps',
  height = 0.2, width = 2, units = 'in'
)

## CATL1
# P24534, EEF1B2, Elongation factor 1-beta
library(tidyverse)

P24534_database_info <- tribble(
  ~ name, ~ start, ~ end,
  'protein', 1, 225,
  'EF1_GNE', 142, 225,
  'EF-1_beta_acid', 103, 130
)

P24534_result <- tibble(
  cleavage_site = c(143, 144, 145)
)

# example plot
P24534_example <- ggplot() +
  geom_rect(
    data = P24534_database_info,
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
      x = P24534_result$cleavage_site, 
      xend = P24534_result$cleavage_site, 
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
      'EF1_GNE' = color_1,
      'EF-1_beta_acid' = color_2
    )
  ) +
  scale_color_manual(
    values = c(
      'protein' = 'black',
      'EF1_GNE' = 'transparent',
      'EF-1_beta_acid' = 'transparent'
    )
  ) +
  theme_void()

ggsave(
  filename = 'figures/figure2/P24534_example.eps',
  height = 0.2, width = 2, units = 'in'
)

# P29692, EEF1D, Elongation factor 1-delta
library(tidyverse)

P29692_database_info <- tribble(
  ~ name, ~ start, ~ end,
  'protein', 1, 281,
  'EF1_GNE', 198, 281,
  'EF-1_beta_acid', 159, 186
)

P29692_result <- tibble(
  cleavage_site = c(60, 61, 62)
)

# example plot
P29692_example <- ggplot() +
  geom_rect(
    data = P29692_database_info,
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
      x = P29692_result$cleavage_site, 
      xend = P29692_result$cleavage_site, 
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
      'EF1_GNE' = color_3,
      'EF-1_beta_acid' = color_4
    )
  ) +
  scale_color_manual(
    values = c(
      'protein' = 'black',
      'EF1_GNE' = 'transparent',
      'EF-1_beta_acid' = 'transparent'
    )
  ) +
  theme_void()

ggsave(
  filename = 'figures/figure2/P29692_example.eps',
  height = 0.2, width = 2, units = 'in'
)

# P31946, YWHAB, 14-3-3 protein beta/alpha
library(tidyverse)

P31946_database_info <- tribble(
  ~ name, ~ start, ~ end,
  'protein', 1, 246,
  '14-3-3', 11, 231
)

P31946_result <- tibble(
  cleavage_site = c(30, 31)
)

# example plot
P31946_example <- ggplot() +
  geom_rect(
    data = P31946_database_info,
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
      x = P31946_result$cleavage_site, 
      xend = P31946_result$cleavage_site, 
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
      '14-3-3' = color_5
    )
  ) +
  scale_color_manual(
    values = c(
      'protein' = 'black',
      '14-3-3' = 'transparent'
    )
  ) +
  theme_void()

ggsave(
  filename = 'figures/figure2/P31946_example.eps',
  height = 0.2, width = 2, units = 'in'
)

### figure 2G, Nterm structure feature
library(reticulate)

# use specific virtual env created by anaconda
use_condaenv(
  condaenv = '/opt/anaconda3/envs/structuremap',
  required = TRUE
)

# execute the python script for Nterm structure
source_python("data_analysis/Nterm_structuremap_common.py")

# generate a table example for Fisher's Exact test
library(gridExtra)

Fisher_Exact_test_table <- tibble(
  'Yes' = c('a', 'b'),
  'No' = c('c', 'd')
)

grid.table(Fisher_Exact_test_table)

## save the result 
# proteoform N-terminus structural information
common_Nterm_alphafold_N_terminus_tb <- tibble(common_Nterm_alphafold_N_terminus) |> 
      filter(common_Nterm != 0) |> 
      select(-protein_number, -UniProt_Accession, -start.position)

write_csv(
  common_Nterm_alphafold_N_terminus_tb,
  file = 'data_source/common_Nterm/common_Nterm_alphafold_N_terminus_tb.csv'
)

# enrichment analysis
write_csv(
  enrichment_N_terminus,
  file = 'data_source/common_Nterm/common_Nterm_enrichment_N_terminus.csv'
)

# calculate the percentage of secondary structure
library(tidyverse)

common_Nterm_alphafold_N_terminus_tb <- read_csv(
  'data_source/common_Nterm/common_Nterm_alphafold_N_terminus_tb.csv'
)

common_Nterm_secondary_structure_percentage <- common_Nterm_alphafold_N_terminus_tb |> 
  count(structure_group) |> 
  mutate(
    percentage = (n / sum(n)) * 100
  )

# circular barplot
circular_barplot_secondary_structure_percentage <- common_Nterm_secondary_structure_percentage |> 
  ggplot() +
  geom_bar(
    aes(
      x = 3, 
      y = percentage, 
      fill = structure_group
    ),
    stat = 'identity', 
    show.legend = FALSE
  ) +
  scale_fill_manual(
    values = c(
      'HELX' = color_1,
      'BEND' = color_2,
      'STRN' = color_3,
      'TURN' = color_4,
      'unstructured' = 'gray70'
    )
  ) +
  xlim(1, 4) +
  labs(x = '', y = '') +
  coord_polar(theta = "y") +
  theme(
    axis.text = element_blank(), 
    axis.ticks = element_blank()
  )

ggsave(
  filename = 'figures/figure2/circular_barplot_secondary_structure_percentage.eps',
  plot = circular_barplot_secondary_structure_percentage,
  height = 2, width = 2, units = 'in'
)

# calculate the percentage of solvent accessiblity
common_Nterm_solvent_accessibility_percentage <- common_Nterm_alphafold_N_terminus_tb |> 
  count(high_acc_5) |> 
  mutate(
    percentage = (n / sum(n)) * 100
  ) |> 
  mutate(
    high_acc_5 = as.character(high_acc_5)
  )

# circular barplot
circular_barplot_solvent_accessibility_percentage <- common_Nterm_solvent_accessibility_percentage |> 
  ggplot() +
  geom_bar(
    aes(
      x = 3, 
      y = percentage, 
      fill = high_acc_5
    ),
    stat = 'identity', 
    show.legend = FALSE
  ) +
  scale_fill_manual(
    values = c(
      '0' = color_3,
      '1' = color_4
    )
  ) +
  xlim(1, 4) +
  labs(x = '', y = '') +
  coord_polar(theta = "y") +
  theme(
    axis.text = element_blank(), 
    axis.ticks = element_blank()
  )

ggsave(
  filename = 'figures/figure2/circular_barplot_solvent_accessibility_percentage.eps',
  plot = circular_barplot_solvent_accessibility_percentage,
  height = 2, width = 2, units = 'in'
)
