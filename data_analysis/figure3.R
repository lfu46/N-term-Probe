# import packages
library(tidyverse)

# generate unique N-terminal proteoform identified in each cell line
# HEK293T
unique_Nterm_HEK293T <- setdiff(
  Nterm_list_comb$HEK293T, c(Nterm_list_comb$Jurkat, Nterm_list_comb$`THP-1`)
) |> 
  as_tibble() |> 
  separate(
    value, into = c('UniProt_Accession', 'start.position'), 
    sep = '_', convert = TRUE, remove = FALSE
  ) |> 
  select(Index = value, UniProt_Accession, start.position)

write_csv(
  unique_Nterm_HEK293T,
  'data_source/unique_Nterm/unique_Nterm_HEK293T.csv'
)

# Jurkat
unique_Nterm_Jurkat <- setdiff(
  Nterm_list_comb$Jurkat, c(Nterm_list_comb$`THP-1`, Nterm_list_comb$HEK293T)
) |> 
  as_tibble() |> 
  separate(
    value, into = c('UniProt_Accession', 'start.position'), 
    sep = '_', convert = TRUE, remove = FALSE
  ) |> 
  select(Index = value, UniProt_Accession, start.position)

write_csv(
  unique_Nterm_Jurkat,
  'data_source/unique_Nterm/unique_Nterm_Jurkat.csv'
)

# THP-1
unique_Nterm_THP1 <- setdiff(
  Nterm_list_comb$`THP-1`, c(Nterm_list_comb$HEK293T, Nterm_list_comb$Jurkat)
) |> 
  as_tibble() |> 
  separate(
    value, into = c('UniProt_Accession', 'start.position'), 
    sep = '_', convert = TRUE, remove = FALSE
  ) |> 
  select(Index = value, UniProt_Accession, start.position)

write_csv(
  unique_Nterm_THP1,
  'data_source/unique_Nterm/unique_Nterm_THP1.csv'
)

# total Nterm protein
total_Nterm_comb <- Nterm_list_comb |> 
  unlist() |> 
  unname() |> 
  unique() |> 
  as_tibble() |> 
  separate(
    value, into = c('UniProt_Accession', 'start.position'), 
    sep = '_', convert = TRUE, remove = FALSE
  ) |> 
  select(Index = value, UniProt_Accession, start.position)

write_csv(
  total_Nterm_comb,
  'data_source/unique_Nterm/total_Nterm_comb.csv'
)

### figure 3A, GO and KEGG analysis of identified unique N-terminal proteoform in HEK293T cells
## GO analysis
library(clusterProfiler)
library(org.Hs.eg.db)

# HEK293T
unique_Nterm_HEK293T_protein_GO <- enrichGO(
  gene = unique_Nterm_HEK293T |> distinct(UniProt_Accession) |> pull(),
  OrgDb = org.Hs.eg.db,
  keyType = 'UNIPROT',
  ont = 'ALL',
  pvalueCutoff = 1,
  qvalueCutoff = 1,
  universe = total_Nterm_comb |> distinct(UniProt_Accession) |> pull()
)

write_csv(
  unique_Nterm_HEK293T_protein_GO@result,
  file = 'data_source/unique_Nterm/unique_Nterm_HEK293T_protein_GO.csv'
)

# Jurkat
unique_Nterm_Jurkat_protein_GO <- enrichGO(
  gene = unique_Nterm_Jurkat |> distinct(UniProt_Accession) |> pull(),
  OrgDb = org.Hs.eg.db,
  keyType = 'UNIPROT',
  ont = 'ALL',
  pvalueCutoff = 1,
  qvalueCutoff = 1,
  universe = total_Nterm_comb |> distinct(UniProt_Accession) |> pull()
)

write_csv(
  unique_Nterm_Jurkat_protein_GO@result,
  file = 'data_source/unique_Nterm/unique_Nterm_Jurkat_protein_GO.csv'
)

# THP1
unique_Nterm_THP1_protein_GO <- enrichGO(
  gene = unique_Nterm_THP1 |> distinct(UniProt_Accession) |> pull(),
  OrgDb = org.Hs.eg.db,
  keyType = 'UNIPROT',
  ont = 'ALL',
  pvalueCutoff = 1,
  qvalueCutoff = 1,
  universe = total_Nterm_comb |> distinct(UniProt_Accession) |> pull()
)

write_csv(
  unique_Nterm_THP1_protein_GO@result,
  file = 'data_source/unique_Nterm/unique_Nterm_THP1_protein_GO.csv'
)

## KEGG analysis
library(clusterProfiler)

# HEK293T
unique_Nterm_HEK293T_protein_KEGG <- enrichKEGG(
  gene = unique_Nterm_HEK293T |> distinct(UniProt_Accession) |> pull(),
  organism = "hsa",
  keyType = 'uniprot',
  pvalueCutoff = 1,
  qvalueCutoff = 1,
  universe = total_Nterm_comb |> distinct(UniProt_Accession) |> pull()
)

write_csv(
  unique_Nterm_HEK293T_protein_KEGG@result,
  'data_source/unique_Nterm/unique_Nterm_HEK293T_protein_KEGG.csv'
)

# Jurkat
unique_Nterm_Jurkat_protein_KEGG <- enrichKEGG(
  gene = unique_Nterm_Jurkat |> distinct(UniProt_Accession) |> pull(),
  organism = "hsa",
  keyType = 'uniprot',
  pvalueCutoff = 1,
  qvalueCutoff = 1,
  universe = total_Nterm_comb |> distinct(UniProt_Accession) |> pull()
)

write_csv(
  unique_Nterm_Jurkat_protein_KEGG@result,
  'data_source/unique_Nterm/unique_Nterm_Jurkat_protein_KEGG.csv'
)

# THP-1
unique_Nterm_THP1_protein_KEGG <- enrichKEGG(
  gene = unique_Nterm_THP1 |> distinct(UniProt_Accession) |> pull(),
  organism = "hsa",
  keyType = 'uniprot',
  pvalueCutoff = 1,
  qvalueCutoff = 1,
  universe = total_Nterm_comb |> distinct(UniProt_Accession) |> pull()
)

write_csv(
  unique_Nterm_THP1_protein_KEGG@result,
  'data_source/unique_Nterm/unique_Nterm_THP1_protein_KEGG.csv'
)

## import results from GO and KEGG enrichment analysis
unique_Nterm_HEK293T_protein_GO <- read_csv(
  'data_source/unique_Nterm/unique_Nterm_HEK293T_protein_GO.csv'
)

unique_Nterm_Jurkat_protein_GO <- read_csv(
  'data_source/unique_Nterm/unique_Nterm_Jurkat_protein_GO.csv'
)

unique_Nterm_THP1_protein_GO <- read_csv(
  'data_source/unique_Nterm/unique_Nterm_THP1_protein_GO.csv'
)

unique_Nterm_HEK293T_protein_KEGG <- read_csv(
  'data_source/unique_Nterm/unique_Nterm_HEK293T_protein_KEGG.csv'
)

unique_Nterm_Jurkat_protein_KEGG <- read_csv(
  'data_source/unique_Nterm/unique_Nterm_Jurkat_protein_KEGG.csv'
)

unique_Nterm_THP1_protein_KEGG <- read_csv(
  'data_source/unique_Nterm/unique_Nterm_THP1_protein_KEGG.csv'
)

## bar plot
# HEK293T
unique_Nterm_HEK293T_GO_KEGG <- bind_rows(
  unique_Nterm_HEK293T_protein_GO |> 
    filter(Description %in% c(
      'chromosome organization',
      'DNA binding',
      'nuclear body',
      'cadherin binding',
      'molecular adaptor activity',
      'supramolecular complex',
      'mRNA splicing, via spliceosome',
      'cytoplasmic stress granule',
      'regulation of cell cycle'
    ))
)

barplot_unique_Nterm_HEK293T_GO_KEGG <- unique_Nterm_HEK293T_GO_KEGG |> 
  ggplot() +
  geom_bar(
    aes(
      x = fct_reorder(Description, -log10(p.adjust)), 
      y = -log10(p.adjust)
    ),
    fill = color_1, color = 'transparent', stat = 'identity'
  ) +
  labs(x = 'GO & KEGG', y = '-log10(adjust P value)') +
  coord_flip() +
  theme(
    axis.text = element_text(size = 8, color = 'black')
  )

ggsave(
  filename = 'figures/figure3/barplot_unique_Nterm_HEK293T_GO_KEGG.eps',
  plot = barplot_unique_Nterm_HEK293T_GO_KEGG,
  height = 2, width = 3, units = 'in'
)

### figure 3B, GO and KEGG analysis of identified unique N-terminal proteoform in HEK293T cells
# THP-1
unique_Nterm_THP1_GO_KEGG <- bind_rows(
  unique_Nterm_THP1_protein_GO |> 
    filter(Description %in% c(
      'secretory granule',
      'vesicle lumen',
      'cytoplasmic vesicle lumen',
      'lysosome',
      'ficolin-1-rich granule',
      'cell surface',
      'lytic vacuole',
      'secretory granule membrane',
      'calcium ion binding'
    ))
)

barplot_unique_Nterm_THP1_GO_KEGG <- unique_Nterm_THP1_GO_KEGG |> 
  ggplot() +
  geom_bar(
    aes(
      x = fct_reorder(Description, -log10(p.adjust)), 
      y = -log10(p.adjust)
    ),
    fill = color_3, color = 'transparent', stat = 'identity'
  ) +
  labs(x = 'GO & KEGG', y = '-log10(adjust P value)') +
  coord_flip() +
  theme(
    axis.text = element_text(size = 8, color = 'black')
  )

ggsave(
  filename = 'figures/figure3/barplot_unique_Nterm_THP1_GO_KEGG.eps',
  plot = barplot_unique_Nterm_THP1_GO_KEGG,
  height = 2, width = 3.5, units = 'in'
)

### figure 3C, enriched cleaving proteases in each cell line
# TopFIND ExploreR retrieves general and position specific information 
# for a list of protein termini as well as protease specific analysis tools. 
# (https://topfind.clip.msl.ubc.ca/topfinder)
# generate Nterm sequence
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

## generate topfinder id
# HEK293T
topfinder_id_unique_Nterm_HEK293T <- unique_Nterm_HEK293T |> 
  left_join(human_fasta_tibble, by = 'UniProt_Accession') |> 
  mutate(
    Nterm_sequence = substr(Sequence, start = start.position, stop = Full_Protein_Length)
  ) |> 
  mutate(
    UniProt_Accession_Nterm_sequence = paste(UniProt_Accession, Nterm_sequence, sep = ' ')
  )

writeLines(
  topfinder_id_unique_Nterm_HEK293T$UniProt_Accession_Nterm_sequence,
  'data_source/unique_Nterm/topfinder_id_unique_Nterm_HEK293T.txt'
)

# Jurkat
topfinder_id_unique_Nterm_Jurkat <- unique_Nterm_Jurkat |> 
  left_join(human_fasta_tibble, by = 'UniProt_Accession') |> 
  mutate(
    Nterm_sequence = substr(Sequence, start = start.position, stop = Full_Protein_Length)
  ) |> 
  mutate(
    UniProt_Accession_Nterm_sequence = paste(UniProt_Accession, Nterm_sequence, sep = ' ')
  )

writeLines(
  topfinder_id_unique_Nterm_Jurkat$UniProt_Accession_Nterm_sequence,
  'data_source/unique_Nterm/topfinder_id_unique_Nterm_Jurkat.txt'
)

# THP-1
topfinder_id_unique_Nterm_THP1 <- unique_Nterm_THP1 |> 
  left_join(human_fasta_tibble, by = 'UniProt_Accession') |> 
  mutate(
    Nterm_sequence = substr(Sequence, start = start.position, stop = Full_Protein_Length)
  ) |> 
  mutate(
    UniProt_Accession_Nterm_sequence = paste(UniProt_Accession, Nterm_sequence, sep = ' ')
  )

writeLines(
  topfinder_id_unique_Nterm_THP1$UniProt_Accession_Nterm_sequence,
  'data_source/unique_Nterm/topfinder_id_unique_Nterm_THP1.txt'
)

## import topfinder results
library(tidyverse)

# HEK293T
unique_Nterm_topfinder_HEK293T <- read_delim(
  'data_source/unique_Nterm/2025_02_27_unique_Nterm_HEK293T_02272025/2025_02_27_unique_Nterm_HEK293T_02272025_Full_Table.txt',
  col_names = TRUE,
  name_repair = 'universal'
)

# Jurkat
unique_Nterm_topfinder_Jurkat <- read_delim(
  'data_source/unique_Nterm/2025_02_27_unique_Nterm_Jurkat_02272025/2025_02_27_unique_Nterm_Jurkat_02272025_Full_Table.txt',
  col_names = TRUE,
  name_repair = 'universal'
)

# THP-1
unique_Nterm_topfinder_THP1 <- read_delim(
  'data_source/unique_Nterm/2025_02_27_unique_Nterm_THP1_02272025/2025_02_27_unique_Nterm_THP1_02272025_Full_Table.txt',
  col_names = TRUE,
  name_repair = 'universal'
)

# enriched cleaving proteases in each cell line
enriched_protease_HEK293T <- unique_Nterm_topfinder_HEK293T |> 
  separate_rows(Cleaving.proteases, sep = ';') |> 
  filter(str_detect(Cleaving.proteases, 'GRAA')) |> 
  count(Cleaving.proteases) |> 
  mutate(
    cell = 'HEK293T'
  )

enriched_protease_Jurkat <- unique_Nterm_topfinder_Jurkat |> 
  separate_rows(Cleaving.proteases, sep = ';') |> 
  filter(str_detect(Cleaving.proteases, 'CATL1|MEP1B|CATS|MEP1A|CATB')) |> 
  count(Cleaving.proteases) |> 
  mutate(
    cell = 'Jurkat'
  )

enriched_protease_THP1 <- unique_Nterm_topfinder_THP1 |> 
  separate_rows(Cleaving.proteases, sep = ';') |> 
  filter(str_detect(Cleaving.proteases, 'HTRA2|CASP1|MMP11|ELNE|CASP7|MMP3')) |> 
  count(Cleaving.proteases) |> 
  mutate(
    cell = 'THP-1'
  )

# output the unique Nterm result for each cleaving protease
enriched_protease_Nterm_comb <- bind_rows(
  unique_Nterm_topfinder_HEK293T |> 
    separate_rows(Cleaving.proteases, sep = ';') |> 
    filter(str_detect(Cleaving.proteases, 'GRAA')) |> 
    select(Accession, Recommended.Protein.Name, Other.Names.and.IDs, P1..Position, Cleaving.proteases) |> 
    mutate(cell = 'HEK293T'),
  
  unique_Nterm_topfinder_Jurkat |> 
    separate_rows(Cleaving.proteases, sep = ';') |> 
    filter(str_detect(Cleaving.proteases, 'CATL1|MEP1B|CATS|MEP1A|CATB')) |> 
    select(Accession, Recommended.Protein.Name, Other.Names.and.IDs, P1..Position, Cleaving.proteases) |> 
    mutate(cell = 'Jurkat'),
  
  unique_Nterm_topfinder_THP1 |> 
    separate_rows(Cleaving.proteases, sep = ';') |> 
    filter(str_detect(Cleaving.proteases, 'HTRA2|CASP1|MMP11|ELNE|CASP7|MMP3')) |> 
    select(Accession, Recommended.Protein.Name, Other.Names.and.IDs, P1..Position, Cleaving.proteases) |> 
    mutate(cell = 'THP-1')
)

write_csv(
  enriched_protease_Nterm_comb,
  file = 'data_source/unique_Nterm/enriched_protease_Nterm_comb.csv'
)

# combine result form HEK293T, Jurkat and THP-1
enriched_protease_comb <- bind_rows(
  enriched_protease_HEK293T,
  enriched_protease_Jurkat,
  enriched_protease_THP1
)

# bar plot
barplot_enriched_protease_comb <- enriched_protease_comb |> 
  ggplot() +
  geom_bar(
    aes(
      x = fct_reorder(Cleaving.proteases, n),
      y = n,
      fill = cell
    ),
    stat = 'identity', show.legend = FALSE
  ) +
  facet_grid(. ~ cell, scales = 'free', space = 'free') +
  labs(x = '', y = 'Cleavages') +
  scale_fill_manual(
    values = c(
      'HEK293T' = color_1,
      'Jurkat' = color_2,
      'THP-1' = color_3
    )
  ) +
  theme(
    axis.text.x = element_text(size = 8, color = 'black', angle = 90, hjust = 1),
    axis.text.y = element_text(size = 8, color = 'black'),
    axis.title = element_text(size = 8, color = 'black'),
    strip.text = element_text(size = 8, color = 'black')
  )

ggsave(
  filename = 'figures/figure3/barplot_enriched_protease_comb.eps',
  plot = barplot_enriched_protease_comb,
  height = 2, width = 2.5, units = 'in'
)

### figure 3D, N-terminus features in each cell line
N_terminus_feature_HEK293T <- unique_Nterm_topfinder_HEK293T |> 
  select(Distance.To.signal.peptide, Distance.to.propeptide.lost, Distance.to.last.transmembrane.domain..shed.) |> 
  pivot_longer(
   cols =  Distance.To.signal.peptide:Distance.to.last.transmembrane.domain..shed., names_to = 'N_terminus_feature', values_to = 'distance'
  ) |> 
  filter(!is.na(distance)) |> 
  count(N_terminus_feature) |> 
  mutate(
    cell = 'HEK293T'
  )

N_terminus_feature_Jurkat <- unique_Nterm_topfinder_Jurkat |> 
  select(Distance.To.signal.peptide, Distance.to.propeptide.lost, Distance.to.last.transmembrane.domain..shed.) |> 
  pivot_longer(
    cols =  Distance.To.signal.peptide:Distance.to.last.transmembrane.domain..shed., names_to = 'N_terminus_feature', values_to = 'distance'
  ) |> 
  filter(!is.na(distance)) |> 
  count(N_terminus_feature) |> 
  mutate(
    cell = 'Jurkat'
  )

N_terminus_feature_THP1 <- unique_Nterm_topfinder_THP1 |> 
  select(Distance.To.signal.peptide, Distance.to.propeptide.lost, Distance.to.last.transmembrane.domain..shed.) |> 
  pivot_longer(
    cols =  Distance.To.signal.peptide:Distance.to.last.transmembrane.domain..shed., names_to = 'N_terminus_feature', values_to = 'distance'
  ) |> 
  filter(!is.na(distance)) |> 
  count(N_terminus_feature) |> 
  mutate(
    cell = 'THP-1'
  )

N_terminus_feature_THP1_noCount <- unique_Nterm_topfinder_THP1 |> 
  select(Accession, Recommended.Protein.Name, Other.Names.and.IDs, P1..Position, Distance.To.signal.peptide, Distance.to.propeptide.lost, Distance.to.last.transmembrane.domain..shed.) |> 
  pivot_longer(
    cols =  Distance.To.signal.peptide:Distance.to.last.transmembrane.domain..shed., names_to = 'N_terminus_feature', values_to = 'distance'
  ) |> 
  filter(!is.na(distance))

write_csv(
  N_terminus_feature_THP1_noCount,
  file = 'data_source/unique_Nterm/N_terminus_feature_THP1_noCount.csv'
)

# combine results from HEK293T, Jurkat and THP-1
N_terminus_feature_comb <- bind_rows(
  N_terminus_feature_HEK293T,
  N_terminus_feature_Jurkat,
  N_terminus_feature_THP1
) |> 
  mutate(
    N_terminus_feature = case_when(
      str_detect(N_terminus_feature, 'signal.peptide') ~ 'Signal.Peptide',
      str_detect(N_terminus_feature, 'transmembrane.domain') ~ 'Transmembrane.Domain',
      str_detect(N_terminus_feature, 'propeptide.lost') ~ 'Propeptide.Lost'
    )
  )

# pairwise chi-square test
N_terminus_feature_comb_matrix_tb <- N_terminus_feature_comb |> 
  pivot_wider(names_from = cell, values_from = n)

N_terminus_feature_comb_matrix <- data.matrix(N_terminus_feature_comb_matrix_tb)
rownames(N_terminus_feature_comb_matrix) <- N_terminus_feature_comb_matrix_tb$N_terminus_feature

# HEK293T vs. Jurkat
chisq.test(N_terminus_feature_comb_matrix[, c('HEK293T', 'Jurkat')])

# Jurkat vs. THP-1
chisq.test(N_terminus_feature_comb_matrix[, c('Jurkat', 'THP-1')])

# THP-1 vs. HEK293T
chisq.test(N_terminus_feature_comb_matrix[, c('THP-1', 'HEK293T')])

# bar plot
barplot_N_terminus_feature_comb <- N_terminus_feature_comb |> 
  ggplot() +
  geom_bar(
    aes(
      x = N_terminus_feature,
      y = n,
      fill = N_terminus_feature
    ),
    stat = 'identity'
  ) +
  facet_wrap(vars(cell), nrow = 1) +
  scale_fill_manual(
    name = '',
    values = c(
      'Signal.Peptide' = color_1,
      'Transmembrane.Domain' = color_2,
      'Propeptide.Lost' = color_3
    )
  ) +
  labs(x = '', y = '') +
  theme(
    axis.text.x = element_blank(),
    axis.text.y = element_text(size = 8, color = 'black'),
    strip.text = element_text(size = 8, color = 'black')
  )

ggsave(
  filename = 'figures/figure3/barplot_N_terminus_feature_comb.eps',
  plot = barplot_N_terminus_feature_comb,
  height = 1.65, width = 4, units = 'in'
)

### figure 3E, protein N-termini cleaving protease examples for HEK293T and Jurkat
library(tidyverse)

## HEK293T, GRAA
# O95793, STAU1, Double-stranded RNA-binding protein Staufen homolog 1
O95793_database_info <- tribble(
  ~ name, ~ start, ~ end,
  'protein', 1, 577,
  'dsrm', 103, 160,
  'dsrm', 186, 249,
  'dsrm', 287, 352,
  'Staufen_C', 447, 557
)

O95793_result <- tibble(
  cleavage_site = c(159, 161)
)

# example plot
O95793_example <- ggplot() +
  geom_rect(
    data = O95793_database_info,
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
      x = O95793_result$cleavage_site, 
      xend = O95793_result$cleavage_site, 
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
      'dsrm' = color_1,
      'Staufen_C' = color_2
    )
  ) +
  scale_color_manual(
    values = c(
      'protein' = 'black',
      'dsrm' = 'transparent',
      'Staufen_C' = 'transparent'
    )
  ) +
  theme_void()

ggsave(
  filename = 'figures/figure3/O95793_example.eps',
  height = 0.2, width = 2, units = 'in'
)

# P10412, H1-4, Histone H1.4
library(tidyverse)

P10412_database_info <- tribble(
  ~ name, ~ start, ~ end,
  'protein', 1, 219,
  'Linker_histone', 38, 108
)

P10412_result <- tibble(
  cleavage_site = c(68)
)

# example plot
P10412_example <- ggplot() +
  geom_rect(
    data = P10412_database_info,
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
      x = P10412_result$cleavage_site, 
      xend = P10412_result$cleavage_site, 
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
      'Linker_histone' = color_3
    )
  ) +
  scale_color_manual(
    values = c(
      'protein' = 'black',
      'Linker_histone' = 'transparent'
    )
  ) +
  theme_void()

ggsave(
  filename = 'figures/figure3/P10412_example.eps',
  height = 0.2, width = 2, units = 'in'
)

# Q93077, H2AC6, Histone H2A type 1-C
library(tidyverse)

Q93077_database_info <- tribble(
  ~ name, ~ start, ~ end,
  'protein', 1, 130,
  'Histone', 5, 89,
  'Histone_H2A_C', 92, 126
)

Q93077_result <- tibble(
  cleavage_site = c(22)
)

# example plot
Q93077_example <- ggplot() +
  geom_rect(
    data = Q93077_database_info,
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
      x = Q93077_result$cleavage_site, 
      xend = Q93077_result$cleavage_site, 
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
      'Histone' = color_4,
      'Histone_H2A_C' = color_5
    )
  ) +
  scale_color_manual(
    values = c(
      'protein' = 'black',
      'Histone' = 'transparent',
      'Histone_H2A_C' = 'transparent'
    )
  ) +
  theme_void()

ggsave(
  filename = 'figures/figure3/Q93077_example.eps',
  height = 0.2, width = 2, units = 'in'
)

## Jurkat, CATB
library(tidyverse)

# O14979, HNRNPDL, Heterogeneous nuclear ribonucleoprotein D-like
O14979_database_info <- tribble(
  ~ name, ~ start, ~ end,
  'protein', 1, 420,
  'RRM_1', 151, 218,
  'RRM_1', 235, 294
)

O14979_result <- tibble(
  cleavage_site = c(43, 236, 237, 278)
)

# example plot
O14979_example <- ggplot() +
  geom_rect(
    data = O14979_database_info,
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
      x = O14979_result$cleavage_site, 
      xend = O14979_result$cleavage_site, 
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
  filename = 'figures/figure3/O14979_example.eps',
  height = 0.2, width = 2, units = 'in'
)

# P07910, HNRNPC, Heterogeneous nuclear ribonucleoproteins C1/C2
P07910_database_info <- tribble(
  ~ name, ~ start, ~ end,
  'protein', 1, 420,
  'RRM_1', 18, 80
)

P07910_result <- tibble(
  cleavage_site = c(76, 80)
)

# example plot
P07910_example <- ggplot() +
  geom_rect(
    data = P07910_database_info,
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
      x = P07910_result$cleavage_site, 
      xend = P07910_result$cleavage_site, 
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
      'RRM_1' = color_2
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
  filename = 'figures/figure3/P07910_example.eps',
  height = 0.2, width = 2, units = 'in'
)

# Q92945, KHSRP, Far upstream element-binding protein 2
Q92945_database_info <- tribble(
  ~ name, ~ start, ~ end,
  'protein', 1, 711,
  'KH_1', 146, 209,
  'KH_1', 236, 301,
  'KH_1', 327, 387,
  'KH_1', 427, 493,
  'FUBP_C', 611, 632,
  'FUBP_C', 670, 685
)

Q92945_result <- tibble(
  cleavage_site = c(269, 271)
)

# example plot
Q92945_example <- ggplot() +
  geom_rect(
    data = Q92945_database_info,
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
      x = Q92945_result$cleavage_site, 
      xend = Q92945_result$cleavage_site, 
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
      'KH_1' = color_3,
      'FUBP_C' = color_4
    )
  ) +
  scale_color_manual(
    values = c(
      'protein' = 'black',
      'KH_1' = 'transparent',
      'FUBP_C' = 'transparent'
    )
  ) +
  theme_void()

ggsave(
  filename = 'figures/figure3/Q92945_example.eps',
  height = 0.2, width = 2, units = 'in'
)

## THP-1, HTRA2
library(tidyverse)

# Q9BQE3, TUBA1C, Tubulin alpha-1C chain
Q9BQE3_database_info <- tribble(
  ~ name, ~ start, ~ end,
  'protein', 1, 449,
  'Tubulin', 3, 213,
  'Tubulin_C', 263, 392
)

Q9BQE3_result <- tibble(
  cleavage_site = c(
    65, 68, 69,
    251, 254,
    294,
    327, 329,
    353, 355, 356,
    405, 406, 407, 408, 409, 411, 414, 415
  )
)

# example plot
Q9BQE3_example <- ggplot() +
  geom_rect(
    data = Q9BQE3_database_info,
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
      x = Q9BQE3_result$cleavage_site, 
      xend = Q9BQE3_result$cleavage_site, 
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
      'Tubulin' = color_1,
      'Tubulin_C' = color_2
    )
  ) +
  scale_color_manual(
    values = c(
      'protein' = 'black',
      'Tubulin' = 'transparent',
      'Tubulin_C' = 'transparent'
    )
  ) +
  theme_void()

ggsave(
  filename = 'figures/figure3/Q9BQE3_example.eps',
  height = 0.2, width = 2, units = 'in'
)

# P07437, TUBB, Tubulin beta chain
library(tidyverse)

P07437_database_info <- tribble(
  ~ name, ~ start, ~ end,
  'protein', 1, 444,
  'Tubulin', 3, 211,
  'Tubulin_C', 261, 382
)

P07437_result <- tibble(
  cleavage_site = c(
    66, 68, 242
  )
)

# example plot
P07437_example <- ggplot() +
  geom_rect(
    data = P07437_database_info,
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
      x = P07437_result$cleavage_site, 
      xend = P07437_result$cleavage_site, 
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
      'Tubulin' = color_3,
      'Tubulin_C' = color_4
    )
  ) +
  scale_color_manual(
    values = c(
      'protein' = 'black',
      'Tubulin' = 'transparent',
      'Tubulin_C' = 'transparent'
    )
  ) +
  theme_void()

ggsave(
  filename = 'figures/figure3/P07437_example.eps',
  height = 0.2, width = 2, units = 'in'
)

## THP-1, CASP1 & CASP7
library(tidyverse)

# P60709, ACTB, Actin, cytoplasmic 1
P60709_database_info <- tribble(
  ~ name, ~ start, ~ end,
  'protein', 1, 375,
  'Actin', 3, 375
)

P60709_result <- tibble(
  cleavage_site = c(
    13, 14, 243, 244
  )
)

# example plot
P60709_example <- ggplot() +
  geom_rect(
    data = P60709_database_info,
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
      x = P60709_result$cleavage_site, 
      xend = P60709_result$cleavage_site, 
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
      'Actin' = color_5
    )
  ) +
  scale_color_manual(
    values = c(
      'protein' = 'black',
      'Actin' = 'transparent'
    )
  ) +
  theme_void()

ggsave(
  filename = 'figures/figure3/P60709_example.eps',
  height = 0.2, width = 2, units = 'in'
)

# P06396, GSN, Gelsolin
library(tidyverse)

P06396_database_info <- tribble(
  ~ name, ~ start, ~ end,
  'protein', 1, 782,
  'Gelsolin', 76, 158,
  'Gelsolin', 197, 270,
  'Gelsolin', 316, 389,
  'Gelsolin', 455, 536,
  'Gelsolin', 577, 642,
  'Gelsolin', 681, 757
)

P06396_result <- tibble(
  cleavage_site = c(
    401, 402, 405, 407, 401, 402, 405, 407
  )
)

# example plot
P06396_example <- ggplot() +
  geom_rect(
    data = P06396_database_info,
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
      x = P06396_result$cleavage_site, 
      xend = P06396_result$cleavage_site, 
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
      'Gelsolin' = color_1
    )
  ) +
  scale_color_manual(
    values = c(
      'protein' = 'black',
      'Gelsolin' = 'transparent'
    )
  ) +
  theme_void()

ggsave(
  filename = 'figures/figure3/P06396_example.eps',
  height = 0.2, width = 2, units = 'in'
)

## THP-1, MMP3 & MMP11
library(tidyverse)

# P60709, ACTB, Actin, cytoplasmic 1
P60709_database_info <- tribble(
  ~ name, ~ start, ~ end,
  'protein', 1, 375,
  'Actin', 3, 375
)

P60709_result <- tibble(
  cleavage_site = c(
    294, 295, 296, 297, 298, 299, 300
  )
)

# example plot
P60709_example_MMP <- ggplot() +
  geom_rect(
    data = P60709_database_info,
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
      x = P60709_result$cleavage_site, 
      xend = P60709_result$cleavage_site, 
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
      'Actin' = color_5
    )
  ) +
  scale_color_manual(
    values = c(
      'protein' = 'black',
      'Actin' = 'transparent'
    )
  ) +
  theme_void()

ggsave(
  filename = 'figures/figure3/P60709_example_MMP.eps',
  height = 0.2, width = 2, units = 'in'
)

# P06396, GSN, Gelsolin
library(tidyverse)

P06396_database_info <- tribble(
  ~ name, ~ start, ~ end,
  'protein', 1, 782,
  'Gelsolin', 76, 158,
  'Gelsolin', 197, 270,
  'Gelsolin', 316, 389,
  'Gelsolin', 455, 536,
  'Gelsolin', 577, 642,
  'Gelsolin', 681, 757
)

P06396_result <- tibble(
  cleavage_site = c(
    51, 52, 53, 415, 420
  )
)

# example plot
P06396_example_MMP <- ggplot() +
  geom_rect(
    data = P06396_database_info,
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
      x = P06396_result$cleavage_site, 
      xend = P06396_result$cleavage_site, 
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
      'Gelsolin' = color_1
    )
  ) +
  scale_color_manual(
    values = c(
      'protein' = 'black',
      'Gelsolin' = 'transparent'
    )
  ) +
  theme_void()

ggsave(
  filename = 'figures/figure3/P06396_example_MMP.eps',
  height = 0.2, width = 2, units = 'in'
)


### figure 3G, Nterm structure features for secondary structure
library(reticulate)

# use specific virtual env created by anaconda
use_condaenv(
  condaenv = '/opt/anaconda3/envs/structuremap',
  required = TRUE
)

# execute the python script for Nterm structure
source_python("data_analysis/Nterm_structuremap_common_unique.py")

## AlphaFold annotation for unique Nterm in each cell line
# HEK293T
unique_Nterm_HEK293T_alphafold_annotation_R <- unique_Nterm_HEK293T |> 
  left_join(
    unique_Nterm_HEK293T_alphafold_accessibility_smooth, by = join_by(
      'UniProt_Accession' == 'protein_id', 'start.position' == 'position'
    )
  )

write_csv(
  unique_Nterm_HEK293T_alphafold_annotation_R,
  'data_source/unique_Nterm/unique_Nterm_HEK293T_alphafold_annotation_R.csv'
)

# Jurkat
unique_Nterm_Jurkat_alphafold_annotation_R <- unique_Nterm_Jurkat |> 
  left_join(
    unique_Nterm_Jurkat_alphafold_accessibility_smooth, by = join_by(
      'UniProt_Accession' == 'protein_id', 'start.position' == 'position'
    )
  )

write_csv(
  unique_Nterm_Jurkat_alphafold_annotation_R,
  'data_source/unique_Nterm/unique_Nterm_Jurkat_alphafold_annotation_R.csv'
)

# THP-1
unique_Nterm_THP1_alphafold_annotation_R <- unique_Nterm_THP1 |> 
  left_join(
    unique_Nterm_THP1_alphafold_accessibility_smooth, by = join_by(
      'UniProt_Accession' == 'protein_id', 'start.position' == 'position'
    )
  )

write_csv(
  unique_Nterm_THP1_alphafold_annotation_R,
  'data_source/unique_Nterm/unique_Nterm_THP1_alphafold_annotation_R.csv'
)

# import the results of AlphaFold annotation from each cell line
library(tidyverse)

# HEK293T 
unique_Nterm_HEK293T_alphafold_annotation_R <- read_csv(
  'data_source/unique_Nterm/unique_Nterm_HEK293T_alphafold_annotation_R.csv'
)

# Jurkat
unique_Nterm_Jurkat_alphafold_annotation_R <- read_csv(
  'data_source/unique_Nterm/unique_Nterm_Jurkat_alphafold_annotation_R.csv'
)

# THP-1
unique_Nterm_THP1_alphafold_annotation_R <- read_csv(
  'data_source/unique_Nterm/unique_Nterm_THP1_alphafold_annotation_R.csv'
)

## calculate the percentage of secondary structure for each cell line
# HEK293T
unique_Nterm_HEK293T_secondary_group_percentage <- unique_Nterm_HEK293T_alphafold_annotation_R |> 
  filter(!is.na(structure_group)) |> 
  count(structure_group) |> 
  mutate(
    percentage = n/sum(n),
    cell = 'HEK293T'
  )

# Jurkat
unique_Nterm_Jurkat_secondary_group_percentage <- unique_Nterm_Jurkat_alphafold_annotation_R |> 
  filter(!is.na(structure_group)) |> 
  count(structure_group) |> 
  mutate(
    percentage = n/sum(n),
    cell = 'Jurkat'
  )

# THP-1
unique_Nterm_THP1_secondary_group_percentage <- unique_Nterm_THP1_alphafold_annotation_R |> 
  filter(!is.na(structure_group)) |> 
  count(structure_group) |> 
  mutate(
    percentage = n/sum(n),
    cell = 'THP-1'
  )

# combine result from HEK293T, Jurkat and THP-1
unique_Nterm_secondary_group_percentage_comb <- bind_rows(
  unique_Nterm_HEK293T_secondary_group_percentage,
  unique_Nterm_Jurkat_secondary_group_percentage,
  unique_Nterm_THP1_secondary_group_percentage
)

write_csv(
  unique_Nterm_secondary_group_percentage_comb,
  file = 'data_source/unique_Nterm/unique_Nterm_secondary_group_percentage_comb.csv'
)

# ANOVA test/Turkey HSD
ANOVA_result <- aov(
  percentage ~ structure_group, data = unique_Nterm_secondary_group_percentage_comb
)
summary(ANOVA_result)
TukeyHSD(ANOVA_result)

# bar plot
barplot_unique_Nterm_secondary_group_percentage_comb <- unique_Nterm_secondary_group_percentage_comb |> 
  ggplot() +
  geom_bar(
    aes(
      x = cell,
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
  coord_polar(theta = "y") +
  theme(
    axis.text.x = element_text(size = 8, color = 'black'),
    axis.text.y = element_text(size = 8, color = 'black'),
    legend.text = element_text(size = 8, color = 'black')
  )

ggsave(
  filename = 'figures/figure3/barplot_unique_Nterm_secondary_group_percentage_comb.eps',
  plot = barplot_unique_Nterm_secondary_group_percentage_comb,
  height = 2, width = 5, units = 'in'
)

### figure 3H, THP-1 signal peptide example
library(tidyverse)

# P49755, TMED10, Transmembrane emp24 domain-containing protein 10
library(tidyverse)

P49755_database_info <- tribble(
  ~ name, ~ start, ~ end,
  'protein', 1, 219,
  'EMP24_GP25L', 31, 213
)

P49755_result <- tibble(
  cleavage_site = c(
    32, 144
  )
)

# example plot
P49755_example <- ggplot() +
  geom_rect(
    data = P49755_database_info,
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
      x = P49755_result$cleavage_site, 
      xend = P49755_result$cleavage_site, 
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
      'EMP24_GP25L' = color_2
    )
  ) +
  scale_color_manual(
    values = c(
      'protein' = 'black',
      'EMP24_GP25L' = 'transparent'
    )
  ) +
  theme_void()

ggsave(
  filename = 'figures/figure3/P49755_example.eps',
  height = 0.2, width = 2, units = 'in'
)

### figure 3I, THP-1 transmembrane domain example
# O75063, FAM20B, Glycosaminoglycan xylosylkinase
library(tidyverse)

O75063_database_info <- tribble(
  ~ name, ~ start, ~ end,
  'protein', 1, 409,
  'Fam20C', 191, 398
)

O75063_result <- tibble(
  cleavage_site = c(
    25
  )
)

# example plot
O75063_example <- ggplot() +
  geom_rect(
    data = O75063_database_info,
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
      x = O75063_result$cleavage_site, 
      xend = O75063_result$cleavage_site, 
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
      'Fam20C' = color_3
    )
  ) +
  scale_color_manual(
    values = c(
      'protein' = 'black',
      'Fam20C' = 'transparent'
    )
  ) +
  theme_void()

ggsave(
  filename = 'figures/figure3/O75063_example.eps',
  height = 0.2, width = 2, units = 'in'
)
