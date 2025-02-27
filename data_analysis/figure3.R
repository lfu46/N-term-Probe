# import packages
library(tidyverse)

# generate unique N-terminal proteoform identified in each cell line
unique_Nterm_HEK293T <- setdiff(
  Nterm_list_comb$HEK293T, c(Nterm_list_comb$Jurkat, Nterm_list_comb$`THP-1`)
) |> 
  as_tibble() |> 
  separate(
    value, into = c('UniProt_Accession', 'start.position'), 
    sep = '_', convert = TRUE, remove = FALSE
  ) |> 
  select(Index = value, UniProt_Accession, start.position)

unique_Nterm_Jurkat <- setdiff(
  Nterm_list_comb$Jurkat, c(Nterm_list_comb$`THP-1`, Nterm_list_comb$HEK293T)
) |> 
  as_tibble() |> 
  separate(
    value, into = c('UniProt_Accession', 'start.position'), 
    sep = '_', convert = TRUE, remove = FALSE
  ) |> 
  select(Index = value, UniProt_Accession, start.position)

unique_Nterm_THP1 <- setdiff(
  Nterm_list_comb$`THP-1`, c(Nterm_list_comb$HEK293T, Nterm_list_comb$Jurkat)
) |> 
  as_tibble() |> 
  separate(
    value, into = c('UniProt_Accession', 'start.position'), 
    sep = '_', convert = TRUE, remove = FALSE
  ) |> 
  select(Index = value, UniProt_Accession, start.position)

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

### figure 3A, N-terminal proteoform topFinder information
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


### figure 3B, GO and KEGG analysis of identified unique N-terminal proteoform in HEK293T, Jurkat and THP-1 cells
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

# Jurkat
unique_Nterm_Jurkat_GO_KEGG <- bind_rows(
  unique_Nterm_Jurkat_protein_GO |> 
    filter(Description %in% c(
      'NADH metabolic process',
      'NAD catabolic process',
      'canonical glycolysis',
      'glucose catabolic process to pyruvate',
      'ATP-dependent protein folding chaperone',
      'NADH regeneration',
      'nucleotide metabolic process'
    )), 
  unique_Nterm_Jurkat_protein_KEGG |> 
    filter(Description %in% c(
      'Spliceosome',
      'Proteasome'
    ))
)

barplot_unique_Nterm_Jurkat_GO_KEGG <- unique_Nterm_Jurkat_GO_KEGG |> 
  ggplot() +
  geom_bar(
    aes(
      x = fct_reorder(Description, -log10(p.adjust)), 
      y = -log10(p.adjust)
    ),
    fill = color_2, color = 'transparent', stat = 'identity'
  ) +
  labs(x = 'GO & KEGG', y = '-log10(adjust P value)') +
  coord_flip() +
  theme(
    axis.text = element_text(size = 8, color = 'black')
  )

ggsave(
  filename = 'figures/figure3/barplot_unique_Nterm_Jurkat_GO_KEGG.eps',
  plot = barplot_unique_Nterm_Jurkat_GO_KEGG,
  height = 2, width = 3.5, units = 'in'
)

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

### figure 3C, protein domain and complex analysis of identified unique N-terminal proteoform in HEK293T, Jurkat and THP-1 cells
# import protein domain information from Pfam 
# (https://ftp.ebi.ac.uk/pub/databases/Pfam/releases/Pfam37.0/proteomes/, 9606.tsv.gz, 2024-05-28 13:25	2.8M)
human_pfam <- read_tsv(
  'data_source/pfam/9606.tsv',
  skip = 3, 
  col_names = FALSE
) |> 
  select(X6, X1)

colnames(human_pfam) <- c('TERM', 'GENE')

# protein domain enrichment analysis
library(clusterProfiler)

# HEK293T
unique_Nterm_HEK293T_protein_pfam <- enricher(
  gene = unique_Nterm_HEK293T |> distinct(UniProt_Accession) |> pull(),
  pvalueCutoff = 1,
  universe = total_Nterm_comb |> distinct(UniProt_Accession) |> pull(),
  qvalueCutoff = 1,
  TERM2GENE = human_pfam
)

write_csv(
  unique_Nterm_HEK293T_protein_pfam@result,
  'data_source/unique_Nterm/unique_Nterm_HEK293T_protein_pfam.csv'
)

# Jurkat
unique_Nterm_Jurkat_protein_pfam <- enricher(
  gene = unique_Nterm_Jurkat |> distinct(UniProt_Accession) |> pull(),
  pvalueCutoff = 1,
  universe = total_Nterm_comb |> distinct(UniProt_Accession) |> pull(),
  qvalueCutoff = 1,
  TERM2GENE = human_pfam
)

write_csv(
  unique_Nterm_Jurkat_protein_pfam@result,
  'data_source/unique_Nterm/unique_Nterm_Jurkat_protein_pfam.csv'
)

# THP-1
unique_Nterm_THP1_protein_pfam <- enricher(
  gene = unique_Nterm_THP1 |> distinct(UniProt_Accession) |> pull(),
  pvalueCutoff = 1,
  universe = total_Nterm_comb |> distinct(UniProt_Accession) |> pull(),
  qvalueCutoff = 1,
  TERM2GENE = human_pfam
)

write_csv(
  unique_Nterm_THP1_protein_pfam@result,
  'data_source/unique_Nterm/unique_Nterm_THP1_protein_pfam.csv'
)

# import protein complex information from CORUM
# (https://mips.helmholtz-muenchen.de/corum/download, UniProt-CORUM Mapping, Corum 5.1 release, 2025-01-21)
human_corum <- read_delim(
  'data_source/corum/corum_uniprotCorumMapping.txt',
  col_names = TRUE,
  name_repair = 'universal'
) |> 
  mutate(
    corum_id = paste('corum_id', corum_id, sep = '_')
  ) |> 
  select(corum_id, UniProtKB_accession_number)

colnames(human_corum) <- c('TERM', 'GENE')

# protein complex enrichment analysis
library(clusterProfiler)

# HEK293T
unique_Nterm_HEK293T_protein_corum <- enricher(
  gene = unique_Nterm_HEK293T |> distinct(UniProt_Accession) |> pull(),
  pvalueCutoff = 1,
  universe = total_Nterm_comb |> distinct(UniProt_Accession) |> pull(),
  qvalueCutoff = 1,
  TERM2GENE = human_corum
)

write_csv(
  unique_Nterm_HEK293T_protein_corum@result,
  'data_source/unique_Nterm/unique_Nterm_HEK293T_protein_corum.csv'
)

# Jurkat
unique_Nterm_Jurkat_protein_corum <- enricher(
  gene = unique_Nterm_Jurkat |> distinct(UniProt_Accession) |> pull(),
  pvalueCutoff = 1,
  universe = total_Nterm_comb |> distinct(UniProt_Accession) |> pull(),
  qvalueCutoff = 1,
  TERM2GENE = human_corum
)

write_csv(
  unique_Nterm_Jurkat_protein_corum@result,
  'data_source/unique_Nterm/unique_Nterm_Jurkat_protein_corum.csv'
)

# THP-1
unique_Nterm_THP1_protein_corum <- enricher(
  gene = unique_Nterm_THP1 |> distinct(UniProt_Accession) |> pull(),
  pvalueCutoff = 1,
  universe = total_Nterm_comb |> distinct(UniProt_Accession) |> pull(),
  qvalueCutoff = 1,
  TERM2GENE = human_corum
)

write_csv(
  unique_Nterm_THP1_protein_corum@result,
  'data_source/unique_Nterm/unique_Nterm_THP1_protein_corum.csv'
)

# import results from protein domain and protein complex enrichment analysis
unique_Nterm_HEK293T_protein_corum <- read_csv(
  'data_source/unique_Nterm/unique_Nterm_HEK293T_protein_corum.csv'
)

unique_Nterm_HEK293T_protein_pfam <- read_csv(
  'data_source/unique_Nterm/unique_Nterm_HEK293T_protein_pfam.csv'
)

unique_Nterm_Jurkat_protein_corum <- read_csv(
  'data_source/unique_Nterm/unique_Nterm_Jurkat_protein_corum.csv'
)

unique_Nterm_Jurkat_protein_pfam <- read_csv(
  'data_source/unique_Nterm/unique_Nterm_Jurkat_protein_pfam.csv'
)

unique_Nterm_THP1_protein_corum <- read_csv(
  'data_source/unique_Nterm/unique_Nterm_THP1_protein_corum.csv'
)

unique_Nterm_THP1_protein_pfam <- read_csv(
  'data_source/unique_Nterm/unique_Nterm_THP1_protein_pfam.csv'
)

# generate combine reults from each cell line
unique_Nterm_comb_pfam_corum <- bind_rows(
  # HEK293T
  unique_Nterm_HEK293T_protein_corum |> 
    filter(Description %in% c(
      'corum_id_5199'
    )) |> 
    mutate(
      cell = 'HEK293T',
      category = 'CORUM'
    ),
  
  unique_Nterm_HEK293T_protein_pfam |> 
    filter(Description %in% c(
      'PF02037',
      'PF01585'
    )) |> 
    mutate(
      cell = 'HEK293T',
      category = 'Pfam'
    ),
  
  # Jurkat
  unique_Nterm_Jurkat_protein_corum |> 
    filter(Description %in% c(
      'corum_id_8850',
      'corum_id_181'
    )) |> 
    mutate(
      cell = 'Jurkat',
      category = 'CORUM'
    ),
  
  unique_Nterm_Jurkat_protein_pfam |> 
    filter(Description %in% c(
      'PF00227'
    )) |> 
    mutate(
      cell = 'Jurkat',
      category = 'Pfam'
    ),
  
  # THP-1
  unique_Nterm_THP1_protein_corum |> 
    filter(Description %in% c(
      'corum_id_193',
      'corum_id_3055'
    )) |> 
    mutate(
      cell = 'THP-1',
      category = 'CORUM'
    ),
  
  unique_Nterm_THP1_protein_pfam |> 
    filter(Description %in% c(
      'PF00012'
    )) |> 
    mutate(
      cell = 'THP-1',
      category = 'Pfam'
    )
)

# bar plot
barplot_unique_Nterm_comb_pfam_corum <- unique_Nterm_comb_pfam_corum |> 
  ggplot() +
  geom_bar(
    aes(
      x = factor(Description, Description), 
      y = -log10(p.adjust),
      fill = cell,
      color = category
    ),
    stat = 'identity'
  ) +
  labs(x = '', y = '-log10(adjust P value)') +
  scale_fill_manual(
    values = c(
      'HEK293T' = color_1,
      'Jurkat' = color_2,
      'THP-1' = color_3
    )
  ) +
  scale_color_manual(
    values = c(
      'Pfam' = color_4,
      'CORUM' = color_6
    )
  ) +
  coord_flip() +
  theme(
    axis.text = element_text(size = 8, color = 'black'),
    legend.text = element_text(size = 8, color = 'black'),
    axis.title = element_text(size = 8, color = 'black')
  )

ggsave(
  filename = 'figures/figure3/barplot_unique_Nterm_comb_pfam_corum.eps',
  plot = barplot_unique_Nterm_comb_pfam_corum,
  height = 2, width = 3.5, units = 'in'
)

### figure 3D, Nterm structure features
library(reticulate)

# use specific virtual env created by anaconda
use_condaenv(
  condaenv = '/opt/anaconda3/envs/structuremap',
  required = TRUE
)

# execute the python script for Nterm structure
source_python("data_analysis/Nterm_structuremap_common_unique.py")


