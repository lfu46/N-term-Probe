# import packages
library(tidyverse)

### figure S2A, GO & KEGG analysis for Jarket cell
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
  filename = 'figures/figureS2/barplot_unique_Nterm_Jurkat_GO_KEGG.eps',
  plot = barplot_unique_Nterm_Jurkat_GO_KEGG,
  height = 2, width = 3.5, units = 'in'
)

### figure S2B, protein domain and complex analysis of identified unique N-terminal proteoform in HEK293T, Jurkat and THP-1 cells
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
  filename = 'figures/figureS2/barplot_unique_Nterm_comb_pfam_corum.eps',
  plot = barplot_unique_Nterm_comb_pfam_corum,
  height = 2, width = 3.5, units = 'in'
)

### figure S2C, N-terminal proteoform topFinder information
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

## calculate the # of N-terminal proteoform in each cell line
# HEK293T
annotate_unique_Nterm_HEK293T <- unique_Nterm_topfinder_HEK293T |> 
  select(
    UniProt.curated.start, 
    Alternative.Spliced.Start,
    Cleaving.proteases,
    Other.experimental.terminus.evidences,
    Alternative.Translation.Start
  ) |> 
  pivot_longer(UniProt.curated.start:Alternative.Translation.Start, names_to = 'category', values_to = 'value') |> 
  filter(!is.na(value))

not_annotate_unique_Nterm_HEK293T <- unique_Nterm_topfinder_HEK293T |> 
  filter(if_all(c(
    UniProt.curated.start, 
    Alternative.Spliced.Start,
    Cleaving.proteases,
    Other.experimental.terminus.evidences,
    Alternative.Translation.Start
  ), is.na))

# unique Nterm HEK293T topfinder annotation
annotation_unique_Nterm_HEK293T <- bind_rows(
  annotate_unique_Nterm_HEK293T |> 
    count(category),
  
  tibble(
    category = 'Not.Annotated',
    n = not_annotate_unique_Nterm_HEK293T |> nrow()
  )
)

colnames(annotation_unique_Nterm_HEK293T) <- c('category', 'HEK293T')

# Jurkat
annotate_unique_Nterm_Jurkat <- unique_Nterm_topfinder_Jurkat |> 
  select(
    UniProt.curated.start, 
    Alternative.Spliced.Start,
    Cleaving.proteases,
    Other.experimental.terminus.evidences,
    Alternative.Translation.Start
  ) |> 
  pivot_longer(UniProt.curated.start:Alternative.Translation.Start, names_to = 'category', values_to = 'value') |> 
  filter(!is.na(value))

not_annotate_unique_Nterm_Jurkat <- unique_Nterm_topfinder_Jurkat |> 
  filter(if_all(c(
    UniProt.curated.start, 
    Alternative.Spliced.Start,
    Cleaving.proteases,
    Other.experimental.terminus.evidences,
    Alternative.Translation.Start
  ), is.na))

# unique Nterm Jurkat topfinder annotation
annotation_unique_Nterm_Jurkat <- bind_rows(
  annotate_unique_Nterm_Jurkat |> 
    count(category),
  
  tibble(
    category = 'Not.Annotated',
    n = not_annotate_unique_Nterm_Jurkat |> nrow()
  )
)

colnames(annotation_unique_Nterm_Jurkat) <- c('category', 'Jurkat')

# THP-1
annotate_unique_Nterm_THP1 <- unique_Nterm_topfinder_THP1 |> 
  select(
    UniProt.curated.start, 
    Alternative.Spliced.Start,
    Cleaving.proteases,
    Other.experimental.terminus.evidences,
    Alternative.Translation.Start
  ) |> 
  pivot_longer(UniProt.curated.start:Alternative.Translation.Start, names_to = 'category', values_to = 'value') |> 
  filter(!is.na(value))

not_annotate_unique_Nterm_THP1 <- unique_Nterm_topfinder_THP1 |> 
  filter(if_all(c(
    UniProt.curated.start, 
    Alternative.Spliced.Start,
    Cleaving.proteases,
    Other.experimental.terminus.evidences,
    Alternative.Translation.Start
  ), is.na))

# unique Nterm Jurkat topfinder annotation
annotation_unique_Nterm_THP1 <- bind_rows(
  annotate_unique_Nterm_THP1 |> 
    count(category),
  
  tibble(
    category = 'Not.Annotated',
    n = not_annotate_unique_Nterm_THP1 |> nrow()
  )
)

colnames(annotation_unique_Nterm_THP1) <- c('category', 'THP1')

# combine results from HEK293T, Jurkat and THP-1
annotation_unique_Nterm_comb <- annotation_unique_Nterm_HEK293T |> 
  left_join(
    annotation_unique_Nterm_Jurkat, by = 'category'
  ) |> 
  left_join(
    annotation_unique_Nterm_THP1, by = 'category'
  )

# chi-square test
annotation_unique_Nterm_comb_matrix <- data.matrix(
  annotation_unique_Nterm_comb |> select(-category)
)

rownames(annotation_unique_Nterm_comb_matrix) <- annotation_unique_Nterm_comb$category

# pairwise chi-square tests
# HEK293T vs. Jurkat
chisq.test(
  annotation_unique_Nterm_comb_matrix[, c('HEK293T', 'Jurkat')]
)

# Jurkat vs. THP-1
chisq.test(
  annotation_unique_Nterm_comb_matrix[, c('Jurkat', 'THP1')]
)

# THP-1 vs. HEK293T
chisq.test(
  annotation_unique_Nterm_comb_matrix[, c('THP1', 'HEK293T')]
)

# bar plot
barplot_annotation_unique_Nterm_comb <- annotation_unique_Nterm_comb |> 
  pivot_longer(
    cols = 'HEK293T':'THP1', names_to = 'cell', values_to = 'count'
  ) |> 
  mutate(
    category = factor(category, levels = c(
      'Not.Annotated',
      'Alternative.Translation.Start',
      'Other.experimental.terminus.evidences',
      'Cleaving.proteases',
      'Alternative.Spliced.Start',
      'UniProt.curated.start'
    )) 
  ) |> 
  ggplot() +
  geom_bar(
    aes(
      x = cell,
      y = count,
      fill = category
    ),
    stat = 'identity', 
    position = 'stack'
  ) +
  labs(x = '', y = '') +
  scale_fill_manual(
    name = '',
    values = c(
      'Not.Annotated' = 'gray70',
      'Alternative.Translation.Start' = color_5,
      'Other.experimental.terminus.evidences' = color_4,
      'Cleaving.proteases' = color_3,
      'Alternative.Spliced.Start' = color_2,
      'UniProt.curated.start' = color_1
    )
  ) +
  theme(
    axis.text.x = element_text(size = 8, color = 'black', angle = 30, hjust = 1),
    axis.text.y = element_text(size = 8, color = 'black'),
    legend.text = element_text(size = 8, color = 'black')
  )

ggsave(
  filename = 'figures/figureS2/barplot_annotation_unique_Nterm_comb.eps',
  plot = barplot_annotation_unique_Nterm_comb,
  height = 2, width = 4, units = 'in'
)

### figure S2D, Nterm structure features for solvent accessibility
## calculate the percentage of solvent accessibility for each cell line
# HEK293T
unique_Nterm_HEK293T_solvent_accessibility_percentage <- unique_Nterm_HEK293T_alphafold_annotation_R |> 
  filter(!is.na(high_acc_5)) |> 
  count(high_acc_5) |> 
  mutate(
    percentage = n/sum(n),
    cell = 'HEK293T'
  )

# Jurkat
unique_Nterm_Jurkat_solvent_accessibility_percentage <- unique_Nterm_Jurkat_alphafold_annotation_R |> 
  filter(!is.na(high_acc_5)) |> 
  count(high_acc_5) |> 
  mutate(
    percentage = n/sum(n),
    cell = 'Jurkat'
  )

# THP-1
unique_Nterm_THP1_solvent_accessibility_percentage <- unique_Nterm_THP1_alphafold_annotation_R |> 
  filter(!is.na(high_acc_5)) |> 
  count(high_acc_5) |> 
  mutate(
    percentage = n/sum(n),
    cell = 'THP-1'
  )

# combine result from HEK293T, Jurkat and THP-1
unique_Nterm_solvent_accessibility_percentage_comb <- bind_rows(
  unique_Nterm_HEK293T_solvent_accessibility_percentage,
  unique_Nterm_Jurkat_solvent_accessibility_percentage,
  unique_Nterm_THP1_solvent_accessibility_percentage
)

# ANOVA test/Turkey HSD
ANOVA_result <- aov(
  percentage ~ cell, data = unique_Nterm_solvent_accessibility_percentage_comb
)
summary(ANOVA_result)
TukeyHSD(ANOVA_result)

# bar plot
barplot_unique_Nterm_solvent_accessibility_percentage_comb <- unique_Nterm_solvent_accessibility_percentage_comb |> 
  mutate(
    high_acc_5 = as.character(high_acc_5)
  ) |> 
  ggplot() +
  geom_bar(
    aes(
      x = cell,
      y = percentage,
      fill = high_acc_5
    ), 
    stat = 'identity',
    position = 'stack'
  ) +
  labs(x = '', y = '') +
  scale_fill_manual(
    name = '',
    values = c(
      '0' = color_3,
      '1' = color_4
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
  filename = 'figures/figureS2/barplot_unique_Nterm_solvent_accessibility_percentage_comb.eps',
  plot = barplot_unique_Nterm_solvent_accessibility_percentage_comb,
  height = 2, width = 5, units = 'in'
)

### figure S2E, Nterm structure features for IDR
## calculate the percentage of IDR for each cell line
# HEK293T
unique_Nterm_HEK293T_IDR_percentage <- unique_Nterm_HEK293T_alphafold_annotation_R |> 
  filter(!is.na(IDR)) |> 
  count(IDR) |> 
  mutate(
    percentage = n/sum(n),
    cell = 'HEK293T'
  )

# Jurkat
unique_Nterm_Jurkat_IDR_percentage <- unique_Nterm_Jurkat_alphafold_annotation_R |> 
  filter(!is.na(IDR)) |> 
  count(IDR) |> 
  mutate(
    percentage = n/sum(n),
    cell = 'Jurkat'
  )

# THP-1
unique_Nterm_THP1_IDR_percentage <- unique_Nterm_THP1_alphafold_annotation_R |> 
  filter(!is.na(IDR)) |> 
  count(IDR) |> 
  mutate(
    percentage = n/sum(n),
    cell = 'THP-1'
  )

# combine result from HEK293T, Jurkat and THP-1
unique_Nterm_IDR_percentage_comb <- bind_rows(
  unique_Nterm_HEK293T_IDR_percentage,
  unique_Nterm_Jurkat_IDR_percentage,
  unique_Nterm_THP1_IDR_percentage
)

# ANOVA test/Turkey HSD
ANOVA_result <- aov(
  percentage ~ cell, data = unique_Nterm_IDR_percentage_comb
)
summary(ANOVA_result)
TukeyHSD(ANOVA_result)

# bar plot
barplot_unique_Nterm_IDR_percentage_comb <- unique_Nterm_IDR_percentage_comb |> 
  mutate(
    IDR = as.character(IDR)
  ) |> 
  ggplot() +
  geom_bar(
    aes(
      x = cell,
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
      '0' = color_2,
      '1' = color_1
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
  filename = 'figures/figureS2/barplot_unique_Nterm_IDR_percentage_comb.eps',
  plot = barplot_unique_Nterm_IDR_percentage_comb,
  height = 2, width = 5, units = 'in'
)

### figure S2F, THP-1 propeptide example
library(tidyverse)

# O14773, TPP1, Tripeptidyl-peptidase 1
O14773_database_info <- tribble(
  ~ name, ~ start, ~ end,
  'protein', 1, 563,
  'Peptidase_S8', 317, 490,
  'Pro-kuma_activ', 33, 176
)

O14773_result <- tibble(
  cleavage_site = c(
    20, 199, 324, 325, 511
  )
)

# example plot
O14773_example <- ggplot() +
  geom_rect(
    data = O14773_database_info,
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
      x = O14773_result$cleavage_site, 
      xend = O14773_result$cleavage_site, 
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
      'Peptidase_S8' = color_1,
      'Pro-kuma_activ' = color_2
    )
  ) +
  scale_color_manual(
    values = c(
      'protein' = 'black',
      'Peptidase_S8' = 'transparent',
      'Pro-kuma_activ' = 'transparent'
    )
  ) +
  theme_void()

ggsave(
  filename = 'figures/figureS2/O14773_example.eps',
  height = 0.2, width = 2, units = 'in'
)

# P25774, CTSS, Cathepsin S
P25774_database_info <- tribble(
  ~ name, ~ start, ~ end,
  'protein', 1, 331,
  'Peptidase_C1', 115, 329,
  'Inhibitor_I29', 28, 87
)

P25774_result <- tibble(
  cleavage_site = c(
    114
  )
)

# example plot
P25774_example <- ggplot() +
  geom_rect(
    data = P25774_database_info,
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
      x = P25774_result$cleavage_site, 
      xend = P25774_result$cleavage_site, 
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
      'Peptidase_C1' = color_3,
      'Inhibitor_I29' = color_4
    )
  ) +
  scale_color_manual(
    values = c(
      'protein' = 'black',
      'Peptidase_C1' = 'transparent',
      'Inhibitor_I29' = 'transparent'
    )
  ) +
  theme_void()

ggsave(
  filename = 'figures/figureS2/P25774_example.eps',
  height = 0.2, width = 2, units = 'in'
)

# Q9UBR2, CTSZ, Cathepsin Z
Q9UBR2_database_info <- tribble(
  ~ name, ~ start, ~ end,
  'protein', 1, 303,
  'Peptidase_C1', 62, 280
)

Q9UBR2_result <- tibble(
  cleavage_site = c(
    61, 70, 79, 192, 246, 260, 261
  )
)

# example plot
Q9UBR2_example <- ggplot() +
  geom_rect(
    data = Q9UBR2_database_info,
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
      x = Q9UBR2_result$cleavage_site, 
      xend = Q9UBR2_result$cleavage_site, 
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
      'Peptidase_C1' = color_3
    )
  ) +
  scale_color_manual(
    values = c(
      'protein' = 'black',
      'Peptidase_C1' = 'transparent'
    )
  ) +
  theme_void()

ggsave(
  filename = 'figures/figureS2/Q9UBR2_example.eps',
  height = 0.2, width = 2, units = 'in'
)

