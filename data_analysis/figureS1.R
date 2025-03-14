# import packages
library(tidyverse)

### figure S1A, venn diagram, HEK293T duplicates
# import result from MSFragger
HEK293T_1_07102024 <- read_tsv(
  'data_source/HEK293T_duplicates/NGlyco_HEK293T_1_2_1_07102024_psm.tsv',
  col_names = TRUE,
  name_repair = "universal"
)

HEK293T_2_07102024 <- read_tsv(
  'data_source/HEK293T_duplicates/NGlyco_HEK293T_1_2_2_07102024_psm.tsv',
  col_names = TRUE,
  name_repair = "universal"
)

## filter out the PSMs which contain N-term modification
# HEK293T_1_07102024
HEK293T_1_07102024_N_term <- HEK293T_1_07102024 |> 
  filter(
    str_detect(Assigned.Modifications, "N-term"),
    str_detect(Entry.Name, 'HUMAN')
  ) |> 
  select(
    UniProt_Accession = Protein.ID,
    Protein.Start,
    Gene,
    Entry.Name
  ) |> 
  mutate(
    Index = paste(UniProt_Accession, Protein.Start, sep = '_'),
    Exp = "Exp1"
  )

HEK293T_1_07102024_N_term_distinct <- HEK293T_1_07102024_N_term |> 
  distinct()

# HEK293T_2_07102024
HEK293T_2_07102024_N_term <- HEK293T_2_07102024 |> 
  filter(
    str_detect(Assigned.Modifications, "N-term"),
    str_detect(Entry.Name, 'HUMAN')
  ) |> 
  select(
    UniProt_Accession = Protein.ID,
    Protein.Start,
    Gene,
    Entry.Name
  ) |> 
  mutate(
    Index = paste(UniProt_Accession, Protein.Start, sep = '_'),
    Exp = "Exp2"
  )

HEK293T_2_07102024_N_term_distinct <- HEK293T_2_07102024_N_term |> 
  distinct()

## remove N-term which only identified in one replicate with one psm
overlap_HEK293T <- intersect(
  HEK293T_1_07102024_N_term_distinct |> pull(Index),
  HEK293T_2_07102024_N_term_distinct |> pull(Index)
)

# HEK293T_1_07102024
HEK293T_1_07102024_N_term_unique_filtered <- HEK293T_1_07102024_N_term |> 
  filter(! Index %in% overlap_HEK293T) |> 
  group_by(Index) |> 
  filter(n() > 1) |> 
  ungroup() |> 
  distinct()

# HEK293T_2_07102024
HEK293T_2_07102024_N_term_unique_filtered <- HEK293T_2_07102024_N_term |> 
  filter(! Index %in% overlap_HEK293T) |> 
  group_by(Index) |> 
  filter(n() > 1) |> 
  ungroup() |> 
  distinct()

## combine overlap and unique N-term in each replicate
# HEK293T_1_07102024
HEK293T_1_07102024_N_term_combined <- bind_rows(
  HEK293T_1_07102024_N_term_distinct |> 
    filter(Index %in% overlap_HEK293T),
  HEK293T_1_07102024_N_term_unique_filtered
)

# HEK293T_2_07102024
HEK293T_2_07102024_N_term_combined <- bind_rows(
  HEK293T_2_07102024_N_term_distinct |> 
    filter(Index %in% overlap_HEK293T),
  HEK293T_2_07102024_N_term_unique_filtered
)

# venn diagram
library(VennDiagram)

venn.diagram(
  x = list(HEK293T_1_07102024_N_term_combined$Index, HEK293T_2_07102024_N_term_combined$Index),
  category.names = c("", ""),
  cex = 0,
  filename = 'figures/figureS1/venn_diagram_HEK293T_identification.tiff', 
  fill = c(color_1, color_2),
  output = FALSE,
  imagetype = "tiff",
  height = 2,
  width = 2,
  units = c("in"),
  resolution = 1200,
  margin = 0.02,
  col = 'transparent'
)

### figure S1B, venn diagram, Jurkat and THP-1
# import result from MSFragger
# Jurkat
Jurkat_07102024 <- read_tsv(
  "data_source/Jurkat/N_term_Jurkat_07102024_psm.tsv",
  col_names = TRUE,
  name_repair = "universal"
)

Jurkat_07102024_N_term <- Jurkat_07102024 |> 
  filter(
    str_detect(Assigned.Modifications, "N-term"),
    str_detect(Entry.Name, 'HUMAN')
  ) |> 
  select(
    UniProt_Accession = Protein.ID,
    Protein.Start,
    Gene,
    Entry.Name
  ) |> 
  mutate(
    Index = paste(UniProt_Accession, Protein.Start, sep = '_')
  ) |> 
  distinct()

# THP-1
THP1_07102024 <- read_tsv(
  "data_source/THP1/N_term_THP1_07102024_psm.tsv",
  col_names = TRUE,
  name_repair = "universal"
)

THP1_07102024_N_term <- THP1_07102024 |> 
  filter(
    str_detect(Assigned.Modifications, "N-term"),
    str_detect(Entry.Name, 'HUMAN')
  ) |> 
  select(
    UniProt_Accession = Protein.ID,
    Protein.Start,
    Gene,
    Entry.Name
  ) |> 
  mutate(
    Index = paste(UniProt_Accession, Protein.Start, sep = '_')
  ) |> 
  distinct()

# check the overlap Nterm
overlap_Jurkat_THP1 <- intersect(
  Jurkat_07102024_N_term$Index, THP1_07102024_N_term$Index
)

# venn diagram
library(VennDiagram)

venn.diagram(
  x = list(Jurkat_07102024_N_term$Index, THP1_07102024_N_term$Index),
  category.names = c("", ""),
  cex = 0,
  filename = 'figures/figureS1/venn_diagram_Jurkat_THP1_identification.tiff', 
  fill = c(color_3, color_4),
  output = FALSE,
  imagetype = "tiff",
  height = 2,
  width = 2,
  units = c("in"),
  resolution = 1200,
  margin = 0.02,
  col = 'transparent'
)

### figure S1C
# CC
barplot_GO_CC_common_Nterm <- GO_common_Nterm |> 
  filter(Description %in% c(
    'cytosolic small ribosomal subunit',
    'protein folding chaperone complex',
    'Sm-like protein family complex',
    'spliceosomal snRNP complex',
    'proteasome complex',
    'Lsm1-7-Pat1 complex',
    'chaperonin-containing T-complex',
    'U2 snRNP',
    'U6 snRNP',
    'proton-transporting two-sector ATPase complex'
  )) |> 
  ggplot() +
  geom_bar(
    aes(
      x = fct_reorder(Description, -log10(p.adjust)), 
      y = -log10(p.adjust)
    ),
    fill = color_2, color = 'transparent', stat = 'identity'
  ) +
  labs(x = 'Cellular Component', y = '-log10(adjust P value)') +
  coord_flip() +
  theme(
    axis.text = element_text(size = 8, color = 'black')
  )

ggsave(
  filename = 'figures/figureS1/barplot_GO_CC_common_Nterm.eps',
  plot = barplot_GO_CC_common_Nterm,
  height = 2, width = 4, units = 'in'
)

### figure S1D
# MF
barplot_GO_MF_common_Nterm <- GO_common_Nterm |> 
  filter(Description %in% c(
    'intramolecular oxidoreductase activity',
    'protein-disulfide reductase activity',
    'ATPase regulator activity',
    'adenyl-nucleotide exchange factor activity',
    'translation initiation factor activity',
    'thioredoxin peroxidase activity',
    'aldehyde-lyase activity',
    'proton channel activity',
    '2 iron, 2 sulfur cluster binding',
    'rRNA binding'
  )) |> 
  ggplot() +
  geom_bar(
    aes(
      x = fct_reorder(Description, -log10(p.adjust)), 
      y = -log10(p.adjust)
    ),
    fill = color_3, color = 'transparent', stat = 'identity'
  ) +
  labs(x = 'Molecular Function', y = '-log10(adjust P value)') +
  coord_flip() +
  theme(
    axis.text = element_text(size = 8, color = 'black')
  )

ggsave(
  filename = 'figures/figureS1/barplot_GO_MF_common_Nterm.eps',
  plot = barplot_GO_MF_common_Nterm,
  height = 2, width = 4, units = 'in'
)

### figure S1E
# calculate the percentage of IDR
common_Nterm_IDR_percentage <- common_Nterm_alphafold_N_terminus_tb |> 
  count(IDR) |> 
  mutate(
    percentage = (n / sum(n)) * 100
  ) |> 
  mutate(
    IDR = as.character(IDR)
  )

# circular barplot
circular_barplot_IDR_percentage <- common_Nterm_IDR_percentage |> 
  ggplot() +
  geom_bar(
    aes(
      x = 3, 
      y = percentage, 
      fill = IDR
    ),
    stat = 'identity', 
    show.legend = FALSE
  ) +
  scale_fill_manual(
    values = c(
      '0' = color_2,
      '1' = color_1
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
  filename = 'figures/figureS1/circular_barplot_IDR_percentage.eps',
  plot = circular_barplot_IDR_percentage,
  height = 2, width = 2, units = 'in'
)

### figure S1F
## Q15459, SF3A1, Splicing factor 3A subunit 1
library(tidyverse)

Q15459_database_info <- tribble(
  ~ name, ~ start, ~ end,
  'protein', 1, 793,
  'ubiquitin', 722, 788,
  'Surp', 50, 101,
  'Surp', 164, 215,
  'PRP21_like_P', 233, 472
)

Q15459_result <- tibble(
  cleavage_site = c(2, 3, 6, 51, 116, 132, 382, 397, 403, 458, 473, 685)
)

# example plot
Q15459_example <- ggplot() +
  geom_rect(
    data = Q15459_database_info,
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
      x = Q15459_result$cleavage_site, 
      xend = Q15459_result$cleavage_site, 
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
      'ubiquitin' = color_4,
      'Surp' = color_5,
      'PRP21_like_P' = color_1
    )
  ) +
  scale_color_manual(
    values = c(
      'protein' = 'black',
      'ubiquitin' = 'transparent',
      'Surp' = 'transparent',
      'PRP21_like_P' = 'transparent'
    )
  ) +
  theme_void()

ggsave(
  filename = 'figures/figureS1/Q15459_example.eps',
  height = 0.2, width = 2, units = 'in'
)
