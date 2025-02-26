# import packages
library(tidyverse)

### figure 2A, venn diagram, HEK293T duplicates
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
  filename = 'figures/figure2/venn_diagram_HEK293T_identification.tiff', 
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

### figure 2B, venn diagram, Jurkat and THP-1
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
  filename = 'figures/figure2/venn_diagram_Jurkat_THP1_identification.tiff', 
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

### figure 2C, UpSet plot for HEK293T, Jurkat and THP-1
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

### figure 2D, GO analysis of commonly identified N-term in HEK293T, Jurkat and THP-1 cells
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
  filename = 'figures/figure2/barplot_GO_CC_common_Nterm.eps',
  plot = barplot_GO_CC_common_Nterm,
  height = 2, width = 4, units = 'in'
)

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
  filename = 'figures/figure2/barplot_GO_MF_common_Nterm.eps',
  plot = barplot_GO_MF_common_Nterm,
  height = 2, width = 4, units = 'in'
)

### figure 2E, KEGG analysis of commonly identified N-term in HEK293T, Jurkat and THP-1 cells
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

### figure 2F, N-termial proteoform topFinder information
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

