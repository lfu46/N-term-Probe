# import packages
library(tidyverse)

# comparison of half-life between N-terminal proteoforms with their parent proteins
HEK_Nterm_WP_delta_half_life <- HEK_Nterm_Kd_half_life_LaminB_Tcomplex |> 
  left_join(
    HEK_WP_Kd_half_life_LaminB_Tcomplex, by = c('UniProt_Accession', 'Gene', 'Entry.Name'),
    suffix = c('.Nt', '.WP')
  ) |> 
  filter(!is.na(half_life.WP)) |> 
  select(Index:Entry.Name, half_life.Nt, RSS.Nt, half_life.WP, RSS.WP) |> 
  mutate(
    delta_half_life = half_life.Nt - half_life.WP,
    Percentile = percent_rank(delta_half_life)
  ) |> 
  arrange(Percentile)

write_csv(
  HEK_Nterm_WP_delta_half_life,
  file = 'data_source/Nterm_WP_comparison/HEK_Nterm_WP_delta_half_life.csv'
)

# total overlap protein
overlap_protein <- HEK_Nterm_WP_delta_half_life |> 
  distinct(UniProt_Accession) |> 
  pull()

# Top 20% protein, Nterm shorter half-life
top20_protein <- HEK_Nterm_WP_delta_half_life |> 
  filter(Percentile < 0.2) |> 
  distinct(UniProt_Accession) |> 
  pull()

# Bottom 20% protein, Nterm longer half-life
bottom20_protein <- HEK_Nterm_WP_delta_half_life |> 
  filter(Percentile > 0.8) |> 
  distinct(UniProt_Accession) |> 
  pull()

## GO and KEGG analysis
library(clusterProfiler)
library(org.Hs.eg.db)

# GO
top20_comparison_GO <- enrichGO(
  gene = top20_protein,
  OrgDb = org.Hs.eg.db,
  universe = overlap_protein,
  keyType = 'UNIPROT',
  ont = 'ALL',
  pvalueCutoff = 1,
  qvalueCutoff = 1
)

write_csv(
  top20_comparison_GO@result, file = 'data_source/Nterm_WP_comparison/top20_comparison_GO.csv'
)

bottom20_comparison_GO <- enrichGO(
  gene = bottom20_protein,
  OrgDb = org.Hs.eg.db,
  universe = overlap_protein,
  keyType = 'UNIPROT',
  ont = 'ALL',
  pvalueCutoff = 1,
  qvalueCutoff = 1
)

write_csv(
  bottom20_comparison_GO@result, file = 'data_source/Nterm_WP_comparison/bottom20_comparison_GO.csv'
)

# KEGG
top20_comparison_KEGG <- enrichKEGG(
  gene = top20_protein,
  organism = 'hsa',
  keyType = 'uniprot',
  universe = overlap_protein,
  pvalueCutoff = 1, 
  qvalueCutoff = 1
)

write_csv(
  top20_comparison_KEGG@result, file = 'data_source/Nterm_WP_comparison/top20_comparison_KEGG.csv'
)

bottom20_comparison_KEGG <- enrichKEGG(
  gene = bottom20_protein,
  organism = 'hsa',
  keyType = 'uniprot',
  universe = overlap_protein,
  pvalueCutoff = 1, 
  qvalueCutoff = 1
)

write_csv(
  bottom20_comparison_KEGG@result, file = 'data_source/Nterm_WP_comparison/bottom20_comparison_KEGG.csv'
)

## protein complex and domain analysis
library(clusterProfiler)

# Pfam
human_pfam <- read_tsv(
  'data_source/pfam/9606.tsv',
  skip = 3, 
  col_names = FALSE
) |> 
  dplyr::select(X6, X1)

colnames(human_pfam) <- c('TERM', 'GENE')

top20_comparison_pfam <- enricher(
  gene = top20_protein,
  pvalueCutoff = 1,
  universe = overlap_protein,
  qvalueCutoff = 1,
  TERM2GENE = human_pfam
)

write_csv(
  top20_comparison_pfam@result,
  'data_source/Nterm_WP_comparison/top20_comparison_pfam.csv'
)

bottom20_comparison_pfam <- enricher(
  gene = bottom20_protein,
  pvalueCutoff = 1,
  universe = overlap_protein,
  qvalueCutoff = 1,
  TERM2GENE = human_pfam
)

write_csv(
  bottom20_comparison_pfam@result,
  'data_source/Nterm_WP_comparison/bottom20_comparison_pfam.csv'
)

# CORUM
human_corum <- read_delim(
  'data_source/corum/corum_uniprotCorumMapping.txt',
  col_names = TRUE,
  name_repair = 'universal'
) |> 
  mutate(
    corum_id = paste('corum_id', corum_id, sep = '_')
  ) |> 
  dplyr::select(corum_id, UniProtKB_accession_number)

colnames(human_corum) <- c('TERM', 'GENE')

top20_comparison_corum <- enricher(
  gene = top20_protein,
  pvalueCutoff = 1,
  universe = overlap_protein,
  qvalueCutoff = 1,
  TERM2GENE = human_corum
)

write_csv(
  top20_comparison_corum@result,
  'data_source/Nterm_WP_comparison/top20_comparison_corum.csv'
)

bottom20_comparison_corum <- enricher(
  gene = bottom20_protein,
  pvalueCutoff = 1,
  universe = overlap_protein,
  qvalueCutoff = 1,
  TERM2GENE = human_corum
)

write_csv(
  bottom20_comparison_corum@result,
  'data_source/Nterm_WP_comparison/bottom20_comparison_corum.csv'
)

## protein structure analysis
library(reticulate)

# use specific virtual env created by anaconda
use_condaenv(
  condaenv = '/opt/anaconda3/envs/structuremap',
  required = TRUE
)

# execute the python script for Nterm structure
source_python("data_analysis/Nterm_WP_comparison_structuremap.py")

library(tidyverse)
# output the result for N-terminus structure
top20_bottom20_protein_structure <- Nterm_WP_delta_half_life_alphafold_N_terminus |> 
  as_tibble() |> 
  filter(top20 == 1 | bottom20 == 1)

write_csv(
  top20_bottom20_protein_structure,
  file = 'data_source/Nterm_WP_comparison/top20_bottom20_protein_structure.csv'
)

# output the result for N-terminus enrichment analysis
write_csv(
  enrichment_structure,
  file = 'data_source/Nterm_WP_comparison/enrichment_structure.csv'
)

# output the result for N-terminus 3D clustering
write_csv(
  enrichment_top20_bottom20_proximity,
  file = 'data_source/Nterm_WP_comparison/enrichment_top20_bottom20_proximity.csv'
)
