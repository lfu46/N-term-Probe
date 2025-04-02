# import packages
library(tidyverse)
library(clusterProfiler)
library(org.Hs.eg.db)

# Nterm, Fast turnover, Top 20%
Nterm_fast_turnover_protein <- HEK_Nterm_Kd_half_life_LaminB_Tcomplex |> 
  filter(Percentile < 0.2) |> 
  distinct(UniProt_Accession) |> 
  pull()

# Nterm, Stable, Bottom 20%
Nterm_stable_protein <- HEK_Nterm_Kd_half_life_LaminB_Tcomplex |> 
  filter(Percentile > 0.8) |> 
  distinct(UniProt_Accession) |> 
  pull()

# Nterm total protein
Nterm_protein_total <- HEK_Nterm_Kd_half_life_LaminB_Tcomplex |> 
  distinct(UniProt_Accession) |> 
  pull()

## Gene Ontology analysis
# Fast turnover, Top 20%
Nterm_fast_turnover_GO <- enrichGO(
  gene = Nterm_fast_turnover_protein,
  OrgDb = org.Hs.eg.db,
  universe = Nterm_protein_total,
  keyType = 'UNIPROT',
  ont = 'ALL',
  pvalueCutoff = 1,
  qvalueCutoff = 1
)

write_csv(
  Nterm_fast_turnover_GO@result, file = 'data_source/GO_KEGG_analysis/Nterm_fast_turnover_GO.csv'
)

# Stable, Bottom 20%
Nterm_stable_GO <- enrichGO(
  gene = Nterm_stable_protein,
  OrgDb = org.Hs.eg.db,
  universe = Nterm_protein_total,
  keyType = 'UNIPROT',
  ont = 'ALL',
  pvalueCutoff = 1,
  qvalueCutoff = 1
)

write_csv(
  Nterm_stable_GO@result, file = 'data_source/GO_KEGG_analysis/Nterm_stable_GO.csv'
)

## KEGG analysis
# Fast turnover, Top 20%
Nterm_fast_turnover_KEGG <- enrichKEGG(
  gene = Nterm_fast_turnover_protein,
  organism = 'hsa',
  keyType = 'uniprot',
  universe = Nterm_protein_total,
  pvalueCutoff = 1, 
  qvalueCutoff = 1
)

write_csv(
  Nterm_fast_turnover_KEGG@result, file = 'data_source/GO_KEGG_analysis/Nterm_fast_turnover_KEGG.csv'
)

# Stable, Bottom 20%
Nterm_stable_KEGG <- enrichKEGG(
  gene = Nterm_stable_protein,
  organism = 'hsa',
  keyType = 'uniprot',
  universe = Nterm_protein_total,
  pvalueCutoff = 1, 
  qvalueCutoff = 1
)

write_csv(
  Nterm_stable_KEGG@result, file = 'data_source/GO_KEGG_analysis/Nterm_stable_KEGG.csv'
)
