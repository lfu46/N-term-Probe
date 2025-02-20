# import packages
packages_names <- c('tidyverse', 'clusterProfiler', 'org.Hs.eg.db')
lapply(packages_names, require, character.only = TRUE)

# Nterm, Fast turnover, half-life < 7 h
Nterm_fast_turnover_protein <- HEK_Nterm_Kd_half_life_LaminB_Tcomplex |> 
  filter(category == 'Fast turnover') |> 
  distinct(UniProt_Accession) |> 
  pull()

# Nterm, Stable, half-life = 200 h
Nterm_stable_protein <- HEK_Nterm_Kd_half_life_LaminB_Tcomplex |> 
  filter(category == 'Stable') |> 
  distinct(UniProt_Accession) |> 
  pull()

# Nterm total protein
Nterm_protein_total <- HEK_Nterm_Kd_half_life_LaminB_Tcomplex |> 
  distinct(UniProt_Accession) |> 
  pull()

## Gene Ontology analysis
# Fast turnover, half-life < 7 h
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

# Stable, half-life = 200 h
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
# Fast turnover, half-life < 7 h
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

# Stable, half-life = 200 h
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
