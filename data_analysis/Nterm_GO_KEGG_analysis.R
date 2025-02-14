#import packages
packages_names <- c('tidyverse', 'clusterProfiler', 'org.Hs.eg.db')
lapply(packages_names, require, character.only = TRUE)

#Nterm Top 25%, short half life
Nterm_short_half_life_protein <- HEK_Nterm_Kd_half_life |> 
  filter(Percentile <= 0.25) |> 
  distinct(UniProt_Accession) |> 
  pull()

#Nterm Bottom 25%, long half life
Nterm_long_half_life_protein <- HEK_Nterm_Kd_half_life |> 
  filter(Percentile >= 0.75) |> 
  distinct(UniProt_Accession) |> 
  pull()

#Nterm total protein
Nterm_protein_total <- HEK_Nterm_Kd_half_life |> 
  distinct(UniProt_Accession) |> 
  pull()

##Gene Ontology analysis
#Top 25%, short half life
Nterm_short_half_life_GO <- enrichGO(
  gene = Nterm_short_half_life_protein,
  OrgDb = org.Hs.eg.db,
  universe = Nterm_protein_total,
  keyType = 'UNIPROT',
  ont = 'ALL',
  pvalueCutoff = 1,
  qvalueCutoff = 1
)

write_csv(
  Nterm_short_half_life_GO@result, file = 'data_source/GO_KEGG_analysis/Nterm_short_half_life_GO.csv'
)

#Bottom 25%, long half life
Nterm_long_half_life_GO <- enrichGO(
  gene = Nterm_long_half_life_protein,
  OrgDb = org.Hs.eg.db,
  universe = Nterm_protein_total,
  keyType = 'UNIPROT',
  ont = 'ALL',
  pvalueCutoff = 1,
  qvalueCutoff = 1
)

write_csv(
  Nterm_long_half_life_GO@result, file = 'data_source/GO_KEGG_analysis/Nterm_long_half_life_GO.csv'
)

##KEGG analysis
#Top 25%, short half life
Nterm_short_half_life_KEGG <- enrichKEGG(
  gene = Nterm_short_half_life_protein,
  organism = 'hsa',
  keyType = 'uniprot',
  universe = Nterm_protein_total,
  pvalueCutoff = 1, 
  qvalueCutoff = 1
)

write_csv(
  Nterm_short_half_life_KEGG@result, file = 'data_source/GO_KEGG_analysis/Nterm_short_half_life_KEGG.csv'
)

#Bottom 25%, long half life
Nterm_long_half_life_KEGG <- enrichKEGG(
  gene = Nterm_long_half_life_protein,
  organism = 'hsa',
  keyType = 'uniprot',
  universe = Nterm_protein_total,
  pvalueCutoff = 1, 
  qvalueCutoff = 1
)

write_csv(
  Nterm_long_half_life_KEGG@result, file = 'data_source/GO_KEGG_analysis/Nterm_long_half_life_KEGG.csv'
)
