# import packages
packages_names <- c("tidyverse", "rstatix", "clusterProfiler", "org.Hs.eg.db")
lapply(packages_names, require, character.only = TRUE)

# get annotation from ENZYME, PFAM and PROSITE
enzyme_pfam_prosite_annotation <- tibble(AnnotationDbi::select(
  org.Hs.eg.db, keys = HEK_Nterm_Kd_half_life_LaminB_Tcomplex$UniProt_Accession, 
  columns = c("ENZYME", "PFAM", "PROSITE"), 
  keytype = "UNIPROT")
) |> 
  separate(ENZYME, sep = '\\.', into = c('ENZYME', 'sub1', 'sub2', 'sub3')) |> 
  dplyr::select(!sub1:sub3) |> 
  mutate(
    ENZYME = ifelse(!is.na(ENZYME), paste('ENZYME', ENZYME, sep = '_'), ENZYME)
  ) |> 
  pivot_longer(cols = ENZYME:PROSITE, names_to = 'Type', values_to = 'TERM') |> 
  dplyr::select(TERM, GENE = UNIPROT) |> 
  filter(!is.na(TERM))

# get annotation from CORUM (https://mips.helmholtz-muenchen.de/corum/)
complex_annotation <- read_delim(
  'data_source/Enzyme_Motif_Domain_Complex_analysis/corum_uniprotCorumMapping.txt',
  col_names = TRUE
) |> 
  mutate(
    CORUM_id = paste('CORUM_id', corum_id, sep = '_')
  ) |> 
  dplyr::select(UniProt_Accession = UniProtKB_accession_number, CORUM_id) |> 
  filter(UniProt_Accession %in% HEK_Nterm_Kd_half_life_LaminB_Tcomplex$UniProt_Accession)
colnames(complex_annotation) <- c('TERM', 'GENE')

# combine ENZYME, PFAM, PROSITE and CORUM
functional_property_database <- bind_rows(
  enzyme_pfam_prosite_annotation,
  complex_annotation
)

# Nterm protein half life and Kd
HEK_Nterm_protein_Kd_half_life <- HEK_Nterm_Kd_half_life_LaminB_Tcomplex |> 
  group_by(UniProt_Accession) |> 
  summarize(
    half_life_median = median(half_life),
    Kd_median = median(Kd)
  )

## Nterm protein functional property GSEA
# descending Kd
geneList_des_Kd <- HEK_Nterm_protein_Kd_half_life$Kd_median
names(geneList_des_Kd) <- HEK_Nterm_protein_Kd_half_life$UniProt_Accession
geneList_des_Kd <- sort(geneList_des_Kd, decreasing = TRUE)

Nterm_protein_functional_property_GSEA_des_Kd <- GSEA(
  geneList = geneList_des_Kd,
  TERM2GENE = functional_property_database,
  pvalueCutoff = 1,
  scoreType = 'pos'
)

write_csv(Nterm_protein_functional_property_GSEA_des_Kd@result, file = 'data_source/Enzyme_Motif_Domain_Complex_analysis/Nterm_protein_functional_property_GSEA_des_Kd.csv')

# descending half life
geneList_des_half_life <- HEK_Nterm_protein_Kd_half_life$half_life_median
names(geneList_des_half_life) <- HEK_Nterm_protein_Kd_half_life$UniProt_Accession
geneList_des_half_life <- sort(geneList_des_half_life, decreasing = TRUE)

Nterm_protein_functional_property_GSEA_des_half_life <- GSEA(
  geneList = geneList_des_half_life,
  TERM2GENE = functional_property_database,
  pvalueCutoff = 1,
  scoreType = 'pos'
)

write_csv(Nterm_protein_functional_property_GSEA_des_half_life@result, file = 'data_source/Enzyme_Motif_Domain_Complex_analysis/Nterm_protein_functional_property_GSEA_des_half_life.csv')
