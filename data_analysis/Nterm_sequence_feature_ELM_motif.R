# import packages
library(tidyverse)

# import ELM classes downloaded from 
# The Eukaryotic Linear Motif resources for Functional Sites in Proteins (http://elm.eu.org/searchdb.html)
ELM_classes <- read_tsv(
  'data_source/ELM_degron/elm_classes.tsv',
  skip = 5,
  col_names = TRUE,
  name_repair = 'universal'
)

# annotate Nterm sequences using motif regex
HEK_Nterm_sequence_ELM_motif <- HEK_Nterm_Kd_half_life_sequence |> 
  rowwise() |> 
  mutate(
    matched_motifs = list(
      ELM_classes %>%
        filter(str_detect(Nterm_sequence, Regex)) %>%
        pull(ELMIdentifier)
    )
  ) %>%
  ungroup() %>%
  mutate(matched_motifs = sapply(matched_motifs, function(x) paste(x, collapse = ";"))) |> 
  separate_rows(matched_motifs, sep = ';') |> 
  distinct()

# calculate medain half life for each ELM motif
library(rstatix)

HEK_Nterm_sequence_ELM_motif_median <- HEK_Nterm_sequence_ELM_motif |> 
  group_by(matched_motifs) |> 
  get_summary_stats(half_life, type = 'median') |> 
  arrange(median)

write_csv(HEK_Nterm_sequence_ELM_motif_median, file = 'data_source/ELM_degron/HEK_Nterm_sequence_ELM_motif_median.csv')

## GSEA for ELM motif
# construct database
ELM_motif_database <- HEK_Nterm_sequence_ELM_motif |> 
  select(matched_motifs, Index) |> 
  distinct()

# genelist descending Kd
geneList_des_Kd <- HEK_Nterm_Kd_half_life_sequence$Kd
names(geneList_des_Kd) <- HEK_Nterm_Kd_half_life_sequence$Index
geneList_des_Kd <- sort(geneList_des_Kd, decreasing = TRUE)

# GSEA
library(clusterProfiler)

Nterm_ELM_motif_GSEA_des_Kd <- GSEA(
  geneList = geneList_des_Kd,
  TERM2GENE = ELM_motif_database,
  pvalueCutoff = 1,
  scoreType = 'pos'
)

write_csv(
  Nterm_ELM_motif_GSEA_des_Kd@result,
  'data_source/ELM_degron/Nterm_ELM_motif_GSEA_des_Kd.csv'
)
