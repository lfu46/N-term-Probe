# import packages
packages_names <- c("tidyverse", 'rstatix', 'clusterProfiler', 'org.Hs.eg.db')
lapply(packages_names, require, character.only = TRUE)

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
        pull(FunctionalSiteName)
    )
  ) %>%
  ungroup() %>%
  mutate(matched_motifs = sapply(matched_motifs, function(x) paste(x, collapse = ";"))) |> 
  separate_rows(matched_motifs, sep = ';') |> 
  distinct()

# calculate medain half life for each ELM motif
HEK_Nterm_sequence_ELM_motif_median <- HEK_Nterm_sequence_ELM_motif |> 
  group_by(matched_motifs) |> 
  get_summary_stats(half_life, type = 'median') |> 
  arrange(median)

write_csv(HEK_Nterm_sequence_ELM_motif_median, file = 'data_source/ELM_degron/HEK_Nterm_sequence_ELM_motif_median.csv')
