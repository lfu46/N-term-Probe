# import packages
library(tidyverse)

# generate data frame for TopFinder (https://topfind.clip.msl.ubc.ca/topfinder) upload
topfinder_id <- HEK_Nterm_Kd_half_life_sequence |> 
  mutate(
    topfinder_id = paste(UniProt_Accession, Nterm_sequence, sep = ' ')
  ) |> 
  select(topfinder_id)

write_delim(topfinder_id, 
            file = 'data_source/Nterm_topfinder/topfinder_id.txt',
            col_names = FALSE,
            quote = 'none')

# import result from TopFinder
Nterm_topfinder_result <- read_delim(
  'data_source/Nterm_topfinder/2025_04_02_Nterm_degradation_04022025_R/2025_04_02_Nterm_degradation_04022025_R_Full_Table.txt',
  col_names = TRUE,
  name_repair = 'universal'
) |> 
  filter(Protein.found == 'YES') |> 
  mutate(Index = paste(Accession, P1..Position, sep = '_'))

# Nterm cleaving proteases
Nterm_degron_topfinder_cleaving_proteases <- Nterm_topfinder_result |> 
  select(Index, Cleaving.proteases) |> 
  left_join(HEK_Nterm_ELM_N_degron, by = 'Index') |> 
  filter(!is.na(Cleaving.proteases)) |> 
  separate_rows(Cleaving.proteases, sep = ';')

# cleaving protease, N-degron
cleaving_protease_N_degron <- Nterm_degron_topfinder_cleaving_proteases |> 
  group_by(Cleaving.proteases) |> 
  count(ELM_N_degron)

write_csv(
  cleaving_protease_N_degron,
  file = 'data_source/ELM_degron/cleaving_protease_N_degron.csv'
)

# N-degron, cleaving protease
N_degron_cleaving_protease <- Nterm_degron_topfinder_cleaving_proteases |> 
  group_by(ELM_N_degron) |> 
  count(Cleaving.proteases)

write_csv(
  N_degron_cleaving_protease,
  file = 'data_source/ELM_degron/N_degron_cleaving_protease.csv'
)

# Wilcoxon rank-sum test
library(rstatix)

Cleaving.Proteases.List <- Nterm_degron_topfinder_cleaving_proteases |> 
  group_by(Cleaving.proteases) |> 
  get_summary_stats(half_life, type = 'median') |> 
  filter(n > 1) |> 
  pull(Cleaving.proteases)

Nterm_degron_topfinder_cleaving_proteases |>
  filter(Cleaving.proteases %in% Cleaving.Proteases.List) |> 
  wilcox_test(half_life ~ Cleaving.proteases) |>
  filter(p < 0.05)

# Nterm signal peptide
Nterm_topfinder_signal_peptide <- Nterm_topfinder_result |> 
  select(Index, signal = Distance.To.signal.peptide) |> 
  filter(!is.na(signal)) |> 
  left_join(HEK_Nterm_Kd_half_life_LaminB_Tcomplex, by = 'Index') |> 
  mutate(
    signal_type = ifelse(signal >= -10 & signal <= 10, 'signal peptide', 'non signal peptide')
  )

# Wilcoxon rank-sum test
Nterm_topfinder_signal_peptide |> 
  wilcox_test(half_life ~ signal_type)

# Nterm propeptide
Nterm_topfinder_propeptide <- Nterm_topfinder_result |> 
  select(Index, propeptide = Distance.to.propeptide.lost) |> 
  filter(!is.na(propeptide)) |> 
  left_join(HEK_Nterm_Kd_half_life_LaminB_Tcomplex, by = 'Index') |> 
  mutate(
    propeptide_type = ifelse(propeptide >= -10 & propeptide <= 10, 'propeptide', 'non propeptide')
  )

# Wilcoxon rank-sum test
Nterm_topfinder_propeptide |> 
  wilcox_test(half_life ~ propeptide_type)

# Nterm transmembrane
Nterm_topfinder_transmembrane <- Nterm_topfinder_result |> 
  select(Index, transmembrane = Distance.to.last.transmembrane.domain..shed.) |> 
  filter(!is.na(transmembrane)) |> 
  left_join(HEK_Nterm_Kd_half_life_LaminB_Tcomplex, by = 'Index') |> 
  mutate(
    transmembrane_type = ifelse(transmembrane >= -10 & transmembrane <= 10, 'transmembrane', 'non transmembrane')
  )

# Wilcoxon rank-sum test
Nterm_topfinder_transmembrane |> 
  wilcox_test(half_life ~ transmembrane_type)

