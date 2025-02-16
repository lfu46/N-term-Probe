#import packages
packages_names <- c("tidyverse", "rstatix")
lapply(packages_names, require, character.only = TRUE)

#generate data frame for TopFinder (https://topfind.clip.msl.ubc.ca/topfinder) upload
topfinder_id <- HEK_Nterm_Kd_half_life_sequence |> 
  mutate(
    topfinder_id = paste(UniProt_Accession, Nterm_sequence, sep = ' ')
  ) |> 
  select(topfinder_id)

write_delim(topfinder_id, 
            file = 'data_source/Nterm_topfinder/topfinder_id.txt',
            col_names = FALSE,
            quote = 'none')

#import result from TopFinder
Nterm_topfinder_result <- read_delim(
  'data_source/Nterm_topfinder/2025_02_15_Nterm_02152025/2025_02_15_Nterm_02152025_Full_Table.txt',
  col_names = TRUE,
  name_repair = 'universal'
) |> 
  filter(Protein.found == 'YES') |> 
  mutate(Index = paste(Accession, P1..Position, sep = '_'))

#Nterm cleaving proteases
Nterm_topfinder_cleaving_proteases <- Nterm_topfinder_result |> 
  select(Index, Cleaving.proteases) |> 
  left_join(HEK_Nterm_Kd_half_life, by = 'Index')

#Wilcoxon rank-sum test
Nterm_topfinder_cleaving_proteases |> 
  wilcox_test(half_life ~ Cleaving.proteases) |> 
  filter(p < 0.05)

#Nterm signal peptide
Nterm_topfinder_signal_peptide <- Nterm_topfinder_result |> 
  select(Index, signal = Distance.To.signal.peptide) |> 
  filter(!is.na(signal)) |> 
  left_join(HEK_Nterm_Kd_half_life, by = 'Index') |> 
  mutate(
    signal_type = ifelse(signal >= -10 & signal <= 10, 'signal peptide', 'non signal peptide')
  )

#Wilcoxon rank-sum test
Nterm_topfinder_signal_peptide |> 
  wilcox_test(half_life ~ signal_type)

#Nterm propeptide
Nterm_topfinder_propeptide <- Nterm_topfinder_result |> 
  select(Index, propeptide = Distance.to.propeptide.lost) |> 
  filter(!is.na(propeptide)) |> 
  left_join(HEK_Nterm_Kd_half_life, by = 'Index') |> 
  mutate(
    propeptide_type = ifelse(propeptide >= -10 & propeptide <= 10, 'propeptide', 'non propeptide')
  )

#Wilcoxon rank-sum test
Nterm_topfinder_propeptide |> 
  wilcox_test(half_life ~ propeptide_type)

#Nterm transmembrane
Nterm_topfinder_transmembrane <- Nterm_topfinder_result |> 
  select(Index, transmembrane = Distance.to.last.transmembrane.domain..shed.) |> 
  filter(!is.na(transmembrane)) |> 
  left_join(HEK_Nterm_Kd_half_life, by = 'Index') |> 
  mutate(
    transmembrane_type = ifelse(transmembrane >= -10 & transmembrane <= 10, 'transmembrane', 'non transmembrane')
  )

#Wilcoxon rank-sum test
Nterm_topfinder_transmembrane |> 
  wilcox_test(half_life ~ transmembrane_type)

