# import packages
library(tidyverse)

# import ubiquitination dataset from PhosphoSitePlus v6.7.7 
# (https://www.phosphosite.org/homeAction, 'Downloads/Ubiquitination_site_dataset.gz')
ubiquitination_site_human <- read_delim(
  'data_source/ubiquitination_site/Ubiquitination_site_dataset',
  skip = 3, 
  col_names = TRUE,
  name_repair = 'universal'
) |> 
  filter(
    ORGANISM == 'human',
    Ambiguous_Site == 0
  ) |> 
  select(UniProt_Accession = ACC_ID, MOD_RSD) |> 
  separate(MOD_RSD, into = c('Site', 'ub'), sep = '-') |> 
  mutate(
    site_position = as.numeric(str_extract(Site, '\\d+'))
  ) |> 
  select(UniProt_Accession, site_position) |> 
  distinct()

# integrate the Nterm half life results with ubiquitination site dataset
Nterm_ubiquitination_site_occurrence <- ubiquitination_site_human |> 
  left_join(HEK_Nterm_Kd_half_life_sequence, by = 'UniProt_Accession', relationship = 'many-to-many') |> 
  mutate(
    start = Protein.Start,
    end = ifelse(Protein.Start + 12 <= Full_Protein_Length, Protein.Start + 12, Full_Protein_Length)
  ) |> 
  select(UniProt_Accession, site_position, start, end, half_life, Kd_adj, Index, Nterm_sequence, Nterm_13mer, category) |> 
  mutate(
    occupancy = ifelse(site_position >= start & site_position <= end, 'occupied', 'unoccupied'),
    Lys_number = str_count(Nterm_13mer, 'K')
  ) |> 
  filter(occupancy == 'occupied') |> 
  group_by(Index, Lys_number ,half_life, category) |> 
  count(occupancy) |> 
  mutate(
    occupancy_percentage = n / Lys_number
  ) |> 
  ungroup() |> 
  filter(occupancy_percentage != Inf)

# Wilcoxon rank-sum test
library(rstatix)

Nterm_ubiquitination_site_occurrence |> 
  wilcox_test(occupancy_percentage ~ category)

# Spearman correlation test
cor.test(
  Nterm_ubiquitination_site_occurrence$half_life,
  Nterm_ubiquitination_site_occurrence$occupancy_percentage,
  method = 'spearman'
)
