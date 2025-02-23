library(tidyverse)

# import phosphorylation dataset from PhosphoSitePlus v6.7.7 
# (https://www.phosphosite.org/homeAction, 'Downloads/Phosphorylation_site_dataset.gz')
phosphorylation_site_human <- read_delim(
  'data_source/phosphorylation_site/Phosphorylation_site_dataset',
  skip = 3, 
  col_names = TRUE,
  name_repair = 'universal'
) |> 
  filter(
    ORGANISM == 'human',
    Ambiguous_Site == 0
  ) |> 
  select(UniProt_Accession = ACC_ID, MOD_RSD) |> 
  separate(MOD_RSD, into = c('Site', 'p'), sep = '-') |> 
  mutate(
    site_position = as.numeric(str_extract(Site, '\\d+'))
  ) |> 
  select(UniProt_Accession, site_position) |> 
  distinct()

# integrate the Nterm half life results with phosphorylation site dataset
Nterm_phosphorylation_site_occurrence <- phosphorylation_site_human |> 
  left_join(HEK_Nterm_Kd_half_life_sequence, by = 'UniProt_Accession', relationship = 'many-to-many') |> 
  mutate(
    start = Protein.Start,
    end = ifelse(Protein.Start + 12 <= Full_Protein_Length, Protein.Start + 12, Full_Protein_Length)
  ) |> 
  select(UniProt_Accession, site_position, start, end, half_life, Kd, Index, Nterm_sequence, Nterm_13mer, category) |> 
  mutate(
    occupancy = ifelse(site_position >= start & site_position <= end, 'occupied', 'unoccupied')
  ) |> 
  filter(occupancy == 'occupied') |> 
  group_by(Index, half_life, category) |> 
  count(occupancy) |> 
  ungroup()

# Wilcoxon rank-sum test
library(rstatix)

Nterm_phosphorylation_site_occurrence |> 
  group_by(category) |> 
  get_summary_stats(n)
  # wilcox_test(n ~ category)

# Spearman correlation test
cor.test(
  Nterm_phosphorylation_site_occurrence$half_life,
  Nterm_phosphorylation_site_occurrence$n,
  method = 'spearman'
)
