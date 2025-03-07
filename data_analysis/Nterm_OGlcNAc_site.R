# import packages
library(tidyverse)
library(readxl)

# OGlcNAc_site_HEK293T <- read_xlsx(
#   'data_source/OGlcNAc_site/OG_glycopeptide_Top_tb_HEK293T_singlesite.xlsx',
#   col_name = TRUE,
#   .name_repair = 'universal'
# ) |> 
#   select(UniProt_Accession = UniprotID, site_position)

OGlcNAc_site_HEK293T <- read_delim(
  'data_source/OGlcNAc_site/O-GlcNAc_site_dataset',
  skip = 3, 
  col_names = TRUE,
  name_repair = 'universal'
) |> 
  filter(
    ORGANISM == 'human',
    Ambiguous_Site == 0
  ) |> 
  select(UniProt_Accession = ACC_ID, MOD_RSD) |> 
  separate(MOD_RSD, into = c('Site', 'gl'), sep = '-') |> 
  mutate(
    site_position = as.numeric(str_extract(Site, '\\d+'))
  ) |> 
  select(UniProt_Accession, site_position) |> 
  distinct()

# integrate the Nterm half life results with O-GlcNAcylation site dataset
Nterm_OGlcNAc_site_occurrence <- OGlcNAc_site_HEK293T |> 
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

Nterm_OGlcNAc_site_occurrence |> 
  wilcox_test(n ~ category)

# Spearman correlation test
cor.test(
  Nterm_OGlcNAc_site_occurrence$half_life,
  Nterm_OGlcNAc_site_occurrence$n,
  method = 'spearman'
)
