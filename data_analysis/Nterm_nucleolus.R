# import packages
packages_names <- c("tidyverse", 'readxl', 'rstatix')
lapply(packages_names, require, character.only = TRUE)

# import data for nucleolus protein sub-localization
nucleolus_localization <- read_xlsx(
  'data_source/nucleolus/41586_2023_5767_MOESM4_ESM.xlsx',
  skip = 1,
  col_names = TRUE,
  .name_repair = 'universal'
) |> 
  select(Gene.Symbol = Gene.Symbol.and.Synonyms, Localization = Localization.in.this.paper) |> 
  separate_rows(Gene.Symbol, Localization, sep = ', ')

# Nterm nucleolus half life
Nterm_nucleolus_localization_half_life <- nucleolus_localization |> 
  left_join(HEK_Nterm_Kd_half_life_LaminB_Tcomplex, by = join_by('Gene.Symbol' == 'Gene'), relationship = 'many-to-many') |> 
  filter(!is.na(half_life))

# Wilcoxon rank-sum test
Nterm_nucleolus_localization_half_life |> 
  wilcox_test(half_life ~ Localization) |> 
  filter(p < 0.05)
