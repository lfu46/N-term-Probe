# import packages
library(tidyverse)

# Nterm half life by cell doubling normalization
library(rstatix)

HEK_Nterm_Kd_half_life_LaminB_Tcomplex <- HEK_Nterm_linear_fitting |> 
  filter(parameters %in% c('Kd', 'RSS')) |> 
  pivot_wider(names_from = `parameters`, values_from = `values`) |> 
  mutate(
    Kd_adj = Kd - 0.03,
    half_life = log(2)/Kd_adj
  ) |> 
  filter(half_life > 0) |>
  mutate(
    half_life = ifelse(half_life < 0 | half_life > 200, 200, half_life)
  ) |> 
  # count(half_life == 200)
  # get_summary_stats(half_life)
  select(Index, UniProt_Accession, Protein.Start, Gene, Entry.Name, Kd_adj, half_life, RSS) |> 
  mutate(
    Percentile = percent_rank(half_life)
  ) |> 
  arrange(Percentile) |> 
  mutate(
    category = case_when(
      Percentile < 0.2 ~ 'Fast turnover',
      Percentile > 0.8 ~ 'Stable',
      .default = 'Median'
    )
  )

write_csv(HEK_Nterm_Kd_half_life_LaminB_Tcomplex, file = 'data_source/Kd_half_life/HEK_Nterm_Kd_half_life_LaminB_Tcomplex.csv')

# Whole Proteome half life by cell doubling normalization
library(rstatix)

HEK_WP_Kd_half_life_LaminB_Tcomplex <- HEK_WP_linear_fitting |> 
  filter(parameters %in% c('Kd', 'RSS')) |> 
  pivot_wider(names_from = `parameters`, values_from = `values`) |> 
  mutate(
    Kd_adj = Kd - 0.03,
    half_life = log(2)/Kd_adj
  ) |> 
  mutate(
    half_life = ifelse(half_life < 0 | half_life > 200, 200, half_life)
  ) |> 
  # count(half_life == 200) |>
  # get_summary_stats(half_life) |> 
  select(UniProt_Accession, Gene, Entry.Name, Kd_adj, half_life, RSS) |> 
  mutate(
    Percentile = percent_rank(half_life)
  ) |> 
  arrange(Percentile)

write_csv(HEK_WP_Kd_half_life_LaminB_Tcomplex, file = 'data_source/Kd_half_life/HEK_WP_Kd_half_life_LaminB_Tcomplex.csv')
