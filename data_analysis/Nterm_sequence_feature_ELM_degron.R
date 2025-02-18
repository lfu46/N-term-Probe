# import packages
packages_names <- c("tidyverse", "rstatix")
lapply(packages_names, require, character.only = TRUE)

# ELM N-degron
HEK_Nterm_ELM_N_degron <- HEK_Nterm_Kd_half_life_sequence |> 
  mutate(
    ELM_N_degron = case_when(
      Nterm_terminus %in% c('F', 'Y', 'L', 'I', 'W') ~ 'FYLIW',
      Nterm_terminus %in% c('R', 'K') ~ 'RK',
      Nterm_terminus %in% c('E', 'D') ~ 'ED',
      Nterm_terminus %in% c('N', 'Q') ~ 'NQ',
      Nterm_terminus %in% c('C') ~ 'C',
      Nterm_terminus %in% c('M') ~ 'M'
    )  
  ) |> 
  mutate(
    ELM_N_degron = ifelse(is.na(ELM_N_degron), 'Others', ELM_N_degron)
  )

# ELM N-degron half life median
ELM_N_degron_half_life_median <- HEK_Nterm_ELM_N_degron |> 
  group_by(ELM_N_degron) |> 
  get_summary_stats(half_life, type = 'median') |> 
  arrange(median)

write_csv(ELM_N_degron_half_life_median, file = 'data_source/ELM_degron/ELM_N_degron_half_life_median.csv')

# Wilcoxon rank-sum test
HEK_Nterm_ELM_N_degron |>
  wilcox_test(half_life ~ ELM_N_degron) |> 
  filter(p < 0.05)

