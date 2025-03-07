# import packages
library(tidyverse)

# linear model
linear_model <- function(df) {
  tryCatch({
    # Fit the model
    model <- lm(
      ln_deg_ratio_avg ~ timepoint,
      data = df
    )
    
    # Extract coefficients
    params <- coef(model)
    
    # Calculate rss
    ss_residual <- sum(residuals(model)^2)
    
    # Calculate R^2
    r_squared <- summary(model)$r.squared
    
    # Return coefficients and R^2
    tibble(lnA = params["(Intercept)"], Kd = params["timepoint"], R2_RSS = ss_residual)
  }, error = function(e) {
    # Return NA if fitting fails
    tibble(lnA = NA, Kd = NA, R2_RSS = NA)
  })
}

# calculate cell doubling rate according to BCA result
Rep1_BCA <- tribble(
  ~ Exp, ~ timepoint, ~ concentration, ~ ratio,
  'Rep1', 0, 3.73, 3.73/3.73,
  'Rep1', 3, 3.40, 3.40/3.73,
  'Rep1', 6, 4.15, 4.15/3.73,
  'Rep1', 9, 4.46, 4.46/3.73,
  'Rep1', 12, 5.15, 5.15/3.73,
  'Rep1', 24, 7.75, 7.75/3.73
)

Rep2_BCA <- tribble(
  ~ Exp, ~ timepoint, ~ concentration, ~ ratio,
  'Rep2', 0, 2.83, 2.83/2.83,
  'Rep2', 3, 3.21, 3.21/2.83,
  'Rep2', 6, 4.41, 4.41/2.83,
  'Rep2', 9, 4.28, 4.28/2.83,
  'Rep2', 12, 4.81, 4.81/2.83,
  'Rep2', 24, 7.58, 7.58/2.83
)

Rep3_BCA <- tribble(
  ~ Exp, ~ timepoint, ~ concentration, ~ ratio,
  'Rep3', 0, 1.63, 1.63/1.63,
  'Rep3', 3, 2.00, 2.00/1.63,
  'Rep3', 6, 2.17, 2.17/1.63,
  'Rep3', 9, 2.49, 2.49/1.63,
  'Rep3', 12, 3.17, 3.17/1.63,
  'Rep3', 24, 4.75, 4.75/1.63
)

cell_double_rate_kinetic <- bind_rows(
  Rep1_BCA, Rep2_BCA, Rep3_BCA
) |> 
  select(Exp, timepoint, ratio) |> 
  mutate(
    ln_deg_ratio_avg = log(ratio)
  ) |> 
  group_by(Exp) |> 
  group_modify(~ linear_model(.x)) |> 
  ungroup()

K_cd <- mean(cell_double_rate_kinetic$Kd)

# Nterm half life by cell doubling normalization
HEK_Nterm_Kd_half_life_cell_doubling <- HEK_Nterm_curve_fitting_combined |> 
  filter(parameters %in% c('Kd', 'RSS')) |> 
  pivot_wider(names_from = `parameters`, values_from = `values`) |> 
  filter(Kd < 5) |> 
  mutate(
    half_life = log(2)/(Kd - K_cd)
  )

# calculate normalization factor using half life of Lamin B and T-complex
HEK_Nterm_LaminB_half_life <- HEK_Nterm_Kd_half_life_cell_doubling |> 
  filter(UniProt_Accession %in% c(
    'P20700', #	Lamin-B1
    'Q03252' # Lamin-B2
  ))

HEK_Nterm_Tcomplex_half_life <- HEK_Nterm_Kd_half_life_cell_doubling |> 
  filter(UniProt_Accession %in% c(
    'P48643', # T-complex protein 1 subunit epsilon
    'P78371', # T-complex protein 1 subunit beta
    'Q99832', # T-complex protein 1 subunit eta
    'P49368', # T-complex protein 1 subunit gamma
    'P17987', # T-complex protein 1 subunit alpha
    'P50991', # T-complex protein 1 subunit delta
    'P50990', # T-complex protein 1 subunit theta
    'P40227', # T-complex protein 1 subunit zeta
    'Q92526' # T-complex protein 1 subunit zeta-2
  ))

Nterm_LaminB_Tcomplex_combined <- bind_rows(
  HEK_Nterm_LaminB_half_life, 
  HEK_Nterm_Tcomplex_half_life
)

library(rstatix)

Nterm_LaminB_Tcomplex_half_life_filter_criteria <- Nterm_LaminB_Tcomplex_combined |> 
  filter(half_life < 0) |> 
  get_summary_stats(half_life, type = 'median')

# Nterm half life normalization by Lamin B and T-complex
HEK_Nterm_Kd_half_life_LaminB_Tcomplex <- HEK_Nterm_Kd_half_life_cell_doubling |> 
  filter(half_life > Nterm_LaminB_Tcomplex_half_life_filter_criteria |> pull(median)) |> 
  mutate(
    half_life = ifelse(half_life < 0 | half_life > 200, 200, half_life)
  ) |> 
  # count(half_life == 200)
  # get_summary_stats(half_life)
  select(Index, UniProt_Accession, Protein.Start, Gene, Entry.Name, Kd, half_life, RSS) |> 
  mutate(
    Percentile = percent_rank(half_life)
  ) |> 
  arrange(Percentile) |> 
  mutate(
    category = case_when(
      half_life < 6 ~ 'Fast turnover',
      half_life == 200 ~ 'Stable',
      .default = 'Median'
    )
  )

write_csv(HEK_Nterm_Kd_half_life_LaminB_Tcomplex, file = 'data_source/Kd_half_life/HEK_Nterm_Kd_half_life_LaminB_Tcomplex.csv')

# Whole Proteome half life by cell doubling normalization
HEK_WP_Kd_half_life_cell_doubling <- HEK_WP_curve_fitting_combined |> 
  filter(parameters %in% c('Kd', 'RSS')) |> 
  pivot_wider(names_from = `parameters`, values_from = `values`) |> 
  filter(Kd < 5) |> 
  mutate(
    half_life = log(2)/(Kd - K_cd)
  )

# calculate normalization factor using half life of Lamin B and T-complex
HEK_WP_LaminB_half_life <- HEK_WP_Kd_half_life_cell_doubling |> 
  filter(UniProt_Accession %in% c(
    'P20700', #	Lamin-B1
    'Q03252' # Lamin-B2
  ))

HEK_WP_Tcomplex_half_life <- HEK_WP_Kd_half_life_cell_doubling |> 
  filter(UniProt_Accession %in% c(
    'P48643', # T-complex protein 1 subunit epsilon
    'P78371', # T-complex protein 1 subunit beta
    'Q99832', # T-complex protein 1 subunit eta
    'P49368', # T-complex protein 1 subunit gamma
    'P17987', # T-complex protein 1 subunit alpha
    'P50991', # T-complex protein 1 subunit delta
    'P50990', # T-complex protein 1 subunit theta
    'P40227', # T-complex protein 1 subunit zeta
    'Q92526' # T-complex protein 1 subunit zeta-2
  ))

Whole_Proteome_LaminB_Tcomplex_combined <- bind_rows(
  HEK_WP_LaminB_half_life, 
  HEK_WP_Tcomplex_half_life
)

library(rstatix)

Whole_Proteome_LaminB_Tcomplex_half_life_filter_criteria <- Whole_Proteome_LaminB_Tcomplex_combined |> 
  filter(half_life < 0) |> 
  get_summary_stats(half_life, type = 'median')

# Whole Proteome half life normalization by Lamin B and T-complex
HEK_WP_Kd_half_life_LaminB_Tcomplex <- HEK_WP_Kd_half_life_cell_doubling |> 
  filter(half_life > Whole_Proteome_LaminB_Tcomplex_half_life_filter_criteria |> pull(median)) |> 
  mutate(
    half_life = ifelse(half_life < 0 | half_life > 200, 200, half_life)
  ) |> 
  # count(half_life == 200)
  # get_summary_stats(half_life)
  select(UniProt_Accession, Gene, Entry.Name, Kd, half_life, RSS) |> 
  mutate(
    Percentile = percent_rank(half_life)
  ) |> 
  arrange(Percentile)

write_csv(HEK_WP_Kd_half_life_LaminB_Tcomplex, file = 'data_source/Kd_half_life/HEK_WP_Kd_half_life_LaminB_Tcomplex.csv')
