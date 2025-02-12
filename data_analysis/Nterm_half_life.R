#import packages
packages_names <- c("tidyverse", "readxl", "writexl", 'rstatix')
lapply(packages_names, require, character.only = TRUE)

#linear model
linear_model <- function(df) {
  tryCatch({
    # Fit the model
    model <- lm(
      ln_ratio ~ timepoint,
      data = df
    )
    
    # Extract coefficients
    params <- coef(model)
    
    # Calculate rss
    ss_residual <- sum(residuals(model)^2)
    
    # Calculate R^2
    r_squared <- summary(model)$r.squared
    
    # Return coefficients and R^2
    tibble(lnA = params["(Intercept)"], Kd = params["timepoint"], RSS = ss_residual)
  }, error = function(e) {
    # Return NA if fitting fails
    tibble(lnA = NA, Kd = NA, RSS = NA)
  })
}

#calculate cell doubling rate according to BCA result
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
  mutate(ln_ratio = log(ratio)) |> 
  group_by(Exp) |> 
  group_modify(~ linear_model(.x)) |> 
  ungroup()

K_cd <- mean(cell_double_rate_kinetic$Kd)

#calculate half life
HEK_Nt_Kd_half_life <- HEK_Nterm_curve_fitting_combined |> 
  filter(parameters %in% c('Kd', 'RSS')) |> 
  pivot_wider(names_from = `parameters`, values_from = `values`) |> 
  mutate(
    half_life = log(2)/(Kd - K_cd)
  ) |> 
  filter(half_life > 0) |> 
  select(Index, UniProt_Accession, Gene, Entry.Name, Kd, half_life, RSS) |> 
  mutate(
    Percentile = percent_rank(half_life)
  ) |> 
  arrange(Percentile)

HEK_Nt_Kd_half_life |> 
  ggplot() +
  geom_histogram(aes(x = half_life), bins = 50) +
  xlim(0, 100)

HEK_Nt_Kd_half_life |> 
  get_summary_stats(half_life)
