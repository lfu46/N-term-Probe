# import packages
library(tidyverse)

## curve fitting for Lamin B and T-complex proteoforms
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
    tibble(lnA = params["(Intercept)"], Kd = -params["timepoint"], RSS = ss_residual)
  }, error = function(e) {
    # Return NA if fitting fails
    tibble(lnA = NA, Kd = NA, RSS = NA)
  })
}

# curve fitting
HEK_WP_curve_fitting_LaminB_Tcomplex <- HEK_WP_deg_ratio |> 
  filter(
    UniProt_Accession %in% c(
      'P20700', #	Lamin-B1
      'Q03252', # Lamin-B2
      'P48643', # T-complex protein 1 subunit epsilon
      'P78371', # T-complex protein 1 subunit beta
      'Q99832', # T-complex protein 1 subunit eta
      'P49368', # T-complex protein 1 subunit gamma
      'P17987', # T-complex protein 1 subunit alpha
      'P50991', # T-complex protein 1 subunit delta
      'P50990', # T-complex protein 1 subunit theta
      'P40227', # T-complex protein 1 subunit zeta
      'Q92526' # T-complex protein 1 subunit zeta-2
    )
  ) |> 
  select(UniProt_Accession, Gene, Entry.Name, timepoint, deg_ratio_avg) |> 
  mutate(
    ln_deg_ratio_avg = log(deg_ratio_avg)
  ) |> 
  group_by(UniProt_Accession, Gene, Entry.Name) |> 
  group_modify(~ linear_model(.x)) |> 
  ungroup()

# filter
HEK_WP_LaminB_Tcomplex_result <- HEK_WP_curve_fitting_LaminB_Tcomplex |> 
  filter(RSS < 0.05) |> 
  select(UniProt_Accession, Gene, Entry.Name, lnA:RSS) |> 
  pivot_longer(cols = lnA:RSS, names_to = 'parameters', values_to = 'values') |> 
  mutate(
    model = 'linear fitting'
  )

write_csv(HEK_WP_LaminB_Tcomplex_result, file = 'data_source/LaminB_Tcomplex_normalization/HEK_WP_LaminB_Tcomplex_result.csv')

## half life calculation for Lamin B and T-complex proteoforms
HEK_WP_LaminB_Tcomplex_half_life <- HEK_WP_LaminB_Tcomplex_result |> 
  pivot_wider(names_from = `parameters`, values_from = `values`) |> 
  mutate(
    Kd_adj = Kd - 0.055,
    half_life = log(2)/Kd_adj
  ) |> 
  filter(half_life < 0)

write_csv(HEK_WP_LaminB_Tcomplex_half_life, file = 'data_source/LaminB_Tcomplex_normalization/HEK_WP_LaminB_Tcomplex_half_life.csv')

# ratio calculation for each time point
timepoint <- c(0, 3, 6, 9, 12, 24)
HEK_WP_LaminB_Tcomplex_adjusted_ratio <- HEK_WP_LaminB_Tcomplex_half_life |> 
  select(UniProt_Accession, Gene, Entry.Name, lnA, Kd, Kd_adj) |>
  crossing(timepoint) |>
  mutate(
    ratio = exp(lnA) * exp(-Kd_adj * timepoint)
  )

write_csv(
  HEK_WP_LaminB_Tcomplex_adjusted_ratio, 
  file = 'data_source/LaminB_Tcomplex_normalization/HEK_WP_LaminB_Tcomplex_adjusted_ratio.csv'
)

# ratio normalization factor
ratio_norm_facs_WP <- HEK_WP_LaminB_Tcomplex_adjusted_ratio |> 
  group_by(timepoint) |> 
  summarize(
    ratio_norm_fac = max(ratio)
  )

write_csv(ratio_norm_facs_WP, file = 'data_source/LaminB_Tcomplex_normalization/ratio_norm_facs_WP.csv')

