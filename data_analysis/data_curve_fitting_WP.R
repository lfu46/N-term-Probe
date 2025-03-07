# import packages
library(tidyverse)

## calculate degradation ratio
# HEK_WP_1
HEK_WP_1_deg_ratio <- HEK_Whole_Proteome_1_irs |> 
  mutate(
    ratio_0 = deg_126_0h_TMTi_irs/deg_126_0h_TMTi_irs,
    ratio_3 = deg_127_3h_TMTi_irs/deg_126_0h_TMTi_irs,
    ratio_6 = deg_128_6h_TMTi_irs/deg_126_0h_TMTi_irs,
    ratio_9 = deg_129_9h_TMTi_irs/deg_126_0h_TMTi_irs,
    ratio_12 = deg_130_12h_TMTi_irs/deg_126_0h_TMTi_irs,
    ratio_24 = deg_131_24h_TMTi_irs/deg_126_0h_TMTi_irs
  ) |> 
  select(UniProt_Accession, Gene, Entry.Name, starts_with("ratio")) |> 
  pivot_longer(cols = starts_with("ratio"), names_to = 'deg_timepoint', values_to = 'deg_ratio') |> 
  separate(col = deg_timepoint, sep = "_", into = c('ratio', 'timepoint'), convert = TRUE) |> 
  select(UniProt_Accession, Gene, Entry.Name, timepoint, deg_ratio)

write_csv(HEK_WP_1_deg_ratio, file = 'data_source/degradation_ratio/HEK_WP_1_deg_ratio.csv')

# HEK_WP_2
HEK_WP_2_deg_ratio <- HEK_Whole_Proteome_2_irs |> 
  mutate(
    ratio_0 = deg_126_0h_TMTi_irs/deg_126_0h_TMTi_irs,
    ratio_3 = deg_127_3h_TMTi_irs/deg_126_0h_TMTi_irs,
    ratio_6 = deg_128_6h_TMTi_irs/deg_126_0h_TMTi_irs,
    ratio_9 = deg_129_9h_TMTi_irs/deg_126_0h_TMTi_irs,
    ratio_12 = deg_130_12h_TMTi_irs/deg_126_0h_TMTi_irs,
    ratio_24 = deg_131_24h_TMTi_irs/deg_126_0h_TMTi_irs
  ) |> 
  select(UniProt_Accession, Gene, Entry.Name, starts_with("ratio")) |> 
  pivot_longer(cols = starts_with("ratio"), names_to = 'deg_timepoint', values_to = 'deg_ratio') |> 
  separate(col = deg_timepoint, sep = "_", into = c('ratio', 'timepoint'), convert = TRUE) |> 
  select(UniProt_Accession, Gene, Entry.Name, timepoint, deg_ratio)

write_csv(HEK_WP_2_deg_ratio, file = 'data_source/degradation_ratio/HEK_WP_2_deg_ratio.csv')

# HEK_WP_3
HEK_WP_3_deg_ratio <- HEK_Whole_Proteome_3_irs |> 
  mutate(
    ratio_0 = deg_126_0h_TMTi_irs/deg_126_0h_TMTi_irs,
    ratio_3 = deg_127_3h_TMTi_irs/deg_126_0h_TMTi_irs,
    ratio_6 = deg_128_6h_TMTi_irs/deg_126_0h_TMTi_irs,
    ratio_9 = deg_129_9h_TMTi_irs/deg_126_0h_TMTi_irs,
    ratio_12 = deg_130_12h_TMTi_irs/deg_126_0h_TMTi_irs,
    ratio_24 = deg_131_24h_TMTi_irs/deg_126_0h_TMTi_irs
  ) |> 
  select(UniProt_Accession, Gene, Entry.Name, starts_with("ratio")) |> 
  pivot_longer(cols = starts_with("ratio"), names_to = 'deg_timepoint', values_to = 'deg_ratio') |> 
  separate(col = deg_timepoint, sep = "_", into = c('ratio', 'timepoint'), convert = TRUE) |> 
  select(UniProt_Accession, Gene, Entry.Name, timepoint, deg_ratio)

write_csv(HEK_WP_3_deg_ratio, file = 'data_source/degradation_ratio/HEK_WP_3_deg_ratio.csv')

## get overlap of quantified protein from each replicate
# overlap in three replicates
overlap_protein <- intersect(
  intersect(HEK_WP_1_deg_ratio$UniProt_Accession, HEK_WP_2_deg_ratio$UniProt_Accession), 
  HEK_WP_3_deg_ratio$UniProt_Accession
)

# overlap in HEK_WP_1 and HEK_WP_2
overlap_1_2 <- setdiff(
  intersect(HEK_WP_1_deg_ratio$UniProt_Accession, HEK_WP_2_deg_ratio$UniProt_Accession), 
  overlap_protein
)

# overlap in HEK_WP_2 and HEK_WP_3
overlap_2_3 <- setdiff(
  intersect(HEK_WP_2_deg_ratio$UniProt_Accession, HEK_WP_3_deg_ratio$UniProt_Accession), 
  overlap_protein
)

# overlap in HEK_WP_3 and HEK_WP_1
overlap_3_1 <- setdiff(
  intersect(HEK_WP_3_deg_ratio$UniProt_Accession, HEK_WP_1_deg_ratio$UniProt_Accession), 
  overlap_protein
)

# unique in HEK_WP_1
unique_1 <- setdiff(
  HEK_WP_1_deg_ratio$UniProt_Accession, 
  union(HEK_WP_2_deg_ratio$UniProt_Accession, HEK_WP_3_deg_ratio$UniProt_Accession)
)

# unique in HEK_WP_2
unique_2 <- setdiff(
  HEK_WP_2_deg_ratio$UniProt_Accession, 
  union(HEK_WP_3_deg_ratio$UniProt_Accession, HEK_WP_1_deg_ratio$UniProt_Accession)
)

# unique in HEK_WP_3
unique_3 <- setdiff(
  HEK_WP_3_deg_ratio$UniProt_Accession, 
  union(HEK_WP_1_deg_ratio$UniProt_Accession, HEK_WP_2_deg_ratio$UniProt_Accession)
)

## combine results from three replicates
# overlap_protein
overlap_WP_deg_ratio <- HEK_WP_1_deg_ratio |> 
  filter(UniProt_Accession %in% overlap_protein) |> 
  left_join(
    HEK_WP_2_deg_ratio |> filter(UniProt_Accession %in% overlap_protein),
    by = c('UniProt_Accession', 'Gene', 'Entry.Name', 'timepoint'),
    suffix = c('.WP_1', '.WP_2')
  ) |> 
  left_join(
    HEK_WP_3_deg_ratio |> filter(UniProt_Accession %in% overlap_protein),
    by = c('UniProt_Accession', 'Gene', 'Entry.Name', 'timepoint')
  ) |> 
  mutate(
    deg_ratio_avg = (deg_ratio.WP_1 + deg_ratio.WP_2 + deg_ratio)/3
  ) |> 
  select(UniProt_Accession, Gene, Entry.Name, timepoint, deg_ratio_avg)

# overlap in HEK_WP_1 and HEK_WP_2
overlap_1_2_deg_ratio <- HEK_WP_1_deg_ratio |> 
  filter(UniProt_Accession %in% overlap_1_2) |> 
  left_join(
    HEK_WP_2_deg_ratio |> filter(UniProt_Accession %in% overlap_1_2),
    by = c('UniProt_Accession', 'Gene', 'Entry.Name', 'timepoint'),
    suffix = c('.WP_1', '.WP_2')
  ) |> 
  mutate(
    deg_ratio_avg = (deg_ratio.WP_1 + deg_ratio.WP_2)/2
  ) |> 
  select(UniProt_Accession, Gene, Entry.Name, timepoint, deg_ratio_avg)

# overlap in HEK_WP_2 and HEK_WP_3
overlap_2_3_deg_ratio <- HEK_WP_2_deg_ratio |> 
  filter(UniProt_Accession %in% overlap_2_3) |> 
  left_join(
    HEK_WP_3_deg_ratio |> filter(UniProt_Accession %in% overlap_2_3),
    by = c('UniProt_Accession', 'Gene', 'Entry.Name', 'timepoint'),
    suffix = c('.WP_2', '.WP_3')
  ) |> 
  mutate(
    deg_ratio_avg = (deg_ratio.WP_2 + deg_ratio.WP_3)/2
  ) |> 
  select(UniProt_Accession, Gene, Entry.Name, timepoint, deg_ratio_avg)

# overlap in HEK_WP_3 and HEK_WP_1
overlap_3_1_deg_ratio <- HEK_WP_3_deg_ratio |> 
  filter(UniProt_Accession %in% overlap_3_1) |> 
  left_join(
    HEK_WP_1_deg_ratio |> filter(UniProt_Accession %in% overlap_3_1),
    by = c('UniProt_Accession', 'Gene', 'Entry.Name', 'timepoint'),
    suffix = c('.WP_3', '.WP_1')
  ) |> 
  mutate(
    deg_ratio_avg = (deg_ratio.WP_3 + deg_ratio.WP_1)/2
  ) |> 
  select(UniProt_Accession, Gene, Entry.Name, timepoint, deg_ratio_avg)

# unique in HEK_WP_1
unique_1_deg_ratio <- HEK_WP_1_deg_ratio |> 
  filter(UniProt_Accession %in% unique_1) |> 
  mutate(
    deg_ratio_avg = deg_ratio
  ) |> 
  select(UniProt_Accession, Gene, Entry.Name, timepoint, deg_ratio_avg)

# unique in HEK_WP_2
unique_2_deg_ratio <- HEK_WP_2_deg_ratio |> 
  filter(UniProt_Accession %in% unique_2) |> 
  mutate(
    deg_ratio_avg = deg_ratio
  ) |> 
  select(UniProt_Accession, Gene, Entry.Name, timepoint, deg_ratio_avg)

# unique in HEK_WP_3
unique_3_deg_ratio <- HEK_WP_3_deg_ratio |> 
  filter(UniProt_Accession %in% unique_3) |> 
  mutate(
    deg_ratio_avg = deg_ratio
  ) |> 
  select(UniProt_Accession, Gene, Entry.Name, timepoint, deg_ratio_avg)

# combine results from three replicates
HEK_WP_deg_ratio <- bind_rows(
  overlap_WP_deg_ratio,
  overlap_1_2_deg_ratio,
  overlap_2_3_deg_ratio,
  overlap_3_1_deg_ratio,
  unique_1_deg_ratio,
  unique_2_deg_ratio,
  unique_3_deg_ratio
)

write_csv(HEK_WP_deg_ratio, file = 'data_source/degradation_ratio/HEK_WP_deg_ratio.csv')

## curve fitting models
# non-linear model
non_linear_model <- function(df) {
  tryCatch({
    # Fit the model
    model <- nls(
      deg_ratio_avg ~ (A - B) * exp(-Kd * timepoint) + B,
      data = df,
      # A refers to the maximum of the curve and should be 1 in an ideal case
      # B accounts for a potential curve offset which ideally should be 0
      # the initial Kd value can be fine-tuned to achieve better results
      start = list(A = 1, B = 0, Kd = 0.015)
    )
    
    # Extract coefficients
    params <- coef(model)
    
    # Calculate R^2
    fitted_values <- predict(model, df)
    ss_total <- sum((df$deg_ratio_avg - mean(df$deg_ratio_avg))^2)
    ss_residual <- sum((df$deg_ratio_avg - fitted_values)^2)
    r_squared <- 1 - (ss_residual / ss_total)
    
    # Return coefficients and R^2
    tibble(A = params["A"], B = params["B"], Kd = params["Kd"], RSS = ss_residual)
  }, error = function(e) {
    # Return NA if fitting fails
    tibble(A = NA, B = NA, Kd = NA, RSS = NA)
  })
}

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

## curve fitting
# non-linear model fitting
HEK_WP_curve_fitting_result_1 <- HEK_WP_deg_ratio %>%
  group_by(UniProt_Accession, Gene, Entry.Name) |> 
  group_modify(~ non_linear_model(.x)) |> 
  ungroup()

HEK_WP_nonlinear_fitting <- HEK_WP_curve_fitting_result_1 |> 
  filter(RSS < 0.05) |> 
  select(UniProt_Accession, Gene, Entry.Name, A:RSS) |> 
  pivot_longer(cols = A:RSS, names_to = 'parameters', values_to = 'values') |> 
  mutate(
    model = 'nonlinear fitting'
  )

write_csv(HEK_WP_nonlinear_fitting, file = 'data_source/curve_fitting/HEK_WP_nonlinear_fitting.csv')

# linear model fitting
HEK_WP_curve_fitting_result_2 <- HEK_WP_curve_fitting_result_1 |> 
  filter(is.na(A) | RSS >= 0.05) |> 
  left_join(HEK_WP_deg_ratio, by = c('UniProt_Accession', 'Gene', 'Entry.Name')) |> 
  select(UniProt_Accession, Gene, Entry.Name, timepoint, deg_ratio_avg) |> 
  mutate(
    ln_deg_ratio_avg = log(deg_ratio_avg)
  ) |> 
  group_by(UniProt_Accession, Gene, Entry.Name) |> 
  group_modify(~ linear_model(.x)) |> 
  ungroup()

HEK_WP_linear_fitting <- HEK_WP_curve_fitting_result_2 |> 
  filter(RSS < 0.1) |> 
  select(UniProt_Accession, Gene, Entry.Name, lnA:RSS) |> 
  pivot_longer(cols = lnA:RSS, names_to = 'parameters', values_to = 'values') |> 
  mutate(
    model = 'linear fitting'
  )

write_csv(HEK_WP_linear_fitting, file = 'data_source/curve_fitting/HEK_WP_linear_fitting.csv')

# combine non-linear fitting and linear fitting result
HEK_WP_curve_fitting_combined <- bind_rows(
  HEK_WP_nonlinear_fitting,
  HEK_WP_linear_fitting
)

write_csv(HEK_WP_curve_fitting_combined, file = 'data_source/curve_fitting/HEK_WP_curve_fitting_combined.csv')
