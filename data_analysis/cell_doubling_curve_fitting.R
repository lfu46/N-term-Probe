# import packages
library(tidyverse)

## cell doubling kinetics
# linear model
linear_model_cell_doubling <- function(df) {
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
  group_modify(~ linear_model_cell_doubling(.x)) |> 
  ungroup()

