# import packages
library(tidyverse)

## calculate degradation ratio
# HEK_Nt_1
HEK_Nterm_1_deg_ratio <- HEK_Nterm_1_irs |> 
  mutate(
    ratio_0 = deg_126_0h_TMTi_irs/deg_126_0h_TMTi_irs,
    ratio_3 = deg_127_3h_TMTi_irs/deg_126_0h_TMTi_irs,
    ratio_6 = deg_128_6h_TMTi_irs/deg_126_0h_TMTi_irs,
    ratio_9 = deg_129_9h_TMTi_irs/deg_126_0h_TMTi_irs,
    ratio_12 = deg_130_12h_TMTi_irs/deg_126_0h_TMTi_irs,
    ratio_24 = deg_131_24h_TMTi_irs/deg_126_0h_TMTi_irs
  ) |> 
  select(Index, UniProt_Accession, Protein.Start, Gene, Entry.Name, starts_with("ratio")) |> 
  pivot_longer(cols = starts_with("ratio"), names_to = 'deg_timepoint', values_to = 'deg_ratio') |> 
  separate(col = deg_timepoint, sep = "_", into = c('ratio', 'timepoint'), convert = TRUE) |> 
  select(Index, UniProt_Accession, Protein.Start, Gene, Entry.Name, timepoint, deg_ratio)

write_csv(HEK_Nterm_1_deg_ratio, file = 'data_source/degradation_ratio/HEK_Nterm_1_deg_ratio.csv')

# HEK_Nt_2
HEK_Nterm_2_deg_ratio <- HEK_Nterm_2_irs |> 
  mutate(
    ratio_0 = deg_126_0h_TMTi_irs/deg_126_0h_TMTi_irs,
    ratio_3 = deg_127_3h_TMTi_irs/deg_126_0h_TMTi_irs,
    ratio_6 = deg_128_6h_TMTi_irs/deg_126_0h_TMTi_irs,
    ratio_9 = deg_129_9h_TMTi_irs/deg_126_0h_TMTi_irs,
    ratio_12 = deg_130_12h_TMTi_irs/deg_126_0h_TMTi_irs,
    ratio_24 = deg_131_24h_TMTi_irs/deg_126_0h_TMTi_irs
  ) |> 
  select(Index, UniProt_Accession, Protein.Start, Gene, Entry.Name, starts_with("ratio")) |> 
  pivot_longer(cols = starts_with("ratio"), names_to = 'deg_timepoint', values_to = 'deg_ratio') |> 
  separate(col = deg_timepoint, sep = "_", into = c('ratio', 'timepoint'), convert = TRUE) |> 
  select(Index, UniProt_Accession, Protein.Start, Gene, Entry.Name, timepoint, deg_ratio)

write_csv(HEK_Nterm_2_deg_ratio, file = 'data_source/degradation_ratio/HEK_Nterm_2_deg_ratio.csv')

# HEK_Nt_3
HEK_Nterm_3_deg_ratio <- HEK_Nterm_3_irs |> 
  mutate(
    ratio_0 = deg_126_0h_TMTi_irs/deg_126_0h_TMTi_irs,
    ratio_3 = deg_127_3h_TMTi_irs/deg_126_0h_TMTi_irs,
    ratio_6 = deg_128_6h_TMTi_irs/deg_126_0h_TMTi_irs,
    ratio_9 = deg_129_9h_TMTi_irs/deg_126_0h_TMTi_irs,
    ratio_12 = deg_130_12h_TMTi_irs/deg_126_0h_TMTi_irs,
    ratio_24 = deg_131_24h_TMTi_irs/deg_126_0h_TMTi_irs
  ) |> 
  select(Index, UniProt_Accession, Protein.Start, Gene, Entry.Name, starts_with("ratio")) |> 
  pivot_longer(cols = starts_with("ratio"), names_to = 'deg_timepoint', values_to = 'deg_ratio') |> 
  separate(col = deg_timepoint, sep = "_", into = c('ratio', 'timepoint'), convert = TRUE) |> 
  select(Index, UniProt_Accession, Protein.Start, Gene, Entry.Name, timepoint, deg_ratio)

write_csv(HEK_Nterm_3_deg_ratio, file = 'data_source/degradation_ratio/HEK_Nterm_3_deg_ratio.csv')

## get overlap quantified N-term from each replicate
# overlap in three replicates
overlap_Nterm <- intersect(
  intersect(HEK_Nterm_1_deg_ratio$Index, HEK_Nterm_2_deg_ratio$Index), HEK_Nterm_3_deg_ratio$Index
)

# overlap in HEK_Nt_1 and HEK_Nt_2
overlap_1_2 <- setdiff(intersect(HEK_Nterm_1_deg_ratio$Index, HEK_Nterm_2_deg_ratio$Index), overlap_Nterm)

#overlap in HEK_Nt_2 and HEK_Nt_3
overlap_2_3 <- setdiff(intersect(HEK_Nterm_2_deg_ratio$Index, HEK_Nterm_3_deg_ratio$Index), overlap_Nterm)

# overlap in HEK_Nt_3 and HEK_Nt_1
overlap_3_1 <- setdiff(intersect(HEK_Nterm_3_deg_ratio$Index, HEK_Nterm_1_deg_ratio$Index), overlap_Nterm)

# unique in HEK_Nt_1
unique_1 <- setdiff(HEK_Nterm_1_deg_ratio$Index, union(HEK_Nterm_2_deg_ratio$Index, HEK_Nterm_3_deg_ratio$Index))

# unique in HEK_Nt_2
unique_2 <- setdiff(HEK_Nterm_2_deg_ratio$Index, union(HEK_Nterm_3_deg_ratio$Index, HEK_Nterm_1_deg_ratio$Index))

# uniquqe in HEK_Nt_3
unique_3 <- setdiff(HEK_Nterm_3_deg_ratio$Index, union(HEK_Nterm_1_deg_ratio$Index, HEK_Nterm_2_deg_ratio$Index))

## combine results from three replicates
# overlap_Nterm
overlap_Nterm_deg_ratio <- HEK_Nterm_1_deg_ratio |> 
  filter(Index %in% overlap_Nterm) |> 
  left_join(
    HEK_Nterm_2_deg_ratio |> filter(Index %in% overlap_Nterm),
    by = c('Index', 'UniProt_Accession', 'Protein.Start', 'Gene', 'Entry.Name', 'timepoint'),
    suffix = c('.Nt_1', '.Nt_2')
  ) |> 
  left_join(
    HEK_Nterm_3_deg_ratio |> filter(Index %in% overlap_Nterm),
    by = c('Index', 'UniProt_Accession', 'Protein.Start', 'Gene', 'Entry.Name', 'timepoint')
  ) |> 
  mutate(
    deg_ratio_avg = (deg_ratio.Nt_1 + deg_ratio.Nt_2 + deg_ratio)/3
  ) |> 
  select(Index, UniProt_Accession, Protein.Start, Gene, Entry.Name, timepoint, deg_ratio_avg)

# overlap in HEK_Nt_1 and HKE_Nt_2
overlap_1_2_deg_ratio <- HEK_Nterm_1_deg_ratio |> 
  filter(Index %in% overlap_1_2) |> 
  left_join(
    HEK_Nterm_2_deg_ratio |> filter(Index %in% overlap_1_2),
    by = c('Index', 'UniProt_Accession', 'Protein.Start', 'Gene', 'Entry.Name', 'timepoint'),
    suffix = c('.Nt_1', '.Nt_2')
  ) |> 
  mutate(
    deg_ratio_avg = (deg_ratio.Nt_1 + deg_ratio.Nt_2)/2
  ) |> 
  select(Index, UniProt_Accession, Protein.Start, Gene, Entry.Name, timepoint, deg_ratio_avg)

# overlap in HEK_Nt_2 and HEK_Nt_3
overlap_2_3_deg_ratio <- HEK_Nterm_2_deg_ratio |> 
  filter(Index %in% overlap_2_3) |> 
  left_join(
    HEK_Nterm_3_deg_ratio |> filter(Index %in% overlap_2_3),
    by = c('Index', 'UniProt_Accession', 'Protein.Start', 'Gene', 'Entry.Name', 'timepoint'),
    suffix = c('.Nt_2', '.Nt_3')
  ) |> 
  mutate(
    deg_ratio_avg = (deg_ratio.Nt_2 + deg_ratio.Nt_3)/2
  ) |> 
  select(Index, UniProt_Accession, Protein.Start, Gene, Entry.Name, timepoint, deg_ratio_avg)

# overlap in HEK_Nt_3 and HEK_Nt_1
overlap_3_1_deg_ratio <- HEK_Nterm_3_deg_ratio |> 
  filter(Index %in% overlap_3_1) |> 
  left_join(
    HEK_Nterm_1_deg_ratio |> filter(Index %in% overlap_3_1),
    by = c('Index', 'UniProt_Accession', 'Protein.Start', 'Gene', 'Entry.Name', 'timepoint'),
    suffix = c('.Nt_3', '.Nt_1')
  ) |> 
  mutate(
    deg_ratio_avg = (deg_ratio.Nt_3 + deg_ratio.Nt_1)/2
  ) |> 
  select(Index, UniProt_Accession, Protein.Start, Gene, Entry.Name, timepoint, deg_ratio_avg)

# unique in HEK_Nt_1
unique_1_deg_ratio <- HEK_Nterm_1_deg_ratio |> 
  filter(Index %in% unique_1) |> 
  mutate(
    deg_ratio_avg = deg_ratio
  ) |> 
  select(Index, UniProt_Accession, Protein.Start, Gene, Entry.Name, timepoint, deg_ratio_avg)

# unique in HEK_Nt_2
unique_2_deg_ratio <- HEK_Nterm_2_deg_ratio |> 
  filter(Index %in% unique_2) |> 
  mutate(
    deg_ratio_avg = deg_ratio
  ) |> 
  select(Index, UniProt_Accession, Protein.Start, Gene, Entry.Name, timepoint, deg_ratio_avg)

# unique in HEK_Nt_3
unique_3_deg_ratio <- HEK_Nterm_3_deg_ratio |> 
  filter(Index %in% unique_3) |> 
  mutate(
    deg_ratio_avg = deg_ratio
  ) |> 
  select(Index, UniProt_Accession, Protein.Start, Gene, Entry.Name, timepoint, deg_ratio_avg)

# combine results from three replicates
HEK_Nterm_deg_ratio <- bind_rows(
  overlap_Nterm_deg_ratio,
  overlap_1_2_deg_ratio,
  overlap_2_3_deg_ratio,
  overlap_3_1_deg_ratio,
  unique_1_deg_ratio,
  unique_2_deg_ratio,
  unique_3_deg_ratio
)

write_csv(HEK_Nterm_deg_ratio, file = 'data_source/degradation_ratio/HEK_Nterm_deg_ratio.csv')

## ratio normalization by LaminB and T-complex proteoforms
HEK_Nterm_deg_ratio_adj <- HEK_Nterm_deg_ratio |> 
  left_join(
    ratio_norm_facs_Nterm, by = 'timepoint'
  ) |> 
  mutate(
    deg_ratio_adj = deg_ratio_avg/ratio_norm_fac
  )

write_csv(HEK_Nterm_deg_ratio_adj, file = 'data_source/degradation_ratio/HEK_Nterm_deg_ratio_adj.csv')

## curve fitting
# linear model
linear_model <- function(df) {
  tryCatch({
    # Fit the model
    model <- lm(
      ln_deg_ratio_adj ~ timepoint,
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
HEK_Nterm_curve_fitting_result <- HEK_Nterm_deg_ratio_adj |> 
  select(Index, UniProt_Accession, Protein.Start, Gene, Entry.Name, timepoint, deg_ratio_adj) |> 
  mutate(
    ln_deg_ratio_adj = log(deg_ratio_adj)
  ) |> 
  group_by(Index, UniProt_Accession, Protein.Start, Gene, Entry.Name) |> 
  group_modify(~ linear_model(.x)) |> 
  ungroup()

HEK_Nterm_linear_fitting <- HEK_Nterm_curve_fitting_result |> 
  filter(RSS < 0.05) |> 
  select(Index, UniProt_Accession, Protein.Start, Gene, Entry.Name, lnA:RSS) |> 
  pivot_longer(cols = lnA:RSS, names_to = 'parameters', values_to = 'values') |> 
  mutate(
    model = 'linear fitting'
  )

write_csv(HEK_Nterm_linear_fitting, file = 'data_source/curve_fitting/HEK_Nterm_linear_fitting.csv')

# # non-linear model fitting
# HEK_Nterm_curve_fitting_result_2 <- HEK_Nterm_curve_fitting_result_1 |> 
#   filter(is.na(A) | RSS >= 0.05) |> 
#   left_join(HEK_Nterm_deg_ratio, by = c('Index', 'UniProt_Accession', 'Protein.Start', 'Gene', 'Entry.Name')) |> 
#   select(Index, UniProt_Accession, Protein.Start, Gene, Entry.Name, timepoint, deg_ratio_avg) |> 
#   mutate(
#     ln_deg_ratio_avg = log(deg_ratio_avg)
#   ) |> 
#   group_by(Index, UniProt_Accession, Protein.Start, Gene, Entry.Name) |> 
#   group_modify(~ linear_model(.x)) |> 
#   ungroup()
# 
# HEK_Nterm_linear_fitting <- HEK_Nterm_curve_fitting_result_2 |> 
#   filter(RSS < 0.05) |> 
#   select(Index, UniProt_Accession, Protein.Start, Gene, Entry.Name, lnA:RSS) |> 
#   pivot_longer(cols = lnA:RSS, names_to = 'parameters', values_to = 'values') |> 
#   mutate(
#     model = 'linear fitting'
#   )
# 
# write_csv(HEK_Nterm_linear_fitting, file = 'data_source/curve_fitting/HEK_Nterm_linear_fitting.csv')
