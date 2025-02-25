# import packages
library(tidyverse)

# define color
color_1 <- '#4b1978'
color_2 <- '#18647c'
color_3 <- '#0164c9'
color_4 <- '#00b32d'
color_5 <- '#c73001'
color_6 <- '#ffb042'

# psm 
HEK_Nt_1_psm <- read_csv(
  'data_source/psm/HEK_Nt_1_psm.csv'
)
HEK_Nt_2_psm <- read_csv(
  'data_source/psm/HEK_Nt_2_psm.csv'
)
HEK_Nt_3_psm <- read_csv(
  'data_source/psm/HEK_Nt_3_psm.csv'
)

# TMT intensity normalized psm
HEK_Nt_1_psm_TMTi <- read_csv(
  'data_source/normalized_data/HEK_Nt_1_psm_TMTi.csv'
)
HEK_Nt_2_psm_TMTi <- read_csv(
  'data_source/normalized_data/HEK_Nt_2_psm_TMTi.csv'
)
HEK_Nt_3_psm_TMTi <- read_csv(
  'data_source/normalized_data/HEK_Nt_3_psm_TMTi.csv'
)

# Remove heavy labeled and not N-term modified psm
HEK_Nt_1_psm_TMTi_light_N_term <- read_csv(
  'data_source/normalized_data/HEK_Nt_1_psm_TMTi_light_N_term.csv'
)
HEK_Nt_2_psm_TMTi_light_N_term <- read_csv(
  'data_source/normalized_data/HEK_Nt_2_psm_TMTi_light_N_term.csv'
)
HEK_Nt_3_psm_TMTi_light_N_term <- read_csv(
  'data_source/normalized_data/HEK_Nt_3_psm_TMTi_light_N_term.csv'
)

# group by N-term Index
HEK_Nterm_1 <- read_csv(
  'data_source/normalized_data/HEK_Nterm_1.csv'
)
HEK_Nterm_2 <- read_csv(
  'data_source/normalized_data/HEK_Nterm_2.csv'
)
HEK_Nterm_3 <- read_csv(
  'data_source/normalized_data/HEK_Nterm_3.csv'
)

#internal reference scaling normalization
HEK_Nterm_1_irs <- read_csv(
  'data_source/normalized_data/HEK_Nterm_1_irs.csv'
)
HEK_Nterm_2_irs <- read_csv(
  'data_source/normalized_data/HEK_Nterm_2_irs.csv'
)
HEK_Nterm_3_irs <- read_csv(
  'data_source/normalized_data/HEK_Nterm_3_irs.csv'
)

# degradation ratio
HEK_Nterm_1_deg_ratio <- read_csv(
  'data_source/degradation_ratio/HEK_Nterm_1_deg_ratio.csv'
)
HEK_Nterm_2_deg_ratio <- read_csv(
  'data_source/degradation_ratio/HEK_Nterm_2_deg_ratio.csv'
)
HEK_Nterm_3_deg_ratio <- read_csv(
  'data_source/degradation_ratio/HEK_Nterm_3_deg_ratio.csv'
)
HEK_Nterm_deg_ratio <- read_csv(
  'data_source/degradation_ratio/HEK_Nterm_deg_ratio.csv'
)

# curve fitting
HEK_Nterm_curve_fitting_combined <- read_csv(
  'data_source/curve_fitting/HEK_Nterm_curve_fitting_combined.csv'
)

# half life for Nterm
HEK_Nterm_Kd_half_life_LaminB_Tcomplex <- read_csv(
  'data_source/Kd_half_life/HEK_Nterm_Kd_half_life_LaminB_Tcomplex.csv'
)

# half life for WP

# Nterm sequence
HEK_Nterm_Kd_half_life_sequence <- read_csv(
  'data_source/Nterm_sequence/HEK_Nterm_Kd_half_life_sequence.csv'
)

