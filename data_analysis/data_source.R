#import packages
packages_names <- c("tidyverse", "readxl", "writexl")
lapply(packages_names, require, character.only = TRUE)

#psm 
HEK_Nt_1_psm <- read_csv(
  'data_source/psm/HEK_Nt_1_psm.csv'
)
HEK_Nt_2_psm <- read_csv(
  'data_source/psm/HEK_Nt_2_psm.csv'
)
HEK_Nt_3_psm <- read_csv(
  'data_source/psm/HEK_Nt_3_psm.csv'
)

#TMT intensity normalized psm
HEK_Nt_1_psm_TMTi <- read_csv(
  'data_source/normalized_data/HEK_Nt_1_psm_TMTi.csv'
)
HEK_Nt_2_psm_TMTi <- read_csv(
  'data_source/normalized_data/HEK_Nt_2_psm_TMTi.csv'
)
HEK_Nt_3_psm_TMTi <- read_csv(
  'data_source/normalized_data/HEK_Nt_3_psm_TMTi.csv'
)

#Remove heavy labeled and not N-term modified psm
HEK_Nt_1_psm_TMTi_light_N_term <- read_csv(
  'data_source/normalized_data/HEK_Nt_1_psm_TMTi_light_N_term.csv'
)
HEK_Nt_2_psm_TMTi_light_N_term <- read_csv(
  'data_source/normalized_data/HEK_Nt_2_psm_TMTi_light_N_term.csv'
)
HEK_Nt_3_psm_TMTi_light_N_term <- read_csv(
  'data_source/normalized_data/HEK_Nt_3_psm_TMTi_light_N_term.csv'
)

#group by N-term Index
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

#degradation ratio
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

