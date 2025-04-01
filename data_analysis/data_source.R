# import packages
library(tidyverse)

# define color
color_1 <- '#4b1978'
color_2 <- '#18647c'
color_3 <- '#0164c9'
color_4 <- '#00b32d'
color_5 <- '#ffb042'

### psm 
# Nterm
HEK_Nt_1_psm <- read_csv(
  'data_source/psm/HEK_Nt_1_psm.csv'
)
HEK_Nt_2_psm <- read_csv(
  'data_source/psm/HEK_Nt_2_psm.csv'
)
HEK_Nt_3_psm <- read_csv(
  'data_source/psm/HEK_Nt_3_psm.csv'
)

# Whole Proteome
HEK_WP_1_psm <- read_csv(
  'data_source/psm/HEK_WP_1_psm.csv'
)
HEK_WP_2_psm <- read_csv(
  'data_source/psm/HEK_WP_2_psm.csv'
)
HEK_WP_3_psm <- read_csv(
  'data_source/psm/HEK_WP_3_psm.csv'
)

### TMT intensity normalized psm
# Nterm
HEK_Nt_1_psm_TMTi <- read_csv(
  'data_source/normalized_data/HEK_Nt_1_psm_TMTi.csv'
)
HEK_Nt_2_psm_TMTi <- read_csv(
  'data_source/normalized_data/HEK_Nt_2_psm_TMTi.csv'
)
HEK_Nt_3_psm_TMTi <- read_csv(
  'data_source/normalized_data/HEK_Nt_3_psm_TMTi.csv'
)

# Whole Proteome
HEK_WP_1_psm_TMTi <- read_csv(
  'data_source/normalized_data/HEK_WP_1_psm_TMTi.csv'
)
HEK_WP_2_psm_TMTi <- read_csv(
  'data_source/normalized_data/HEK_WP_2_psm_TMTi.csv'
)
HEK_WP_3_psm_TMTi <- read_csv(
  'data_source/normalized_data/HEK_WP_3_psm_TMTi.csv'
)

# Nterm, remove heavy labeled and not N-term modified psm
HEK_Nt_1_psm_TMTi_light_N_term <- read_csv(
  'data_source/normalized_data/HEK_Nt_1_psm_TMTi_light_N_term.csv'
)
HEK_Nt_2_psm_TMTi_light_N_term <- read_csv(
  'data_source/normalized_data/HEK_Nt_2_psm_TMTi_light_N_term.csv'
)
HEK_Nt_3_psm_TMTi_light_N_term <- read_csv(
  'data_source/normalized_data/HEK_Nt_3_psm_TMTi_light_N_term.csv'
)
# Whole Proteome, remove heavy labeled and not N-term modified psm
HEK_WP_1_psm_TMTi_light <- read_csv(
  'data_source/normalized_data/HEK_WP_1_psm_TMTi_light.csv'
)
HEK_WP_2_psm_TMTi_light <- read_csv(
  'data_source/normalized_data/HEK_WP_2_psm_TMTi_light.csv'
)
HEK_WP_3_psm_TMTi_light <- read_csv(
  'data_source/normalized_data/HEK_WP_3_psm_TMTi_light.csv'
)

# Nterm, group by N-term Index
HEK_Nterm_1 <- read_csv(
  'data_source/normalized_data/HEK_Nterm_1.csv'
)
HEK_Nterm_2 <- read_csv(
  'data_source/normalized_data/HEK_Nterm_2.csv'
)
HEK_Nterm_3 <- read_csv(
  'data_source/normalized_data/HEK_Nterm_3.csv'
)

# Whole Proteome, group by UniProt_Accession
HEK_Whole_Proteome_1 <- read_csv(
  'data_source/normalized_data/HEK_Whole_Proteome_1.csv'
)
HEK_Whole_Proteome_2 <- read_csv(
  'data_source/normalized_data/HEK_Whole_Proteome_2.csv'
)
HEK_Whole_Proteome_3 <- read_csv(
  'data_source/normalized_data/HEK_Whole_Proteome_3.csv'
)

### internal reference scaling normalization
# Nterm
HEK_Nterm_1_irs <- read_csv(
  'data_source/normalized_data/HEK_Nterm_1_irs.csv'
)
HEK_Nterm_2_irs <- read_csv(
  'data_source/normalized_data/HEK_Nterm_2_irs.csv'
)
HEK_Nterm_3_irs <- read_csv(
  'data_source/normalized_data/HEK_Nterm_3_irs.csv'
)
# Whole Proteome
HEK_Whole_Proteome_1_irs <- read_csv(
  'data_source/normalized_data/HEK_Whole_Proteome_1_irs.csv'
)
HEK_Whole_Proteome_2_irs <- read_csv(
  'data_source/normalized_data/HEK_Whole_Proteome_2_irs.csv'
)
HEK_Whole_Proteome_3_irs <- read_csv(
  'data_source/normalized_data/HEK_Whole_Proteome_3_irs.csv'
)

### degradation ratio
# Nterm
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
HEK_Nterm_deg_ratio_adj <- read_csv(
  'data_source/degradation_ratio/HEK_Nterm_deg_ratio_adj.csv'
)
# Whole Proteome
HEK_WP_1_deg_ratio <- read_csv(
  'data_source/degradation_ratio/HEK_WP_1_deg_ratio.csv'
)
HEK_WP_2_deg_ratio <- read_csv(
  'data_source/degradation_ratio/HEK_WP_2_deg_ratio.csv'
)
HEK_WP_3_deg_ratio <- read_csv(
  'data_source/degradation_ratio/HEK_WP_3_deg_ratio.csv'
)
HEK_WP_deg_ratio <- read_csv(
  'data_source/degradation_ratio/HEK_WP_deg_ratio.csv'
)
HEK_WP_deg_ratio_adj <- read_csv(
  'data_source/degradation_ratio/HEK_WP_deg_ratio_adj.csv'
)

### cell doubling kinetics
K_cd_1 <- 0.0338
K_cd_2 <- 0.0397
K_cd_3 <- 0.0443

### Lamin B and T-complex normalization
## Nterm
# curve fitting
HEK_Nterm_LaminB_Tcomplex_result <- read_csv(
  'data_source/LaminB_Tcomplex_normalization/HEK_Nterm_LaminB_Tcomplex_result.csv'
)
# half life
HEK_Nterm_LaminB_Tcomplex_half_life <- read_csv(
  'data_source/LaminB_Tcomplex_normalization/HEK_Nterm_LaminB_Tcomplex.csv'
)
# adjusted ratio
HEK_Nterm_LaminB_Tcomplex_adjusted_ratio <- read_csv(
  'data_source/LaminB_Tcomplex_normalization/HEK_Nterm_LaminB_Tcomplex_adjusted_ratio.csv'
)
# normalizatio factor
ratio_norm_facs_Nterm <- read_csv(
  'data_source/LaminB_Tcomplex_normalization/ratio_norm_facs_Nterm.csv'
)

## Whole Proteome
# curve fitting
HEK_WP_LaminB_Tcomplex_result <- read_csv(
  'data_source/LaminB_Tcomplex_normalization/HEK_WP_LaminB_Tcomplex_result.csv'
)
# half life
HEK_WP_LaminB_Tcomplex_half_life <- read_csv(
  'data_source/LaminB_Tcomplex_normalization/HEK_WP_LaminB_Tcomplex_half_life.csv'
)
# adjusted ratio
HEK_WP_LaminB_Tcomplex_adjusted_ratio <- read_csv(
  'data_source/LaminB_Tcomplex_normalization/HEK_WP_LaminB_Tcomplex_adjusted_ratio.csv'
)
# normalizatio factor
ratio_norm_facs_WP <- read_csv(
  'data_source/LaminB_Tcomplex_normalization/ratio_norm_facs_WP.csv'
)

### half-life calculation
# half-life for Nterm
HEK_Nterm_Kd_half_life_LaminB_Tcomplex <- read_csv(
  'data_source/Kd_half_life/HEK_Nterm_Kd_half_life_LaminB_Tcomplex.csv'
)
# half-life for WP
HEK_WP_Kd_half_life_LaminB_Tcomplex <- read_csv(
  'data_source/Kd_half_life/HEK_WP_Kd_half_life_LaminB_Tcomplex.csv'
)

### Nterm sequence
HEK_Nterm_Kd_half_life_sequence <- read_csv(
  'data_source/Nterm_sequence/HEK_Nterm_Kd_half_life_sequence.csv'
)

### Nterm ELM N-degron
HEK_Nterm_ELM_N_degron <- read_csv(
  'data_source/ELM_degron/HEK_Nterm_ELM_N_degron.csv'
)

### half-life comparison between Nterm and Whole Proteome
HEK_Nterm_WP_delta_half_life <- read_csv(
  'data_source/Nterm_WP_comparison/HEK_Nterm_WP_delta_half_life.csv'
)
