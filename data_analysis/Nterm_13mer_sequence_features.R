library(tidyverse)

# for executing python script in Rstudio
library(reticulate)

# use specific virtual env created by anaconda
use_condaenv(
  condaenv = '/opt/anaconda3/envs/Nterm_probe',
  required = TRUE
)

# execute the python script for calculating the sequence features
source_python("data_analysis/Nterm_sequence_features.py")

# import the result from localCIDER (data_analysis/Nterm_sequence_features.py)
Nterm_13mer_sequence_features <- read_csv(
  'data_source/Nterm_13mer_sequence_features/HEK_Nterm_Kd_half_life_sequence_features.csv'
)

## Wilcoxon rank-sum test
library(rstatix)

# hydropathy
Nterm_13mer_sequence_features |> 
  wilcox_test(hydropathy ~ category)

# FCR
Nterm_13mer_sequence_features |> 
  wilcox_test(FCR ~ category)

# NCPR
Nterm_13mer_sequence_features |> 
  wilcox_test(NCPR ~ category)

# isoelectric point
Nterm_13mer_sequence_features |> 
  wilcox_test(isoelectric_point ~ category)

# kappa
Nterm_13mer_sequence_features |> 
  wilcox_test(kappa ~ category)
