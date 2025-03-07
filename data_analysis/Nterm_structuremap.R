# import packages
library(tidyverse)

# for executing python script in Rstudio
library(reticulate)

# use specific virtual env created by anaconda
use_condaenv(
  condaenv = '/opt/anaconda3/envs/structuremap',
  required = TRUE
)

# execute the python script for Nterm structure
source_python("data_analysis/Nterm_structuremap.py")

# combine Nterm AlphaFold information with Nterm half life result
Nterm_alphafold_half_life <- tibble(Nterm_alphafold_accessibility_smooth) |> 
  mutate(
    Index = str_c(protein_id, position, sep = '_')
  ) |> 
  select(
    Index, UniProt_Accession = protein_id, secondary_structure:IDR
  ) |> 
  filter(Index %in% HEK_Nterm_Kd_half_life_sequence$Index) |> 
  left_join(HEK_Nterm_Kd_half_life_sequence, by = c('Index', 'UniProt_Accession'))

# Wilcoxon rank-sum test
library(rstatix)

# secondary_structure
Nterm_alphafold_half_life |> 
  wilcox_test(half_life ~ secondary_structure) |> 
  filter(p < 0.05)

# structure_group
Nterm_alphafold_half_life |> 
  wilcox_test(half_life ~ structure_group) |> 
  filter(p < 0.05)

ks.test(
  Nterm_alphafold_half_life |> filter(structure_group == 'unstructured') |> pull(half_life),
  Nterm_alphafold_half_life |> filter(structure_group == 'BEND') |> pull(half_life)
)

# accessibility
Nterm_alphafold_half_life |> 
  wilcox_test(half_life ~ accessibility) |> 
  filter(p < 0.05)

ks.test(
  Nterm_alphafold_half_life |> filter(accessibility == 'high exposure') |> pull(half_life),
  Nterm_alphafold_half_life |> filter(accessibility == 'low exposure') |> pull(half_life)
)

# IDR
Nterm_alphafold_half_life |> 
  wilcox_test(half_life ~ IDR) |> 
  filter(p < 0.05)

ks.test(
  Nterm_alphafold_half_life |> filter(IDR == 'IDR') |> pull(half_life),
  Nterm_alphafold_half_life |> filter(IDR == 'Structural Region') |> pull(half_life)
)
