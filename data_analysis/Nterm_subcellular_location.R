# import packages
library(tidyverse)

# import sub-cellular location data from Human Protein Atlas (https://www.proteinatlas.org/)
hpa_subcellular_location <- read_tsv(
  'data_source/HPA_subcellular_location/subcellular_location.tsv',
  col_names = TRUE,
  name_repair = 'universal'
) |> 
  filter(Reliability != 'Uncertain') |> 
  separate_rows(Main.location, sep = ';') |> 
  dplyr::select(Gene.name, Main.location) |> 
  group_by(Gene.name) |> 
  dplyr::filter(n() == 1) |> 
  ungroup()

# Nterm sub-cellular location half life
HEK_Nterm_Kd_half_life_subcellular <- HEK_Nterm_Kd_half_life_LaminB_Tcomplex |> 
  left_join(hpa_subcellular_location, by = join_by('Gene' == 'Gene.name')) |> 
  filter(!is.na(Main.location))

# sub-cellular location median half life
library(rstatix)

Nterm_subcellular_median_half_life <- HEK_Nterm_Kd_half_life_subcellular |> 
  group_by(Main.location) |> 
  get_summary_stats(half_life, type = 'median') |> 
  arrange(desc(median))

write_csv(Nterm_subcellular_median_half_life, file = 'data_source/HPA_subcellular_location/Nterm_subcellular_median_half_life.csv')

# Wilcoxon rank-sum test
HEK_Nterm_Kd_half_life_subcellular |> 
  wilcox_test(half_life ~ Main.location) |> 
  filter(p.adj.signif != 'ns')


