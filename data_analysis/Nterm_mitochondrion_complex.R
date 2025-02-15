#import packages
packages_names <- c("tidyverse", "rstatix")
lapply(packages_names, require, character.only = TRUE)

#get annotation from CORUM (https://mips.helmholtz-muenchen.de/corum/)
complex_annotation_id <- read_delim(
  'data_source/Enzyme_Motif_Domain_Complex_analysis/corum_uniprotCorumMapping.txt',
  col_names = TRUE
) |> 
  mutate(
    CORUM_id = paste('CORUM_id', corum_id, sep = '_')
  ) |> 
  dplyr::select(UniProt_Accession = UniProtKB_accession_number, CORUM_id) |> 
  filter(UniProt_Accession %in% HEK_Nterm_Kd_half_life$UniProt_Accession)
colnames(complex_annotation_id) <- c('TERM', 'GENE')

complex_name_human <- read_delim(
  'data_source/Enzyme_Motif_Domain_Complex_analysis/corum_humanComplexes.txt',
  col_names = TRUE
) |> 
  mutate(
    CORUM_id = paste('CORUM_id', complex_id, sep = '_')
  ) |> 
  select(CORUM_id, complex_name)

complex_annotation_name <- complex_annotation_id |> 
  left_join(complex_name_human, by = join_by('GENE' == 'CORUM_id'))

#filter out mitochondrial complexes
mitochondrial_complex <- complex_annotation_name |> 
  filter(str_detect(complex_name, 'mitochondrial'))

#Nterm mitochondrial complexes half life
Nterm_mito_half_life <- mitochondrial_complex |> 
  left_join(HEK_Nterm_Kd_half_life, by = join_by('TERM' == 'UniProt_Accession'), relationship = 'many-to-many')

#get median value of Nterm mitochondiral complexes half life
Nterm_mito_half_life_median <- Nterm_mito_half_life |> 
  group_by(complex_name) |> 
  get_summary_stats(half_life, type = 'median') |> 
  arrange(median)

write_csv(Nterm_mito_half_life_median, file = 'data_source/mitochondrial_complexes/Nterm_mito_half_life_median.csv')

#Wilcoxon rank-sum test
Nterm_mito_half_life |> 
  filter(
    str_detect(complex_name, 'MICOS|ribosom|TIM23|Respiratory')
  ) |> 
  mutate(
    complex_name = case_when(
      str_detect(complex_name, 'MICOS') ~ 'MICOS complex related',
      str_detect(complex_name, 'ribosom') ~ 'Mitochondrial ribosome related',
      str_detect(complex_name, 'TIM23') ~ 'TIM23 complex related',
      str_detect(complex_name, 'Respiratory chain complex') ~ 'Respiratory chain complex I'
    )
  ) |> 
  wilcox_test(half_life ~ complex_name) |> 
  filter(p < 0.05)
