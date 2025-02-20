# import packages
packages_names <- c("tidyverse", 'readxl', 'rstatix')
lapply(packages_names, require, character.only = TRUE)

# import result from literature 
# (https://www.sciencedirect.com/science/article/pii/S1535947620350209?via%3Dihub, Supplemental Table S1)
cell_cycle_marker <- read_xlsx(
  'data_source/cell_cycle/158182_1_supp_471831_q5c4ws.xlsx',
  sheet = 'Proteomics and comparisons',
  skip = 4, 
  col_names = TRUE,
  .name_repair = 'universal'
) |> 
  select(UniProt_Accession = Unique.Uniprot.ID, CC_HighestLogFC, CC_Group...10) |> 
  filter(CC_HighestLogFC > log2(1.5)) |> 
  mutate(
    cell_cycle_phase = case_when(
      CC_Group...10 == 1 ~ 'G1',
      CC_Group...10 == 2 ~ 'G1+S',
      CC_Group...10 == 3 ~ 'S',
      CC_Group...10 == 4 ~ 'S+G2/M',
      CC_Group...10 == 5 ~ 'G2/M',
      CC_Group...10 == 6 ~ 'G2/M+G1'
    )
  ) |> 
  select(UniProt_Accession, cell_cycle_phase)

# combine cell cycle marker result with Nterm half-life result
Nterm_cell_cycle_marker_half_life <- HEK_Nterm_Kd_half_life_LaminB_Tcomplex |> 
  left_join(cell_cycle_marker, by = 'UniProt_Accession') |> 
  filter(!is.na(cell_cycle_phase))

write_csv(Nterm_cell_cycle_marker_half_life, file = 'data_source/cell_cycle/Nterm_cell_cycle_marker_half_life.csv')

