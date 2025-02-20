# import packages
packages_names <- c('tidyverse', 'rstatix', 'readxl')
lapply(packages_names, require, character.only = TRUE)

### import stress granule, processing body and cajal body database
## stress granule 
## (https://rnagranuledb.lunenfeld.ca/, Species == 'Homo sapiens', Type = 'Stress Granules', Gold standard (only tier 1))
stress_granule_list <- read_csv(
  'data_source/cytoplasmic_body/user_gold_report_stress_granule.csv',
  skip = 2,
  col_names = TRUE,
  name_repair = 'universal'
) |> 
  select(Gene.Name) |> 
  distinct()

## processing body 
## (https://rnagranuledb.lunenfeld.ca/, Species == 'Homo sapiens', Type = 'P-bodies', Gold standard (only tier 1))
processing_body_list <- read_csv(
  'data_source/cytoplasmic_body/user_gold_report_P_body.csv',
  skip = 2,
  col_names = TRUE,
  name_repair = 'universal'
) |> 
  select(Gene.Name) |> 
  distinct()

## cajal body 
## (https://www.sciencedirect.com/science/article/pii/S2211124724010854?via%3Dihub, Table S4 and Table S6)
# Table S4
cajal_body_list_1 <- read_xlsx(
  'data_source/cytoplasmic_body/mmc6.xlsx',
  sheet = 'HEK_insitu_FC>10',
  col_names = TRUE,
  .name_repair = 'universal'
) |> 
  select(UniProt_Accession = Protein.IDs) |> 
  separate_rows(UniProt_Accession, sep = ';') |> 
  distinct()

# Table S6
cajal_body_list_2 <- read_xlsx(
  'data_source/cytoplasmic_body/mmc7.xlsx',
  sheet = 'HEK_turbo-Coilin_FC>10',
  col_names = TRUE,
  .name_repair = 'universal'
) |> 
  select(UniProt_Accession = Protein.IDs) |> 
  separate_rows(UniProt_Accession, sep = ';') |> 
  distinct()

# combine cajal body list 1 and list 2
cajal_body_list_combined <- bind_rows(
  cajal_body_list_1,
  cajal_body_list_2
) |> 
  distinct()

## combine cytoplasmic body protein list with Nterm half-life data
Nterm_cytoplasmic_body_half_life <- bind_rows(
  HEK_Nterm_Kd_half_life_LaminB_Tcomplex |> 
    filter(Gene %in% stress_granule_list$Gene.Name) |> 
    mutate(
      category = 'stress granule'
    ),
  
  HEK_Nterm_Kd_half_life_LaminB_Tcomplex |> 
    filter(Gene %in% processing_body_list$Gene.Name) |> 
    mutate(
      category = 'processing body'
    ),
  
  HEK_Nterm_Kd_half_life_LaminB_Tcomplex |> 
    filter(UniProt_Accession %in% cajal_body_list_combined$UniProt_Accession) |> 
    mutate(
      category = 'cajal body'
    )
)

# Wilcoxon rank-sum test
Nterm_cytoplasmic_body |> 
  wilcox_test(half_life ~ category)

