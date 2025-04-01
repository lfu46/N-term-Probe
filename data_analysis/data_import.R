# import packages
library(tidyverse)

### import degradation data from MSFragger
## Nterm
# HEK_Nt_1
HEK_Nt_1_psm <- read_tsv(
  'data_source/raw_data/HEK_Nt_1_psm.tsv',
  col_names = TRUE,
  name_repair = "universal"
) |> 
  filter(Is.Unique == TRUE) |>
  select(Peptide, Modified.Peptide, Observed.M.Z, Hyperscore, 
         Protein.Start, Protein.End, Is.Unique,
         Assigned.Modifications, UniProt_Accession = Protein.ID, Gene, Entry.Name,
         'deg_126_0h' = ..126_0h, 'deg_127_3h' = ..127_3h, 'deg_128_6h' = ..128_6h,
         'deg_129_9h' = ..129_9h, 'deg_130_12h' = ..130_12h, 'deg_131_24h' = ..131_24h) |> 
  filter(str_detect(Entry.Name, 'HUMAN')) |> 
  mutate(
    Index = str_c(UniProt_Accession, Protein.Start, sep = '_')
  ) |> 
  select(Index, everything())

write_csv(HEK_Nt_1_psm, file = 'data_source/psm/HEK_Nt_1_psm.csv')

# HEK_Nt_2
HEK_Nt_2_psm <- read_tsv(
  'data_source/raw_data/HEK_Nt_2_psm.tsv',
  col_names = TRUE,
  name_repair = "universal"
) |> 
  filter(Is.Unique == TRUE) |> 
  select(Peptide, Modified.Peptide, Observed.M.Z, Hyperscore, 
         Protein.Start, Protein.End, Is.Unique,
         Assigned.Modifications, UniProt_Accession = Protein.ID, Gene, Entry.Name,
         'deg_126_0h' = ..126_0h, 'deg_127_3h' = ..127_3h, 'deg_128_6h' = ..128_6h,
         'deg_129_9h' = ..129_9h, 'deg_130_12h' = ..130_12h, 'deg_131_24h' = ..131_24h) |> 
  filter(str_detect(Entry.Name, 'HUMAN')) |> 
  mutate(
    Index = str_c(UniProt_Accession, Protein.Start, sep = '_')
  ) |> 
  select(Index, everything())

write_csv(HEK_Nt_2_psm, file = 'data_source/psm/HEK_Nt_2_psm.csv')

# HEK_Nt_3
HEK_Nt_3_psm <- read_tsv(
  'data_source/raw_data/HEK_Nt_3_psm.tsv',
  col_names = TRUE,
  name_repair = "universal"
) |> 
  filter(Is.Unique == TRUE) |> 
  select(Peptide, Modified.Peptide, Observed.M.Z, Hyperscore, 
         Protein.Start, Protein.End, Is.Unique,
         Assigned.Modifications, UniProt_Accession = Protein.ID, Gene, Entry.Name,
         'deg_126_0h' = ..126_0h, 'deg_127_3h' = ..127_3h, 'deg_128_6h' = ..128_6h,
         'deg_129_9h' = ..129_9h, 'deg_130_12h' = ..130_12h, 'deg_131_24h' = ..131_24h) |> 
  filter(str_detect(Entry.Name, 'HUMAN')) |> 
  mutate(
    Index = str_c(UniProt_Accession, Protein.Start, sep = '_')
  ) |> 
  select(Index, everything())

write_csv(HEK_Nt_3_psm, file = 'data_source/psm/HEK_Nt_3_psm.csv')

## Whole Proteome
# HEK_WP_1
HEK_WP_1_psm <- read_tsv(
  'data_source/raw_data/HEK_WP_1_psm.tsv',
  col_names = TRUE,
  name_repair = "universal"
) |> 
  filter(Is.Unique == TRUE) |> 
  select(Peptide, Modified.Peptide, Observed.M.Z, Hyperscore, 
         Protein.Start, Protein.End, Is.Unique,
         Assigned.Modifications, UniProt_Accession = Protein.ID, Gene, Entry.Name,
         'deg_126_0h' = ..126_0h, 'deg_127_3h' = ..127_3h, 'deg_128_6h' = ..128_6h,
         'deg_129_9h' = ..129_9h, 'deg_130_12h' = ..130_12h, 'deg_131_24h' = ..131_24h) |> 
  filter(str_detect(Entry.Name, 'HUMAN')) |> 
  mutate(
    Index = str_c(UniProt_Accession, Protein.Start, sep = '_')
  ) |> 
  select(Index, everything())

write_csv(HEK_WP_1_psm, file = 'data_source/psm/HEK_WP_1_psm.csv')

# HEK_WP_2
HEK_WP_2_psm <- read_tsv(
  'data_source/raw_data/HEK_WP_2_psm.tsv',
  col_names = TRUE,
  name_repair = "universal"
) |> 
  filter(Is.Unique == TRUE) |> 
  select(Peptide, Modified.Peptide, Observed.M.Z, Hyperscore, 
         Protein.Start, Protein.End, Is.Unique,
         Assigned.Modifications, UniProt_Accession = Protein.ID, Gene, Entry.Name,
         'deg_126_0h' = ..126_0h, 'deg_127_3h' = ..127_3h, 'deg_128_6h' = ..128_6h,
         'deg_129_9h' = ..129_9h, 'deg_130_12h' = ..130_12h, 'deg_131_24h' = ..131_24h) |> 
  filter(str_detect(Entry.Name, 'HUMAN')) |> 
  mutate(
    Index = str_c(UniProt_Accession, Protein.Start, sep = '_')
  ) |> 
  select(Index, everything())

write_csv(HEK_WP_2_psm, file = 'data_source/psm/HEK_WP_2_psm.csv')

# HEK_WP_3
HEK_WP_3_psm <- read_tsv(
  'data_source/raw_data/HEK_WP_3_psm.tsv',
  col_names = TRUE,
  name_repair = "universal"
) |> 
  filter(Is.Unique == TRUE) |> 
  select(Peptide, Modified.Peptide, Observed.M.Z, Hyperscore, 
         Protein.Start, Protein.End, Is.Unique,
         Assigned.Modifications, UniProt_Accession = Protein.ID, Gene, Entry.Name,
         'deg_126_0h' = ..126_0h, 'deg_127_3h' = ..127_3h, 'deg_128_6h' = ..128_6h,
         'deg_129_9h' = ..129_9h, 'deg_130_12h' = ..130_12h, 'deg_131_24h' = ..131_24h) |> 
  filter(str_detect(Entry.Name, 'HUMAN')) |> 
  mutate(
    Index = str_c(UniProt_Accession, Protein.Start, sep = '_')
  ) |> 
  select(Index, everything())

write_csv(HEK_WP_3_psm, file = 'data_source/psm/HEK_WP_3_psm.csv')
