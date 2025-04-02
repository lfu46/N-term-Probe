# import packages
library(tidyverse)
library(Biostrings)

# import human fasta downloaded from UniProt (https://www.uniprot.org/)
human_fasta <- readAAStringSet(
  'data_source/fasta_file/uniprotkb_reviewed_true_AND_model_organ_2025_02_15.fasta'
)

# build tibble using human fasta
human_fasta_tibble <- tibble(
  Name = names(human_fasta),
  Sequence = as.character(human_fasta),
  Length = width(human_fasta)
) |> 
  mutate(
    Name = sub(' .*', '', Name),
    Full_Protein_Length = as.numeric(Length)
  ) |> 
  separate(Name, into = c('sp', 'UniProt_Accession', 'name'), sep = '\\|') |> 
  select(UniProt_Accession, Sequence, Full_Protein_Length)

# generate Nterm sequence
HEK_Nterm_Kd_half_life_sequence <- HEK_Nterm_Kd_half_life_LaminB_Tcomplex |> 
  left_join(human_fasta_tibble, by = 'UniProt_Accession') |> 
  mutate(
    Nterm_sequence = substr(Sequence, start = Protein.Start, stop = Full_Protein_Length),
    Nterm_13mer = substr(Sequence, start = Protein.Start, stop = Protein.Start + 12),
    Nterm_terminus = substr(Sequence, start = Protein.Start, stop = Protein.Start),
    Nterm_15mer = substr(Sequence, start = Protein.Start - 7, stop = Protein.Start + 7)
  )

write_csv(HEK_Nterm_Kd_half_life_sequence, file = 'data_source/Nterm_sequence/HEK_Nterm_Kd_half_life_sequence.csv')
