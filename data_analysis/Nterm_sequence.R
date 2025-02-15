#import packages
packages_names <- c("tidyverse", 'Biostrings')
lapply(packages_names, require, character.only = TRUE)

#import human fasta downloaded from UniProt (https://www.uniprot.org/)
human_fasta <- readAAStringSet(
  'data_source/uniprotkb_reviewed_true_AND_model_organ_2025_02_15.fasta'
)

