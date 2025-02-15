#import packages
packages_names <- c("tidyverse", "rstatix")
lapply(packages_names, require, character.only = TRUE)

#generate data frame for TopFinder (https://topfind.clip.msl.ubc.ca/topfinder) upload
topfinder_id <- HEK_Nterm_Kd_half_life_sequence |> 
  mutate(
    topfinder_id = paste(UniProt_Accession, Nterm_sequence, sep = ' ')
  ) |> 
  select(topfinder_id)

write_delim(topfinder_id, 
            file = 'data_source/Nterm_topfinder/topfinder_id.txt',
            col_names = FALSE,
            quote = 'none')
