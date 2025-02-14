#import packages
packages_names <- c("tidyverse", "readxl", "writexl", "rstatix")
lapply(packages_names, require, character.only = TRUE)

##TMT channel intensity normalization
#HEK_Nt_1
target_mean_deg_TMTi_Nt_1 <- mean(colSums(HEK_Nt_1_psm |> select(starts_with('deg'))))
norm_facs_Nt_1 <- target_mean_deg_TMTi_Nt_1/colSums(HEK_Nt_1_psm |> select(starts_with('deg')))
deg_TMTi_Nt_1 <- tibble(sweep(HEK_Nt_1_psm |> select(starts_with('deg')), 2, norm_facs_Nt_1, FUN = '*'))
colnames(deg_TMTi_Nt_1) <- c(
  'deg_126_0h_TMTi', 'deg_127_3h_TMTi', 'deg_128_6h_TMTi', 'deg_129_9h_TMTi', 'deg_130_12h_TMTi', 'deg_131_24h_TMTi'
)

HEK_Nt_1_psm_TMTi <- bind_cols(HEK_Nt_1_psm, deg_TMTi_Nt_1)
write_csv(HEK_Nt_1_psm_TMTi, file = 'data_source/normalized_data/HEK_Nt_1_psm_TMTi.csv')

#HEK_Nt_2
target_mean_deg_TMTi_Nt_2 <- mean(colSums(HEK_Nt_2_psm |> select(starts_with('deg'))))
norm_facs_Nt_2 <- target_mean_deg_TMTi_Nt_2/colSums(HEK_Nt_2_psm |> select(starts_with('deg')))
deg_TMTi_Nt_2 <- tibble(sweep(HEK_Nt_2_psm |> select(starts_with('deg')), 2, norm_facs_Nt_2, FUN = '*'))
colnames(deg_TMTi_Nt_2) <- c(
  'deg_126_0h_TMTi', 'deg_127_3h_TMTi', 'deg_128_6h_TMTi', 'deg_129_9h_TMTi', 'deg_130_12h_TMTi', 'deg_131_24h_TMTi'
)

HEK_Nt_2_psm_TMTi <- bind_cols(HEK_Nt_2_psm, deg_TMTi_Nt_2)
write_csv(HEK_Nt_2_psm_TMTi, file = 'data_source/normalized_data/HEK_Nt_2_psm_TMTi.csv')

#HEK_Nt_3
target_mean_deg_TMTi_Nt_3 <- mean(colSums(HEK_Nt_3_psm |> select(starts_with('deg'))))
norm_facs_Nt_3 <- target_mean_deg_TMTi_Nt_3/colSums(HEK_Nt_3_psm |> select(starts_with('deg')))
deg_TMTi_Nt_3 <- tibble(sweep(HEK_Nt_3_psm |> select(starts_with('deg')), 2, norm_facs_Nt_3, FUN = '*'))
colnames(deg_TMTi_Nt_3) <- c(
  'deg_126_0h_TMTi', 'deg_127_3h_TMTi', 'deg_128_6h_TMTi', 'deg_129_9h_TMTi', 'deg_130_12h_TMTi', 'deg_131_24h_TMTi'
)

HEK_Nt_3_psm_TMTi <- bind_cols(HEK_Nt_3_psm, deg_TMTi_Nt_3)
write_csv(HEK_Nt_3_psm_TMTi, file = 'data_source/normalized_data/HEK_Nt_3_psm_TMTi.csv')

##Remove all the heavy labeled and not N-term modified psm
#HEK_Nt_1
HEK_Nt_1_psm_TMTi_light_N_term <- HEK_Nt_1_psm_TMTi |> 
  filter(str_detect(Assigned.Modifications, "N-term")) |> 
  filter(!str_detect(Assigned.Modifications, "K\\(8.0142\\)")) |> 
  filter(!str_detect(Assigned.Modifications, "K\\(237.1771\\)")) |> 
  filter(!str_detect(Assigned.Modifications, "R\\(10.0082\\)"))

write_csv(HEK_Nt_1_psm_TMTi_light_N_term, file = 'data_source/normalized_data/HEK_Nt_1_psm_TMTi_light_N_term.csv')

#HEK_Nt_2
HEK_Nt_2_psm_TMTi_light_N_term <- HEK_Nt_2_psm_TMTi |> 
  filter(str_detect(Assigned.Modifications, "N-term")) |> 
  filter(!str_detect(Assigned.Modifications, "K\\(8.0142\\)")) |> 
  filter(!str_detect(Assigned.Modifications, "K\\(237.1771\\)")) |> 
  filter(!str_detect(Assigned.Modifications, "R\\(10.0082\\)"))

write_csv(HEK_Nt_2_psm_TMTi_light_N_term, file = 'data_source/normalized_data/HEK_Nt_2_psm_TMTi_light_N_term.csv')

#HEK_Nt_3
HEK_Nt_3_psm_TMTi_light_N_term <- HEK_Nt_3_psm_TMTi |> 
  filter(str_detect(Assigned.Modifications, "N-term")) |> 
  filter(!str_detect(Assigned.Modifications, "K\\(8.0142\\)")) |> 
  filter(!str_detect(Assigned.Modifications, "K\\(237.1771\\)")) |> 
  filter(!str_detect(Assigned.Modifications, "R\\(10.0082\\)"))

write_csv(HEK_Nt_3_psm_TMTi_light_N_term, file = 'data_source/normalized_data/HEK_Nt_3_psm_TMTi_light_N_term.csv')

#group by N-term index
#HEK_Nt_1
HEK_Nterm_1 <- HEK_Nt_1_psm_TMTi_light_N_term |> 
  group_by(Index, UniProt_Accession, Protein.Start, Gene, Entry.Name) |> 
  summarise(
    deg_126_0h_TMTi = sum(deg_126_0h_TMTi),
    deg_127_3h_TMTi = sum(deg_127_3h_TMTi),
    deg_128_6h_TMTi = sum(deg_128_6h_TMTi),
    deg_129_9h_TMTi = sum(deg_129_9h_TMTi),
    deg_130_12h_TMTi = sum(deg_130_12h_TMTi),
    deg_131_24h_TMTi = sum(deg_131_24h_TMTi)
  ) |> 
  ungroup()

write_csv(HEK_Nterm_1, file = 'data_source/normalized_data/HEK_Nterm_1.csv')

#HEK_Nt_2
HEK_Nterm_2 <- HEK_Nt_2_psm_TMTi_light_N_term |> 
  group_by(Index, UniProt_Accession, Protein.Start, Gene, Entry.Name) |> 
  summarise(
    deg_126_0h_TMTi = sum(deg_126_0h_TMTi),
    deg_127_3h_TMTi = sum(deg_127_3h_TMTi),
    deg_128_6h_TMTi = sum(deg_128_6h_TMTi),
    deg_129_9h_TMTi = sum(deg_129_9h_TMTi),
    deg_130_12h_TMTi = sum(deg_130_12h_TMTi),
    deg_131_24h_TMTi = sum(deg_131_24h_TMTi)
  ) |> 
  ungroup()

write_csv(HEK_Nterm_2, file = 'data_source/normalized_data/HEK_Nterm_2.csv')

#HEK_Nt_3
HEK_Nterm_3 <- HEK_Nt_3_psm_TMTi_light_N_term |> 
  group_by(Index, UniProt_Accession, Protein.Start, Gene, Entry.Name) |> 
  summarise(
    deg_126_0h_TMTi = sum(deg_126_0h_TMTi),
    deg_127_3h_TMTi = sum(deg_127_3h_TMTi),
    deg_128_6h_TMTi = sum(deg_128_6h_TMTi),
    deg_129_9h_TMTi = sum(deg_129_9h_TMTi),
    deg_130_12h_TMTi = sum(deg_130_12h_TMTi),
    deg_131_24h_TMTi = sum(deg_131_24h_TMTi)
  ) |> 
  ungroup()

write_csv(HEK_Nterm_3, file = 'data_source/normalized_data/HEK_Nterm_3.csv')

##internal reference scaling normalization
#overlap of quantified N-term in three replicates
overlap_Nterm <- Reduce(intersect, list(
  HEK_Nterm_1$Index,
  HEK_Nterm_2$Index,
  HEK_Nterm_3$Index
))

HEK_Nterm_1_overlap <- HEK_Nterm_1 |> filter(Index %in% overlap_Nterm)
HEK_Nterm_2_overlap <- HEK_Nterm_2 |> filter(Index %in% overlap_Nterm)
HEK_Nterm_3_overlap <- HEK_Nterm_3 |> filter(Index %in% overlap_Nterm)

#calculate irs normalization factors
irs_factor <- tribble(
  ~ timepoint, ~ HEK_Nt_1, ~ HEK_Nt_2, ~ HEK_Nt_3,
  '126_0h', sum(HEK_Nterm_1_overlap$deg_126_0h_TMTi), sum(HEK_Nterm_2_overlap$deg_126_0h_TMTi), sum(HEK_Nterm_3_overlap$deg_126_0h_TMTi),
  '127_3h', sum(HEK_Nterm_1_overlap$deg_127_3h_TMTi), sum(HEK_Nterm_2_overlap$deg_127_3h_TMTi), sum(HEK_Nterm_3_overlap$deg_127_3h_TMTi),
  '128_6h', sum(HEK_Nterm_1_overlap$deg_128_6h_TMTi), sum(HEK_Nterm_2_overlap$deg_128_6h_TMTi), sum(HEK_Nterm_3_overlap$deg_128_6h_TMTi),
  '129_9h', sum(HEK_Nterm_1_overlap$deg_129_9h_TMTi), sum(HEK_Nterm_2_overlap$deg_129_9h_TMTi), sum(HEK_Nterm_3_overlap$deg_129_9h_TMTi),
  '130_12h', sum(HEK_Nterm_1_overlap$deg_130_12h_TMTi), sum(HEK_Nterm_2_overlap$deg_130_12h_TMTi), sum(HEK_Nterm_3_overlap$deg_130_12h_TMTi),
  '131_24h', sum(HEK_Nterm_1_overlap$deg_131_24h_TMTi), sum(HEK_Nterm_2_overlap$deg_131_24h_TMTi), sum(HEK_Nterm_3_overlap$deg_131_24h_TMTi)
) |> 
  mutate(
    average = (HEK_Nt_1 + HEK_Nt_2 + HEK_Nt_3)/3,
    irs_factor_Nt_1 = average/HEK_Nt_1,
    irs_factor_Nt_2 = average/HEK_Nt_2, 
    irs_factor_Nt_3 = average/HEK_Nt_3
  )

#apply irs normalization factors to data
#HEK_Nt_1
HEK_Nterm_1_irs <- HEK_Nterm_1 |> 
  mutate(
    deg_126_0h_TMTi_irs = deg_126_0h_TMTi * irs_factor |> filter(timepoint == "126_0h") |> pull(irs_factor_Nt_1),
    deg_127_3h_TMTi_irs = deg_127_3h_TMTi * irs_factor |> filter(timepoint == "127_3h") |> pull(irs_factor_Nt_1),
    deg_128_6h_TMTi_irs = deg_128_6h_TMTi * irs_factor |> filter(timepoint == "128_6h") |> pull(irs_factor_Nt_1),
    deg_129_9h_TMTi_irs = deg_129_9h_TMTi * irs_factor |> filter(timepoint == "129_9h") |> pull(irs_factor_Nt_1),
    deg_130_12h_TMTi_irs = deg_130_12h_TMTi * irs_factor |> filter(timepoint == "130_12h") |> pull(irs_factor_Nt_1),
    deg_131_24h_TMTi_irs = deg_131_24h_TMTi * irs_factor |> filter(timepoint == "131_24h") |> pull(irs_factor_Nt_1)
  )

write_csv(HEK_Nterm_1_irs, file = 'data_source/normalized_data/HEK_Nterm_1_irs.csv')

#HEK_Nt_2
HEK_Nterm_2_irs <- HEK_Nterm_2 |> 
  mutate(
    deg_126_0h_TMTi_irs = deg_126_0h_TMTi * irs_factor |> filter(timepoint == "126_0h") |> pull(irs_factor_Nt_2),
    deg_127_3h_TMTi_irs = deg_127_3h_TMTi * irs_factor |> filter(timepoint == "127_3h") |> pull(irs_factor_Nt_2),
    deg_128_6h_TMTi_irs = deg_128_6h_TMTi * irs_factor |> filter(timepoint == "128_6h") |> pull(irs_factor_Nt_2),
    deg_129_9h_TMTi_irs = deg_129_9h_TMTi * irs_factor |> filter(timepoint == "129_9h") |> pull(irs_factor_Nt_2),
    deg_130_12h_TMTi_irs = deg_130_12h_TMTi * irs_factor |> filter(timepoint == "130_12h") |> pull(irs_factor_Nt_2),
    deg_131_24h_TMTi_irs = deg_131_24h_TMTi * irs_factor |> filter(timepoint == "131_24h") |> pull(irs_factor_Nt_2)
  )

write_csv(HEK_Nterm_2_irs, file = 'data_source/normalized_data/HEK_Nterm_2_irs.csv')

#HEK_Nt_3
HEK_Nterm_3_irs <- HEK_Nterm_3 |> 
  mutate(
    deg_126_0h_TMTi_irs = deg_126_0h_TMTi * irs_factor |> filter(timepoint == "126_0h") |> pull(irs_factor_Nt_3),
    deg_127_3h_TMTi_irs = deg_127_3h_TMTi * irs_factor |> filter(timepoint == "127_3h") |> pull(irs_factor_Nt_3),
    deg_128_6h_TMTi_irs = deg_128_6h_TMTi * irs_factor |> filter(timepoint == "128_6h") |> pull(irs_factor_Nt_3),
    deg_129_9h_TMTi_irs = deg_129_9h_TMTi * irs_factor |> filter(timepoint == "129_9h") |> pull(irs_factor_Nt_3),
    deg_130_12h_TMTi_irs = deg_130_12h_TMTi * irs_factor |> filter(timepoint == "130_12h") |> pull(irs_factor_Nt_3),
    deg_131_24h_TMTi_irs = deg_131_24h_TMTi * irs_factor |> filter(timepoint == "131_24h") |> pull(irs_factor_Nt_3)
  )

write_csv(HEK_Nterm_3_irs, file = 'data_source/normalized_data/HEK_Nterm_3_irs.csv')

##coefficient of variance
#overlap of quantified N-term in three replicates
overlap_Nterm <- Reduce(intersect, list(
  HEK_Nterm_1_irs$Index,
  HEK_Nterm_2_irs$Index,
  HEK_Nterm_3_irs$Index
))

#extract data from each result table
HEK_Nterm_1_irs_overlap <- HEK_Nterm_1_irs |> 
  filter(Index %in% overlap_Nterm)

HEK_Nterm_2_irs_overlap <- HEK_Nterm_2_irs |> 
  filter(Index %in% overlap_Nterm)

HEK_Nterm_3_irs_overlap <- HEK_Nterm_3_irs |> 
  filter(Index %in% overlap_Nterm)

#define make_CV function
make_CV <- function(df1, df2, df3) {
  
  #extract data
  df_126_0h <- log10(
    bind_cols(
      df1 |> select(deg_126_0h_TMTi_irs),
      df2 |> select(deg_126_0h_TMTi_irs),
      df3 |> select(deg_126_0h_TMTi_irs)
    )
  )
  
  df_127_3h <- log10(
    bind_cols(
      df1 |> select(deg_127_3h_TMTi_irs),
      df2 |> select(deg_127_3h_TMTi_irs),
      df3 |> select(deg_127_3h_TMTi_irs)
    )
  )
  
  df_128_6h <- log10(
    bind_cols(
      df1 |> select(deg_128_6h_TMTi_irs),
      df2 |> select(deg_128_6h_TMTi_irs),
      df3 |> select(deg_128_6h_TMTi_irs)
    )
  )
  
  df_129_9h <- log10(
    bind_cols(
      df1 |> select(deg_129_9h_TMTi_irs),
      df2 |> select(deg_129_9h_TMTi_irs),
      df3 |> select(deg_129_9h_TMTi_irs)
    )
  )
  
  df_130_12h <- log10(
    bind_cols(
      df1 |> select(deg_130_12h_TMTi_irs),
      df2 |> select(deg_130_12h_TMTi_irs),
      df3 |> select(deg_130_12h_TMTi_irs)
    )
  )
  
  df_131_24h <- log10(
    bind_cols(
      df1 |> select(deg_131_24h_TMTi_irs),
      df2 |> select(deg_131_24h_TMTi_irs),
      df3 |> select(deg_131_24h_TMTi_irs)
    )
  )
  
  #calculate average, sd and cv
  df_126_0h$ave <- rowMeans(df_126_0h)
  df_126_0h$sd <- apply(df_126_0h[1:3], 1, sd)
  df_126_0h$cv <- 100 * df_126_0h$sd / df_126_0h$ave
  
  df_127_3h$ave <- rowMeans(df_127_3h)
  df_127_3h$sd <- apply(df_127_3h[1:3], 1, sd)
  df_127_3h$cv <- 100 * df_127_3h$sd / df_127_3h$ave
  
  df_128_6h$ave <- rowMeans(df_128_6h)
  df_128_6h$sd <- apply(df_128_6h[1:3], 1, sd)
  df_128_6h$cv <- 100 * df_128_6h$sd / df_128_6h$ave
  
  df_129_9h$ave <- rowMeans(df_129_9h)
  df_129_9h$sd <- apply(df_129_9h[1:3], 1, sd)
  df_129_9h$cv <- 100 * df_129_9h$sd / df_129_9h$ave
  
  df_130_12h$ave <- rowMeans(df_130_12h)
  df_130_12h$sd <- apply(df_130_12h[1:3], 1, sd)
  df_130_12h$cv <- 100 * df_130_12h$sd / df_130_12h$ave
  
  df_131_24h$ave <- rowMeans(df_131_24h)
  df_131_24h$sd <- apply(df_131_24h[1:3], 1, sd)
  df_131_24h$cv <- 100 * df_131_24h$sd / df_131_24h$ave
  
  #generate output
  ave_df <- data.frame(df_126_0h$ave, df_127_3h$ave, df_128_6h$ave, df_129_9h$ave, df_130_12h$ave, df_131_24h$ave)
  sd_df <- data.frame(df_126_0h$sd, df_127_3h$sd, df_128_6h$sd, df_129_9h$sd, df_130_12h$sd, df_131_24h$sd)
  cv_df <- data.frame(df_126_0h$cv, df_127_3h$cv, df_128_6h$cv, df_129_9h$cv, df_130_12h$cv, df_131_24h$cv)
  return(list(ave_df, sd_df, cv_df))
  
}

#calculate cv for the overlap quantified N-term in three replicates
list_ave_sd_cv_overlap_Nterm <- make_CV(
  df1 = HEK_Nterm_1_irs_overlap,
  df2 = HEK_Nterm_2_irs_overlap,
  df3 = HEK_Nterm_3_irs_overlap
)

#boxplot for cv
cv_overlap_Nterm <- tibble(list_ave_sd_cv_overlap_Nterm[[3]]) |> 
  pivot_longer(ends_with("cv"), names_to = "Exp", values_to = "CV")

boxplot_cv_overlap_Nterm <- cv_overlap_Nterm |> 
  ggplot() +
  geom_boxplot(aes(x = Exp, y = CV)) +
  labs(x = "", y = "CV (%)") +
  coord_cartesian(ylim = c(0, 25)) +
  theme(
    axis.text.x = element_text(angle = 30, hjust = 1)
  )

ggsave(
  filename = 'figures/boxplot_cv_overlap_Nterm.png',
  plot = boxplot_cv_overlap_Nterm,
  height = 640, width = 640, dpi = 300, units = 'px'
)
