# import packages
library(tidyverse)

## TMT channel intensity normalization
# HEK_WP_1
target_mean_deg_TMTi_WP_1 <- mean(colSums(HEK_WP_1_psm |> select(starts_with('deg'))))
norm_facs_WP_1 <- target_mean_deg_TMTi_WP_1/colSums(HEK_WP_1_psm |> select(starts_with('deg')))
deg_TMTi_WP_1 <- tibble(sweep(HEK_WP_1_psm |> select(starts_with('deg')), 2, norm_facs_WP_1, FUN = '*'))
colnames(deg_TMTi_WP_1) <- c(
  'deg_126_0h_TMTi', 'deg_127_3h_TMTi', 'deg_128_6h_TMTi', 'deg_129_9h_TMTi', 'deg_130_12h_TMTi', 'deg_131_24h_TMTi'
)

HEK_WP_1_psm_TMTi <- bind_cols(HEK_WP_1_psm, deg_TMTi_WP_1)
write_csv(HEK_WP_1_psm_TMTi, file = 'data_source/normalized_data/HEK_WP_1_psm_TMTi.csv')

# HEK_WP_2
target_mean_deg_TMTi_WP_2 <- mean(colSums(HEK_WP_2_psm |> select(starts_with('deg'))))
norm_facs_WP_2 <- target_mean_deg_TMTi_WP_2/colSums(HEK_WP_2_psm |> select(starts_with('deg')))
deg_TMTi_WP_2 <- tibble(sweep(HEK_WP_2_psm |> select(starts_with('deg')), 2, norm_facs_WP_2, FUN = '*'))
colnames(deg_TMTi_WP_2) <- c(
  'deg_126_0h_TMTi', 'deg_127_3h_TMTi', 'deg_128_6h_TMTi', 'deg_129_9h_TMTi', 'deg_130_12h_TMTi', 'deg_131_24h_TMTi'
)

HEK_WP_2_psm_TMTi <- bind_cols(HEK_WP_2_psm, deg_TMTi_WP_2)
write_csv(HEK_WP_2_psm_TMTi, file = 'data_source/normalized_data/HEK_WP_2_psm_TMTi.csv')

# HEK_WP_3
target_mean_deg_TMTi_WP_3 <- mean(colSums(HEK_WP_3_psm |> select(starts_with('deg'))))
norm_facs_WP_3 <- target_mean_deg_TMTi_WP_3/colSums(HEK_WP_3_psm |> select(starts_with('deg')))
deg_TMTi_WP_3 <- tibble(sweep(HEK_WP_3_psm |> select(starts_with('deg')), 2, norm_facs_WP_3, FUN = '*'))
colnames(deg_TMTi_WP_3) <- c(
  'deg_126_0h_TMTi', 'deg_127_3h_TMTi', 'deg_128_6h_TMTi', 'deg_129_9h_TMTi', 'deg_130_12h_TMTi', 'deg_131_24h_TMTi'
)

HEK_WP_3_psm_TMTi <- bind_cols(HEK_WP_3_psm, deg_TMTi_WP_3)
write_csv(HEK_WP_3_psm_TMTi, file = 'data_source/normalized_data/HEK_WP_3_psm_TMTi.csv')

## Remove all the heavy labeled and not N-term modified psm
# HEK_WP_1
HEK_WP_1_psm_TMTi_light <- HEK_WP_1_psm_TMTi |> 
  filter(!str_detect(Assigned.Modifications, "K\\(8.0142\\)")) |> 
  filter(!str_detect(Assigned.Modifications, "K\\(237.1771\\)")) |> 
  filter(!str_detect(Assigned.Modifications, "R\\(10.0082\\)"))

write_csv(HEK_WP_1_psm_TMTi_light, file = 'data_source/normalized_data/HEK_WP_1_psm_TMTi_light.csv')

# HEK_WP_2
HEK_WP_2_psm_TMTi_light <- HEK_WP_2_psm_TMTi |> 
  filter(!str_detect(Assigned.Modifications, "K\\(8.0142\\)")) |> 
  filter(!str_detect(Assigned.Modifications, "K\\(237.1771\\)")) |> 
  filter(!str_detect(Assigned.Modifications, "R\\(10.0082\\)"))

write_csv(HEK_WP_2_psm_TMTi_light, file = 'data_source/normalized_data/HEK_WP_2_psm_TMTi_light.csv')

# HEK_WP_3
HEK_WP_3_psm_TMTi_light <- HEK_WP_3_psm_TMTi |> 
  filter(!str_detect(Assigned.Modifications, "K\\(8.0142\\)")) |> 
  filter(!str_detect(Assigned.Modifications, "K\\(237.1771\\)")) |> 
  filter(!str_detect(Assigned.Modifications, "R\\(10.0082\\)"))

write_csv(HEK_WP_3_psm_TMTi_light, file = 'data_source/normalized_data/HEK_WP_3_psm_TMTi_light.csv')

## group by UniProt_Accession
# HEK_WP_1
HEK_Whole_Proteome_1 <- HEK_WP_1_psm_TMTi_light |> 
  group_by(UniProt_Accession, Gene, Entry.Name) |> 
  summarise(
    deg_126_0h_TMTi = sum(deg_126_0h_TMTi),
    deg_127_3h_TMTi = sum(deg_127_3h_TMTi),
    deg_128_6h_TMTi = sum(deg_128_6h_TMTi),
    deg_129_9h_TMTi = sum(deg_129_9h_TMTi),
    deg_130_12h_TMTi = sum(deg_130_12h_TMTi),
    deg_131_24h_TMTi = sum(deg_131_24h_TMTi)
  ) |> 
  ungroup()

write_csv(HEK_Whole_Proteome_1, file = 'data_source/normalized_data/HEK_Whole_Proteome_1.csv')

# HEK_WP_2
HEK_Whole_Proteome_2 <- HEK_WP_2_psm_TMTi_light |> 
  group_by(UniProt_Accession, Gene, Entry.Name) |> 
  summarise(
    deg_126_0h_TMTi = sum(deg_126_0h_TMTi),
    deg_127_3h_TMTi = sum(deg_127_3h_TMTi),
    deg_128_6h_TMTi = sum(deg_128_6h_TMTi),
    deg_129_9h_TMTi = sum(deg_129_9h_TMTi),
    deg_130_12h_TMTi = sum(deg_130_12h_TMTi),
    deg_131_24h_TMTi = sum(deg_131_24h_TMTi)
  ) |> 
  ungroup()

write_csv(HEK_Whole_Proteome_2, file = 'data_source/normalized_data/HEK_Whole_Proteome_2.csv')

# HEK_WP_3
HEK_Whole_Proteome_3 <- HEK_WP_3_psm_TMTi_light |> 
  group_by(UniProt_Accession, Gene, Entry.Name) |> 
  summarise(
    deg_126_0h_TMTi = sum(deg_126_0h_TMTi),
    deg_127_3h_TMTi = sum(deg_127_3h_TMTi),
    deg_128_6h_TMTi = sum(deg_128_6h_TMTi),
    deg_129_9h_TMTi = sum(deg_129_9h_TMTi),
    deg_130_12h_TMTi = sum(deg_130_12h_TMTi),
    deg_131_24h_TMTi = sum(deg_131_24h_TMTi)
  ) |> 
  ungroup()

write_csv(HEK_Whole_Proteome_3, file = 'data_source/normalized_data/HEK_Whole_Proteome_3.csv')

## internal reference scaling normalization
# overlap of quantified proteins in three replicates
overlap_protein <- Reduce(intersect, list(
  HEK_Whole_Proteome_1$UniProt_Accession,
  HEK_Whole_Proteome_2$UniProt_Accession,
  HEK_Whole_Proteome_3$UniProt_Accession
))

HEK_Whole_Proteome_1_overlap <- HEK_Whole_Proteome_1 |> filter(UniProt_Accession %in% overlap_protein)
HEK_Whole_Proteome_2_overlap <- HEK_Whole_Proteome_2 |> filter(UniProt_Accession %in% overlap_protein)
HEK_Whole_Proteome_3_overlap <- HEK_Whole_Proteome_3 |> filter(UniProt_Accession %in% overlap_protein)

# calculate irs normalization factors
irs_factor <- tribble(
  ~ timepoint, ~ HEK_WP_1, ~ HEK_WP_2, ~ HEK_WP_3,
  '126_0h', sum(HEK_Whole_Proteome_1_overlap$deg_126_0h_TMTi), sum(HEK_Whole_Proteome_2_overlap$deg_126_0h_TMTi), sum(HEK_Whole_Proteome_3_overlap$deg_126_0h_TMTi),
  '127_3h', sum(HEK_Whole_Proteome_1_overlap$deg_127_3h_TMTi), sum(HEK_Whole_Proteome_2_overlap$deg_127_3h_TMTi), sum(HEK_Whole_Proteome_3_overlap$deg_127_3h_TMTi),
  '128_6h', sum(HEK_Whole_Proteome_1_overlap$deg_128_6h_TMTi), sum(HEK_Whole_Proteome_2_overlap$deg_128_6h_TMTi), sum(HEK_Whole_Proteome_3_overlap$deg_128_6h_TMTi),
  '129_9h', sum(HEK_Whole_Proteome_1_overlap$deg_129_9h_TMTi), sum(HEK_Whole_Proteome_2_overlap$deg_129_9h_TMTi), sum(HEK_Whole_Proteome_3_overlap$deg_129_9h_TMTi),
  '130_12h', sum(HEK_Whole_Proteome_1_overlap$deg_130_12h_TMTi), sum(HEK_Whole_Proteome_2_overlap$deg_130_12h_TMTi), sum(HEK_Whole_Proteome_3_overlap$deg_130_12h_TMTi),
  '131_24h', sum(HEK_Whole_Proteome_1_overlap$deg_131_24h_TMTi), sum(HEK_Whole_Proteome_2_overlap$deg_131_24h_TMTi), sum(HEK_Whole_Proteome_3_overlap$deg_131_24h_TMTi)
) |> 
  mutate(
    average = (HEK_WP_1 + HEK_WP_2 + HEK_WP_3)/3,
    irs_factor_WP_1 = average/HEK_WP_1,
    irs_factor_WP_2 = average/HEK_WP_2, 
    irs_factor_WP_3 = average/HEK_WP_3
  )

## apply irs normalization factors to data
# HEK_WP_1
HEK_Whole_Proteome_1_irs <- HEK_Whole_Proteome_1 |> 
  mutate(
    deg_126_0h_TMTi_irs = deg_126_0h_TMTi * irs_factor |> filter(timepoint == "126_0h") |> pull(irs_factor_WP_1),
    deg_127_3h_TMTi_irs = deg_127_3h_TMTi * irs_factor |> filter(timepoint == "127_3h") |> pull(irs_factor_WP_1),
    deg_128_6h_TMTi_irs = deg_128_6h_TMTi * irs_factor |> filter(timepoint == "128_6h") |> pull(irs_factor_WP_1),
    deg_129_9h_TMTi_irs = deg_129_9h_TMTi * irs_factor |> filter(timepoint == "129_9h") |> pull(irs_factor_WP_1),
    deg_130_12h_TMTi_irs = deg_130_12h_TMTi * irs_factor |> filter(timepoint == "130_12h") |> pull(irs_factor_WP_1),
    deg_131_24h_TMTi_irs = deg_131_24h_TMTi * irs_factor |> filter(timepoint == "131_24h") |> pull(irs_factor_WP_1)
  )

write_csv(HEK_Whole_Proteome_1_irs, file = 'data_source/normalized_data/HEK_Whole_Proteome_1_irs.csv')

# HEK_WP_2
HEK_Whole_Proteome_2_irs <- HEK_Whole_Proteome_2 |> 
  mutate(
    deg_126_0h_TMTi_irs = deg_126_0h_TMTi * irs_factor |> filter(timepoint == "126_0h") |> pull(irs_factor_WP_2),
    deg_127_3h_TMTi_irs = deg_127_3h_TMTi * irs_factor |> filter(timepoint == "127_3h") |> pull(irs_factor_WP_2),
    deg_128_6h_TMTi_irs = deg_128_6h_TMTi * irs_factor |> filter(timepoint == "128_6h") |> pull(irs_factor_WP_2),
    deg_129_9h_TMTi_irs = deg_129_9h_TMTi * irs_factor |> filter(timepoint == "129_9h") |> pull(irs_factor_WP_2),
    deg_130_12h_TMTi_irs = deg_130_12h_TMTi * irs_factor |> filter(timepoint == "130_12h") |> pull(irs_factor_WP_2),
    deg_131_24h_TMTi_irs = deg_131_24h_TMTi * irs_factor |> filter(timepoint == "131_24h") |> pull(irs_factor_WP_2)
  )

write_csv(HEK_Whole_Proteome_2_irs, file = 'data_source/normalized_data/HEK_Whole_Proteome_2_irs.csv')

# HEK_WP_3
HEK_Whole_Proteome_3_irs <- HEK_Whole_Proteome_3 |> 
  mutate(
    deg_126_0h_TMTi_irs = deg_126_0h_TMTi * irs_factor |> filter(timepoint == "126_0h") |> pull(irs_factor_WP_3),
    deg_127_3h_TMTi_irs = deg_127_3h_TMTi * irs_factor |> filter(timepoint == "127_3h") |> pull(irs_factor_WP_3),
    deg_128_6h_TMTi_irs = deg_128_6h_TMTi * irs_factor |> filter(timepoint == "128_6h") |> pull(irs_factor_WP_3),
    deg_129_9h_TMTi_irs = deg_129_9h_TMTi * irs_factor |> filter(timepoint == "129_9h") |> pull(irs_factor_WP_3),
    deg_130_12h_TMTi_irs = deg_130_12h_TMTi * irs_factor |> filter(timepoint == "130_12h") |> pull(irs_factor_WP_3),
    deg_131_24h_TMTi_irs = deg_131_24h_TMTi * irs_factor |> filter(timepoint == "131_24h") |> pull(irs_factor_WP_3)
  )

write_csv(HEK_Whole_Proteome_3_irs, file = 'data_source/normalized_data/HEK_Whole_Proteome_3_irs.csv')

## coefficient of variance
# overlap of quantified proteins in three replicates
overlap_protein <- Reduce(intersect, list(
  HEK_Whole_Proteome_1_irs$UniProt_Accession,
  HEK_Whole_Proteome_2_irs$UniProt_Accession,
  HEK_Whole_Proteome_3_irs$UniProt_Accession
))

# extract data from each result table
HEK_Whole_Proteome_1_irs_overlap <- HEK_Whole_Proteome_1_irs |> 
  filter(UniProt_Accession %in% overlap_protein)

HEK_Whole_Proteome_2_irs_overlap <- HEK_Whole_Proteome_2_irs |> 
  filter(UniProt_Accession %in% overlap_protein)

HEK_Whole_Proteome_3_irs_overlap <- HEK_Whole_Proteome_3_irs |> 
  filter(UniProt_Accession %in% overlap_protein)

# define make_CV function
make_CV <- function(df1, df2, df3) {
  
  # extract data
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
  
  # calculate average, sd and cv
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
  
  # generate output
  ave_df <- data.frame(df_126_0h$ave, df_127_3h$ave, df_128_6h$ave, df_129_9h$ave, df_130_12h$ave, df_131_24h$ave)
  sd_df <- data.frame(df_126_0h$sd, df_127_3h$sd, df_128_6h$sd, df_129_9h$sd, df_130_12h$sd, df_131_24h$sd)
  cv_df <- data.frame(df_126_0h$cv, df_127_3h$cv, df_128_6h$cv, df_129_9h$cv, df_130_12h$cv, df_131_24h$cv)
  return(list(ave_df, sd_df, cv_df))
  
}

# calculate cv for the overlap quantified proteins in three replicates
list_ave_sd_cv_overlap_Whole_Proteome <- make_CV(
  df1 = HEK_Whole_Proteome_1_irs_overlap,
  df2 = HEK_Whole_Proteome_2_irs_overlap,
  df3 = HEK_Whole_Proteome_3_irs_overlap
)

# boxplot for cv
cv_overlap_protein <- tibble(list_ave_sd_cv_overlap_Whole_Proteome[[3]]) |> 
  pivot_longer(ends_with("cv"), names_to = "Exp", values_to = "CV")

boxplot_cv_overlap_Whole_Proteome <- cv_overlap_protein |> 
  ggplot() +
  geom_boxplot(aes(x = Exp, y = CV)) +
  labs(x = "", y = "CV (%)") +
  coord_cartesian(ylim = c(0, 20)) +
  theme(
    axis.text.x = element_text(size = 8, angle = 30, hjust = 1, family = 'arial'),
    axis.text.y = element_text(size = 8, family = 'arial'),
  )

ggsave(
  filename = 'figures/boxplot_cv_overlap_Whole_Proteome.png',
  plot = boxplot_cv_overlap_Whole_Proteome,
  height = 640, width = 640, dpi = 300, units = 'px'
)
