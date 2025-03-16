# import packages
library(tidyverse)

### figure S3A, CV for Nterm
## coefficient of variance
# overlap of quantified N-term in three replicates
overlap_Nterm <- Reduce(intersect, list(
  HEK_Nterm_1_irs$Index,
  HEK_Nterm_2_irs$Index,
  HEK_Nterm_3_irs$Index
))

# extract data from each result table
HEK_Nterm_1_irs_overlap <- HEK_Nterm_1_irs |> 
  filter(Index %in% overlap_Nterm)

HEK_Nterm_2_irs_overlap <- HEK_Nterm_2_irs |> 
  filter(Index %in% overlap_Nterm)

HEK_Nterm_3_irs_overlap <- HEK_Nterm_3_irs |> 
  filter(Index %in% overlap_Nterm)

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

# calculate cv for the overlap quantified N-term in three replicates
list_ave_sd_cv_overlap_Nterm <- make_CV(
  df1 = HEK_Nterm_1_irs_overlap,
  df2 = HEK_Nterm_2_irs_overlap,
  df3 = HEK_Nterm_3_irs_overlap
)

# boxplot for cv
cv_overlap_Nterm <- tibble(list_ave_sd_cv_overlap_Nterm[[3]]) |> 
  pivot_longer(ends_with("cv"), names_to = "Exp", values_to = "CV")

boxplot_cv_overlap_Nterm <- cv_overlap_Nterm |> 
  ggplot() +
  geom_boxplot(aes(x = Exp, y = CV)) +
  labs(x = "", y = "CV (%)") +
  coord_cartesian(ylim = c(0, 25)) +
  theme(
    axis.text.x = element_text(size = 8, angle = 30, hjust = 1, family = 'arial'),
    axis.text.y = element_text(size = 8, family = 'arial')
  )

ggsave(
  filename = 'figures/figureS3/boxplot_cv_overlap_Nterm.eps',
  plot = boxplot_cv_overlap_Nterm,
  height = 2, width = 2, units = 'in'
)

### figure S3B, CV for Whole Proteome
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
    axis.text.y = element_text(size = 8, family = 'arial')
  )

ggsave(
  filename = 'figures/figureS3/boxplot_cv_overlap_Whole_Proteome.eps',
  plot = boxplot_cv_overlap_Whole_Proteome,
  height = 2, width = 2, units = 'in'
)

