# import packages
library(tidyverse)

### figure 1C, protein level workflow
# import open search result
# source data NGlyco_9_05012024
opensearch_result <- read_tsv(
  'data_source/open_search/global.modsummary.tsv',
  col_names = TRUE,
  name_repair = 'universal'
)

# summarize the Top5 detected modification PSMs
Top5_modification_psm <- opensearch_result |> 
  select(Modification:E_LF_NGlyco_9_3_05012024_3_PSMs) |> 
  slice(1:5) |> 
  pivot_longer(!Modification:Mass.Shift, names_to = 'Exp', values_to = 'PSMs') |> 
  group_by(Modification) |> 
  mutate(
    avg_PSMs = mean(PSMs)
  )

# bar plot
barplot_opensearch_psm <- Top5_modification_psm |> 
  ggplot() +
  geom_bar(
    aes(
      x = Mass.Shift, 
      y = avg_PSMs/3
    ), 
    stat = 'identity'
  ) +
  labs(x = 'Mass Shift', y = 'Number of PSMs') +
  scale_y_continuous(expand = c(0, 0)) +
  # theme_classic() +
  theme(
    panel.grid.major = element_line(color = "gray", linewidth = 0),
    panel.grid.minor = element_line(color = "gray", linewidth = 0),
    axis.title = element_text(size = 8),
    axis.text.x = element_text(size = 8, color = "black"),
    axis.text.y = element_text(size = 8, color = "black")
  )

ggsave(
  filename = 'figures/figure1/barplot_opensearch_psm.eps',
  device = 'eps',
  height = 2, width = 2.5, units = 'in'
)

### figure 1D, yeast peptide, specificity test
# import yeast peptide experiment result
# source data NGlyco_yeast_05012024
Yeast_Test <- tribble(
  ~ Exp, ~ Total, ~ N_term, ~ K,
  "Rep 1", 15738, 7598, 27,
  "Rep 2", 11602, 5293, 16,
  "Rep 3", 12409, 5771, 15
)

Yeast_Test <- Yeast_Test |> 
  pivot_longer(cols = c(N_term, K), names_to = "N_term_K", values_to = "Value") |> 
  mutate(Percentage_N_term_K = Value/Total, 
         Total_Percentage = Total/Total)

# bar plot
barplot_yeast_peptide_test <- ggplot() +
  geom_bar(data = Yeast_Test, 
           aes(x = Exp, y = Percentage_N_term_K, fill = N_term_K), stat = "identity") +
  geom_bar(data = Yeast_Test |> filter(N_term_K == "N_term"), 
           aes(x = Exp, y = Total_Percentage), stat = "identity", fill = NA, color = "black") +
  labs(x = "", y = "Modified Percentage", fill = "Residues") +
  scale_fill_manual(values = c("N_term" = color_1, "K" = color_2)) +
  theme(
    axis.title = element_text(size = 8),
    axis.text.x = element_text(size = 8, color = "black", angle = 30, hjust = 1, family = 'arial'),
    axis.text.y = element_text(size = 8, color = "black", family = 'arial'),
    legend.title = element_text(size = 8, color = "black", family = 'arial'),
    legend.text = element_text(size = 8, color = "black", family = 'arial')
  )

ggsave(
  filename = 'figures/figure1/barplot_yeast_peptide_test.eps', 
  plot = barplot_yeast_peptide_test, 
  height = 1.5, width = 2.5, units = "in"
)

### figure 1E, experimental condition optimization
# import optimization result
# source data NGlyco_*_05012024 (* represents 6, 7, 8 or 9)
Different_pH <- tribble(
  ~ Exp, ~ pH, ~ Value,
  "Exp_1", "pH 6.5", 772,
  "Exp_2", "pH 6.5", 827,
  "Exp_3", "pH 6.5", 754,
  "Exp_1", "pH 7.4", 2364,
  "Exp_2", "pH 7.4", 2390,
  "Exp_3", "pH 7.4", 2533,
  "Exp_1", "pH 8.0", 3177,
  "Exp_2", "pH 8.0", 3329,
  "Exp_3", "pH 8.0", 3575,
  "Exp_1", "pH 8.6", 3270,
  "Exp_2", "pH 8.6", 2696,
  "Exp_3", "pH 8.6", 2537,
)

Different_pH_average <- Different_pH |> 
  group_by(pH) |> 
  mutate(avg = mean(Value))

# bar plot
barplot_different_pH <- ggplot() +
  geom_bar(data = Different_pH_average |> filter(Exp == "Exp_1"),
           aes(x = pH, y = avg), stat = "identity", fill = color_1) +
  geom_point(data = Different_pH, 
             aes(x = pH, y = Value)) +
  labs(x = "", y = "Identified Unique N-term") +
  theme(
    axis.title = element_text(size = 8),
    axis.text.x = element_text(size = 8, color = "black", angle = 30, hjust = 1),
    axis.text.y = element_text(size = 8, color = "black"),
  )

ggsave(
  filename = 'figures/figure1/barplot_different_pH.eps', 
  plot = barplot_different_pH, 
  height = 2, width = 1.5, units = "in"
)

### figure 1F, 2PCA alkyne comparison experiment
# import comparison experiment result
# source data 2PCA_05292024 and NGlyco_05292024
Glyco_2PCA_Comparison <- tribble(
  ~ Condition, ~ Exp, ~ Value,
  "Glyco, pH 8.0", "Exp 1", 3464,
  "Glyco, pH 8.0", "Exp 2", 3566,
  "Glyco, pH 8.0", "Exp 3", 3367,
  "2PCA, pH 7.4", "Exp 1", 1992,
  "2PCA, pH 7.4", "Exp 2", 1874,
  "2PCA, pH 7.4", "Exp 3", 2129
)

Glyco_2PCA_Comparison_Average <- Glyco_2PCA_Comparison |> 
  group_by(Condition) |> 
  mutate(avg = mean(Value))

# bar plot
barplot_Glyco_2PCA_comparison <- ggplot() +
  geom_bar(data = Glyco_2PCA_Comparison_Average |> filter(Exp == "Exp 1"),
           aes(x = Condition, y = avg, fill = Condition), stat = "identity") +
  geom_point(data = Glyco_2PCA_Comparison, 
             aes(x = Condition, y = Value)) +
  labs(x = "", y = "Identified Unique N-term") +
  scale_fill_manual(values = c(
    '2PCA, pH 7.4' = color_2,
    'Glyco, pH 8.0' = color_1
  )) +
  theme(
    axis.title = element_text(size = 8),
    axis.text.x = element_text(size = 8, color = "black", angle = 15, hjust = 1),
    axis.text.y = element_text(size = 8, color = "black"),
    legend.position = 'none'
  )

ggsave(
  filename = 'figures/figure1/barplot_Glyco_2PCA_comparison.eps', 
  plot = barplot_Glyco_2PCA_comparison, 
  height = 2, width = 1.5, units = "in"
)

