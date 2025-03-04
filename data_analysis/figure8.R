# import packages
library(tidyverse)

### figure 8A, GO and KEGG analysis for top20% and bottom 20% proteins
## import GO and KEGG enrichment results
top20_comparison_GO <- read_csv(
  'data_source/Nterm_WP_comparison/top20_comparison_GO.csv'
)

bottom20_comparison_GO <- read_csv(
  'data_source/Nterm_WP_comparison/bottom20_comparison_GO.csv'
)

top20_comparison_KEGG <- read_csv(
  'data_source/Nterm_WP_comparison/top20_comparison_KEGG.csv'
)

bottom20_comparison_KEGG <- read_csv(
  'data_source/Nterm_WP_comparison/bottom20_comparison_KEGG.csv'
)

## dot plot
top20_bottom20_comb <- bind_rows(
  top20_comparison_GO |> 
    filter(
      Description %in% c(
        'protein folding',
        'endoplasmic reticulum lumen',
        'vesicle lumen',
        'secretory granule lumen',
        'extracellular organelle',
        'proteasome complex',
        'proteasome regulatory particle, base subcomplex'
      )
    ) |> 
    mutate(
      category = 'top20'
    ),
  
  top20_comparison_KEGG |> 
    filter(
      Description %in% c(
        'Proteasome',
        'Protein processing in endoplasmic reticulum'
      )
    ) |> 
    mutate(
      category = 'top20'
    ),
  
  bottom20_comparison_GO |> 
    filter(
      Description %in% c(
        'nucleotide binding',
        'anion binding',
        'ribonucleotide binding',
        'carbohydrate derivative binding',
        'ATP binding',
        'helicase activity',
        'ATP hydrolysis activity',
        'pyrophosphatase activity',
        'structural constituent of ribosome'
      )
    ) |> 
    mutate(
      category = 'bottom20'
    )
)

write_csv(
  top20_bottom20_comb,
  file = 'data_source/Nterm_WP_comparison/top20_bottom20_comb.csv'
)

dotplot_top20_bottom20_GO_KEGG <- top20_bottom20_comb |> 
  mutate(
    FoldEnrichment = ifelse(category == 'bottom20', -FoldEnrichment, FoldEnrichment)
  ) |> 
  ggplot() +
  geom_point(
    aes(
      x = FoldEnrichment,
      y = -log10(p.adjust),
      size = Count,
      fill = category
    ),
    shape = 21
  ) +
  geom_vline(xintercept = 0, color = 'black', linetype = 'dashed') +
  theme(
    axis.text = element_text(size = 8, color = 'black', family = 'arial'),
    axis.title = element_text(size = 8, color = 'black', family = 'arial'),
    legend.title = element_text(size = 8, color = 'black', family = 'arial'),
    legend.text = element_text(size = 8, color = 'black', family = 'arial')
  )

ggsave(
  filename = 'figures/figure8/dotplot_top20_bottom20_GO_KEGG.eps',
  plot = dotplot_top20_bottom20_GO_KEGG,
  height = 2, width = 3.5, units = 'in'
)

### figure 8B, protein complex analysis for top20% and bottom 20% proteins
## import corum enrichment results
top20_comparison_corum <- read_csv(
  'data_source/Nterm_WP_comparison/top20_comparison_corum.csv'
)

bottom20_comparison_corum <- read_csv(
  'data_source/Nterm_WP_comparison/bottom20_comparison_corum.csv'
)

## bar plot
# top20
barplot_top20_corum <- top20_comparison_corum |> 
  filter(
    Description %in% c(
      'corum_id_181',
      'corum_id_8372',
      'corum_id_193',
      'corum_id_5615',
      'corum_id_8391'
    )
  ) |> 
  ggplot() +
  geom_bar(
    aes(
      x = fct_reorder(Description, -pvalue),
      y = -log10(pvalue)
    ),
    stat = 'identity', fill = color_3
  ) +
  labs(x = '', y = '') +
  coord_flip() +
  theme(
    axis.text = element_text(size = 8, color = 'black')
  )

ggsave(
  filename = 'figures/figure8/barplot_top20_corum.eps',
  plot = barplot_top20_corum,
  height = 2, width = 2.5, units = 'in'
)

# bottom20
barplot_bottom20_corum <- bottom20_comparison_corum |> 
  filter(
    Description %in% c(
      'corum_id_3055',
      'corum_id_306',
      'corum_id_305',
      'corum_id_3040',
      'corum_id_7627'
    )
  ) |> 
  ggplot() +
  geom_bar(
    aes(
      x = fct_reorder(Description, -pvalue),
      y = -log10(pvalue)
    ),
    stat = 'identity', fill = color_4
  ) +
  labs(x = '', y = '') +
  coord_flip() +
  theme(
    axis.text = element_text(size = 8, color = 'black')
  )

ggsave(
  filename = 'figures/figure8/barplot_bottom20_corum.eps',
  plot = barplot_bottom20_corum,
  height = 2, width = 2.5, units = 'in'
)
