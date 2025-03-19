# import packages
import pandas as pd
import numpy as np
import tqdm as tqdm

# import strucutremap functions as instructed by tutorial 
# (https://github.com/MannLabs/structuremap/blob/main/nbs/tutorial.ipynb)
import structuremap.utils
structuremap.utils.set_logger()
from structuremap.processing import download_alphafold_cif, download_alphafold_pae, format_alphafold_data, annotate_accessibility, get_smooth_score, annotate_proteins_with_idr_pattern, get_extended_flexible_pattern, get_proximity_pvals, perform_enrichment_analysis, perform_enrichment_analysis_per_protein, evaluate_ptm_colocalization, extract_motifs_in_proteome

# import Nterm Whole Proteome comparison results
HEK_Nterm_WP_delta_half_life = pd.read_csv(
  "data_source/Nterm_WP_comparison/HEK_Nterm_WP_delta_half_life.csv",
  header = 0,
  index_col = None
)

# top 20%, Nterm shorter half-life
HEK_Nterm_WP_delta_half_life["top20"] = np.where(
  HEK_Nterm_WP_delta_half_life.Percentile < 0.2, 1, 0
)

# bottom 20%, Nterm longer half-life
HEK_Nterm_WP_delta_half_life["bottom20"] = np.where(
  HEK_Nterm_WP_delta_half_life.Percentile > 0.8, 1, 0
)

# get overlap protein
overlap_protein_list = HEK_Nterm_WP_delta_half_life["UniProt_Accession"].unique().tolist()

# format AlphaFold data input
HEK_Nterm_WP_delta_half_life_alphafold_annotation = format_alphafold_data(
  directory = "data_source/Nterm_structuremap/Nterm_degradation_cif",
  protein_ids = overlap_protein_list
)

# annotate prediction-aware part-sphere exposure (pPSE) values
HEK_Nterm_WP_delta_half_life_full_sphere_exposure = annotate_accessibility(
    df = HEK_Nterm_WP_delta_half_life_alphafold_annotation,
    max_dist = 24,
    max_angle = 180,
    error_dir = "data_source/Nterm_structuremap/Nterm_degradation_pae"
)

HEK_Nterm_WP_delta_half_life_alphafold_accessibility = HEK_Nterm_WP_delta_half_life_alphafold_annotation.merge(
  HEK_Nterm_WP_delta_half_life_full_sphere_exposure,
  how = "left",
  on = ["protein_id", "AA", "position"]
)

HEK_Nterm_WP_delta_half_life_part_sphere_exposure = annotate_accessibility(
    df = HEK_Nterm_WP_delta_half_life_alphafold_annotation,
    max_dist = 12,
    max_angle = 70,
    error_dir = "data_source/Nterm_structuremap/Nterm_degradation_pae"
)

HEK_Nterm_WP_delta_half_life_alphafold_accessibility = HEK_Nterm_WP_delta_half_life_alphafold_accessibility.merge(
  HEK_Nterm_WP_delta_half_life_part_sphere_exposure,
  how = "left",
  on = ["protein_id", "AA", "position"]
)

HEK_Nterm_WP_delta_half_life_alphafold_accessibility["high_acc_5"] = np.where(
  HEK_Nterm_WP_delta_half_life_alphafold_accessibility.nAA_12_70_pae <= 5, 1, 0
)

HEK_Nterm_WP_delta_half_life_alphafold_accessibility["low_acc_5"] = np.where(
  HEK_Nterm_WP_delta_half_life_alphafold_accessibility.nAA_12_70_pae > 5, 1, 0
)

# annotate instriscly disorder region (IDR)
HEK_Nterm_WP_delta_half_life_alphafold_accessibility_smooth = get_smooth_score(
  HEK_Nterm_WP_delta_half_life_alphafold_accessibility,
  np.array(['nAA_24_180_pae']),
  [10]
).reset_index(drop=True)

HEK_Nterm_WP_delta_half_life_alphafold_accessibility_smooth["IDR"] = np.where(
  HEK_Nterm_WP_delta_half_life_alphafold_accessibility_smooth.nAA_24_180_pae_smooth10 <= 34.27, 1, 0
)

# annotate short IDRs
HEK_Nterm_WP_delta_half_life_alphafold_accessibility_smooth_pattern = annotate_proteins_with_idr_pattern(
  HEK_Nterm_WP_delta_half_life_alphafold_accessibility_smooth,
  min_structured_length = 80,
  max_unstructured_length = 20
).reset_index(drop=True)

HEK_Nterm_WP_delta_half_life_alphafold_accessibility_smooth_pattern_extended = get_extended_flexible_pattern(
  HEK_Nterm_WP_delta_half_life_alphafold_accessibility_smooth_pattern,
  ["flexible_pattern"],
  [5]
).reset_index(drop=True)

Nterm_site_dict = {
  "top20" : [
    "A",
    "R",
    "N",
    "D",
    "C",
    "E",
    "Q",
    "G",
    "H",
    "I",
    "L",
    "K",
    "M",
    "F",
    "P",
    "S",
    "T",
    "W",
    "Y",
    "V"
  ],
  "bottom20" : [
    "A",
    "R",
    "N",
    "D",
    "C",
    "E",
    "Q",
    "G",
    "H",
    "I",
    "L",
    "K",
    "M",
    "F",
    "P",
    "S",
    "T",
    "W",
    "Y",
    "V"
  ]
}

# annotate top20 and bottom20 N-terminus
Nterm_WP_delta_half_life_alphafold_N_terminus = HEK_Nterm_WP_delta_half_life_alphafold_accessibility_smooth_pattern_extended.merge(
  HEK_Nterm_WP_delta_half_life,
  how = "left",
  left_on = ["protein_id", "position"],
  right_on = ["UniProt_Accession", "Protein.Start"]
)

Nterm_WP_delta_half_life_alphafold_N_terminus = Nterm_WP_delta_half_life_alphafold_N_terminus.fillna(0)

# N-terminus enrichment analysis
rois_list = ["BEND", "HELX", "STRN", "TURN", "IDR", "high_acc_5", "low_acc_5", "flexible_pattern"]

enrichment_results = {}

for roi in tqdm(rois_list):
    enrichment_results[f"enrichment_top20_bottom20_{roi}"] = perform_enrichment_analysis(
        df=Nterm_WP_delta_half_life_alphafold_N_terminus,
        ptm_types=["top20", "bottom20"],
        rois=[roi],
        ptm_site_dict=Nterm_site_dict,
        quality_cutoffs=[0]
    )

# N-terminus 3D proximity
enrichment_top20_bottom20_proximity = get_proximity_pvals(
  df = Nterm_WP_delta_half_life_alphafold_N_terminus,
  ptm_types = ["top20", "bottom20"],
  ptm_site_dict = Nterm_site_dict,
  error_dir = "data_source/Nterm_structuremap/Nterm_degradation_pae", 
  per_site_metric = 'mean',
  error_operation = 'plus',
  n_random = 10000,
  random_seed = 44
)
