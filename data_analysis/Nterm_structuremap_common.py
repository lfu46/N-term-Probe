# import packages
import pandas as pd
import numpy as np

# import strucutremap functions as instructed by tutorial 
# (https://github.com/MannLabs/structuremap/blob/main/nbs/tutorial.ipynb)
import structuremap.utils
structuremap.utils.set_logger()
from structuremap.processing import download_alphafold_cif, download_alphafold_pae, format_alphafold_data, annotate_accessibility, get_smooth_score, annotate_proteins_with_idr_pattern, get_extended_flexible_pattern, get_proximity_pvals, perform_enrichment_analysis, perform_enrichment_analysis_per_protein, evaluate_ptm_colocalization, extract_motifs_in_proteome

# import common Nterm data from results
common_Nterm_protein = pd.read_csv(
  'data_source/common_Nterm/common_Nterm_protein.csv',
  header = 0,
  index_col = None
)

common_Nterm_protein_list = common_Nterm_protein["UniProt_Accession"].unique().tolist()

## download AlphaFold data
# crystallographic information file
# valid_proteins_cif, invalid_proteins_cif, existing_proteins_cif = download_alphafold_cif(
#     proteins = total_Nterm_protein_list,
#     out_folder = "data_source/Nterm_structuremap/Nterm_degradation_cif"
# )

# predicted aligned error
# valid_proteins_pae, invalid_proteins_pae, existing_proteins_pae = download_alphafold_pae(
#     proteins = total_Nterm_protein_list,
#     out_folder = "data_source/Nterm_structuremap/Nterm_degradation_pae"
# )

# format AlphaFold data input
common_Nterm_alphafold_annotation = format_alphafold_data(
  directory = "data_source/Nterm_structuremap/Nterm_degradation_cif",
  protein_ids = common_Nterm_protein_list
)

# annotate prediction-aware part-sphere exposure (pPSE) values
common_Nterm_full_sphere_exposure = annotate_accessibility(
    df = common_Nterm_alphafold_annotation,
    max_dist = 24,
    max_angle = 180,
    error_dir = "data_source/Nterm_structuremap/Nterm_degradation_pae"
)

common_Nterm_alphafold_accessibility = common_Nterm_alphafold_annotation.merge(
  common_Nterm_full_sphere_exposure,
  how = "left",
  on = ["protein_id", "AA", "position"]
)

common_Nterm_part_sphere_exposure = annotate_accessibility(
    df = common_Nterm_alphafold_annotation,
    max_dist = 12,
    max_angle = 70,
    error_dir = "data_source/Nterm_structuremap/Nterm_degradation_pae"
)

common_Nterm_alphafold_accessibility = common_Nterm_alphafold_accessibility.merge(
  common_Nterm_part_sphere_exposure,
  how = "left",
  on = ["protein_id", "AA", "position"]
)

common_Nterm_alphafold_accessibility["high_acc_5"] = np.where(
  common_Nterm_alphafold_accessibility.nAA_12_70_pae <= 5, 1, 0
)

common_Nterm_alphafold_accessibility["low_acc_5"] = np.where(
  common_Nterm_alphafold_accessibility.nAA_12_70_pae > 5, 1, 0
)

# annotate instrisicly disorder region (IDR)
common_Nterm_alphafold_accessibility_smooth = get_smooth_score(
  common_Nterm_alphafold_accessibility,
  np.array(['nAA_24_180_pae']),
  [10]
).reset_index(drop=True)

common_Nterm_alphafold_accessibility_smooth["IDR"] = np.where(
  common_Nterm_alphafold_accessibility_smooth.nAA_24_180_pae_smooth10 <= 34.27, 1, 0
)

# anntate short IDRs
# common_Nterm_alphafold_accessibility_smooth_pattern = annotate_proteins_with_idr_pattern(
#   common_Nterm_alphafold_accessibility_smooth,
#   min_structured_length = 80,
#   max_unstructured_length = 20
# )

# common_Nterm_alphafold_accessibility_smooth_pattern_extended = get_extended_flexible_pattern(
#   common_Nterm_alphafold_accessibility_smooth_pattern,
#   ["flexible_pattern"],
#   [5]
# )

# annotate common Nterm data
common_Nterm_alphafold_N_terminus = common_Nterm_alphafold_accessibility_smooth.merge(
  common_Nterm_protein,
  how = "left",
  left_on = ["protein_id", "position"],
  right_on = ["UniProt_Accession", "start.position"]
)

common_Nterm_alphafold_N_terminus = common_Nterm_alphafold_N_terminus.fillna(0)

common_Nterm_alphafold_N_terminus["common_Nterm"] = common_Nterm_alphafold_N_terminus["common_Nterm"].apply(lambda x: 1 if x != 0 else x)

Nterm_site_dict = {
  "common_Nterm" : [
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

# perform N-terminus enrichement analysis
enrichment_N_terminus = perform_enrichment_analysis(
  df = common_Nterm_alphafold_N_terminus,
  ptm_types = ["common_Nterm"],
  rois = ["BEND", "HELX", "STRN", "TURN", "IDR", "high_acc_5", "low_acc_5"],
  ptm_site_dict = Nterm_site_dict,
  quality_cutoffs = [0]
)

# HELX
# enrichment_N_terminus_HELX = perform_enrichment_analysis(
#   df = common_Nterm_alphafold_N_terminus,
#   ptm_types = ["common_Nterm"],
#   rois = ["HELX"],
#   ptm_site_dict = Nterm_site_dict,
#   quality_cutoffs = [0]
# )

# STRN
# enrichment_N_terminus_STRN = perform_enrichment_analysis(
#   df = common_Nterm_alphafold_N_terminus,
#   ptm_types = ["common_Nterm"],
#   rois = ["STRN"],
#   ptm_site_dict = Nterm_site_dict,
#   quality_cutoffs = [0]
# )

# TURN
# enrichment_N_terminus_TURN = perform_enrichment_analysis(
#   df = common_Nterm_alphafold_N_terminus,
#   ptm_types = ["common_Nterm"],
#   rois = ["TURN"],
#   ptm_site_dict = Nterm_site_dict,
#   quality_cutoffs = [0]
# )

# IDR
# enrichment_N_terminus_IDR = perform_enrichment_analysis(
#   df = common_Nterm_alphafold_N_terminus,
#   ptm_types = ["common_Nterm"],
#   rois = ["IDR"],
#   ptm_site_dict = Nterm_site_dict,
#   quality_cutoffs = [0]
# )

# high_acc_5
# enrichment_N_terminus_high_acc_5 = perform_enrichment_analysis(
#   df = common_Nterm_alphafold_N_terminus,
#   ptm_types = ["common_Nterm"],
#   rois = ["high_acc_5"],
#   ptm_site_dict = Nterm_site_dict,
#   quality_cutoffs = [0]
# )

# low_acc_5
# enrichment_N_terminus_low_acc_5 = perform_enrichment_analysis(
#   df = common_Nterm_alphafold_N_terminus,
#   ptm_types = ["common_Nterm"],
#   rois = ["low_acc_5"],
#   ptm_site_dict = Nterm_site_dict,
#   quality_cutoffs = [0]
# )

# flexible_pattern
# enrichment_N_terminus_flexible_pattern = perform_enrichment_analysis(
#   df = common_Nterm_alphafold_N_terminus,
#   ptm_types = ["common_Nterm"],
#   rois = ["flexible_pattern"],
#   ptm_site_dict = Nterm_site_dict,
#   quality_cutoffs = [0]
# )

# N-terminus colocalization
# common_Nterm_colocaliztion = evaluate_ptm_colocalization(
#     df = common_Nterm_alphafold_N_terminus,
#     ptm_target = 'self',
#     ptm_types = ['common_Nterm'],
#     ptm_dict = Nterm_site_dict,
#     pae_dir = "data_source/Nterm_structuremap/Nterm_degradation_pae",
#     min_dist = 1,
#     max_dist = 35,
#     dist_step = 5
# )

# N-terminus 3D structures
common_Nterm_proximity = get_proximity_pvals(
    df = common_Nterm_alphafold_N_terminus,
    ptm_types = ['common_Nterm'],
    ptm_site_dict = Nterm_site_dict,
    error_dir = "data_source/Nterm_structuremap/Nterm_degradation_pae",
    per_site_metric = 'mean',
    error_operation = 'plus',
    n_random = 10000,
    random_seed = 44
)
