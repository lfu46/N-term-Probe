import pandas as pd
import numpy as np

# import strucutremap functions as instructed by tutorial 
# (https://github.com/MannLabs/structuremap/blob/main/nbs/tutorial.ipynb)
# original paper: The structural context of posttranslational modifications at a proteome-wide scale
# (https://journals.plos.org/plosbiology/article?id=10.1371/journal.pbio.3001636#sec009)
import structuremap.utils
structuremap.utils.set_logger()
from structuremap.processing import download_alphafold_cif, download_alphafold_pae, format_alphafold_data, annotate_accessibility, get_smooth_score, annotate_proteins_with_idr_pattern, get_extended_flexible_pattern, get_proximity_pvals, perform_enrichment_analysis, perform_enrichment_analysis_per_protein, evaluate_ptm_colocalization, extract_motifs_in_proteome
from structuremap.plotting import plot_enrichment, plot_ptm_colocalization

# import the Nterm protein data from results
HEK_Nterm_Kd_half_life_sequence = pd.read_csv(
  'data_source/Nterm_sequence/HEK_Nterm_Kd_half_life_sequence.csv',
  header = 0,
  index_col = None
)

Nterm_degradation_protein_list = HEK_Nterm_Kd_half_life_sequence["UniProt_Accession"].unique().tolist()

# download AlphaFold data
crystallographic information file
valid_proteins_cif, invalid_proteins_cif, existing_proteins_cif = download_alphafold_cif(
    proteins = Nterm_degradation_protein_list,
    out_folder = "data_source/Nterm_structuremap/Nterm_degradation_cif"
)

predicted aligned error
valid_proteins_pae, invalid_proteins_pae, existing_proteins_pae = download_alphafold_pae(
    proteins = Nterm_degradation_protein_list,
    out_folder = "data_source/Nterm_structuremap/Nterm_degradation_pae"
)

# format AlphaFold data input
Nterm_degradation_alphafold_annotation = format_alphafold_data(
  directory = "data_source/Nterm_structuremap/Nterm_degradation_cif",
  protein_ids = Nterm_degradation_protein_list
)

# annotate prediction-aware part-sphere exposure (pPSE) values
Nterm_degradation_full_sphere_exposure = annotate_accessibility(
    df = Nterm_degradation_alphafold_annotation, 
    max_dist = 24, 
    max_angle = 180, 
    error_dir = "data_source/Nterm_structuremap/Nterm_degradation_pae"
)

Nterm_degradation_alphafold_accessibility = Nterm_degradation_alphafold_annotation.merge(
  Nterm_degradation_full_sphere_exposure,
  how = "left",
  on = ["protein_id", "AA", "position"]
)

Nterm_degradation_part_sphere_exposure = annotate_accessibility(
    df = Nterm_degradation_alphafold_annotation, 
    max_dist = 12, 
    max_angle = 70, 
    error_dir = "data_source/Nterm_structuremap/Nterm_degradation_pae"
)

Nterm_degradation_alphafold_accessibility = Nterm_degradation_alphafold_accessibility.merge(
  Nterm_degradation_part_sphere_exposure,
  how = "left",
  on = ["protein_id", "AA", "position"]
)

Nterm_degradation_alphafold_accessibility["accessibility"] = np.where(
  Nterm_degradation_alphafold_accessibility.nAA_12_70_pae <= 5, "high exposure", "low exposure"
)

# annotate intrisicly disorder region (IDR)
Nterm_degradation_alphafold_accessibility_smooth = get_smooth_score(
  Nterm_degradation_alphafold_accessibility,
  np.array(["nAA_24_180_pae"]),
  [10]
).reset_index(drop=True)

Nterm_degradation_alphafold_accessibility_smooth["IDR"] = np.where(
  Nterm_degradation_alphafold_accessibility_smooth.nAA_24_180_pae_smooth10 <= 34.27, "IDR", "Structural Region"
)
