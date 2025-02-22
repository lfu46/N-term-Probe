import pandas as pd
import numpy as np

# import strucutremap functions as instructed by tutorial (https://github.com/MannLabs/structuremap/blob/main/nbs/tutorial.ipynb)
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

Nterm_protein_list = HEK_Nterm_Kd_half_life_sequence["UniProt_Accession"].unique().tolist()

## download AlphaFold data
# crystallographic information file
# valid_proteins_cif, invalid_proteins_cif, existing_proteins_cif = download_alphafold_cif(
#     proteins = Nterm_protein_list,
#     out_folder = "data_source/Nterm_structuremap/Nterm_degradation_cif"
# )

# predicted aligned error
# valid_proteins_pae, invalid_proteins_pae, existing_proteins_pae = download_alphafold_pae(
#     proteins = Nterm_protein_list,
#     out_folder = "data_source/Nterm_structuremap/Nterm_degradation_pae"
# )

# format AlphaFold data input
alphafold_annotation = format_alphafold_data(
  directory = "data_source/Nterm_structuremap/Nterm_degradation_cif",
  protein_ids = Nterm_protein_list
)

