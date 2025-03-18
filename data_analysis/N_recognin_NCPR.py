# import packages
import pandas as pd
from tqdm import tqdm

# import localCIDER package (https://pappulab.github.io/localCIDER/, Version 0.1.20)
from localcider.sequenceParameters import SequenceParameters

# import N-recognin binding region sequence
N_recognin_binding_region_sequence = pd.read_csv(
  "data_source/ELM_degron/N_recognin_binding_region_sequence.csv"
)

# get sequence parameters for binding region sequence
N_recognin_binding_region_sequence.loc[:, "sequence_parameter"] = N_recognin_binding_region_sequence["binding_region"].astype(pd.StringDtype()).apply(SequenceParameters)

# Get all unique UniProt_Accession IDs
protein_ids = N_recognin_binding_region_sequence["UniProt_Accession"].unique()

# Dictionary to store results
NCPR_result = {}

for protein_id in tqdm(protein_ids):
    # Filter the data for the current protein ID
    filtered_data = N_recognin_binding_region_sequence[N_recognin_binding_region_sequence["UniProt_Accession"] == protein_id]

    # Compute linear NCPR for the filtered sequences
    NCPR_result[protein_id] = filtered_data["sequence_parameter"].apply(
      lambda obj: obj.get_linear_NCPR() if obj is not None else None
    )

    
