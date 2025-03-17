# import packages
import pandas as pd
from tqdm import tqdm

# import localCIDER package (https://pappulab.github.io/localCIDER/, Version 0.1.20)
# using for calculate various properties of a peptide sequence
from localcider.sequenceParameters import SequenceParameters

# import Nterm 13mer sequence
HEK_Nterm_Kd_half_life_sequence = pd.read_csv(
  'data_source/Nterm_sequence/HEK_Nterm_Kd_half_life_sequence.csv',
  header = 0,
  index_col = None
)

# get sequence parameters for Nterm_13mer
HEK_Nterm_Kd_half_life_sequence.loc[:, "Nterm_13mers_sequence_parameter"]= HEK_Nterm_Kd_half_life_sequence["Nterm_13mer"].astype(pd.StringDtype()).apply(SequenceParameters)

## calculate hydropathy, FCR, NCPR, isoelectric point and kappa values
methods = {
    "hydropathy": "get_mean_hydropathy",
    "FCR": "get_FCR",
    "NCPR": "get_NCPR",
    "isoelectric_point": "get_isoelectric_point",
    "kappa": "get_kappa"
}

for col, method in tqdm(methods.items()):
    HEK_Nterm_Kd_half_life_sequence.loc[:, col] = HEK_Nterm_Kd_half_life_sequence["Nterm_13mers_sequence_parameter"].apply(
        lambda obj: getattr(obj, method)() if obj is not None else None
    )

# export result
HEK_Nterm_Kd_half_life_sequence.to_csv(
  "data_source/Nterm_13mer_sequence_features/HEK_Nterm_13mer_Kd_half_life_sequence_features.csv", 
  index = False
)
