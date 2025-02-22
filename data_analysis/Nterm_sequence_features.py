import pandas as pd

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
# hydropathy
HEK_Nterm_Kd_half_life_sequence.loc[:, "hydropathy"] = HEK_Nterm_Kd_half_life_sequence["Nterm_13mers_sequence_parameter"].apply(lambda obj: obj.get_mean_hydropathy() if obj is not None else None)

# fraction of charged residue (FCR)
HEK_Nterm_Kd_half_life_sequence.loc[:, "FCR"] = HEK_Nterm_Kd_half_life_sequence["Nterm_13mers_sequence_parameter"].apply(lambda obj: obj.get_FCR() if obj is not None else None)

# net charge per residue (NCPR)
HEK_Nterm_Kd_half_life_sequence.loc[:, "NCPR"] = HEK_Nterm_Kd_half_life_sequence["Nterm_13mers_sequence_parameter"].apply(lambda obj: obj.get_NCPR() if obj is not None else None)

# isoelectric point
HEK_Nterm_Kd_half_life_sequence.loc[:, "isoelectric_point"] = HEK_Nterm_Kd_half_life_sequence["Nterm_13mers_sequence_parameter"].apply(lambda obj: obj.get_isoelectric_point() if obj is not None else None)

# kappa
HEK_Nterm_Kd_half_life_sequence.loc[:, "kappa"] = HEK_Nterm_Kd_half_life_sequence["Nterm_13mers_sequence_parameter"].apply(lambda obj: obj.get_kappa() if obj is not None else None)

# export result
HEK_Nterm_Kd_half_life_sequence.to_csv(
  "data_source/Nterm_13mer_sequence_features/HEK_Nterm_Kd_half_life_sequence_features.csv", 
  index = False
)
