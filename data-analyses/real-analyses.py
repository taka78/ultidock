import scipy as sp
import numpy as np
import pandas as pd

OUTPUT_PATH = "output.csv"

print(f"Loading: {OUTPUT_PATH}")
data = pd.read_csv(OUTPUT_PATH)

#choosing the affinity values less than -6 kcal/mol and Dist from RMSD l.b. less than 2Å and Best mode RMSD u.b. less than 3Å and creating a new dataframe with them
great_ones = data[(data['Affinity (kcal/mol)'] < -7) & (data['Dist from RMSD l.b.'] < 5) & (data['Best Mode RMSD u.b.'] < 10)]

#choosing the only if a ZINC_ID is duplicated from great_ones and creating a new dataframe with them
fabolous_ones = great_ones[great_ones["ZINC_ID"].duplicated(keep=False)]

we_need_this = fabolous_ones[(fabolous_ones['MODEL'] > 2)]

print(we_need_this)
# Optional: export filtered data
we_need_this.to_csv("selected_candidates.csv", index=False)
