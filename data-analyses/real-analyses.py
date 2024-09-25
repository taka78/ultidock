import scipy as sp
import numpy as np
import pandas as pd

 ###those ara basically output.csv files but because i had so much data, i needed to do it in 2 parts.
data_part1 = pd.read_csv("docking-outputs-not-wrecked-ep1-part-1.csv")
data_part2 = pd.read_csv("docking-outputs-not-wrecked-ep1-part-2.csv")


data = pd.concat([data_part1, data_part2], ignore_index=True)

#choosing the affinity values less than -6 kcal/mol and Dist from RMSD l.b. less than 2Å and Best mode RMSD u.b. less than 3Å and creating a new dataframe with them
great_ones = data[(data['Affinity (kcal/mol)'] < -8) & (data['Dist from RMSD l.b.'] < 5) & (data['Best Mode RMSD u.b.'] < 10)]

#choosing the only if a ZINC_ID is duplicated from great_ones and creating a new dataframe with them
fabolous_ones = great_ones[great_ones["ZINC_ID"].duplicated()]

we_need_this = fabolous_ones[(fabolous_ones['MODEL'] > 2)]


print(we_need_this)


