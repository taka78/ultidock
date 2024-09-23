import sys
import glob
import re
import os
import pandas as pd

files = glob.glob('/mnt/g/docking-outputs-not-wrecked-ep1-part-2/*.pdbqt')
main_df = pd.DataFrame()
print(files)

pdbqt_data_for_this_file = []
name = []
model = []


for file in files:
    with open(file) as fp:
        # read all lines in a list
        lines = fp.readlines()
    for line in lines:
        if line.startswith("REMARK  Name"):
            name_splited = line.split()
            name.append(name_splited)
        elif line.startswith("REMARK VINA RESULT:"):
            line_splited = line.split()
            pdbqt_data_for_this_file.append(line_splited)
        elif line.startswith("MODEL"):
            model_splited = line.split()
            model.append(model_splited)

    pdbqt_df = pd.DataFrame(pdbqt_data_for_this_file, columns=['REMARK', 'VINA', 'RESULT', 'Affinity (kcal/mol)', 'Dist from RMSD l.b.', 'Best Mode RMSD u.b.'])
    name_df = pd.DataFrame(name, columns=['REMARK', 'Name','= because im lazy','ZINC ID'])
    model_df = pd.DataFrame(model, columns=['MODEL','Number'])

    pdbqt_data = {'ZINC_ID': name_df['ZINC ID'], 'MODEL' : model_df['Number'], 'Affinity (kcal/mol)': pdbqt_df['Affinity (kcal/mol)'], 'Dist from RMSD l.b.': pdbqt_df['Dist from RMSD l.b.'], 'Best Mode RMSD u.b.': pdbqt_df['Best Mode RMSD u.b.']}
    pdbqt_data_frame = pd.DataFrame(pdbqt_data)

    print(pdbqt_data_frame)

main_df = main_df._append(pdbqt_data_frame)
print(main_df)


main_df.to_csv("output.csv", index=False)
main_df.to_json("output.json", orient='records')

