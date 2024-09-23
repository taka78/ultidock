import sys
import glob
import re
import os
import pandas as pd

files = glob.glob('/mnt/g/docking-outputs/*.pdbqt')
main_df = pd.DataFrame()
print(files)

target_zinc_ids = ["ZINC000924554774", "ZINC000870235067", "ZINC000575957385", "ZINC000811426569","ZINC001168733702","ZINC000488675918","ZINC000286190168"] ##for nerds
target_zinc_ids_2 = ["ZINC000004471415", "ZINC000000808027", "ZINC000000134139", "ZINC000043591551"] ##for geeks
pdbqt_data_for_this_file = []
name = []
model = []
name_df = pd.DataFrame(columns=["REMARK", "Name", "ZINC ID", "File"])  # Add a "File" column


for file in files:
    print(file)
    with open(file) as fp:
        lines = fp.readlines()
    processed_zinc_ids = []
    for line in lines:
        if line.startswith('REMARK  Name'):
            name_splited = line.split()
            name = name_splited[3]
            zinc_id = name
            if zinc_id not in processed_zinc_ids:  # Check if ZINC ID is already processed
                processed_zinc_ids.append(zinc_id)  # Add to processed set
                name_df = name_df._append({"ZINC ID": name, "File": file}, ignore_index=True)
            
    pdbqt_data_for_this_file = pd.DataFrame(name_df)

main_df = main_df._append(pdbqt_data_for_this_file)
print(main_df)
filtered_df_1 = name_df[name_df["ZINC ID"].isin(target_zinc_ids)]
print("Files containing target nerds:")
for index, row in filtered_df_1.iterrows():
    print(row["File"])
filtered_df_1.to_csv("nerds.csv", index=False)

filtered_df_2 = name_df[name_df["ZINC ID"].isin(target_zinc_ids_2)]
print("Files containing target geeks:")
for index, row in filtered_df_2.iterrows():
    print(row["File"])
filtered_df_2.to_csv("geeks.csv", index=False)