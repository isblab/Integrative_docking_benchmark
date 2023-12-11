import glob
import os
import sys
import pandas as pd

input_folder = sys.argv[1]

for file in glob.glob(os.path.join(input_folder, '*.txt')):
    df = pd.read_csv(file, delim_whitespace=True, skiprows=1, header=None)
    df.columns = ["Index", "Model", "Atom1", "Atom2", "SASD"]

    filtered_df = df[df["Atom1"].str.split('-').str[2] != df["Atom2"].str.split('-').str[2]]
    df2 = pd.concat([
        filtered_df["Atom1"].str.split('-').str[1].rename('res1'),
        filtered_df["Atom1"].str.split('-').str[2].rename('prot1'),
        filtered_df["Atom2"].str.split('-').str[1].rename('res2'),
        filtered_df["Atom2"].str.split('-').str[2].rename('prot2'),
    ], axis=1)

    output_filename = os.path.splitext(os.path.basename(file))[0] + '_processed_crosslinks.csv'
    output_path = os.path.join(input_folder, output_filename)

    df2.to_csv(output_path, index=False)
    print(df2)
    df2.to_csv(f'{file}_crosslinks.csv', index = False)
