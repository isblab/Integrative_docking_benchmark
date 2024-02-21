import pandas as pd
import os
import sys
from Bio import PDB

def calculate_distance(res1, prot1, res2, prot2, structure):
    model = structure[0]
    chain1 = model[prot1]
    chain2 = model[prot2]

    ca1 = chain1[res1]['CA']
    ca2 = chain2[res2]['CA']
    distance = ca1 - ca2

    return distance

def check_distances(selected_rows, pdb_file):
    parser = PDB.PDBParser(QUIET=True)
    structure = parser.get_structure('protein', pdb_file)

    for idx, row in selected_rows.iterrows():
        distance = calculate_distance(row['res1'], row['prot1'], row['res2'], row['prot2'], structure)

    return distance

if __name__ == "__main__":
    input_file = sys.argv[1]
    num_xlinks = int(sys.argv[2])
    xlinker = sys.argv[3]
    pdb_file = sys.argv[4]

    threshold = 32 if xlinker == 'DSSO' else 20

    output_directory = os.path.dirname(input_file)
    output_file_path = os.path.join(output_directory, f"{os.path.splitext(os.path.basename(input_file))[0].split('_')[0]}_{xlinker}_{num_xlinks}_filtered.csv")

    df = pd.read_csv(input_file)

    # Consider the files which have less crosslinks, don't need to "sample"
    if num_xlinks >= len(df):
        df.to_csv(output_file_path, index=False)

    else:
        filtered_df = pd.DataFrame()
        while len(filtered_df) < num_xlinks:
            selected_rows = df.sample(1)
            distance = check_distances(selected_rows, pdb_file)

            # Check for duplicate crosslinks and crosslinks which have ca distance more than given cutoff in the pdb file
            if not (filtered_df.isin(selected_rows.values[0]).all(axis=1).any() or distance > threshold):
                filtered_df = pd.concat([filtered_df, selected_rows], ignore_index=True)
            else:
                print('Sample again')
        filtered_df.to_csv(output_file_path, index=False)
