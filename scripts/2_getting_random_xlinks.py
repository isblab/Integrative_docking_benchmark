import pandas as pd
import os, sys

input_file = sys.argv[1]
num_xlinks = int(sys.argv[2])
xlinker = sys.argv[3]

output_directory = os.path.dirname(input_file)
output_file_path = os.path.join(output_directory, f"{os.path.splitext(os.path.basename(input_file))[0].split('_')[0]}_{xlinker}_{num_xlinks}.csv")
df = pd.read_csv(input_file)

selected_rows = df.sample(num_xlinks)
selected_rows.to_csv(output_file_path, index=False)
