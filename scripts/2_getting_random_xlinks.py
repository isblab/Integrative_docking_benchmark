import pandas as pd
import os, sys

input_file = sys.argv[1]
num_xlinks = int(sys.argv[2])

# Get the directory of the input file
output_directory = os.path.dirname(input_file)

# Create the output file path
output_file_path = os.path.join(output_directory, f"{os.path.splitext(os.path.basename(input_file))[0].split('_')[0]}_{num_xlinks}.csv")
print(output_file_path)

# Read the input CSV file
df = pd.read_csv(input_file)

# Sample rows and save to the output CSV file
chosen_rows = df.sample(num_xlinks)
chosen_rows.to_csv(output_file_path, index=False)

print(f"File saved successfully to: {output_file_path}")
