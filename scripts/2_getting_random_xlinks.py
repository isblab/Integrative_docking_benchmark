import pandas as pd

def read_and_choose_lines(csv_file, output_file, num_lines=6):
    df = pd.read_csv(csv_file)
    chosen_rows = df.sample(n=num_lines)
    chosen_rows.to_csv(output_file, index=False)

read_and_choose_lines(sys.argv[1],sys.argv[2],int(sys.argv[3]))
