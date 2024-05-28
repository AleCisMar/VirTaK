import argparse
import pandas as pd
import numpy as np

# Script used to get VirTaK test statistics
# To get a summary:
# python get_stats.py -f dsDNA_list.txt | cut -d',' -f 1,3,9 | sort -u > dsDNA_short_stats.csv

def process_file(file_path):
    # Read input CSV file
    df = pd.read_csv(file_path, index_col=0)

    # Identify row names without '|'
    rows_without_pipe = [row for row in df.index if '|' not in row]
    rows_with_pipe = [row for row in df.index if '|' in row]

    # Replace values associated with rows without '|'
    df.loc[rows_without_pipe, :] = np.nan
    df.loc[:, rows_without_pipe] = np.nan

    # Extract values below the diagonal
    values_below_diagonal = df.where(np.tril(np.ones(df.shape), k=-1).astype(bool))

    # Create a list of non-NaN values
    non_nan_values = values_below_diagonal.stack().dropna().tolist()

    # Calculate statistics
    min_value = min(non_nan_values)
    q1_value = np.percentile(non_nan_values, 25)
    q2_value = np.percentile(non_nan_values, 50)
    q3_value = np.percentile(non_nan_values, 75)
    max_value = max(non_nan_values)
    mean_value = np.mean(non_nan_values)
    std_dev_value = np.std(non_nan_values)

    # Print statistics in a format suitable for R
    file_name = file_path.split('/')[1]
    for value in non_nan_values:
        print(f"{file_name},{value},{len(rows_with_pipe)},{min_value},{q1_value},{q2_value},{q3_value},{max_value},{mean_value},{std_dev_value}")

def main():
    parser = argparse.ArgumentParser(description="Process files and calculate statistics.")
    parser.add_argument('-f', '--file_list', required=True, help='File containing a list of files to process')
    args = parser.parse_args()

    # Read the list of files from the input file
    with open(args.file_list, 'r') as file:
        files_to_process = file.read().splitlines()

    # Process each file in the list
    print("Family,Distance,N,Min,Q1,Q2,Q3,Max,Mean,SD")
    for file_path in files_to_process:
        process_file(file_path)

if __name__ == "__main__":
    main()
