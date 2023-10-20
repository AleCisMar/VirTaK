#!/usr/bin/env python
import argparse
from scipy.spatial import distance
import os
import tempfile
import shutil

def main():
    parser = argparse.ArgumentParser(description="Knowing the taxa associated to the least dissimilar viruses (in kmer_search.py output files), the user can calculate the threshold for a specific taxon (for example Coronaviridae). This code will search for the kmer profiles that match the user defined string (Coronaviridae, Nidovirales or whatever) in the kmer profile database, calculate all vs all Bray-Curtis dissimilarities and get the maximum value (threshold) for that specific taxon. Finally it will update the input file adding the columns 'taxon, threshold, and number of viruses' to the lines that match the string in the input file")
    parser.add_argument('-i', '--input', required=True, help='Input file. Is the output of kmer_search.py. Example: file_matches.txt')
    parser.add_argument('-d', '--database', required=True, help='Kmer profile database. Example: db/VMR_MSL38_v1_complete_genomes.kmers')
    parser.add_argument('-s', '--string', required=True, help='String to search in input file. Example: Coronaviridae, or Nidovirales, etc.')
    args = parser.parse_args()

    matching_lines = []
    with open(args.database, 'r') as input_database:
        print(f"Searching for lines in {args.database} that match {args.string}")
        for line in input_database:
            if args.string in line:
                matching_lines.append(line)

    num_matching_lines = len(matching_lines)
    print(f"Found {num_matching_lines} lines matching {args.string}")

    max_dissimilarity = 0
    for i, line1 in enumerate(matching_lines):
        for j, line2 in enumerate(matching_lines):
            if i != j:
                kmer_profile1 = [float(value) for value in line1.split('\t')[1:]]
                kmer_profile2 = [float(value) for value in line2.split('\t')[1:]]
                dissimilarity = distance.braycurtis(kmer_profile1, kmer_profile2)
                max_dissimilarity = max(max_dissimilarity, dissimilarity)
    print(f"Threshold for {args.string} = {max_dissimilarity}")

    with tempfile.NamedTemporaryFile(mode='w', delete=False) as temp_file:
        with open(args.input, 'r') as input_file:
            for line in input_file:
                if "Profile" in line:
                    header = line.strip() + '\tTaxon\tThreshold\tNumber of Viruses'
                elif args.string in line:
                    new_columns = f'{args.string}\t{max_dissimilarity}\t{num_matching_lines}'
                    updated_line = line.strip() + '\t' + new_columns
                    temp_file.write(updated_line + '\n')
                else:
                    temp_file.write(line)

    shutil.move(temp_file.name, args.input)

if __name__ == "__main__":
    main()