import argparse
import re
import os
from collections import defaultdict
from itertools import combinations
from scipy.spatial import distance
from scipy.cluster import hierarchy
from Bio import Phylo
import numpy as np

def read_kmer_database(database_file):
    kmer_profiles = {}
    with open(database_file, 'r') as db_file:
        next(db_file)  # Skip the header line
        for line in db_file:
            parts = line.strip().split('\t')
            profile_name = parts[0]
            profile_values = list(map(int, parts[1:]))
            kmer_profiles[profile_name] = profile_values
    return kmer_profiles

def merge_profiles(query_string, kmer_files, kmer_profiles):
    merged_profiles = {}
    # Find lines (profiles) matching the input string in the kmer database
    matching_profiles = [profile for profile in kmer_profiles.keys() if query_string in profile]
    for profile in matching_profiles:
        merged_profiles[profile] = kmer_profiles[profile]
    # Merge kmer profiles from .kmer files
    for kmer_file in kmer_files:
        file_name = os.path.splitext(os.path.basename(kmer_file))[0]
        with open(kmer_file, 'r') as file:
            for line in file:
                parts = line.strip().split('\t')
                profile_name = parts[0]
                profile_values = list(map(int, parts[1:]))
                if profile_name not in merged_profiles:
                    merged_profiles[profile_name] = []
                merged_profiles[profile_name].extend(profile_values)
    return merged_profiles

def save_merged_profiles_to_file(merged_profiles, output_file):
    with open(output_file, "w") as f:
        for profile, counts in merged_profiles.items():
            f.write(f"{profile}\t" + '\t'.join(map(str, counts)) + "\n")

def main():
    parser = argparse.ArgumentParser(description="Merge kmer profiles from database and .kmer files and save to an output file")
    parser.add_argument('-d', '--database', required=True, help='Kmer profile database')
    parser.add_argument('-l', '--list', required=True, help='List of .kmer files to process')
    parser.add_argument('-s', '--string', required=True, help='String to search in the kmer database')
    parser.add_argument('-o', '--output', required=True, help='Output file for merged profiles')
    args = parser.parse_args()

    # Step 1: Read the input kmer profile database
    kmer_profiles = read_kmer_database(args.database)

    # Step 2: Find lines (profiles) matching the input string in kmer database
    print(f"Matching string '{args.string}' in the kmer database...")
    matching_profiles = [profile for profile in kmer_profiles.keys() if args.string in profile]
    #print(f"Matching profiles: {matching_profiles}")

    # Step 3: Merge kmer profiles from .kmer files and those containing the input string
    kmer_files = [line.strip() for line in open(args.list) if line.strip().endswith(".kmer")]
    merged_profiles = merge_profiles(args.string, kmer_files, kmer_profiles)
    #print("Merged profiles:")
    #for profile, counts in merged_profiles.items():
    #    print(f"{profile}: {counts}")

    # Step 4: Save the merged profiles to the output file
    print(f"Saving merged profiles to {args.output}...")
    save_merged_profiles_to_file(merged_profiles, args.output)
    print(f"Merged profiles saved to {args.output}")

if __name__ == "__main__":
    main()
