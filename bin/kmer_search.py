#!/usr/bin/env python
import os
import argparse
import itertools
from collections import defaultdict
from scipy.spatial import distance

def generate_kmers(k):
    # Define the letters
    letters = ['A', 'T', 'G', 'C']

    # Generate all possible k-mers
    all_kmers = [''.join(p) for p in itertools.product(letters, repeat=k)]

    return all_kmers

def get_kmer_counts(sequence, kmers):
    sequence_length = len(sequence)
    kmer_counts = {kmer: 0 for kmer in kmers}

    # Count the occurrences of each k-mer in the sequence
    for i in range(sequence_length - len(kmers[0]) + 1):
        kmer = sequence[i:i+len(kmers[0])]
        if kmer in kmers:
            kmer_counts[kmer] += 1

    return kmer_counts

def calculate_bray_curtis(profile1, profile2):
    numerator = sum(min(profile1[k], profile2[k]) for k in profile1.keys())
    denominator1 = sum(profile1[k] for k in profile1.keys())
    denominator2 = sum(profile2[k] for k in profile2.keys())
    
    return 1.0 - (2.0 * numerator / (denominator1 + denominator2))

def main():
    parser = argparse.ArgumentParser(description="Calculates k-mer profiles for every file.fasta file in the current directory (creates file.kmer files) and searches for the n (default=10) viruses with the least dissimilar (Bray-Curtis) k-mer profiles in the input database. The results are written to files named file_matches.txt with 'profile names' in the first column, and 'dissimilarity values' in the second column")
    parser.add_argument('-d', '--database', required=True, help='Input k-mer profile database. Example: db/VMR_MSL38_v1_complete_genomes.kmers')
    parser.add_argument('-k', '--kmer_length', type=int, default=4, help='OPTIONAL: Length of k-mers. Default = 4')
    parser.add_argument('-n', '--n_matches', type=int, default=10, help='OPTIONAL: Number of top matches to display in output file. Default = 10')
    args = parser.parse_args()

    # Read the input k-mer profile database
    print(f"Reading input k-mer profiles database: {args.database}")
    kmer_profiles = {}
    with open(args.database, 'r') as db_file:
        header = db_file.readline().strip().split('\t')
        kmers = header[1:]
        for line in db_file:
            parts = line.strip().split('\t')
            profile_name = parts[0]
            profile_values = list(map(int, parts[1:]))
            kmer_profiles[profile_name] = dict(zip(kmers, profile_values))

    # Generate all possible k-mers
    kmers = generate_kmers(args.kmer_length)

    # Calculate k-mer profiles for all .fasta files in the current directory    
    for file_name in os.listdir('.'):
        if file_name.endswith('.fasta'):
            print(f"Creating {args.kmer_length}-mer profile for {file_name}")
            kmer_profile = defaultdict(int)
            with open(file_name, 'r') as fasta_file:
                for line in fasta_file:
                    if line.startswith('>'):
                        continue
                    sequence = line.strip()
                    kmer_counts = get_kmer_counts(sequence, kmers)
                    for kmer, count in kmer_counts.items():
                        kmer_profile[kmer] += count

            # Save the new k-mer profile to a .kmer file with the same name as the input fasta file
            kmer_file_name = os.path.splitext(file_name)[0] + '.kmer'
            with open(kmer_file_name, 'w') as kmer_file:
                kmer_file.write(f"{file_name}\t")
                kmer_file.write('\t'.join(str(kmer_profile[kmer]) for kmer in kmers))
                kmer_file.write('\n')
            print(f"{file_name} {args.kmer_length}-mer profile written to {kmer_file_name}")

            # Calculate Bray-Curtis dissimilarity to every k-mer profile in the database
            print(f"Comparing {kmer_file_name} profile with profiles in {args.database}")
            n = args.n_matches
            top_n_matches = [(None, float('inf'))] * n

            for profile_name, profile_values in kmer_profiles.items():
                dissimilarity = calculate_bray_curtis(kmer_profile, profile_values)

                for i, (match_name, match_value) in enumerate(top_n_matches):
                    if dissimilarity < match_value:
                        top_n_matches.insert(i, (profile_name, dissimilarity))
                        break
            top_n_matches = top_n_matches[:n]

            result_file_name = os.path.splitext(file_name)[0] + '_matches.txt'
            with open(result_file_name, 'w') as output_file:
                output_file.write(f"Top {n} matches for {kmer_file_name}:\n")
                print(f"Top {n} matches for {kmer_file_name}:")
                print("Profile Name\tDissimilarity Value")
                for match in top_n_matches:
                    profile_name, dissimilarity = match                    
                    output_file.write(f"{profile_name}\t{dissimilarity}\n")
                    print(f"{profile_name}\t{dissimilarity}")

            print(f"Results written to {result_file_name}")

if __name__ == "__main__":
    main()
