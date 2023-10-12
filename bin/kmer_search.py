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
    parser = argparse.ArgumentParser(description="Calculate k-mer profiles for every .fasta file in the current directory and search for the virus with the least dissimilar (Bray-Curtis) k-mer profile in the input database")
    parser.add_argument('-d', '--database', required=True, help='Input k-mer profile database')
    parser.add_argument('-k', '--kmer_length', type=int, default=4, help='OPTIONAL: Length of k-mers. Default = 4')
    parser.add_argument('-o', '--output', required=True, help='Output file for results')
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
            best_match_name = None
            best_match_value = float('inf')
            for profile_name, profile_values in kmer_profiles.items():                
                dissimilarity = calculate_bray_curtis(kmer_profile, profile_values)
                if dissimilarity < best_match_value:
                    best_match_name = profile_name
                    best_match_value = dissimilarity

            # Write results to the output file
            with open(args.output, 'a') as output_file:
                output_file.write(f"{kmer_file_name}\t{best_match_name}\t{best_match_value}\n")
            print(f"Best match for {kmer_file_name}: {best_match_name} with a distance of {best_match_value}")
    print(f"Results written to {args.output}")

if __name__ == "__main__":
    main()
