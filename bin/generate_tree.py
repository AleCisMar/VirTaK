#!/usr/bin/env python
import argparse
import re
import os
import pandas as pd
from scipy.spatial import distance
from scipy.cluster import hierarchy
from scipy.cluster.hierarchy import to_tree
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

def bray_curtis_distance(v1, v2):
    numerator = sum(abs(x - y) for x, y in zip(v1, v2))
    denominator = sum(v1) + sum(v2)
    return numerator / denominator

def calculate_bray_curtis_matrix(ids, data):
    num_samples = len(ids)
    matrix = np.zeros((num_samples, num_samples))

    for i in range(num_samples):
        for j in range(i, num_samples):
            distance_value = bray_curtis_distance(data[i], data[j])
            matrix[i][j] = distance_value
            matrix[j][i] = distance_value  # Matrix is symmetric, so set the mirrored value

    return matrix

def save_matrix(output_file, ids, matrix):
    with open(output_file, 'w') as file:
        file.write('\t' + '\t'.join(ids) + '\n')
        for i, row in enumerate(matrix):
            file.write(ids[i] + '\t' + '\t'.join(map(str, row)) + '\n')

def compute_distance_tree(distance_matrix):
    condensed_distance = distance.squareform(distance_matrix.values)
    linkage = hierarchy.linkage(condensed_distance, method='average')
    dendrogram = hierarchy.dendrogram(linkage, labels=distance_matrix.index, orientation='right', leaf_font_size=8)
    return linkage, dendrogram

def generate_newick_from_tree(tree, leaf_names):
    if tree.is_leaf():
        return leaf_names[tree.id]
    else:
        left_node = generate_newick_from_tree(tree.left, leaf_names)
        right_node = generate_newick_from_tree(tree.right, leaf_names)
        branch_length = tree.dist
        return f"({left_node}:{branch_length},{right_node}:{branch_length})"

def main():
    parser = argparse.ArgumentParser(description="Merge input kmer profiles with specified kmer profiles in the kmer database, calculate Bray-Curtis dissimilarity matrix, and compute a distance tree")
    parser.add_argument('-d', '--database', required=True, help='Kmer profile database')
    parser.add_argument('-l', '--list', required=True, help='File listing input ".kmer" files to process')
    parser.add_argument('-s', '--string', required=True, help='String to search in the kmer database. For example: Coronaviridae')
    parser.add_argument('-n', '--newick', required=True, help='Output Newick file')
    args = parser.parse_args()

    merged_profiles_output = f"{args.string}_merged.kmers"
    matrix_output = f"{args.string}_merged.mat"

    # Step 1: Merge kmer profiles
    kmer_profiles = read_kmer_database(args.database)
    kmer_files = [line.strip() for line in open(args.list) if line.strip().endswith(".kmer")]
    print(f"Merging {kmer_files} profiles with profiles in {args.database} whose header match the word {args.string}")
    merged_profiles = merge_profiles(args.string, kmer_files, kmer_profiles)
    save_merged_profiles_to_file(merged_profiles, merged_profiles_output)
    print(f"Merged profiles written to {merged_profiles_output}")

    # Step 2: Calculate Bray-Curtis dissimilarity matrix
    print(f"Calculating Bray-Curtis dissimilarity matrix")
    ids = list(merged_profiles.keys())
    data = [merged_profiles[profile] for profile in ids]
    matrix = calculate_bray_curtis_matrix(ids, data)
    save_matrix(matrix_output, ids, matrix)
    print(f"Dissimilarity matrix written to {matrix_output}")
    distance_matrix = pd.DataFrame(matrix, index=ids, columns=ids)
    leaf_names = list(distance_matrix.index)

    # Step 2: Compute the distance tree
    print(f"Computing distance tree")
    linkage, dendrogram = compute_distance_tree(distance_matrix)

    # Step 3: Convert linkage matrix to a tree structure
    tree = to_tree(linkage)

    # Step 4: Generate Newick format from the tree structure
    newick_tree = generate_newick_from_tree(tree, leaf_names)

    # Step 5: Save the Newick tree to the output file
    newick_otuput = args.newick
    with open(newick_otuput, "w") as f:
        f.write(newick_tree + ";")
    print(f"Tree written to {newick_otuput}")

if __name__ == "__main__":
    main()
