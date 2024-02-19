#!/usr/bin/env python
import argparse
import os
import re
import shutil
import subprocess
import sys
import csv
import pandas as pd
import numpy as np
from scipy.spatial.distance import euclidean
from Bio import Phylo
from Bio.Phylo.TreeConstruction import DistanceMatrix, DistanceTreeConstructor

def create_directory(directory_name):
    if not os.path.exists(directory_name):
        os.makedirs(directory_name)

def search_strings_in_pfamscan_database(database, strings):
    matching_lines = []
    found = False

    with open(database, 'r') as db:
        for line in db:
            if any(target_string in line for target_string in strings):
                line = re.sub(r' +', '\t', line)
                matching_lines.append(line.strip())

    return matching_lines

def get_query_domain_profiles(files):
    domain_profiles = []
    for file in files:
        with open(file, 'r') as f:
            for line in f:
                # Remove lines starting with '#' and empty lines
                if not line.startswith('#') and line.strip():
                    line = re.sub(r' +', '\t', line)
                    domain_profiles.append(line.strip())

    return domain_profiles

def get_unannotated_proteins(faa_files, pfamscan_files, out_dir):
    output_files = []
    for faa_file, file in zip(faa_files, pfamscan_files):
        faa_headers = set()
        with open(faa_file, 'r') as faa:
            for line in faa:
                if line.startswith('>'):
                    faa_headers.add(line.strip()[1:])  # Remove the ">" character

        found_headers = set()
        for header in faa_headers:
            with open(file, 'r') as pf:
                for line in pf:
                    if not line.startswith('#') and line.strip():
                        if header in line:
                            found_headers.add(header)

        missing_headers = faa_headers - found_headers
        print(f"found {len(missing_headers)} unnanotated proteins in {faa_file}")
        
        if missing_headers:
            # Extract and save unannotated proteins
            for header in missing_headers:
                output_file_path = os.path.join(out_dir, f"{header}.faa")
                output_files.append(output_file_path)
                with open(faa_file, 'r') as faa, open(output_file_path, 'w') as out_file:
                    record_lines = []
                    write_record = False
                    for line in faa:
                        if line.startswith('>'):
                            if record_lines:
                                out_file.write(''.join(record_lines))
                                record_lines = []
                            if line.strip()[1:] == header:
                                write_record = True
                                out_file.write(line)
                            else:
                                write_record = False
                        elif write_record:
                            record_lines.append(line)
                    if record_lines:
                        out_file.write(''.join(record_lines))
    return output_files

def create_temporal_database(original_database, faa_files, temporal_database):
    shutil.copyfile(original_database, temporal_database)

    # Open the new database file and modify headers
    with open(temporal_database, 'r') as td_file:
        lines = td_file.readlines()

    # Open the new database file again for writing modified content
    with open(temporal_database, 'w') as td_file:
        for line in lines:
            if line.startswith('>'):
                # Modify the header here based on your desired format
                header_info = line.split(';')[0].split('=')[1]
                header_parts = header_info.split('_')[1]
                merged_header = f"{line.split('|')[0]}_{header_parts}"
                header = merged_header.strip()[1:]

                # Write the modified header to the new database file
                td_file.write(f">{header}\n")
            else:
                # Write other lines without modification
                td_file.write(line)

    if faa_files:
        # Append the content of each .faa file to the new database file
        with open(temporal_database, 'a') as td:
            for faa_file in faa_files:
                with open(faa_file, 'r') as faa:
                    td.write(faa.read())

def get_s_unnanotated_proteins(faa_database, strings, s_domain_profiles, out_dir):
    output_files = []
    matching_lines = []
    found = False

    with open(faa_database, 'r') as db:
        for line in db:
            if line.startswith('>'):
                header = line.strip()[1:].split(maxsplit=1)[0]
                if any(target_string in header for target_string in strings):
                    matching_lines.append(header)
    
    annotated_proteins = set(column.split('\t')[0] for column in s_domain_profiles)
    unannotated_proteins = set(matching_lines) - annotated_proteins
    print(f"Found {len(unannotated_proteins)} unnanotated proteins")

    if unannotated_proteins:
        print(f"Writting unnanotated proteins to {out_dir}")
        # Extract and save unannotated proteins
        for protein in unannotated_proteins:
            header_info = None
            header_parts = None

            with open(faa_database, 'r') as db:
                record_lines = []
                write_record = False
                output_file_path = None
                out_file = None

                for line in db:
                    if line.startswith('>'):
                        header_info = line.split(';')[0].split('=')[1]
                        header_parts = header_info.split('_')[1]
                        merged_header = f"{line.split('|')[0]}_{header_parts}"
                        header = merged_header.strip()[1:]

                        if line.strip()[1:].split(maxsplit=1)[0] == protein:
                            output_file_name = f"{header}.faa"
                            output_file_path = os.path.join(out_dir, output_file_name)
                            print(f"Creating {output_file_path}")
                            out_file = open(output_file_path, 'w')
                            out_file.write(f">{header}\n")
                            write_record = True
                        else:
                            write_record = False
                    elif write_record:
                        record_lines.append(line)
                        out_file.write(line)

                if write_record and record_lines:
                    out_file.write(''.join(record_lines))

                if out_file:
                    out_file.close()

                if output_file_path:
                    output_files.append(output_file_path)

    return output_files

def run_jackhmmer(query_files, database, out_dir, threads):
    cluster_counter = 1  # Initialize the cluster counter

    for query in query_files:
        query_filename = os.path.splitext(os.path.basename(query))[0].strip()
        cluster_name = f"cluster{cluster_counter}"
        jh_out = os.path.join(out_dir, f"{cluster_name}.out")
        jh_aln = os.path.join(out_dir, f"{cluster_name}.msa")
        jh_tbl = os.path.join(out_dir, f"{cluster_name}.tbl")

        # Check if any .msa file exists in the output directory
        existing_msa_files = [f for f in os.listdir(out_dir) if f.endswith(".msa")]
        if existing_msa_files:
            print(f"Checking existing .msa files...")
            matching_msa = next((msa for msa in existing_msa_files if query_filename in get_msa_content(os.path.join(out_dir, msa))), None)
            if matching_msa:
                print(f"{query} already represented in {matching_msa}. Skipping jackhmmer.")
                continue

        print(f"Running jackhmmer for {query}...\nCreating {jh_aln}")
        subprocess.run(['jackhmmer', '-o', jh_out, '-A', jh_aln, '--tblout', jh_tbl, '--noali', '--cpu', threads, '-N', '10', query, database])
        cluster_counter += 1

def get_msa_content(file_path):
    with open(file_path, 'r') as msa_file:
        return msa_file.read()

def process_intersecting_clusters(out_dir):
    msa_files = [f for f in os.listdir(out_dir) if f.endswith(".msa")]
    clusters = {}
    for file in msa_files:
        cluster = os.path.splitext(file)[0]
        with open(os.path.join(out_dir, file), 'r') as f:
            content = f.read()
            ids = set(re.findall(r"#=GS (\S+)/", content))
        clusters[cluster] = ids
    
    intersecting_clusters = set()
    for cluster_i in clusters.keys():
        ids_i = clusters[cluster_i]
        for cluster_j in [c for c in clusters.keys() if c > cluster_i]:
            ids_j = clusters[cluster_j]
            if any(id_i in ids_j for id_i in ids_i):
                intersecting_clusters.add((cluster_i, cluster_j))
    # Print number of intersecting clusters
    print(f"Number of intersecting clusters: {len(intersecting_clusters)}\n")
    # Print intersecting clusters
    print("Intersecting clusters:")
    for cluster_pair in intersecting_clusters:
        print(f"{cluster_pair[0]} and {cluster_pair[1]}")
    # Initialize sets for representative and dropped clusters
    representative_clusters = set()
    dropped_clusters = set()
    # Compare lengths and choose representatives
    for cluster_i, cluster_j in intersecting_clusters:
        if len(clusters[cluster_i]) == len(clusters[cluster_j]):
            representative_clusters.add(cluster_i)
        else:
            max_len_cluster = max(cluster_i, cluster_j, key=lambda c: len(clusters[c]))
            min_len_cluster = min(cluster_i, cluster_j, key=lambda c: len(clusters[c]))
            representative_clusters.add(max_len_cluster)
            dropped_clusters.add(min_len_cluster)
    # Print representative clusters
    print("\nRepresentative clusters:")
    for cluster in representative_clusters:
        print(cluster)
    # Print dropped clusters
    print("\nDropped clusters:")
    for cluster in dropped_clusters:
        print(cluster)
        # Rename corresponding .msa files to .dropped
        msa_file_path = os.path.join(out_dir, f"{cluster}.msa")
        dropped_file_path = os.path.join(out_dir, f"{cluster}.dropped")
        if os.path.exists(msa_file_path):
            os.rename(msa_file_path, dropped_file_path)

def get_cluster_counts(ids_set, out_dir):
    existing_msa_files = [f for f in os.listdir(out_dir) if f.endswith(".msa")]
    counts_dict = {id.split('|')[0]: {os.path.splitext(cluster)[0]: 0 for cluster in existing_msa_files} for id in ids_set}

    for msa_file in existing_msa_files:
        msa_path = os.path.join(out_dir, msa_file)
        with open(msa_path, 'r') as msa:
            for line in msa:
                if line.startswith('#=GS'):
                    parts = line.split()
                    if len(parts) >= 3:
                        # Extracting ID and cluster from the line
                        msa_id = parts[1].split('_')[0]
                        cluster = os.path.splitext(msa_file)[0]
                        # Updating the counts if id is in ids_set
                        if msa_id in counts_dict:
                            counts_dict[msa_id][cluster] += 1
                            #for id, counts in counts_dict.items():
                            #    print(f"ID: {id}, Counts: {counts}")

    return counts_dict

def rename_ids(counts_dict, ids_mapping):
    renamed_dict = {}
    for key_in_counts, counts in counts_dict.items():
        for original_id, truncated_id in ids_mapping.items():
            if key_in_counts == truncated_id:
                renamed_dict[original_id] = counts
                break
    
    return renamed_dict

def filter_rows_with_zero_sum(input_file_path, output_file_path):
    # Read the file into a DataFrame
    data = pd.read_csv(input_file_path, sep='\t', index_col=0)
    # Filter out rows with a sum of zero
    data_filtered = data[data.sum(axis=1) != 0]
    # Write the modified DataFrame back to a file
    data_filtered.to_csv(output_file_path, sep='\t')

def read_data(input_file):
    data = {}
    ids = []
    with open(input_file, 'r') as file:
        lines = file.readlines()
        for line in lines[1:]:
            fields = line.strip().split('\t')
            id = fields[0]
            scores = list(map(float, fields[1:]))
            data[id] = scores
            ids.append(id)
    return ids, data

def bray_curtis_distance(v1, v2):
    numerator = sum(abs(x - y) for x, y in zip(v1, v2))
    denominator = sum(v1) + sum(v2)
    return numerator / denominator

def calculate_bray_curtis_matrix(ids, data):
    num_samples = len(ids)
    matrix = []
    matrix_sym = np.zeros((num_samples, num_samples))

    for i in range(num_samples):
        row = []
        for j in range(i + 1):
            if i == j:
                row.append(0)  # Diagonal elements are zero
            else:
                distance = bray_curtis_distance(data[ids[i]], data[ids[j]])
                row.append(distance)
        for j in range(i, num_samples): # For symmetric matrix
            distance = bray_curtis_distance(data[ids[i]], data[ids[j]])
            matrix_sym[i][j] = distance
            matrix_sym[j][i] = distance

        matrix.append(row)

    return matrix, matrix_sym

def calculate_euclidean_distance_matrix(ids, data):
    num_samples = len(ids)
    matrix = []
    matrix_sym = np.zeros((num_samples, num_samples))

    for i in range(num_samples):
        row = []
        for j in range(i + 1):
            if i == j:
                row.append(0)  # Diagonal elements are zero
            else:
                distance = euclidean(data[ids[i]], data[ids[j]])
                row.append(distance)
        for j in range(i, num_samples): # For symmetric matrix
            distance = euclidean(data[ids[i]], data[ids[j]])
            matrix_sym[i][j] = distance
            matrix_sym[j][i] = distance

        matrix.append(row)

    return matrix, matrix_sym

def save_matrix_to_csv(matrix, ids, filename):
    with open(filename, 'w', newline='') as csvfile:
        csvwriter = csv.writer(csvfile)
        # Write header (IDs corresponding to columns)
        csvwriter.writerow([''] + ids)
        # Write rows with row IDs and corresponding data
        for i, row_id in enumerate(ids):
            csvwriter.writerow([row_id] + list(matrix[i]))

def reorder_matrix_sym_with_tree(tree, matrix_file_path, output_file_path):
    # Read data from the CSV file
    matrix_data = pd.read_csv(matrix_file_path, index_col=0)
    # Get the leaves of the tree in hierarchical order
    leaves = [leaf.name for leaf in tree.get_terminals()]
    # Reorder rows and columns based on the tree leaves
    reordered_matrix = matrix_data.loc[leaves, leaves]
    # Save the reordered matrix to a new file
    reordered_matrix.to_csv(output_file_path)

def compute_distance_tree(distance_matrix, ids):
    # Create a DistanceMatrix object
    dist_matrix = DistanceMatrix(ids, matrix=distance_matrix)
    #print(dist_matrix)
    # Use the NeighborJoining algorithm
    constructor = DistanceTreeConstructor()
    tree = constructor.nj(dist_matrix)
    tree.root_at_midpoint()

    return tree

def reorder_rows_based_on_tree(tree, data_file_path, output_file_path):
    # Read data from the file
    data = pd.read_csv(data_file_path, sep='\t', index_col=0)
    # Get the leaves of the tree in hierarchical order
    leaves = [leaf.name for leaf in tree.get_terminals()]
    # Reorder rows based on the tree leaves
    reordered_data = data.loc[leaves]
    # Save the reordered data to a new file
    reordered_data.to_csv(output_file_path)


def main():
    parser = argparse.ArgumentParser(description="PanPhylo (Pangenomic and Phylogenomic analysis) is intended to reconstruct the phylogenomic relationships between a group of assembled viruses (metagenomic) and a group of known viruses selected from previous knowledge of taxonomic membership. It takes the same input list.txt as VirTaK.py, aswell as the VirTaK database, and some files created during the taxonomic classification process. It creates a domain (and proteins without domain annotation) content table which is used to calculate a Bray-Curtis distance matrix. Distances are represented by a distance tree built with the neighbor joining algorithm")
    parser.add_argument('-l', '--list', required=True, help='Input list of query .fasta files')
    parser.add_argument('-d', '--database', required=True, help='Path to directory where VirTaK database is located. Must include database name without extension. Example: db/VMR_MSL38_v1_complete_genomes')
    parser.add_argument('-s', '--strings_file', required=False, help='Optional. File containing a list of strings to search in the VirTaK database. For example: Coronaviridae')
    parser.add_argument('-o', '--output_dir', default='PanPhylo', help='Optional. Output directory name')
    parser.add_argument('-t', '--threads', default='2', help='Number of parallel CPU workers to use for HMMER (default=2)')
    args = parser.parse_args()

    output_directory = args.output_dir
    print(f"Creating {output_directory}/ directory")
    create_directory(output_directory)

    list_file = args.list
    print(f"Searching for domain profiles of query viruses specified in {list_file}:")
    with open(list_file, 'r') as lf:
        files = [line.strip().replace('.fasta', '.pfamscan') for line in lf.readlines()]
        search_files = ", ".join(files)
        print(f"\t{search_files}\n")
    q_domain_profiles = get_query_domain_profiles(files)

    hmm_directory = f"{output_directory}/HMMER"
    print(f"Creating {hmm_directory}/ directory")
    create_directory(hmm_directory)

    print(f"Looking for proteins without domain annotation. faa files:")
    with open(list_file, 'r') as lf:
        faa_files = [line.strip().replace('.fasta', '.faa') for line in lf.readlines()]
        search_faa_files = ", ".join(faa_files)
        print(f"\t{search_faa_files}\n")

    faa_database = f"{args.database}.faa"
    unannotated_files = get_unannotated_proteins(faa_files, files, hmm_directory)
    if unannotated_files:
        u_files = ", ".join(unannotated_files)
        print(f"\nFiles created:\n\t{u_files}")
        tmp_db = os.path.join(hmm_directory, "tmp_db.faa")
        print(f"\nCreating temporal database:\n\t{tmp_db}")
        create_temporal_database(faa_database, unannotated_files, tmp_db)

################### If strings file is provided ################### 
    if args.strings_file:
        strings_file = args.strings_file
        print(f"\nStrings file provided: \n\t{strings_file}\n")
        with open(strings_file, 'r') as sf:
            strings = [line.strip() for line in sf.readlines()]
            search_strings = ", ".join(strings)
        pfamscan_database = f"{args.database}.pfamscan"
        print(f"Will search domain profiles of all viruses matching: \n\t{search_strings} \nin: \n\t{pfamscan_database}\n")
        s_domain_profiles = search_strings_in_pfamscan_database(pfamscan_database, strings)
        domain_profiles = q_domain_profiles + s_domain_profiles

        print(f"Searching for unannotatd proteins of viruses matching \n\t{search_strings} \nin: \n\t{faa_database} and {pfamscan_database}\n")
        s_unannotated_files = get_s_unnanotated_proteins(faa_database, strings, s_domain_profiles, hmm_directory)
        if s_unannotated_files:
            tmp_db = os.path.join(hmm_directory, "tmp_db.faa")
            if os.path.exists(tmp_db):
                print(f"Temporal database {tmp_db} already created")
            else:
                print(f"\nCreating temporal database:\n\t{tmp_db}")
                create_temporal_database(faa_database, unannotated_files, tmp_db)
        all_unannotated_files = unannotated_files + s_unannotated_files
    else:
        domain_profiles = q_domain_profiles
        all_unannotated_files = unannotated_files
###################################################################
    
    domain_set = set(column.split('\t')[6] for column in domain_profiles)
    domains = ", ".join(domain_set)
    print(f"\n{len(domain_set)} domains identified in the dataset:\n\n{domains}")
    ids_set = set(column.split('\t')[0].rsplit('_', 1)[0] for column in domain_profiles)

    print(f"\nGetting domain counts for every genome ({len(ids_set)} genomes)\n")
    domain_counts = {id: {domain: 0 for domain in domain_set} for id in ids_set}
    for entry in domain_profiles:
        id, domain = entry.split('\t')[0].rsplit('_', 1)[0], entry.split('\t')[6]
        domain_counts[id][domain] += 1

    #print("Domain Counts:")
    #for id, counts in domain_counts.items():
    #    print(f"ID: {id}, Counts: {counts}")

    if all_unannotated_files:
        run_jackhmmer(all_unannotated_files, tmp_db, hmm_directory, args.threads)
        print("\nSearching for intersecting clusters...")
        process_intersecting_clusters(hmm_directory)
        counts = get_cluster_counts(ids_set, hmm_directory)
        ids_mapping = {original_id: original_id.split('|')[0] for original_id in ids_set}
        cluster_counts = rename_ids(counts, ids_mapping)

        #print("Cluster Counts:")
        #for id, counts in cluster_counts.items():
        #    print(f"ID: {id}, Counts: {counts}")

        merged_counts = {}
        # Iterate over keys in both domain_counts and cluster_counts
        for id_key in set(domain_counts.keys()) | set(cluster_counts.keys()):
            merged_counts[id_key] = {}
            # Update counts for domains
            for domain, count in domain_counts.get(id_key, {}).items():
                merged_counts[id_key][domain] = count
            # Update counts for clusters
            for cluster, count in cluster_counts.get(id_key, {}).items():
                merged_counts[id_key][cluster] = count
        # Print the merged counts
        #print("Merged Counts:")
        #for id, counts in merged_counts.items():
        #    print(f"ID: {id}, Counts: {counts}")
    else:
        merged_counts = domain_counts
    
    print("\nCreating counts table")
    counts_file_path = f"{output_directory}/counts.txt"
    # Open the file for writing
    with open(counts_file_path, 'w', newline='') as csvfile:
        # Create a CSV writer with tab delimiter
        csv_writer = csv.writer(csvfile, delimiter='\t')
        # Write header row with column names (domains/clusters)
        header_row = ['\t'] + list(merged_counts[next(iter(merged_counts))].keys())
        csv_writer.writerow(header_row)
        # Write data rows
        for id, counts in merged_counts.items():
            # Convert counts to list in the order of header
            data_row = [id] + [counts.get(column, 0) for column in header_row[1:]]
            csv_writer.writerow(data_row)
    print(f"Counts table written to {counts_file_path}")

######################### Create results files #########################
    tree_output_path = os.path.join(output_directory, "results_nj.tree")
    output_dist_path = os.path.join(output_directory, "distance_matrix.csv")
    out_reordered_dist = os.path.join(output_directory, "reordered_distance_matrix.csv")
    reordered_counts_path = os.path.join(output_directory, "reordered_counts.csv")

    filter_rows_with_zero_sum(counts_file_path, counts_file_path)
    print(f"Computing Bray-Curtis matrix from {counts_file_path}")
    ids, data = read_data(counts_file_path)
    matrix, matrix_sym = calculate_bray_curtis_matrix(ids, data) 

    save_matrix_to_csv(matrix_sym, ids, output_dist_path)

    # Compute the distance tree
    print(f"Computing tree from distance matrix")
    nj_tree = compute_distance_tree(matrix, ids)
    Phylo.write(nj_tree, tree_output_path, "newick")
    print(f"Tree written to {tree_output_path}")
    
    print("Reordering distance matrix based on tree clusterings")
    reorder_matrix_sym_with_tree(nj_tree, output_dist_path, out_reordered_dist)
    print(f"Reordered distance matrix written to {out_reordered_dist}")

    print("Reordering counts file based on tree clusterings")
    filter_rows_with_zero_sum(counts_file_path, counts_file_path) 
    reorder_rows_based_on_tree(nj_tree, counts_file_path, reordered_counts_path)  
    print(f"Reordered counts file written to {reordered_counts_path}")  

    print("Creating transposed counts file based")
    df = pd.read_csv(reordered_counts_path, index_col=0)
    df_transposed = df.transpose()
    transposed_counts_path = os.path.join(output_directory, "transposed_counts.txt")
    df_transposed.to_csv(transposed_counts_path, sep='\t')
    print(f"Transposed counts file written to {transposed_counts_path}")
    print(f"Reordering transposed counts file rows based on euclidean distance")
    filter_rows_with_zero_sum(transposed_counts_path, transposed_counts_path)
    ids2, data2 = read_data(transposed_counts_path)
    matrix2, matrix_sym2 = calculate_euclidean_distance_matrix(ids2, data2)
    nj_tree2 = compute_distance_tree(matrix2, ids2)
    reordered_t_counts_path = os.path.join(output_directory, "reordered_transposed_counts.csv")
    reorder_rows_based_on_tree(nj_tree2, transposed_counts_path, reordered_t_counts_path)
    print(f"Reordered transposed counts file written to {reordered_t_counts_path}")

if __name__ == "__main__":
    main()
