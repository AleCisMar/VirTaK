#!/usr/bin/env python
import os
import argparse
import itertools
from collections import defaultdict, Counter
from scipy.spatial import distance
import subprocess
import re

def generate_kmers(k):
    # Define the letters
    letters = ['A', 'T', 'G', 'C']

    # Generate all possible k-mers
    all_kmers = [''.join(p) for p in itertools.product(letters, repeat=k)]

    return all_kmers

def get_kmer_counts(sequence, kmers):
    sequence = sequence.upper()
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

def run_prodigal(file_name):
    prodigal_out = "prodigal.out"
    with open(prodigal_out, 'w') as prodigal_out_file:
        fasta_file = file_name
        faa_file = os.path.splitext(file_name)[0] + '.faa'
        subprocess.run(['prodigal', '-i', fasta_file, '-a', faa_file, '-p', 'meta', '-q'], stdout=prodigal_out_file)
        modify_faa_header(faa_file)

def modify_faa_header(faa_file):
    with open(faa_file, 'r') as file:
        lines = file.readlines()
    header_count = 0
    modified_lines = []
    for line in lines:
        if line.startswith(">"):
            header_count += 1
            heder_fields = line.strip().split("_")
            header = heder_fields[0].lstrip(">")

            if "|" in header:
                accession = header.split("|")[0]
                modified_lines.append(f">{accession}_{header_count}\n")
            else:
                modified_lines.append(f">{header}_{header_count}\n")
        else:
            modified_lines.append(line)
    with open(faa_file, 'w') as file:
        file.writelines(modified_lines)

def run_pfamscan(in_file, database):
    out_file = os.path.splitext(in_file)[0] + '.pfamscan'
    subprocess.run(['pfam_scan.pl', '-fasta', in_file, '-dir', database, '-outfile', out_file])

def get_pfams(file_path):
    pfams = set()

    with open(file_path, 'r') as file:
        skip_lines = True
        for line in file:
            line = line.strip()
            if line.startswith('#'):
                continue
            if skip_lines:
                skip_lines = False
                continue
            line = re.sub(r' +', '\t', line)
            columns = line.split('\t')
            
            # Check if the line has at least 7 columns
            if len(columns) > 6:
                pfam = columns[6]
                pfams.add(pfam)

    return list(pfams)

def main():
    parser = argparse.ArgumentParser(description="Calculates k-mer profiles for every file.fasta file in the current directory (creates file.kmer files) and searches for the n (default=10) viruses with the least dissimilar (Bray-Curtis) k-mer profiles in the input database. The results are written to files named file_matches.txt with 'profile names' in the first column, and 'dissimilarity values' in the second column")
    parser.add_argument('-l', '--list', required=True, help='Input list of .fasta files to process')
    parser.add_argument('-d', '--database', required=True, help='Input k-mer profile database. Example: db/VMR_MSL38_v1_complete_genomes.kmers')
    parser.add_argument('-f', '--fasta_database', required=True, help='Input fasta database. Example: db/VMR_MSL38_v1_complete_genomes.fasta')
    parser.add_argument('-p', '--pfam_path', required=True, help='Path to directory where Pfam-A.hmm database is located. Example: db/')
    parser.add_argument('-o', '--output_file', required=True, help='Output file name. Example: results.txt')
    parser.add_argument('-k', '--kmer_length', type=int, default=4, help='OPTIONAL: Length of k-mers. Default = 4')
    parser.add_argument('-n', '--n_matches', type=int, default=10, help='OPTIONAL: Number of top matches to display in output file. Default = 10')
    args = parser.parse_args()

    # Read the input k-mer profile database
    #print(f"Reading input k-mer profiles database: {args.database}")
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

    with open(args.output_file, 'w') as output_file:
        #all_results = []
        # Calculate k-mer profiles for all .fasta files in the current directory
        with open(args.list, 'r') as input_list:

            for file_name in input_list:
                file_name = file_name.strip()
            #if file_name.endswith('.fasta'):
                all_results = []
                print(f"Creating {args.kmer_length}-mer profile for {file_name}")
                #output_file.write(f"Query: {file_name}\n")
                kmer_profile = defaultdict(int)
                kmer_profile_rc = defaultdict(int)
                with open(file_name, 'r') as fasta_file:
                    sequence_l=""
                    for line in fasta_file:
                        if line.startswith('>'):
                            continue
                        sequence = line.strip()
                        reversed_sequence = sequence[::-1]
                        complement_dict = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
                        complement_sequence = ''.join([complement_dict[base] for base in reversed_sequence])
                        sequence_l += line.strip()
                        sequence_length = len(sequence_l)
                        #print(f"{file_name} has {sequence_length} bp")
                        kmer_counts = get_kmer_counts(sequence, kmers)
                        for kmer, count in kmer_counts.items():
                            kmer_profile[kmer] += count

                        kmer_counts_rc = get_kmer_counts(complement_sequence, kmers)
                        for kmer, count in kmer_counts_rc.items():
                            kmer_profile_rc[kmer] += count
                    #print(f"{file_name} has {sequence_length} bp")
                    #output_file.write(f"{file_name} has {sequence_length} bp\n")
    
                # Save the new k-mer profile to a .kmer file with the same name as the input fasta file
                #kmer_file_name = os.path.splitext(file_name)[0] + '.kmer'
                #with open(kmer_file_name, 'w') as kmer_file:
                #    kmer_file.write(f"{file_name}\t")
                #    kmer_file.write('\t'.join(str(kmer_profile[kmer]) for kmer in kmers))
                #    kmer_file.write('\n')
                #print(f"{file_name} {args.kmer_length}-mer profile written to {kmer_file_name}")
    
                # Calculate Bray-Curtis dissimilarity to every k-mer profile in the database
                print(f"Comparing {file_name} k-mer profile with k-mer profiles in {args.database}")
                n = args.n_matches
                top_n_matches = [(None, float('inf'))] * n
                top_n_matches_rc = [(None, float('inf'))] * n
    
                for profile_name, profile_values in kmer_profiles.items():
                    dissimilarity = calculate_bray_curtis(kmer_profile, profile_values)
                    for i, (match_name, match_value) in enumerate(top_n_matches):
                        if dissimilarity < match_value:
                            top_n_matches.insert(i, (profile_name, dissimilarity))
                            break
                    top_n_matches = top_n_matches[:n]

                for profile_name, profile_values in kmer_profiles.items():
                    dissimilarity = calculate_bray_curtis(kmer_profile_rc, profile_values)
                    for i, (match_name, match_value) in enumerate(top_n_matches_rc):
                        if dissimilarity < match_value:
                            top_n_matches_rc.insert(i, (profile_name, dissimilarity))
                            break
                    top_n_matches_rc = top_n_matches_rc[:n]

                combined_top_matches = top_n_matches + top_n_matches_rc
                combined_top_matches.sort(key=lambda x: x[1])
                top_matches = combined_top_matches[:n]
                
                print(f"Running prodigal for {file_name}")
                run_prodigal(file_name)
                faa_file_name = os.path.splitext(file_name)[0] + '.faa'
                with open(faa_file_name, 'r') as faa_file:
                    n_faa = 0
                    for line in faa_file:
                        if line.startswith(">"):
                            n_faa += 1

                database_path = args.pfam_path

                print(f"Running pfamscan for {file_name}")
                run_pfamscan(faa_file_name, database_path)
                pfamscan_out = os.path.splitext(faa_file_name)[0] + '.pfamscan'
                unique_pfams = get_pfams(pfamscan_out)

                #result_file_name = os.path.splitext(file_name)[0] + '_matches.txt'
                #with open(args.output_file, 'w') as output_file:
                query_pfams = ", ".join(unique_pfams)
                
                output_file.write(f"# Query: {file_name}\n# Domains: {query_pfams}\n# Top {n} matches for {file_name} ({sequence_length} bp, {n_faa} CDS):\n")
                print(f"Top {n} matches for {file_name}:")
                output_file.write(f"# Profile Name\tKmer Dissimilarity\tKmer Similarity\tLength\tLength Proportion\tCDS count\tCDS proportion\tUnion Domains\tQuery Unique Domains\tSubject Unique Domains\tIntersection Domains\tDomain count dissimilarity\tDomain count similarity\tOverall score\n")
                #print("Profile Name\tDissimilarity\tSimilarity\tLength\tLength Proportion")
                match_count = 0
                for match in top_matches:
                    profile_name, dissimilarity = match
                    similarity = 1 - dissimilarity
                    fasta_length = 0
                    fasta_header = None
                    match_count += 1

                    print(f"{match_count}: {profile_name}")
                    accession = profile_name.split("|")[0]
                    output_fasta_name = f"FILES/{accession}.fasta"
                    #output_fasta_name = os.path.splitext(file_name)[0] + "_" + str(match_count) + ".fasta"
                    if os.path.exists(output_fasta_name):
                        print(f"{output_fasta_name} already exists. Moving forward...")
                        with open(output_fasta_name, 'r') as existing_file:
                            for line in existing_file:
                                if line.startswith('>'):
                                    continue
                                else:
                                    sequence = line.strip()
                                    fasta_length += len(sequence)
                    else:
                        print(f"Extracting {accession} sequence from {args.fasta_database} database")
                        os.makedirs("FILES", exist_ok=True)
                        with open(output_fasta_name, "w") as output_fasta:
                            with open(args.fasta_database, 'r') as fasta_db:
                                sequence_found = False
                                for line in fasta_db:
                                    if line.startswith('>'):
                                        if sequence_found:
                                            break
                                        fasta_header = line[1:].strip()                                
                                        if fasta_header == profile_name:
                                            sequence_found = True
                                    elif sequence_found:
                                        sequence = line.strip()
                                        fasta_length += len(sequence)
                                        output_fasta.write(">" + fasta_header + "\n")
                                        output_fasta.write(sequence + "\n")
                                        break

                    faa_file_name = f"{os.path.splitext(output_fasta_name)[0]}.faa"
                    if os.path.exists(faa_file_name):
                        print(f"{faa_file_name} already exists. Moving forward...")
                    else:
                        print(f"Running prodigal for {accession}")
                        run_prodigal(output_fasta_name)
                    
                    with open(faa_file_name, 'r') as faa_file:
                        faa_count = 0
                        for line in faa_file:
                            if line.startswith(">"):
                                faa_count += 1

                    if faa_count == n_faa:
                        faa_proportion = 1
                    else:
                        minimum_value = min(faa_count, n_faa)
                        maximum_value = max(faa_count, n_faa)
                        faa_proportion = minimum_value / maximum_value                         

                    if fasta_length == sequence_length:
                        length_proportion = 1
                    else:
                        minimum_value = min(fasta_length, sequence_length)
                        maximum_value = max(fasta_length, sequence_length)
                        length_proportion = minimum_value / maximum_value

                    pfamscan_out_file = f"{os.path.splitext(faa_file_name)[0]}.pfamscan"
                    if os.path.exists(pfamscan_out_file):
                        print(f"{pfamscan_out_file} already exists. Moving forward...")
                    else:
                        print(f"Running pfamscan for {accession}")
                        run_pfamscan(faa_file_name, database_path)
                    
                    unique_pfams_match = get_pfams(pfamscan_out_file)
                    combined_pfams = set(unique_pfams + unique_pfams_match)
                    unique_query_pfams = set(unique_pfams) - set(unique_pfams_match)
                    unique_subject_pfams = set(unique_pfams_match) - set(unique_pfams)
                    intersection_pfams = set(unique_pfams).intersection(unique_pfams_match)
                    pfam_counts_out = Counter(unique_pfams)
                    pfam_counts_out_file = Counter(unique_pfams_match)
                    all_pfams = list(combined_pfams)
                    union_domains = ", ".join(all_pfams)
                    query_only_domains = ", ".join(unique_query_pfams)
                    subject_only_domains = ", ".join(unique_subject_pfams)
                    intersection_domains = ", ".join(intersection_pfams)
                    pfam_count_array_out = [pfam_counts_out[pfam] for pfam in all_pfams]
                    pfam_count_array_out_file = [pfam_counts_out_file[pfam] for pfam in all_pfams]
                    bray_curtis_dissimilarity = distance.braycurtis(pfam_count_array_out, pfam_count_array_out_file)
                    bray_curtis_similarity = 1 - bray_curtis_dissimilarity
                    average_score = (similarity + length_proportion + faa_proportion + bray_curtis_similarity)/4

                    result = (profile_name, dissimilarity, similarity, fasta_length, length_proportion, faa_count, faa_proportion,  union_domains, query_only_domains, subject_only_domains, intersection_domains, bray_curtis_dissimilarity, bray_curtis_similarity, average_score)
                    all_results.append(result)
                sorted_results = sorted(all_results, key=lambda x: x[-1], reverse=True)
                for result in sorted_results:
                    output_file.write('\t'.join(map(str, result)) + '\n')
                    #print(f"{profile_name}\t{dissimilarity}\t{similarity}\t{fasta_length}\t{length_proportion}")

    print(f"Results written to {args.output_file}")

if __name__ == "__main__":
    main()
