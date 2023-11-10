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

def get_pfams(file_path, accession_filter=None):
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
                if accession_filter is None or (accession_filter is not None and accession_filter in line):
                    pfams.add(pfam)

    return list(pfams)

def run_blastn(query, database, output):
    err_file = 'blastn.err'
    out_file = output
    evalue = '1e-20'
    with open(out_file, 'w') as out_blast_file, open(err_file, 'w') as err_blast_file:
        subprocess.run(['blastn', '-query', query, '-db', database, '-outfmt', '6 sacc pident', '-evalue', evalue], stdout=out_blast_file, stderr=err_blast_file)

def run_blastp(query, database, output):
    err_file = 'blastp.err'
    out_file = output
    evalue = '1e-10'
    with open(out_file, 'w') as out_blast_file, open(err_file, 'w') as err_blast_file:
        subprocess.run(['blastp', '-query', query, '-db', database, '-outfmt', '6 sacc pident', '-evalue', evalue], stdout=out_blast_file, stderr=err_blast_file)

def main():
    parser = argparse.ArgumentParser(description="Calculates k-mer profiles for every fasta file in input list, and compares them against the k-mer profiles in the input database. It also evaluates nucleotide, amino acid, domain content, number of protein coding sequences, and genome length similarities to compute a relatedness score")
    parser.add_argument('-l', '--list', required=True, help='Input list of .fasta files to process')
    parser.add_argument('-d', '--database', required=True, help='Path to directory where VirTaK database is located. Must include database name without extension. Example: db/VMR_MSL38_v1_complete_genomes')
    parser.add_argument('-p', '--pfam_path', required=True, help='Path to directory where Pfam-A.hmm database is located. Example: db/')
    parser.add_argument('-o', '--output_file', required=True, help='Output file name. Example: results.txt')
    parser.add_argument('-k', '--kmer_length', type=int, default=4, help='OPTIONAL: Length of k-mers. Default = 4')
    parser.add_argument('-n', '--n_matches', type=int, default=10, help='OPTIONAL: Number of top matches to display in output file. Default = 10')
    args = parser.parse_args()

    with open(args.list, 'r') as input_list:
        n_files = sum(1 for line in input_list)
        print(f"==================================================\nReading {args.list}\n{n_files} files in {args.list}")
    # Open output file
    with open(args.output_file, 'w') as output_file:
        # Open input list
        with open(args.list, 'r') as input_list:
            # Open input file
            query_number = 1
            for file_name in input_list:
                file_name = file_name.strip()
                print(f"==================================================\nQuery {query_number}: {file_name}")
                query_number += 1
                all_results = []

                # Run blastn to detect closely related genomes
                print("Searching for closely related genomes")
                nucl_database = f"{args.database}.fasta"
                blastn_output = os.path.splitext(file_name)[0] + '_blastn.out'
                if os.path.exists(blastn_output):
                    print(f"{blastn_output} already exists. Moving forward...")
                else:
                    print(f"Running blastn for {file_name} against {nucl_database}")
                    run_blastn(file_name, nucl_database, blastn_output)
                #blastn_matches = set()
                blastn_dict = {}
                with open(blastn_output, 'r') as file:
                    for line in file:
                        fields = line.strip().split('\t')
                        if len(fields) == 2:
                            accessions = fields[0].split('|')[0]
                            pident = float(fields[1]) / 100
                            #blastn_matches.add(accessions)
                            if accessions in blastn_dict:
                                old_pident, count = blastn_dict[accessions]
                                new_pident = (old_pident * count + pident) / (count + 1)
                                blastn_dict[accessions] = (new_pident, count + 1)
                            else:
                                blastn_dict[accessions] = (pident, 1)
                #unique_blastn_matches = len(blastn_matches)
                unique_blastn_matches = len(blastn_dict)
                print(f"Found {unique_blastn_matches} closely related genomes")
                #print(blastn_matches)

                # Run prodigal to predict protein coding genes
                print(f"Predicting protein coding genes")
                faa_file_name = os.path.splitext(file_name)[0] + '.faa'
                if os.path.exists(faa_file_name):
                    print(f"{faa_file_name} already exists. Moving forward...")
                else:
                    print(f"Running prodigal for {file_name}")
                    run_prodigal(file_name)
                    print(f"{faa_file_name} created")
                with open(faa_file_name, 'r') as faa_file:
                    n_faa = 0
                    for line in faa_file:
                        if line.startswith(">"):
                            n_faa += 1
                print(f"There are {n_faa} protein coding genes in {file_name}")

                # Run blastp to detect genomes with homologous proteins
                print("Searching for genomes with homologous proteins")
                prot_database = f"{args.database}.faa"
                blastp_output = os.path.splitext(file_name)[0] + '_blastp.out'
                if os.path.exists(blastp_output):
                    print(f"{blastp_output} already exists. Moving forward...")
                else:                
                    print(f"Running blastp for {faa_file_name} against {prot_database}")
                    run_blastp(faa_file_name, prot_database, blastp_output)
                #blastp_matches = set()
                blastp_dict = {}
                with open(blastp_output, 'r') as file:
                    for line in file:
                        fields = line.strip().split('\t')
                        if len(fields) == 2:
                            accessions = fields[0].split('|')[0]
                            pident = float(fields[1]) / 100
                            #blastp_matches.add(accessions)
                            if accessions in blastp_dict:
                                old_pident, count = blastp_dict[accessions]
                                new_pident = (old_pident * count + pident) / (count + 1)
                                blastp_dict[accessions] = (new_pident, count + 1)
                            else:
                                blastp_dict[accessions] = (pident, 1)
                #unique_blastp_matches = len(blastp_matches)
                max_count = max(count for pident, count in blastp_dict.values())                
                normalized_blastp_dict = {}
                for accession, (pident, count) in blastp_dict.items():
                    normalized_pident = (pident * count) / max_count
                    normalized_blastp_dict[accession] = normalized_pident
                unique_blastp_matches = len(normalized_blastp_dict)
                print(f"Found {unique_blastp_matches} genomes with homologous proteins")
                #print(blastp_matches)

                # Run pfamscan to detect distantly related genomes
                print("Searching for distantly related genomes")
                database_path = args.pfam_path
                pfamscan_out = os.path.splitext(faa_file_name)[0] + '.pfamscan'
                if os.path.exists(pfamscan_out):
                    print(f"{pfamscan_out} already exists. Moving forward...")
                else:
                    print(f"Running pfamscan for {file_name}")
                    run_pfamscan(faa_file_name, database_path)
                unique_pfams = get_pfams(pfamscan_out)
                #print(unique_pfams)
                query_pfams = ", ".join(unique_pfams)
                # Get list of genomes that share at least one domain
                pfamscan_database = f"{args.database}.pfamscan"
                matching_accessions = set()
                with open(pfamscan_database, 'r') as file:
                    for line in file:
                        line = line.strip()
                        if not line.startswith('#'):
                            line = re.sub(r' +', '\t', line)
                            columns = line.split('\t')
                            if len(columns) > 6:
                                pfam = columns[6]
                                accession = columns[0].split('|')[0]
                                if pfam in unique_pfams:
                                    matching_accessions.add(accession)
                non_redundant_accessions = list(matching_accessions)
                n_matching_accessions = len(non_redundant_accessions)
                print(f"Found {n_matching_accessions} genomes that share at least one domain")
                #print(non_redundant_accessions)
                # Get all domains of genomes that share at leas one domain
                all_domains = {}
                for accession in non_redundant_accessions:
                    all_domains[accession] = set()
                with open(pfamscan_database, 'r') as file:
                    for line in file:
                        line = line.strip()
                        if not line.startswith('#'):
                            line = re.sub(r' +', '\t', line)
                            columns = line.split('\t')
                            if len(columns) > 6:
                                pfam = columns[6]
                                accession = columns[0].split('|')[0]
                                if accession in non_redundant_accessions:
                                    all_domains[accession].add(pfam)
                # Calculate domain count similarity
                pfam_count_similarities = {}
                for accession in all_domains:
                    query_counts = Counter(unique_pfams)
                    subject_counts = Counter(all_domains[accession])
                    all_domains[accession] |= set(unique_pfams)
                    all_pfams = list(all_domains[accession])
                    query_array = [query_counts[pfam] for pfam in all_pfams]
                    subject_array = [subject_counts[pfam] for pfam in all_pfams]
                    pfam_count_dissimilarity = distance.braycurtis(query_array, subject_array)
                    pfam_count_similarity = 1 - pfam_count_dissimilarity
                    pfam_count_similarities[accession] = pfam_count_similarity

                # Run kmer_search to detect genomes with similar kmer profile
                print("Searching for genomes with similar kmer profiles")
                # Read the input k-mer profile database
                kmer_profiles = {}
                kmer_database = f"{args.database}.kmers"
                with open(kmer_database, 'r') as db_file:
                    header = db_file.readline().strip().split('\t')
                    kmers = header[1:]
                    for line in db_file:
                        parts = line.strip().split('\t')
                        profile_name = parts[0]
                        profile_values = list(map(int, parts[1:]))
                        kmer_profiles[profile_name] = dict(zip(kmers, profile_values))
                # Generate all possible k-mers
                kmers = generate_kmers(args.kmer_length)
                # Generate k-mer profile for input file
                print(f"Generating {args.kmer_length}-mer profile for {file_name}")
                kmer_profile = defaultdict(int)
                kmer_profile_rc = defaultdict(int)
                with open(file_name, 'r') as fasta_file:
                    sequence_l=""
                    for line in fasta_file:
                        if line.startswith('>'):
                            continue
                        sequence = line.strip()
                        sequence = sequence.upper()
                        reversed_sequence = sequence[::-1]
                        complement_dict = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'X': 'X', 'V': 'V', 'N': 'N', 'R': 'Y', 'Y': 'R', 'W': 'W', 'S': 'S'}
                        complement_sequence = ''.join([complement_dict[base] for base in reversed_sequence])
                        sequence_l += line.strip()
                        sequence_length = len(sequence_l)
                        kmer_counts = get_kmer_counts(sequence, kmers)
                        kmer_counts_rc = get_kmer_counts(complement_sequence, kmers)
                        combined_kmer_counts = {kmer: kmer_counts[kmer] + kmer_counts_rc[kmer] for kmer in kmers}
                        for kmer, count in combined_kmer_counts.items():
                            kmer_profile[kmer] += count
                print(f"Comparing {file_name} k-mer profile with k-mer profiles in {kmer_database}")
                all_matches = []
                for profile_name, profile_values in kmer_profiles.items():
                    dissimilarity = calculate_bray_curtis(kmer_profile, profile_values)
                    all_matches.append((profile_name, dissimilarity))

                # Get faa_proportion
                print("Calculating CDS proportions")
                accession_counts = defaultdict(int)
                with open(prot_database, 'r') as faa_file:
                    for line in faa_file:
                        if line.startswith('>'):
                            header = line.strip()
                            accession = header.split('|')[0][1:]
                            accession_counts[accession] += 1
                faa_proportion_dict = {}
                for accession in accession_counts:
                    faa_count = accession_counts[accession]
                    if faa_count == n_faa:
                        faa_proportion = 1
                    else:
                        minimum_value = min(faa_count, n_faa)
                        maximum_value = max(faa_count, n_faa)
                        faa_proportion = minimum_value / maximum_value
                    faa_proportion_dict[accession] = faa_proportion

                # Get length_proportion
                print("Calculating length proportions")
                sequence_lengths = {}
                with open(nucl_database, 'r') as fasta_db:
                    header = None
                    sequence = []
                    for line in fasta_db:
                        line = line.strip()
                        if line.startswith('>'):
                            if header is not None and sequence:
                                fasta_length = len(''.join(sequence))
                                accessions = header.split('|')[0][1:]
                                sequence_lengths[accessions] = fasta_length
                            header = line
                            sequence = []
                        else:
                            sequence.append(line)
                    if header is not None and sequence:
                        fasta_length = len(''.join(sequence))
                        accessions = header.split('|')[0][1:]
                        sequence_lengths[accessions] = fasta_length
                length_proportions = {}
                for accession in sequence_lengths:
                    fasta_length = sequence_lengths[accession]
                    if fasta_length == sequence_length:
                        length_proportion = 1
                    else:
                        minimum_value = min(fasta_length, sequence_length)
                        maximum_value = max(fasta_length, sequence_length)
                        length_proportion = minimum_value / maximum_value
                    length_proportions[accession] = length_proportion

                # Get statistics
                for profile_name, dissimilarity in all_matches:
                    similarity = 1 - dissimilarity
                    genome = profile_name.split('|')[0]
                    blastn_proportion, _ = blastn_dict.get(genome, (0, 0))
                    blastp_proportion = normalized_blastp_dict.get(genome, 0)
                    domain_count_similarity = pfam_count_similarities.get(genome, 0)
                    faa_proportion_value = faa_proportion_dict.get(genome, 0)
                    length_proportion_value = length_proportions.get(genome, 0)
                    w1 = 0.285714286 * blastn_proportion
                    w2 = 0.238095238 * blastp_proportion
                    w3 = 0.19047619 * domain_count_similarity
                    w4 = 0.142857143 * similarity
                    w5 = 0.095238095 * faa_proportion_value
                    w6 = 0.047619048 * length_proportion_value
                    score = w1 + w2 + w3 + w4 + w5 + w6
                    if score >= 0.714285714:
                        prediction = "Closely related"
                    elif score >= 0.476190476:
                        prediction = "Related"
                    elif score >= 0.285714286:
                        prediction = "Distantly related"
                    elif score >= 0.142857143:
                        prediction = "Similar genomes"
                    else:
                        prediction = "Not related"

                    #full report vs short report... length, blastn matches, cds count, blastp matches, union domains, query unique domains, subject unique domains, intersection, raw score, weighted score
                    result = (profile_name, blastn_proportion, blastp_proportion, domain_count_similarity, similarity, faa_proportion_value, length_proportion_value, score, prediction)
                    all_results.append(result)

                print(f"Top {args.n_matches} matches for {file_name}:")
                output_file.write(f"# Query: {file_name}\n# Domains: {query_pfams}\n# Top {args.n_matches} matches for {file_name} ({sequence_length} bp, {n_faa} CDS):\n")
                output_file.write("Virus\tBlastn ANI\tBlastp AAI\tDomain Count Similarity\tKmer Similarity\tCDS proportion\tLength Proportion\tScore\tPrediction\n")
                sorted_results = sorted(all_results, key=lambda x: x[-2], reverse=True)
                top_n_results = sorted_results[:args.n_matches]
                result_count = 1
                for result in top_n_results:
                    result_str = '\t'.join(map(str, result))
                    output_file.write(result_str + '\n')
                    print(f"{result_count}: {result_str}")
                    result_count += 1

    print(f"==================================================\nResults written to {args.output_file}")

if __name__ == "__main__":
    main()
