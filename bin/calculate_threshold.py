#!/usr/bin/env python
import argparse
import re
from collections import defaultdict
from itertools import combinations
from scipy.spatial import distance

def read_kmer_database(database_file):
    kmer_profiles = {}
    with open(database_file, 'r') as db_file:
        # Skip header
        next(db_file)
        for line in db_file:
            parts = line.strip().split('\t')
            profile_name = parts[0]
            profile_values = list(map(int, parts[1:]))
            kmer_profiles[profile_name] = profile_values
    return kmer_profiles

def read_input_file(input_file):
    entries = []
    with open(input_file, 'r') as in_file:
        for line in in_file:
            parts = line.strip().split('\t')
            fasta_name = parts[0]
            virus_name = parts[1]
            dissimilarity = float(parts[2])
            entries.append((fasta_name, virus_name, dissimilarity))
    return entries

def find_matching_suffix(entries, suffix):
    matching_words = set()
    for _, virus_name, _ in entries:
        # Step 1: Remove everything after the suffix
        virus_name = re.sub(suffix + r'.*$', suffix, virus_name)

        # Step 2: Extract the last field after the last "_"
        last_field = virus_name.rsplit('_', 1)[-1]
        matching_words.add(last_field)
    return matching_words

def remove_duplicates(matching_words):
    unique_words = set()
    non_redundant_list = []
    for word in matching_words:
        if word not in unique_words:
            non_redundant_list.append(word)
            unique_words.add(word)
    return non_redundant_list

def calculate_threshold(kmer_profiles, element):
    profiles_to_compare = []
    for profile_name in kmer_profiles.keys():
        if element in profile_name:
            profiles_to_compare.append(profile_name)

    max_dissimilarity = 0.0
    for pair in combinations(profiles_to_compare, 2):
        profile1 = kmer_profiles[pair[0]]
        profile2 = kmer_profiles[pair[1]]
        dissimilarity = 1.0 - (2.0 * sum(min(a, b) for a, b in zip(profile1, profile2)) /
                               (sum(profile1) + sum(profile2)))
        if dissimilarity > max_dissimilarity:
            max_dissimilarity = dissimilarity

    return max_dissimilarity

def count_matching_lines(database_file, element):
    count = 0
    with open(database_file, 'r') as db_file:
        # Skip header
        next(db_file)
        for line in db_file:
            if element in line:
                count += 1
    return count

def main():
    parser = argparse.ArgumentParser(description="Given the taxonomy of the virus with the best-matching profile, calculate the Bray-Curtis threshold for the specified taxonomic rank")
    parser.add_argument('-i', '--input', required=True, help='Input file. Is the output of kmer_search.py')
    parser.add_argument('-o', '--output', required=True, help='Output file')
    parser.add_argument('-d', '--database', required=True, help='Kmer profile database')
    parser.add_argument('-t', '--taxonomy', required=True, choices=["family", "order", "class", "phylum", "kingdom", "realm"], help='Suffix to search for (e.g., family, order, class, phylum, kingdom, realm)')
    args = parser.parse_args()

    suffixes = {
        "family": "viridae",
        "order": "virales",
        "class": "viricetes",
        "phylum": "viricota",
        "kingdom": "virae",
        "realm": "viria"
    }

    # Step 1: Read the input kmer profile database
    kmer_profiles = read_kmer_database(args.database)

    # Step 2: Read the input file
    entries = read_input_file(args.input)
    print(f"Read {len(entries)} entries from {args.input}")

    # Step 3: Find matching words for the given suffix
    print(f"Finding matching words with suffix '{suffixes[args.taxonomy]}'...")
    matching_words = find_matching_suffix(entries, suffixes[args.taxonomy])
    print(f"Matching words: {', '.join(matching_words)}")

    # Step 4: Remove duplicates from the matching words
    non_redundant_list = remove_duplicates(matching_words)
    print(f"Non-redundant list: {', '.join(non_redundant_list)}")

    # Step 5: Calculate the threshold for each element in the non-redundant list
    thresholds = {}
    print("Calculating thresholds...")
    for element in non_redundant_list:
        threshold = calculate_threshold(kmer_profiles, element)
        thresholds[element] = threshold
        print(f"Threshold for {element}: {threshold}")

        # Step 6: Print the number of matching lines in the input database
        matching_line_count = count_matching_lines(args.database, element)
        print(f"Number of matching lines in the database for {element}: {matching_line_count}")

    # Step 7: For each entry in the input file, print the threshold
    #print("\nResults:")
    with open(args.output, 'w') as output_file:
        for entry in entries:
            fasta_name, virus_name, dissimilarity = entry
            matching_word = re.sub(suffixes[args.taxonomy] + r'.*$', suffixes[args.taxonomy], virus_name)
            last_field = matching_word.rsplit('_', 1)[-1]
            threshold = thresholds[last_field]
            output_file.write(f"{fasta_name}\t{virus_name}\t{dissimilarity}\t{threshold}\n")

if __name__ == "__main__":
    main()