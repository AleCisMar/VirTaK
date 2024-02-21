#!/usr/bin/env python
import argparse

def extract_pfam_info(input_file, pfam_seed_file, output_file):
    # Read the Pfam-A.seed file in binary mode and decode using 'latin-1'
    with open(pfam_seed_file, 'rb') as pfam_file:
        pfam_content = pfam_file.read().decode('latin-1')

    # Read the list of Pfam IDs from the input file
    with open(input_file, 'r') as input_file:
        pfam_ids = input_file.read().splitlines()

    # Split the Pfam-A.seed content into entries
    pfam_entries = pfam_content.split("# STOCKHOLM 1.0")[1:]

    # Create a dictionary to store Pfam information
    pfam_info_dict = {}

    # Populate the dictionary with Pfam information
    for entry in pfam_entries:
        lines = entry.split('\n')
        pfam_id = None
        accession = None
        description = None
        comments = []

        for line in lines:
            if line.startswith("#=GF ID"):
                pfam_id = line.split()[2]
            elif line.startswith("#=GF AC"):
                accession = line.split()[2].split('.')[0]
            elif line.startswith("#=GF DE"):
                description = ' '.join(line.split()[2:])
            elif line.startswith("#=GF CC"):
                comments.append(' '.join(line.split()[2:]))

        if pfam_id is not None:
            pfam_info_dict[pfam_id] = [accession, description, comments]

    # Write the formatted output to a tab-delimited file
    with open(output_file, 'w') as output_file:
        output_file.write("Pfam_ID\tAccession\tDescription\tComments\n")
        for pfam_id in pfam_ids:
            if pfam_id in pfam_info_dict:
                accession, description, comments = pfam_info_dict[pfam_id]
                formatted_comments = ' '.join(comments)
                output_file.write(f"{pfam_id}\t{accession}\t{description}\t{formatted_comments}\n")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Extract Pfam information based on Pfam IDs.")
    parser.add_argument("input_file", help="Path to the file containing Pfam IDs")
    parser.add_argument("pfam_seed_file", help="Path to Pfam-A.seed file")

    output_file = "domain_info.tsv"

    args = parser.parse_args()

    extract_pfam_info(args.input_file, args.pfam_seed_file, output_file)
