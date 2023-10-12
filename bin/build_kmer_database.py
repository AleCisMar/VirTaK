#!/usr/bin/env python
import argparse
from Bio import Entrez, SeqIO
import itertools
import os

def download_genbank_record(accession):
    Entrez.email = args.email
    handle = Entrez.efetch(db="nucleotide", id=accession, rettype="gb", retmode="text")
    return SeqIO.read(handle, "genbank")

def combine_sequences(sequences):
    return "X".join(sequences)

def format_accession(accession):
    return accession.split(".")[0]

def format_header(accessions, organism, full_taxonomy):
    accession_str = "|".join(accessions)
    taxonomy_str = "_".join(reversed(full_taxonomy))
    header = ">{0}|{1}_{2}".format(accession_str, organism, taxonomy_str)
    header = header.replace("(", "_").replace(")", "_").replace(",", "_").replace(":", "_").replace(";", "_").replace(" ", "_").replace("/", "_")
    return header

def process_line(line):
    accessions = line.strip().split(";")
    records = []
    for accession in accessions:
        record = download_genbank_record(accession)
        records.append(record)
    if len(records) == 1:
        record = records[0]
        organism = record.annotations.get("organism", "")
        full_taxonomy = record.annotations.get("taxonomy", [])
        header = format_header([format_accession(record.id)], organism, full_taxonomy)
        return header, str(record.seq)
    else:
        combined_sequences = [str(record.seq) for record in records]
        accessions = [format_accession(record.id) for record in records]
        organism = records[0].annotations.get("organism", "")
        full_taxonomy = records[0].annotations.get("taxonomy", [])
        header = format_header(accessions, organism, full_taxonomy)
        combined_sequence = combine_sequences(combined_sequences)
        return header, combined_sequence

def generate_kmers(k):
    letters = ['A', 'T', 'G', 'C']
    all_kmers = [''.join(p) for p in itertools.product(letters, repeat=k)]
    return all_kmers

def get_kmer_counts(sequence, kmers):
    sequence_length = len(sequence)
    kmer_counts = {kmer: 0 for kmer in kmers}

    for i in range(sequence_length - len(kmers[0]) + 1):
        kmer = sequence[i:i+len(kmers[0])]
        if kmer in kmers:
            kmer_counts[kmer] += 1

    return kmer_counts

def append_to_existing_fasta(fasta_file_path, new_records):
    with open(fasta_file_path, "a") as fasta_file:
        for header, sequence in new_records:
            fasta_file.write(f"{header}\n{sequence}\n")

def main(args):
    output_fasta_path = f"{args.output_prefix}.fasta"

    if not os.path.isfile(output_fasta_path):
        # Output file doesn't exist, perform the normal download and processing
        with open(args.input_list, "r") as input_file, open(output_fasta_path, "w") as fasta_file:
            for line in input_file:
                accession = line.strip()
                print(f"Downloading DNA sequence for: {accession}")
                header, sequence = process_line(line)
                fasta_file.write(f"{header}\n{sequence}\n")
    else:
        # Output file exists, only download and append missing records
        print(f"{args.output_prefix}.fasta exists, will only download and append missing records")
        existing_accessions = set()
        with open(output_fasta_path, "r") as fasta_file:
            for line in fasta_file:
                if line.startswith(">"):
                    accession = line.split("|")[0][1:]
                    existing_accessions.add(accession)

        missing_records = []
        with open(args.input_list, "r") as input_file:
            for line in input_file:
                accessions = line.strip().split(";")
                first_field_accession = accessions[0].split("|")[0]  # Take the first field before "|"
                if first_field_accession in existing_accessions:
                    print(f"Skipping {accessions} as it is already downloaded.")
                else:
                    print(f"Downloading missing DNA sequence: {accessions}")
                    header, sequence = process_line(line)
                    if header and sequence:
                        missing_records.append((header, sequence))

        if missing_records:
            append_to_existing_fasta(output_fasta_path, missing_records)
            print("Missing records downloaded and appended to the existing .fasta file.")

    kmers = generate_kmers(args.kmer_length)

    with open(f"{args.output_prefix}.kmers", "w") as kmer_database:
        kmer_database.write("Fasta Record Name\t" + "\t".join(kmers) + "\n")
        for record in SeqIO.parse(output_fasta_path, "fasta"):
            identifier = record.description
            print(f"Generating kmer profile for: {identifier}")
            kmer_counts = get_kmer_counts(str(record.seq), kmers)
            kmer_database.write(identifier + "\t" + "\t".join(str(kmer_counts[kmer]) for kmer in kmers) + "\n")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Download and process GenBank records from a list of accessions to obtain complete genomes in fasta format and create a kmer profile database.")
    parser.add_argument("-i", "--input_list", required=True, help="Input file containing a list of GenBank accessions. Each line may have one or more accessions (segmented genomes) separated by ;")
    parser.add_argument("-e", "--email", required=True, help="Your email address for NCBI Entrez utilities")
    parser.add_argument("-o", "--output_prefix", required=True, help="Prefix for the output files")
    parser.add_argument('-k', '--kmer_length', default=4, type=int, help='OPTIONAL: Length of k-mers. Default = 4')
    args = parser.parse_args()
    main(args)
