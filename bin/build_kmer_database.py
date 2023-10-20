#!/usr/bin/env python
import argparse
from Bio import Entrez, SeqIO
import itertools
import os
import string

def download_genbank_record(accession):
    Entrez.email = args.email
    handle = Entrez.efetch(db="nucleotide", id=accession, rettype="gb", retmode="text")
    return SeqIO.read(handle, "genbank")

def combine_sequences(sequences):
    return "X".join(sequences)

def format_accession(accession):
    return accession.split(".")[0]

def format_header(accessions, annotation_file):
    full_accessions = "|".join(accessions)
    accession_str = full_accessions.split("|")[0]
    #print(f'Looking for {accession_str}') 
    translator = str.maketrans(",;:.-()[]/ ", "___________")
    matching_line = None
    with open(annotation_file, 'r', encoding='latin-1') as file1:
        for file1_line in file1:
            file1_line = file1_line.translate(translator)
            file1_columns = file1_line.strip().split('\t')
            file1_accessions = file1_columns[21]
            if accession_str in file1_accessions:
                #print(f'in {file1_accessions}')
                matching_line = file1_columns
                break
    if matching_line is not None:
        header = f'>{full_accessions}|{matching_line[18]}-{matching_line[16]}-{matching_line[14]}-{matching_line[12]}-{matching_line[10]}-{matching_line[8]}-{matching_line[6]}-{matching_line[4]}-{matching_line[2]}-Viruses-{matching_line[24]}'
    return header

def process_line(line):
    accessions = line.strip().split(";")
    records = []
    for accession in accessions:
        record = download_genbank_record(accession)
        records.append(record)
    if len(records) == 1:
        record = records[0]
        header = format_header([format_accession(record.id)], args.annotation_file)
        return header, str(record.seq)
    else:
        combined_sequences = [str(record.seq) for record in records]
        accessions = [format_accession(record.id) for record in records]
        header = format_header(accessions, args.annotation_file)
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
            print(f"DNA sequences downloaded and written to {args.output_prefix}.fasta file")
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
        print(f"kmer profiles written to {args.output_prefix}.kmers file")

    print(f"Creating DNA.kmers and RNA.kmers files")
    rna_filename = 'RNA.kmers'
    dna_filename = 'DNA.kmers'
    target_strings = ["RNA", "Riboviria", "dsDNA_RT"]
    with open(f"{args.output_prefix}.kmers", "r") as kmer_database, open(rna_filename, 'w') as rna_file, open(dna_filename, 'w') as dna_file:
        header = next(kmer_database)
        rna_file.write(header)
        dna_file.write(header)
        for line in kmer_database:
            if any(target in line for target in target_strings):
                rna_file.write(line)
            else:
                dna_file.write(line)    
    print("Lines containing 'RNA', 'Riboviria', or 'dsDNA_RT' saved in RNA.kmers")
    print("Other lines saved in DNA.kmers")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Download and process GenBank records from a list of accessions to obtain complete genomes in fasta format and create a kmer profile database. Will also generate a fasta file with complete genomes and separate .kmers files for DNA and RNA viruses.")
    parser.add_argument("-i", "--input_list", required=True, help="Input file containing a list of GenBank accessions. Each line may have one or more accessions (segmented genomes) separated by ;")
    parser.add_argument("-e", "--email", required=True, help="Your email address for NCBI Entrez utilities")
    parser.add_argument("-o", "--output_prefix", required=True, help="Prefix for the output files. Example: VMR_MSL38_v1_complete_genomes")
    parser.add_argument("-a", "--annotation_file", required=True, help="Annotation file needed to name fasta headers. Example: lists/VMR_MSL38_v1_complete_genomes.txt")
    parser.add_argument('-k', '--kmer_length', default=4, type=int, help='OPTIONAL: Length of k-mers. Default = 4')
    args = parser.parse_args()
    main(args)
