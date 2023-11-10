#!/usr/bin/env python
import argparse
from Bio import Entrez, SeqIO
import itertools
import os
import string
from Bio.Seq import Seq
import subprocess
from urllib.error import HTTPError

def download_genbank_record(accession):
    try:
        Entrez.email = args.email
        handle = Entrez.efetch(db="nucleotide", id=accession, rettype="gb", retmode="text")
        return SeqIO.read(handle, "genbank")
    except HTTPError as e:
        print(f"HTTPError: {e}")

def combine_sequences(sequences):
    return "X".join(sequences)

def format_accession(accession):
    return accession.split(".")[0]

def format_header(accessions, annotation_file):
    full_accessions = "|".join(accessions)
    accession_str = full_accessions.split("|")[0]
    translator = str.maketrans(",;:.-()[]/ ", "___________")
    matching_line = None
    with open(annotation_file, 'r', encoding='latin-1') as file1:
        for file1_line in file1:
            file1_line = file1_line.translate(translator)
            file1_columns = file1_line.strip().split('\t')
            file1_accessions = file1_columns[21]
            if accession_str in file1_accessions:
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

def run_prodigal(file_name):
    prodigal_out = "prodigal.out"
    with open(prodigal_out, 'w') as prodigal_out_file:
        fasta_file = file_name
        faa_file = os.path.splitext(file_name)[0] + '.faa'
        subprocess.run(['prodigal', '-i', fasta_file, '-a', faa_file, '-p', 'meta', '-q'], stdout=prodigal_out_file)

def create_blast_db(input_file, db_type):
    out_file = f"{db_type}_db.out"
    err_file = f"{db_type}_db.err"
    with open(out_file, 'w') as out_file_handle, open(err_file, 'w') as err_file_handle:
        subprocess.run(['makeblastdb', '-in', input_file, '-dbtype', db_type], stdout=out_file_handle, stderr=err_file_handle)

def run_pfamscan(in_file, database):
    out_file = os.path.splitext(in_file)[0] + '.pfamscan'
    subprocess.run(['pfam_scan.pl', '-fasta', in_file, '-dir', database, '-outfile', out_file])        

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
                first_field_accession = accessions[0].split("|")[0]
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

    dbtype = "nucl"
    if os.path.exists(f"{output_fasta_path}.nhr") or os.path.exists(f"{output_fasta_path}.nin") or os.path.exists(f"{output_fasta_path}.nsq"):
        print("Nucleic acid database already created")
    else:
        create_blast_db(output_fasta_path, dbtype)
        print("Nucleic acid database successfully created")

    kmers = generate_kmers(args.kmer_length)

    output_kmers_path = f"{args.output_prefix}.kmers"
    if not os.path.isfile(output_kmers_path):
        with open(f"{args.output_prefix}.kmers", "w") as kmer_database:
            kmer_database.write("Fasta Record Name\t" + "\t".join(kmers) + "\n")
            for record in SeqIO.parse(output_fasta_path, "fasta"):
                identifier = record.description
                print(f"Generating kmer profile for: {identifier}")
                forward_seq = str(record.seq)
                reverse_seq = str(Seq(forward_seq).reverse_complement())
                forward_kmer_counts = get_kmer_counts(forward_seq, kmers)
                reverse_kmer_counts = get_kmer_counts(reverse_seq, kmers)
                combined_kmer_counts = {kmer: forward_kmer_counts[kmer] + reverse_kmer_counts[kmer] for kmer in kmers}
                kmer_database.write(identifier + "\t" + "\t".join(str(combined_kmer_counts[kmer]) for kmer in kmers) + "\n")
            print(f"kmer profiles written to {args.output_prefix}.kmers file")
    else:
        print(f"{args.output_prefix}.kmers file already exists. Moving forward...")

    faa_file = f"{args.output_prefix}.faa"
    if not (os.path.exists(faa_file)):
        print(f"Running prodigal for {output_fasta_path}")
        run_prodigal(output_fasta_path)
        print(f"Protein translations written to {faa_file}")
    else:
        print(f"{faa_file} already exists. Moving forward...")

    dbtype = "prot"
    if os.path.exists(f"{faa_file}.phr") or os.path.exists(f"{faa_file}.pin") or os.path.exists(f"{faa_file}.psq"):
        print("Protein database already created")
    else:
        create_blast_db(faa_file, dbtype)
        print("Protein database successfully created")

    pfam_database_path = args.pfam_path
    pfamscan_file = f"{args.output_prefix}.pfamscan"
    if not (os.path.exists(pfamscan_file)):
        print(f"Running pfamscan for {faa_file}")
        run_pfamscan(faa_file, pfam_database_path)
        print(f"pfamscan results written to {pfamscan_file}")
    else:
        print(f"{pfamscan_file} already exists. Moving forward...")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Create VirTaK database. From input file containing a list of GenBank accessions downloads and process fasta genomes to: create blastn database, predict protein coding sequences, create blastp database, create kmer profile database, and create a protein domain database")
    parser.add_argument("-i", "--input_list", required=True, help="Input file containing a list of GenBank accessions. Each line may have one or more accessions (segmented genomes) separated by ;")    
    parser.add_argument("-o", "--output_prefix", required=True, help="Prefix for the output files. Example: VMR_MSL38_v1_complete_genomes")
    parser.add_argument("-a", "--annotation_file", required=True, help="Annotation file needed to name fasta headers. Example: lists/VMR_MSL38_v1_complete_genomes.txt")
    parser.add_argument("-p", "--pfam_path", required=True, help="Path to directory where Pfam-A.hmm database is located. Example: db/")
    parser.add_argument("-e", "--email", required=True, help="Your email address for NCBI Entrez utilities")
    parser.add_argument('-k', '--kmer_length', default=4, type=int, help='OPTIONAL: Length of k-mers. Default = 4')
    args = parser.parse_args()
    main(args)
