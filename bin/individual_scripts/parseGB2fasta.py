import argparse
from Bio import SeqIO

def parse_genbank(input_file, output_file):
    with open(input_file, "r") as genbank_file, open(output_file, "w") as fasta_file:
        for record in SeqIO.parse(genbank_file, "genbank"):
            accession = record.annotations["accessions"][0]
            organism = record.annotations.get("organism", "Unknown_Organism")
            taxonomy = "_".join(reversed(record.annotations.get("taxonomy", ["Unknown_Taxonomy"])))
            sequence = str(record.seq)

            # Combine accession, source, and taxonomy into a custom header
            header = f"{accession}|{organism}_{taxonomy}"

            # Write the header and sequence to the FASTA file
            fasta_file.write(f">{header}\n")
            fasta_file.write(f"{sequence}\n")

def main():
    parser = argparse.ArgumentParser(description="Parse GenBank file and create a FASTA file with custom headers.")
    parser.add_argument("-i","--input", required=True, help="Input GenBank file")
    parser.add_argument("-o","--output", required=True, help="Output FASTA file")
    args = parser.parse_args()

    input_file = args.input
    output_file = args.output

    parse_genbank(input_file, output_file)

if __name__ == "__main__":
    main()
