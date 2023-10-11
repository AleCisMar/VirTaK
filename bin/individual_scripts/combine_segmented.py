import argparse
from Bio import SeqIO

def combine_sequences(input_file, output_file):
    sequences_dict = {}
    
    # Read the input FASTA file and store sequences in a dictionary
    with open(input_file, "r") as fasta_file:
        for record in SeqIO.parse(fasta_file, "fasta"):
            header = record.description
            sequence = str(record.seq)
            
            # Extract the matching part after "|"
            matching_part = header.split("|")[1]
            
            # Extract the first field before "|" for this sequence
            first_field = header.split("|")[0]
            
            # Add the sequence and the first field to the dictionary under the matching part key
            if matching_part not in sequences_dict:
                sequences_dict[matching_part] = {"sequences": [], "first_fields": []}
            sequences_dict[matching_part]["sequences"].append(sequence)
            sequences_dict[matching_part]["first_fields"].append(first_field)

    # Write combined sequences to the output file
    with open(output_file, "w") as output_fasta:
        for matching_part, data in sequences_dict.items():
            sequences = data["sequences"]
            first_fields = data["first_fields"]
            
            # Combine sequences with "X" separator
            combined_sequence = "X".join(sequences)
            
            # Combine first fields into a single string
            combined_first_fields = "|".join(first_fields)
            
            # Create the combined header
            combined_header = f"{combined_first_fields}|{matching_part}"
            
            # Write combined sequence and original descriptions to the output file
            output_fasta.write(f">{combined_header}\n")
            output_fasta.write(f"{combined_sequence}\n")

def main():
    parser = argparse.ArgumentParser(description="Combine sequences of segmented genomes. Combined sequences are separated by X")
    parser.add_argument("-i","--input", required=True, help="input FASTA file")
    parser.add_argument("-o","--output", required=True, help="output FASTA file")
    args = parser.parse_args()

    input_file = args.input
    output_file = args.output
    combine_sequences(input_file, output_file)

if __name__ == "__main__":
    main()
