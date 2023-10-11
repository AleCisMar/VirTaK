import argparse
import itertools
from Bio import SeqIO

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

def main():
    parser = argparse.ArgumentParser(description="Calculate k-mer counts for each DNA sequence in a multi-FASTA file")
    parser.add_argument("-i", "--input", required=True, help="Input multi-FASTA file")
    parser.add_argument("-k", "--kmer-length", type=int, required=True, help="Length of k-mers")
    parser.add_argument("-o", "--output", required=True, help="Output table file")

    args = parser.parse_args()

    # Generate all possible k-mers
    kmers = generate_kmers(args.kmer_length)

    with open(args.output, "w") as output_file:
        # Write the header row with k-mers
        output_file.write("Fasta Record Name\t" + "\t".join(kmers) + "\n")

        # Process the input multi-FASTA file
        for record in SeqIO.parse(args.input, "fasta"):
            # Get the first field of the identifier (before the first space)
            identifier = record.description

            # Calculate k-mer counts for the sequence
            kmer_counts = get_kmer_counts(str(record.seq), kmers)

            # Write the results to the output file
            output_file.write(identifier + "\t" + "\t".join(str(kmer_counts[kmer]) for kmer in kmers) + "\n")

    print(f"Number of k-mers: {len(kmers)}")

if __name__ == "__main__":
    main()
