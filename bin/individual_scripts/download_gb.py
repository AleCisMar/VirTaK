from Bio import Entrez
import argparse

def download_genbank_records(id_list, out_file, e_mail, batch_size):
    OpenIDlist = open(id_list)
    ReadIDlist = OpenIDlist.read()
    SeqIDLists = ", ".join(ReadIDlist.split("\n"))

    GenomeSeqFile = out_file

    Entrez.email = e_mail
    print("The introduced E-mail is: " + Entrez.email)

    batchsize = batch_size
    n_ids = sum(1 for line in open(id_list))
    print("Fetching genbank records...")
    for start in range(0, n_ids, batchsize):
        handle = Entrez.efetch(db="nucleotide", id=SeqIDLists, rettype="gb", retmode="text", retstart=start, retmax=batchsize)
        print('Starting at {} Fetching {} records'.format(start, batchsize))
        print("Writing genbank records to: " + GenomeSeqFile)
        with open(GenomeSeqFile, "a") as GenomeSeqFile_handle:
            GenomeSeqFile_handle.write(handle.read())
        handle.close()

def main():
    parser = argparse.ArgumentParser(description="A script to download genbank records from a list of accessions")
    parser.add_argument('-i', '--id_list', required=True, help='A text file containing the list of IDs')
    parser.add_argument('-o', '--out_file', required=True, help='Output genbank.gb file name')
    parser.add_argument('-e', '--e_mail', required=True, help='Your E-mail address')
    parser.add_argument('-b', '--batch_size', required=True, type=int, help='Batch size')
    args = parser.parse_args()

    download_genbank_records(args.id_list, args.out_file, args.e_mail, args.batch_size)

if __name__ == "__main__":
    main()
