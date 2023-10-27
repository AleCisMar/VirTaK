# VirTaK (Virus Taxonomy from Kmers)
 
VirTaK is intended to quickly relate a virus with known taxonomy to a metagenomic assembled virus. It is based on a k-mer (default = 4-mer) profile database build from complete genomes within the ICTV Virus Metadata Resource (https://ictv.global/vmr). For each query genome in fasta format listed in the input file, it will calculate a k-mer (default = 4-mer) profile and search for the n (default = 10) viruses with the least dissimilar profiles (Bray-Curtis dissimilarity) in the kmer profile database. Then, it will compare the genome length, the number of coding sequences, and the similarity in domain content between each query genome and its corresponding most similar viruses to calculate an overall score. In general, the virus with the highest overall score will be more likely to actually be related to the query virus.

## Dependencies

VirTaK codes where developed and tested with Python 3.9.18 using the following libraries:

* Biopython 1.78
* Matplotlib 3.7.2
* NumPy 1.26.0
* Pandas 2.0.3
* SciPy 1.11.3
* Prodigal 2.6.3
* Hmmer 3.3.2
* Pfam_scan 1.6

## Installation and setup

### Download the repository as a ZIP archive or clone the repository with:

```{bash, eval=FALSE, echo=TRUE}
git clone https://github.com/AleCisMar/VirTaK.git
```

Once unpacked, within the VirTak directory: 

### Create VirTaK environment with all dependencies from environment.yml file:

```{bash, eval=FALSE, echo=TRUE}
conda env create -f environment.yml
```
### Download databases:

* The VMR complete genomes database is available at https://drive.google.com/file/d/1C-zm-d9eTlarlSKQ_19eQTm-V7m_SNh_/view?usp=sharing.
* Pfam-A database. Download and decompress both Pfam-A.hmm.gz and Pfam-A.hmm.dat.gz from https://www.ebi.ac.uk/interpro/download/pfam/.

### Create index files for HMM database

```{bash, eval=FALSE, echo=TRUE}
hmmpress Pfam-A.hmm
```

## Execution

Make sure to activate VirTaK environment before executing its codes:

```{bash, eval=FALSE, echo=TRUE}
conda activate VirTaK
```

### VirTaK.py:

Create a list of files.fasta to process. Example:

```{bash, eval=FALSE, echo=TRUE}
ls *.fasta > list.txt
```

This code helps answering the question: to whom our assembled viruses are related to? 

```{bash, eval=FALSE, echo=TRUE}
usage: VirTaK.py [-h] -l LIST -d DATABASE -f FASTA_DATABASE -p PFAM_PATH -o OUTPUT_FILE [-k KMER_LENGTH] [-n N_MATCHES]

Calculates k-mer profiles for every file.fasta listed in an input list and searches for the n (default=10) viruses with the least dissimilar (Bray-Curtis) k-mer profiles in the input database. Finally, to test if the most similar viruses are actually related to the query virus it compares the genome length, number of coding sequences, and the similarity in domain content between each query genome and its corresponding most similar viruses.

optional arguments:
  -h, --help            show this help message and exit
  -l LIST, --list LIST  Input list of .fasta files to process
  -d DATABASE, --database DATABASE
                        Input k-mer profile database. Example: db/VMR_MSL38_v1_complete_genomes.kmers
  -f FASTA_DATABASE, --fasta_database FASTA_DATABASE
                        Input fasta database. Example: db/VMR_MSL38_v1_complete_genomes.fasta
  -p PFAM_PATH, --pfam_path PFAM_PATH
                        Path to directory where Pfam-A.hmm database is located. Example: db/
  -o OUTPUT_FILE, --output_file OUTPUT_FILE
                        Output file name. Example: results.txt
  -k KMER_LENGTH, --kmer_length KMER_LENGTH
                        OPTIONAL: Length of k-mers. Default = 4
  -n N_MATCHES, --n_matches N_MATCHES
                        OPTIONAL: Number of top matches to display in output file. Default = 10
```

### build_virtak_database.py (optional):

A 4-mer database is already provided in db/VMR_MSL38_v1_complete_genomes.kmers, along with db/DNA.kmers and db/RNA.kmers. Execute this code only if you want to build an updated database (from an updated VMR) or if you want to build a database with different k-mer length. If the download process is unexpectedly interrupted we can resume the process by executing the same command.

```{bash, eval=FALSE, echo=TRUE}
usage: build_virtak_database.py [-h] -i INPUT_LIST -e EMAIL -o OUTPUT_PREFIX -a ANNOTATION_FILE [-k KMER_LENGTH]

Download and process GenBank records from a list of accessions to obtain complete genomes in fasta format and create a kmer profile database. Will also generate a fasta file with complete genomes and
separate .kmers files for DNA and RNA viruses.

optional arguments:
  -h, --help            show this help message and exit
  -i INPUT_LIST, --input_list INPUT_LIST
                        Input file containing a list of GenBank accessions. Each line may have one or more accessions (segmented genomes) separated by ;
  -e EMAIL, --email EMAIL
                        Your email address for NCBI Entrez utilities
  -o OUTPUT_PREFIX, --output_prefix OUTPUT_PREFIX
                        Prefix for the output files. Example: VMR_MSL38_v1_complete_genomes
  -a ANNOTATION_FILE, --annotation_file ANNOTATION_FILE
                        Annotation file needed to name fasta headers. Example: lists/VMR_MSL38_v1_complete_genomes.txt
  -k KMER_LENGTH, --kmer_length KMER_LENGTH
                        OPTIONAL: Length of k-mers. Default = 4
```

The output consists of two main files (database_name.fasta and database_name.kmers) and two additional files (DNA.kmers and RNA.kmers).
