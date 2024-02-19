# VirTaK (Virus Taxonomy from Kmers)
 
VirTaK is intended to relate a virus with known taxonomy to a metagenomic assembled virus. It is based on a k-mer (default = 4-mer) profile database build from complete genomes within the ICTV Virus Metadata Resource (https://ictv.global/vmr). For each query genome in fasta format listed in the input file, it will calculate a k-mer (default = 4-mer) profile and compare it to the k-mer profile of viruses in the k-mer profile database. Then, it will evaluate nucleotide, amino acid, domain content, number of coding sequences, and genome length similarities to calculate an overall score. In general, the virus with the highest overall score will be more likely to actually be related to the query virus.

## Dependencies

VirTaK codes where developed and tested with Python 3.9.18 using the following libraries:

* Biopython 1.78
* Blast 2.12.0
* Hmmer 3.3.2
* NumPy 1.26.0
* Pandas 2.1.1
* Pfam_scan 1.6
* Prodigal 2.6.3
* SciPy 1.11.3

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
Installing pfam_scan may be problematic in some systems. If this is the case try removing pfam_scan from the environment.yml file and execute the above code. Then, activate the created environment and try:

```{bash, eval=FALSE, echo=TRUE}
conda install pfam_scan --channel Bioconda
```

### Download VirTaK database:

* The VMR complete genomes database is available at https://drive.google.com/file/d/1TqNQsjCyDSvLUAGJufLoHvBxwd91n04T/view?usp=sharing.
This is a compressed directory that contains all files needed to execute VirTaK.py and build_virtak_database.py. It includes the Pfam-A database.

## Execution

Make sure to activate VirTaK environment before executing its codes:

```{bash, eval=FALSE, echo=TRUE}
conda activate VirTaK
```

Also, it may be necessary to activate execution permissions:

```{bash, eval=FALSE, echo=TRUE}
chmod 755 /path/to/VirTaK.py
```

Additionally it may be confortable to create a symbolic link within the bin directory:

If bin directory does not exist:
```{bash, eval=FALSE, echo=TRUE}
cd
mkdir bin
```
Then:
```{bash, eval=FALSE, echo=TRUE}
cd bin
ln -s /path/to/VirTaK.py .
```

### VirTaK.py:

Create a list of files.fasta to process. Example:

```{bash, eval=FALSE, echo=TRUE}
ls *.fasta > list.txt
```

This code helps answering the question: to whom our assembled viruses are related to? 

```{bash, eval=FALSE, echo=TRUE}
usage: VirTaK.py [-h] -l LIST -d DATABASE -p PFAM_PATH -o OUTPUT_FILE [-k KMER_LENGTH] [-n N_MATCHES]

Calculates k-mer profiles for every fasta file in input list, and compares them against the k-mer profiles in the input database. It also evaluates nucleotide, amino acid, domain content, number of
protein coding sequences, and genome length similarities to compute a relatedness score

optional arguments:
  -h, --help            show this help message and exit
  -l LIST, --list LIST  Input list of .fasta files to process
  -d DATABASE, --database DATABASE
                        Path to directory where VirTaK database is located. Must include database name without extension. Example: db/VMR_MSL38_v1_complete_genomes
  -p PFAM_PATH, --pfam_path PFAM_PATH
                        Path to directory where Pfam-A.hmm database is located. Example: db/
  -o OUTPUT_FILE, --output_file OUTPUT_FILE
                        Output file name. Example: results.txt
  -k KMER_LENGTH, --kmer_length KMER_LENGTH
                        OPTIONAL: Length of k-mers. Default = 4
  -n N_MATCHES, --n_matches N_MATCHES
                        OPTIONAL: Number of top matches to display in output file. Default = 10
```

### PanPhylo.py:

This code helps answering: how are query viruses related to other known viruses?, are they novel endemic viruses?, how many genes they share with their
closest relatives? which genes support such relationships?

```{bash, eval=FALSE, echo=TRUE}
usage: PanPhylo.py [-h] -l LIST -d DATABASE [-s STRINGS_FILE] [-o OUTPUT_DIR] [-t THREADS]

PanPhylo (Pangenomic and Phylogenomic analysis) is intended to reconstruct the phylogenomic relationships between a group of assembled viruses (metagenomic) and a group of known viruses selected from previous knowledge
of taxonomic membership. It takes the same input list.txt as VirTaK.py, aswell as the VirTaK database, and some files created during the taxonomic classification process. It creates a domain (and proteins without domain
annotation) content table which is used to calculate a Bray-Curtis distance matrix. Distances are represented by a distance tree built with the neighbor joining algorithm

optional arguments:
  -h, --help            show this help message and exit
  -l LIST, --list LIST  Input list of query .fasta files
  -d DATABASE, --database DATABASE
                        Path to directory where VirTaK database is located. Must include database name without extension. Example: db/VMR_MSL38_v1_complete_genomes
  -s STRINGS_FILE, --strings_file STRINGS_FILE
                        Optional. File containing a list of strings to search in the VirTaK database. For example: Coronaviridae
  -o OUTPUT_DIR, --output_dir OUTPUT_DIR
                        Optional. Output directory name
  -t THREADS, --threads THREADS
                        Number of parallel CPU workers to use for HMMER (default=2)
```

### build_virtak_database.py (optional):

Execute this code only if you want to build an updated database (from an updated VMR) or if you want to build a database with different k-mer length. If the download process is unexpectedly interrupted we can resume the process by executing the same command.

```{bash, eval=FALSE, echo=TRUE}
usage: build_virtak_database.py [-h] -i INPUT_LIST -o OUTPUT_PREFIX -a ANNOTATION_FILE -p PFAM_PATH -e EMAIL [-k KMER_LENGTH]

Create VirTaK database. From input file containing a list of GenBank accessions downloads and process fasta genomes to: create blastn database, predict protein coding sequences, create blastp database,
create kmer profile database, and create a protein domain database

optional arguments:
  -h, --help            show this help message and exit
  -i INPUT_LIST, --input_list INPUT_LIST
                        Input file containing a list of GenBank accessions. Each line may have one or more accessions (segmented genomes) separated by ;
  -o OUTPUT_PREFIX, --output_prefix OUTPUT_PREFIX
                        Prefix for the output files. Example: VMR_MSL38_v1_complete_genomes
  -a ANNOTATION_FILE, --annotation_file ANNOTATION_FILE
                        Annotation file needed to name fasta headers. Example: lists/VMR_MSL38_v1_complete_genomes.txt
  -p PFAM_PATH, --pfam_path PFAM_PATH
                        Path to directory where Pfam-A.hmm database is located. Example: db/
  -e EMAIL, --email EMAIL
                        Your email address for NCBI Entrez utilities
  -k KMER_LENGTH, --kmer_length KMER_LENGTH
                        OPTIONAL: Length of k-mers. Default = 4
```

The output consists of four main files (database_name.fasta, database_name.faa, database_name.kmers, and database_name.pfamscan) and associated blast database files.
