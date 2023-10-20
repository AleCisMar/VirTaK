# VirTaK (Virus Taxonomy from Kmers)
 
VirTaK is intended to quickly relate a virus with known taxonomy to a metagenomic assembled virus. It is based on a k-mer (default = 4-mer) profile database build from complete genomes within the ICTV Virus Metadata Resource (https://ictv.global/vmr). For each query assembly in fasta format it will calculate a k-mer (default = 4-mer) profile and search for the n (default = 10) viruses with the least dissimilar profiles (Bray-Curtis dissimilarity) in the kmer profile database.

## Dependencies

VirTaK codes where developed and tested with Python 3.9.18 using the following libraries:
* Biopython 1.78
* Matplotlib 3.7.2
* NumPy 1.26.0
* Pandas 2.0.3
* SciPy 1.11.3

## Installation

You can download the repository as a ZIP archive or clone the repository with:

```{bash, eval=FALSE, echo=TRUE}
git clone https://github.com/AleCisMar/VirTaK.git
``` 
Once unpacked, within the VirTak directory:
### Option 1 (recommended): create environment with all dependencies from environment.yml file:
```{bash, eval=FALSE, echo=TRUE}
conda env create -f environment.yml
```

### Option 2: create environment with all dependencies from a single command:

```{bash, eval=FALSE, echo=TRUE}
conda create -n VirTaK python=3.9 biopython matplotlib numpy pandas scipy
```

### Option 3: create a conda environment and install further dependencies one at a time:

```{bash, eval=FALSE, echo=TRUE}
conda create -n VirTaK python=3.9
```

```{bash, eval=FALSE, echo=TRUE}
conda activate VirTaK
```
Example:
```{bash, eval=FALSE, echo=TRUE}
conda install biopython=1.78
```

## Execution

Make sure to activate VirTaK environment before executing its codes:

```{bash, eval=FALSE, echo=TRUE}
conda activate VirTaK
```

VirTaK is made up of four main codes: build_kmer_database.py, kmer_search.py, calculate_threshold.py and generate_tree.py

### build_kmer_database.py (optional):

A 4-mer database is already provided in db/VMR_MSL38_v1_complete_genomes.kmers, along with db/DNA.kmers and db/RNA.kmers. Execute this code only if you want to build an updated database or if you want to build a database with different k-mer length.

```{bash, eval=FALSE, echo=TRUE}
usage: build_kmer_database.py [-h] -i INPUT_LIST -e EMAIL -o OUTPUT_PREFIX -a ANNOTATION_FILE [-k KMER_LENGTH]

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

The output consists of two files: database_name.fasta and database_name.kmers
NOTE: a copy of database_name.fasta is used by PanPhylo (available at https://drive.google.com/file/d/1C-zm-d9eTlarlSKQ_19eQTm-V7m_SNh_/view?usp=sharing). So, if we want to update the database used by PanPhylo we need to execute this code.

### kmer_search.py:

This is the code that will quickly relate a virus with known taxonomy to a metagenomic assemblied virus. It helps answering the question: to whom our assembled viruses are related to?

```{bash, eval=FALSE, echo=TRUE}
usage: kmer_search.py [-h] -d DATABASE [-k KMER_LENGTH] [-n N_MATCHES]

Calculates k-mer profiles for every file.fasta file in the current directory (creates file.kmer files) and searches for the n (default=10) viruses with the least dissimilar (Bray-Curtis) k-mer profiles
in the input database. The results are written to files named file_matches.txt with 'profile names' in the first column, and 'dissimilarity values' in the second column

optional arguments:
  -h, --help            show this help message and exit
  -d DATABASE, --database DATABASE
                        Input k-mer profile database. Example: db/VMR_MSL38_v1_complete_genomes.kmers
  -k KMER_LENGTH, --kmer_length KMER_LENGTH
                        OPTIONAL: Length of k-mers. Default = 4
  -n N_MATCHES, --n_matches N_MATCHES
                        OPTIONAL: Number of top matches to display in output file. Default = 10
```

Additionally to the output results file it will create .kmer files for each .fasta file in the current directory.

### calculate_threshold.py:

This is not compulsory but it will calculate the maximum distance between profiles in a specified taxon. If kmer_search.py found that the least dissimilar viruses to the query virus are coronaviruses, we can calculate the maximum distance between all profiles belonging to the family Coronaviridae or to the order Nidovirales, etc. which will help us decide whether the query virus actually belongs to the given taxon or not. It helps answering the question: do assembled viruses belong to the taxon of their least dissimilar virus?, or how close may the relationship be?

```{bash, eval=FALSE, echo=TRUE}
usage: calculate_threshold.py [-h] -i INPUT -d DATABASE -s STRING

Knowing the taxa associated to the least dissimilar viruses (in kmer_search.py output files), the user can calculate the threshold for a specific taxon (for example Coronaviridae). This code will search
for the kmer profiles that match the user defined string (Coronaviridae, Nidovirales or whatever) in the kmer profile database, calculate all vs all Bray-Curtis dissimilarities and get the maximum value
(threshold) for that specific taxon. Finally it will update the input file adding the columns 'taxon, threshold, and number of viruses' to the lines that match the string in the input file

optional arguments:
  -h, --help            show this help message and exit
  -i INPUT, --input INPUT
                        Input file. Is the output of kmer_search.py. Example: file_matches.txt
  -d DATABASE, --database DATABASE
                        Kmer profile database. Example: db/VMR_MSL38_v1_complete_genomes.kmers
  -s STRING, --string STRING
                        String to search in input file. Example: Coronaviridae, or Nidovirales, etc.
```

### generate_tree.py:

This is not compulsory but it will merge the profile of viruses listed in an input list with the profiles of all viruses belonging to a specified taxon, perform an all vs all comparison to obtain a distance matrix and finally a distance tree. Such tree will be a proxy of the phylogenetic relationships within the taxon and how is the query virus related to other member viruses. It helps answering the questions: how are our assembled viruses related to other viruses belonging to the taxon of the least dissimilar virus?, do they represent a novel lineage?

```{bash, eval=FALSE, echo=TRUE}
usage: generate_tree.py [-h] -d DATABASE -l LIST -s STRING -n NEWICK

Merge input kmer profiles with specified kmer profiles in the kmer database, calculate Bray-Curtis dissimilarity matrix, and compute a distance tree

optional arguments:
  -h, --help            show this help message and exit
  -d DATABASE, --database DATABASE
                        Kmer profile database
  -l LIST, --list LIST  File listing input ".kmer" files to process
  -s STRING, --string STRING
                        String to search in the kmer database. For example: Coronaviridae
  -n NEWICK, --newick NEWICK
                        Output Newick file```
