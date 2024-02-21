Seven files can be found:
* counts.txt - Number of times each domain or unnanotated cluster is found in every genome
* distance_matrix.csv - Distance matrix. Built from counts.txt. Used for distance tree construction
* results_nj.tree - Neighbor joining tree in newick format. Built from distance_matrix.csv. 
* reordered_counts.csv - rows in counts.txt are reordered based on results_nj.tree
* reordered_distance_matrix.csv - rows and columns in distance_matrix.csv are reordered based on results_nj.tree. Can be used to see the distance between pairs of genomes
* transposed_counts.txt - transposed version of reordered_counts.csv. Used to cluster domains with euclidean distance
* reordered_transposed_counts.csv - rows in transposed_counts.txt are reordered based on euclidean distance. Can be used to visualize the similarities in domain content
* domains.txt - list of domains identified in the dataset. Can be used with get_domain_info.py and get_gos.sh

To get annotations related to the domains identified in the dataset we can execute the following command:
```{bash, eval=FALSE, echo=TRUE}
get_domain_info.py domains.txt ~/db/virtak_db/Pfam-A.seed
```
Make sure to use the correct paths to Pfam-A.seed. This will create a tab delimited file named domain_info.tsv

To get Gene Ontologies associated to the Pfam accessions we can execute the following command:
```{bash, eval=FALSE, echo=TRUE}
get_gos.sh domain_info.tsv ~/db/virtak_db/pfam_gomf_most_specific.txt ~/db/virtak_db/go.txt
```
This will create a tab delimited file named domain_gos.tsv
