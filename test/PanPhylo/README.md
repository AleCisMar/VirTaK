Seven files can be found:
* counts.txt - Number of times each domain or unnanotated cluster is found in every genome
* distance_matrix.csv - Distance matrix. Built from counts.txt. Used for distance tree construction
* results_nj.tree - Neighbor joining tree in newick format. Built from distance_matrix.csv. 
* reordered_counts.csv - rows in counts.txt are reordered based on results_nj.tree
* reordered_distance_matrix.csv - rows and columns in distance_matrix.csv are reordered based on results_nj.tree. Can be used to see the distance between pairs of genomes
* transposed_counts.txt - transposed version of reordered_counts.csv. Used to cluster domains with euclidean distance
* reordered_transposed_counts.csv - rows in transposed_counts.txt are reordered based on euclidean distance. Can be used to visualize the similarities in domain content
