import argparse
import pandas as pd
from scipy.spatial import distance
from scipy.cluster import hierarchy

def parse_args():
    parser = argparse.ArgumentParser(description="Compute a distance tree from a symmetric distance matrix.")
    parser.add_argument("-i", "--input", required=True, help="Input distance matrix file")
    parser.add_argument("-o", "--output", required=True, help="Output tree file in Newick format")
    return parser.parse_args()

def read_distance_matrix(input_file):
    df = pd.read_csv(input_file, sep='\t', index_col=0)
    return df

def compute_distance_tree(distance_matrix):
    condensed_distance = distance.squareform(distance_matrix.values)
    linkage = hierarchy.linkage(condensed_distance, method='average')
    dendrogram = hierarchy.dendrogram(linkage, labels=distance_matrix.index, orientation='right', leaf_font_size=8)
    return linkage, dendrogram

def save_dendrogram_to_newick(linkage, dendrogram, output_file):
    def generate_newick(node_id):
        if node_id < len(dendrogram['leaves']):
            return str(dendrogram['ivl'][node_id])
        else:
            left_node = generate_newick(int(linkage[node_id - len(dendrogram['leaves']), 0]))
            right_node = generate_newick(int(linkage[node_id - len(dendrogram['leaves']), 1]))
            return f"({left_node},{right_node})"

    tree_str = generate_newick(len(linkage) + len(dendrogram['leaves']) - 1)
    
    with open(output_file, "w") as f:
        f.write(tree_str + ";")

def main():
    args = parse_args()
    input_file = args.input
    output_file = args.output

    # Step 1: Read the symmetric distance matrix
    distance_matrix = read_distance_matrix(input_file)

    # Step 2: Compute the distance tree
    linkage, dendrogram = compute_distance_tree(distance_matrix)

    # Step 3: Save the dendrogram to Newick format
    save_dendrogram_to_newick(linkage, dendrogram, output_file)

if __name__ == "__main__":
    main()
