import argparse
import numpy as np

def read_data(input_file):
    data = {}
    ids = []
    with open(input_file, 'r') as file:
        lines = file.readlines()
        header = lines[0].strip().split('\t')[1:]  # Extract column headers
        for line in lines[1:]:
            fields = line.strip().split('\t')
            id = fields[0]
            scores = list(map(float, fields[1:]))
            data[id] = scores
            ids.append(id)
    return ids, header, data

def bray_curtis_distance(v1, v2):
    numerator = sum(abs(x - y) for x, y in zip(v1, v2))
    denominator = sum(v1) + sum(v2)
    return numerator / denominator

def calculate_bray_curtis_matrix(ids, data):
    num_samples = len(ids)
    matrix = np.zeros((num_samples, num_samples))

    for i in range(num_samples):
        for j in range(i, num_samples):
            distance = bray_curtis_distance(data[ids[i]], data[ids[j]])
            matrix[i][j] = distance
            matrix[j][i] = distance  # Matrix is symmetric, so set the mirrored value

    return matrix

def save_matrix(output_file, ids, matrix):
    with open(output_file, 'w') as file:
        file.write('\t' + '\t'.join(ids) + '\n')
        for i, row in enumerate(matrix):
            file.write(ids[i] + '\t' + '\t'.join(map(str, row)) + '\n')

def main():
    parser = argparse.ArgumentParser(description="Calculate Bray-Curtis dissimilarity matrix from a tab-delimited file")
    parser.add_argument("-i", "--input", required=True, help="Input tab-delimited file with scores")
    parser.add_argument("-o", "--output", required=True, help="Output file to save the Bray-Curtis matrix")

    args = parser.parse_args()

    ids, header, data = read_data(args.input)
    matrix = calculate_bray_curtis_matrix(ids, data)
    save_matrix(args.output, ids, matrix)

if __name__ == "__main__":
    main()
