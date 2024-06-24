from Bio.Phylo.TreeConstruction import DistanceMatrix, DistanceTreeConstructor

# Convert our custom distance matrix to Biopython's format
def convert_to_lower_triangle(matrix, names):
    """Explicitly create the full lower triangle matrix."""
    n = len(names)
    lower_triangle = []
    for i in range(n):
        row = []
        for j in range(i + 1):
            row.append(matrix.matrix[i][j])
        lower_triangle.append(row)


    # Print sizes for debugging
    #print("Number of names:", len(names))
    #print("Size of biopy_matrix:", len(lower_triangle))
    #print("Biopy matrix content:")
    #for row in lower_triangle:
    #    print(row)

    return DistanceMatrix(names, lower_triangle)


# Generate guide tree using Biopython
def generate_guide_tree(matrix, sequence_names):
    biopy_matrix = convert_to_lower_triangle(matrix, sequence_names)
    constructor = DistanceTreeConstructor()
    tree = constructor.upgma(biopy_matrix)
    return tree