# GuideTree.py

from Bio.Phylo.TreeConstruction import DistanceMatrix, DistanceTreeConstructor
import numpy as np  # Make sure numpy is imported, as it's needed for the isinstance check

def construct_guide_tree(distance_matrix, names):
    """
    Construct a guide tree using UPGMA algorithm from a distance matrix.

    :param distance_matrix: 2D numpy array, the matrix of pairwise distances.
    :param names: list of strings, names of the sequences corresponding to the distance matrix.
    :return: Bio.Phylo.BaseTree.Tree, the guide tree.
    """

    # Check if distance_matrix is a numpy array and convert it to a list of lists if it is
    if isinstance(distance_matrix, np.ndarray):
        distance_matrix = distance_matrix.tolist()

    # Add a print statement here to inspect the distance matrix
    print(distance_matrix)  # This will print the matrix to your console/terminal

    # Convert the list of lists to a Bio.Phylo.TreeConstruction.DistanceMatrix
    bio_distance_matrix = DistanceMatrix(names=names, matrix=distance_matrix)

    # Initialize a tree constructor.
    constructor = DistanceTreeConstructor()

    # Build a guide tree using the UPGMA algorithm.
    guide_tree = constructor.upgma(bio_distance_matrix)

    return guide_tree
