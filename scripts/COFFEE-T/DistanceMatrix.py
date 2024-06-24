import numpy as np
from PairwiseAlignmentScore import  needleman_wunsch, create_substitution_matrix

def lower_triangle_matrix(full_matrix):
    """
    Convert a full distance matrix to a lower triangular matrix.

    :param full_matrix: 2D list (full distance matrix).
    :return: 2D list (lower triangular matrix).
    """
    lower_matrix = []
    for i in range(len(full_matrix)):
        row = full_matrix[i][:i+1]
        lower_matrix.append(row)
    return lower_matrix

def compute_distance(seq1, seq2, substitution_matrix, gap_penalty):
    # Get the score and aligned sequences from the Needleman-Wunsch algorithm
    score, align1, align2 = needleman_wunsch(seq1, seq2, substitution_matrix, gap_penalty)

    # Calculate the number of mismatches and gaps
    mismatches = sum(1 for a, b in zip(align1, align2) if a != b)
    gaps = align1.count('-') + align2.count('-')

    # Define the distance as the number of mismatches and gaps (or some other function)
    distance = mismatches + gaps

    return distance

def build_distance_matrix(sequences, substitution_matrix, gap_penalty):
    """
    Calculate the distance matrix for a set of sequences.

    :param sequences: list of sequences to be compared.
    :param substitution_matrix: dictionary representing the substitution matrix.
    :param gap_penalty: integer representing the gap penalty.
    :return: 2D list (distance matrix).
    """
    num_sequences = len(sequences)
    distance_matrix = []

    for i in range(num_sequences):
        row = []
        for j in range(num_sequences):
            if i == j:
                row.append(0)  # The distance between a sequence and itself is zero
            else:
                # Now we are passing all required arguments to the function
                distance = compute_distance(sequences[i], sequences[j], substitution_matrix, gap_penalty)
                row.append(distance)
        distance_matrix.append(row)

    return distance_matrix

