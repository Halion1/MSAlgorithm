# Import .py files

from classes import DistanceMatrixOwn
from itertools import combinations
from Bio import pairwise2


def compute_kmer_distance(seq1, seq2, k=3):
    """Compute distance based on shared k-mers."""

    def generate_kmers(seq, k):
        """Generate set of k-mers for a sequence."""
        return {seq[i:i + k] for i in range(len(seq) - k + 1)}

    kmers_seq1 = generate_kmers(seq1.data, k)
    kmers_seq2 = generate_kmers(seq2.data, k)

    shared_kmers = kmers_seq1.intersection(kmers_seq2)

    # Return the inverse of shared k-mer count; add a small constant to avoid division by zero
    return 1.0 / (len(shared_kmers) + 0.01)


def compute_distance_matrix(sequences):
    """Compute the distance matrix for a set of sequences using pairwise alignment scores."""
    n = len(sequences)
    matrix = DistanceMatrixOwn(n)

    # Loop over all unique pairs of sequences
    for i, seqA in enumerate(sequences):
        for j, seqB in enumerate(sequences):
            if j <= i:
                continue
            score = pairwise2.align.globalxx(seqA.seq, seqB.seq, one_alignment_only=True, score_only=True)

            # Use the negative of the score as the distance, as higher scores indicate more similarity
            distance = -score

            matrix.set_distance(i, j, distance)

    return matrix

