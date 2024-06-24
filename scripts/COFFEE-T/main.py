from Bio import SeqIO, Phylo
from Bio.Align import MultipleSeqAlignment
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from PairwiseAlignmentScore import pairwise_alignment_score
from GuideTree import construct_guide_tree
from DistanceMatrix import build_distance_matrix, lower_triangle_matrix
from ProgressiveAlign import progressive_alignment
from checkseq import check_sequences, pad_sequences


substitution_matrix = {
    ('A', 'A'): 1, ('A', 'T'): -1, ('A', 'C'): -1, ('A', 'G'): -1, ('A', '-'): -2,
    ('T', 'A'): -1, ('T', 'T'): 1, ('T', 'C'): -1, ('T', 'G'): -1, ('T', '-'): -2,
    # ... and so on for all pairs
}

if __name__ == "__main__":
    # Sample sequences; in practice, these could be read from a file.
    sequences = ["ACGGGT", "ACGTTT", "ACGG", "ACGT", "ACGTTTT", "ACGGGG"]

    # Extract unique characters from the sequences for the substitution matrix.
    characters = set(''.join(sequences))

    # Define a substitution matrix.
    substitution_matrix = {(x, y): +1 if x == y else -1 for x in characters for y in characters}

    # Define a gap penalty.
    gap_penalty = -1  # Adjust based on needs.

    # Build distance matrix.
    distance_matrix = build_distance_matrix(sequences, substitution_matrix, gap_penalty)

    # Convert to a lower triangular matrix
    lower_distance_matrix = lower_triangle_matrix(distance_matrix)
    print(lower_distance_matrix)

    # Specify names for sequences.
    names = [f"Seq{i}" for i in range(1, len(sequences) + 1)]

    guide_tree = construct_guide_tree(lower_distance_matrix, names)
    Phylo.draw_ascii(guide_tree)

    # Align sequences using the guide tree
    result_alignment = progressive_alignment(guide_tree, substitution_matrix, gap_penalty, sequences, names=names)

    print("Aligned sequences:", result_alignment)

    # Check if all sequences are the same length and pad them if necessary
    print("Before padding:", result_alignment)
    aligned_seqs = pad_sequences(result_alignment)
    print("After padding:", aligned_seqs)

    # Create SeqRecord objects from the aligned sequences
    msa_records = [SeqRecord(Seq(seq.seq), id=seq.id) for seq in aligned_seqs]

    msa = MultipleSeqAlignment(msa_records)
    print(msa)
